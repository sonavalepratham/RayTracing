#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Light { // Light source
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    Vec3f position; // position in space
    float intensity; // intensity
};

struct Material { // body material, contains color
    Material(const Vec2f &a, const Vec3f &color, const float &spec) : albedo(a), diffuse_color(color), specular_exponent(spec) {}
    Material() : albedo(1,0), diffuse_color(), specular_exponent() {}
    // albedo is diffused reflectance or reflectivity. it's value is from 0 to 1.
    Vec2f albedo; // array of 2 elements. 1st diffused albedo, 2nd specular albedo
    //Light is scattered in many directions when it strikes a rough surface
    Vec3f diffuse_color;
    // specular_exponent(shininess), infinite for mirror. 0 for diffused material.
    float specular_exponent;
};

Vec3f reflect(const Vec3f &I, const Vec3f &N) { // I is incoming ray on surface(invert), N is normal on surface.
/*

returns reflection vector of Incident ray from surface in reverse direction.
R+I = (2*N.I)*N
R = (2*N.I)*N - I // reflection vector

*/
    return I - N*2.f*(I*N);
}

struct Sphere { // 3D sphere
    Vec3f center;
    float radius;
    Material material;

    Sphere(const Vec3f &c, const float &r, const Material &m) : center(c), radius(r), material(m) {}

/*
    returns true if intersecting else false.
    t0 gives distance between origin to intersection.
*/
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const { // orig - from which ray is coming(camera), dir - direction of ray, t0 - distance between center and ray
        Vec3f L = center - orig; // vector on origin to center.
        float tca = L*dir;  // dot product, distance between origin & center to direction intersection.
        float d2 = L*L - tca*tca; // distance between center to direction intersection
        if (d2 > radius*radius) return false; // if distance between direction & center is greater than radius return false (not intersecting)
        float thc = sqrtf(radius*radius - d2); // circle edge to center-direction intersection distance.
        t0       = tca - thc;  // t0 & t1 is distance between circle edge to origin.
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false; // source is inside sphere
        return true;
    }
};

/*
    returns true if ray intersect to any of sphere.
    output variables:
    hit: point where ray hit to sphere.
    N: normal from center to hit.
    material: material of intersecting sphere.
*/
bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, Vec3f &hit, Vec3f &N, Material &material) {
    float spheres_dist = std::numeric_limits<float>::max();
    for (size_t i=0; i < spheres.size(); i++) {
        float dist_i;
        if (spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist) { // if distance is less than previous sphere intersection it means current sphere is near to camera
            spheres_dist = dist_i;
            hit = orig + dir*dist_i;
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    return spheres_dist!=std::numeric_limits<float>::max(); // it means ray intersect to any of the sphere.
}

/*
    return color based on intersect or not.
*/
Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    Vec3f point, N; // point: intersection of ray on sphere. N: normal on sphere at intersection point.
    Material material; // material of sphere.
    if (!scene_intersect(orig, dir, spheres, point, N, material)) {
        return Vec3f(0.2, 0.7, 0.8); // background color, if not intersect any of sphere.
    }

    float diffuse_light_intensity = 0, specular_light_intensity = 0;
    for (size_t i=0; i<lights.size(); i++) {
        Vec3f light_dir      = (lights[i].position - point).normalize();
        diffuse_light_intensity  += lights[i].intensity * std::max(0.f, light_dir*N); // light intensity will be higher if N and light_dir are parallel.
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N)*dir), material.specular_exponent)*lights[i].intensity;
        // dot product of parallel vector is multiplication of magnitude, and perpendicular is 0.
    }
    return material.diffuse_color * diffuse_light_intensity * material.albedo[0] + Vec3f(1., 1., 1.)*specular_light_intensity * material.albedo[1];
}

void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width    = 1024;
    const int height   = 768;
    const int fov      = M_PI/2.;
    std::vector<Vec3f> framebuffer(width*height);

    #pragma omp parallel for
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height; // x co ordinate with aspect ratio.
            float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.); // y co ordinate, negative sign applied as height is decreasing with increase of angle
            Vec3f dir = Vec3f(x, y, -1).normalize(); // z is -1 as frame is at 1 unit distance from camera & on opposite side(inside verticle plane).
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, spheres, lights); // cast ray and store background in frame buffer.
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary); // open in binary mode
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max>1) c = c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j]))); // converting in range 0-255
        }
    }
    ofs.close();
}

int main() {
    Material      ivory(Vec2f(0.6,  0.3), Vec3f(0.4, 0.4, 0.3),   50.);
    Material red_rubber(Vec2f(0.9,  0.1), Vec3f(0.3, 0.1, 0.1),   10.);

    std::vector<Sphere> spheres;

    spheres.push_back(Sphere(Vec3f(-3,    0,   -16), 2,      ivory));
    spheres.push_back(Sphere(Vec3f(-1.0, -1.5, -12), 2, red_rubber));
    spheres.push_back(Sphere(Vec3f( 1.5, -0.5, -18), 3, red_rubber));
    spheres.push_back(Sphere(Vec3f( 7,    5,   -18), 4,      ivory));
    
    std::vector<Light>  lights;
    lights.push_back(Light(Vec3f(-20, 20,  20), 1.5));
    lights.push_back(Light(Vec3f( 30, 50, -25), 1.8));
    lights.push_back(Light(Vec3f( 30, 20,  30), 1.7));

    render(spheres, lights);

    return 0;
}
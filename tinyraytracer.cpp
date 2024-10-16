#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "geometry.h"

struct Sphere { // 3D sphere
    Vec3f center;
    float radius;

    Sphere(const Vec3f &c, const float &r) : center(c), radius(r) {}

/*
    returns true if intersecting else false.
*/
    bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const { // orig - from which ray is coming(camera), dir - direction of ray, t0 - distance between center and ray
        Vec3f L = center - orig; // vector on origin to center.
        float tca = L*dir;  // dot product, distance between origin & center to direction intersection.
        float d2 = L*L - tca*tca; // distance between center to direction intersection
        if (d2 > radius*radius) return false; // if distance between direction & center is greater than radius return false (not intersecting)
        float thc = sqrtf(radius*radius - d2); // circle edge to center-direction intersection distance.
        t0       = tca - thc; 
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }
};

/*
    return color based on intersect or not.
*/

Vec3f cast_ray(const Vec3f &orig, const Vec3f &dir, const Sphere &sphere) {
    float sphere_dist = std::numeric_limits<float>::max();
    if (!sphere.ray_intersect(orig, dir, sphere_dist)) {
        return Vec3f(0.2, 0.7, 0.8); // background color
    }
    return Vec3f(0.4, 0.4, 0.3);
}

void render(const Sphere &sphere) {
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
            framebuffer[i+j*width] = cast_ray(Vec3f(0,0,0), dir, sphere); // cast ray and store background in frame buffer.
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm", std::ofstream::out | std::ofstream::binary); // open in binary mode
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j]))); // converting in range 0-255
        }
    }
    ofs.close();
}

int main() {
    Sphere sphere(Vec3f(-3, 0, -16), 2);
    render(sphere);

    return 0;
}
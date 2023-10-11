#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>

struct Vec3 {
    double x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}
    Vec3 operator+(const Vec3 &b) const { return Vec3(x + b.x, y + b.y, z + b.z); }
    Vec3 operator-(const Vec3 &b) const { return Vec3(x - b.x, y - b.y, z - b.z); }
    Vec3 operator*(double b) const { return Vec3(x * b, y * b, z * b); }
    double dot(const Vec3 &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec3 cross(const Vec3 &b) const { return Vec3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
    Vec3 normalize() const {
        double norm = std::sqrt(dot(*this));
        return Vec3(x / norm, y / norm, z / norm);
    }
};

struct Ray {
    Vec3 origin, direction;
    Ray(const Vec3 &origin, const Vec3 &direction) : origin(origin), direction(direction) {}
};

struct Sphere {
    Vec3 center;
    double radius;
    Vec3 color;
    Sphere(const Vec3 &center, double radius, const Vec3 &color) : center(center), radius(radius), color(color) {}
    
    bool intersect(const Ray &ray, double &t) const {
        Vec3 oc = ray.origin - center;
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * c;
        if (discriminant > 0) {
            t = (-b - sqrt(discriminant)) / 2.0;
            if (t > 0.0) return true;
            t = (-b + sqrt(discriminant)) / 2.0;
            return t > 0.0;
        }
        return false;
    }
};

Vec3 trace(const Ray &ray, const Sphere &sphere) {
    double t = 0;
    if (!sphere.intersect(ray, t)) return Vec3(0, 0, 0);  // Missed sphere, return black color.
    return sphere.color;  // Hit sphere, return its color.
}

int main() {
    const int width = 512;
    const int height = 512;
    Sphere sphere(Vec3(0, 0, -5), 1, Vec3(0.5, 0.5, 0.5));

    std::vector<Vec3> framebuffer(width * height);

    #pragma omp parallel for
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            double x = (2.0 * i - width) / height;
            double y = (2.0 * j - height) / height;
            Ray ray(Vec3(0, 0, 0), Vec3(x, y, -1).normalize());
            framebuffer[i + j * width] = trace(ray, sphere);
        }
    }

    std::cout << "P3\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < width * height; ++i) {
        std::cout << (int)(255 * framebuffer[i].x) << " " << (int)(255 * framebuffer[i].y) << " " << (int)(255 * framebuffer[i].z) << "\n";
    }

    return 0;
}

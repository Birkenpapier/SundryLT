#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(0.0, 1.0);

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

double clamp(double val, double min_val, double max_val) {
    return std::max(min_val, std::min(val, max_val));
}

struct Ray {
    Vec3 origin, direction;
    Ray(const Vec3 &origin, const Vec3 &direction) : origin(origin), direction(direction) {}
};

struct Material {
    Vec3 color;
    double albedo;    
    bool isMirror;    
    Material(const Vec3 &color, double albedo, bool isMirror) : color(color), albedo(albedo), isMirror(isMirror) {}
};

struct Sphere {
    Vec3 center;
    double radius;
    Material material;

    Sphere(const Vec3 &center, double radius, const Material &material) : center(center), radius(radius), material(material) {}

    bool intersect(const Ray &ray, double &t) const {
        Vec3 oc = ray.origin - center;
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * c;
        if (discriminant > 0) {
            t = (-b - std::sqrt(discriminant)) / 2.0;
            if (t > 0.0) return true;
            t = (-b + std::sqrt(discriminant)) / 2.0;
            return t > 0.0;
        }
        return false;
    }
};

Vec3 randomPointInHemisphere(const Vec3& normal) {
    Vec3 inSphere;
    do {
        inSphere = Vec3(distribution(generator) * 2 - 1, distribution(generator) * 2 - 1, distribution(generator) * 2 - 1);
    } while (inSphere.dot(inSphere) > 1.0 || inSphere.dot(normal) < 0);
    return inSphere.normalize();
}

Vec3 trace(const Ray &ray, const Sphere &sphere, int depth = 0) {
    const Vec3 lightPos(5, 5, 0);
    const Vec3 lightColor(1, 1, 1);
    const double lightIntensity = 5.0;

    double t = 0;
    if (!sphere.intersect(ray, t) || depth > 2) {
        return Vec3(0, 0, 0); 
    }

    Vec3 hitPoint = ray.origin + ray.direction * t;
    Vec3 normal = (hitPoint - sphere.center).normalize();

    Vec3 toLight = (lightPos - hitPoint).normalize();
    double diff = std::max(0.0, normal.dot(toLight));
    Vec3 diffuseColor = sphere.material.color * (diff * lightIntensity * sphere.material.albedo);

    Vec3 reflectedColor(0, 0, 0);
    if (sphere.material.isMirror) {
        Vec3 reflectedDir = ray.direction - normal * 2.0 * ray.direction.dot(normal);
        Ray reflectedRay(hitPoint + reflectedDir * 1e-3, reflectedDir);
        reflectedColor = trace(reflectedRay, sphere, depth + 1);
    }

    return diffuseColor + reflectedColor;
}

int main() {
    const int width = 512;
    const int height = 512;
    Sphere sphere(Vec3(0, 0, -5), 1, Material(Vec3(0.5, 0.5, 0.5), 0.8, true));

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
        std::cout << (int)(255 * clamp(framebuffer[i].x, 0.0, 1.0)) << " "
                  << (int)(255 * clamp(framebuffer[i].y, 0.0, 1.0)) << " "
                  << (int)(255 * clamp(framebuffer[i].z, 0.0, 1.0)) << "\n";
    }

    return 0;
}

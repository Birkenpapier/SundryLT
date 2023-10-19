#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

// Vec3 class to represent vectors and points
struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0) {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z) {}

    Vec3 operator+(const Vec3 &other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator-(const Vec3 &other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    Vec3 operator/(double scalar) const {
        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    double dot(const Vec3 &other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    Vec3 cross(const Vec3 &other) const {
        return Vec3(y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x);
    }

    Vec3 normalize() const {
        double magnitude = std::sqrt(x * x + y * y + z * z);
        return Vec3(x / magnitude, y / magnitude, z / magnitude);
    }
};

// Ray class to represent rays
struct Ray {
    Vec3 origin, direction;
    Ray(const Vec3 &origin, const Vec3 &direction) : origin(origin), direction(direction.normalize()) {}
};

// Material class to define material properties of objects
struct Material {
    Vec3 color;
    double albedo;    
    bool isMirror;    
    bool isTransparent;
    double refractiveIndex;
    
    Material(const Vec3 &color, double albedo, bool isMirror, bool isTransparent = false, double refractiveIndex = 1)
        : color(color), albedo(albedo), isMirror(isMirror), isTransparent(isTransparent), refractiveIndex(refractiveIndex) {}
};

// Sphere class to represent spherical objects
struct Sphere {
    Vec3 center;
    double radius;
    Material material;

    Sphere(const Vec3 &center, double radius, const Material &material)
        : center(center), radius(radius), material(material) {}

    bool intersect(const Ray &ray, double &t) const {
        Vec3 oc = ray.origin - center;
        double a = ray.direction.dot(ray.direction);
        double b = 2.0 * oc.dot(ray.direction);
        double c = oc.dot(oc) - radius * radius;
        double discriminant = b * b - 4 * a * c;
        if (discriminant < 0) return false;
        t = (-b - sqrt(discriminant)) / (2.0 * a);
        return t > 0;
    }
};

double clamp(double value, double min, double max) {
    if (value < min) return min;
    if (value > max) return max;
    return value;
}

double fresnel(const Vec3 &I, const Vec3 &N, const double &ior) {
    double cosi = clamp(I.dot(N), -1.0, 1.0);
    double etai = 1, etat = ior;
    if (cosi > 0) { std::swap(etai, etat); }
    double sint = etai / etat * sqrt(std::max(0.0, 1 - cosi * cosi));
    if (sint >= 1) {
        return 1;
    } else {
        double cost = sqrt(std::max(0.0, 1 - sint * sint));
        cosi = fabs(cosi);
        double Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
        double Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
        return (Rs * Rs + Rp * Rp) / 2;
    }
}

Vec3 refract(const Vec3 &I, const Vec3 &N, const double &ior) {
    double cosi = clamp(I.dot(N), -1.0, 1.0);
    double etai = 1, etat = ior;
    Vec3 n = N;
    if (cosi < 0) {
        cosi = -cosi;
    } else {
        std::swap(etai, etat);
        n = N * -1;
    }
    double eta = etai / etat;
    double k = 1 - eta * eta * (1 - cosi * cosi);
    return k < 0 ? Vec3(0, 0, 0) : I * eta + n * (eta * cosi - sqrt(k));
}

Vec3 trace(const Ray &ray, const Sphere &sphere, int depth = 0) {
    if (depth > 4) return Vec3(); // prevent deep recursion

    double t;
    if (sphere.intersect(ray, t)) {
        Vec3 hitPoint = ray.origin + ray.direction * t;
        Vec3 normal = (hitPoint - sphere.center).normalize();
        Vec3 toLight = (Vec3(0, 0, 10) - hitPoint).normalize();

        double lightIntensity = std::max(0.0, normal.dot(toLight));

        if (sphere.material.isTransparent) {
            double fresnelEffect = fresnel(ray.direction, normal, sphere.material.refractiveIndex);
            bool inside = ray.direction.dot(normal) > 0;
            Vec3 refractionColor(0, 0, 0);

            if (fresnelEffect < 1) {
                Vec3 refractionDir = refract(ray.direction, normal, sphere.material.refractiveIndex);
                Ray refractionRay(hitPoint - normal * 1e-3, refractionDir);
                refractionColor = trace(refractionRay, sphere, depth + 1);
            }

            Vec3 reflectionDir = ray.direction - normal * 2.0 * ray.direction.dot(normal);
            Ray reflectionRay(hitPoint + normal * 1e-3, reflectionDir);
            Vec3 reflectionColor = trace(reflectionRay, sphere, depth + 1);

            return reflectionColor * fresnelEffect + refractionColor * (1 - fresnelEffect);
        }

        Vec3 color = sphere.material.color * lightIntensity * sphere.material.albedo;
        return color;
    }

    return Vec3(); // black background
}

int main() {
    const int width = 800;
    const int height = 600;
    const double aspectRatio = width / static_cast<double>(height);
    const Vec3 cameraPos(0, 0, 0);
    std::vector<Vec3> framebuffer(width * height);

    Sphere sphere(Vec3(0, 0, -5), 1, Material(Vec3(0.8, 0.8, 0.8), 0.8, false, true, 1.5));

    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            double x = (2 * (i + 0.5) / (double)width - 1) * tan(30.0 / 2.0) * aspectRatio;
            double y = -(2 * (j + 0.5) / (double)height - 1) * tan(30.0 / 2.0);
            Vec3 dir = Vec3(x, y, -1).normalize();
            Ray ray(cameraPos, dir);
            framebuffer[i + j * width] = trace(ray, sphere);
        }
    }

    std::cout << "P3\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height * width; i++) {
        Vec3 &c = framebuffer[i];
        std::cout << (int)(255 * clamp(c.x, 0.0, 1.0)) << " "
                  << (int)(255 * clamp(c.y, 0.0, 1.0)) << " "
                  << (int)(255 * clamp(c.z, 0.0, 1.0)) << "\n";
    }

    return 0;
}

// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

// Function to get the intersection point of a ray with a sphere
Vector3d getPointSphere(const Vector3d& sphere_center, const Vector3d& ray_origin, const Vector3d& ray_direction, double sphere_radius)
{
    double A = ray_direction.dot(ray_direction);
    double B = 2.0 * ray_direction.dot(ray_origin - sphere_center);
    double C = (ray_origin - sphere_center).dot(ray_origin - sphere_center) - (sphere_radius * sphere_radius);

    double discriminant = (B * B) - (4 * A * C);
    double t;

    if (discriminant >= 0)
    {
        double t1 = (-B + sqrt(discriminant)) / (2.0 * A);
        double t2 = (-B - sqrt(discriminant)) / (2.0 * A);
        t = (t1 < t2 && t1 >= 0) ? t1 : t2;
    }

    return ray_origin + t * ray_direction;
}

// Function to check if a ray intersects with a sphere
bool raysphere(const Vector3d& sphere_center, const Vector3d& ray_origin, const Vector3d& ray_direction, double sphere_radius)
{
    double A = ray_direction.dot(ray_direction);
    double B = 2.0 * ray_direction.dot(ray_origin - sphere_center);
    double C = (ray_origin - sphere_center).dot(ray_origin - sphere_center) - (sphere_radius * sphere_radius);

    double discriminant = (B * B) - (4 * A * C);

    return discriminant >= 0;
}

// Function to get the intersection point of a ray with a parallelogram
Vector3d getPoint(const Vector3d& u_vector, const Vector3d& v_vector, const Vector3d& d_vector, const Vector3d& a_vector, const Vector3d& e_vector)
{
    Matrix3d A;
    A << -u_vector, -v_vector, d_vector;
    Vector3d ae_vector = a_vector - e_vector;
    Vector3d solution_vector = A.colPivHouseholderQr().solve(ae_vector);
    return a_vector + solution_vector(0) * u_vector + solution_vector(1) * v_vector;
}

// Function to check if a ray intersects with a parallelogram
bool raytri(const Vector3d& u_vector, const Vector3d& v_vector, const Vector3d& d_vector, const Vector3d& a_vector, const Vector3d& e_vector)
{
    Matrix3d A;
    A << -u_vector, -v_vector, d_vector;
    Vector3d ae_vector = a_vector - e_vector;
    Vector3d solution_vector = A.colPivHouseholderQr().solve(ae_vector);

    return solution_vector(2) >= 0 && solution_vector(0) >= 0 && solution_vector(0) <= 1 && solution_vector(1) >= 0 && solution_vector(1) <= 1;
}

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Check if the ray intersects with the parallelogram
            if (raytri(pgram_u, pgram_v, ray_direction, pgram_origin, ray_origin))
            {
                // The ray hit the parallelogram, compute the exact intersection point
                Vector3d ray_intersection = getPoint(pgram_u, pgram_v, ray_direction, pgram_origin, ray_origin);

                // Compute normal at the intersection point
                Vector3d ray_normal = pgram_v.cross(pgram_u).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Check if the ray intersects with the parallelogram
            if (raytri(pgram_u, pgram_v, ray_direction, pgram_origin, ray_origin))
            {
                // The ray hit the parallelogram, compute the exact intersection point
                Vector3d ray_intersection = getPoint(pgram_u, pgram_v, ray_direction, pgram_origin, ray_origin);

                // Compute normal at the intersection point
                Vector3d ray_normal = pgram_v.cross(pgram_u).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    MatrixXd R = MatrixXd::Zero(800, 800); // Store the red channel
    MatrixXd G = MatrixXd::Zero(800, 800); // Store the green channel
    MatrixXd B = MatrixXd::Zero(800, 800); // Store the blue channel
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    // Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    // Material parameters
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < R.cols(); ++i)
    {
        for (unsigned j = 0; j < R.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray (origin point and direction)
            const Vector3d ray_origin = camera_origin;
            const Vector3d ray_direction = pixel_center - camera_origin;

            // Intersect with the sphere
            if (raysphere(sphere_center, ray_origin, ray_direction, sphere_radius))
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection = getPointSphere(sphere_center, ray_origin, ray_direction, sphere_radius);

                // Compute normal at the intersection point
                Vector3d ray_normal = ((sphere_center - ray_intersection) / sphere_radius).normalized();

                // Add shading parameters
                Vector3d l = (ray_intersection - light_position).normalized();
                Vector3d v = (ray_intersection - ray_origin).normalized();
                Vector3d h = ((v + l) / (v + l).norm()).normalized();

                Vector3d diffuse = diffuse_color * fmax(0, ray_normal.dot(l));
                Vector3d specular = specular_color * pow(fmax(0, ray_normal.dot(h)), specular_exponent);

                // Combine shading components
                R(i, j) = ambient + diffuse(0) + specular(0);
                G(i, j) = ambient + diffuse(1) + specular(1);
                B(i, j) = ambient + diffuse(2) + specular(2);

                // Clamp to zero
                R(i, j) = std::max(R(i, j), 0.);
                G(i, j) = std::max(G(i, j), 0.);
                B(i, j) = std::max(B(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(R, G, B, A, filename);
}

int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}

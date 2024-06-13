////////////////////////////////////////////////////////////////////////////////
#include <algorithm>
#include <complex>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include <Eigen/Dense>
// Shortcut to avoid  everywhere, DO NOT USE IN .h
using namespace Eigen;
////////////////////////////////////////////////////////////////////////////////

// Problem message "identifier "DATA_DIR" is undefined"
#ifndef DATA_DIR
#define DATA_DIR ""
#endif

const std::string root_path = DATA_DIR;

// Computes the determinant of the matrix whose columns are the vector u and v
double inline det(const Vector2d &u, const Vector2d &v)
{
    Matrix2d mat;
    mat << u, v;
    double determinant = mat.determinant();
    return determinant;
}

// Return true iff [a,b] intersects [c,d]
bool intersect_segment(const Vector2d &a, const Vector2d &b, const Vector2d &c, const Vector2d &d)
{
    auto orientation = [](const Vector2d &p, const Vector2d &q, const Vector2d &r) {
        double val = det(q - p, r - p);
        if (val == 0) return 0;  // 0: Collinear
        return (val > 0) ? 1 : 2; // 1: Clockwise 2: Counterclock wise
    };

    int o1 = orientation(a, b, c);
    int o2 = orientation(a, b, d);
    int o3 = orientation(c, d, a);
    int o4 = orientation(c, d, b);

    // General case
    if (o1 != o2 && o3 != o4)
        return true;

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool is_inside(const std::vector<Vector2d> &poly, const Vector2d &query)
{
    // 1. Compute bounding box and set coordinate of a point outside the polygon
    Vector2d outside(0, 0);
    double min_x = poly[0].x();
    double min_y = poly[0].y();
    double max_x = poly[0].x();
    double max_y = poly[0].y();

    for (const auto& point : poly) {
        if (point.x() < min_x) min_x = point.x();
        if (point.y() < min_y) min_y = point.y();
        if (point.x() > max_x) max_x = point.x();
        if (point.y() > max_y) max_y = point.y();
    }
    outside = Vector2d(max_x + 1, max_y + 1);

    // 2. Cast a ray from the query point to the 'outside' point, count number of intersections
    int count = 0;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Vector2d &a = poly[i];
        const Vector2d &b = poly[(i + 1) % poly.size()];
        if (intersect_segment(query, outside, a, b))
            count++;
    }

    return count % 2 == 1;
    //return true;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<Vector2d> load_xyz(const std::string &filename)
{
    std::vector<Vector2d> points;
    std::ifstream in(filename);
    // When it's undable to open file
    if (!in) {
        std::cerr << "Unable to open file " << filename << std::endl;
        return points;
    }
    // Read the number of points and read each point's coordinates
    int n;
    in >> n;
    points.reserve(n);
    for (int i = 0; i < n; ++i) {
        double x, y, z;
        in >> x >> y >> z;
        points.push_back(Vector2d(x, y));
    }
    return points;
}

void save_xyz(const std::string &filename, const std::vector<Vector2d> &points)
{
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Unable to open file " << filename << std::endl;
        return;
    }
    // Saves points to the file
    // Writes the number of points and each point's coordinates
    out << points.size() << "\n";
    for (const auto &point : points) {
        out << point.x() << " " << point.y() << " 0\n";
    }
}

std::vector<Vector2d> load_obj(const std::string &filename)
{
    std::ifstream in(filename);
    std::vector<Vector2d> points;
    std::vector<Vector2d> poly;
    char key;
    while (in >> key)
    {
        if (key == 'v')
        {
            double x, y, z;
            in >> x >> y >> z;
            points.push_back(Vector2d(x, y));
        }
        else if (key == 'f')
        {
            std::string line;
            std::getline(in, line);
            std::istringstream ss(line);
            int id;
            while (ss >> id)
            {
                poly.push_back(points[id - 1]);
            }
        }
    }
    return poly;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{
    const std::string points_path = root_path + "/points.xyz";
    const std::string poly_path = root_path + "/polygon.obj";

    std::vector<Vector2d> points = load_xyz(points_path);

    ////////////////////////////////////////////////////////////////////////////////
    //Point in polygon
    std::vector<Vector2d> poly = load_obj(poly_path);
    std::vector<Vector2d> result;
    for (size_t i = 0; i < points.size(); ++i)
    {
        if (is_inside(poly, points[i]))
        {
            result.push_back(points[i]);
        }
    }
    save_xyz("output.xyz", result);

    return 0;
}

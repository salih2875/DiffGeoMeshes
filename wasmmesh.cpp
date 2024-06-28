#include <vector>
#include <iterator>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>
#include "colormap.h"


class Vertix
{
public:
    double x;
    double y;
    double z;
    Vertix() : x(0), y(0), z(0) {}
    Vertix(double x, double y, double z) : x(x), y(y), z(z) {}
};

class Triangle
{
public:
    int x;
    int y;
    int z;
    Triangle() : x(0), y(0), z(0) {}
    Triangle(int x, int y, int z) : x(x), y(y), z(z) {}
};

void read_file(std::ifstream &file, std::vector<Vertix> &vertices, std::vector<Triangle> &triangles)
{
    std::string current_line;
    std::getline(file, current_line);
    if (current_line != "OFF")
    {
        std::runtime_error("Why no OFF header?");
    }

    std::getline(file, current_line);
    std::istringstream iss(current_line);
    int number_of_vertices, number_of_faces, edges;

    iss >> number_of_vertices >> number_of_faces >> edges;

    for (int ith_vert = 0; ith_vert < number_of_vertices; ++ith_vert)
    {
        std::getline(file, current_line);
        iss.clear();
        iss.str(current_line);
        // I assume vertix has 3 coordinates
        double x1, y1, z1;
        iss >> x1 >> y1 >> z1;
        Vertix v;
        v.x = x1;
        v.y = y1;
        v.z = z1;
        vertices.push_back(v);
    }

    for (int ith_triangle = 0; ith_triangle < number_of_faces; ++ith_triangle)
    {
        std::getline(file, current_line);
        iss.clear();
        iss.str(current_line);
        int x1, y1, z1;
        int dummy;
        iss >> dummy >> x1 >> y1 >> z1;
        Triangle t;
        t.x = x1;
        t.y = y1;
        t.z = z1;
        triangles.push_back(t);
    }
}

// find 1 ring neighborhood of x
std::vector<Triangle> one_ring_neighborhood(int index, std::vector<Triangle> triangles)
{
    std::vector<Triangle> neighborhood;
    for (auto cur_triangle : triangles)
    {
        if (index == cur_triangle.x || index == cur_triangle.y || index == cur_triangle.z)
        {
            neighborhood.push_back(cur_triangle);
        }
    }

    return neighborhood;
}

double Square(float x)
{
    double ans = x * x;
    return ans;
}

double norm_of_vector(Vertix v)
{
    double norm = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return norm;
}

double dot_product(Vertix a, Vertix b)
{
    double dotp = a.x * b.x + a.y * b.y + a.z * b.z;
    return dotp;
}

double angleof(Vertix a, Vertix b)
{
    // sometimes angle>1 so acos is nan
    // maybe i should check it
    double norm_of_a = norm_of_vector(a);
    double norm_of_b = norm_of_vector(b);
    double dotp = dot_product(a, b);

    double angle = acos(dotp / (norm_of_a * norm_of_b));
    return angle;
}

Vertix vector_from_points(Vertix b, Vertix e)
{
    Vertix c;
    c.x = b.x - e.x;
    c.y = b.y - e.y;
    c.z = b.z - e.z;
    return c;
}

int check_triangle_obtuse(std::vector<Vertix> vertices, Triangle triangle, int index)
{
    int type_of_obtuse = 0;

    std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};
    dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());
    int x1 = index;
    int y1 = dummyvec[0];
    int z1 = dummyvec[1];

    Vertix a = vector_from_points(vertices[x1], vertices[y1]);
    Vertix b = vector_from_points(vertices[x1], vertices[z1]);
    double alpha = angleof(a, b);

    double pi = M_PI / 2;

    if (alpha > pi)
    {
        type_of_obtuse = 1;
    }
    Vertix a1 = vector_from_points(vertices[y1], vertices[x1]);
    Vertix b1 = vector_from_points(vertices[y1], vertices[z1]);
    double beta = angleof(a1, b1);

    if (beta > pi)
    {
        type_of_obtuse = 2;
    }

    Vertix a2 = vector_from_points(vertices[z1], vertices[x1]);
    Vertix b2 = vector_from_points(vertices[z1], vertices[y1]);
    double theta = angleof(a2, b2);

    if (theta > pi)
    {
        type_of_obtuse = 2;
    }

    return type_of_obtuse;
}

double voronoi_area_of_non_obtuse(std::vector<Vertix> vertices, Triangle triangle, int index)
{
    double a_mixed = 0;
    // int x1{triangle.x}; int y1{triangle.y}; int z1{triangle.z};

    std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};

    dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());

    int x1 = index;

    int y1 = dummyvec[0];

    int z1 = dummyvec[1];

    // our index is point P in triangle PRQ
    // Vectors to P
    Vertix a = vector_from_points(vertices[x1], vertices[y1]);
    Vertix b = vector_from_points(vertices[x1], vertices[z1]);

    // angle Q
    Vertix v1 = vector_from_points(vertices[y1], vertices[x1]);
    Vertix v2 = vector_from_points(vertices[y1], vertices[z1]);

    // angle R
    Vertix v3 = vector_from_points(vertices[z1], vertices[x1]);
    Vertix v4 = vector_from_points(vertices[z1], vertices[y1]);

    double norm1 = norm_of_vector(a);
    double norm2 = norm_of_vector(b);
    // double angleq = angleof(v1,v2);
    // double angler = angleof(v3,v4);
    double cot_q = 1 / tan(angleof(v1, v2));
    double cot_r = 1 / tan(angleof(v3, v4));

    a_mixed = (cot_q * norm1 * norm1 + cot_r * norm2 * norm2) / 8;
    return a_mixed;
}

double voronoi_area_of_obtuse(std::vector<Vertix> vertices, Triangle triangle, int index)
{
    double a_mixed = 0;

    // testing deletion
    std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};
    dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());
    int x1 = index;
    int y1 = dummyvec[0];
    int z1 = dummyvec[1];

    Vertix const a1 = vector_from_points(vertices[x1], vertices[y1]);
    Vertix const b1 = vector_from_points(vertices[x1], vertices[z1]);
    Vertix const c1 = vector_from_points(vertices[y1], vertices[z1]);

    const double a = norm_of_vector(a1);
    const double b = norm_of_vector(b1);
    const double c = norm_of_vector(c1);
    const double p = (a + b + c) / 2;
    a_mixed = sqrt(p * (p - a) * (p - b) * (p - c));

    return a_mixed;
}

double A_mixed(std::vector<Vertix> vertices, std::vector<Triangle> triangles, int index)
{
    double amixed = 0;
    std::vector<Triangle> neighborhood = one_ring_neighborhood(index, triangles);
    double sum = 0;
    for (auto triangle : neighborhood)
    {
        int flag = check_triangle_obtuse(vertices, triangle, index);
        if (flag == 0)
        {
            sum += voronoi_area_of_non_obtuse(vertices, triangle, index);
        }
        else if (flag == 1)
        {
            double area = voronoi_area_of_obtuse(vertices, triangle, index);
            sum += area / 2.0f;
        }
        else if (flag == 2)
        {
            double area = voronoi_area_of_obtuse(vertices, triangle, index);
            sum += area / 4.0f;
        }
    }
    amixed += sum;
    return amixed;
}

double mean_curvature(std::vector<Vertix> vertices, std::vector<Triangle> triangles, int index, double A_mixed_value)
{
    // our index is P opposite angles are q and r
    Vertix meancurv{0, 0, 0};

    for (auto triangle : triangles)
    {
        std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};
        dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());
        int x1 = index;
        int y1 = dummyvec[0];
        int z1 = dummyvec[1];

        // our index is point P in triangle PRQ
        // Vectors to P
        Vertix a = vector_from_points(vertices[x1], vertices[y1]);
        Vertix b = vector_from_points(vertices[x1], vertices[z1]);

        // angle Q
        Vertix v1 = vector_from_points(vertices[y1], vertices[x1]);
        Vertix v2 = vector_from_points(vertices[y1], vertices[z1]);

        // angle R
        Vertix v3 = vector_from_points(vertices[z1], vertices[x1]);
        Vertix v4 = vector_from_points(vertices[z1], vertices[y1]);

        double const cot_q = 1 / tan(angleof(v1, v2));
        double cot_r = 1 / tan(angleof(v3, v4));
        meancurv.x += cot_q * b.x + cot_r * a.x;
        meancurv.y += cot_q * b.y + cot_r * a.y;
        meancurv.z += cot_q * b.z + cot_r * a.z;
    }

    double denom = 1 / (2 * A_mixed_value);

    Vertix mean;
    mean.x = meancurv.x * denom;
    mean.y = meancurv.y * denom;
    mean.z = meancurv.z * denom;
    double Kh = norm_of_vector(mean) / 2;
    return Kh;
}

double gauss_curvature(std::vector<Vertix> vertices, std::vector<Triangle> triangles, int index, double A_mixed_value)
{
    double sum = 0;
    for (auto triangle : triangles)
    {
        std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};
        dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());
        int x1 = index;
        int y1 = dummyvec[0];
        int z1 = dummyvec[1];

        Vertix const a = vector_from_points(vertices[x1], vertices[y1]);
        Vertix const b = vector_from_points(vertices[x1], vertices[z1]);

        double theta = angleof(a, b);

        sum += theta;
    }

    double gauss = (2 * M_PI - sum) / A_mixed_value;
    return gauss;
}

double first_principal(double meancurv, double gausscurv)
{
    double delta = 0;
    if (meancurv * meancurv > gausscurv)
    {
        delta = meancurv * meancurv - gausscurv;
    }

    double k = meancurv + sqrt(delta);
    return k;
}

double second_principal(double meancurv, double gausscurv)
{
    double delta = 0;
    if (meancurv * meancurv > gausscurv)
    {
        delta = meancurv * meancurv - gausscurv;
    }

    double k = meancurv - sqrt(delta);
    return k;
}

void normalize(double *scalars, int sv)
{
    for (int i = 0; i < sv; ++i)
    {
        if (std::isinf(scalars[i]))
        {
            scalars[i] = 0;
        }
        if (std::isnan(scalars[i]))
        {
            scalars[i] = 0;
        }
    }
    double min = *std::min_element(scalars, scalars + sv);
    double max = *std::max_element(scalars, scalars + sv);
    double delta = max - min;
    for (int i = 0; i < sv; ++i)
    {
        scalars[i] = (scalars[i] - min) / delta;
    }
}

void interpolate(double *scalars, int sv, std::vector<double> &colors)
{
    for (int i = 0; i < sv; ++i)
    {
        double x = scalars[i];
        if (x < 0)
        {
            x = 0;
        }
        if (x > 1)
        {
            x = 1;
        }
        const int lo = std::floor(x * (N_VIRIDIS - 1));
        const int hi = std::ceil(x * (N_VIRIDIS - 1));
        const double r = (viridis[lo].r + viridis[hi].r) / 2;
        const double g = (viridis[lo].g + viridis[hi].g) / 2;
        const double b = (viridis[lo].b + viridis[hi].b) / 2;
        colors.push_back(r);
        colors.push_back(g);
        colors.push_back(b);
    }
}

// WASM export function
extern "C"
{
    void printgauss(int sv, int st, double *before_vertices, int *before_triangles, double *gauss, double *mean, double *first, double *second, double *gcolors, double *mcolors, double *fcolors, double *scolors)
    {
        std::vector<Vertix> vertices;
        for (size_t i = 0; i < 3 * sv; i += 3)
        {
            Vertix v;
            v.x = before_vertices[i];
            v.y = before_vertices[i + 1];
            v.z = before_vertices[i + 2];
            vertices.push_back(v);
        }
        std::vector<Triangle> triangles;
        for (size_t i = 0; i < 3 * st; i += 3)
        {
            Triangle v;
            v.x = before_triangles[i];
            v.y = before_triangles[i + 1];
            v.z = before_triangles[i + 2];
            triangles.push_back(v);
        }

        for (size_t i = 0; i < vertices.size(); ++i)
        {
            double amixed = A_mixed(vertices, triangles, i);
            std::vector<Triangle> neighborhood = one_ring_neighborhood(i, triangles);
            gauss[i] = gauss_curvature(vertices, neighborhood, i, amixed);
            mean[i] = mean_curvature(vertices, neighborhood, i, amixed);
            first[i] = first_principal(mean[i], gauss[i]);
            second[i] = second_principal(mean[i], gauss[i]);
            // printf("Try number %d Gauss is: %f\n", i, m);
        }

        normalize(gauss, sv);
        normalize(mean, sv);
        normalize(first, sv);
        normalize(second, sv);

        std::vector<double> vgcolors;
        std::vector<double> vmcolors;
        std::vector<double> vfcolors;
        std::vector<double> vscolors;
        interpolate(gauss, sv, vgcolors);
        interpolate(mean, sv, vmcolors);
        interpolate(first, sv, vfcolors);
        interpolate(second, sv, vscolors);

        for (int i = 0; i < vgcolors.size(); ++i)
        {
            gcolors[i] = vgcolors[i];
        }
        for (int i = 0; i < vmcolors.size(); ++i)
        {
            mcolors[i] = vmcolors[i];
        }
        for (int i = 0; i < vfcolors.size(); ++i)
        {
            fcolors[i] = vfcolors[i];
        }
        for (int i = 0; i < vscolors.size(); ++i)
        {
            scolors[i] = vscolors[i];
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc > 2)
    {
        throw std::runtime_error("Only one file path is required");
    }

    std::string Filename;
    Filename = argv[1];

    std::vector<Vertix> vertices;
    std::vector<Triangle> triangles;
    std::ifstream file(Filename);
    if (!file)
    {
        throw std::runtime_error("Provide valid file path");
    }
    read_file(file, vertices, triangles);

    // TODO: make this all in one file and just read this one file in python
    std::ofstream gaussfile;
    gaussfile.open("gauss.txt");

    std::ofstream meanfile;
    meanfile.open("mean.txt");

    std::ofstream frfile;
    frfile.open("first.txt");

    std::ofstream secondfile;
    secondfile.open("second.txt");

    for (size_t index = 0; index < vertices.size(); ++index)
    {

        std::vector<Triangle> neighbourhood = one_ring_neighborhood(index, triangles);

        double a_mixed = A_mixed(vertices, neighbourhood, index);
        if (std::isnan(a_mixed) || std::isinf(a_mixed))
        {
            a_mixed = 0;
        }

        double gauss = gauss_curvature(vertices, neighbourhood, index, a_mixed);
        if (std::isnan(gauss) || std::isinf(gauss))
        {
            gauss = 0;
        }
        if (index != vertices.size() - 1)
            gaussfile << gauss << ",";
        if (index == vertices.size() - 1)
            gaussfile << gauss;

        double meancurv = mean_curvature(vertices, neighbourhood, index, a_mixed);
        if (std::isnan(meancurv) || std::isinf(meancurv))
        {
            meancurv = 0;
        }
        if (index != vertices.size() - 1)
            meanfile << meancurv << ",";
        if (index == vertices.size() - 1)
            meanfile << meancurv;

        double k1 = first_principal(meancurv, gauss);
        if (std::isnan(k1) || std::isinf(k1))
        {
            k1 = 0;
        }
        if (index != vertices.size() - 1)
            frfile << k1 << ",";
        if (index == vertices.size() - 1)
            frfile << k1;

        double k2 = second_principal(meancurv, gauss);
        if (std::isnan(k2) || std::isinf(k2))
        {
            k2 = 0;
        }
        if (index != vertices.size() - 1)
            secondfile << k2 << ",";
        if (index == vertices.size() - 1)
            secondfile << k2;
    }

    gaussfile.close();
    meanfile.close();
    frfile.close();
    secondfile.close();
}
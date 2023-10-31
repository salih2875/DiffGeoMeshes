#include <vector>
#include <iterator>
#include <fstream>
#include <string>
#include <string.h>
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm>

// TODO: create scalar_product function

class Vertix
{
public:
    double x;
    double y;
    double z;
    Vertix() : x(0) , y(0), z(0) {}
    Vertix(double x, double y, double z) : x(x), y(y), z(z) { }
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

class Color
{
public:
    double r;
    double g;
    double b;
    Color() : r(0), g(0), b(0) {}
    Color(double r, double g, double b) : r(r), g(g), b(b) { }
};

// declare true_vertices vector in main and in read_file make it equal to vertices
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
    // std::vector<Triangle> triangles;
    // for (size_t i=0; i<before_triangles.size(); i+=3) {
    //     Triangle v;
    //     v.x = before_triangles[i]; 
    //     v.y = before_triangles[i+1];
    //     v.z = before_triangles[i+2];
    //     triangles.push_back(v);
    // }

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
    // sometimes angle>1 so acos is nan i dont know why
    double norm_of_a = norm_of_vector(a);
    double norm_of_b = norm_of_vector(b);
    double dotp = dot_product(a, b);

    // double checkangle = dotp / (norm_of_a*norm_of_b);
    // std::cout << checkangle << std::endl;
    // if (checkangle>1) {
    //     std::cout << "angle between vectors is lower than -1" << "\n";
    // }

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

Triangle change_indexes(Triangle triangle, int index)
{
    std::vector<int> dummyvec{triangle.x, triangle.y, triangle.z};

    dummyvec.erase(std::remove(dummyvec.begin(), dummyvec.end(), index), dummyvec.end());

    Triangle t;

    t.x = index;

    t.y = dummyvec[0];

    t.z = dummyvec[1];

    return t;
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
    // std::vector<Vertix> vertices;
    // for (size_t i=0; i<before_vertices.size(); i+=3) {
    //     Vertix v;
    //     v.x = before_vertices[i]; 
    //     v.y = before_vertices[i+1];
    //     v.z = before_vertices[i+2];
    //     vertices.push_back(v);
    // }
    
    // std::vector<Triangle> triangles;
    // for (size_t i=0; i<before_triangles.size(); i+=3) {
    //     Triangle v;
    //     v.x = before_triangles[i]; 
    //     v.y = before_triangles[i+1];
    //     v.z = before_triangles[i+2];
    //     triangles.push_back(v);
    // }

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
    // std::vector<Vertix> vertices;
    // for (size_t i=0; i<before_vertices.size(); i+=3) {
    //     Vertix v;
    //     v.x = before_vertices[i]; 
    //     v.y = before_vertices[i+1];
    //     v.z = before_vertices[i+2];
    //     vertices.push_back(v);
    // }
    
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


void normalize(double *scalars,int sv) {
    for (int i=0; i<sv; ++i) {
        if (std::isinf(scalars[i])) {
            scalars[i]=0;
        }
        if (std::isnan(scalars[i])) {
            scalars[i]=0;
        }
    }
    double min = *std::min_element(scalars,scalars+sv);
    double max = *std::max_element(scalars,scalars+sv);
    double delta = max - min;
    for (int i=0; i<sv; ++i) {
        scalars[i] = (scalars[i]-min)/delta;
    }
}

const int N_VIRIDIS = 255;

Color viridis[N_VIRIDIS+6] = {
    { 0.267004, 0.004874, 0.329415 },
    { 0.268510, 0.009605, 0.335427 },
    { 0.269944, 0.014625, 0.341379 },
    { 0.271305, 0.019942, 0.347269 },
    { 0.272594, 0.025563, 0.353093 },
    { 0.273809, 0.031497, 0.358853 },
    { 0.274952, 0.037752, 0.364543 },
    { 0.276022, 0.044167, 0.370164 },
    { 0.277018, 0.050344, 0.375715 },
    { 0.277941, 0.056324, 0.381191 },
    { 0.278791, 0.062145, 0.386592 },
    { 0.279566, 0.067836, 0.391917 },
    { 0.280267, 0.073417, 0.397163 },
    { 0.280894, 0.078907, 0.402329 },
    { 0.281446, 0.084320, 0.407414 },
    { 0.281924, 0.089666, 0.412415 },
    { 0.282327, 0.094955, 0.417331 },
    { 0.282656, 0.100196, 0.422160 },
    { 0.282910, 0.105393, 0.426902 },
    { 0.283091, 0.110553, 0.431554 },
    { 0.283197, 0.115680, 0.436115 },
    { 0.283229, 0.120777, 0.440584 },
    { 0.283187, 0.125848, 0.444960 },
    { 0.283072, 0.130895, 0.449241 },
    { 0.282884, 0.135920, 0.453427 },
    { 0.282623, 0.140926, 0.457517 },
    { 0.282290, 0.145912, 0.461510 },
    { 0.281887, 0.150881, 0.465405 },
    { 0.281412, 0.155834, 0.469201 },
    { 0.280868, 0.160771, 0.472899 },
    { 0.280255, 0.165693, 0.476498 },
    { 0.279574, 0.170599, 0.479997 },
    { 0.278826, 0.175490, 0.483397 },
    { 0.278012, 0.180367, 0.486697 },
    { 0.277134, 0.185228, 0.489898 },
    { 0.276194, 0.190074, 0.493001 },
    { 0.275191, 0.194905, 0.496005 },
    { 0.274128, 0.199721, 0.498911 },
    { 0.273006, 0.204520, 0.501721 },
    { 0.271828, 0.209303, 0.504434 },
    { 0.270595, 0.214069, 0.507052 },
    { 0.269308, 0.218818, 0.509577 },
    { 0.267968, 0.223549, 0.512008 },
    { 0.266580, 0.228262, 0.514349 },
    { 0.265145, 0.232956, 0.516599 },
    { 0.263663, 0.237631, 0.518762 },
    { 0.262138, 0.242286, 0.520837 },
    { 0.260571, 0.246922, 0.522828 },
    { 0.258965, 0.251537, 0.524736 },
    { 0.257322, 0.256130, 0.526563 },
    { 0.255645, 0.260703, 0.528312 },
    { 0.253935, 0.265254, 0.529983 },
    { 0.252194, 0.269783, 0.531579 },
    { 0.250425, 0.274290, 0.533103 },
    { 0.248629, 0.278775, 0.534556 },
    { 0.246811, 0.283237, 0.535941 },
    { 0.244972, 0.287675, 0.537260 },
    { 0.243113, 0.292092, 0.538516 },
    { 0.241237, 0.296485, 0.539709 },
    { 0.239346, 0.300855, 0.540844 },
    { 0.237441, 0.305202, 0.541921 },
    { 0.235526, 0.309527, 0.542944 },
    { 0.233603, 0.313828, 0.543914 },
    { 0.231674, 0.318106, 0.544834 },
    { 0.229739, 0.322361, 0.545706 },
    { 0.227802, 0.326594, 0.546532 },
    { 0.225863, 0.330805, 0.547314 },
    { 0.223925, 0.334994, 0.548053 },
    { 0.221989, 0.339161, 0.548752 },
    { 0.220057, 0.343307, 0.549413 },
    { 0.218130, 0.347432, 0.550038 },
    { 0.216210, 0.351535, 0.550627 },
    { 0.214298, 0.355619, 0.551184 },
    { 0.212395, 0.359683, 0.551710 },
    { 0.210503, 0.363727, 0.552206 },
    { 0.208623, 0.367752, 0.552675 },
    { 0.206756, 0.371758, 0.553117 },
    { 0.204903, 0.375746, 0.553533 },
    { 0.203063, 0.379716, 0.553925 },
    { 0.201239, 0.383670, 0.554294 },
    { 0.199430, 0.387607, 0.554642 },
    { 0.197636, 0.391528, 0.554969 },
    { 0.195860, 0.395433, 0.555276 },
    { 0.194100, 0.399323, 0.555565 },
    { 0.192357, 0.403199, 0.555836 },
    { 0.190631, 0.407061, 0.556089 },
    { 0.188923, 0.410910, 0.556326 },
    { 0.187231, 0.414746, 0.556547 },
    { 0.185556, 0.418570, 0.556753 },
    { 0.183898, 0.422383, 0.556944 },
    { 0.182256, 0.426184, 0.557120 },
    { 0.180629, 0.429975, 0.557282 },
    { 0.179019, 0.433756, 0.557430 },
    { 0.177423, 0.437527, 0.557565 },
    { 0.175841, 0.441290, 0.557685 },
    { 0.174274, 0.445044, 0.557792 },
    { 0.172719, 0.448791, 0.557885 },
    { 0.171176, 0.452530, 0.557965 },
    { 0.169646, 0.456262, 0.558030 },
    { 0.168126, 0.459988, 0.558082 },
    { 0.166617, 0.463708, 0.558119 },
    { 0.165117, 0.467423, 0.558141 },
    { 0.163625, 0.471133, 0.558148 },
    { 0.162142, 0.474838, 0.558140 },
    { 0.160665, 0.478540, 0.558115 },
    { 0.159194, 0.482237, 0.558073 },
    { 0.157729, 0.485932, 0.558013 },
    { 0.156270, 0.489624, 0.557936 },
    { 0.154815, 0.493313, 0.557840 },
    { 0.153364, 0.497000, 0.557724 },
    { 0.151918, 0.500685, 0.557587 },
    { 0.150476, 0.504369, 0.557430 },
    { 0.149039, 0.508051, 0.557250 },
    { 0.147607, 0.511733, 0.557049 },
    { 0.146180, 0.515413, 0.556823 },
    { 0.144759, 0.519093, 0.556572 },
    { 0.143343, 0.522773, 0.556295 },
    { 0.141935, 0.526453, 0.555991 },
    { 0.140536, 0.530132, 0.555659 },
    { 0.139147, 0.533812, 0.555298 },
    { 0.137770, 0.537492, 0.554906 },
    { 0.136408, 0.541173, 0.554483 },
    { 0.135066, 0.544853, 0.554029 },
    { 0.133743, 0.548535, 0.553541 },
    { 0.132444, 0.552216, 0.553018 },
    { 0.131172, 0.555899, 0.552459 },
    { 0.129933, 0.559582, 0.551864 },
    { 0.128729, 0.563265, 0.551229 },
    { 0.127568, 0.566949, 0.550556 },
    { 0.126453, 0.570633, 0.549841 },
    { 0.125394, 0.574318, 0.549086 },
    { 0.124395, 0.578002, 0.548287 },
    { 0.123463, 0.581687, 0.547445 },
    { 0.122606, 0.585371, 0.546557 },
    { 0.121831, 0.589055, 0.545623 },
    { 0.121148, 0.592739, 0.544641 },
    { 0.120565, 0.596422, 0.543611 },
    { 0.120092, 0.600104, 0.542530 },
    { 0.119738, 0.603785, 0.541400 },
    { 0.119512, 0.607464, 0.540218 },
    { 0.119423, 0.611141, 0.538982 },
    { 0.119483, 0.614817, 0.537692 },
    { 0.119699, 0.618490, 0.536347 },
    { 0.120081, 0.622161, 0.534946 },
    { 0.120638, 0.625828, 0.533488 },
    { 0.121380, 0.629492, 0.531973 },
    { 0.122312, 0.633153, 0.530398 },
    { 0.123444, 0.636809, 0.528763 },
    { 0.124780, 0.640461, 0.527068 },
    { 0.126326, 0.644107, 0.525311 },
    { 0.128087, 0.647749, 0.523491 },
    { 0.130067, 0.651384, 0.521608 },
    { 0.132268, 0.655014, 0.519661 },
    { 0.134692, 0.658636, 0.517649 },
    { 0.137339, 0.662252, 0.515571 },
    { 0.140210, 0.665859, 0.513427 },
    { 0.143303, 0.669459, 0.511215 },
    { 0.146616, 0.673050, 0.508936 },
    { 0.150148, 0.676631, 0.506589 },
    { 0.153894, 0.680203, 0.504172 },
    { 0.157851, 0.683765, 0.501686 },
    { 0.162016, 0.687316, 0.499129 },
    { 0.166383, 0.690856, 0.496502 },
    { 0.170948, 0.694384, 0.493803 },
    { 0.175707, 0.697900, 0.491033 },
    { 0.180653, 0.701402, 0.488189 },
    { 0.185783, 0.704891, 0.485273 },
    { 0.191090, 0.708366, 0.482284 },
    { 0.196571, 0.711827, 0.479221 },
    { 0.202219, 0.715272, 0.476084 },
    { 0.208030, 0.718701, 0.472873 },
    { 0.214000, 0.722114, 0.469588 },
    { 0.220124, 0.725509, 0.466226 },
    { 0.226397, 0.728888, 0.462789 },
    { 0.232815, 0.732247, 0.459277 },
    { 0.239374, 0.735588, 0.455688 },
    { 0.246070, 0.738910, 0.452024 },
    { 0.252899, 0.742211, 0.448284 },
    { 0.259857, 0.745492, 0.444467 },
    { 0.266941, 0.748751, 0.440573 },
    { 0.274149, 0.751988, 0.436601 },
    { 0.281477, 0.755203, 0.432552 },
    { 0.288921, 0.758394, 0.428426 },
    { 0.296479, 0.761561, 0.424223 },
    { 0.304148, 0.764704, 0.419943 },
    { 0.311925, 0.767822, 0.415586 },
    { 0.319809, 0.770914, 0.411152 },
    { 0.327796, 0.773980, 0.406640 },
    { 0.335885, 0.777018, 0.402049 },
    { 0.344074, 0.780029, 0.397381 },
    { 0.352360, 0.783011, 0.392636 },
    { 0.360741, 0.785964, 0.387814 },
    { 0.369214, 0.788888, 0.382914 },
    { 0.377779, 0.791781, 0.377939 },
    { 0.386433, 0.794644, 0.372886 },
    { 0.395174, 0.797475, 0.367757 },
    { 0.404001, 0.800275, 0.362552 },
    { 0.412913, 0.803041, 0.357269 },
    { 0.421908, 0.805774, 0.351910 },
    { 0.430983, 0.808473, 0.346476 },
    { 0.440137, 0.811138, 0.340967 },
    { 0.449368, 0.813768, 0.335384 },
    { 0.458674, 0.816363, 0.329727 },
    { 0.468053, 0.818921, 0.323998 },
    { 0.477504, 0.821444, 0.318195 },
    { 0.487026, 0.823929, 0.312321 },
    { 0.496615, 0.826376, 0.306377 },
    { 0.506271, 0.828786, 0.300362 },
    { 0.515992, 0.831158, 0.294279 },
    { 0.525776, 0.833491, 0.288127 },
    { 0.535621, 0.835785, 0.281908 },
    { 0.545524, 0.838039, 0.275626 },
    { 0.555484, 0.840254, 0.269281 },
    { 0.565498, 0.842430, 0.262877 },
    { 0.575563, 0.844566, 0.256415 },
    { 0.585678, 0.846661, 0.249897 },
    { 0.595839, 0.848717, 0.243329 },
    { 0.606045, 0.850733, 0.236712 },
    { 0.616293, 0.852709, 0.230052 },
    { 0.626579, 0.854645, 0.223353 },
    { 0.636902, 0.856542, 0.216620 },
    { 0.647257, 0.858400, 0.209861 },
    { 0.657642, 0.860219, 0.203082 },
    { 0.668054, 0.861999, 0.196293 },
    { 0.678489, 0.863742, 0.189503 },
    { 0.688944, 0.865448, 0.182725 },
    { 0.699415, 0.867117, 0.175971 },
    { 0.709898, 0.868751, 0.169257 },
    { 0.720391, 0.870350, 0.162603 },
    { 0.730889, 0.871916, 0.156029 },
    { 0.741388, 0.873449, 0.149561 },
    { 0.751884, 0.874951, 0.143228 },
    { 0.762373, 0.876424, 0.137064 },
    { 0.772852, 0.877868, 0.131109 },
    { 0.783315, 0.879285, 0.125405 },
    { 0.793760, 0.880678, 0.120005 },
    { 0.804182, 0.882046, 0.114965 },
    { 0.814576, 0.883393, 0.110347 },
    { 0.824940, 0.884720, 0.106217 },
    { 0.835270, 0.886029, 0.102646 },
    { 0.845561, 0.887322, 0.099702 },
    { 0.855810, 0.888601, 0.097452 },
    { 0.866013, 0.889868, 0.095953 },
    { 0.876168, 0.891125, 0.095250 },
    { 0.886271, 0.892374, 0.095374 },
    { 0.896320, 0.893616, 0.096335 },
    { 0.906311, 0.894855, 0.098125 },
    { 0.916242, 0.896091, 0.100717 },
    { 0.926106, 0.897330, 0.104071 },
    { 0.935904, 0.898570, 0.108131 },
    { 0.945636, 0.899815, 0.112838 },
    { 0.955300, 0.901065, 0.118128 },
    { 0.964894, 0.902323, 0.123941 },
    { 0.974417, 0.903590, 0.130215 },
    { 0.983868, 0.904867, 0.136897 },
    { 0.993248, 0.906157, 0.143936 }
};


void interpolate(double* scalars, int sv, std::vector<double> &colors) {
    for (int i=0; i<sv; ++i) {
        double x = scalars[i];
        if (x<0) {
            x = 0;
        }
        if (x > 1) {
            x = 1;
        }
        const int lo = std::floor(x*(N_VIRIDIS-1));
        const int hi = std::ceil(x*(N_VIRIDIS-1));
        const double r =  (viridis[lo].r + viridis[hi].r) / 2;
        const double g = (viridis[lo].g + viridis[hi].g) / 2;
        const double b = (viridis[lo].b + viridis[hi].b) / 2;
        colors.push_back(r);
        colors.push_back(g);
        colors.push_back(b);
    }
}

extern "C" {
void printgauss(int sv, int st, double* before_vertices, int* before_triangles, double* gauss, double* mean, double* first, double* second, double* gcolors, double* mcolors, double* fcolors, double* scolors) {
    std::vector<Vertix> vertices;
    for (size_t i=0; i<3*sv; i+=3) {
        Vertix v;
        v.x = before_vertices[i]; 
        v.y = before_vertices[i+1];
        v.z = before_vertices[i+2];
        vertices.push_back(v);
    }
    std::vector<Triangle> triangles;
    for (size_t i=0; i<3*st; i+=3) {
        Triangle v;
        v.x = before_triangles[i]; 
        v.y = before_triangles[i+1];
        v.z = before_triangles[i+2];
        triangles.push_back(v);
    }
    
    for (size_t i=0; i<vertices.size(); ++i) {
        double amixed = A_mixed(vertices,triangles,i);
        std::vector<Triangle> neighborhood = one_ring_neighborhood(i, triangles);
        gauss[i] = gauss_curvature(vertices, neighborhood, i, amixed);
        mean[i] = mean_curvature(vertices,neighborhood, i, amixed);
        first[i] = first_principal(mean[i], gauss[i]);
        second[i] = second_principal(mean[i], gauss[i]);
        // printf("Try number %d Gauss is: %f\n", i, m);
    }

    normalize(gauss,sv);
    normalize(mean,sv);
    normalize(first,sv);
    normalize(second,sv);
    
    std::vector<double> vgcolors;
    std::vector<double> vmcolors;
    std::vector<double> vfcolors;
    std::vector<double> vscolors;
    interpolate(gauss, sv, vgcolors);
    interpolate(mean, sv, vmcolors);
    interpolate(first, sv, vfcolors);
    interpolate(second, sv, vscolors);
   
    for (int i=0; i<vgcolors.size(); ++i) {
        gcolors[i] = vgcolors[i];
    }
    for (int i=0; i<vmcolors.size(); ++i) {
        mcolors[i] = vmcolors[i];
    }
    for (int i=0; i<vfcolors.size(); ++i) {
        fcolors[i] = vfcolors[i];
    }
    for (int i=0; i<vscolors.size(); ++i) {
        scolors[i] = vscolors[i];
    }


    // for (int i=0; i<sv; ++i ) printf("%d %f\n",i, gcolors[i]);
}

}




int main()
{


    int v;
    int f;
    std::cin>>v>>f;
    // std::vector<int> triangles;
    // std::vector<double> vertices;
    const int sv = 3*v;
    const int st = 3*f;
    double vertices[3000];
    int triangles[5000];

    
    for (int i=0; i<3*v; ++i) {
        // double x;
        std::cin >> vertices[i];
        // vertices.push_back(x);
    }

    for (int i=0; i<3*f; ++i) {
        // int x;
        std::cin >> triangles[i] ;
        // triangles.push_back(x);
    }

    double gauss[sv];
    double mean[sv];
    double first[sv];
    double second[sv];


    // for (int i=0; i<v; ++i) {
    //     double x,y,z;
    //     std::cin >> x >> y >> z;
    //     vertices.push_back(Vertix(x,y,z));
    // }

    // for (int i=0; i<f; ++i) {
    //     int x,y,z;
    //     std::cin >> x >> y >> z; 
    //     triangles.push_back(Triangle(x,y,z));
    // }

     // printgauss(v, f, vertices, triangles);
    // for (int i=0; i<v; ++i) {
    //     std::vector<Triangle> neighborhood = one_ring_neighborhood(i, triangles);
    //     double amixed = A_mixed(vertices, triangles, i);
    //     double gauss = gauss_curvature(vertices, neighborhood, i, amixed);
    //     std::cout << amixed << ", ";
    // }
}


/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Robb Tran
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits>
#include <random>

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100
#define MAX_REFLECTION 3

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define EPSILON 1e-5
#define PI 3.14159
#define SOFT_SHADOW true
#define ANTI_ALIASING true
#define SUB_LIGHTS 30
#define REFLECT_RATIO 0.1

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct Ray
{
    double origin[3];
    double direction[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
int num_rays = 0;

double a, rad, x_min, y_min, x_max, y_max, screen_width, cell_width, screen_height, cell_height;

void clamp(double(&color)[3]) {
    if (color[0] > 1.0) color[0] = 1.0;
    if (color[1] > 1.0) color[1] = 1.0;
    if (color[2] > 1.0) color[2] = 1.0;
}

inline double rand_val()
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<double> dist(0.f, 1);

    return dist(rng);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

// --- helper --- //
void normalize(double& x, double& y, double& z) {
    double magnitude = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    if (magnitude >= 0)
    {
        x = x / magnitude;
        y = y / magnitude;
        z = z / magnitude;
    }
    else
    {
        x = 0;
        y = 0;
        z = 0;
    }
}

double dot_product(const double v1[3], const double v2[3]) {
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

void cross_product(const double a[1], const double b[3], double (&c)[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void reflect_dir(const double L[3], const double N[3], double (&R)[3]) {
    // R = 2 * (L DOT N) * N - L
    double factor = 2 * dot_product(N, L);
    R[0] = factor * N[0] - L[0];
    R[1] = factor * N[1] - L[1];
    R[2] = factor * N[2] - L[2];
}

// --- init --- //
void init_screen() {
    // set up
    a = static_cast<double> (WIDTH) / static_cast<double> (HEIGHT);
    rad = fov * PI / 180.0;
    x_min = -a * tan(rad / 2.0);
    y_min = -tan(rad / 2.0);
    x_max = a * tan(rad / 2.0);
    y_max = tan(rad / 2.0);
    screen_width = x_max - x_min;
    cell_width = screen_width / static_cast<double> (WIDTH);
    screen_height = y_max - y_min;
    cell_height = screen_height / static_cast<double> (HEIGHT);
}

void init_ray(int row, int col, Ray (&r)[4]) {
    // calc
    double x = x_min + (2 * static_cast<double>(row) + 1) / 2.0 * cell_width;
    double y = y_min + (2 * static_cast<double>(col) + 1) / 2.0 * cell_height;
    double z = -1;

    // normalize
    //printf("B4 Origin: %f %f %f - Dest: %f %f %f\n", 0.0, 0.0, 0.0, x, y, z);
    normalize(x, y, z);

    // set Ray
    r[0].origin[0] = r[0].origin[1] = r[0].origin[2] = 0;
    r[0].direction[0] = (x - cell_width) / 4.0;
    r[0].direction[1] = y;
    r[0].direction[2] = z;
    normalize(r[0].direction[0], r[0].direction[1], r[0].direction[2]);

    r[1].origin[0] = r[1].origin[1] = r[1].origin[2] = 0;
    r[1].direction[0] = (x + cell_width) / 4.0;
    r[1].direction[1] = y;
    r[1].direction[2] = z;
    normalize(r[1].direction[0], r[1].direction[1], r[1].direction[2]);

    r[2].origin[0] = r[2].origin[1] = r[2].origin[2] = 0;
    r[2].direction[0] = x;
    r[2].direction[1] = (y - cell_height) / 4.0;
    r[2].direction[2] = z;
    normalize(r[2].direction[0], r[2].direction[1], r[2].direction[2]);

    r[3].origin[0] = r[3].origin[1] = r[3].origin[2] = 0;
    r[3].direction[0] = x;
    r[3].direction[1] = (y + cell_height) / 4.0;
    r[3].direction[2] = z;
    normalize(r[3].direction[0], r[3].direction[1], r[3].direction[2]);
}

// soft shadow feature
void init_light() {
    if (!SOFT_SHADOW) return;
    int initial_num_lights = num_lights;
    double color[3], center[3];

    for (int i = 0; i < initial_num_lights; i++) {
        color[0] = lights[i].color[0] / SUB_LIGHTS;
        color[1] = lights[i].color[1] / SUB_LIGHTS;
        color[2] = lights[i].color[2] / SUB_LIGHTS;

        center[0] = lights[i].position[0];
        center[1] = lights[i].position[1];
        center[2] = lights[i].position[2];

        lights[i].color[0] = color[0];
        lights[i].color[1] = color[1];
        lights[i].color[2] = color[2];

        for (int j = 0; j < SUB_LIGHTS - 1; j++) {
            lights[num_lights].color[0] = color[0];
            lights[num_lights].color[1] = color[1];
            lights[num_lights].color[2] = color[2];
            
            lights[num_lights].position[0] = center[0] + rand_val();
            lights[num_lights].position[1] = center[1] + rand_val();
            lights[num_lights].position[2] = center[2] + rand_val();

            num_lights++;
        }
    }
}

// --- intersect --- //
bool intersect_sphere(const Ray& r, double& t_min, int& idx) {
    t_min = (std::numeric_limits<double>::max)();
    idx = -1;
    double x0 = r.origin[0], y0 = r.origin[1], z0 = r.origin[2];
    double xd = r.direction[0], yd = r.direction[1], zd = r.direction[2];
    double xc, yc, zc, radius;
    double a, b, c;
    double t0, t1, t_temp;

    for (int i = 0; i < num_spheres; i++) {
        // get other sphere variables
        xc = spheres[i].position[0], yc = spheres[i].position[1], zc = spheres[i].position[2];
        radius = spheres[i].radius;
        // get values
        a = 1;
        b = 2 * (xd * (x0 - xc) + yd * (y0 - yc) + zd * (z0 - zc));
        c = pow((x0 - xc), 2) + pow((y0 - yc), 2) + pow((z0 - zc), 2) - pow(radius, 2);
        // get t0, t1
        t0 = -b + sqrt(b * b - 4 * c) / 2.0;
        t1 = -b - sqrt(b * b - 4 * c) / 2.0;
        t_temp = min(t0, t1);
        // get nearest t
        if (t_temp > EPSILON && t_temp < t_min) {
            idx = i;
            t_min = t_temp;
        }
    }
    return (idx != -1);
}

// moller-trumbore algorithm, no need for projection and plane precomputation
bool intersect_triangle(const Ray& r, double& t_min, int& idx) {
    double edge1[3], edge2[3], h[3], s[3], q[3];
    double a, f, u, v;
    t_min = (std::numeric_limits<double>::max)();
    idx = -1;
    for (int i = 0; i < num_triangles; i++) {
        Triangle triangle = triangles[i];
        // B - A edge1 = v1 - v0
        edge1[0] = triangle.v[1].position[0] - triangle.v[0].position[0];
        edge1[1] = triangle.v[1].position[1] - triangle.v[0].position[1];
        edge1[2] = triangle.v[1].position[2] - triangle.v[0].position[2];

        // C - A edge2 = v2 - v0
        edge2[0] = triangle.v[2].position[0] - triangle.v[0].position[0];
        edge2[1] = triangle.v[2].position[1] - triangle.v[0].position[1];
        edge2[2] = triangle.v[2].position[2] - triangle.v[0].position[2];

        // some more math
        cross_product(r.direction, edge2, h); // h = edge2 CROSS r.direction
        a = dot_product(edge1, h); // a = edge1 DOT h
        if (abs(a) < EPSILON) continue; // parallel

        // more math
        f = 1.0 / a;
        // s = r.origin - v0
        s[0] = r.origin[0] - triangle.v[0].position[0];
        s[1] = r.origin[1] - triangle.v[0].position[1];
        s[2] = r.origin[2] - triangle.v[0].position[2];
        // u = f * s DOT h
        u = f * dot_product(s, h);
        if (u < 0.0 || u > 1.0) continue;

        // q = edge1 CROSS s
        cross_product(s, edge1, q);
        // v = f * ray.direction DOT q
        v = f * dot_product(r.direction, q);
        if (v < 0.0 || u + v > 1.0) continue;
        
        // compute t = f * edge2 DOT q
        double t = f * dot_product(edge2, q);
        if (t > EPSILON && t < t_min) {
            t_min = t;
            idx = i;
        }
    }
    return (idx != -1);
}

// --- ray tracing --- //
void calc_color_ratio(int& idx, double pos[3], double(&ratio)[3]) {
    double edgeAB[3], edgeAC[3], normal[3], edgePA[3], edgePB[3], edgePC[3], normalPBC[3], normalPCA[3], normalPAB[3];

    // B - A edge1 = v1 - v0
    edgeAB[0] = triangles[idx].v[1].position[0] - triangles[idx].v[0].position[0];
    edgeAB[1] = triangles[idx].v[1].position[1] - triangles[idx].v[0].position[1];
    edgeAB[2] = triangles[idx].v[1].position[2] - triangles[idx].v[0].position[2];

    // C - A edge2 = v2 - v0
    edgeAC[0] = triangles[idx].v[2].position[0] - triangles[idx].v[0].position[0];
    edgeAC[1] = triangles[idx].v[2].position[1] - triangles[idx].v[0].position[1];
    edgeAC[2] = triangles[idx].v[2].position[2] - triangles[idx].v[0].position[2];

    // ABC area
    cross_product(edgeAB, edgeAC, normal);
    double areaABC = sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2)) * 0.5;

    // get edgePA
    edgePA[0] = triangles[idx].v[0].position[0] - pos[0];
    edgePA[1] = triangles[idx].v[0].position[1] - pos[1];
    edgePA[2] = triangles[idx].v[0].position[2] - pos[2];
    // get edgePA
    edgePB[0] = triangles[idx].v[1].position[0] - pos[0];
    edgePB[1] = triangles[idx].v[1].position[1] - pos[1];
    edgePB[2] = triangles[idx].v[1].position[2] - pos[2];
    // get edgePA
    edgePC[0] = triangles[idx].v[2].position[0] - pos[0];
    edgePC[1] = triangles[idx].v[2].position[1] - pos[1];
    edgePC[2] = triangles[idx].v[2].position[2] - pos[2];

    cross_product(edgePB, edgePC, normalPBC);
    double areaPBC = sqrt(pow(normalPBC[0], 2) + pow(normalPBC[1], 2) + pow(normalPBC[2], 2)) * 0.5;

    cross_product(edgePC, edgePA, normalPCA);
    double areaPCA = sqrt(pow(normalPCA[0], 2) + pow(normalPCA[1], 2) + pow(normalPCA[2], 2)) * 0.5;

    cross_product(edgePA, edgePB, normalPAB);
    double areaPAB = sqrt(pow(normalPAB[0], 2) + pow(normalPAB[1], 2) + pow(normalPAB[2], 2)) * 0.5;

    ratio[0] = areaPBC / areaABC;
    ratio[1] = areaPCA / areaABC;
    ratio[2] = areaPAB / areaABC;
}

bool in_shadow(const Light& light, const Vertex& point) {
    return true;
}

void calc_shading(Vertex& hitPoint, Light& light, double(&result)[3]) {
    double diffuse[3], specular[3], light_dir[3];
    // get light direction, normalize light_dir = light_pos - hitPoint.Pos
    light_dir[0] = light.position[0] - hitPoint.position[0];
    light_dir[1] = light.position[1] - hitPoint.position[1];
    light_dir[2] = light.position[2] - hitPoint.position[2];
    normalize(light_dir[0], light_dir[1], light_dir[2]);

    // get diffuse diff = max(light_dir DOT hitpoint.normal, 0)
    double diff = max(dot_product(light_dir, hitPoint.normal), 0.0);
    // diffuse = light.color * (hitPoint.color_diff * diff)
    diffuse[0] = light.color[0] * hitPoint.color_diffuse[0] * diff;
    diffuse[1] = light.color[1] * hitPoint.color_diffuse[1] * diff;
    diffuse[2] = light.color[2] * hitPoint.color_diffuse[2] * diff;

    // get view dir, normalize view_dir = -hitPoint.pos
    double view_dir[3];
    view_dir[0] = -hitPoint.position[0];
    view_dir[1] = -hitPoint.position[1];
    view_dir[2] = -hitPoint.position[2];
    normalize(view_dir[0], view_dir[1], view_dir[2]);

    // get reflect dir r_dir = 2 * (light_dir DOT hitpoint.Normal) * hitPoint.normal - light_dir
    double r_dir[3];
    reflect_dir(light_dir, hitPoint.normal, r_dir);

    // get specular spec = max(view_dir DOT r_dir, shininess) 
    double spec = pow(max(dot_product(view_dir, r_dir), 0.0), hitPoint.shininess);
    // specular = light.color * (hitPointt.color_specular * spec)
    specular[0] = light.color[0] * hitPoint.color_specular[0] * spec;
    specular[1] = light.color[1] * hitPoint.color_specular[1] * spec;
    specular[2] = light.color[2] * hitPoint.color_specular[2] * spec;

    // result = diffuse + specular
    result[0] = diffuse[0] + specular[0];
    result[1] = diffuse[1] + specular[1];
    result[2] = diffuse[2] + specular[2];
}

void calc_radiance(const Ray& ray, int times, double (&radiance)[3]) {
    double t1, t2;
    int idx1, idx2;
    bool hit_sphere = intersect_sphere(ray, t1, idx1);
    bool hit_triangle = intersect_sphere(ray, t2, idx2);

    // nothing hit
    if (!hit_sphere && !hit_triangle) {
        radiance[0] = radiance[1] = radiance[2] = 1.0;
        return;
    }

    Vertex hit_point;
    // if ray hits only sphere or hits sphere before triangle
    if ((hit_sphere && !hit_triangle) && (hit_sphere && hit_triangle && t1 < t2)) {
        printf("HIT SPHERE\n");
        // get hit point pos
        double hit_point_pos[3] = { ray.origin[0] + ray.direction[0] * t1, ray.origin[1] + ray.direction[1] * t1, ray.origin[2] + ray.direction[2] * t1 };
        // get normal for hit point and normalize
        double normal[3] = { hit_point_pos[0] - spheres[idx1].position[0], hit_point_pos[1] - spheres[idx1].position[1], hit_point_pos[2] - spheres[idx1].position[2] };
        normalize(normal[0], normal[1], normal[2]);

        // create hitpoint
        hit_point = {
            { hit_point_pos[0], hit_point_pos[1], hit_point_pos[2] },
            {spheres[idx1].color_diffuse[0], spheres[idx1].color_diffuse[1], spheres[idx1].color_diffuse[2] },
            {spheres[idx1].color_specular[0], spheres[idx1].color_specular[1], spheres[idx1].color_specular[2] },
            { normal[0], normal[1], normal[2] },
            spheres[idx1].shininess
        };
    } // hit triangle
    else {
        double ratio[3];
        double hit_point_pos[3] = { ray.origin[0] + ray.direction[0] * t2, ray.origin[1] + ray.direction[1] * t2, ray.origin[2] + ray.direction[2] * t2 };
        calc_color_ratio(idx2, hit_point_pos, ratio);
        Triangle triangle = triangles[idx2];
        
        // get diffuse
        double diffuse[3];
        diffuse[0] = triangle.v[0].color_diffuse[0] * ratio[0] + triangle.v[1].color_diffuse[0] * ratio[1] + triangle.v[2].color_diffuse[0] * ratio[2];
        diffuse[1] = triangle.v[0].color_diffuse[1] * ratio[0] + triangle.v[1].color_diffuse[1] * ratio[1] + triangle.v[2].color_diffuse[1] * ratio[2];
        diffuse[2] = triangle.v[0].color_diffuse[2] * ratio[0] + triangle.v[1].color_diffuse[2] * ratio[1] + triangle.v[2].color_diffuse[2] * ratio[2];

        // get specular
        double specular[3];
        specular[0] = triangle.v[0].color_specular[0] * ratio[0] + triangle.v[1].color_specular[0] * ratio[1] + triangle.v[2].color_specular[0] * ratio[2];
        specular[1] = triangle.v[0].color_specular[1] * ratio[0] + triangle.v[1].color_specular[1] * ratio[1] + triangle.v[2].color_specular[1] * ratio[2];
        specular[2] = triangle.v[0].color_specular[2] * ratio[0] + triangle.v[1].color_specular[2] * ratio[1] + triangle.v[2].color_specular[2] * ratio[2];

        // get normal
        double normal[3];
        normal[0] = triangle.v[0].normal[0] * ratio[0] + triangle.v[1].normal[0] * ratio[1] + triangle.v[2].normal[0] * ratio[2];
        normal[1] = triangle.v[0].normal[1] * ratio[0] + triangle.v[1].normal[1] * ratio[1] + triangle.v[2].normal[1] * ratio[2];
        normal[2] = triangle.v[0].normal[2] * ratio[0] + triangle.v[1].normal[2] * ratio[1] + triangle.v[2].normal[2] * ratio[2];

        // shininess
        double shininess = triangle.v[0].shininess * ratio[0] + triangle.v[1].shininess * ratio[1] + triangle.v[2].shininess * ratio[2];

        // create hitpoint
        hit_point = {
            { hit_point_pos[0], hit_point_pos[1], hit_point_pos[2] },
            { diffuse[0], diffuse[1], diffuse[2] },
            { specular[0], specular[1], specular[2] },
            { normal[0], normal[1], normal[2] },
            shininess
        };
    }

    double current_ray_color[3] = { 0, 0, 0 };
    // for each light
    for (int i = 0; i < num_lights; i++) {
        if (!in_shadow(lights[i], hit_point)) {
            double result[3];
            calc_shading(hit_point, lights[i], result);
            current_ray_color[0] += result[0];
            current_ray_color[1] += result[1];
            current_ray_color[2] += result[2];
        }
    }

    // base case - max reflections reached, update radiance with only current_ray_color
    if (times >= MAX_REFLECTION) {
        radiance[0] = current_ray_color[0];
        radiance[1] = current_ray_color[1];
        radiance[2] = current_ray_color[2];
        return;
    }

    times++;

    // reflect ray, recurse
    double negate_ray_dir[3] = { -ray.direction[0] , -ray.direction[1] , -ray.direction[2] };
    double r_dir[3];
    reflect_dir(negate_ray_dir, hit_point.normal, r_dir);
    Ray reflect_ray = {
        { hit_point.position[0], hit_point.position[1], hit_point.position[2] },
        { r_dir[0], r_dir[1], r_dir[2] },
    };
    double r_color[3];
    calc_radiance(reflect_ray, times, r_color);

    // done with recursion, update radiance
    radiance[0] = current_ray_color[0] * (1 - REFLECT_RATIO) + r_color[0] * REFLECT_RATIO;
    radiance[1] = current_ray_color[1] * (1 - REFLECT_RATIO) + r_color[1] * REFLECT_RATIO;
    radiance[2] = current_ray_color[2] * (1 - REFLECT_RATIO) + r_color[2] * REFLECT_RATIO;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    init_screen();
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      Ray r[4];
      init_ray(x, y, r);
      double color[3];
      for (int k = 0; k < 4; k++) {
          double color_temp[3];
          calc_radiance(r[k], 0, color_temp);
          color[0] += color_temp[0];
          color[1] += color_temp[1];
          color[2] += color_temp[2];
      }
      color[0] /= 4.0;
      color[1] /= 4.0;
      color[2] /= 4.0;
      plot_pixel(x, y, static_cast<int>(color[0] * 255), static_cast<int>(color[1] * 255), static_cast<int>(color[2] * 255));
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}


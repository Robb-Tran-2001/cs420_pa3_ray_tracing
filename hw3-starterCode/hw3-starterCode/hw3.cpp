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

double aspect_ratio, rad, x_min, y_min, x_max, y_max, screen_width, cell_width, screen_height, cell_height;

// helper struct to help with all the annoying double[3] math
struct Vector {
    double x, y, z;
    Vector(double x1 = 0, double y1 = 0, double z1 = 0) { x = x1; y = y1; z = z1; }
    Vector operator+(const Vector& v) const { return Vector(x + v.x, y + v.y, z + v.z); }
    Vector operator-(const Vector& v) const { return Vector(x - v.x, y - v.y, z - v.z); }
    Vector operator-() const { return Vector(-x, -y, -z); }
    Vector operator*(double a) const { return Vector(x * a, y * a, z * a); }
    Vector operator/(double a) const { return Vector(x / a, y / a, z / a); }
    Vector mult(const Vector& a) const { return Vector(x * a.x, y * a.y, z * a.z); }
    Vector& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vector& a) const { return x * a.x + y * a.y + z * a.z; }
    double to_norm() { return sqrt(x * x + y * y + z * z); }
    Vector cross(const Vector& a) const { return Vector(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x); }
};

struct Vertex
{
  Vector position;
  Vector color_diffuse;
  Vector color_specular;
  Vector normal;
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  Vector position;
  Vector color_diffuse;
  Vector color_specular;
  double shininess;
  double radius;
};

struct Light
{
  Vector position;
  Vector color;
};

struct Ray
{
    Vector o;
    Vector d;
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];

Vector ambient_light;
Vector background_color = Vector(1.0, 1.0, 1.0);

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;
int num_rays = 0;

// --- helper --- //
void clamp(Vector color) {
    if (color.x > 1.0) color.x = 1.0;
    if (color.y > 1.0) color.y = 1.0;
    if (color.z > 1.0) color.z = 1.0;
}

bool in_range(double v, double low, double high) {
    return (v >= low && v <= high);
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

Vector reflect_dir(const Vector& I, const Vector& N) //I and N should be normalized
{
    Vector reflect = I - N * 2 * N.dot(I);
    return reflect.norm();
}

Vector refract_dir(const Vector& I, const Vector& N, const double& ratio) {
    double k = 1.0 - ratio * ratio * (1.0 - N.dot(I) * N.dot(I));
    if (k < 0.0)
        return Vector(0, 0, 0);
    return I * ratio - N * (ratio * N.dot(I) + sqrt(k));
}

// --- init --- //
void init_screen() {
    // set up
    aspect_ratio = static_cast<double> (WIDTH) / static_cast<double> (HEIGHT);
    rad = fov * PI / 180.0;
    x_min = -aspect_ratio * tan(rad / 2.0);
    y_min = -tan(rad / 2.0);
    x_max = aspect_ratio * tan(rad / 2.0);
    y_max = tan(rad / 2.0);
    screen_width = x_max - x_min;
    cell_width = screen_width / static_cast<double> (WIDTH);
    screen_height = y_max - y_min;
    cell_height = screen_height / static_cast<double> (HEIGHT);
}

void init_ray(int row, int col, Ray (&r)[4]) {
    // calc
    double x = x_min + (2 * row + 1) / 2.0 * cell_width;
    double y = y_min + (2 * col + 1) / 2.0 * cell_height;
    double z = -1;

    // normalize
    //printf("B4 Origin: %f %f %f - Dest: %f %f %f\n", 0.0, 0.0, 0.0, x, y, z);
    r[0] = { Vector(0, 0, 0), Vector(x - cell_width/4.0, y, z).norm() };
    r[1] = { Vector(0, 0, 0), Vector(x + cell_width / 4.0, y, z).norm() };
    r[2] = { Vector(0, 0, 0), Vector(x, y - cell_height / 4.0, z).norm() };
    r[3] = { Vector(0, 0, 0), Vector(x, y + cell_height / 4.0, z).norm() };
}

// soft shadow feature
void init_light() {
    if (!SOFT_SHADOW) return;
    int initial_num_lights = num_lights;

    for (int i = 0; i < initial_num_lights; i++) {
        Vector color = lights[i].color / SUB_LIGHTS;
        Vector center = lights[i].position;

        lights[i].color = color;

        for (int j = 0; j < SUB_LIGHTS - 1; j++) {
            lights[num_lights].color = color;
            lights[num_lights].position = Vector(center.x + rand_val(), center.y + rand_val(), center.z + rand_val());
            num_lights++;
        }
    }
}

// --- intersect --- //
bool intersect_sphere(const Ray& r, double& t_min, int& idx) {
    t_min = (std::numeric_limits<double>::max)();
    idx = -1;

    for (int i = 0; i < num_spheres; i++) {
        // get other sphere variables
        Vector oc = Vector(r.o - spheres[i].position);
        // get a, b, c
        double a = 1;
        double b = 2 * oc.dot(r.d);
        double c = oc.dot(oc) - pow(spheres[i].radius, 2);
        double delta = b * b - 4 * a * c;
        if (delta < 0) continue;
        // get t0, t1
        double t0 = (-b + sqrt(delta)) / 2;
        double t1 = (-b - sqrt(delta)) / 2;
        double t_temp = min(t0, t1);
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
    t_min = (std::numeric_limits<double>::max)();
    idx = -1;
    for (int i = 0; i < num_triangles; i++) {
        Triangle triangle = triangles[i];
        // B - A edge1 = v1 - v0
        Vector edge1 = triangle.v[1].position - triangle.v[0].position;

        // C - A edge2 = v2 - v0
        Vector edge2 = triangle.v[2].position - triangle.v[0].position;

        // some more math
        Vector s = r.o - triangle.v[0].position;
        Vector s1 = r.d.cross(edge2);
        Vector s2 = s.cross(edge1);
        double k = s1.dot(edge1);
        if (abs(k) < EPSILON) continue;
       
        // barycentric math
        double t_temp = 1 / k * s2.dot(edge2);
        double b1 = 1 / k * s1.dot(s);
        double b2 = 1 / k * s2.dot(r.d);
        double b3 = 1 - b1 - b2;

        if (t_temp < EPSILON || !in_range(b1, 0, 1) || !in_range(b2, 0, 1) || !in_range(b3, 0, 1))
            continue;
        if (t_temp < t_min) {
            t_min = t_temp;
            idx = i;
        }
    }
    return (idx != -1);
}

// --- ray tracing --- //
void calc_color_ratio(int& idx, Vector pos, double(&ratio)[3]) {
    // get vertices
    Vertex vertex0 = triangles[idx].v[0];
    Vertex vertex1 = triangles[idx].v[1];
    Vertex vertex2 = triangles[idx].v[1];

    // get 2 edges
    Vector ab = vertex1.position - vertex0.position;

    // C - A edge2 = v2 - v0
    Vector ac = vertex2.position - vertex0.position;

    // ABC area
    double areaABC = ab.cross(ac).to_norm() * 0.5;

    // get edgePA
    Vector pa = vertex0.position - pos;

    // get edgePB
    Vector pb = vertex1.position - pos;

    // get edgePC
    Vector pc = vertex2.position - pos;

    double areaPBC = pb.cross(pc).to_norm() * 0.5;
    double areaPCA = pc.cross(pa).to_norm() * 0.5;
    double areaPAB = pa.cross(pb).to_norm() * 0.5;

    ratio[0] = areaPBC / areaABC;
    ratio[1] = areaPCA / areaABC;
    ratio[2] = areaPAB / areaABC;
}

bool in_shadow(const Light& light, const Vertex& point) {
    Vector dir = (light.position - point.position).norm();
    Ray shadow_ray = { point.position + dir * 5 * EPSILON, dir };

    // check if shadow ray intersects with any shapes, if it does, no shadow
    double t1, t2; int i1, i2;
    bool hit_sphere = intersect_sphere(shadow_ray, t1, i1);
    bool hit_triangle = intersect_triangle(shadow_ray, t2, i2);

    // hit nothing
    if (!hit_sphere && !hit_triangle)
        return false;

    // hit point should not exceed light position
    Vector hit_pos;
    if ((hit_sphere && !hit_triangle) || (hit_sphere && hit_triangle && t1 < t2))
        hit_pos = shadow_ray.o + shadow_ray.d * t1;
    else
        hit_pos = shadow_ray.o + shadow_ray.d * t2;

    double point_to_light_len = (point.position - light.position).to_norm();
    double point_to_hit_len = (point.position - hit_pos).to_norm();
    if (point_to_hit_len - point_to_light_len > EPSILON)
        return false;

    return true;
}

Vector calc_shading(Vertex& hit_point, Light& light) {
    Vector diffuse, specular;
    // get light direction, normalize light_dir = light_pos - hitPoint.Pos
    Vector light_dir = (light.position - hit_point.position).norm();
    // get diffuse diff = max(light_dir DOT hitpoint.normal, 0)
    double diff = max(light_dir.dot(hit_point.normal), 0.0);
    // diffuse = light.color * (hitPoint.color_diff * diff)
    diffuse = light.color.mult(hit_point.color_diffuse * diff);

    // get view dir, normalize view_dir = -hitPoint.pos
    Vector view_dir = (-hit_point.position).norm();
    // get reflect dir r_dir = 2 * (light_dir DOT hitpoint.Normal) * hitPoint.normal - light_dir
    Vector r_dir = reflect_dir(-light_dir, hit_point.normal);
    // get specular spec = max(view_dir DOT r_dir, shininess) 
    double spec = pow(max(view_dir.dot(r_dir), 0.0), hit_point.shininess);
    // specular = light.color * (hitPointt.color_specular * spec)
    specular = light.color.mult(hit_point.color_specular * spec);

    // result = diffuse + specular
    return diffuse + specular;
}

Vector calc_radiance(const Ray& r, int times) {
    double t1, t2;
    int idx1, idx2;
    bool hit_sphere = intersect_sphere(r, t1, idx1);
    bool hit_triangle = intersect_sphere(r, t2, idx2);

    // nothing hit
    if (!hit_sphere && !hit_triangle) {
        return background_color;
    }

    Vertex hit_point;
    // if ray hits only sphere or hits sphere before triangle
    if ((hit_sphere && !hit_triangle) && (hit_sphere && hit_triangle && t1 < t2)) {
        Sphere sphere = spheres[idx1];
        hit_point = {
            r.o + r.d * t1,
            sphere.color_diffuse,
            sphere.color_specular,
            (r.o + r.d * t1 - sphere.position).norm(),
            sphere.shininess
        };
    } // hit triangle
    else {
        double ratio[3];
        Vector pos = r.o + r.d * t2;
        calc_color_ratio(idx2, pos, ratio);
        Triangle triangle = triangles[idx2];
        Vector diffuse = triangle.v[0].color_diffuse * ratio[0] + triangle.v[1].color_diffuse * ratio[1] + triangle.v[2].color_diffuse * ratio[2];

        Vector specular = triangle.v[0].color_specular * ratio[0] + triangle.v[1].color_specular * ratio[1] + triangle.v[2].color_specular * ratio[2];

        Vector normal = triangle.v[0].normal * ratio[0] + triangle.v[1].normal * ratio[1] + triangle.v[2].normal * ratio[2];

        double shininess = triangle.v[0].shininess * ratio[0] + triangle.v[1].shininess * ratio[1] + triangle.v[2].shininess * ratio[2];
        hit_point = { pos, diffuse, specular, normal.norm(), shininess };
    }

    Vector current_ray_color = Vector(0, 0, 0);
    // for each light
    for (int i = 0; i < num_lights; i++) {
        if (!in_shadow(lights[i], hit_point)) {
            current_ray_color = current_ray_color + calc_shading(hit_point, lights[i]);
        }
    }

    // base case - max reflections reached, update radiance with only current_ray_color
    if (times >= MAX_REFLECTION) {
        return current_ray_color;
    }

    times++;

    // reflect ray, recurse
    Vector r_dir = reflect_dir(r.d, hit_point.normal);
    Ray reflect_ray = { hit_point.position, r_dir };
    Vector reflect_color = calc_radiance(reflect_ray, times);
    
    return current_ray_color * (1 - REFLECT_RATIO) + reflect_color * REFLECT_RATIO;
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
        if (ANTI_ALIASING)
        {
            Ray rays[4];
            Vector color;
            init_ray(x, y, rays);
            for (int k = 0; k < 4; k++)
            {
                color = color + calc_radiance(rays[k], 0);
            }
            color = color / 4;
            plot_pixel(x, y, (int)(color.x * 255), (int)(color.y * 255), (int)(color.z * 255));
        }
        else
        {
            Ray ray;
            Vector color;
            double xx = x_min + (2 * x + 1) * cell_width / 2.0f;
            double yy = y_min + (2 * y + 1) * cell_height / 2.0f;
            double zz = -1;

            ray = { Vector(0,0,0),  Vector(xx,yy,zz).norm() };
            color = calc_radiance(ray, 1) + ambient_light;
            clamp(color);
            plot_pixel(x, y, (int)(color.x * 255), (int)(color.y * 255), (int)(color.z * 255));
        }
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

void parse_doubles(FILE* file, const char *check, Vector& p)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf", &p.x, &p.y, &p.z);
  printf("%s %lf %lf %lf\n",check, p.x, p.y, p.z);
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
    init_light();
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


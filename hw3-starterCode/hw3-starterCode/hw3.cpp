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

#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

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

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

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

void cross_product(const double a[1], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

void init_ray(int row, int col, Ray &r) {
    // set up
    double a = static_cast<double> (WIDTH) / static_cast<double> (HEIGHT);
    double rad = fov * PI / 180.0;
    double x_min = -a * tan(rad / 2.0);
    double y_min = -tan(rad / 2.0);
    double x_max = a * tan(rad / 2.0);
    double y_max = tan(rad / 2.0);
    
    // calc
    double x = x_min + static_cast<double> (row) * (x_max - x_min) / static_cast<double> (WIDTH);
    double y = y_min + static_cast<double> (col) * (y_max - y_min) / static_cast<double> (HEIGHT);
    double z = -1;

    // normalize
    //printf("B4 Origin: %f %f %f - Dest: %f %f %f\n", 0.0, 0.0, 0.0, x, y, z);
    normalize(x, y, z);

    // set Ray
    r.origin[0] = r.origin[1] = r.origin[2] = 0;
    r.direction[0] = x;
    r.direction[1] = y;
    r.direction[2] = z;
}

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
    double edge1[3], edge2[3], normal[3], h[3], s[3], q[3];
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

//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      Ray r;
      init_ray(x, y, r);
      //printf("Origin: %f %f %f - Dest: %f %f %f\n", r.origin[0], r.origin[1], r.origin[2], r.destination[0], r.destination[1], r.destination[2]);
      double ts, tt;
      int idx_s, idx_t;
      bool intersect = intersect_sphere(r, ts, idx_s);
      bool intersect2 = intersect_triangle(r, tt, idx_t);
      if (intersect2) printf("Origin: %f %f %f - Dest: %f %f %f\n Intersects triangle %d at %f \n ----------- \n", r.origin[0], r.origin[1], r.origin[2], r.direction[0], r.direction[1], r.direction[2], idx_t, tt);
      // else printf("DOESN'T INTERSECT ANY\n");
      plot_pixel(x, y, x % 256, y % 256, (x+y) % 256);
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


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <time.h> // for seeding rand()
#ifdef USEGLEW
#include <GL/glew.h>
#endif
//  OpenGL with prototypes for glext
#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/glut.h>
// Tell Xcode IDE to not gripe about OpenGL deprecation
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#else
#include <GL/glut.h>
#endif

//#include "Clarkson-Delaunay.h"
int *BuildTriangleIndexList (void *pointList, float factor, int numberOfInputPoints, int numDimensions, int clockwise, int *numTriangleVertices );

//  Default resolution
//  For Retina displays compile with -DRES=2
#ifndef RES
#define RES 1
#endif

//  Macro for sin & cos in degrees
#define PI 3.14159265
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))
#define Acos(th) acos(th)*180/PI
#define Asin(th) asin(th)*180/PI


/*
 *  Convenience routine to output raster text
 *  Use VARARGS to make this more flexible
 */
#define LEN 8192  //  Maximum length of text string
void Print(const char* format , ...)
{
   char    buf[LEN];
   char*   ch=buf;
   va_list args;
   //  Turn the parameters into a character string
   va_start(args,format);
   vsnprintf(buf,LEN,format,args);
   va_end(args);
   //  Display the characters one at a time at the current raster position
   while (*ch)
      glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*ch++);
}

/*
 *  Check for OpenGL errors
 */
void ErrCheck(const char* where)
{
   int err = glGetError();
   if (err) fprintf(stderr,"ERROR: %s [%s]\n",gluErrorString(err),where);
}

/*
 *  Print message to stderr and exit
 */
void Fatal(const char* format , ...)
{
   va_list args;
   va_start(args,format);
   vfprintf(stderr,format,args);
   va_end(args);
   exit(1);
}


/*
 *  Draw a tree
 */
void tree()
{
   //  Save transformation
   glPushMatrix();
   for (int th=0;th<=360;th+=1){  
      //  Save transformation
      glPushMatrix();
      glRotatef(th,0,1,0);
      // tree
      //glColor3f(0,.5,.5);
      glBegin(GL_LINE_STRIP);
      glNormal3d(0,-1,0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(0.4,0.0,0.0);
      glNormal3d(1/1.414,1/1.414,0); // diag
      glVertex3f(0.2,0.2,0.0);
      glNormal3d(0,-1,0);
      glVertex3f(0.3,0.2,0.0);
      glNormal3d(1/1.414,1/1.414,0); // diag
      glVertex3f(0.1,0.4,0.0);
      glNormal3d(0,-1,0);
      glVertex3f(0.2,0.4,0.0);
      glNormal3d(1/1.414,1/1.414,0); // diag
      glVertex3f(0.0,0.6,0.0);
      glEnd();
      //  Undo transformations
      glPopMatrix();
   }
   //  Undo transformations
   glPopMatrix();
}

/*
 *  Draw vertex in polar coordinates with normal
 *  Used for ball function
 */
void Vertex(double th,double ph)
{
   double x = Sin(th)*Cos(ph);
   double y = Cos(th)*Cos(ph);
   double z =         Sin(ph);
   //  For a sphere at the origin, the position
   //  and normal vectors are the same
   glNormal3d(x,y,z);
   glVertex3d(x,y,z);
}

/*
 *  Draw a ball
 *     at (x,y,z)
 *     radius (r)
 */
void ball(double x,double y,double z,double r)
{
   //  Save transformation
   glPushMatrix();
   //  Offset, scale and rotate
   glTranslated(x,y,z);
   glScaled(r,r,r);
   //  White ball with yellow specular
   int emission = 0;
   float yellow[]   = {1.0,1.0,0.0,1.0};
   float Emission[] = {0.0,0.0,0.01*emission,1.0};
   glColor3f(1,1,1);
   int shiny = 1;
   glMaterialf(GL_FRONT,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT,GL_SPECULAR,yellow);
   glMaterialfv(GL_FRONT,GL_EMISSION,Emission);
   //  Bands of latitude
   int inc = 10;
   for (int ph=-90;ph<90;ph+=inc)
   {
      glBegin(GL_QUAD_STRIP);
      for (int th=0;th<=360;th+=2*inc)
      {
         Vertex(th,ph);
         Vertex(th,ph+inc);
      }
      glEnd();
   }
   //  Undo transofrmations
   glPopMatrix();
}


void displayScene_Tree_test(){
   glColor3f(0,.5,.2);
   tree();

   glPushMatrix();
   glTranslatef(1 ,0,0);
   glScalef(1,1.25,.75);
   glColor3f(0,.5,.2);
   tree();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(1 ,0,0);
   glScalef(1.75,1,1);
   glColor3f(0,.5,.2);
   tree();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(1 ,0,1);
   glScalef(1,1,1);
   glColor3f(0,.5,.2);
   tree();
   glPopMatrix();

   glPushMatrix();
   glTranslatef(1 ,0,-1);
   glScalef(1,1,1);
   glColor3f(0,.5,.2);
   tree();
   glPopMatrix();
}

//
//  Reverse n bytes
//
static void Reverse(void* x,const int n)
{
   int k;
   char* ch = (char*)x;
   for (k=0;k<n/2;k++)
   {
      char tmp = ch[k];
      ch[k] = ch[n-1-k];
      ch[n-1-k] = tmp;
   }
}

//
//  Load texture from BMP file
//
unsigned int LoadTexBMP(const char* file)
{
   unsigned int   texture;    // Texture name
   FILE*          f;          // File pointer
   unsigned short magic;      // Image magic
   unsigned int   dx,dy,size; // Image dimensions
   unsigned short nbp,bpp;    // Planes and bits per pixel
   unsigned char* image;      // Image data
   unsigned int   off;        // Image offset
   unsigned int   k;          // Counter
   unsigned int   max;        // Maximum texture dimensions

   //  Open file
   f = fopen(file,"rb");
   if (!f) Fatal("Cannot open file %s\n",file);
   //  Check image magic
   if (fread(&magic,2,1,f)!=1) Fatal("Cannot read magic from %s\n",file);
   if (magic!=0x4D42 && magic!=0x424D) Fatal("Image magic not BMP in %s\n",file);
   //  Read header
   if (fseek(f,8,SEEK_CUR) || fread(&off,4,1,f)!=1 ||
       fseek(f,4,SEEK_CUR) || fread(&dx,4,1,f)!=1 || fread(&dy,4,1,f)!=1 ||
       fread(&nbp,2,1,f)!=1 || fread(&bpp,2,1,f)!=1 || fread(&k,4,1,f)!=1)
     Fatal("Cannot read header from %s\n",file);
   //  Reverse bytes on big endian hardware (detected by backwards magic)
   if (magic==0x424D)
   {
      Reverse(&off,4);
      Reverse(&dx,4);
      Reverse(&dy,4);
      Reverse(&nbp,2);
      Reverse(&bpp,2);
      Reverse(&k,4);
   }
   //  Check image parameters
   glGetIntegerv(GL_MAX_TEXTURE_SIZE,(int*)&max);
   if (dx<1 || dx>max) Fatal("%s image width %d out of range 1-%d\n",file,dx,max);
   if (dy<1 || dy>max) Fatal("%s image height %d out of range 1-%d\n",file,dy,max);
   if (nbp!=1)  Fatal("%s bit planes is not 1: %d\n",file,nbp);
   if (bpp!=24) Fatal("%s bits per pixel is not 24: %d\n",file,bpp);
   if (k!=0)    Fatal("%s compressed files not supported\n",file);
#ifndef GL_VERSION_2_0
   //  OpenGL 2.0 lifts the restriction that texture size must be a power of two
   for (k=1;k<dx;k*=2);
   if (k!=dx) Fatal("%s image width not a power of two: %d\n",file,dx);
   for (k=1;k<dy;k*=2);
   if (k!=dy) Fatal("%s image height not a power of two: %d\n",file,dy);
#endif

   //  Allocate image memory
   size = 3*dx*dy;
   image = (unsigned char*) malloc(size);
   if (!image) Fatal("Cannot allocate %d bytes of memory for image %s\n",size,file);
   //  Seek to and read image
   if (fseek(f,off,SEEK_SET) || fread(image,size,1,f)!=1) Fatal("Error reading data from image %s\n",file);
   fclose(f);
   //  Reverse colors (BGR -> RGB)
   for (k=0;k<size;k+=3)
   {
      unsigned char temp = image[k];
      image[k]   = image[k+2];
      image[k+2] = temp;
   }

   //  Sanity check
   ErrCheck("LoadTexBMP");
   //  Generate 2D texture
   glGenTextures(1,&texture);
   glBindTexture(GL_TEXTURE_2D,texture);
   //  Copy image
   glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,dx,dy,0,GL_RGB,GL_UNSIGNED_BYTE,image);
   if (glGetError()) Fatal("Error in glTexImage2D %s %dx%d\n",file,dx,dy);
   //  Scale linearly when image size doesn't match
   glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);

   //  Free image memory
   free(image);
   //  Return texture name
   return texture;
}




int mode=0;       //  Projection mode
int move=1;       //  Move light
int th=0;         //  Azimuth of view angle
int ph=15;        //  Elevation of view angle
int fov=55;       //  Field of view (for perspective)
double asp=1;     //  Aspect ratio
double dim=5.0;   //  Size of world
int view_pts=0;   //  Flag to render points

// for First Person
// eye position
double Ex_fp=-1;
double Ey_fp=.5;
double Ez_fp=1;

// Light values
int light     =   1; 
int one       =   1;  // Unit value
int distance  =  20;  // Light distance
int inc       =  10;  // Ball increment
int emission  =   0;  // Emission intensity (%)
int ambient   =  40;  // Ambient intensity (%)
int diffuse   =  70;  // Diffuse intensity (%)
int specular  =  40;  // Specular intensity (%)
int shininess =   4;  // Shininess (power of two)
float shiny   =   1;  // Shininess (value)

int omega = 90; // argument of perihelion (angle where min=0) 
int incl  = 90; // inclination angle of plane
int Omega = 0; // Longitude of ascending node
int zh        =  90;  // Light azimuth
int ylight    =  45;  // Elevation of light
int minute    =   0;  // There are 1440 minutes in a day



// crater variables
int numCraters = 0;
int numPts = 0;
int num_vtx;

unsigned int texture; // Texture names

#define CANVAS_LEN 10
#define MAX_DIM 5000 // km or something
#define N 2097152 // 2^21 pts ~ 25mb / ~455 craters 
int pts[N*2]; // 2 planar coordinates per pt
int pts_z[N]; // 1 z value per pt
int color[N]; // 1 color per pt

#define CRATER_ANGULAR_RESOLUTION 5 //int

typedef struct {double r; double z;} rz; // pt for cross-section profile

rz cross0[11] = { // just one for now
   {.04, -.9},
   {.06, -1.0},
   {.23, -1.0},
   {.40, -.93},
   {.56, -.80},
   {.60, -.67},
   {.60, -.53},
   {.80,  .20},
   {.84,  .27},
   {.88,  .27},
   {1.0,  0.0}
}; 
rz cross1[9] = { 
   {.23, -1.0},
   {.40, -.93},
   {.56, -.80},
   {.60, -.67},
   {.60, -.53},
   {.80,  .05},
   {.84,  .1},
   {.88,  .05},
   {1.0,  0.0}
}; 

rz* cross[2] = {cross0, cross1};
double cross_ctr[2] = {-.8, -1};
int cross_len[2] = {11, 9};


typedef struct {int x; int y;} xy; // pt for center

xy next_center = {0, 0};
int next_angle = 0;

typedef struct {int min; int max;} range;

range r_diameter = {500, 1500};
range r_eccentricity = {50, 100};
int r_eccentricity_lim = 100;
range r_depth = {100, 500};
range r_dim = {-MAX_DIM, MAX_DIM};
range r_color = {40, 50};
int r_color_lim = 100;
int r_num_cross = 2;

#define DEFAULT_GRAY 50 

int eccen_angle_incr = 5;

typedef struct {
   // holds crater parameters
   // shape generation
   // ellipse (polar): r(theta) = a* b_n / sqrt(1-(1-b_n^2)cos(theta)^2)
   int a;     // scalar / length of semimajor axis "diameter" (0,50] in km
   double b_n;   // normalized semiminor axis (b/a) "eccentricity" (0,1]
   int depth; // depth scalar in km
   int cross; // cross-sectional shape - later expand to chebyshev
   //int text;  // texture index
   // location / orientation
   xy* center;
   int az;  // azimuthal orientation of semimajor axis
} crater;

void remove_point(int index) {
   int i;
   for(i = index; i < numPts - 2; i++) {
      pts[2*i] = pts[2*i+2];
      pts[2*i+1] = pts[2*i+3];
      pts_z[i] = pts_z[i+1];

      color[i] = color[i+1];
   }
   numPts--;
}

void find_overlaps(crater* c){
   int xp,yp, // pt
       xc,yc;
       //xe,ye; // center

   double d,r, theta;

   xc = c->center->x;
   yc = c->center->y;

   for(int i = 0; i < numPts; i++) {
      xp = pts[2*i];
      yp = pts[2*i + 1];

      // distance between center & pt
      d = sqrt(((xp-xc)*(xp-xc) + (yp-yc)*(yp-yc)));
      // angle
      theta = Asin(abs(yp-yc)/d) + 180;

      // distance to edge of footprint
      r = c->a*c->b_n/sqrt(Sin(theta)*Sin(theta) + c->b_n*c->b_n*Cos(theta)*Cos(theta));

      // xe = (int)floor(r*Cos(theta)) + c->center->x;
      // ye = (int)floor(r*Sin(theta)) + c->center->y;

      if(d <= r){
         //printf("Deletion candidate: Edge: (%d, %d) Point: (%d, %d) Center: (%d, %d)\n", xe,ye, xp,yp, xc,yc);
         //printf("       Theta: %f, D: %f, R: %f\n", theta, d, r);
         remove_point(i);
         i--;
      }
   }
}

int get_random_in_range(range r){
   return (r.min + (r.max - r.min)*(float)rand()/(float)(RAND_MAX));
}

void draw_crater(crater* c){
   // Draw crater here
   int p,Theta,theta,x,y,z;
   double r,x_d,y_d,z_d;
   
   find_overlaps(c);

   for(Theta=0; Theta<360; Theta+=CRATER_ANGULAR_RESOLUTION){
      //printf("  theta=%d\n", theta);
      theta = (Theta + c->az)%360;
      for(p=0; p<cross_len[c->cross]; p++){
         // cylindrical coordinates r, theta, z
         // scale the template radius pt w/ the given radius
         // in our case a = 100 km, b_r = 1 (no eccentricity)
         r = cross[c->cross][p].r * c->a*c->b_n/sqrt(Sin(Theta)*Sin(Theta) + c->b_n*c->b_n*Cos(Theta)*Cos(Theta));

         // cyl to cartesian
         x_d = r*Cos(theta);
         y_d = r*Sin(theta);
         z_d = cross[c->cross][p].z * c->depth;

         x = (int)floor(x_d) + c->center->x;
         y = (int)floor(y_d) + c->center->y;
         z = (int)floor(z_d);

         if(abs(x) <= MAX_DIM && abs(y) <= MAX_DIM){
            pts[2*numPts] = x;
            pts[2*numPts + 1] = y;
            pts_z[numPts] = z;
            color[numPts] = get_random_in_range(r_color); 
            numPts++;
         }         
      }
   }
   pts[2*numPts] = c->center->x;
   pts[2*numPts + 1] = c->center->y;
   pts_z[numPts] = cross_ctr[c->cross]*c->depth;
   color[numPts] = get_random_in_range(r_color); 
   numPts++;
   numCraters++;
}

void draw_random_crater(){
   next_center.x = get_random_in_range(r_dim);//(MAX_DIM*2)*(float)rand()/(float)(RAND_MAX) - MAX_DIM;
   next_center.y = get_random_in_range(r_dim);//(MAX_DIM*2)*(float)rand()/(float)(RAND_MAX) - MAX_DIM;


   crater c;

   c.a = get_random_in_range(r_diameter); //r_diameter.min + (r_diameter.max - r_diameter.min)*(float)rand()/(float)(RAND_MAX);
   c.b_n = (float)get_random_in_range(r_eccentricity)/(float)r_eccentricity_lim; //(r_eccentricity.min + (r_eccentricity.max - r_eccentricity.min)*(float)rand()/(float)(RAND_MAX))/100.0; 
   c.depth = get_random_in_range(r_depth); //r_depth.min + (r_depth.max - r_depth.min)*(float)rand()/(float)(RAND_MAX); // i guess lol idk
   c.cross = rand()%r_num_cross;
   c.center = &next_center;
   c.az = rand()%360;
   printf("Crater Parameters: a=%d b_n=%f d=%d center=(%d, %d) az=%d\n\n", c.a, c.b_n, c.depth, c.center->x, c.center->y, c.az); 
   draw_crater(&c);


}

void draw_centered_crater(){
   next_center.x = 0;
   next_center.y = 0;
   crater c;

   c.a = get_random_in_range(r_diameter); //r_diameter.min + (r_diameter.max - r_diameter.min)*(float)rand()/(float)(RAND_MAX);
   c.b_n = (float)get_random_in_range(r_eccentricity)/(float)r_eccentricity_lim; //(r_eccentricity.min + (r_eccentricity.max - r_eccentricity.min)*(float)rand()/(float)(RAND_MAX))/100.0; 
   c.depth = get_random_in_range(r_depth); //r_depth.min + (r_depth.max - r_depth.min)*(float)rand()/(float)(RAND_MAX); // i guess lol idk
   c.cross = 1;
   c.center = &next_center;
   c.az = rand()%360;
   printf("Crater Parameters: a=%d b_n=%f d=%d center=(%d, %d) az=%d\n\n", c.a, c.b_n, c.depth, c.center->x, c.center->y, c.az); 
   draw_crater(&c);

}


void reset_canvas(){
   memset(pts,0,sizeof(pts));
   memset(pts_z,0,sizeof(pts_z));
   memset(color,0,sizeof(color));

   // add corners
      
   pts[0] = MAX_DIM;
   pts[1] = MAX_DIM;
   pts_z[0] = 0;

   pts[2] = -MAX_DIM;
   pts[3] = MAX_DIM;
   pts_z[1] = 0;

   pts[4] = -MAX_DIM;
   pts[5] = -MAX_DIM;
   pts_z[2] = 0;

   pts[6] = MAX_DIM;
   pts[7] = -MAX_DIM;
   pts_z[3] = 0;

   // initialize default colors
   color[0] = DEFAULT_GRAY;
   color[1] = DEFAULT_GRAY;
   color[2] = DEFAULT_GRAY;
   color[3] = DEFAULT_GRAY;

   //zero out crater counter
   numCraters = 0;
   // reset location
   next_center.x = 0;
   next_center.y = 0;

   numPts = 4; // always starts w/ 4 corners of canvas
}


void displayScene(){
   int i,v0,v1,v2;
   float x1,y1,z1,
         x2,y2,z2,
         x3,y3,z3,
         g1, g2, g3; 
   
   // find triangles
   int* vtx = BuildTriangleIndexList(
      (void*)pts,
      0, // integer points
      numPts,
      2, // 2 dims
      1, // ccw
      &num_vtx
   );
   float white[] = {1,1,1,1};
   float Emission[]  = {0.0,0.0,0.01*emission,1.0};
   glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,shiny);
   glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,white);
   glMaterialfv(GL_FRONT_AND_BACK,GL_EMISSION,Emission);
   glPushMatrix();
   glEnable(GL_TEXTURE_2D);
   glEnable(GL_DEPTH_TEST);
   glEnable(GL_CULL_FACE);
   glTexEnvi(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
   glColor3f(1,1,1);
   glBindTexture(GL_TEXTURE_2D,texture);
   
   for(i=0; i<num_vtx; i+=3){
      // draw triangle
      v0 = vtx[i];
      x1 = ((float)pts[2*v0])/MAX_DIM*CANVAS_LEN;
      y1 = ((float)pts[2*v0+1])/MAX_DIM*CANVAS_LEN;
      z1 = ((float)pts_z[v0])/MAX_DIM*CANVAS_LEN;

      v1 = vtx[i+1];
      x2 = ((float)pts[2*v1])/MAX_DIM*CANVAS_LEN;
      y2 = ((float)pts[2*v1+1])/MAX_DIM*CANVAS_LEN;
      z2 = ((float)pts_z[v1])/MAX_DIM*CANVAS_LEN;

      v2 = vtx[i+2];
      x3 = ((float)pts[2*v2])/MAX_DIM*CANVAS_LEN;
      y3 = ((float)pts[2*v2+1])/MAX_DIM*CANVAS_LEN;
      z3 = ((float)pts_z[v2])/MAX_DIM*CANVAS_LEN;

      // normal vector math
      //  Planar vector 0
      float dx0 = x1-x2;
      float dy0 = y1-y2;
      float dz0 = z1-z2;

      //  Planar vector 1
      float dx1 = x3-x1;
      float dy1 = y3-y1;
      float dz1 = z3-z1;

      //  Normal
      float Nx = dy0*dz1 - dy1*dz0;
      float Ny = dz0*dx1 - dz1*dx0;
      float Nz = dx0*dy1 - dx1*dy0;
      glBindTexture(GL_TEXTURE_2D,texture);
      glBegin(GL_TRIANGLES);
      glNormal3f(Nx,Nz,Ny);

      if(z1 == 0 && z2 == 0 && z3 == 0){
         g1 = g2 = g3 = (float)DEFAULT_GRAY/(float)r_color_lim;;
      } else {
         g1 = (float)color[v0]/(float)r_color_lim;
         g2 = (float)color[v1]/(float)r_color_lim;
         g3 = (float)color[v2]/(float)r_color_lim;

      }

      glColor3f(g1,g1,g1); 
      glTexCoord2f(((x1+CANVAS_LEN)/(float)(CANVAS_LEN*2)),((y1+CANVAS_LEN)/(float)(CANVAS_LEN*2)));
      glVertex3f(x1,z1,y1);

      glColor3f(g2,g2,g2);
      glTexCoord2f(((x2+CANVAS_LEN)/(float)(CANVAS_LEN*2)),((y2+CANVAS_LEN)/(float)(CANVAS_LEN*2)));
      glVertex3f(x2,z2,y2);

      glColor3f(g3,g3,g3);
      glTexCoord2f(((x3+CANVAS_LEN)/(float)(CANVAS_LEN*2)),((y3+CANVAS_LEN)/(float)(CANVAS_LEN*2)));
      glVertex3f(x3,z3,y3);

      glEnd();
      
      if(view_pts){
         glPointSize(10);
         glColor3f(.5,0,0);
         glBegin(GL_POINTS);

         glNormal3f(Nx,Nz,Ny);
         
         glVertex3f(x1,z1,y1);
         glVertex3f(x2,z2,y2);
         glVertex3f(x3,z3,y3);
         glEnd();
      }      
   }
   glDisable(GL_CULL_FACE);
   glDisable(GL_DEPTH_TEST);
   glDisable(GL_TEXTURE_2D);
   glPopMatrix();

   free(vtx); // no mem leaks >:(

}

void displayParams(){
   glColor3f(1,1,1);
   glWindowPos2i(5,5);
   Print("Viewport: Angle=%d,%d Dim=%.1f FOV=%d View=%s",th,ph,dim,fov,mode?"First Person":"Perspective");
   glWindowPos2i(5,45);
   Print("Sun Location: Altitude=%d Azimuth=%d",ylight,zh);
   glWindowPos2i(5,25);
   Print("NumCraters=%d",numCraters);
   //  Render the scene and make it visible
   ErrCheck("display");
}

void configureLighting(){
   // LIGHT AND SHADING
   glShadeModel(GL_SMOOTH);
   if (light) {
      //  Enable Lighting
      //  Translate intensity to color vectors
      float Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0};
      float Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0};
      float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
      //double angle = (minute*360.0)/1440.0 + omega;
      //  Light position
      float Position[]  = {distance*Sin(zh)*Cos(ylight) ,  distance*Sin(zh)*Sin(ylight),  distance*Cos(ylight),1.0};
      //  Draw light position as ball (still no lighting here)
      glColor3f(1,1,1);
      ball(Position[0],Position[2],Position[1] , 0.1);
      //  OpenGL should normalize normal vectors
      glEnable(GL_NORMALIZE);
      //  Enable lighting
      glEnable(GL_LIGHTING);
      //  Location of viewer for specular calculations
      glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,0);
      //  glColor sets ambient and diffuse color materials
      glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
      glEnable(GL_COLOR_MATERIAL);
      //  Enable light 0
      glEnable(GL_LIGHT0);
      //  Set ambient, diffuse, specular components and position of light 0
      glLightfv(GL_LIGHT0,GL_AMBIENT ,Ambient);
      glLightfv(GL_LIGHT0,GL_DIFFUSE ,Diffuse);
      glLightfv(GL_LIGHT0,GL_SPECULAR,Specular);
      glLightfv(GL_LIGHT0,GL_POSITION,Position);
   }
   else glDisable(GL_LIGHTING);
}

/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
   //  Erase the window and the depth buffer
   glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
   //  Enable Z-buffering in OpenGL
   glEnable(GL_DEPTH_TEST);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective - set eye position
   if (mode == 0)
   {
      double Ex = -2*dim*Sin(th)*Cos(ph);
      double Ey = +2*dim        *Sin(ph);
      double Ez = +2*dim*Cos(th)*Cos(ph);
      gluLookAt(Ex,Ey,Ez , 0,0,0 , 0,Cos(ph),0);
   }
   else { // mode == 1 - POV
      // Determine Eye Location
      // Decide where you are looking
      // Decide up (y-axis for now)
      double dx_fp = Cos(th);
      double dz_fp = Sin(th);
      double dy_fp = .1;
      double Cx = Ex_fp + dx_fp;
      double Cy = Ey_fp + dy_fp;
      double Cz = Ez_fp + dz_fp;

      gluLookAt(Ex_fp,Ey_fp,Ez_fp, Cx,Cy,Cz, 0,1,0);
   }
   

   configureLighting();
   
   displayScene(); //  Draw scene

   glDisable(GL_LIGHTING);

   displayParams(); //  Display parameters
   
   glFlush();
   glutSwapBuffers();
}

/*
 *  GLUT calls this routine when there is nothing else to do
 */
void idle()
{
   double t = glutGet(GLUT_ELAPSED_TIME)/1000.0;
   //minute = (minute+1)%1440;
   ylight = fmod(90*t,360.0); // moves sun TODO: make this do something less stupid, somehow changing alt/az like a real sun would
   glutPostRedisplay();
}

/*
 *  Set projection
 */
static void Project()
{
   //  Tell OpenGL we want to manipulate the projection matrix
   glMatrixMode(GL_PROJECTION);
   //  Undo previous transformations
   glLoadIdentity();
   //  Perspective transformation
   gluPerspective(fov,asp,dim/4,4*dim);
   //  Switch to manipulating the model matrix
   glMatrixMode(GL_MODELVIEW);
   //  Undo previous transformations
   glLoadIdentity();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width,int height)
{
   //  Ratio of the width to the height of the window
   asp = (height>0) ? (double)width/height : 1;
   //  Set the viewport to the entire window
   glViewport(0,0, RES*width,RES*height);
   //  Set projection
   Project();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key,int x,int y)
{
   if(mode==1) {
      // First Person 
      if (key == GLUT_KEY_RIGHT)
         th += 5;
      else if (key == GLUT_KEY_LEFT)
         th -= 5;
      else if (key == GLUT_KEY_UP){
         double dx_fp = Cos(th);
         double dz_fp = Sin(th);
         double dy_fp = 0;
         double dt = .1; // idk???
         Ex_fp += dt*dx_fp;
         Ey_fp += dt*dy_fp;
         Ez_fp += dt*dz_fp;
      }
      else if (key == GLUT_KEY_DOWN){
         double dx_fp = Cos(th);
         double dz_fp = Sin(th);
         double dy_fp = 0;
         double dt = .1; // idk???
         Ex_fp -= dt*dx_fp;
         Ey_fp -= dt*dy_fp;
         Ez_fp -= dt*dz_fp;
      }

   } else {
      // Perspective
      //  Right arrow key - increase angle by 5 degrees
      if (key == GLUT_KEY_RIGHT)
         th += 5;
      //  Left arrow key - decrease angle by 5 degrees
      else if (key == GLUT_KEY_LEFT)
         th -= 5;
      //  Up arrow key - increase elevation by 5 degrees
      else if (key == GLUT_KEY_UP)
         ph += 1;
      //  Down arrow key - decrease elevation by 5 degrees
      else if (key == GLUT_KEY_DOWN)
         ph -= 1;
   }
   // for any mode
   //  PageUp key - increase dim
   if (key == GLUT_KEY_PAGE_UP)
      dim += 0.1;
   //  PageDown key - decrease dim
   else if (key == GLUT_KEY_PAGE_DOWN && dim>1)
      dim -= 0.1;
   //  Keep angles to +/-360 degrees
   th = th%360;
   ph = ph%360;
   zh = zh%360;
   ylight = ylight%360; // altitude from 0 (on horizon) to 90 (overhead)
   //  Update projection
   Project();
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch,int x,int y)
{
   //  Exit on ESC
   if (ch == 27)
      exit(0);
   //  Reset view angle
   else if (ch == '0')
      th = ph = 0;
   //  Toggle light movement
   else if (ch == 'm')
      move = 1-move;
   //  Move light
   else if (ch == 'a')
      zh += 1;
   else if (ch == 'd')
      zh -= 1;
   //  Toggle lighting
   else if (ch == 'l' || ch == 'L')
      light = 1-light;
   //  Light elevation
   else if (ch=='s')
      ylight -= 1;
   else if (ch=='w')
      ylight += 1;
   //  Switch display mode
   else if (ch == 'v')
      mode = (mode+1)%2;
   else if (ch == 'V')
      mode = (mode-1)%2;
   //  Change field of view angle
   else if (ch == '-' && ch>1)
      fov--;
   else if (ch == '+' && ch<179)
      fov++;
   //  Add or Clear Craters
   else if (ch == 'c')
      draw_random_crater();  
   else if (ch == 'x')
      reset_canvas();
   else if (ch == '1') 
      draw_centered_crater();
   else if (ch == 'p')
      view_pts = (view_pts+1)%2;
   

   //  Keep angles to +/-360 degrees
   // viewport angles
   th = th%360;
   ph = ph%360;
   // light pos angles
   zh = zh%360;
   ylight = ylight%360;
   //  Reproject
   Project();
   //  Animate if requested
   glutIdleFunc(move?idle:NULL);
   //  Tell GLUT it is necessary to redisplay the scene
   glutPostRedisplay();
}

/*
 *  Start up GLUT and tell it what to do
 */
int main(int argc,char* argv[])
{
   //  Initialize GLUT
   glutInit(&argc,argv);
   //  Request double buffered, true color window with Z buffering at 600x600
   glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
   glutInitWindowSize(600,800);
   glutCreateWindow("Giselle Koo");
#ifdef USEGLEW
   //  Initialize GLEW
   if (glewInit()!=GLEW_OK) Fatal("Error initializing GLEW\n");
#endif
   //  Set callbacks
   glutIdleFunc(idle);
   glutDisplayFunc(display);
   glutReshapeFunc(reshape);
   glutSpecialFunc(special);
   glutKeyboardFunc(key);
   // import texture
   texture = LoadTexBMP("dirt.bmp");
   //  Pass control to GLUT so it can interact with the user
   reset_canvas();
   srand(time(NULL));
   ErrCheck("init");
   glutMainLoop();
   return 0;
}

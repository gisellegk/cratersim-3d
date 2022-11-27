#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
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
//  Default resolution
//  For Retina displays compile with -DRES=2
#ifndef RES
#define RES 1
#endif


#define Cos(th) cos(3.14159265/180*(th))
#define Sin(th) sin(3.14159265/180*(th))


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

// CRATER DEV HERE
// crater variables
int numCraters = 0;

struct Crater {
   // holds crater parameters
   // shape generation
   // ellipse (polar): r(theta) = a* b_n / sqrt(1-(1-b_n^2)cos(theta)^2)
   int a;     // scalar / length of semimajor axis "diameter" (0,50]?
   int b_n;   // normalized semiminor axis (b/a) "eccentricity" (0,1]
   int cross; // cross-sectional shape - later expand to chebyshev
   int text;  // texture index
   // location / orientation
   int x;   // origin x
   int y;   // origin y 
   int az;  // azimuthal orientation of semimajor axis
} crater;

void crater(){
   // Draw crater here
}

void placeCrater(){
   //place crater here
}


int mode=0;       //  Projection mode
int move=1;       //  Move light
int th=0;         //  Azimuth of view angle
int ph=15;        //  Elevation of view angle
int fov=55;       //  Field of view (for perspective)
double asp=1;     //  Aspect ratio
double dim=5.0;   //  Size of world

// for First Person
// eye position
double Ex_fp=-1;
double Ey_fp=.5;
double Ez_fp=1;

// Light values
int one       =   1;  // Unit value
int distance  =   500;  // Light distance
int inc       =  10;  // Ball increment
int emission  =   0;  // Emission intensity (%)
int ambient   =  10;  // Ambient intensity (%)
int diffuse   =  50;  // Diffuse intensity (%)
int specular  =   0;  // Specular intensity (%)
int shininess =   0;  // Shininess (power of two)
float shiny   =   1;  // Shininess (value)
int zh        =  90;  // Light azimuth
int ylight    =  45;  // Elevation of light
typedef struct {float x,y,z;} vtx;
typedef struct {int A,B,C;} tri;
#define n 500
vtx is[n];

//  Macro for sin & cos in degrees
#define Cos(th) cos(3.14159265/180*(th))
#define Sin(th) sin(3.14159265/180*(th))

void displayScene(){

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

void displayParams(){
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

   //  Enable Lighting
   //  Translate intensity to color vectors
   float Ambient[]   = {0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0};
   float Diffuse[]   = {0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0};
   float Specular[]  = {0.01*specular,0.01*specular,0.01*specular,1.0};
   //  Light position
   float Position[]  = {distance*Cos(zh)*Sin(ylight),distance*Cos(ylight),distance*Sin(zh)*Cos(ylight),1.0};
   //  Draw light position as ball (still no lighting here)
   glColor3f(1,1,1);
   ball(Position[0],Position[1],Position[2] , 0.1);
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
      double dy_fp = 0;
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
   zh = fmod(90*t,360.0); // moves sun TODO: make this do something less stupid, somehow changing alt/az like a real sun would
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
         ph += 5;
      //  Down arrow key - decrease elevation by 5 degrees
      else if (key == GLUT_KEY_DOWN)
         ph -= 5;
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
   ylight = ylight%90; // altitude from 0 (on horizon) to 90 (overhead)
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
   else if (ch == 'c'){
      // TODO: call add crater here
      numCraters++;
   }
   else if (ch == 'x'){
      // TODO: call reset here
      numCraters = 0;
   }
   //  Keep angles to +/-360 degrees
   // viewport angles
   th = th%360;
   ph = ph%360;
   // light pos angles
   zh = zh%360;
   ylight = ylight%90;
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
   glutInitWindowSize(600,600);
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
   //  Pass control to GLUT so it can interact with the user
   ErrCheck("init");
   glutMainLoop();
   return 0;
}

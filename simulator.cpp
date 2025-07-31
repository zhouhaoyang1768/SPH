///////////////////////////////////////////////////////////////////////////////
// This is an extensively reworked version of the "demo.c" file
// from Jos Stam's original "Stable Fluids" code:
//
// http://www.dgp.toronto.edu/people/stam/reality/Research/zip/CDROM_GDC03.zip
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdlib>
#include <cstdio>
#include <cstring>
#ifdef USING_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#include "SPH2.h"
SPH fluid;

int mouseDown[3];
int oldMouseX, oldMouseY, mouseX, mouseY;

int   windowX      = 1024;
int   windowY      = 1024;
bool  drawVelocity = false;
int   res          = 128;
float force        = 5.0;
float source       = 100.0;

bool animate = true;
bool step = false;

// the fluid simulation object
//FLUID_2D* fluid;

///////////////////////////////////////////////////////////////////////////////
// Keyboard command processing function
///////////////////////////////////////////////////////////////////////////////
void keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{

		case 'q':
		case 'Q':
			exit(0);
			break;

		case 'v':
		case 'V':
			drawVelocity = !drawVelocity;
			break;

    case 'a':
    case 'A':
      animate = !animate;
      break;
    case 'r':
    case 'R':
      fluid.reset();
      break;
    case ' ':
      step = true;
      animate = true;
      break;
	}
}

///////////////////////////////////////////////////////////////////////////////
// Mouse position and click processing function
///////////////////////////////////////////////////////////////////////////////
void mouseCallback(int button, int state, int x, int y)
{
	oldMouseX = mouseX = x;
	oldMouseX = mouseY = y;

	mouseDown[button] = state == GLUT_DOWN;

}

///////////////////////////////////////////////////////////////////////////////
// Mouse movement processing function
///////////////////////////////////////////////////////////////////////////////
void motionCallback(int x, int y)
{
	mouseX = x;
	mouseY = y;
}

///////////////////////////////////////////////////////////////////////////////
// Window shaping function
///////////////////////////////////////////////////////////////////////////////
void reshapeCallback(int width, int height)
{
	glutReshapeWindow(width, height);

	windowX = width;
	windowY = height;
}

///////////////////////////////////////////////////////////////////////////////
// Idle command processing function
///////////////////////////////////////////////////////////////////////////////
void idleCallback()
{
  if (animate)
  {
	// static int count = 0;
	// ++count;
	// if (count > 500) animate = false;
	if (mouseDown[0])
	{
		fluid.add((float)mouseX/windowX, ((float)(windowY - mouseY))/windowY);
	}
	if (mouseDown[2])
	{
		fluid.removeDirt((float)mouseX/windowX, ((float)(windowY - mouseY))/windowY);
	}
	for (int i = 0; i < 5; ++i)
		fluid.update();
  }
  if (step) { step = false; animate = false; }

  glutPostRedisplay();
}

void timedIdleCallBack(int)
{
  idleCallback();
  glutTimerFunc(1, &timedIdleCallBack, 0);
}

///////////////////////////////////////////////////////////////////////////////
// The drawing function
///////////////////////////////////////////////////////////////////////////////
void displayCallback()
{
	//glViewport(0, 0, windowX, windowY);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(234/256.0, 208/256.0, 168/256.0, 1.0f);
//	glClearColor(0, 0, 0, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);


    

  fluid.render();

	glutSwapBuffers();
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char ** argv)
{
	glutInit(&argc, argv);
	
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowSize(windowX, windowY);
	glutCreateWindow(argv[0]);

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers ();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers ();

	glutKeyboardFunc(keyboardCallback);
	glutMouseFunc(mouseCallback);
	glutMotionFunc(motionCallback);
	glutReshapeFunc(reshapeCallback);
	//glutIdleFunc(idleCallback);
  glutTimerFunc(1, &timedIdleCallBack, 0);
	glutDisplayFunc(displayCallback);
	glViewport(0, 0, windowX, windowY);
  
  //fluid = new FLUID_2D_BOUNDED(res, res, 0.1, mode);
  //fluid = new FLUID_2D_PERIODIC(res, res, 0.1);

	glutMainLoop ();

  return 0;
}

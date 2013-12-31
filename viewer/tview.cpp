// Rotating shape viewer
// Adapted from Edward Angel's Moving Viewer example for OpenGL

// We use the Lookat function in the display callback to point the viewer,
// whose position can be altered by the x, X, y, Y, z, and Z keys.

#include <stdlib.h>
#include <GL/glut.h>
#include <assert.h>
#include "plylib.h"
#include "modelseg.h" // includes "viewer.h" already

static GLfloat theta[] = {0.0, 0.0, 0.0}; // Amount rotated in x, y, z
static GLint axis = 2;
static GLdouble viewer[] = {0.0, 0.0, -15.0}; // initial viewer location
Model model; // The actual tooth model read in from a PLY file
int fvlist[MAX_NUMBER_OF_SEGMENTS][MAX_NUMBER_OF_TRIANGLES * 3 + 1];
int numSegments; // Actual number of segments in segments array
float minx, miny, minz, maxx, maxy, maxz; // Bounding box of model

void InitMain ()
{
	// Define light parameters
	float Light_ambient[] = {0.5, 0.5, 0.5, 1.0};
	float Light_diffuse[] = {0.5, 0.5, 0.5, 1.0};
	float Light_specular[] = {1.0, 1.0, 1.0, 1.0};
	float Light_position[] = {0.0, 0.0, 100.0, 1.0};

	// Set light parameters
	glLightfv (GL_LIGHT0, GL_AMBIENT, Light_ambient);
	glLightfv (GL_LIGHT0, GL_DIFFUSE, Light_diffuse);
	glLightfv (GL_LIGHT0, GL_SPECULAR, Light_specular);
	glLightfv (GL_LIGHT0, GL_POSITION, Light_position);

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glShadeModel (GL_SMOOTH);
	glPixelStorei (GL_UNPACK_ALIGNMENT, 1);

	glEnable (GL_DEPTH_TEST);

	// Uncomment these two lines to enable light
	//glEnable (GL_LIGHTING);
	//glEnable (GL_LIGHT0);

	// Load the model
	int i;
	char filename[] = "../data/Case1600_Tooth_19.ply";
	LoadModel (filename, &model);

	// Check model's fvlist for validity
	for (i = 0; i < model.nfaces; i++)
	{
		assert (0 <= model.fvlist[3 * i]);
		assert (model.fvlist[3 * i] < model.nverts);
		assert (0 <= model.fvlist[3 * i + 1]);
		assert (model.fvlist[3 * i + 1] < model.nverts);
		assert (0 <= model.fvlist[3 * i + 2]);
		assert (model.fvlist[3 * i + 2] < model.nverts);
	}

	// Compute rectangular prism containing model
	minx = model.vertices[0].coord[0];
	miny = model.vertices[0].coord[1];
	minz = model.vertices[0].coord[2];
	maxx = model.vertices[0].coord[0];
	maxy = model.vertices[0].coord[1];
	maxz = model.vertices[0].coord[2];
	for (i = 1; i < model.nverts; i++)
	{
		if (model.vertices[i].coord[0] < minx)
			minx = model.vertices[i].coord[0];
		if (model.vertices[i].coord[1] < miny)
			miny = model.vertices[i].coord[1];
		if (model.vertices[i].coord[2] < minz)
			minz = model.vertices[i].coord[2];
		if (model.vertices[i].coord[0] > maxx)
			maxx = model.vertices[i].coord[0];
		if (model.vertices[i].coord[1] > maxy)
			maxy = model.vertices[i].coord[1];
		if (model.vertices[i].coord[2] > maxz)
			maxz = model.vertices[i].coord[2];
	}

	// Compute center of the bounding prism
	float cx = 0.5 * (minx + maxx);
	float cy = 0.5 * (miny + maxy);
	float cz = 0.5 * (minz + maxz);

	// Translate model so prism center becomes origin
	for (i = 0; i < model.nverts; i++)
	{
		model.vertices[i].coord[0] -= cx;
		model.vertices[i].coord[1] -= cy;
		model.vertices[i].coord[2] -= cz;
		model.vertices[i].coord[0] *= 1000;
		model.vertices[i].coord[1] *= 1000;
		model.vertices[i].coord[2] *= 1000;
	}
	minx -= cx; miny -= cy; minz -= cz;
	maxx -= cx; maxy -= cy; maxz -= cz;
	minx *= 1000;
	miny *= 1000;
	minz *= 1000;
	maxx *= 1000;
	maxy *= 1000;
	maxz *= 1000;
	cx *= 1000;
	cy *= 1000;
	cz *= 1000;

	// Print prism information to standard output
	printf ("Box containing model:\n");
	printf ("Min = (%f, %f, %f)\n", minx, miny, minz);
	printf ("Max = (%f, %f, %f)\n", maxx, maxy, maxz);
	printf ("Center = (%f, %f, %f)\n\n", cx, cy, cz);

	// Segment model
	char segfilename[256] = "../jdata/gwh99.out";
	ModelSegmenter modseg (model);
	modseg.ReadSegmentation (segfilename, fvlist, numSegments);

	// Print segmentation information to standard output
	printf ("Segmentation:\n");
	printf ("Number of segments = %d\n", numSegments);
}

void DrawBoundingBox ()
{
	glBegin (GL_LINE_LOOP);
	glVertex3f (minx, miny, minz);
	glVertex3f (maxx, miny, minz);
	glVertex3f (maxx, maxy, minz);
	glVertex3f (minx, maxy, minz);
	glEnd ();
	glBegin (GL_LINE_LOOP);
	glVertex3f (minx, miny, maxz);
	glVertex3f (maxx, miny, maxz);
	glVertex3f (maxx, maxy, maxz);
	glVertex3f (minx, maxy, maxz);
	glEnd ();
	glBegin (GL_LINE_LOOP);
	glVertex3f (minx, miny, minz);
	glVertex3f (minx, maxy, minz);
	glVertex3f (minx, maxy, maxz);
	glVertex3f (minx, miny, maxz);
	glEnd ();
	glBegin (GL_LINE_LOOP);
	glVertex3f (maxx, miny, minz);
	glVertex3f (maxx, maxy, minz);
	glVertex3f (maxx, maxy, maxz);
	glVertex3f (maxx, miny, maxz);
	glEnd ();
	glBegin (GL_LINE_LOOP);
	glVertex3f (minx, miny, minz);
	glVertex3f (maxx, miny, minz);
	glVertex3f (maxx, miny, maxz);
	glVertex3f (minx, miny, maxz);
	glEnd ();
	glBegin (GL_LINE_LOOP);
	glVertex3f (minx, maxy, minz);
	glVertex3f (maxx, maxy, minz);
	glVertex3f (maxx, maxy, maxz);
	glVertex3f (minx, maxy, maxz);
	glEnd ();
}

void DrawSegments ()
{
	// Draw each segment using model, fvlist, and numSegments
	// fvlist[i][] is the fvlist for segment #i
	// Should do basically the same thing as DrawModel in viewer.cpp,
	//   for each segment, but color each segment a bit differently

	if (model.has_normals)
		glInterleavedArrays (GL_N3F_V3F, 0, model.vertices);
	else
	{
		glDisableClientState (GL_NORMAL_ARRAY);
		glEnableClientState (GL_VERTEX_ARRAY);
		glVertexPointer (3, GL_FLOAT, 6 * sizeof (float),
		                 (char *)model.vertices + 3 * sizeof (float));
	}

	int i;
	for (i = 0; i < numSegments; i++)
	{
		// Draw ith segment
		glColor4f ((float)(i+1)/numSegments, 0.5, 0.5, 1.0);
		/*if (i == 1)
			glColor4f (0, 0, 1.0, 1.0);*/
		glDrawElements (GL_TRIANGLES,
		                fvlist[i][MAX_NUMBER_OF_TRIANGLES * 3] * 3,
		                GL_UNSIGNED_INT,
		                fvlist[i]);
	}
}

void display ()
{
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Update viewer position in modelview matrix
	glLoadIdentity ();
	gluLookAt (viewer[0], viewer[1], viewer[2], 0.0, 0.0, 0.0,
	                                            0.0, 1.0, 0.0);

	// Rotate model
	glRotatef (theta[0], 1.0, 0.0, 0.0);
	glRotatef (theta[1], 0.0, 1.0, 0.0);
	glRotatef (theta[2], 0.0, 0.0, 1.0);

	// Draw model and its bounding box in normal coordinates
	glColor4f (1.0, 1.0, 1.0, 1.0);
	DrawBoundingBox ();
	DrawSegments ();
	//DrawModel (&model);

	glFlush ();
	glutSwapBuffers ();
}

void mouse (int btn, int state, int x, int y)
{
	if (btn == GLUT_LEFT_BUTTON && state == GLUT_DOWN) axis = 0;
	if (btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) axis = 1;
	if (btn == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) axis = 2;
	theta[axis] += 2.0;
	if (theta[axis] > 360.0) theta[axis] -= 360.0;
	display ();
}

void keys (unsigned char key, int x, int y)
{
	// Use x, X, y, Y, z, and Z keys to move viewer

	if (key == 'x') viewer[0] -= 1.0;
	if (key == 'X') viewer[0] += 1.0;
	if (key == 'y') viewer[1] -= 1.0;
	if (key == 'Y') viewer[1] += 1.0;
	if (key == 'z') viewer[2] -= 1.0;
	if (key == 'Z') viewer[2] += 1.0;
	display ();
}

void myReshape (int w, int h)
{
	glViewport (0, 0, w, h);

	// Use a perspective view
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	if (w <= h) glFrustum (-2.0, 2.0, -2.0 * (GLfloat)h / (GLfloat)w,
	                       2.0 * (GLfloat)h / (GLfloat)w, 2.0, 20.0);
	else glFrustum (-2.0, 2.0, -2.0 * (GLfloat)w / (GLfloat) h,
	                2.0 * (GLfloat)w / (GLfloat)h, 2.0, 20.0);

	// Or we can use gluPerspective
	// gluPerspective (45.0, w/h, 1.0, 10.0);

	glMatrixMode (GL_MODELVIEW);
}

int main (int argc, char **argv)
{
	glutInit (&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_ALPHA | GLUT_DEPTH);
	glutInitWindowSize (500, 500);
	glutCreateWindow ("Tooth Segmentation Viewer");

	InitMain ();
	glutReshapeFunc (myReshape);
	glutDisplayFunc (display);
	glutMouseFunc (mouse);
	glutKeyboardFunc (keys);

	glutMainLoop ();

	return 0;
}





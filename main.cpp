#include <GL/glew.h>
#include <glut/glut.h>
#include <stdio.h>
#include <stdlib.h>

#include "treerenderer.h"
#include "treegenerator.h"

tree_model *tm;

GLdouble rotation_x;
GLdouble rotation_y;
GLdouble position_z;
GLdouble scale;

int      click_button;
GLdouble click_rotation_x;
GLdouble click_rotation_y;
GLdouble click_scale;                                                          
GLdouble click_nx;
GLdouble click_ny;


int init()
{

	rotation_x = 0.0;
    rotation_y = 0.0;
    position_z = 5.0;
    scale = 1.0;  

	//glRotatef(270.f,0.f,1.f,0.f);
	init_tree(tm);
	glEnable(GL_DEPTH_TEST);

	glEnable(GL_NORMALIZE);                                                    
    glEnable(GL_DEPTH_TEST);
   // glEnable(GL_LIGHTING);
   // glEnable(GL_LIGHT0);

	return 1;
}

void draw()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslated(0.0, 0.0, -position_z);
    glRotated(rotation_x, 1.0, 0.0, 0.0);
    glRotated(rotation_y, 0.0, 1.0, 0.0);
    glScaled(scale, scale, scale); 


	//glMatrixMode(GL_MODELVIEW);
    //glLoadIdentity();
   // glRotatef(270.f,0.f,0.f,1.f);

	render_tree(tm);

	glutSwapBuffers();
}

void reshape(int w, int h)
{
    GLdouble x = 0.5 * (GLdouble) w / (GLdouble) h;
    GLdouble y = 0.5;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-x, x, -y, y, 1.0, 10.0);

    glViewport(0, 0, w, h);
}


void motion(int x, int y)
{
    GLdouble nx = (GLdouble) x / glutGet(GLUT_WINDOW_WIDTH);
    GLdouble ny = (GLdouble) y / glutGet(GLUT_WINDOW_HEIGHT);

    GLdouble dx = nx - click_nx;
    GLdouble dy = ny - click_ny;

    if (click_button == GLUT_LEFT_BUTTON)
    {
        rotation_x = click_rotation_x +  90.0 * dy;
        rotation_y = click_rotation_y + 180.0 * dx;

        if (rotation_x >   90.0) rotation_x  =  90.0;
        if (rotation_x <  -90.0) rotation_x  = -90.0;
        if (rotation_y >  180.0) rotation_y -= 360.0;
        if (rotation_y < -180.0) rotation_y += 360.0;
    }
    if (click_button == GLUT_RIGHT_BUTTON)
    {
        scale = click_scale - dy;                                              \
    }

    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    click_nx = (GLdouble) x / glutGet(GLUT_WINDOW_WIDTH);
    click_ny = (GLdouble) y / glutGet(GLUT_WINDOW_HEIGHT);

    click_button     = button;
    click_rotation_x = rotation_x;
    click_rotation_y = rotation_y;
    click_scale      = scale;
}

int main(int argc, char *argv[])
{
    int num_steps = 0;
    if(argc > 1) num_steps = atoi(argv[1]);
    else num_steps = 2;
	TreeGenerator tg;
    //tg.init_shadow_grid();
    //tg.propogate_shadow(50, 50, 50);
	tg.grow_tree(num_steps);
	tg.print_tree();


	tm = get_tree_model(tg.branches);

	

    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize(640, 480);
    glutInit(&argc, argv);

    glutCreateWindow(argv[0]);

    glutReshapeFunc(reshape);
    glutDisplayFunc(draw);
    glutMotionFunc(motion);
    glutMouseFunc(mouse);

    if (glewInit() == GLEW_OK)
    {
        init();
        glutMainLoop();
    }
    return 0;
}
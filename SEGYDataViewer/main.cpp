#include <GL/glut.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include "segy_reader.h"

int num_batch;

int num_x;
int num_t;

int timer = 0;

int data_offset = 0;

float scalar = 0;

float * data = NULL;

void drawStuff(void)
{
    float val;
    glPointSize(20);
    glBegin(GL_POINTS);
    for ( int k(0)
        , x(0)
        ; x < num_batch
        ; ++x
        )
    for ( int t(0)
        ; t < num_t
        ; ++t
        , ++k
        )
    {
        val = 0.5f + data[k+data_offset] * scalar;
        glColor3f
        ( val
        , val
        , val
        );
        glVertex3f 
        ( -1 + 2.0f * x / num_batch
        , 0
        , -1 + 2.0f * t / num_t
        );
    }
    glEnd();
    timer++;
    if ( timer > 10 )
    {
        timer = 0;
        data_offset += num_t * num_batch;
        if ( data_offset >= num_t * num_x )
        {
            data_offset = 0;
        }
    }
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff();
  glutSwapBuffers();
}

void init(void)
{
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
    /* aspect ratio */ 1.0,
    /* Z near */ 1.0, /* Z far */ 1000.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, 4.75,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, +2.0);
  glRotatef(90, 1.0, 0.0, 0.0);
  //glRotatef(-20, 0.0, 0.0, 1.0);
}

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        default : break ;
    }
}

void idle(void)
{
  glutPostRedisplay();
}

float find_abs_max ( float const * data , int n )
{
    float ret = 0;
    float temp;
    for ( int k(0)
        ; k < n
        ; ++k
        )
    {
        temp = fabs ( data[k] );
        if ( temp > ret )
        {
            ret = temp;
        }
    }
    return ret;
}

int main( int argc , char ** argv )
{
    // if ( argc == 2 )
    // {
    //     SEPReader reader ( argv[1] );
    //     num_t = reader . n1;
    //     num_batch = reader . n2;
    //     num_x = reader . n2 
    //           * reader . n3
    //           * reader . n4
    //           * reader . n5
    //           * reader . n6
    //           * reader . n7
    //           * reader . n8
    //           ;
    //     std::cout << num_x << std::endl;
    //     data = new float [ num_x * num_t ];
    //     memset ( &data[0] , 0 , num_x * num_t );
    //     reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );
    //     scalar = 0.45f / find_abs_max ( & data [ 0 ] , num_x * num_t );
    //     std::cout << scalar << std::endl;
    // }
    // else
    // if ( argc > 2 )
    // {
    //     SEPReader reader ( argv[1] , false );
    //     num_t = reader . n1;
    //     num_batch = reader . n2;
    //     num_x = reader . n2 
    //           * reader . n3
    //           * reader . n4
    //           * reader . n5
    //           * reader . n6
    //           * reader . n7
    //           * reader . n8
    //           ;
    //     std::cout << num_x << std::endl;
    //     data = new float [ num_x * num_t ];
    //     memset ( &data[0] , 0 , num_x * num_t );
    //     for ( int k(2)
    //         ; k < argc
    //         ; ++k
    //         )
    //     {
    //         std::cout << "k=" << k-1 << "     ";
    //         reader . OpenDataFile ( argv[k] );
    //         reader . read_sepval ( & data [ (k-1) * num_t * num_batch ] , reader . o1 , reader . o2 , reader . o3 , num_batch * num_t );
    //     }
    //     scalar = 0.45f / find_abs_max ( & data [ 0 ] , num_x * num_t );
    //     std::cout << scalar << std::endl;
    // }
    // else
    // {
    //     std::cout << "try ./a.out sep_file.h" << std::endl;
    //     return 0;
    // }

    if ( argc != 2 )
    {
        std::cout << "try ./a.out sep_file.h" << std::endl;
        return 0;
    }

    SEGYReader reader ( argv[1] );

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    init();
    glutMainLoop();

    return 0;
}


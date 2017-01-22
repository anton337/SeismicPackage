#include "viewer.h"

float * data = NULL;

float * wave_data = NULL;

float * reverse_data = NULL;

float timer = 0;

int num_x = 0;
int num_z = 0;

float scalar = 1;

float reverse_scalar = 1;

float wave_scalar = 1;

float scalar_1 = 1;
float scalar_2 = 1;
float scalar_3 = 1;

void set_timer ( float _timer ) { timer = _timer; }
void set_num_x ( int _num_x ) { num_x = _num_x; }
void set_num_z ( int _num_z ) { num_z = _num_z; }
void set_data_ptr ( float * _data ) { data = _data; }
void set_wave_ptr ( float * _wave_data ) { wave_data = _wave_data; }
void set_reverse_ptr ( float * _reverse_data ) { reverse_data = _reverse_data; }
void set_scalar ( float _scalar ) { scalar = _scalar; }
void set_wave_scalar ( float _wave_scalar ) { wave_scalar = _wave_scalar; }
void set_reverse_scalar ( float _reverse_scalar ) { reverse_scalar = _reverse_scalar; }

void drawString (void * font, char const *s, float x, float y, float z)
{
     unsigned int i;
     glRasterPos3f(x, y, z);
     for (i = 0; i < strlen (s); i++)
     {
         glutBitmapCharacter (font, s[i]);
     }
}

void drawStuff(void)
{
    if ( timer > 0 )
    {
        glColor3f ( 1 , 1 , 1 );
        std::stringstream str;
        str << std::showpoint << std::setprecision(6) << timer << "ms" << std::endl;
        drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), -0.275, 0.1, -0.275);
    }
    if ( data && wave_data && reverse_data )
    {
        float TH, PHI;
        glPointSize(7);
        glBegin(GL_POINTS);
        for ( int k(0)
            , x(0)
            ; x < num_x
            ; ++x
            )
        {
            for ( int z(0)
                ; z < num_z
                ; ++z
                , ++k
                )
            {
                glColor3f
                ( 0.5 + scalar_1 *         scalar *         data [k] 
                , 0.5 + scalar_2 * reverse_scalar * reverse_data [k] 
                , 0.5 + scalar_3 *    wave_scalar *    wave_data [k] 
                );
                TH = 2*M_PI*x/num_x;
                PHI = M_PI*z/num_z;
                glVertex3f 
                ( sinf(TH)*sinf(PHI)
                , -cosf(TH)*sinf(PHI)
                , cosf(PHI)
                );
            }
        }
        glEnd();

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
    /* Z near */ 0.1, /* Z far */ 1000.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, 1.2* /*0.6125*/4.75,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, +2.0);
  glRotatef(90, 1.0, 0.0, 0.0);
  //glRotatef(-20, 0.0, 0.0, 1.0);
}

void idle(void)
{
  glutPostRedisplay();
}

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        case 'a': scalar_1 /= 1.1; break;
        case 'q': scalar_1 *= 1.1; break;
        case 's': scalar_2 /= 1.1; break;
        case 'w': scalar_2 *= 1.1; break;
        case 'd': scalar_3 /= 1.1; break;
        case 'e': scalar_3 *= 1.1; break;
        default : break ;
    }
}

void worker ( int argc , char ** argv )
{

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    init();
    glutMainLoop();


}

boost::thread       * show ( int argc , char ** argv )
{

    return new boost::thread ( worker , argc , argv );

}


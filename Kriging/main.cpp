
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <math.h>

int width ;
int height;

float * data = new float [ width * height ];

float scalar = .01;

void drawStuff(void)
{
    glPointSize(2);
    glBegin(GL_POINTS);
    float val;
    for ( int k(0)
        , x(0)
        ; x < width
        ; ++x
        )
    for ( int t(0)
        ; t < height
        ; ++t
        , ++k
        )
    {
        val = 0.5f + data[k] * scalar;
        glColor3f
        ( val
        , val
        , val
        );
        glVertex3f 
        ( -1 + 2.0f * x / width
        , 0
        , -1 + 2.0f * t / height
        );
    }
    glEnd();
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
        case 'a': scalar *= 1.1 ; break ;
        case 'q': scalar /= 1.1 ; break ;
        default : break ;
    }
}

void idle(void)
{
  glutPostRedisplay();
}

struct selection_window
{

    std::vector < std::pair < float , float > > win;

    selection_window ( std::vector < float > center 
                     , float width 
                     )
    {
        for ( int k(0)
            ; k < center . size ()
            ; ++k
            )
        {
            win . push_back ( std::pair < float , float > ( center[k]-width , center[k]+width ) );
        }
    }

    selection_window ()
    {

    }

};

struct point
{

    std::vector < float > pos;

    std::vector < float > val;

    point ( std::vector < float > const & _pos
          , std::vector < float > const & _val
          )
    : pos ( _pos )
    , val ( _val )
    {

    }

};

struct collection 
{

    std::vector < point > data;

    void insert ( point const & p )
    {
        data . push_back ( p );
    }

};

struct point_bin
{

    selection_window win;

    collection data;

    point_bin ( std::vector < float > const & center
              , float width
              )
    : win ( center 
          , width
          )
    {

    }

    point_bin ()
    {

    }

    void insert ( point const & p )
    {
        data . insert ( p );
    }

};

inline float sqr ( float const & x )
{
    return x*x;
}

void interpolate ( point_bin const & bin
                 , std::vector < float > const & pos
                 , std::vector < float >       & val
                 )
{
    float sigma_2 = 50;
    for ( std::size_t j(0)
        ; j < bin . data . data . size ()
        ; ++j
        )
    {
        for ( std::size_t k(0)
            ; k < val . size ()
            ; ++k
            )
        {
            float dist = 0;
            for ( std::size_t x(0)
                ; x < pos . size ()
                ; ++x
                )
            {
                dist += sqr ( pos[x] - bin . data . data[j] . pos[x] );
            }
            val[k] += bin . data . data[j] . val[0] * exp ( - dist / sigma_2 );
        }
    }
}



int main( int argc , char ** argv )
{

    int num_pts = 1000000;

    std::vector < point > pt;

    for ( int k(0)
        ; k < num_pts
        ; ++k
        )
    {
        std::vector < float > pos ( 2 );
        pos[0] = rand()%1000;
        pos[1] = rand()%1000;
        std::vector < float > val ( 1 );
        val[0] = sin(pos[0]/50)*cos(pos[1]/50)*sin((pos[0]+pos[1])/200);
        pt . push_back ( point ( pos , val ) );
    }

    int num_x = 200;
    int num_y = 200;
    width = 1000;
    height = 1000;
    int scal_x = width / num_x;
    int scal_y = height / num_y;
    std::vector < std::vector < point_bin > > win ( num_x , std::vector < point_bin > ( num_y ) );
std::cout << num_x << " " << num_y << std::endl;
std::cout << scal_x << " " << scal_y << std::endl;
std::cout << "p1" << std::endl;
    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int y(0)
            ; y < num_y
            ; ++y
            )
        {
            std::vector < float > cen (2);
            cen[0] = x*scal_x;
            cen[1] = y*scal_y;
            win[x][y] = point_bin ( cen , scal_x );
        }
    }
std::cout << "p2" << std::endl;

    for ( int k(0)
        ; k < pt . size ()
        ; ++k
        )
    {
        int ind_x = pt[k].pos[0]/scal_x;
        int ind_y = pt[k].pos[1]/scal_y;
        for ( int dx=-1; dx <= 1; ++dx )
        for ( int dy=-1; dy <= 1; ++dy )
        if ( ind_x+dx >= 0 && ind_x+dx < num_x )
        if ( ind_y+dy >= 0 && ind_y+dy < num_x )
        win[ind_x+dx][ind_y+dy] . insert ( pt[k] );
    }
std::cout << "p3" << std::endl;

    std::vector < point > interp_pt;
    data = new float [ width * height ];

    for ( int k(0)
        , x(0)
        ; x < width
        ; ++x
        )
    {
        for ( int y(0)
            ; y < height
            ; ++y
            , ++k
            )
        {
            std::vector < float > pos ( 2 );
            pos[0] = x;
            pos[1] = y;
            std::vector < float > val ( 1 );
            interpolate ( win[(int)(x/scal_x)][(int)(y/scal_y)] 
                        , pos 
                        , val 
                        );
            interp_pt . push_back ( point ( pos , val ) );
            data[k] = val[0];
        }
    }
std::cout << "p4" << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width,height);
    glutCreateWindow(argv[1]);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    init();
    glutMainLoop();

    return 0;

}


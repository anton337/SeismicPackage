#include <GL/glut.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"
#include <complex>

int width  = 1000;
int height = 1000;

int num_th = 3000;
int num_r  = 1500 ;

float * transfer_data = NULL;

float scalar = 1;

void drawStuff(void)
{
  glLineWidth(1);
  float val;
  float x_pos;
  float y_pos;
  float R;
  float TH;
  glColor3f(1,0,0);
  glBegin(GL_LINES);
  glVertex3f(-1,0,0);
  glVertex3f( 1,0,0);
  glVertex3f(0,0,-1);
  glVertex3f(0,0, 1);
  for ( int th(0); th-1 < 120; ++th )
  {
  TH = 2*M_PI*(float)th/120;
  R = 1;
  x_pos = R * cos ( TH );
  y_pos = R * sin ( TH );
  if ( th > 0 )
  glVertex3f(x_pos/2,0,y_pos/2);
  glVertex3f(x_pos/2,0,y_pos/2);
  }
  glEnd();
  glBegin(GL_POINTS);
  for ( int th(0), k(0); th < num_th; ++th )
  for ( int r(0); r < num_r ; ++r, ++k )
  {
  TH = 2*M_PI*(float)th/num_th;
  R = 2*(float)r/num_r;
  val = scalar * transfer_data[k];
  glColor3f(val,val,val);
  x_pos = R * cos ( TH );
  y_pos = R * sin ( TH );
  glVertex3f(x_pos/2,0,y_pos/2);
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

void idle(void)
{
  glutPostRedisplay();
}

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 'a': scalar *= 1.1 ; break ;
        case 's': scalar /= 1.1 ; break ;
        case 27 : exit(0) ; break ;
        default : break ;
    }
}

float get_argument ( int &   arg 
                   , int     argc
                   , char ** argv
                   )
{
    ++arg;
    if ( arg > argc )
    {
        std::cout << "parameter error" << std::endl;
        exit(1);
    }
    return atof ( argv[arg] );
}

template < typename T >
inline void mem_cpy ( T * out , T * in , int num )
{
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        out[k] = in[k];
    }
}

template < typename T >
inline void mem_cpy_r2c ( fftwf_complex * out , T * in , int num )
{
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        out[k][0] = in[k];
        out[k][1] = 0;
    }
}

std::complex < float > butterworth_transfer ( std::complex < float > z )
{
    std::complex < float > c1 ( 1 , 0 );
    std::complex < float > c2 ( 2 , 0 );
    return c1/(c1+c2*z+c2*z*z+z*z*z);
}

void apply_butterworth_transfer ( fftwf_complex * data , int num_t )
{
    for ( int k(0) ; k < num_t ; ++k )
    {
        std::complex < float > out ( data[k][0] , data[k][1] );
        float arg ( 2 * M_PI * k / (float)num_t );
        std::complex < float > z ( cos ( arg ) , - sin ( arg ) );
        out *= butterworth_transfer ( z );
        data[k][0] = std::real ( out );
        data[k][1] = std::imag ( out );
    }
}

enum filterType { BUTTERWORTH = 1
                };

int main( int argc , char ** argv )
{

    if (argc == 1)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string file_name( argv[arg] );

    arg++;

    std::string output_file_name( argv[arg] );

    SEPReader reader ( file_name . c_str () );
    int num_t = reader . n1;
    int num_x = reader . n2 
              * reader . n3
              * reader . n4
              * reader . n5
              * reader . n6
              * reader . n7
              * reader . n8
              ;
    std::cout << num_t << std::endl;
    std::cout << num_x << std::endl;
    if ( num_t % 2 != 0 )
    {
        std::cout << "num_t should be even" << std::endl;
        exit(1);
    }
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    int num_f = num_t;
    fftwf_complex * fdata = new fftwf_complex [ num_x * num_f ];

    fftwf_complex * tmp  = new fftwf_complex [ num_t ];
    fftwf_complex * ftmp = new fftwf_complex [ num_f ];

    fftwf_plan plan ( fftwf_plan_dft_1d ( num_t
                                        , & tmp[0]
                                        , &ftmp[0]
                                        , FFTW_FORWARD
                                        , FFTW_ESTIMATE
                                        )
                    );

    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    {
        mem_cpy_r2c ( & tmp[0] , & data[x*num_t] , num_t );
        fftwf_execute ( plan );
        mem_cpy < float > ( (float*)(&fdata[x*num_f]) , (float*)(&ftmp[0]) , 2*num_f );
    }

    filterType filter_type = BUTTERWORTH;

    transfer_data = new float [ num_th * num_r ];

    float R;
    float TH;
    for ( int th(0), k(0); th < num_th; ++th )
    for ( int r(0); r < num_r ; ++r, ++k )
    {
    TH = 2*M_PI*(float)th/num_th;
    R = 2*(float)r/num_r;
    std::complex < float > c ( R * cos ( TH ) , R * sin ( TH ) );
    transfer_data[k] = log ( std::abs ( butterworth_transfer ( c ) ) ) ;
    }
    float max_norm ( 0 );
    for ( int k(0) ; k < num_th * num_r ; ++k )
    {
        if ( transfer_data[k] > max_norm )
        {
            max_norm = transfer_data[k];
        }
    }
    for ( int k(0) ; k < num_th * num_r ; ++k )
    {
        transfer_data[k] /= 1e-2*max_norm;
    }

    switch ( filter_type )
    {

        case BUTTERWORTH:
        {
            for ( int x(0)
                ; x < num_x
                ; ++x
                )
            {
                apply_butterworth_transfer ( &fdata[x*num_f] , num_t );
            }
            //mem_cpy ( & tmp[0] , &filter[0] , num_t );
            //fftwf_execute ( plan );
            break;
        }

        default:
        {
            std::cout << "undefined filter" << std::endl;
            exit(1);
        }

    }

    fftwf_destroy_plan ( plan );

    // SEPWriter writer ( output_file_name . c_str () 
    //                  , reader . o1 , reader . d1 , 2 * num_f
    //                  , reader . o2 , reader . d2 , reader . n2
    //                  , reader . o3 , reader . d3 , reader . n3
    //                  , reader . get_header_labels ()
    //                  , reader . get_sort_order ()
    //                  , (output_file_name + std::string("@")) . c_str()
    //                  );

    // writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    // writer . write_sepval ( (float*)fdata , reader . o1 , reader . o2 , reader . o3 , num_x * (2 * num_f) );

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


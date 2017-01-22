#include <GL/glut.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include <algorithm>
#include "sep_writer.h"
#include "sep_reader.h"

void drawStuff(void);

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

float scalar_1 = 1;
float scalar_2 = 1;
float scalar_3 = 1;

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        case 'd': scalar_1 *= 2 ; break ;
        case 'a': scalar_1 /= 2 ; break ;
        case 'w': scalar_2 *= 2 ; break ;
        case 's': scalar_2 /= 2 ; break ;
        case 'q': scalar_3 *= 2 ; break ;
        case 'z': scalar_3 /= 2 ; break ;
        default : break ;
    }
    std::cout << scalar_1 << " " << scalar_2 << " " << scalar_3 << std::endl;
}

void idle(void)
{
  glutPostRedisplay();
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

void propagate ( float   h       // dt
               , int     num_x 
               , int     num_t
               , float   scale_x
               , float   scale_t
               , float * p2 
               , float * p1 
               , float * c 
               , float * vel )
{
    float fact ( h * h / ( scale_x * scale_x ) );
    for ( int t(1)
        ; t+1 < num_t
        ; ++t
        )
    {
        for ( int x(1)
            ; x+1 < num_x
            ; ++x
            )
        {
            if ( fabs ( p2[t+x*num_t] ) > 1e-40 
              || fabs ( p1[t+1+x*num_t] ) > 1e-40 
              || fabs ( p1[t-1+x*num_t] ) > 1e-40 
              || fabs ( p1[t+(x+1)*num_t] ) > 1e-40 
              || fabs ( p1[t+(x-1)*num_t] ) > 1e-40 
               )
            {
                c[t+x*num_t] = fact * vel[t+x*num_t] * vel[t+x*num_t] * ( -4*p1[t+x*num_t] 
                                                                        +    p1[t+(x+1)*num_t]
                                                                        +    p1[t+(x-1)*num_t]
                                                                        +    p1[(t+1)+x*num_t]
                                                                        +    p1[(t-1)+x*num_t]
                                                                        )
                             + 2 * p1[t+x*num_t]
                             -     p2[t+x*num_t]
                             ;
            }
        }
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p2[k] = p1[k];
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p1[k] = c[k];
    }
}

inline float pow_2 ( float const & x )
{
    return x*x;
}

float find_max ( float const * x , int num )
{
    float m ( 0 ) , t;
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        t = fabs(x[k]);
        if ( t > m )
        {
            m = t;
        }
    }
    return m;
}

    float h = 0.00040;

    float time_frac = h / 0.004;

    float scale_x = 4;

    int num_x;
    int num_t;
    int num_t2;

    int s_num_t;
    int s_num_b;

    float * vel     ;

    float * s_data  ;

    float * data_p2 ; 

    float * data_p1 ;

    float * data_c  ;

    float * rdata_p2;

    float * rdata_p1;

    float * rdata_c ;

    float * odata   ;

struct recorder
{

    int ind;

    std::vector < float > trace;

};

std::vector < recorder > recorders;

void drawStuff(void)
{
    static int k(0);
    for ( int _k(0) ; _k < 5 ; ++_k , ++k )
    {
        //if ( k*time_frac > 50 )
        {
            if ( k*time_frac < 1.0*s_num_t )
            {
                for ( int t(0) ; t < recorders.size() ; ++t )
                {
                    rdata_p2[recorders[t].ind] = (k>0)?recorders[t].trace[k-1]:0;
                    rdata_p1[recorders[t].ind] = recorders[t].trace[k];
                    //data_p1[(int)(num_x*x_frac)*num_t2+(int)(num_t+num_t*t_frac)] = s_data [ (int)((s_num_t-1)-(int)k*time_frac) + s_num_t * (int)((s_num_b-1)-t) ];
                }
                if ( k % 100 == 0 ) std::cout << k << "    " << k*time_frac*0.004 << std::endl;
                propagate ( h , num_x , num_t2 , scale_x , scale_x , data_p2 , data_p1 , data_c , vel  );
                //if ( k*time_frac*0.004 > 0.02 )
                propagate ( h , num_x , num_t2 , scale_x , scale_x ,rdata_p2 ,rdata_p1 ,rdata_c , vel  );
                for ( int x(0)
                    ; x < num_x
                    ; ++x
                    )
                for ( int t(0)
                    ; t < num_t
                    ; ++t
                    )
                odata[t+x*num_t] += data_c[t+num_t+x*num_t2] * rdata_p1[t+num_t+x*num_t2];
                // if ( k*time_frac*0.004 > 1.5 ) break;
            }
            else
            {
            }
        }
    }
    
    glPointSize(20);
    glBegin(GL_POINTS);
    for ( int k(0)
        , x(0)
        ; x < num_x
        ; ++x
        )
    for ( int t(0)
        ; t < num_t
        ; ++t
        , ++k
        )
    {
        glColor3f
        ( 0.5 + data_c[t+num_t+num_t2*x] * scalar_1
        , 0.5 + rdata_c[t+num_t+num_t2*x] * scalar_2
        , 0.5 + odata[t+x*num_t] * scalar_3
        );
        glVertex3f 
        ( -1 + 2.0f * x / num_x
        , 0
        , -1 + 2.0f * t / num_t
        );
    }
    glEnd();
    
}

int main( int argc , char ** argv )
{

    if (argc < 3)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string file_name( argv[arg] );

    arg++;

    std::string seismic_file_name( argv[arg] );

    arg++;

    std::string seismic_data_file_name( argv[arg] );

    arg++;

    std::string output_file_name( argv[arg] );

    int num_iter = get_argument ( arg , argc , argv );

    float x_frac = get_argument ( arg , argc , argv );

    float t_frac = get_argument ( arg , argc , argv );

    SEPReader s_reader ( seismic_file_name . c_str () , false );
    s_num_t = s_reader . n1;
    s_num_b = s_reader . n2;
    s_data = new float [ s_num_b * s_num_t ];
    memset ( &s_data[0] , 0 , s_num_b * s_num_t );
    s_reader . OpenDataFile ( seismic_data_file_name . c_str () );
    s_reader . read_sepval ( & s_data [ 0 ] , s_reader . o1 , s_reader . o2 , s_reader . o3 , s_num_b * s_num_t );

    SEPReader reader ( file_name . c_str () );
    num_t = reader . n1;
    num_x = reader . n2 
          * reader . n3
          * reader . n4
          * reader . n5
          * reader . n6
          * reader . n7
          * reader . n8
          ;
    std::cout << num_t << std::endl;
    std::cout << num_x << std::endl;
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    num_t2 = num_t*3;
    vel = new float [ num_x * num_t2 ];

    for ( int x(0)
        , T
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t2
            ; ++t
            , ++k
            )
        {
            T = t - num_t;
            if ( T < 0 ) T = 0;
            if ( T >= num_t ) T = num_t - 1;
            vel[k] = data[T+x*num_t];
        }
    }

    data_p2 = new float [ num_x * num_t2 ];
    memset ( &data_p2[0] , 0 , num_x * num_t2 );

    data_p1 = new float [ num_x * num_t2 ];
    memset ( &data_p1[0] , 0 , num_x * num_t2 );

    data_c  = new float [ num_x * num_t2 ];
    memset ( &data_c [0] , 0 , num_x * num_t2 );

    rdata_p2 = new float [ num_x * num_t2 ];
    memset ( &rdata_p2[0] , 0 , num_x * num_t2 );

    rdata_p1 = new float [ num_x * num_t2 ];
    memset ( &rdata_p1[0] , 0 , num_x * num_t2 );

    rdata_c  = new float [ num_x * num_t2 ];
    memset ( &rdata_c [0] , 0 , num_x * num_t2 );

    odata  = new float [ num_x * num_t ];
    memset ( &odata [0] , 0 , num_x * num_t );

    x_frac = (3000.0+80*25.0)/9200.00;

    float width = 0.5;

    for ( int x(0)
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t2
            ; ++t
            , ++k
            )
        {
            float disp (sqrtf(pow_2(x-num_x*x_frac)+pow_2(t-num_t-num_t*t_frac)));
            data_p2 [ k ] = 0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
        }
    }

    float disp_x ( 25.0 / 4.0 );

    float near_offset ( 425.0 / 25.0 );

    for ( float r(5);r<1000;r+=5 )
    for ( float th(-M_PI);th<M_PI;th+=10.0/r )
    {
        recorder R;
        int x = (int)(num_x*x_frac+r*sin(th));
        int y = (int)(num_t+num_t*t_frac+r*cos(th));
        R . ind = x*num_t2+y;
        if ( x > 0 && x < num_x && y > 0 && y < num_t2 )
        recorders . push_back ( R );
    }

    for ( int k(0) ; 1 /*k < num_iter*/ ; ++k )
    {
        if ( k*time_frac > 50 )
        {
            if ( k*time_frac < 1.0*s_num_t )
            {
                for ( int t(s_num_b-1) ; t < s_num_b ; ++t )
                {
                    //std::cout << "x:" << (int)(num_x*x_frac-disp_x*(t+near_offset)) << std::endl;
                    //std::cout << "t:" << (int)(num_t+num_t*t_frac) << std::endl;
                    //std::cout << "ind:" << (int)(num_x*x_frac-disp_x*(t+near_offset))*num_t2+(int)(num_t+num_t*t_frac) << std::endl;
                    rdata_p2[(int)(num_x*x_frac-disp_x*(t+near_offset))*num_t2+(int)(num_t+num_t*t_frac)] = s_data [ (int)((s_num_t-1)-(int)(k*time_frac)+1) + s_num_t * (int)((s_num_b-1)-t) ];
                    rdata_p1[(int)(num_x*x_frac-disp_x*(t+near_offset))*num_t2+(int)(num_t+num_t*t_frac)] = s_data [ (int)((s_num_t-1)-(int)(k*time_frac)) + s_num_t * (int)((s_num_b-1)-t) ];
                }
            }
            else
            {
                break;
            }
        }
        if ( k % 100 == 0 ) std::cout << k << "    " << k*time_frac*0.004 << std::endl;
        //propagate ( h , num_x , num_t2 , scale_x , scale_x , data_p2 , data_p1 , data_c , vel  );
        propagate ( h , num_x , num_t2 , scale_x , scale_x ,rdata_p2 ,rdata_p1 ,rdata_c , vel  );

        for ( int k(0) ; k < recorders.size() ; ++k ) recorders[k].trace.push_back(rdata_c[recorders[k].ind]);
        //for ( int x(0)
        //    ; x < num_x
        //    ; ++x
        //    )
        //for ( int t(0)
        //    ; t < num_t
        //    ; ++t
        //    )
        //odata[t+x*num_t] += data_c[t+num_t+x*num_t2] * rdata_c[t+num_t+x*num_t2];
        // if ( k*time_frac*0.004 > 1.5 ) break;
    }

    for ( int k(0) ; k < recorders.size() ; ++k ) std::reverse ( recorders[k].trace.begin() , recorders[k].trace.end() );

    rdata_p2 = new float [ num_x * num_t2 ];
    memset ( &rdata_p2[0] , 0 , num_x * num_t2 );

    rdata_p1 = new float [ num_x * num_t2 ];
    memset ( &rdata_p1[0] , 0 , num_x * num_t2 );

    rdata_c  = new float [ num_x * num_t2 ];
    memset ( &rdata_c [0] , 0 , num_x * num_t2 );

#if 0
    for ( int k(0) ; 1 /*k < num_iter*/ ; ++k )
    {
        if ( k*time_frac > 50 )
        {
            if ( k*time_frac < 0.5*s_num_t )
            {
                for ( int t(0) ; t < s_num_b ; ++t )
                {
                    //rdata_p1[(int)(num_x*x_frac-disp_x*(t+near_offset))*num_t2+(int)(num_t+num_t*t_frac)] = s_data [ (int)((s_num_t-1)-(int)k*time_frac) + s_num_t * (int)((s_num_b-1)-t) ];
                    //data_p1[(int)(num_x*x_frac)*num_t2+(int)(num_t+num_t*t_frac)] = s_data [ (int)((s_num_t-1)-(int)k*time_frac) + s_num_t * (int)((s_num_b-1)-t) ];
                }
            }
            else
            {
                break;
            }
        }
        if ( k % 100 == 0 ) std::cout << k << "    " << k*time_frac*0.004 << std::endl;
        propagate ( h , num_x , num_t2 , scale_x , scale_x , data_p2 , data_p1 , data_c , vel  );
        propagate ( h , num_x , num_t2 , scale_x , scale_x ,rdata_p2 ,rdata_p1 ,rdata_c , vel  );
        for ( int x(0)
            ; x < num_x
            ; ++x
            )
        for ( int t(0)
            ; t < num_t
            ; ++t
            )
        odata[t+x*num_t] += data_c[t+num_t+x*num_t2] * rdata_c[t+num_t+x*num_t2];
        // if ( k*time_frac*0.004 > 1.5 ) break;
    }
#endif
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    init();
    glutMainLoop();

    //for ( int x(0)
    //    ; x < num_x
    //    ; ++x
    //    )
    //for ( int t(0)
    //    ; t < num_t
    //    ; ++t
    //    )
    //odata[t+x*num_t] += data_c[t+num_t+x*num_t2] * rdata_c[t+num_t+x*num_t2];

    float inv_max_w = 1 / find_max ( odata , num_x*num_t );
    float max_v = find_max ( data   , num_x*num_t  );

    for ( int x(0)
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t
            ; ++t
            , ++k
            )
        {
            data[k] -= 100*0.5 * odata[t+x*num_t] * max_v * inv_max_w;
        }
    }

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)&data[0] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    return 0;

}


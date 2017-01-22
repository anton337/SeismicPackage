#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"
#include "./ray_tracing.h"

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

void propagate ( float   h 
               , int     num_x 
               , int     num_t 
               , float * p2 
               , float * p1 
               , float * c 
               , float * vel )
{
    float fact ( h * h );
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

struct ray
{
    float  x ;
    float  y ;
    float  z ;
    float vx ;
    float vy ;
    float vz ;
    float  t ;
    float th_init;
    int done ;
    ray(){}
    ray ( float  _x 
        , float  _y 
        , float  _z 
        , float _vx 
        , float _vy 
        , float _vz 
        )
    :  x (  _x )
    ,  y (  _y )
    ,  z (  _z )
    , vx ( _vx )
    , vy ( _vy )
    , vz ( _vz )
    ,  t (   0 )
    , done ( 0 )
    {

    }
};

int main( int argc , char ** argv )
{

    if (argc < 3)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string smooth_file_name( argv[arg] );

    arg++;

    std::string exact_file_name( argv[arg] );

    arg++;

    std::string output_file_name( argv[arg] );

    int num_iter = get_argument ( arg , argc , argv );

    float x_frac = get_argument ( arg , argc , argv );

    float t_frac = get_argument ( arg , argc , argv );

    float h = 4*0.00005;

    SEPReader reader   ( smooth_file_name . c_str () );
    SEPReader e_reader (  exact_file_name . c_str () );
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
    float * travel_time = new float [ num_x * num_t * (240+96) ];
    memset ( &travel_time[0] , 0 , num_x * num_t * (240+96) );
    for ( int k(0) ; k < num_x * num_t * (240+96) ; ++k )
    {
        travel_time[k] = 3000000;
    }
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    float * e_data = new float [ num_x * num_t ];
    memset ( &e_data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );
    e_reader . read_sepval ( & e_data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    float scale_x ( 1 );
    float scale_y ( 1 );
    float scale_z ( 1 );

    velocity_model < float > model ( 1 
                                   , num_x 
                                   , num_t 
                                   , scale_x
                                   , scale_y 
                                   , scale_z 
                                   , data
                                   , e_data
                                   );



    raytracing :: RK4_solver < float 
                             , raytracing :: isotropic_exact_ray_step_functor_rk4 < float >
                             > solver;

    for ( int pos(0)
        ; pos < 240 + 96
        ; ++pos
        )
    {

        std::cout << "pos=" << pos << std::endl;

        x_frac = (float)(175+pos*25)/9200;

        int factor = 1;

        std::vector < ray > rays;

        for ( float th(-M_PI/2); th <= M_PI/2; th += 0.001 )
        {
            ray r ( scale_x*0 
                  , scale_y*num_x*x_frac 
                  , scale_z*num_t*t_frac
                  , 0 
                  , sin(th) 
                  , cos(th) 
                  );
            r . th_init = th;
            rays . push_back ( r );
        }

        int num_frac = 1;
        for ( int fr(1)
            ; fr <= num_frac
            ; fr += 1
            )
        {
            float frac ( fr / (float)num_frac );
            for ( int trace_iter(0) ; trace_iter < 5 ; ++trace_iter )
            {
                std::cout << "iter:" << trace_iter << " size:" << rays.size() << std::endl;

                std::vector < ray > new_rays;
                for ( std::size_t k(0) ; k < rays . size () ; ++k )
                if ( !rays[k] . done )
                {   
                    if ( k%200==0)
                    std::cout << "k=" << k << std::endl;
                    ray & r ( rays[k] );
                    float th = r . th_init;
                    r . x = scale_x*0;
                    r . y = scale_y*num_x*x_frac;
                    r . z = scale_z*num_t*t_frac;
                    r .vx = 0;
                    r .vy = sin(th);
                    r .vz = cos(th);
                    float slow ( 1 / model ( r . x 
                                           , r . y
                                           , r . z
                                           , r . vx
                                           , r . vy
                                           , r . vz
                                           ) 
                               );
                    r . vx *= slow;
                    r . vy *= slow;
                    r . vz *= slow;
                    for ( int iter(0)
                        ; iter < frac*num_iter*factor*scale_x
                        ; ++iter
                        )
                    {
                        solver . RK4 ( h/factor 
                                     , r . x 
                                     , r . vx
                                     , r . y 
                                     , r . vy
                                     , r . z 
                                     , r . vz
                                     , 1
                                     , r . t
                                     , model
                                     );
                        {
                            int x = r . x / scale_x;
                            int y = r . y / scale_y;
                            int z = r . z / scale_z;
                            if ( x < 0 ) {x = 0;break;}
                            if ( y < 0 ) {y = 0;break;}
                            if ( z < 0 ) {z = 0;break;}
                            if ( x >= 1 ) x = 0;
                            if ( y >= num_x ) {y = num_x - 1;break;}
                            if ( z >= num_t ) {z = num_t - 1;break;}
                            data [ z + num_t * y ] = std::max(0.0f,data[z+num_t*y]-10);
                            for ( int ind_y(-3)
                                ; ind_y <= 3
                                ; ++ind_y
                                )
                            {
                                for ( int ind_z(-3)
                                    ; ind_z <= 3
                                    ; ++ind_z
                                    )
                                {
                                    if ( y+ind_y >= 0 && y+ind_y < num_x )
                                    {
                                        if ( z+ind_z >= 0 && z+ind_z < num_t )
                                        {
                                            travel_time [ pos * (num_x*num_t) + z+ind_z + num_t * (y+ind_y) ] = std::min ( (float)r.t , travel_time [ pos * (num_x*num_t) + z+ind_z + num_t * (y+ind_y) ] );
                                        }
                                    }
                                }
                            }
                        }
                    }
                    r . done = true;
                }

                new_rays . clear ();
                for ( std::size_t k(0)
                    ; k + 1 < rays . size ()
                    ; ++k
                    )
                {
                    new_rays . push_back ( rays[k] );
                    if ( raytracing :: norm ( rays[k].x - rays[k+1].x
                                            , rays[k].y - rays[k+1].y
                                            , rays[k].z - rays[k+1].z
                                            ) 
                       > 5
                       )
                    {
                        if ( rays[k+1].th_init > rays[k].th_init + 1e-5 )
                        for ( float th ( rays[k].th_init )
                            ; th < rays[k+1].th_init
                            ; th += 0.5*(rays[k+1].th_init-rays[k].th_init)
                            )
                        {
                            ray r;
                            r . th_init = th;
                            new_rays . push_back (r);
                        }
                    }
                }
                new_rays . push_back ( rays[rays.size()-1] );

                rays . clear ();
                std::vector < ray > ( rays ) . swap ( rays );
                rays = new_rays;
                new_rays . clear ();
                std::vector < ray > ( new_rays ) . swap ( new_rays );

            }
        }

    }

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , 240+96
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)&travel_time[0] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * (240+96) );

    return 0;

}


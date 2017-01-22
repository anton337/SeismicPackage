#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"


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

    std::string output_file_name( argv[arg] );

    int num_iter = get_argument ( arg , argc , argv );

    float x_frac = get_argument ( arg , argc , argv );

    float t_frac = get_argument ( arg , argc , argv );

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
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    int num_t2 = num_t*3;
    float * vel = new float [ num_x * num_t2 ];

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
            data[T+x*num_t] = 3000.0;
            vel[k] = 3000.0;
        }
    }

    float * data_p2 = new float [ num_x * num_t2 ];
    memset ( &data_p2[0] , 0 , num_x * num_t2 );

    float * data_p1 = new float [ num_x * num_t2 ];
    memset ( &data_p1[0] , 0 , num_x * num_t2 );

    float * data_c  = new float [ num_x * num_t2 ];
    memset ( &data_c [0] , 0 , num_x * num_t2 );


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
            float disp (sqrtf(pow_2(x-num_x*x_frac)+pow_2(t-num_t2/3-num_t*t_frac)));
            data_p2 [ k ] = 0.0005*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
        }
    }

    float h = 0.00020;

    for ( int k(0) ; k < num_iter ; ++k )
    {
        if ( k % 100 == 0 ) std::cout << k << std::endl;
        propagate ( h , num_x , num_t2 , data_p2 , data_p1 , data_c , vel  );
    }

    float inv_max_w = 1 / find_max ( data_c , num_x*num_t2 );
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
            data[k] -= 0.5 * data_c[(t+num_t)+x*num_t2] * max_v * inv_max_w;
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


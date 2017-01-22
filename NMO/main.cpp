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

struct point
{
    float t;
    float vel;
    point ( float _t 
          , float _vel
          )
    : t   ( _t   )
    , vel ( _vel )
    {

    }
};

std::vector < point > vel_selection;

float get_velocity ( float t 
                   , std::vector < point > const & vel 
                   )
{
    if ( vel.size() == 0 ) 
    {
        std::cout << "error, empty velocity vector" << std::endl;
        exit(2);
    }
    float a;
    for ( int k(0)
        ; k+1 < vel . size ()
        ; ++k
        )
    {
        if ( t*1000 > vel[k].t && t*1000 <= vel[k+1].t )
        {
            {
                a = t*1000 - vel[k].t;
                a /= vel[k+1].t - vel[k].t;
                return (1-a) * vel[k] . vel + a * vel[k+1] . vel;
            }
        }
    }
    if ( t*1000 > vel[vel.size()-1] . t )
    {
        return vel[vel.size()-1] . vel;
    }
    return vel[0] . vel;
}

inline
float sqr ( float x )
{
    return x*x;
}

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

    std::string data_file_name( argv[arg] );

    arg++;

    std::string time_file_name( argv[arg] );

    arg++;

    std::string time_data_file_name( argv[arg] );

    arg++;

    std::string output_file_name( argv[arg] );

    arg++;

    std::string output_data_file_name( argv[arg] );

    float sample_rate= 0.001*get_argument( arg , argc , argv );
    float min_offset= get_argument( arg , argc , argv );
    float max_offset= get_argument( arg , argc , argv );

    int   side      = get_argument( arg , argc , argv );

    while ( arg+1 < argc )
    {
        float t     = get_argument( arg , argc , argv );
        float vel   = get_argument( arg , argc , argv );
        vel_selection . push_back ( point ( t , vel ) );
    }

    SEPReader reader ( file_name . c_str () , false );
    int num_t = reader . n1;
    int num_batch = reader . n2;
    int num_x = num_batch;
    int num_batches = reader . n3
                    * reader . n4
                    * reader . n5
                    * reader . n6
                    * reader . n7
                    * reader . n8
                    ;
    std::cout << "num_t = " << num_t << std::endl;
    std::cout << "num_x = " << num_x << std::endl;
    float * data = new float [ num_x * num_t * num_batches ];
    memset ( &data[0] , 0 , num_x * num_t * num_batches );

    float * data_out = new float [ num_x * num_t * num_batches ];
    memset ( &data_out[0] , 0 , num_x * num_t * num_batches );

    float * time = new float [ num_x * num_t * num_batches ];
    memset ( &time[0] , 0 , num_x * num_t * num_batches );

    float * weight = new float [ num_x * num_t * num_batches ];
    memset ( &weight[0] , 0 , num_x * num_t * num_batches );

    reader . OpenDataFile ( (data_file_name) . c_str() );

    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * num_batches );

    if ( side == -1 )
    {
        SEPReader reader ( time_file_name . c_str () , false );
        reader . OpenDataFile ( (time_data_file_name) . c_str() );
        reader . read_sepval ( & time [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * num_batches );
    }

    if ( side == 1 ) // forward
    {
        float shift;
        float w;
        float u;
        float t0;
        float offset;
        float vel_inv;
        float vel;
        float arg;
        for ( int n(0)
            , k(0)
            ; n < num_batches
            ; ++n
            )
        {
            for ( int x(0)
                ; x < num_x
                ; ++x
                )
            {
                offset = min_offset + x*(max_offset - min_offset)/num_x;
                for ( int t(0)
                    ; t < num_t
                    ; ++t
                    , ++k
                    )
                {
                    if ( 0.2f * fabs ( offset ) < t )
                    {
                        t0 = t*sample_rate;
                        vel = get_velocity ( t0 , vel_selection );
                        vel_inv = 1 / vel;
                        arg = ( sqr ( t0 ) + sqr ( offset * vel_inv ) );
                        if ( arg > 0 )
                        {
                            shift = sqrtf ( arg ) / sample_rate;
                            w = shift - (int)shift;
                            w *= w;
                            w *= w;
                            u = 1-w;
                            if ( shift >= 0 && shift+1 < num_t )
                            {
                                time[k] = shift;
                                data_out[k] = u*data[x*num_t+(int)shift] + w*data[x*num_t+(int)shift+1];
                            }
                            else
                            {
                                time[k] = 0;
                                data_out[k] = 0;
                            }
                        }
                        else
                        {
                            time[k] = 0;
                            data_out[k] = 0;
                        }
                    }
                }
            }
        }
    }
    else
    if ( side == -1 ) // reverse
    {
        float shift;
        float w;
        float u;
        float t0;
        float offset;
        float vel_inv;
        float vel;
        float arg;
        for ( int n(0)
            , k(0)
            ; n < num_batches
            ; ++n
            )
        {
            for ( int x(0)
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
                    {
                        shift = time[k]; 
                        w = shift - (int)shift;
                        w *= w;
                        w *= w;
                        u = 1-w;
                        if ( shift >= 0 && shift+1 < num_t )
                        {
                            data_out[x*num_t+(int)shift]   += u*data[k];
                            data_out[x*num_t+(int)shift+1] += w*data[k];
                            weight[x*num_t+(int)shift  ] += u;
                            weight[x*num_t+(int)shift+1] += w;
                        }
                    }
                }
                for ( int t(0)
                    ; t < num_t
                    ; ++t
                    )
                {
                    if ( weight[x*num_t+t] > 1e-5 )
                    {
                        data_out[x*num_t+t] /= weight[x*num_t+t];
                    }
                }
            }
        }
    }

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_data_file_name) . c_str()
                     );

    writer . OpenDataFile ( (output_data_file_name) . c_str() );

    writer . write_sepval ( (float*)data_out , reader . o1 , reader . o2 , reader . o3 , num_x * (num_t) * num_batches );

    if ( side == 1 )
    {
        SEPWriter writer ( time_file_name . c_str () 
                         , reader . o1 , reader . d1 , num_t
                         , reader . o2 , reader . d2 , reader . n2
                         , reader . o3 , reader . d3 , reader . n3
                         , reader . get_header_labels ()
                         , reader . get_sort_order ()
                         , (time_data_file_name) . c_str()
                         );

        writer . OpenDataFile ( (time_data_file_name) . c_str() );

        writer . write_sepval ( (float*)time , reader . o1 , reader . o2 , reader . o3 , num_x * (num_t) * num_batches );
    }

    return 0;

}


#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include <algorithm>
#include "sep_writer.h"
#include "sep_reader.h"
#include "viewer.h"

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
               , float   scale_x
               , float   scale_t
               , float * p2 
               , float * p1 
               , float * c 
               , float * vel 
               , bool reflecting_top_boundary = true
               , bool absorbing_boundary_condition = true
               )
{
    h /= scale_x;
    float fact ( h * h );
    for ( int t(1)
        ; t+1 < num_t
        ; ++t
        )
    {
        for ( int x(1)
            , ind_t
            ; x+1 < num_x
            ; ++x
            )
        {
            ind_t = t+x*num_t;
            if ( fabs ( p2[ind_t] ) > 1e-40 
              || fabs ( p1[ind_t+1] ) > 1e-40 
              || fabs ( p1[ind_t-1] ) > 1e-40 
              || fabs ( p1[ind_t+num_t] ) > 1e-40 
              || fabs ( p1[ind_t-num_t] ) > 1e-40 
               )
            {
                c[ind_t] = fact * vel[ind_t] * vel[ind_t] * ( -4*p1[ind_t] 
                                                                        +    p1[ind_t+num_t]
                                                                        +    p1[ind_t-num_t]
                                                                        +    p1[ind_t+1]
                                                                        +    p1[ind_t-1]
                                                                        )
                             + 2 * p1[ind_t]
                             -     p2[ind_t]
                             ;
            }
        }
    }
    if ( absorbing_boundary_condition )
    {
        // absorbing boundary t = 0
        if ( !reflecting_top_boundary )
        for ( int x(0)
            , t(0)
            ; x < num_x
            ; ++x
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+2+x*num_t] 
                           + 2 * p1[t+2+x*num_t] 
                               - p2[t+2+x*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+2+x*num_t] 
                                                      - p2[t+2+x*num_t] 
                                                      + p2[t+x*num_t] 
                                                      ) 
                                - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+2+x*num_t] 
                                                                            - 2 * p1[t+1+x*num_t] 
                                                                                + p1[t+x*num_t] 
                                                                                ) 
                            ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary t = num_t-1
        for ( int x(0)
            , t(num_t-1)
            ; x < num_x
            ; ++x
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t-2+x*num_t] 
                           + 2 * p1[t-2+x*num_t] 
                               - p2[t-2+x*num_t] 
                               + vel[t+x*num_t] * h * ( c[t-2+x*num_t] 
                                                     - p2[t-2+x*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t-2+x*num_t] 
                                                                           - 2 * p1[t-1+x*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary x = 0
        for ( int x(0)
            , t(0)
            ; t < num_t
            ; ++t
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+(x+2)*num_t] 
                           + 2 * p1[t+(x+2)*num_t] 
                               - p2[t+(x+2)*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+(x+2)*num_t] 
                                                     - p2[t+(x+2)*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+(x+2)*num_t] 
                                                                           - 2 * p1[t+(x+1)*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary x = num_x-1
        for ( int x(num_x-1)
            , t(0)
            ; t < num_t
            ; ++t
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+(x-2)*num_t] 
                           + 2 * p1[t+(x-2)*num_t] 
                               - p2[t+(x-2)*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+(x-2)*num_t] 
                                                     - p2[t+(x-2)*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+(x-2)*num_t] 
                                                                           - 2 * p1[t+(x-1)*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
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

/*
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
*/

void calculate_poyinting_vector ( vect        * p 
                                , float const * u
                                , int           num_x
                                , int           num_z
                                )
{
    int disp = 1;
    for ( int x(disp)
        ; x+disp < num_x
        ; ++x
        )
    {
        for ( int z(disp)
            ; z+disp < num_z
            ; ++z
            )
        {
            p[z+x*num_z].vx += (u[z+(x+disp)*num_z] - u[z+(x-disp)*num_z]);
            p[z+x*num_z].vz += (u[z+disp  +x*num_z] - u[z-disp+  x*num_z]);
            p[z+x*num_z].vy = 0;
        }
    }
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

float find_min ( float const * x , int num )
{
    float m ( 1e10 ) , t;
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        t = fabs(x[k]);
        if ( t < m )
        {
            m = t;
        }
    }
    return m;
}

struct source
{

    float pos_x; // m
    float pos_z; // m

    source ( int _pos_x , int _pos_z )
    : pos_x ( _pos_x )
    , pos_z ( _pos_z )
    {

    }

};

struct recorder
{

    float pos_x; // m
    float pos_z; // m

    recorder ( int _pos_x , int _pos_z )
    : pos_x ( _pos_x )
    , pos_z ( _pos_z )
    {

    }

};


inline float sqr ( float x )
{
    return x*x;
}

inline float weight ( float vz , float vx )
{
    float th ( atan2f ( vz , vx ) );
    return (fabs(th)>M_PI/2)?0:sqr(cosf(th));
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  wave propagation based synthetic trace data generator :
//
//  required fields:
//
//  input: 
//      - velocity SEP header file
//      - velocity model dimensions     // get from SEP file
//
//      - grid width x
//      - grid width z
//
//      - wavelet file
//
//      - sampling rate
//      - num_t                         // number of time samples to generate
//      - propagation step size
//
//      - num_rec
//      - num_sou
//      - rec_disp
//      - sou_disp
//      - rec_init
//      - sou_init
//
//     
//  output:   
//      - seismic SEP header file
//      - seismic SEP data files
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc , char ** argv )
{

    boost::thread       * thr = show ( argc , argv );


    if (argc < 3)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

//      - velocity SEP header file
    arg++;
    std::string velocity_file_name( argv[arg] );

//      - velocity model dimensions     // get from SEP file
    SEPReader vel_reader ( velocity_file_name . c_str () , false );
    int vel_num_z = vel_reader . n1;
    int vel_num_x = vel_reader . n2;
    float * vel_data = new float [ vel_num_x * vel_num_z ];
    memset ( &vel_data[0] , 0 , vel_num_x * vel_num_z );
    vel_reader . OpenDataFile ( (velocity_file_name + std::string("@")) . c_str () );
    vel_reader . read_sepval ( & vel_data [ 0 ] , vel_reader . o1 , vel_reader . o2 , vel_reader . o3 , vel_num_x * vel_num_z );

//      - grid width x
    float grid_width_x = get_argument ( arg , argc , argv );
//      - grid width z
    float grid_width_z = get_argument ( arg , argc , argv );

//      - wavelet file
    arg++;
    std::string wavelet_file_name( argv[arg] );
    SEPReader wav_reader ( wavelet_file_name . c_str () , false );
    int wav_num_t = wav_reader . n1;
    float * wav_data = new float [ wav_num_t ];
    memset ( &wav_data[0] , 0 , wav_num_t );
    wav_reader . OpenDataFile ( (wavelet_file_name + std::string("@")) . c_str () );
    wav_reader . read_sepval ( & wav_data [ 0 ] , wav_reader . o1 , wav_reader . o2 , wav_reader . o3 , wav_num_t );

//      - sampling rate
    float sample_rate = 0.001 * get_argument ( arg , argc , argv );
//      - num_t                         // number of time samples to generate
    float max_t = 0.001 * get_argument ( arg , argc , argv );
    int num_t = max_t / sample_rate;
//      - propagation step size
    float h = 0.001 * get_argument ( arg , argc , argv );

//      - num_rec
    int num_rec = get_argument ( arg , argc , argv );
//      - num_sou
    int num_sou = get_argument ( arg , argc , argv );
//      - rec_disp
    float rec_disp = get_argument ( arg , argc , argv );
//      - sou_disp
    float sou_disp = get_argument ( arg , argc , argv );
//      - rec_init
    float rec_init = get_argument ( arg , argc , argv );
//      - sou_init
    float sou_init = get_argument ( arg , argc , argv );

//      - seismic output SEP header file
    arg++;
    std::string seismic_file_header_name( argv[arg] );

//      - seismic output SEP data file
    arg++;
    std::string seismic_file_data_name( argv[arg] );



// pad velocity volume
    int vel_num_x2 = vel_num_x*1;
    int vel_num_z2 = vel_num_z*1;
    float * vel = new float [ vel_num_x2 * vel_num_z2 ];

    bool water_reflection = true;

    for ( int x(0)
        , _x(0)
        , k(0)
        ; x < vel_num_x2
        ; ++x
        )
    {
        for ( int z(0)
            , _z(0)
            ; z < vel_num_z2
            ; ++z
            , ++k
            )
        {
            vel[k] = vel_data[z+x*vel_num_z];
        }
    }

    float * data_p2 = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_p2[0] , 0 , vel_num_x2 * vel_num_z2 );

    float * data_p1 = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_p1[0] , 0 , vel_num_x2 * vel_num_z2 );

    float * data_c  = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_c [0] , 0 , vel_num_x2 * vel_num_z2 );

    float * recorder_data  = new float [ num_rec * num_sou * num_t ];
    memset ( &recorder_data [0] , 0 , num_rec * num_sou * num_t );

    vect * point_f = new vect [ vel_num_x2 * vel_num_z2 ];

    vect * point_o = new vect [ vel_num_x2 * vel_num_z2 ];

    set_num_x ( vel_num_x2 );
    set_num_z ( vel_num_z2 );
    set_data_ptr ( vel );
    set_reverse_ptr ( vel );
    set_scalar ( 0.4 / find_max ( vel , vel_num_x2 * vel_num_z2 ) );
    set_reverse_scalar ( 0.4 / find_max ( vel , vel_num_x2 * vel_num_z2 ) );
    set_wave_ptr ( data_c );
    // set_data_vec_ptr ( point_o );

    for ( int sou(0)
        ; sou < num_sou
        ; ++sou
        )
    {

        std::cout << "sou=" << sou << " pos=" << (int)((sou_init + sou_disp * sou) / grid_width_x) << std::endl;
        std::cout << sou_init << " " << sou_disp << std::endl;

        for ( int k(0)
            ; k < vel_num_x2 * vel_num_z2
            ; ++k
            )
        {
            point_f[k].vx = 0;
            point_f[k].vy = 0;
            point_f[k].vz = 0;
        }

        source S ( (int)((sou_init + sou_disp * sou) / grid_width_x)
                 , (int)(0)
                 );

        std::vector < recorder > recorders;

        for ( int rec(0)
            ; rec < num_rec
            ; ++rec
            )
        {
            recorder R ( (int)((sou_init + sou_disp * sou + rec_init + rec_disp * rec) / grid_width_x)
                       , (int)(1)
                       );
            std::cout << rec << " " << R . pos_x << std::endl;
            recorders . push_back ( R );
        }

        float sample_time ( sample_rate / h );

        for ( int k(0) 
            ; 1 
            ; ++k 
            )
        {
            if ( k < num_t * sample_time )
            {
                // feed wavelet 
                if ( (int)(k*sample_time) < wav_num_t )
                {
                    if ( (int)(k*sample_time) >= 1 )
                    {
                        data_p2[(int)(S.pos_z  )+(int)(S.pos_x)*vel_num_z2] =  wav_data [ (int)(k*sample_time)-1 ];
                        data_p2[(int)(S.pos_z+1)+(int)(S.pos_x)*vel_num_z2] = -wav_data [ (int)(k*sample_time)-1 ];
                    }
                    data_p1[(int)(S.pos_z  )+(int)(S.pos_x)*vel_num_z2] =  wav_data [ (int)(k*sample_time) ];
                    data_p1[(int)(S.pos_z+1)+(int)(S.pos_x)*vel_num_z2] = -wav_data [ (int)(k*sample_time) ];
                }
            }
            else
            {
                break;
            }
            propagate ( h 
                      , vel_num_x2 
                      , vel_num_z2 
                      , grid_width_x 
                      , grid_width_z 
                      , data_p2 
                      , data_p1 
                      , data_c
                      , vel  
                      );
            set_wave_scalar ( 0.4 / find_max ( data_c , vel_num_x2 * vel_num_z2 ) );
            set_timer ( 1000 * k * h );
            // calculate_poyinting_vector ( point_f
            //                            , data_p1
            //                            , vel_num_x2
            //                            , vel_num_z2
            //                            );
            // for ( int x(1)
            //     ; x+1 < vel_num_x2
            //     ; ++x
            //     )
            // for ( int z(1)
            //     ; z+1 < vel_num_z2
            //     ; ++z
            //     )
            // {
            //     {
            //         point_o[z+vel_num_z2*x].vx = -k*(point_f[z+vel_num_z2*(x)].vx*data_p1[z+vel_num_z2*(x)]);
            //         point_o[z+vel_num_z2*x].vz = -k*(point_f[z+vel_num_z2*(x)].vz*data_p1[z+vel_num_z2*(x)]);
            //     }
            // }
            for ( std::size_t rec(0)
                ; rec < recorders . size ()
                ; ++rec
                )
            {
                int t ( k * sample_time );
                int ind ( (int)(recorders[rec].pos_z) + (int)(recorders[rec].pos_x)*vel_num_z2 );
                recorder_data[ sou*num_rec*num_t + rec*num_t + t ] = /* weight ( -point_o[ind].vx , -point_o[ind].vz ) */ ( data_p1[(int)(recorders[rec].pos_z)+(int)(recorders[rec].pos_x)*vel_num_z2] );
            }
        }

        for ( int k(0)
            , size = vel_num_x2 * vel_num_z2
            ; k < size
            ; ++k
            )
        {
            data_p2[k] = 0;
            data_p1[k] = 0;
            data_c [k] = 0;
        }

    }

    SEPWriter writer ( seismic_file_header_name . c_str () 
                     , 0 , 1000*sample_rate , num_t
                     , 1 , 1 , num_rec
                     , 1 , 1 , num_sou
                     , wav_reader . get_header_labels ()
                     , wav_reader . get_sort_order ()
                     , seismic_file_data_name . c_str()
                     );

    writer . OpenDataFile ( (seismic_file_data_name) . c_str() );

    writer . write_sepval ( (float*)&recorder_data[0] , writer . o1 , writer . o2 , writer . o3 , num_sou*num_rec*num_t );



    thr -> join ();

    return 0;

}


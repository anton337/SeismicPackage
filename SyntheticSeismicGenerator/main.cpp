#include <iostream>
#include <string.h>
#include <math.h>
#include <map>
#include "sep_writer.h"

/////////////////////////////////////////////////////////////////////////////////////////
//
// ORMSBY
//
// f1  = low  - cut  frequency
// f2  = low  - pass frequency
// f3  = high - pass frequency
// f4  = high - cut  frequency
// tau = time shift
//
// A(t) = pi [ ( f4^2 sinc^2 (pi f4 (t - tau)) - f3^2 sinc^2 (pi f3 (t - tau)) ) 
//             / (f4 - f3) 
//           - ( f2^2 sinc^2 (pi f2 (t - tau)) - f1^2 sinc^2 (pi f1 (t - tau)) ) 
//             / (f2 - f1) 
//           ]
//
/////////////////////////////////////////////////////////////////////////////////////////
//
// RICKER
//
// f   = peak frequency
// tau = time shift
//
// A(t) = [ 1 - 2 pi^2 f^2 (t - tau)^2 ] exp ( - pi^2 f^2 (t - tau)^2 )
//
/////////////////////////////////////////////////////////////////////////////////////////
inline float sqr ( float const & x )
{
    return x*x;
}
inline float sinc_sqr ( float const & x )
{
    return (x!=0)?sqr(sinf(x)/x):1;
}
void generate_ormsby_data ( int     num_t
                          , float   sample_rate
                          , float   f1
                          , float   f2
                          , float   f3
                          , float   f4
                          , float   A
                          , float   tau
                          , float * data  
                          )
{
    float u4 = M_PI * f4 * f4 * A * sample_rate / ( f4 - f3 );
    float u3 = M_PI * f3 * f3 * A * sample_rate / ( f4 - f3 );
    float u2 = M_PI * f2 * f2 * A * sample_rate / ( f2 - f1 );
    float u1 = M_PI * f1 * f1 * A * sample_rate / ( f2 - f1 );
    float w4 = M_PI * f4;
    float w3 = M_PI * f3;
    float w2 = M_PI * f2;
    float w1 = M_PI * f1;

    for ( int t(0)
        ; t < num_t
        ; ++t
        )
    {
        data[t] += ( u4 * sinc_sqr ( w4 * ( sample_rate * t - tau ) ) - u3 * sinc_sqr ( w3 * ( sample_rate * t - tau ) ) )
                 - ( u2 * sinc_sqr ( w2 * ( sample_rate * t - tau ) ) - u1 * sinc_sqr ( w1 * ( sample_rate * t - tau ) ) )
                 ;
    }

}
void generate_ricker_data ( int     num_t
                          , float   sample_rate
                          , float   f
                          , float   A
                          , float   tau
                          , float * data  
                          )
{
    float u = M_PI * M_PI * f * f;
    float d;

    for ( int t(0)
        ; t < num_t
        ; ++t
        )
    {
        d = sqr(sample_rate * t - tau);
        data[t] += A * ( 1 - 2 * u * d ) * expf ( - u * d );
    }

}

void generate_hyperbolic_ormsby_event ( int                              num_x
                                      , int                              num_t
                                      , float                            sample_rate
                                      , float                            f1
                                      , float                            f2
                                      , float                            f3
                                      , float                            f4
                                      , float                            A
                                      , float                            t0
                                      , float                            apex
                                      , float                            vel
                                      , float                            multiple_coefficient
                                      , std::string                    & offset_str
                                      , float                          * data  
                                      , std::map < std::string , int > & header_ind
                                      , float                          * header
                                      )
{
    int offset_ind ( header_ind [ offset_str ] );
    float tau;
    float offset;
    float vel_inv = 1 / vel;
    for ( int mult(1)
        ; mult < 15
        ; ++mult
        , A *= multiple_coefficient
        )
    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    {
        offset = header [ offset_ind + x * header_ind . size () ];
        tau = sqrtf ( sqr ( mult * t0 ) + sqr ( ( offset - apex ) * vel_inv ) );
        generate_ormsby_data ( num_t
                             , sample_rate
                             , f1
                             , f2
                             , f3
                             , f4
                             , A
                             , tau
                             , & data [ x * num_t ]
                             );
    }
}

void generate_hyperbolic_ricker_event ( int                              num_x
                                      , int                              num_t
                                      , float                            sample_rate
                                      , float                            f
                                      , float                            A
                                      , float                            t0
                                      , float                            apex
                                      , float                            vel
                                      , float                            multiple_coefficient
                                      , std::string                    & offset_str
                                      , float                          * data  
                                      , std::map < std::string , int > & header_ind
                                      , float                          * header
                                      )
{
    int offset_ind ( header_ind [ offset_str ] );
    float tau;
    float offset;
    float vel_inv = 1 / vel;
    for ( int mult(1)
        ; mult < 15
        ; ++mult
        , A *= multiple_coefficient
        )
    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    {
        offset = header [ offset_ind + x * header_ind . size () ];
        tau = sqrtf ( sqr ( mult * t0 ) + sqr ( ( offset - apex ) * vel_inv ) );
        generate_ricker_data ( num_t
                             , sample_rate
                             , f
                             , A
                             , tau
                             , & data [ x * num_t ]
                             );
    }
}

void generate_header ( int                              num_x
                     , int                              num_ffid
                     , int                              num_chan
                     , float                            delta_x
                     , std::map < std::string , int > & header_ind
                     , float                          * header
                     )
{
    int ffid = 1;
    for ( int x(0)
        ; ffid <= num_ffid
        ; ++ffid
        )
    {
        int chan = 1;
        for ( 
            ; chan <= num_chan
            ; ++chan
            , ++x
            )
        {
            header [ x * header_ind . size () + header_ind [ "FFID"   ] ] = ffid;
            header [ x * header_ind . size () + header_ind [ "CHAN"   ] ] = chan;
            header [ x * header_ind . size () + header_ind [ "OFFSET" ] ] = delta_x * chan;
        }
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

int main( int argc , char ** argv )
{

    // argument format
    /*
       ./a.out seismic_output_file_name.SEP                             \
               headers_output_file_name.SEP                             \
               num_ffid                                                 \
               num_chan                                                 \
               delta_x     [meters]                                     \
               measurement_time [ms]                                    \
               sample_rate [ms]                                         \
               event_type_1 event_shape_1 [shape_args_1] [event_args_1] \
               event_type_2 event_shape_2 [shape_args_2] [event_args_2] \
               event_type_3 event_shape_3 [shape_args_3] [event_args_3] \
               .
               .
               .
               additional events
      
      
      
       event_shape 1 : linear
       event_shape 2 : parabolic
       event_shape 3 : hyperboic
      
       apex [meters]
      
       time shift [ms]
      
       vel [m/s]
      
       amplitude
       
       multiple_coefficient [%]                                 
      
      
       event_type 1 : Ricker
      
       frequency [Hz]
      
      
      
       event_type 1 : ORMSBY
      
       f1 [Hz] 
       f2 [Hz]
       f3 [Hz]
       f4 [Hz]
      
    */

    if (argc == 1)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;
    std::string file_name( argv[arg] );

    arg++;
    std::string hdrs_name( argv[arg] );

    int num_ffid = get_argument ( arg , argc , argv ) ;
    std::cout << "num_ffid : " << num_ffid << std::endl;
    int num_chan = get_argument ( arg , argc , argv ) ;
    std::cout << "num_chan : " << num_chan << std::endl;

    float delta_x = get_argument ( arg , argc , argv ) ;
    std::cout << "delta_x : " << delta_x << std::endl;

    int num_x = num_ffid * num_chan;
    float measurement_time = 0.001 * get_argument ( arg , argc , argv ) ;
    std::cout << "measurement_time : " << measurement_time << std::endl;

    float sample_rate = 0.001 * get_argument ( arg , argc , argv ) ;
    std::cout << "sample_rate : " << sample_rate << std::endl;

    int num_t = measurement_time / sample_rate;
    if ( num_t % 2 == 1 ) num_t += 1;
    
    float * data = new float [ num_x * num_t ];

    std::string offset_str = "OFFSET";

    std::map < std::string , int > header_ind;
    header_ind [ "FFID"   ] = 0;
    header_ind [ "CHAN"   ] = 1;
    header_ind [ "OFFSET" ] = 2;
    header_ind [ "SOU_X"  ] = 3;
    header_ind [ "REC_X"  ] = 4;
    int num_header = header_ind . size ();

    std::vector < std::string > headers;
    headers . push_back ( "FFID"   );
    headers . push_back ( "CHAN"   );
    headers . push_back ( "OFFSET" );
    headers . push_back ( "SOU_X"  );
    headers . push_back ( "REC_X"  );

    std::vector < std::string > sort_order;
    sort_order . push_back ( "FFID"   );
    sort_order . push_back ( "CHAN"   );

    memset ( &data   [0] , 0 ,   num_x * num_t      );
    float * header = new float [ num_x * num_header ];
    memset ( &header [0] , 0 ,   num_x * num_header );

    generate_header ( num_x
                    , num_ffid
                    , num_chan
                    , delta_x
                    , header_ind
                    , header
                    );

    int event = 0;

    while ( arg+1 < argc )
    {
        std::cout << "event : " << ++event << std::endl;

        int type  = get_argument ( arg , argc , argv );

        int shape = get_argument ( arg , argc , argv );

        switch ( shape )
        {
            case 1: // linear
            {
                std::cout << "linear mode currently undefined" << std::endl;
                exit(1);
                break;
            }
            case 2: // parabolic
            {
                std::cout << "parabolic mode currently undefined" << std::endl;
                exit(1);
                break;
            }
            case 3: // hyperbolic
            {

                break;
            }
            default:
            {
                std::cout << "undefined shape" << std::endl;
                exit(1);
            }
        }
        float apex = get_argument ( arg , argc , argv );
        float t0   = get_argument ( arg , argc , argv );
        float vel  = get_argument ( arg , argc , argv );
        float A    = get_argument ( arg , argc , argv );
        float multiple_coefficient = get_argument ( arg , argc , argv );

        switch ( type )
        {
            case 1: // ricker
            {
                float f  = get_argument ( arg , argc , argv ); 
                generate_hyperbolic_ricker_event ( num_x
                                                 , num_t
                                                 , sample_rate
                                                 , f
                                                 , A
                                                 , t0 * 0.001
                                                 , apex
                                                 , vel
                                                 , multiple_coefficient
                                                 , offset_str
                                                 , data  
                                                 , header_ind
                                                 , header
                                                 );
                break;
            }
            case 2: // ormsby
            {
                float f1 = get_argument ( arg , argc , argv );
                std::cout << "f1:" << f1 << std::endl;
                float f2 = get_argument ( arg , argc , argv ); 
                std::cout << "f2:" << f2 << std::endl;
                float f3 = get_argument ( arg , argc , argv ); 
                std::cout << "f3:" << f3 << std::endl;
                float f4 = get_argument ( arg , argc , argv ); 
                std::cout << "f4:" << f4 << std::endl;
                generate_hyperbolic_ormsby_event ( num_x
                                                 , num_t
                                                 , sample_rate
                                                 , f1
                                                 , f2
                                                 , f3
                                                 , f4
                                                 , A
                                                 , t0 * 0.001
                                                 , apex
                                                 , vel
                                                 , multiple_coefficient
                                                 , offset_str
                                                 , data  
                                                 , header_ind
                                                 , header
                                                 );
                break;
            }
            default:
            {
                std::cout << "undefined type" << std::endl;
                exit(1);
            }
        }

    }

    if ( arg+1 != argc )
    {
        std::cout << "not enough arguments" << std::endl;
        exit(1);
    }

    SEPWriter writer ( file_name . c_str () 
                     , 0 , 1 , num_t
                     , 1 , 1 , num_chan
                     , 1 , 1 , num_ffid
                     , headers
                     , sort_order
                     , (file_name + std::string("@")) . c_str()
                     , (hdrs_name) . c_str()
                     );

    writer . OpenDataFile ( (file_name + std::string("@")) . c_str() );

    writer . write_sepval ( data , 0 , 1 , 1 , num_x * num_t );

    SEPWriter hdr_writer ( hdrs_name . c_str () 
                         , 0 , 1 , header_ind . size ()
                         , 1 , 1 , num_chan
                         , 1 , 1 , num_ffid
                         , headers
                         , sort_order
                         , (hdrs_name + std::string("@")) . c_str()
                         );

    hdr_writer . OpenDataFile ( (hdrs_name + std::string("@")) . c_str() );

    hdr_writer . write_sepval ( header , 0 , 1 , 1 , num_x * header_ind . size () );

    return 0;
}


#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"
#include "matrix.h"
#include "solvers.h"

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

void construct_forward_operator_linear ( int num_h 
                                       , int num_p 
                                       , int num_t
                                       , float f
                                       , float min_H
                                       , float delta_H
                                       , float min_P
                                       , float delta_P
                                       , fftwf_complex * a 
                                       )
{
    float P = min_P , H;
    float arg;
    float scal ( 2*M_PI/num_t ) ;
    for ( int p(0)
        , k(0)
        ; p < num_p
        ; ++p
        , P += delta_P
        )
    {
        H = min_H;
        for ( int h(0)
            ; h < num_h
            ; ++h
            , H += delta_H
            , ++k
            )
        {
            arg = scal * f * H * P;
            a[k][0] = cosf ( arg );
            a[k][1] = sinf ( arg );
        }
    }
}

void construct_gram_operator_parabolic ( int num_h 
                                       , int num_p 
                                       , int num_t
                                       , float f
                                       , float min_H
                                       , float delta_H
                                       , float min_P
                                       , float delta_P
                                       , fftwf_complex * T 
                                       )
{
    float P = min_P , H;
    float arg;
    float factor = 2*f*M_PI/(num_t);
    float qx = factor*delta_P , wx;
    float r1 , i1;
    float rd , id;
    float r1_;
    H = min_H;
    for ( int p = 0
        ; p < num_p
        ; ++p
        )
    {
        T[p][0] = 0;
        T[p][1] = 0;
    }
    for ( int h = 0
        ; h < num_h
        ; ++h
        , H += delta_H
        ) 
    {
        r1 = 1;
        i1 = 0;
        wx = qx*H*H;
        rd =  cos(wx);
        id =  sin(wx);
        T[0][0] += r1;
        for ( int p = 1
            ; p < num_p
            ; ++p
            ) 
        {
            r1_ = r1;
            r1 = r1 * rd - i1 * id;
            i1 = r1_* id + i1 * rd;
            T[p][0] += r1;
            T[p][1] += i1;
        }
    }
    T[0][0] += 0.1*num_p;
}

void construct_gram_operator_linear ( int num_h 
                                    , int num_p 
                                    , int num_t
                                    , float f
                                    , float min_H
                                    , float delta_H
                                    , float min_P
                                    , float delta_P
                                    , fftwf_complex * T 
                                    )
{
    float P = min_P , H;
    float arg;
    float factor      = 2*f*M_PI/(num_t);
    float qx = factor*delta_P , wx;
    float r1 , i1;
    float rd , id;
    float r1_;
    H = min_H;
    for ( int h = 0
        ; h < num_h
        ; ++h
        , H += delta_H
        ) 
    {
        r1 = 1;
        i1 = 0;
        wx = qx*H;
        rd =  cos(wx);
        id = -sin(wx);
        T[0][0] += r1;
        T[0][1] += i1;
        for ( int p = 1
            ; p < num_p
            ; ++p
            ) 
        {
            r1_ = r1;
            r1 = r1 * rd - i1 * id;
            i1 = r1_* id + i1 * rd;
            T[p][0] += r1;
            T[p][1] += i1;
        }
    }
}

void construct_forward_operator_parabolic ( int num_h 
                                          , int num_p 
                                          , int num_t
                                          , float f
                                          , float min_H
                                          , float delta_H
                                          , float min_P
                                          , float delta_P
                                          , fftwf_complex * a 
                                          )
{
    float P = min_P , H;
    float arg;
    float scal ( 2*M_PI/num_t ) ;
    for ( int p(0)
        , k(0)
        ; p < num_p
        ; ++p
        , P += delta_P
        )
    {
        H = min_H;
        for ( int h(0)
            ; h < num_h
            ; ++h
            , H += delta_H
            , ++k
            )
        {
            arg = scal * f * H * H * P;
            a[k][0] = cosf ( arg );
            a[k][1] = sinf ( arg );
        }
    }
}



template < typename T >
inline T taylor_hyperbolic_expansion_3 ( T const & h 
                                       , T const & p
                                       )
{
    T q ( p * h * h );
    return           q 
           - 0.5f*   q*q      
           + 0.5f*   q*q*q    
           - 0.625f* q*q*q*q  
           + 0.875f* q*q*q*q*q
           ;
}

void construct_forward_operator_hyperbolic ( int num_h 
                                           , int num_p 
                                           , int num_t
                                           , float f
                                           , float min_H
                                           , float delta_H
                                           , float min_P
                                           , float delta_P
                                           , fftwf_complex * a 
                                           )
{
    float P = min_P , H;
    float arg;
    float scal ( 2*M_PI/num_t ) ;
    for ( int p(0)
        , k(0)
        ; p < num_p
        ; ++p
        , P += delta_P
        )
    {
        H = min_H;
        for ( int h(0)
            ; h < num_h
            ; ++h
            , H += delta_H
            , ++k
            )
        {
            arg = scal * f * taylor_hyperbolic_expansion_3 ( H , P );
            a[k][0] = cosf ( arg );
            a[k][1] = sinf ( arg );
        }
    }
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

    std::string output_file_name( argv[arg] );

    int num_p ( get_argument ( arg , argc , argv ) );

    float min_h ( get_argument ( arg , argc , argv ) );
    float max_h ( get_argument ( arg , argc , argv ) );

    float min_p ( get_argument ( arg , argc , argv ) );
    float max_p ( get_argument ( arg , argc , argv ) );

    int mode ( get_argument ( arg , argc , argv ) );

    SEPReader reader ( file_name . c_str () );
    int num_t = reader . n1;
    int num_batch = reader . n2;
    int num_x = reader . n2 
              * reader . n3
              * reader . n4
              * reader . n5
              * reader . n6
              * reader . n7
              * reader . n8
              ;
    std::cout << "num_t : " << num_t << std::endl;
    std::cout << "num_x : " << num_batch << std::endl;
    std::cout << "num_p : " << num_p << std::endl;
    int num_batches = reader . n3
                    * reader . n4
                    * reader . n5
                    * reader . n6
                    * reader . n7
                    * reader . n8
                    ;

    float delta_h ( ( max_h - min_h ) / num_batch );
    float delta_p ( ( max_p - min_p ) / num_p );

    if ( num_t % 2 != 0 )
    {
        std::cout << "ERROR : num_t should be even" << std::endl;
        exit(1);
    }
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    float * odata = new float [ num_p * num_t * num_batches ];
    memset ( &odata[0] , 0 , num_p * num_t * num_batches );

    int num_f = num_t / 2 + 1;

    float * tmp = new float [ num_t ];
    fftwf_complex * ftmp = new fftwf_complex [ num_f ];

    float * tmp2 = new float [ num_t ];
    fftwf_complex * ftmp2 = new fftwf_complex [ num_f ];

    fftwf_complex * ftmp_total  = new fftwf_complex [ num_f * num_batch ];
    fftwf_complex * ftmp2_total = new fftwf_complex [ num_f * num_p ];

    fftwf_complex * forward_array  = new fftwf_complex [ num_batch * num_p  ];
    fftwf_complex * toeplitz_array = new fftwf_complex [ num_p              ];
    fftwf_complex * freq           = new fftwf_complex [ num_batch          ];
    fftwf_complex * freq2          = new fftwf_complex [ num_p              ];
    fftwf_complex * freq3          = new fftwf_complex [ num_p              ];

    fftwf_plan plan ( fftwf_plan_dft_r2c_1d ( num_t
                                            , & tmp[0]
                                            , &ftmp[0]
                                            , FFTW_ESTIMATE
                                            )
                    );

    fftwf_plan inv_plan ( fftwf_plan_dft_c2r_1d ( num_t
                                                , &ftmp2[0]
                                                , & tmp2[0]
                                                , FFTW_ESTIMATE
                                                )
                        );

    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        std::cout << "pass 1" << std::endl;
        for ( int h(0)
            ; h < num_batch
            ; ++h
            )
        {
            mem_cpy ( & tmp[0] , & data[x*num_t*num_batch+h*num_t] , num_t );
            fftwf_execute ( plan );
            mem_cpy ( & ((float*)ftmp_total)[2*h*num_f] , & ((float*)ftmp)[0] , 2*num_f );
        }
        std::cout << "pass 2" << std::endl;
        for ( int f(0)
            ; f < num_f
            ; ++f
            )
        {
            std::cout << "f=" << f << std::endl;
            for ( int h(0)
                ; h < num_batch
                ; ++h
                )
            {
                freq[h][0] = ftmp_total[ num_f*h + f ][0];
                freq[h][1] = ftmp_total[ num_f*h + f ][1];
            }
            switch ( mode )
            {
            case 1:
            {
            construct_forward_operator_linear ( num_batch 
                                              , num_p 
                                              , num_t 
                                              , f
                                              , min_h 
                                              , delta_h 
                                              , min_p 
                                              , delta_p 
                                              , forward_array 
                                              );
            matrix < fftwf_complex > L ( num_p 
                                       , num_batch 
                                       , forward_array 
                                       );
            L ( freq , freq2 );
            }
            break;
            case 2:
            {
            construct_forward_operator_parabolic ( num_batch 
                                                 , num_p 
                                                 , num_t 
                                                 , f
                                                 , min_h 
                                                 , delta_h 
                                                 , min_p 
                                                 , delta_p 
                                                 , forward_array 
                                                 );
            matrix < fftwf_complex > L ( num_p 
                                       , num_batch 
                                       , forward_array 
                                       );
            L ( freq , freq2 );
            }
            break;
            case 10:
            {
            construct_forward_operator_linear ( num_batch 
                                              , num_p 
                                              , num_t 
                                              , f
                                              , min_h 
                                              , delta_h 
                                              , min_p 
                                              , delta_p 
                                              , forward_array 
                                              );
            construct_gram_operator_linear ( num_batch 
                                           , num_p 
                                           , num_t 
                                           , f
                                           , min_h 
                                           , delta_h 
                                           , min_p 
                                           , delta_p 
                                           , toeplitz_array 
                                           );
            matrix < fftwf_complex > L ( num_p 
                                       , num_batch 
                                       , forward_array 
                                       );
            L ( freq , freq2 );
            }
            break;
            case 20:
            {
            construct_forward_operator_parabolic ( num_batch 
                                                 , num_p 
                                                 , num_t 
                                                 , f
                                                 , min_h 
                                                 , delta_h 
                                                 , min_p 
                                                 , delta_p 
                                                 , forward_array 
                                                 );
            construct_gram_operator_parabolic ( num_batch 
                                              , num_p 
                                              , num_t 
                                              , f
                                              , min_h 
                                              , delta_h 
                                              , min_p 
                                              , delta_p 
                                              , toeplitz_array 
                                              );
            matrix < fftwf_complex > L ( num_p 
                                       , num_batch 
                                       , forward_array 
                                       );
            L ( freq , freq3 );
            SymmetricConjugateGradientSolver < toeplitz_matrix < fftwf_complex > , fftwf_complex > solver;
            toeplitz_matrix < fftwf_complex > M ( num_p , num_p , toeplitz_array );
            solver ( M , freq2 , freq3 , 20 );
            }
            break;
            case -20:
            {
            std::cout << num_batch << " " << num_p << std::endl;
            construct_forward_operator_parabolic ( num_batch 
                                                 , num_p 
                                                 , num_t 
                                                 , f
                                                 , min_h 
                                                 , delta_h 
                                                 , min_p 
                                                 , delta_p 
                                                 , forward_array 
                                                 );
            matrix < fftwf_complex > L ( num_p 
                                       , num_batch 
                                       , forward_array 
                                       );
            L . adjoint( freq , freq2 );
            }
            break;
            }
            for ( int p(0)
                ; p < num_p
                ; ++p
                )
            {
                ftmp2_total[ num_f*p + f ][0] = freq2[p][0];
                ftmp2_total[ num_f*p + f ][1] = freq2[p][1];
            }
        }
        std::cout << "pass 3" << std::endl;
        for ( int p(0)
            ; p < num_p
            ; ++p
            )
        {
            mem_cpy ( & ((float*)ftmp2)[0] , & ((float*)ftmp2_total)[2*p*num_f] , 2*num_f );
            fftwf_execute ( inv_plan );
            mem_cpy ( & odata[x*num_t*num_p+p*num_t] , & tmp2[0] , num_t );
        }
        std::cout << "pass 4" << std::endl;
    }

    fftwf_destroy_plan ( plan );

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , num_p
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)odata , reader . o1 , reader . o2 , reader . o3 , num_batches * num_p * num_t );

    return 0;

}


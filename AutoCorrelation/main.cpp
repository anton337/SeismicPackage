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

inline void mem_mult ( fftwf_complex * out , fftwf_complex * in , int num )
{
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        float br( out[k][0] );
        float bi( out[k][1] );
        float ar(  in[k][0] );
        float ai( -in[k][1] );
        out[k][0] = ar*br - ai*bi;
        out[k][1] = ar*bi + ai*br;
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
    int num_batches = reader . n3
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

    int num_f = num_t / 2 + 1;
    fftwf_complex * fdata = new fftwf_complex [ num_x * num_f ];

    float * tmp = new float [ num_t * 1 ];
    fftwf_complex * ftmp = new fftwf_complex [ num_f * 1 ];

    fftwf_plan plan ( fftwf_plan_dft_r2c_2d ( 1
                                            , num_t
                                            , & tmp[0]
                                            , &ftmp[0]
                                            , FFTW_ESTIMATE
                                            )
                    );
std::cout << "pass 1" << std::endl;
    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        for ( int n(0)
            ; n < num_batch
            ; ++n
            )
        {
            mem_cpy ( & tmp[0] , & data[x*num_t*num_batch+num_t*n] , num_t*1 );
            fftwf_execute ( plan );
            mem_cpy < float > ( (float*)(&fdata[x*num_f*num_batch+num_f*n]) , (float*)(&ftmp[0]) , 2*num_f*1 );
        }
    }
std::cout << "pass 2" << std::endl;

    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    {
        mem_mult ( (&fdata[x*num_f*1]) , (&fdata[x*num_f*1]) , num_f*1 );
    }
std::cout << "pass 3" << std::endl;

    fftwf_destroy_plan ( plan );

    fftwf_plan iplan ( fftwf_plan_dft_c2r_2d ( 1
                                             , num_t
                                             , &ftmp[0]
                                             , & tmp[0]
                                             , FFTW_ESTIMATE
                                             )
                     );

std::cout << "pass 4" << std::endl;
    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        for ( int n(0)
            ; n < num_batch
            ; ++n
            )
        {
            mem_cpy < float > ( (float*)(&ftmp[0]) , (float*)(&fdata[x*num_f*num_batch+num_f*n]) , 2*num_f*1 );
            fftwf_execute ( iplan );
            mem_cpy ( & data[x*num_t*num_batch+num_t*n] , & tmp[0] , num_t*1 );
            for ( int k(0)
                , size ( num_t * 1 )
                ; k < size
                ; ++k
                )
            {
                data[x*size + k] /= num_t;
            }
        }
    }
std::cout << "pass 5" << std::endl;

    fftwf_destroy_plan ( iplan );

    {
        float maximum ( 0 );
        for ( int k(0) ; k < num_x * num_t ; ++k )
        {
            if ( fabs(data[k]) > maximum )
            {
                maximum = fabs(data[k]);
            }
        }
        std::cout << "max=" << maximum << std::endl;
    }
std::cout << "pass 6" << std::endl;

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)data , reader . o1 , reader . o2 , reader . o3 , num_x * (num_t) );

    return 0;

}


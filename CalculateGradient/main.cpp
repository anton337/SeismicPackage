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
        float r( out[k][0] );
        float i( out[k][1] );
        out[k][0] = in[k][0]*r - in[k][1]*i;
        out[k][1] = in[k][0]*i + in[k][1]*r;
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

    int width = get_argument ( arg , argc , argv );
    int inv = get_argument ( arg , argc , argv );
    int dir = get_argument ( arg , argc , argv );

    std::cout << "inv=" << inv << std::endl;

    std::cout << "dir=" << dir << std::endl;

    if ( inv )
    {
        for ( int k(0) ; k < num_x * num_t ; ++k )
        {
            data[k] = 1/(data[k]+1);
        }
    }

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

    float * filt = new float [ num_batch * num_t ];
    memset ( &filt[0] , 0 , num_batch * num_t );

    float total_val = 0 , val , r , fact ;
    for ( int fx(0) ; fx < width ; ++fx )
    {
        for ( int ft(0) ; ft < width ; ++ft )
        {
            r = sqrt ( fx*fx + ft*ft ) / (0.5*width);
            if ( dir == 1 ) // z dir
            {
                val = fx * exp ( - r * r );
            }
            else
            if ( dir == 2 ) // x dir
            {
                val = ft * exp ( - r * r );
            }
            fact = 1;
            filt[ft+fx*num_t] = fact * val;
            total_val += fact * val;
            std::cout << val << " \t";
            if ( num_t - ft < num_t )
            {
                fact = ((dir==2)?-1:1);
                filt[num_t-ft+fx*num_t] = fact * val;
                total_val += fact * val;
            }
            if ( num_batch - fx < num_batch )
            {
                fact = ((dir==1)?-1:1);
                filt[ft+(num_batch-fx)*num_t] = fact * val;
                total_val += fact * val;
            }
            if ( num_t - ft < num_t && num_batch - fx < num_batch )
            {
                fact = ((dir==2)?-1:1)*((dir==1)?-1:1);
                filt[num_t-ft+(num_batch-fx)*num_t] = fact * val;
                total_val += fact * val;
            }
        }
        std::cout << std::endl;
    }
    //for ( int fx(0) ; fx < width ; ++fx )
    //{
    //    for ( int ft(0) ; ft < width ; ++ft )
    //    {
    //        filt[ft+fx*num_t] /= total_val;
    //        if ( num_t - ft < num_t )
    //        {
    //            filt[num_t-ft+fx*num_t] /= total_val;
    //        }
    //        if ( num_batch - fx < num_batch )
    //        {
    //            filt[ft+(num_batch-fx)*num_t] /= total_val;
    //        }
    //        if ( num_t - ft < num_t && num_batch - fx < num_batch )
    //        {
    //            filt[num_t-ft+(num_batch-fx)*num_t] /= total_val;
    //        }
    //    }
    //}
    for ( int fx(-width) ; fx < width ; ++fx )
    {
        for ( int ft(-width) ; ft < width ; ++ft )
        {
            std::cout << filt[(num_t+ft)%num_t+((num_batch+fx)%num_batch)*num_t] << " \t";
        }
        std::cout << std::endl;
    }
    float total = 0;
    for ( int k(0) ; k < num_batch * num_t ; ++k )
    {
        total += filt[k];
    }
    std::cout << "total=" << total << std::endl;
    std::cout << "total_val=" << total_val << std::endl;

    int num_f = num_t / 2 + 1;
    fftwf_complex * fdata = new fftwf_complex [ num_x * num_f ];

    float * tmp = new float [ num_t * num_batch ];
    fftwf_complex * ftmp = new fftwf_complex [ num_f * num_batch ];

    fftwf_plan plan ( fftwf_plan_dft_r2c_2d ( num_batch
                                            , num_t
                                            , & tmp[0]
                                            , &ftmp[0]
                                            , FFTW_ESTIMATE
                                            )
                    );

    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        mem_cpy ( & tmp[0] , & data[x*num_t*num_batch] , num_t*num_batch );
        fftwf_execute ( plan );
        mem_cpy < float > ( (float*)(&fdata[x*num_f*num_batch]) , (float*)(&ftmp[0]) , 2*num_f*num_batch );
    }

    for ( int x(0)
        ; x < 1
        ; ++x
        )
    {
        mem_cpy ( & tmp[0] , & filt[0] , num_t*num_batch );
        fftwf_execute ( plan );
    }
    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        mem_mult ( (&fdata[x*num_f*num_batch]) , (&ftmp[0]) , num_f*num_batch );
    }

    fftwf_destroy_plan ( plan );

    fftwf_plan iplan ( fftwf_plan_dft_c2r_2d ( num_batch
                                             , num_t
                                             , &ftmp[0]
                                             , & tmp[0]
                                             , FFTW_ESTIMATE
                                             )
                     );

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

    for ( int x(0)
        ; x < num_batches
        ; ++x
        )
    {
        mem_cpy < float > ( (float*)(&ftmp[0]) , (float*)(&fdata[x*num_f*num_batch]) , 2*num_f*num_batch );
        fftwf_execute ( iplan );
        mem_cpy ( & data[x*num_t*num_batch] , & tmp[0] , num_t*num_batch );
        for ( int k(0)
            , size ( num_t * num_batch )
            ; k < size
            ; ++k
            )
        {
            data[x*size + k] /= num_t*num_batch;
            if ( inv )
            data[x*size + k] = 1/data[x*size + k]-1;
        }
    }

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


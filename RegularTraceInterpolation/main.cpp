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

    int interpolation_factor = get_argument ( arg , argc , argv );

    int num_iters = get_argument ( arg , argc , argv );

    int num_amp_iters = get_argument ( arg , argc , argv );

    int num_actual_amp_iters = get_argument ( arg , argc , argv );

    float slope_frac = get_argument ( arg , argc , argv );

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
    if ( num_t % 2 != 0 )
    {
        std::cout << "num_t should be even" << std::endl;
        exit(1);
    }
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    float * odata = new float [ interpolation_factor * num_x * num_t ];
    memset ( &odata[0] , 0 , interpolation_factor * num_x * num_t );

    int num_f = num_t / 2 + 1;

    float * tmp = new float [ interpolation_factor * num_t * num_batch ];
    fftwf_complex * ftmp = new fftwf_complex [ interpolation_factor * num_f * num_batch ];
    fftwf_complex * ftmp_orig = new fftwf_complex [ interpolation_factor * num_f * num_batch ];

    fftwf_plan plan ( fftwf_plan_dft_r2c_2d ( num_batch * interpolation_factor
                                            , num_t
                                            , & tmp[0]
                                            , &ftmp[0]
                                            , FFTW_ESTIMATE
                                            )
                    );

    fftwf_plan inv_plan ( fftwf_plan_dft_c2r_2d ( num_batch * interpolation_factor
                                                , num_t
                                                , &ftmp[0]
                                                , & tmp[0]
                                                , FFTW_ESTIMATE
                                                )
                        );

    float delta_frac ( 1.0 / num_iters );
    float delta_amp_frac ( 1.0 / num_amp_iters );
    float frac_amp = delta_amp_frac;
    float normalization_scalar = 1.0 / ( num_t * num_batch * interpolation_factor );
    for ( int amp_iter(0)
        ; amp_iter < num_actual_amp_iters
        ; ++amp_iter
        , frac_amp += delta_amp_frac
        )
    {
        float frac = delta_frac;
        for ( int iter(0)
            ; iter < num_iters
            ; ++iter
            , frac += delta_frac
            )
        {
            std::cout << "iter=" << iter << "   frac:" << frac << std::endl;
            for ( int x(0)
                ; x < num_batches
                ; ++x
                )
            {
                float max_amp ( 0 ) , tmp_amp , cur_amp ;
                mem_cpy ( & tmp[0] , & odata[x*num_t*num_batch*interpolation_factor] , num_t*num_batch*interpolation_factor );
                // reset live traces
                for ( int tr(0)
                    ; tr < num_batch
                    ; ++tr
                    )
                {
                    mem_cpy ( & tmp[tr*interpolation_factor*num_t] , & data[x*num_t*num_batch+tr*num_t] , num_t );
                }
                
                fftwf_execute ( plan );
                mem_cpy ( &((float*)ftmp_orig)[0] , &((float*)ftmp)[0] , num_t * num_batch * interpolation_factor );
                // threshold in FK domain
                int fk_tr(0);
                
                for ( int size_tr ( interpolation_factor*num_batch*(slope_frac*frac) )
                    ; fk_tr < size_tr
                    ; ++fk_tr
                    )
                {
                    int f(0);
                    for ( int size_f( (int)(frac * num_f) * 2 )
                        , size_f_max( num_f * 2 )
                        ; f < size_f && f < size_f_max
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f_max+f  ] *= normalization_scalar;
                        ((float*)ftmp)[fk_tr*size_f_max+f+1] *= normalization_scalar;
                        cur_amp = 
                        sqrtf ( 
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f  ] *
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f  ] 
                        +
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f+1] *
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f+1] 
                        );
                        for ( int inter(0)
                            ; inter < interpolation_factor
                            ; ++inter
                            )
                        {
                            tmp_amp = 
                            sqrtf ( 
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f  ] *
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f  ] 
                            +
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f+1] *
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f+1] 
                            );
                            if ( tmp_amp > max_amp )
                            {
                                max_amp = tmp_amp;
                            }
                        }
                        if ( cur_amp < (1-frac_amp) * max_amp )
                        {
                            ((float*)ftmp)[fk_tr*size_f_max+f  ] = 0;
                            ((float*)ftmp)[fk_tr*size_f_max+f+1] = 0;
                        }
                    }
                    for ( int size_f ( num_f * 2 )
                        ; f < size_f
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f+f  ] = 0;
                        ((float*)ftmp)[fk_tr*size_f+f+1] = 0;
                    }
                }
                
                for ( int size_tr ( interpolation_factor*num_batch*(1-(1-slope_frac)*frac) )
                    ; fk_tr < size_tr
                    ; ++fk_tr
                    )
                {
                    for ( int f(0)
                        , size_f ( num_f * 2 )
                        ; f < size_f
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f+f  ] = 0;
                        ((float*)ftmp)[fk_tr*size_f+f+1] = 0;
                    }
                }
                for ( int size_tr ( interpolation_factor*num_batch )
                    ; fk_tr < size_tr
                    ; ++fk_tr
                    )
                {
                    int f(0);
                    for ( int size_f( (int)(frac * num_f) * 2 )
                        , size_f_max( num_f * 2 )
                        ; f < size_f && f < size_f_max
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f_max+f  ] *= normalization_scalar;
                        ((float*)ftmp)[fk_tr*size_f_max+f+1] *= normalization_scalar;
                        cur_amp = 
                        sqrtf ( 
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f  ] *
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f  ] 
                        +
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f+1] *
                        ((float*)ftmp_orig)[fk_tr*size_f_max+f+1] 
                        );
                        for ( int inter(0)
                            ; inter < interpolation_factor
                            ; ++inter
                            )
                        {
                            tmp_amp = 
                            sqrtf ( 
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f  ] *
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f  ] 
                            +
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f+1] *
                            ((float*)ftmp_orig)[inter*num_batch*size_f_max+fk_tr*size_f_max+f+1] 
                            );
                            if ( tmp_amp > max_amp )
                            {
                                max_amp = tmp_amp;
                            }
                        }
                        if ( cur_amp < (1-frac_amp) * max_amp )
                        {
                            ((float*)ftmp)[fk_tr*size_f_max+f  ] = 0;
                            ((float*)ftmp)[fk_tr*size_f_max+f+1] = 0;
                        }
                    }
                    for ( int size_f ( num_f * 2 )
                        ; f < size_f
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f+f  ] = 0;
                        ((float*)ftmp)[fk_tr*size_f+f+1] = 0;
                    }
                    
                    for ( int size_f ( num_f * 2 )
                        ; f < size_f
                        ; f += 2
                        )
                    {
                        ((float*)ftmp)[fk_tr*size_f+f  ] = 0;
                        ((float*)ftmp)[fk_tr*size_f+f+1] = 0;
                    }
                }
                fftwf_execute ( inv_plan );
                
                mem_cpy ( & odata[x*num_t*num_batch*interpolation_factor] , & tmp[0] , num_t*num_batch*interpolation_factor );
            }
        }
    }

    fftwf_destroy_plan ( plan );

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 / interpolation_factor , interpolation_factor * reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)odata , reader . o1 , reader . o2 , reader . o3 , interpolation_factor * num_x * num_t );

    return 0;

}


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

    std::string output_file_name( argv[arg] );

    arg++;

    std::string output_data_file_name( argv[arg] );

    float prev_sample_rate = get_argument( arg , argc , argv );
    float curr_sample_rate = get_argument( arg , argc , argv );

    std::cout << "pass 1" << std::endl;

    SEPReader reader ( file_name . c_str () , false );
    int num_t = reader . n1;
    int new_num_t = reader . n1 * ( prev_sample_rate / curr_sample_rate ) + 1;
    if ( new_num_t % 2 == 1 )
    {
        new_num_t += 1;
    }
    int num_batch = reader . n2;
    int num_x = num_batch;
    int num_batches = reader . n3
                    * reader . n4
                    * reader . n5
                    * reader . n6
                    * reader . n7
                    * reader . n8
                    ;
    std::cout << num_t << std::endl;
    std::cout << num_x << std::endl;
    float * data = new float [ num_x * num_t * num_batches ];
    memset ( &data[0] , 0 , num_x * num_t * num_batches );

    reader . OpenDataFile ( (data_file_name) . c_str() );

    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * num_batches );

    float * out_data = new float [ num_x * new_num_t * num_batches ];
    memset ( &out_data[0] , 0 , num_x * new_num_t * num_batches );

    std::cout << "pass 2" << std::endl;

    int new_t;
    for ( int n(0)
        , k(0)
        ; n < num_batches
        ; ++n
        )
    for ( int x(0)
        ; x < num_x
        ; ++x
        )
    for ( int t(0)
        ; t < new_num_t
        ; ++t
        , ++k
        )
    {
        new_t = t*curr_sample_rate/prev_sample_rate;
        if ( new_t >= num_t )
        {
            new_t = num_t - 1;
        }
        out_data[k] = data[n*num_t*num_x+x*num_t+new_t];
    }

    std::cout << "pass 3" << std::endl;

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , new_num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_data_file_name) . c_str()
                     );

    writer . OpenDataFile ( (output_data_file_name) . c_str() );

    writer . write_sepval ( (float*)out_data , reader . o1 , reader . o2 , reader . o3 , num_x * (new_num_t) * num_batches );

    std::cout << "pass 4" << std::endl;

    return 0;

}


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

    std::cout << "pass 0" << std::endl;

    int pad_top    = get_argument( arg , argc , argv );
    std::cout << "pad top:" << pad_top << std::endl;
    int pad_bottom = get_argument( arg , argc , argv );
    std::cout << "pad bottom:" << pad_bottom << std::endl;
    int pad_left   = get_argument( arg , argc , argv );
    std::cout << "pad left:" << pad_left << std::endl;
    int pad_right  = get_argument( arg , argc , argv );
    std::cout << "pad right:" << pad_right << std::endl;

    std::cout << "pass 1" << std::endl;

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
    std::cout << num_t << std::endl;
    std::cout << num_x << std::endl;
    float * data = new float [ num_x * num_t * num_batches ];
    memset ( &data[0] , 0 , num_x * num_t * num_batches );

    reader . OpenDataFile ( (data_file_name) . c_str() );

    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * num_batches );

    std::cout << "pass 2" << std::endl;

    int new_num_x = num_x + pad_left + pad_right;

    int new_num_t = num_t + pad_top + pad_bottom;

    float * new_data = new float [ new_num_x * new_num_t * num_batches ];
    memset ( &new_data[0] , 0 , new_num_x * new_num_t * num_batches );

    for ( int n(0)
        , k(0)
        ; n < num_batches
        ; ++n
        )
    for ( int x(0) , X
        ; x < new_num_x
        ; ++x
        )
    for ( int t(0) , T
        ; t < new_num_t
        ; ++t
        , ++k
        )
    {
        if ( t < pad_top ) T = 0;
        else if ( t < pad_top + num_t ) T = t - pad_top;
        else if ( t < new_num_t ) T = num_t-1;
        if ( x < pad_left ) X = 0;
        else if ( x < pad_left + num_x ) X = x - pad_left;
        else if ( x < new_num_x ) X = num_x-1;
        new_data[k] = (T==0)?1500:data[T+X*num_t];
    }

    std::cout << "pass 3" << std::endl;

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , new_num_t
                     , reader . o2 , reader . d2 , new_num_x
                     , reader . o3 , reader . d3 , num_batches
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_data_file_name) . c_str()
                     );

    writer . OpenDataFile ( (output_data_file_name) . c_str() );

    writer . write_sepval ( (float*)new_data , reader . o1 , reader . o2 , reader . o3 , new_num_x * (new_num_t) * num_batches );

    std::cout << "pass 4" << std::endl;

    return 0;

}


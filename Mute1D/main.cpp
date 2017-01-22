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
    float x_pos;
    float y_pos;
    point ( float x , float y )
    : x_pos ( x )
    , y_pos ( y )
    {

    }
};

std::vector < point > selection;

bool mute_toggle ( std::vector < point > const & sel
                 , int x_pos 
                 , int y_pos
                 , int count_init
                 )
{
    float m;
    int count = count_init;
    for ( std::size_t k(0)
        ; k + 1 < sel . size ()
        ; ++k
        )
    {
        if ( x_pos > sel[k].x_pos && x_pos <= sel[k+1].x_pos )
        {
            m = (sel[k+1].y_pos - sel[k].y_pos) / (sel[k+1].x_pos - sel[k].x_pos);
            if ( y_pos < sel[k].y_pos + m * ( x_pos - sel[k].x_pos ) )
            {
                count += 1;
            }
        }
        if ( x_pos > sel[k+1].x_pos && x_pos <= sel[k].x_pos )
        {
            m = (sel[k+1].y_pos - sel[k].y_pos) / (sel[k+1].x_pos - sel[k].x_pos);
            if ( y_pos < sel[k].y_pos + m * ( x_pos - sel[k].x_pos ) )
            {
                count += 1;
            }
        }
    }
    return ( count % 2 == 0 );
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

    int   side      = get_argument( arg , argc , argv );
    while ( arg+1 < argc )
    {
        float x_pos     = get_argument( arg , argc , argv );
        float y_pos     = get_argument( arg , argc , argv );
        selection . push_back ( point ( x_pos , y_pos ) );
    }

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
        ; t < num_t
        ; ++t
        , ++k
        )
    {
        data[k] = (mute_toggle(selection,x,t,side>0))?data[k]:0;
    }

    std::cout << "pass 3" << std::endl;

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , reader . n2
                     , reader . o3 , reader . d3 , reader . n3
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_data_file_name) . c_str()
                     );

    writer . OpenDataFile ( (output_data_file_name) . c_str() );

    writer . write_sepval ( (float*)data , reader . o1 , reader . o2 , reader . o3 , num_x * (num_t) * num_batches );

    std::cout << "pass 4" << std::endl;

    return 0;

}


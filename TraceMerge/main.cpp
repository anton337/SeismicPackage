
#include <limits>
#include <iostream>
#include <string.h>
#include <math.h>
#include <map>
#include "sep_reader.h"
#include "sep_writer.h"

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

    if (argc < 3)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string output_file_name( argv[arg] );

    if ( argc == 3 )
    {
        arg++;
        std::string input_file_name( argv[arg] );
        SEPReader reader ( input_file_name . c_str () );
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
        std::cout << num_x << std::endl;
        float * data = new float [ num_x * num_t ];
        memset ( &data[0] , 0 , num_x * num_t );
        reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );
    }
    else
    if ( argc > 3 )
    {
        arg++;
        std::string input_file_name( argv[arg] );
        SEPReader reader ( input_file_name . c_str () , false );
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
        std::cout << num_x << std::endl;
        std::vector < std::string > header_labels;
        header_labels . push_back ( "D2" );
        header_labels . push_back ( "D3" );
        std::vector < std::string > sort_order;
        sort_order . push_back ( "D2" );
        sort_order . push_back ( "D3" );
        float * data = new float [ num_x * num_t ];
        memset ( &data[0] , 0 , num_x * num_t );
        float * hdrs = new float [ num_x * 2 ];
        memset ( &hdrs[0] , 0 , num_x * 2 );
        int num_layers = 0;
        for ( int k(3)
            , x(0)
            , i(0)
            ; k < argc
            ; ++k
            , ++x
            , ++num_layers
            )
        {
            arg++;
            std::cout << "x=" << x << "     ";
            reader . OpenDataFile ( argv[k] );
            reader . read_sepval ( & data [ x * num_t * num_batch ] , reader . o1 , reader . o2 , reader . o3 , num_batch * num_t );
            for ( int j(0)
                ; j < num_batch
                ; ++j
                , ++i
                )
            {
                hdrs[i*2  ] = x;
                hdrs[i*2+1] = j;
            }
        }
        SEPWriter hdr_writer ( (output_file_name + std::string(".hdrs")) . c_str()
                             , 0 , 1 , 2
                             , 0 , 1 , num_batch
                             , 0 , 1 , num_layers
                             , header_labels
                             , sort_order
                             , (output_file_name + std::string(".hdrs") + std::string("@")) . c_str()
                             );
        hdr_writer . OpenDataFile ( (output_file_name + std::string(".hdrs") + std::string("@")) . c_str() );
        hdr_writer . write_sepval ( &hdrs[0] , 0 , 0 , 0 , num_x*2 );
        SEPWriter writer ( (output_file_name) . c_str()
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , num_layers
                         , header_labels
                         , sort_order
                         , (output_file_name + std::string("@")) . c_str()
                         , (output_file_name + std::string(".hdrs")) . c_str()
                         );
        writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );
        writer . write_sepval ( &data[0] , 0 , 0 , 0 , num_x*num_t );
    }

    return 0;

}


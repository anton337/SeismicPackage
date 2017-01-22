
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

struct trace_hdrs
{
    std::vector < float > hdrs;
    int index;
};

std::ostream & operator << ( std::ostream & stream , trace_hdrs & tmp )
{
    stream << "ind=" << tmp.index << "   ";
    for ( std::size_t k(0)
        ; k < tmp.hdrs.size()
        ; ++k
        )
    {
        stream << tmp.hdrs[k] << " ";
    }
    std::cout << std::endl;
    return stream;
}

struct partition
{
    std::vector < float > hdrs;
};

std::ostream & operator << ( std::ostream & stream , partition & tmp )
{
    for ( std::size_t k(0)
        ; k < tmp.hdrs.size()
        ; ++k
        )
    {
        stream << tmp.hdrs[k] << " ";
    }
    std::cout << std::endl;
    return stream;
}

class MergeSort 
{

public:

    MergeSort ( int                                hdr_index
              , std::vector < trace_hdrs >       & x 
              )
    {
        std::vector < std::vector < trace_hdrs > > init;
        for ( std::size_t k(0)
            ; k < x . size ()
            ; ++k
            )
        {
            std::vector < trace_hdrs > a;
            a . push_back ( x[k] );
            init . push_back ( a );
        }
        while(1)
        {
            std::vector < std::vector < trace_hdrs > > next;
            for ( std::size_t k(0)
                ; k < init . size ()
                ; k += 2
                )
            {
                if ( k + 1 < init . size () )
                {
                    std::vector < trace_hdrs > b;
                    this -> operator () ( hdr_index
                                        , init[k] 
                                        , init[k+1]
                                        , b
                                        );
                    next . push_back ( b );
                }
                else
                {
                    next . push_back ( init[k] );
                }
            }
            init . clear ();
            init = next;
            next . clear ();
            if ( init . size () == 1 ) 
            {
                x = init[0];
                break;
            }
        }
    }

    void operator () ( int                                hdr_index
                     , std::vector < trace_hdrs > const & in1
                     , std::vector < trace_hdrs > const & in2
                     , std::vector < trace_hdrs >       & out
                     )
    {
        out . clear ();
        std::size_t index1 = 0;
        std::size_t index2 = 0;
        {
            if ( in1[index1].hdrs[hdr_index] == in2[index2].hdrs[hdr_index] )
            {
                out . push_back ( in1[index1] );
                ++index1;
            }
            else
            if ( in1[index1].hdrs[hdr_index] < in2[index2].hdrs[hdr_index] )
            {
                out . push_back ( in1[index1] );
                ++index1;
            }
            else
            {
                out . push_back ( in2[index2] );
                ++index2;
            }
        }
        while ( 1 )
        {
            if ( index1 >= in1.size() && index2 >= in2.size() )
            {
                break;
            }
            if ( index1 < in1.size () && index2 < in2.size () )
            {
                if ( in1[index1].hdrs[hdr_index] == in2[index2].hdrs[hdr_index] )
                {
                    out . push_back ( in1[index1] );
                    ++index1;
                }
                else
                if ( in1[index1].hdrs[hdr_index] < in2[index2].hdrs[hdr_index] )
                {
                    out . push_back ( in1[index1] );
                    ++index1;
                }
                else
                {
                    out . push_back ( in2[index2] );
                    ++index2;
                }
            }
            else
            if ( index1 < in1 . size() )
            {
                out . push_back ( in1[index1] );
                ++index1;
            }
            else
            {
                out . push_back ( in2[index2] );
                ++index2;
            }
        }
    }

};


void get_extends ( int                                hdr_index 
                 , float                            & min_elem
                 , float                            & max_elem
                 , float                            & delta
                 , int                              & num_unique
                 , std::vector < trace_hdrs > const & in
                 )
{
    int _min_elem ( std::numeric_limits < int > :: max () );
    int _max_elem ( std::numeric_limits < int > :: min () );
    int _delta = 1;
    int _num = 0;
    int _curr = 0;
    for ( std::size_t k(0)
        ; k < in . size ()
        ; ++k
        )
    {
        _curr = in [ k ] . hdrs [ hdr_index ];
        if ( _curr > _max_elem )
        {
            _delta = _curr - _max_elem;
            _max_elem = _curr;
        }
        if ( _curr < _min_elem )
        {
            _min_elem = _curr;
        }
    }
    _num = 1 + (_max_elem - _min_elem) / _delta;
    max_elem = _max_elem;
    min_elem = _min_elem;
    delta = _delta;
    num_unique = _num;
}


struct extends
{
    float min_elem;
    float max_elem;
    float delta;
    int num;
};


struct partition_comparator
{
    bool operator () ( partition const & a 
                     , partition const & b
                     )
    {
        if ( a . hdrs . size () != b . hdrs . size () )
        {
            std::cout << "comparing wrong things." << std::endl;
            exit(1);
        }
        for ( std::size_t k ( 0 )
            ; k < a . hdrs . size ()
            ; ++k
            )
        {
            if ( (int)a . hdrs [k] < (int)b . hdrs [k] )
            {
                return true;
            }
        }
        return false;
    }
};


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

    SEPReader reader ( file_name.c_str() );

    std::string header_name( reader . get_hdrs_file_name () );

    if ( header_name . size () == 0 ) // empty string
    {
        std::cout << "couldn't find headers file" << std::endl;
        exit(1);
    }

    SEPReader hdr_reader ( header_name.c_str() );

    std::vector < trace_hdrs > traces;

    int num_x = hdr_reader . n2 
              * hdr_reader . n3 
              * hdr_reader . n4 
              * hdr_reader . n5 
              * hdr_reader . n6 
              * hdr_reader . n7 
              * hdr_reader . n8 
              ;

    int num_hdrs = hdr_reader . n1;

    float * hdr_array = new float [ num_hdrs * num_x ];



    hdr_reader . read_sepval ( &hdr_array[0]
                             , hdr_reader . o1 
                             , hdr_reader . o2 
                             , hdr_reader . o3
                             , num_hdrs * num_x
                             );


    for ( int k(0)
        , t(0)
        ; k < num_x
        ; ++k
        )
    {
        trace_hdrs trace;
        trace . index = k;
        for ( int j(0)
            ; j < num_hdrs
            ; ++j
            , ++t
            )
        {
            trace . hdrs . push_back ( hdr_array [ t ] );
        }
        traces . push_back ( trace );
    }

    std::vector < std::string > const header_labels ( hdr_reader . get_header_labels () );

    std::map < std::string , int > header_index;

    for ( std::size_t k(0)
        ; k < header_labels . size ()
        ; ++k
        )
    {
        header_index [ header_labels [ k ] ] = k;
    }

    std::vector < int > partition_keys;
    std::vector < std::string > partition_keys_wrd;

    int num_partition_keys = get_argument ( arg , argc , argv );

    for ( int k(0)
        ; k < num_partition_keys
        ; ++k
        )
    {
        arg++;
        partition_keys . push_back ( header_index [ std::string ( argv[ arg ] ) ] );
        partition_keys_wrd . push_back ( std::string ( argv[ arg ] ) );
        std::cout << "partition key : " << argv[ arg ] << std::endl;
    }

    std::vector < int > requested_sort_order;
    std::vector < std::string > requested_sort_order_wrd;

    int num_sort_keys = get_argument ( arg , argc , argv );

    for ( int k(0)
        ; k < num_sort_keys
        ; ++k
        )
    {
        arg++;
        requested_sort_order . push_back ( header_index [ std::string ( argv[ arg ] ) ] );
        requested_sort_order_wrd . push_back ( std::string ( argv[ arg ] ) );
        std::cout << "sort key : " << argv[ arg ] << std::endl;
    }

    std::vector < extends > extnd ( num_sort_keys );

    for ( int k(0)
        ; k < num_sort_keys
        ; ++k
        )
    {
        int hdr_index = requested_sort_order[k];
        get_extends ( hdr_index 
                    , extnd[k] . min_elem
                    , extnd[k] . max_elem
                    , extnd[k] . delta
                    , extnd[k] . num
                    , traces
                    );
    }

    std::map < partition , std::vector < trace_hdrs > , partition_comparator > index_vector;
    std::map < partition , std::vector < trace_hdrs > , partition_comparator > :: iterator it;

    for ( std::size_t k(0)
        ; k < traces . size ()
        ; ++k
        )
    {
        partition p;
        for ( std::size_t j(0)
            ; j < partition_keys . size ()
            ; ++j
            )
        {
            p . hdrs . push_back ( traces[k] . hdrs [ partition_keys [ j ] ] );
        }
        it = index_vector . find ( p );
        if ( it == index_vector . end () )
        {
            std::vector < trace_hdrs > vec;
            vec . push_back ( traces[k] );
            index_vector . insert ( std :: pair < partition , std::vector < trace_hdrs > > ( p , vec ) );
        }
        else
        {
            it -> second . push_back ( traces[k] );
        }
    }

    std::map < partition , std::vector < trace_hdrs > , partition_comparator > orig_index_vector = index_vector;
    std::map < partition , std::vector < trace_hdrs > , partition_comparator > :: iterator orig_it;

    float * trace_data = new float [ reader . n1 ];
         it =      index_vector . begin ();
    orig_it = orig_index_vector . begin ();
    for ( int count
        ;       it !=      index_vector . end ()
        && orig_it != orig_index_vector . end ()
        ; ++it
        , ++orig_it
        , ++count
        )
    {

        //for ( std::size_t k(0)
        //    , size ( requested_sort_order . size () ) 
        //    ; k < size
        //    ; ++k
        //    )
        //{
        //    MergeSort ( requested_sort_order [ (size-1) - k ]
        //              , it -> second 
        //              );
        //}

        SEPWriter hdr_writer ( (output_file_name + std::string(".hdrs")) . c_str()
                             , 0 , 1 , num_hdrs
                             , extnd[1] . min_elem , extnd[1] . delta , extnd[1] . num
                             , extnd[0] . min_elem , extnd[0] . delta , extnd[0] . num
                             , header_labels
                             , requested_sort_order_wrd
                             , (output_file_name + std::string(".hdrs") + std::string("@")) . c_str()
                             );

        hdr_writer . OpenDataFile ( (output_file_name + std::string(".hdrs") + std::string("@")) . c_str() );

        for ( std::size_t j(0)
            ; j < it -> second . size ()
            ; ++j
            )
        {
            hdr_writer . write_sepval ( &(it -> second[j] . hdrs[0]) 
                                      , 0 
                                      , it -> second[j] . hdrs [ requested_sort_order[1] ] 
                                      , it -> second[j] . hdrs [ requested_sort_order[0] ]
                                      , num_hdrs 
                                      );
        }

        SEPWriter writer ( (output_file_name) . c_str()
                         , reader . o1 , reader . d1 , reader . n1
                         , extnd[1] . min_elem , extnd[1] . delta , extnd[1] . num
                         , extnd[0] . min_elem , extnd[0] . delta , extnd[0] . num
                         , header_labels
                         , requested_sort_order_wrd
                         , (output_file_name + std::string("@")) . c_str()
                         , (output_file_name + std::string(".hdrs")) . c_str()
                         );

        writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

        reader . OpenDataFile ( (file_name + std::string("@")).c_str() );

        for ( std::size_t j(0)
            ; j < it -> second . size ()
            ; ++j
            )
        {
            reader .  read_sepval ( &trace_data[0]
                                  , reader . o1
                                  , reader . o2 + reader . d2 * (j % reader . n2)
                                  , reader . o3 + reader . d3 * ((j / reader . n2) % reader . n3)
                                  , reader . n1
                                  );
            writer . write_sepval ( &trace_data[0] 
                                  , reader . o1 
                                  , it -> second[j] . hdrs [ requested_sort_order[1] ] 
                                  , it -> second[j] . hdrs [ requested_sort_order[0] ]
                                  , reader . n1 
                                  );
        }

    }

    return 0;

}


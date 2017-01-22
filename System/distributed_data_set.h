#ifndef DISTRIBUTED_DATA_SET_H
#define DISTRIBUTED_DATA_SET_H

#include <vector>
#include <list>
#include <string>
#include <sstream>

#include "sep_reader.h"
#include "sep_writer.h"

template < typename type , long N=3 >
class DistributedDataSet
{

    type * buffer;

    std::vector < std::pair < long , std::string > > file_names;

    long latest_writer_z;

    long latest_writer_index;

    SEPWriter * latest_writer;

    long latest_reader_z;

    long latest_reader_index;

    SEPReader * latest_reader;

    std::list < std::pair < long , SEPWriter * > > list_writer;

    std::list < std::pair < long , SEPReader * > > list_reader;

    int num_tiles;

    std::string file_name_prefix;

public:

    long size[N];
    long tile[N];
    long ind[N];
    long pad[N];

    DistributedDataSet ( long * _size 
                       , long * _tile
                       , long * _pad
                       , std::string _file_name_prefix
                       , bool read_only = true
                       )
    {
        std::string file_name_prefix = _file_name_prefix;
        for ( int k(0)
            ; k < N
            ; ++k
            )
        {
            size[k] = _size[k];
            tile[k] = _tile[k];
            pad [k] = _pad [k];
            ind [k] =  size[k] / tile[k];
        }
        latest_writer = NULL;
        latest_reader = NULL;
        latest_writer_z = -1;
        latest_writer_index = -1;
        latest_reader_z = -1;
        latest_reader_index = -1;
        buffer = new type [ size[2] ];
        num_tiles = 1;
        for ( int k(0)
            ; k < N
            ; ++k
            )
        {
            num_tiles *= ( size[k] / tile[k] );
        }
        for ( long x(0)
            , num_x(size[0]/tile[0])
            , k(0)
            ; x < num_x
            ; ++x
            )
        {
            for ( long y(0)
                , num_y(size[1]/tile[1])
                ; y < num_y
                ; ++y
                )
            {
                for ( long z(0)
                    , num_z(size[2]/tile[2])
                    ; z < num_z
                    ; ++z
                    , ++k
                    )
                {
                    std::vector < std::string > header_labels;
                    std::vector < std::string > sort_order;
                    std::stringstream ss;
                    ss << file_name_prefix << "-" << k << ".SEP";
                    file_names . push_back ( std::pair < long , std::string > ( k 
                                                                              , std::string ( ss.str() )
                                                                              )
                                           );
                }
            }
        }
        // create files
        if ( read_only != true )
        {
            long total_size = (tile[0]+pad[0]) * (tile[1]+pad[1]) * (tile[2]+pad[1]);
            type * data = new type [ total_size ];
            for ( long x(0)
                , num_x(size[0]/tile[0])
                , k(0)
                ; x < num_x
                ; ++x
                )
            {
                for ( long y(0)
                    , num_y(size[1]/tile[1])
                    ; y < num_y
                    ; ++y
                    )
                {
                    for ( long z(0)
                        , num_z(size[2]/tile[2])
                        ; z < num_z
                        ; ++z
                        , ++k
                        )
                    {
                        std::vector < std::string > header_labels;
                        std::vector < std::string > sort_order;
                        std::stringstream ss;
                        ss << file_name_prefix << "-" << k << ".SEP";
                        SEPWriter writer ( ss.str().c_str()
                                         , tile[2]*z , 1 , tile[2]+pad[2]
                                         , tile[1]*y , 1 , tile[1]+pad[1]
                                         , tile[0]*x , 1 , tile[0]+pad[0]
                                         , header_labels
                                         , sort_order
                                         , std::string ( ss.str() + "@" ) . c_str()
                                         , std::string ( ss.str() + ".hdrs" ) . c_str()
                                         );
                        writer . OpenDataFile ( std::string ( ss.str() + "@" ) . c_str() );
                        writer . write_sepval ( &data[0] 
                                              , tile[2]*z 
                                              , tile[1]*y 
                                              , tile[0]*x 
                                              , total_size
                                              );
                        writer . CloseDataFile ();
                    }
                }
            }
            delete [] data;
        }
    }

    ~DistributedDataSet ( )
    {

    }

    inline 
    void set_tile ( long k
                  , long n
                  , type * data
                  )
    {
        long z = k%ind[2];
        long y = (k/ind[2])%ind[1];
        long x = (k/(ind[2]*ind[1]))%ind[0];
        std::vector < std::string > header_labels;
        std::vector < std::string > sort_order;
        std::stringstream ss;
        ss << file_names[k].second;
        SEPWriter writer ( ss.str().c_str()
                         , tile[2]*z , 1 , tile[2]+pad[2]
                         , tile[1]*y , 1 , tile[1]+pad[1]
                         , tile[0]*x , 1 , tile[0]+pad[0]
                         , header_labels
                         , sort_order
                         , std::string ( ss.str() + "@" ) . c_str()
                         , std::string ( ss.str() + ".hdrs" ) . c_str()
                         );
        writer . OpenDataFile ( std::string ( ss.str() + "@" ) . c_str() );
        writer . write_sepval ( &data[0] 
                              , tile[2]*z 
                              , tile[1]*y 
                              , tile[0]*x 
                              , n
                              );
        writer . CloseDataFile ();
    }

    inline
    void set ( long x
             , long y
             , long z
             , type val
             )
    {
        long index = (z/tile[2]) + ind[1]*((y/tile[1]) + ind[0]*(x/tile[0]));
        if ( index == latest_writer_index && z%tile[2] == latest_writer_z+1 )
        {
            buffer[z%tile[2]] = val;
        }
        else
        {
            if ( index == latest_writer_index )
            {
                latest_writer -> write_sepval ( &buffer[0]
                                              , latest_writer -> o1 + 0
                                              , latest_writer -> o2 + (y)%latest_writer->n2
                                              , latest_writer -> o3 + (x)%latest_writer->n3
                                              , latest_writer -> n1
                                              );
                buffer[z%tile[2]] = val;
            }
            else
            {
                if ( latest_writer == NULL )
                {
                    std::vector < std::string > header_labels;
                    std::vector < std::string > sort_order;
                    latest_writer_index = index;
                    latest_writer = new SEPWriter ( file_names[index].second.c_str()
                                                  , tile[2]*(z/tile[2]) , 1 , tile[2]
                                                  , tile[1]*(y/tile[1]) , 1 , tile[1]
                                                  , tile[0]*(x/tile[0]) , 1 , tile[0]
                                                  , header_labels
                                                  , sort_order
                                                  , std::string ( file_names[index].second + "@" ) . c_str()
                                                  , std::string ( file_names[index].second + ".hdrs" ) . c_str()
                                                  );
                    latest_writer -> OpenDataFile ( std::string ( file_names[index].second + "@" ) . c_str() );
                }
                else
                {
                    bool found = false;
                    std::list < std::pair < long , SEPWriter * > > :: iterator it = list_writer . begin ();
                    while ( it != list_writer . end () )
                    {
                        if ( it->first == latest_writer_index )
                        {
                            found = true;
                            break;
                        }
                        ++it;
                    }
                    if ( !found )
                    {
                        list_writer . push_front ( std::pair < long , SEPWriter * > ( latest_writer_index , latest_writer ) );
                        if ( list_writer . size () > ind[1]*ind[2] )
                        {
                            list_writer . back() . second -> CloseDataFile ();
                            list_writer . pop_back ();
                        }
                    }
                    found = false;
                    it = list_writer . begin ();
                    while ( it != list_writer . end () )
                    {
                        if ( it->first == index )
                        {
                            found = true;
                            latest_writer = it->second;
                            latest_writer_index = index;
                            break;
                        }
                        ++it;
                    }
                    if ( !found )
                    {
                        std::vector < std::string > header_labels;
                        std::vector < std::string > sort_order;
                        latest_writer_index = index;
                        latest_writer = new SEPWriter ( file_names[index].second.c_str()
                                                      , tile[2]*(z/tile[2]) , 1 , tile[2]
                                                      , tile[1]*(y/tile[1]) , 1 , tile[1]
                                                      , tile[0]*(x/tile[0]) , 1 , tile[0]
                                                      , header_labels
                                                      , sort_order
                                                      , std::string ( file_names[index].second + "@" ) . c_str()
                                                      , std::string ( file_names[index].second + ".hdrs" ) . c_str()
                                                      );
                        latest_writer -> OpenDataFile ( std::string ( file_names[index].second + "@" ) . c_str() );
                    }
                }
                buffer[z%tile[2]] = val;
            }
        }
    }

    type read()
    {

    }

    type * read_tile()
    {

    }

};

#endif


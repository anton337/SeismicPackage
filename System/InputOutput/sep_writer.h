#ifndef SEP_WRITER_H
#define SEP_WRITER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>

class SEPWriter
{

    FILE * sep_DATA;
    FILE * sep_head;

public:

    int o1 , d1 , n1;
    int o2 , d2 , n2;
    int o3 , d3 , n3;
    int o4 , d4 , n4;
    int o5 , d5 , n5;
    int o6 , d6 , n6;
    int o7 , d7 , n7;
    int o8 , d8 , n8;

    std::string hdrs_name;

    std::vector < std::string > header_labels;

    std::vector < std::string > sort_order;

    SEPWriter ( const char * sepname 
              , int _o1 , int _d1 , int _n1
              , int _o2 , int _d2 , int _n2
              , int _o3 , int _d3 , int _n3
              , std::vector < std::string > const & _header_labels
              , std::vector < std::string > const & _sort_order
              , const char * dataname = NULL
              , const char * hdrsname = NULL
              )
    : o1 ( _o1 ) , d1 ( _d1 ) , n1 ( _n1 )
    , o2 ( _o2 ) , d2 ( _d3 ) , n2 ( _n2 )
    , o3 ( _o3 ) , d3 ( _d2 ) , n3 ( _n3 )
    , header_labels ( _header_labels )
    , sort_order    ( _sort_order    )
    {
        if ((sep_head = fopen(sepname,"w")) == NULL) {
          std::cout << "Error opening SEP header: " << sepname << std::endl;
          exit(1);
        }  
        WriteHeaderLabels ( header_labels );
        WriteSortOrder ( sort_order );
        SetHeaderInfo ( sepname 
                      , dataname 
                      , hdrsname
                      );
    }

    void WriteHeaderLabels ( std::vector < std::string > const & header_labels )
    {
        for ( std::size_t k(0)
            ; k < header_labels . size ()
            ; ++k
            )
        {
            fprintf ( sep_head , "hdr_label=" );
            fprintf ( sep_head , "%s\n" , header_labels[k].c_str() );
        }
    }

    void WriteSortOrder ( std::vector < std::string > const & sort_order )
    {
        for ( std::size_t k(0)
            ; k < sort_order . size ()
            ; ++k
            )
        {
            fprintf ( sep_head , "sort_order=" );
            fprintf ( sep_head , "%s\n" , sort_order[k].c_str() );
        }
    }

    void SetHeaderInfo ( const char * sepname
                       , const char * dataname 
                       , const char * hdrsname 
                       )
    {
        fprintf ( sep_head , "o1=%d\n" , o1 );
        fprintf ( sep_head , "d1=%d\n" , d1 );
        fprintf ( sep_head , "n1=%d\n" , n1 );
        fprintf ( sep_head , "o2=%d\n" , o2 );
        fprintf ( sep_head , "d2=%d\n" , d2 );
        fprintf ( sep_head , "n2=%d\n" , n2 );
        fprintf ( sep_head , "o3=%d\n" , o3 );
        fprintf ( sep_head , "d3=%d\n" , d3 );
        fprintf ( sep_head , "n3=%d\n" , n3 );
        fprintf ( sep_head , "data_format=xdr_float\n" );
        fprintf ( sep_head , "little_endian=yes\n" );
        char sep_dname[1024];
        if ( dataname )
        {
            sprintf(sep_dname, "%s",dataname);
        }
        else
        {
            sprintf(sep_dname, "%s@", sepname );
        }
        fprintf ( sep_head , "in=%s\n" , sep_dname );
        if ( hdrsname )
        {
            hdrs_name = std::string ( hdrsname );
            fprintf ( sep_head , "hdrs=%s\n" , hdrsname );
        }
        fclose(sep_head);
        if ((sep_head = fopen(sepname,"r")) == NULL) {
          std::cout << "Error opening " << sepname << std::endl;
        }
        if ((sep_DATA = fopen(sep_dname,"w")) == NULL) {
          std::cout << "Error opening " << sep_dname << std::endl;
        }
        setvbuf(sep_DATA, NULL, _IONBF, 0);
        fclose(sep_DATA);
        fclose(sep_head);
    }

    void OpenDataFile ( char const * in )
    {
        std::cout << "Attempting to open data file: \"" << in << "\"" << std::endl;
        if ( ( sep_DATA = fopen ( in , "rw+" ) ) == NULL )
        {
            std::cout << "Error opening data file: \"" << in << "\"" << std::endl;
        }
    }
    
    void CloseDataFile ()
    {
        fclose(sep_DATA);
    }

    void _impl_byte_swap_( float * p_arg ) const
    {
        union { char a[ sizeof( float ) ]; float b; } temp; temp.b = *p_arg;
        for( unsigned d ( 0 ), end( sizeof(float)-1 ); d < (sizeof(float)>>1); ++d, --end ) {
            temp.a[   d ] ^= temp.a[ end ];
            temp.a[ end ] ^= temp.a[   d ];
            temp.a[   d ] ^= temp.a[ end ];
        }
    
        *p_arg = temp.b;
    }

    int write_sepval(float *val, int x,int y,int o,int n,bool byteswap = false) const 
    {
      off64_t foffset, o64, x64, y64, n164, n264;
      int status;
    
      status=1;
      if (x<o1 || x>=o1+n1) {
        status=0;
      }
      if (y<o2 || y>=o2+n2) {
        status=0;
      }
      if (o<o3 || o>=o3+n3) {
        status=0;
      }
      if (status == 1) {
        x64=(off64_t) x - (off64_t) o1;
        y64=(off64_t) y - (off64_t) o2;
        o64=(off64_t) o - (off64_t) o3;
        n164=n1;
        n264=n2;
        foffset = (o64*n164*n264 + y64*n164 + x64) * sizeof(float);
        fseek(sep_DATA, foffset, SEEK_SET);
        if(fwrite(val, sizeof(float), n, sep_DATA) == 0) {
          std::cout << "Error writing sep value at x=" << x << " y=" << y << " o=" << o << std::endl;
          return 0;
        }
      }
      else {
        // error do nothing
        std::cout << "write index out of range" << std::endl;
        std::cout << "x=" << x << " " << o1 << " " << o1+n1 << std::endl;
        std::cout << "y=" << y << " " << o2 << " " << o2+n2 << std::endl;
        std::cout << "z=" << o << " " << o3 << " " << o3+n3 << std::endl;
      }
      return 1;
    }

};

#endif


#ifndef SEP_READER_H
#define SEP_READER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <string>

class SEPReader
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

    SEPReader ( const char * sepname , bool open_data_file = true )
    : o1 ( 0 ) , d1 ( 0 ) , n1 ( 1 )
    , o2 ( 0 ) , d2 ( 0 ) , n2 ( 1 )
    , o3 ( 0 ) , d3 ( 0 ) , n3 ( 1 )
    , o4 ( 0 ) , d4 ( 0 ) , n4 ( 1 )
    , o5 ( 0 ) , d5 ( 0 ) , n5 ( 1 )
    , o6 ( 0 ) , d6 ( 0 ) , n6 ( 1 )
    , o7 ( 0 ) , d7 ( 0 ) , n7 ( 1 )
    , o8 ( 0 ) , d8 ( 0 ) , n8 ( 1 )
    {
        if ((sep_head = fopen(sepname,"r")) == NULL) {
          std::cout << "Error opening SEP header: " << sepname << std::endl;
          exit(1);
        }  
        GetHeaderInfo ( open_data_file );
    }

    std::string const & get_hdrs_file_name ()
    {
        return hdrs_name;
    }

    std::vector < std::string > const & get_header_labels ()
    {
        return header_labels;
    }

    std::vector < std::string > const & get_sort_order ()
    {
        return sort_order;
    }

    void check_duplicates ( std::vector < std::string > & vec , const char * in )
    {
        for ( std::size_t k(0)
            ; k < vec . size ()
            ; ++k
            )
        {
            if ( strncmp ( vec[k].c_str() , in , vec[k].size() ) == 0 )
            {
                std::cout << "duplicate header label, check SEP file." << std::endl;
                exit(1);
            }
        }
    }

    void GetHeaderInfo ( bool open_data_file = true )
    {
        header_labels . clear ();
        sort_order    . clear ();
        char junk1[1024];
        char *junk = junk1;
        char *junk2;
        char in1[1024];
        char *in=in1;
        while (!feof(sep_head)) {
          int ret = fscanf(sep_head, "%s", junk);
          if ( !ret )
          {
            std::cout << "fscan returns " << ret << std::endl;
            exit(1);
          }
          if (strncmp(junk, "o1=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &o1); }
          if (strncmp(junk, "d1=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &d1); }
          if (strncmp(junk, "n1=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &n1); }
          if (strncmp(junk, "o2=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &o2); }
          if (strncmp(junk, "d2=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &d2); }
          if (strncmp(junk, "n2=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &n2); }
          if (strncmp(junk, "o3=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &o3); }
          if (strncmp(junk, "d3=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &d3); }
          if (strncmp(junk, "n3=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &n3); }
          if (strncmp(junk, "o4=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &o4); }
          if (strncmp(junk, "d4=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &d4); }
          if (strncmp(junk, "n4=",3)==0) { junk2 = junk + 3; sscanf(junk2, "%d", &n4); }
          if (strncmp(junk, "hdr_label=",10)==0) { junk2 = junk + 10; sscanf(junk2, "%s", in); check_duplicates (header_labels , in); header_labels . push_back ( std::string ( in ) ); }
          if (strncmp(junk, "sort_order=",11)==0) { junk2 = junk + 11; sscanf(junk2, "%s", in); check_duplicates (sort_order , in); sort_order . push_back ( std::string ( in ) ); }
          if (strncmp(junk, "hdrs=",5)==0) { 
            junk2 = junk + 5; 
            sscanf(junk2, "%s", in); 
            std::cout << "Headers File: " << in << std::endl;
            hdrs_name = std::string ( in );
          }
          if ( open_data_file )
          {
            if (strncmp(junk, "in=\"",4)==0) { 
              junk2 = junk + 4; 
              strtok(junk2, "\"");
              sscanf(junk2, "%s", in); 
              std::cout << "Attempting to open data file: " << in << std::endl;
              if ( ( sep_DATA = fopen ( in , "rw+" ) ) == NULL )
              {
                  std::cout << "Error opening data file: " << in << std::endl;
              }
            }
            else
            if (strncmp(junk, "in=",3)==0) { 
              junk2 = junk + 3; 
              sscanf(junk2, "%s", in); 
              std::cout << "Attempting to open data file: " << in << std::endl;
              if ( ( sep_DATA = fopen ( in , "rw+" ) ) == NULL )
              {
                  std::cout << "Error opening data file: " << in << std::endl;
              }
            }
          }
        }
        std::cout << n1 << " " << n2 << " " << n3 << std::endl;
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

    int read_sepval(float *val, int x,int y,int o,int n,bool byteswap = false) const 
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
        if(fread(val, sizeof(float), n, sep_DATA) == 0) {
          std::cout << "Error reading sep value at x=" << x << " y=" << y << " o=" << o << std::endl;
          return 0;
        }
        if (byteswap)
        {
          for ( int k(0); k < n; ++k )
          {
            _impl_byte_swap_(&val[k]);
          }
        }
      }
      else {
        // error do nothing
        std::cout << "read index out of range" << std::endl;
        std::cout << "x=" << x << " " << o1 << " " << o1+n1 << std::endl;
        std::cout << "y=" << y << " " << o2 << " " << o2+n2 << std::endl;
        std::cout << "z=" << o << " " << o3 << " " << o3+n3 << std::endl;
      }
      return 1;
    }

};

#endif


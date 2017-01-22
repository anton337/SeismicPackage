#ifndef SEGY_READER_H
#define SEGY_READER_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>

class SEGYReader
{

    FILE * segy_file;

    unsigned short num_x;
    unsigned short num_t;
    unsigned short format;
    unsigned short sampr;


public:

    SEGYReader ( const char * segyname )
    {
        if ((segy_file = fopen(segyname,"rb")) == NULL) {
          std::cout << "Error opening SEGY file: " << segyname << std::endl;
          exit(1);
        }

        off64_t foffset;

        unsigned char binary_header [400];
        foffset = 3200;
        fseek ( segy_file , foffset , SEEK_SET );
        if ( fread ( &binary_header , sizeof(unsigned char) , 400 , segy_file ) == 0 )
        {
            std::cout << "error getting binary header" << std::endl;
            return;
        }
        num_x = * reinterpret_cast < unsigned short * > ( & binary_header [ 13 ] );
        std::cout << "good num_x = " << num_x << std::endl;
        // num_x = * reinterpret_cast < unsigned short * > ( & binary_header [ 15 ] );
        // std::cout << "num_x = " << num_x << std::endl;
        sampr = * reinterpret_cast < unsigned short * > ( & binary_header [ 17 ] );
        std::cout << "sampr = " << sampr << std::endl;
        num_t = * reinterpret_cast < unsigned short * > ( & binary_header [ 21 ] );
        std::cout << "num_t = " << num_t << std::endl;
        format = * reinterpret_cast < unsigned short * > ( & binary_header [ 25 ] );
        std::cout << "format = " << format << std::endl;

        // unsigned short x = 258;
        // _impl_byte_swap_ ( &x );
        // std::cout << x << std::endl;

        unsigned char trace_header [240];
        foffset = 3600;
        fseek ( segy_file , foffset , SEEK_SET );
        if ( fread ( &trace_header , sizeof(unsigned char) , 240 , segy_file ) == 0 )
        {
            std::cout << "error getting binary header" << std::endl;
            return;
        }
        unsigned int index;
        index = * reinterpret_cast < unsigned int * > ( & trace_header [ 0 ] );
        _impl_byte_swap_ ( &index );
        std::cout << "index = " << index << std::endl;
        sampr = * reinterpret_cast < unsigned short * > ( & trace_header [ 114 ] );
        _impl_byte_swap_ ( &sampr );
        std::cout << "sampr = " << sampr << std::endl;
        num_t = * reinterpret_cast < unsigned short * > ( & trace_header [ 116 ] );
        _impl_byte_swap_ ( &num_t );
        std::cout << "num_t = " << num_t << std::endl;


        unsigned char trace_data [100000];
        foffset = 3600+240;
        fseek ( segy_file , foffset , SEEK_SET );
        if ( fread ( &trace_data , sizeof(unsigned char) , 100000 , segy_file ) == 0 )
        {
            std::cout << "error getting binary header" << std::endl;
            return;
        }
        index = * reinterpret_cast < unsigned int * > ( & trace_data [ 207*4 ] );
        _impl_byte_swap_ ( &index );
        std::cout << "index = " << index << std::endl;

        int _k = 0;
        while ( 1 )
        {
            for ( int k(0)
                ; k < 50
                ; ++k
                , ++_k
                )
            {
                float val = * reinterpret_cast < float * > ( & trace_data[_k*4] );
                _impl_byte_swap_ ( &val );
                std::cout << "k=" << _k << " \"" << val << "\"" << std::endl;
            }
            char ch;
            std::cin >> ch;
        }

    }

    void _impl_byte_swap_( unsigned short * x ) const
    {
        unsigned short y = *x;
        unsigned short z = (y << 8) | (y >> 8);
        *x = z;
    }

    void _impl_byte_swap_( unsigned int * x ) const
    {
        unsigned int num = *x;
        unsigned int swapped = ((num>>24)&0xff) | // move byte 3 to byte 0
                               ((num<<8)&0xff0000) | // move byte 1 to byte 2
                               ((num>>8)&0xff00) | // move byte 2 to byte 1
                               ((num<<24)&0xff000000); // byte 0 to byte 3
        *x = swapped;
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

    int read_segyval(float *val, int x,int y,int o,int n,bool byteswap = false) const 
    {

        return 0;
    }

};

#endif


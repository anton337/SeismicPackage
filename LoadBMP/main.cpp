#include <iostream>
#include <string.h>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <stdint.h>
#include <stdlib.h>
#include <limits>
#include "sep_writer.h"

int width = 0;
int height = 0;
short BitsPerPixel = 0;
std::vector<unsigned char> Pixels;

void LoadBitmap(const char* FilePath)
{
    std::fstream hFile(FilePath, std::ios::in | std::ios::binary);
    if (!hFile.is_open()) throw std::invalid_argument("Error: File Not Found.");

    hFile.seekg(0, std::ios::end);
    int Length = hFile.tellg();
    hFile.seekg(0, std::ios::beg);
    std::vector<uint8_t> FileInfo(Length);
    hFile.read(reinterpret_cast<char*>(FileInfo.data()), 54);

    if(FileInfo[0] != 'B' && FileInfo[1] != 'M')
    {
        hFile.close();
        throw std::invalid_argument("Error: Invalid File Format. Bitmap Required.");
    }

    if (FileInfo[28] != 24 && FileInfo[28] != 32)
    {
        hFile.close();
        throw std::invalid_argument("Error: Invalid File Format. 24 or 32 bit Image Required.");
    }

    BitsPerPixel = FileInfo[28];
    width = FileInfo[18] + (FileInfo[19] << 8);
    height = FileInfo[22] + (FileInfo[23] << 8);
    uint32_t PixelsOffset = FileInfo[10] + (FileInfo[11] << 8);
    uint32_t size = ((width * BitsPerPixel + 31) / 32) * 4 * height;
    Pixels.resize(size);

    hFile.seekg (PixelsOffset, std::ios::beg);
    hFile.read(reinterpret_cast<char*>(Pixels.data()), size);
    hFile.close();
}

int main(int argc , char ** argv)
{
    std::string bmp_file ( argv[1] );
    std::string out_file ( argv[2] );
    float scalar ( atof ( argv[3] ) );
    LoadBitmap(argv[1]);
    float * data = new float [ width * height ];
    float min_elem ( std::numeric_limits<float>::max() );
    float max_elem ( std::numeric_limits<float>::min() );
    float elem;
    for ( int h(0)
        , k(0)
        ; h < height
        ; ++h
        )
    {
        for ( int w(0)
            ; w < width
            ; ++w
            , ++k
            )
        {
            elem = (float)(Pixels[3*k])*scalar;
            if ( elem < min_elem )
            {
                min_elem = elem;
            }
            if ( elem > max_elem )
            {
                max_elem = elem;
            }
            data[height-1-h+height*w] = elem;
        }
    }
    std::cout << "min:" << min_elem << std::endl;
    std::cout << "max:" << max_elem << std::endl;
    std::vector < std::string > header_labels;
    std::vector < std::string > sort_order;
    SEPWriter writer ( out_file . c_str () 
                     , 0 , 1 , height
                     , 0 , 1 , width
                     , 0 , 1 , 1
                     , header_labels
                     , sort_order
                     , (out_file + std::string("@")) . c_str()
                     );
    writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
    writer . write_sepval ( (float*)data , writer . o1 , writer . o2 , writer . o3 , width * height );
}


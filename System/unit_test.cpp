#include <iostream>
#include <math.h>
#include "distributed_data_set.h"

int main()
{
    std::cout << "hello from unit_test" << std::endl;
    long size[3] = { 1024 , 1024 , 1024 };
    long tile[3] = { 256  , 256  , 256  };
    long pad [3] = { 8    , 8    , 8    };
    DistributedDataSet < float , 3 > data_set ( size 
                                              , tile 
                                              , pad 
                                              , "file" 
                                              , false
                                              );
    clock_t start = clock();
    long tile_size ( (data_set.tile[0]+data_set.pad[0])
                   * (data_set.tile[1]+data_set.pad[1])
                   * (data_set.tile[2]+data_set.pad[2])
                   );
    float * val = new float [ tile_size ];
    for ( int k(0)
        ; k < tile_size
        ; ++k
        )
    {
        val[k] = 2000;
    }
    for ( int k(0)
        , x(0)
        ; x < data_set.ind[0]
        ; ++x
        )
    {
        for ( int y(0)
            ; y < data_set.ind[1]
            ; ++y
            )
        {
            for ( int z(0)
                ; z < data_set.ind[2]
                ; ++z
                , ++k
                )
            {
                data_set . set_tile ( k , tile_size , val );
            }
        }
    }
    clock_t end = clock();
    std::cout << "runtime=" << end - start << std::endl;
    return 0;
}


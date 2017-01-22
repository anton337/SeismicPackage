#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include <algorithm>
#include "sep_writer.h"
#include "sep_reader.h"
#include "viewer.h"
#include "producer_consumer.h"
#include "networking.h"

long total_ind_x = 3;//12;
long total_ind_y = 3;//12;
long total_ind_z = 3;//12;

template < typename T , int N = 3 >
struct Edge
{
    long pos[N]; // +++ corner ; local
    long neg[N]; // --- corner ; local
    long pos_out[N]; // +++ corner ; local
    long neg_out[N]; // --- corner ; local
    long size[N];
    T * data_1;
    Edge ( long * _pos_in // +++ corner ; local
         , long * _neg_in // --- corner ; local
         , long * _pos_out// +++ corner ; local
         , long * _neg_out// --- corner ; local
         , T * _data_1
         , long * _size
         )
    {
        long total_size = 1;
        for ( int k(0)
            ; k < N
            ; ++k
            )
        {
            total_size *= _pos_in[k] - _neg_in[k];
            size[k] = _size[k];
            pos[k] = _pos_in[k];
            neg[k] = _neg_in[k];
            pos_out[k] = _pos_out[k];
            neg_out[k] = _neg_out[k];
        }
        if ( total_size <= 0 )
        {
            std::cout << "edge boundaries error." << std::endl;
            exit(1);
        }
        data_1 = new T [ total_size ];
        update ( _data_1 );
    }
    void update ( T * _data_1 )
    {
        for ( int x(neg[0])
            , k(0)
            ; x < pos[0]
            ; ++x
            )
        {
            for ( int y(neg[1])
                , ind2
                ; y < pos[1]
                ; ++y
                )
            {
                ind2 = size[2] * ( y + size[1] * x );
                for ( int z(neg[2])
                    ; z < pos[2]
                    ; ++z
                    , ++k
                    )
                {
                    data_1[k] = _data_1[ind2+z];
                }
            }
        }
    }
    void apply ( T * _data_1 )
    {
        for ( int x(neg_out[0])
            , k(0)
            ; x < pos_out[0]
            ; ++x
            )
        {
            for ( int y(neg_out[1])
                , ind2
                ; y < pos_out[1]
                ; ++y
                )
            {
                ind2 = size[2] * ( y + size[1] * x );
                for ( int z(neg_out[2])
                    ; z < pos_out[2]
                    ; ++z
                    , ++k
                    )
                {
                    _data_1[ind2+z] = data_1[k];
                }
            }
        }
    }
};

Message & operator << ( Message & m , Edge<float> const & e )
{
    long ind = 0;
    long chunk = 0;
    long * tmp;
    float * tmp_f;
    Chunk c;
    m . data . push_back ( c );
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos[0]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos[1]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos[2]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg[0]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg[1]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg[2]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos_out[0]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos_out[1]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . pos_out[2]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg_out[0]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg_out[1]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . neg_out[2]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . size[0]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . size[1]; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    *tmp = e . size[2]; ind += 4;
    for ( int k(0)
        , total_size ( e . size[0] * e . size[1] * e . size[2] )
        ; k < total_size
        ; ++k
        )
    {
        tmp_f = reinterpret_cast < float * > ( m . data[chunk] . buffer[ind] );
        *tmp_f = e . data_1[k]; ind += 4; if ( ind % 256 == 0 ) { ind = 0; chunk++; Chunk c; m . data . push_back ( c ); }
    }
    return m;
}

Message & operator >> ( Message & m , Edge<float> & e )
{
    long ind = 0;
    long chunk = 0;
    long * tmp;
    float * tmp_f;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos[0] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos[1] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos[2] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg[0] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg[1] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg[2] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos_out[0] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos_out[1] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . pos_out[2] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg_out[0] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg_out[1] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . neg_out[2] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . size[0] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . size[1] = *tmp; ind += 4;
    tmp = reinterpret_cast < long * > ( m . data[chunk] . buffer[ind] );
    e . size[2] = *tmp; ind += 4;
    for ( int k(0)
        , total_size ( e . size[0] * e . size[1] * e . size[2] )
        ; k < total_size
        ; ++k
        )
    {
        tmp_f = reinterpret_cast < float * > ( m . data[chunk] . buffer[ind] );
        e . data_1[k] = *tmp_f; ind += 4; if ( ind % 256 == 0 ) { ind = 0; chunk++; }
    }
    return m;
}

template < typename T , int N , long size_x , long size_y , long size_z , long len >
struct WavePropagationTile;

bool toggle = false;

//ProducerConsumer < WavePropagationTile<float,3,64,64,64,262144> * > propagation_queue(total_ind_x*total_ind_y*total_ind_z);
std::vector      < WavePropagationTile<float,3,64,64,64,262144> * > propagation_vect ;

template < typename T , int N = 3 , long size_x = 64 , long size_y = 64 , long size_z = 64 , long len = 262144 >
struct WavePropagationTile
{
    T * data;
    T * data_1;
    T * data_2;
    T * vel_2;
    long long boundary;
    std::string file_name;
    long pos[N]; // +++ corner ; global
    long neg[N]; // --- corner ; global
    float scale[N];
    float h;
    long index[N];
    Edge<T,N> ** side;
    bool data_present[6];
    long num_sides;
    ProducerConsumer<int> * propagation_barrier;
    WavePropagationTile<T,N,size_x,size_y,size_z,len> ** neighbors;
    long long timer;
    WavePropagationTile( float _h
                       , float * _scale
                       , long * _index
                       , long * _boundary_index
                       )
    {
        timer = 0;
        h = _h;
        data   = new T[len];
        data_1 = new T[len];
        data_2 = new T[len];
        vel_2  = new T[len];
        for ( int k(0)
            ; k < len
            ; ++k
            )
        {
            vel_2[k] = 2000*2000;
        }
        for ( int k(0)
            ; k < N
            ; ++k
            )
        {
            scale[k] = _scale[k];
            index[k] = _index[k];
        }
        {
            side = new Edge<T,N>*[2*N];

            for ( int k(0)
                ; k < N
                ; ++k
                )
            {
                data_present[2*k  ] = index[k]!=0;
                data_present[2*k+1] = index[k]!=_boundary_index[k]-1;
                long size[N];
                size[0] = size_x;
                size[1] = size_y;
                size[2] = size_z;
                long _neg[N];
                long _pos[N];
                long _neg_out[N];
                long _pos_out[N];
                if ( data_present[2*k] )
                {
                    _neg[0] = 1;
                    _neg[1] = 1;
                    _neg[2] = 1;
                    _pos[0] = size_x-1;
                    _pos[1] = size_y-1;
                    _pos[2] = size_z-1;
                    _pos[k] = 2;
                    _neg_out[0] = 1;
                    _neg_out[1] = 1;
                    _neg_out[2] = 1;
                    _pos_out[0] = size_x-1;
                    _pos_out[1] = size_y-1;
                    _pos_out[2] = size_z-1;
                    _neg_out[k] = size[k]-1;
                    _pos_out[k] = size[k];
                    side[2*k] = new Edge<T,N> ( _pos
                                              , _neg
                                              , _pos_out
                                              , _neg_out
                                              , data_1
                                              , size
                                              );
                }
                if ( data_present[2*k+1] )
                {
                    _neg[0] = 1;
                    _neg[1] = 1;
                    _neg[2] = 1;
                    _pos[0] = size_x-1;
                    _pos[1] = size_y-1;
                    _pos[2] = size_z-1;
                    _neg[k] = size[k]-2;
                    _neg_out[0] = 1;
                    _neg_out[1] = 1;
                    _neg_out[2] = 1;
                    _pos_out[0] = size_x-1;
                    _pos_out[1] = size_y-1;
                    _pos_out[2] = size_z-1;
                    _neg_out[k] = 0;
                    _pos_out[k] = 1;
                    side[2*k+1] = new Edge<T,N> ( _pos
                                                , _neg
                                                , _pos_out
                                                , _neg_out
                                                , data_1
                                                , size
                                                );
                }
            }

        }
        num_sides = 0;
        for ( int k(0)
            ; k < 6
            ; ++k
            )
        {
            if ( data_present[k] )
            {
                num_sides ++;
            }
        }
        // propagation_barrier = new ProducerConsumer<int>(num_sides);
    }
    void propagate()
    {
        timer++;
        float fact ( h * h * h / ( scale[0] * scale[1] * scale[2] ) );
        long size_2 = size_y*size_z;
        for ( long x(1)
            ; x+1 < size_x
            ; ++x
            )
        {
            for ( long y(1) , ind_pref
                ; y+1 < size_y
                ; ++y
                )
            {
                ind_pref = size_z*(y+size_y*x);
                for ( long z(1) , ind
                    ; z+1 < size_z
                    ; ++z
                    )
                {
                    ind = ind_pref + z;
                    data[ind] = fact * vel_2[ind] * ( -6*data_1[ind] 
                                                    +    data_1[ind+1]
                                                    +    data_1[ind-1]
                                                    +    data_1[ind+size_z]
                                                    +    data_1[ind-size_z]
                                                    +    data_1[ind+size_2]
                                                    +    data_1[ind-size_2]
                                                    )
                                 + 2 * data_1[ind]
                                 -     data_2[ind]
                                 ;
                }
            }
        }
        for ( long k(0)
            ; k < len
            ; ++k
            )
        {
            data_2[k] = data_1[k];
        }
        for ( long k(0)
            ; k < len
            ; ++k
            )
        {
            data_1[k] = data[k];
        }
        update_edges();
    }
    void load()
    {

    }
    void send_edges()
    {
        for ( int k(0)
            ; k < 6
            ; ++k
            )
        {
            if ( data_present[k] )
            {
                //neighbors[k] -> propagation_barrier -> produce ( k );
                side[k] -> apply ( neighbors[k] -> data_1 );
                //if ( neighbors[k] -> propagation_barrier -> size () == neighbors[k] -> num_sides )
                //{
                //    for ( int j(0)
                //        ; j < 6
                //        ; ++j
                //        )
                //    {
                //        if ( neighbors[k] -> data_present[j] )
                //        {
                //            neighbors[k] -> propagation_barrier -> consume ();
                //        }
                //    }
                //    //propagation_queue . produce ( neighbors[k] );
                //}
            }    
        }

    }
    void update_edges()
    {
        for ( long k(0)
            ; k < 6
            ; ++k
            )
        {
            if ( data_present[k] )
            {
                side[k] -> update ( data_1 );
                side[k] -> update ( data_1 );
            }
        }

    }
};


void propagate_worker ( int k , int n )
{
    for ( int i(k) ; i < propagation_vect . size () ; i+=n )
    {
        propagation_vect [i] -> propagate ( ) ;
    }
}

void send_edges_worker ( int k , int n )
{
    for ( int i(k) ; i < propagation_vect . size () ; i+=n )
    {
        propagation_vect [i] -> send_edges ( ) ;
    }
}

float get_argument ( int &   arg 
                   , int     argc
                   , char ** argv
                   )
{
    ++arg;
    if ( arg > argc )
    {
        //std::cout << "parameter error" << std::endl;
        exit(1);
    }
    return atof ( argv[arg] );
}

template < typename T >
inline void mem_cpy ( T * out , T * in , int num )
{
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        out[k] = in[k];
    }
}

void propagate ( float   h       // dt
               , int     num_x 
               , int     num_t
               , float   scale_x
               , float   scale_t
               , float * p2 
               , float * p1 
               , float * c 
               , float * vel )
{
    float fact ( h * h / ( scale_x * scale_x ) );
    for ( int t(1)
        ; t+1 < num_t
        ; ++t
        )
    {
        for ( int x(1)
            ; x+1 < num_x
            ; ++x
            )
        {
            if ( fabs ( p2[t+x*num_t] ) > 1e-40 
              || fabs ( p1[t+1+x*num_t] ) > 1e-40 
              || fabs ( p1[t-1+x*num_t] ) > 1e-40 
              || fabs ( p1[t+(x+1)*num_t] ) > 1e-40 
              || fabs ( p1[t+(x-1)*num_t] ) > 1e-40 
               )
            {
                c[t+x*num_t] = fact * vel[t+x*num_t] * vel[t+x*num_t] * ( -4*p1[t+x*num_t] 
                                                                        +    p1[t+(x+1)*num_t]
                                                                        +    p1[t+(x-1)*num_t]
                                                                        +    p1[(t+1)+x*num_t]
                                                                        +    p1[(t-1)+x*num_t]
                                                                        )
                             + 2 * p1[t+x*num_t]
                             -     p2[t+x*num_t]
                             ;
            }
        }
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p2[k] = p1[k];
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p1[k] = c[k];
    }
}

float find_max ( float const * x , int num )
{
    float m ( 0 ) , t;
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        t = fabs(x[k]);
        if ( t > m )
        {
            m = t;
        }
    }
    return m;
}

float find_min ( float const * x , int num )
{
    float m ( 1e10 ) , t;
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        t = fabs(x[k]);
        if ( t < m )
        {
            m = t;
        }
    }
    return m;
}

struct source
{

    float pos_x; // m
    float pos_z; // m

    source ( int _pos_x , int _pos_z )
    : pos_x ( _pos_x )
    , pos_z ( _pos_z )
    {

    }

};

struct recorder
{

    float pos_x; // m
    float pos_z; // m

    recorder ( int _pos_x , int _pos_z )
    : pos_x ( _pos_x )
    , pos_z ( _pos_z )
    {

    }

};

void calculate_poyinting_vector ( vect        * p 
                                , float const * u
                                , int           num_x
                                , int           num_z
                                )
{
    int disp = 1;
    for ( int x(disp)
        ; x+disp < num_x
        ; ++x
        )
    {
        for ( int z(disp)
            ; z+disp < num_z
            ; ++z
            )
        {
            p[z+x*num_z].vx += (u[z+(x+disp)*num_z] - u[z+(x-disp)*num_z]);
            p[z+x*num_z].vz += (u[z+disp  +x*num_z] - u[z-disp+  x*num_z]);
            p[z+x*num_z].vy = 0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  wave propagation based synthetic trace data generator :
//
//  required fields:
//
//  input: 
//      - velocity SEP header file
//      - velocity model dimensions     // get from SEP file
//
//      - grid width x
//      - grid width z
//
//      - wavelet file
//
//      - sampling rate
//      - num_t                         // number of time samples to generate
//      - propagation step size
//
//      - num_rec
//      - num_sou
//      - rec_disp
//      - sou_disp
//      - rec_init
//      - sou_init
//
//     
//  output:   
//      - seismic SEP header file
//      - seismic SEP data files
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main( int argc , char ** argv )
{

    float scale[3];
    scale[0] = 25;
    scale[1] = 25;
    scale[2] = 25;
    WavePropagationTile<float,3,64,64,64,262144> **** neigh;
    neigh = new
    WavePropagationTile<float,3,64,64,64,262144> *** [total_ind_x];
    for ( long x(0)
        ; x < total_ind_x
        ; ++x
        )
    {
    neigh[x] = new
    WavePropagationTile<float,3,64,64,64,262144> ** [total_ind_y];
    for ( long y(0)
        ; y < total_ind_y
        ; ++y
        )
    {
    neigh[x][y] = new
    WavePropagationTile<float,3,64,64,64,262144> * [total_ind_z];
    for ( long z(0)
        ; z < total_ind_z
        ; ++z
        )
    {
        long index[3];
        index[0] = x;
        index[1] = y;
        index[2] = z;
        long total_index[3];
        total_index[0] = total_ind_x;
        total_index[1] = total_ind_y;
        total_index[2] = total_ind_z;
        WavePropagationTile<float,3,64,64,64,262144> * v =
        new WavePropagationTile<float,3,64,64,64,262144> ( 1e-5
                                                         , scale
                                                         , index
                                                         , total_index
                                                         );
        neigh[x][y][z] = v;
        propagation_vect  . push_back ( v );
        // propagation_queue . produce ( v
        //                             );
    }
    }
    }
    
    for ( long x(0)
        ; x < total_ind_x
        ; ++x
        )
    {
    for ( long y(0)
        ; y < total_ind_y
        ; ++y
        )
    {
    for ( long z(0)
        ; z < total_ind_z
        ; ++z
        )
    {
        WavePropagationTile<float,3,64,64,64,262144> ** n = new WavePropagationTile<float,3,64,64,64,262144> * [6];
        if ( x == 0 )
        {
            n[0] = NULL;
        }
        else
        {
            n[0] = neigh[x-1][y][z];
        }
        if ( x == total_ind_x-1 )
        {
            n[1] = NULL;
        }
        else
        {
            n[1] = neigh[x+1][y][z];
        }
        if ( y == 0 )
        {
            n[2] = NULL;
        }
        else
        {
            n[2] = neigh[x][y-1][z];
        }
        if ( y == total_ind_y-1 )
        {
            n[3] = NULL;
        }
        else
        {
            n[3] = neigh[x][y+1][z];
        }
        if ( z == 0 )
        {
            n[4] = NULL;
        }
        else
        {
            n[4] = neigh[x][y][z-1];
        }
        if ( z == total_ind_z-1 )
        {
            n[5] = NULL;
        }
        else
        {
            n[5] = neigh[x][y][z+1];
        }
        neigh[x][y][z] -> neighbors = n;
    }
    }
    }
    
    clock_t start = clock();
    int n_thread = 4;
    for ( int it_k(0) ; it_k < 10 ; ++it_k  )
    {
        std::cout << it_k << std::endl;
        std::vector < boost::thread * > thread_th;
        for ( int k(0)
            ; k < n_thread
            ; ++k
            )
        {
            boost::thread * thr1 = new boost::thread ( propagate_worker , k , n_thread );
            thread_th . push_back ( thr1 );
        }
        for ( int k(0)
            ; k < thread_th.size()
            ; ++k
            )
        {
            thread_th[k] -> join ();
            delete thread_th[k];
        }
        std::vector < boost::thread * > thread_send_th;
        for ( int k(0)
            ; k < n_thread
            ; ++k
            )
        {
            boost::thread * thr1 = new boost::thread ( send_edges_worker , k , n_thread );
            thread_send_th . push_back ( thr1 );
        }
        for ( int k(0)
            ; k < thread_send_th.size()
            ; ++k
            )
        {
            thread_send_th[k] -> join ();
            delete thread_send_th[k];
        }
        // for ( int k(0)
        //     ; k < propagation_vect . size ()
        //     ; ++k
        //     )
        // {
        //     propagation_queue . produce ( propagation_vect [ k ] );
        // }
    }
    clock_t end = clock();
    //std::cout << "elapsed time = " << (end - start)/4 << std::endl;
    clock_t start_1 = clock();
    //for ( int k(0) ; k < total_num ; ++k )
    //tile_1 . propagate ();
    clock_t end_1 = clock();
    //std::cout << "elapsed time = " << end_1 - start_1 << std::endl;

/*
    boost::thread       * thr = show ( argc , argv );


    if (argc < 3)
    {
        //std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

//      - velocity SEP header file
    arg++;
    std::string velocity_file_name( argv[arg] );

//      - velocity model dimensions     // get from SEP file
    SEPReader vel_reader ( velocity_file_name . c_str () , false );
    int vel_num_z = vel_reader . n1;
    int vel_num_x = vel_reader . n2;
    float * vel_data = new float [ vel_num_x * vel_num_z ];
    memset ( &vel_data[0] , 0 , vel_num_x * vel_num_z );
    vel_reader . OpenDataFile ( (velocity_file_name + std::string("@")) . c_str () );
    vel_reader . read_sepval ( & vel_data [ 0 ] , vel_reader . o1 , vel_reader . o2 , vel_reader . o3 , vel_num_x * vel_num_z );

//      - grid width x
    float grid_width_x = get_argument ( arg , argc , argv );
//      - grid width z
    float grid_width_z = get_argument ( arg , argc , argv );

//      - wavelet file
    arg++;
    std::string wavelet_file_name( argv[arg] );
    SEPReader wav_reader ( wavelet_file_name . c_str () , false );
    int wav_num_t = wav_reader . n1;
    float * wav_data = new float [ wav_num_t ];
    memset ( &wav_data[0] , 0 , wav_num_t );
    wav_reader . OpenDataFile ( (wavelet_file_name + std::string("@")) . c_str () );
    wav_reader . read_sepval ( & wav_data [ 0 ] , wav_reader . o1 , wav_reader . o2 , wav_reader . o3 , wav_num_t );

//      - sampling rate
    float sample_rate = 0.001 * get_argument ( arg , argc , argv );
//      - num_t                         // number of time samples to generate
    float max_t = 0.001 * get_argument ( arg , argc , argv );
    int num_t = max_t / sample_rate;
//      - propagation step size
    float h = 0.001 * get_argument ( arg , argc , argv );

//      - num_rec
    int num_rec = get_argument ( arg , argc , argv );
//      - num_sou
    int num_sou = get_argument ( arg , argc , argv );
//      - rec_disp
    float rec_disp = get_argument ( arg , argc , argv );
//      - sou_disp
    float sou_disp = get_argument ( arg , argc , argv );
//      - rec_init
    float rec_init = get_argument ( arg , argc , argv );
//      - sou_init
    float sou_init = get_argument ( arg , argc , argv );

//      - seismic output SEP header file
    arg++;
    std::string seismic_file_header_name( argv[arg] );

//      - seismic output SEP data file
    arg++;
    std::string seismic_file_data_name( argv[arg] );



// pad velocity volume
    int vel_num_x2 = vel_num_x*3;
    int vel_num_z2 = vel_num_z*3;
    float * vel = new float [ vel_num_x2 * vel_num_z2 ];

    bool water_reflection = false;

    for ( int x(0)
        , _x(0)
        , k(0)
        ; x < vel_num_x2
        ; ++x
        )
    {
        for ( int z(0)
            , _z(0)
            ; z < vel_num_z2
            ; ++z
            , ++k
            )
        {
            _z = z-vel_num_z;
            if ( water_reflection && _z < 0 )
            {
                vel[k] = 0;
            }
            else
            {
                if ( _z < 0 ) _z = 0;
                if ( _z >= vel_num_z ) _z = vel_num_z - 1;
                _x = x-vel_num_x;
                if ( _x < 0 ) _x = 0;
                if ( _x >= vel_num_x ) _x = vel_num_x - 1;
                vel[k] = vel_data[_z+_x*vel_num_z];
            }
        }
    }

    float * data_p2 = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_p2[0] , 0 , vel_num_x2 * vel_num_z2 );

    float * data_p1 = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_p1[0] , 0 , vel_num_x2 * vel_num_z2 );

    float * data_c  = new float [ vel_num_x2 * vel_num_z2 ];
    memset ( &data_c [0] , 0 , vel_num_x2 * vel_num_z2 );

    float * recorder_data  = new float [ num_rec * num_sou * num_t ];
    memset ( &recorder_data [0] , 0 , num_rec * num_sou * num_t );

    vect * point_f = new vect [ vel_num_x2 * vel_num_z2 ];

    vect * point_o = new vect [ vel_num_x2 * vel_num_z2 ];

    set_num_x ( vel_num_x2 );
    set_num_z ( vel_num_z2 );
    set_data_ptr ( vel );
    set_reverse_ptr ( vel );
    set_scalar ( 0.4 / find_max ( vel , vel_num_x2 * vel_num_z2 ) );
    set_reverse_scalar ( 0.4 / find_max ( vel , vel_num_x2 * vel_num_z2 ) );
    set_wave_ptr ( data_c );
    set_data_vec_ptr ( point_o );


    for ( int sou(0)
        ; sou < num_sou
        ; ++sou
        )
    {

        //std::cout << "sou=" << sou << " pos=" << (int)((sou_init + sou_disp * sou) / grid_width_x) << std::endl;
        //std::cout << sou_init << " " << sou_disp << std::endl;

        //source S ( (int)((sou_init + sou_disp * sou) / grid_width_x)
        //         , (int)(0)
        //         );
        
        source S ( (int)(vel_num_x/2)
                 , (int)(vel_num_z/2)
                 );

        std::vector < recorder > recorders;

        for ( int rec(0)
            ; rec < num_rec
            ; ++rec
            )
        {
            recorder R ( (int)((sou_init + sou_disp * sou + rec_init + rec_disp * rec) / grid_width_x)
                       , (int)(0)
                       );
            recorders . push_back ( R );
        }

        float sample_time ( sample_rate / h );

        for ( int k(0) 
            ; 1 
            ; ++k 
            )
        {
            if ( k < num_t * sample_time )
            {
                // feed wavelet 
                if ( (int)(k*sample_time) < wav_num_t )
                {
                    if ( (int)(k*sample_time) >= 1 )
                    {
                        data_p2[(int)(vel_num_z+S.pos_z)+(int)(vel_num_x+S.pos_x)*vel_num_z2] = wav_data [ (int)(k*sample_time)-1 ];
                    }
                    data_p1[(int)(vel_num_z+S.pos_z)+(int)(vel_num_x+S.pos_x)*vel_num_z2] = wav_data [ (int)(k*sample_time) ];
                }
            }
            else
            {
                break;
            }
            propagate ( h 
                      , vel_num_x2 
                      , vel_num_z2 
                      , grid_width_x 
                      , grid_width_z 
                      , data_p2 
                      , data_p1 
                      , data_c
                      , vel  
                      );
            set_wave_scalar ( 0.4 / find_max ( data_c , vel_num_x2 * vel_num_z2 ) );
            set_timer ( 1000 * k * h );
            calculate_poyinting_vector ( point_f
                                       , data_p1
                                       , vel_num_x2
                                       , vel_num_z2
                                       );
            for ( int x(1)
                ; x+1 < vel_num_x2
                ; ++x
                )
            for ( int z(1)
                ; z+1 < vel_num_z2
                ; ++z
                )
            {
                {
                    point_o[z+vel_num_z2*x].vx = -k*(point_f[z+vel_num_z2*(x)].vx*data_p1[z+vel_num_z2*(x)]);
                    point_o[z+vel_num_z2*x].vz = -k*(point_f[z+vel_num_z2*(x)].vz*data_p1[z+vel_num_z2*(x)]);
                }
            }
            for ( std::size_t rec(0)
                ; rec < recorders . size ()
                ; ++rec
                )
            {
                int t ( k * sample_time );
                // int ind ( (int)(vel_num_z+recorders[rec].pos_z) + (int)(vel_num_x+recorders[rec].pos_x)*vel_num_z2 );
                recorder_data[ sou*num_rec*num_t + rec*num_t + t ] = ( data_p1[(int)(vel_num_z+recorders[rec].pos_z)+(int)(vel_num_x+recorders[rec].pos_x)*vel_num_z2] );// * fabs ( point_f[ind].vx ) / ( sqrtf ( point_f[ind].vx*point_f[ind].vx + point_f[ind].vz*point_f[ind].vz ) + 1e-5 );
            }
        }

        for ( int k(0)
            , size = vel_num_x2 * vel_num_z2
            ; k < size
            ; ++k
            )
        {
            data_p2[k] = 0;
            data_p1[k] = 0;
            data_c [k] = 0;
        }

    }

    SEPWriter writer ( seismic_file_header_name . c_str () 
                     , 0 , 1000*sample_rate , num_t
                     , 1 , 1 , num_rec
                     , 1 , 1 , num_sou
                     , wav_reader . get_header_labels ()
                     , wav_reader . get_sort_order ()
                     , seismic_file_data_name . c_str()
                     );

    writer . OpenDataFile ( (seismic_file_data_name) . c_str() );

    writer . write_sepval ( (float*)&recorder_data[0] , writer . o1 , writer . o2 , writer . o3 , num_sou*num_rec*num_t );



    thr -> join ();

*/

    std::cout << "exit" << std::endl;

    return 0;

}


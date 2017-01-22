#ifndef MATRIX_H
#define MATRIX_H

template < typename T >
class matrix
{
public:
    int m;
    int n;
    T * array;
    matrix ()
    {

    }
    matrix ( int _m
           , int _n
           , T * _array )
    : m ( _m )
    , n ( _n )
    , array ( _array )
    {

    }
    void operator () ( T const * x 
                     , T       * b 
                     ) const
    {
        for ( int _m(0)
            , k(0)
            ; _m < m
            ; ++_m
            )
        {
            b[_m] = 0;
            for ( int _n(0)
                ; _n < n
                ; ++_n
                , ++k
                )
            {
                b[_m] += array[k]*x[_n];
            }
        }
    }
};

template <  >
class matrix < fftwf_complex >
{
public:
    typedef fftwf_complex T;
    int m;
    int n;
    T * array;
    matrix ()
    {

    }
    matrix ( int _m
           , int _n
           , T * _array )
    : m ( _m )
    , n ( _n )
    , array ( _array )
    {

    }
    void operator () ( T const * x 
                     , T       * b 
                     ) const
    {
        for ( int _m(0)
            , k(0)
            ; _m < m
            ; ++_m
            )
        {
            b[_m][0] = 0;
            b[_m][1] = 0;
            for ( int _n(0)
                ; _n < n
                ; ++_n
                , ++k
                )
            {
                b[_m][0] += array[k][0]*x[_n][0] - array[k][1]*x[_n][1];
                b[_m][1] += array[k][0]*x[_n][1] + array[k][1]*x[_n][0];
            }
        }
    }
    void adjoint ( T const * x 
                 , T       * b 
                 ) const
    {
        for ( int _n(0)
            ; _n < n
            ; ++_n
            )
        {
            b[_n][0] = 0;
            b[_n][1] = 0;
        }
        for ( int _m(0)
            , k(0)
            ; _m < m
            ; ++_m
            )
        {
            for ( int _n(0)
                ; _n < n
                ; ++_n
                , ++k
                )
            {
                b[_n][0] += array[k][0]*x[_m][0] + array[k][1]*x[_m][1];
                b[_n][1] += array[k][0]*x[_m][1] - array[k][1]*x[_m][0];
            }
        }
    }
};

template < typename T >
class toeplitz_matrix
{
public:
    int m;
    int n;
    T * array;
    toeplitz_matrix ()
    {

    }
    toeplitz_matrix ( int _m
                    , int _n
                    , T * _array 
                    )
    : m ( _m )
    , n ( _n )
    , array ( _array )
    {

    }
    void operator () ( T const * x 
                     , T       * b 
                     ) const
    {
        for ( int _m(0)
            , k(0)
            ; _m < m
            ; ++_m
            )
        {
            b[_m] = 0;
            for ( int _n(0)
                ; _n < n
                ; ++_n
                , ++k
                )
            {
                if ( _m >= _n )
                {
                    b[_m] += array[(_n+_m)%m]*x[_n];
                }
                else
                {
                    b[_m] += array[(_n+_m)%m]*x[_n];
                }
            }
        }
    }
};

template <  >
class toeplitz_matrix < fftwf_complex >
{
public:
    typedef fftwf_complex T;
    int m;
    int n;
    T * array;
    toeplitz_matrix ()
    {

    }
    toeplitz_matrix ( int _m
                    , int _n
                    , T * _array 
                    )
    : m ( _m )
    , n ( _n )
    , array ( _array )
    {

    }
    void operator () ( T const * x 
                     , T       * b 
                     ) const
    {
        for ( int _m(0)
            , k(0)
            ; _m < m
            ; ++_m
            )
        {
            b[_m][0] = 0;
            b[_m][1] = 0;
            for ( int _n(0)
                ; _n < n
                ; ++_n
                , ++k
                )
            {
                if ( _m >= _n )
                {
                    b[_m][0] += array[(-_n+m+_m)%m][0]*x[_n][0] - array[(-_n+m+_m)%m][1]*x[_n][1];
                    b[_m][1] += array[(-_n+m+_m)%m][0]*x[_n][1] + array[(-_n+m+_m)%m][1]*x[_n][0];
                }
                else
                {
                    b[_m][0] += array[(_n+m-_m)%m][0]*x[_n][0] + array[(_n+m-_m)%m][1]*x[_n][1];
                    b[_m][1] += array[(_n+m-_m)%m][0]*x[_n][1] - array[(_n+m-_m)%m][1]*x[_n][0];
                }
            }
        }
    }
};
#endif


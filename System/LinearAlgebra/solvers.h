#ifndef SOLVERS_H
#define SOLVERS_H

#include "matrix.h"

template < typename T , typename functor >
class BisectionSearch
{

    typedef T value_type;

public:
    value_type operator () ( value_type const & x
                           , functor    const & apply
                           ) const
    {
        value_type disp ( 1e-4 );
        int dir = 1;
        value_type val1 ( apply ( x ) );
        value_type val2 ( apply ( x+disp ) );
        value_type a , b , c ;
        a = x;
        if ( val2 > val1 )
        {
            dir = -1;
        }
        for ( int iter = 0
            ; iter < 100
            ; ++iter
            )
        {
            b = x + dir * disp;
            value_type val2 ( apply ( b ) );
            if ( val2 < val1 )
            {
                disp *= 2;
            }
            else
            {
                break;
            }
        }
        for ( int iter = 0
            ; iter < 10
            ; ++iter
            )
        {
            c = 0.5f * ( a + b );
            value_type val_mid ( apply ( c ) );
            if ( val1 < val2 )
            {
                b = c;
                val2 = val_mid;
            }
            else
            {
                a = c;
                val1 = val_mid;
            }
        }
        return a;
    }

};

template < typename T , typename functor >
class ZeroFinder
{

    typedef T value_type;

public:
    value_type operator () ( value_type const & _s
                           , value_type const & _t
                           , functor    const & apply
                           ) const
    {
        value_type r,fr;
        int n, side=0;
        /* starting values at endpoints of interval */
        value_type  s = _s;
        value_type  t = _t;
        value_type fs = apply(s);
        value_type ft = apply(t);
        int m = 1000000;

        for (n = 0; n < m; n++)
        {
            r = (fs*t - ft*s) / (fs - ft);
            if (fabs(t-s) < 1e-8*fabs(t+s)) 
            {
                return r;
            }
            fr = apply(r);
            if (fr * ft > 0)
            {
                /* fr and ft have same sign, copy r to t */
                t = r; ft = fr;
                if (side==-1) fs /= 2;
                side = -1;
            }
            else if (fs * fr > 0)
            {
                /* fr and fs have same sign, copy r to s */
                s = r;  fs = fr;
                if (side==+1) ft /= 2;
                side = +1;
            }
            else if (fs * fr < 0 && fr * ft < 0)
            {
                // t = r; ft = fr;
                // if (side==-1) fs /= 2;
                // side = -1;

                s = r;  fs = fr;
                if (side==+1) ft /= 2;
                side = +1;
            }
            else
            {
                /* fr * f_ very small (looks like zero) */
                break;
            } 
        }
        if ( r > s+0.0001*fabs(fs) && r < t-0.0001*fabs(ft) )
        {
            return r;
        }
        else
        {
            std::cout << "zero not found" << std::endl;
            return -1;
        }
    }

};

template < typename T >
class LinearSystemEnergyL2
{
    T * b_candidate;
    T const * b;
    matrix < T > const A;
public:
    LinearSystemEnergyL2 ( matrix < T > const & _A
                         , T            const * _b
                         )
    : A ( _A )
    , b ( _b )
    {
        b_candidate = new T [ A . m ];
    }
    ~LinearSystemEnergyL2()
    {
        delete [] b_candidate;
    }
    T operator () ( T const * x ) const
    {
        T ret = 0;
        A ( x , b_candidate );
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            ret += ( b_candidate[k] - b[k] ) * ( b_candidate[k] - b[k] );
        }
        return ret;
    }
};

template < typename T >
class LinearSystemEnergyL1
{
    T * b_candidate;
    T const * b;
    matrix < T > const A;
public:
    LinearSystemEnergyL1 ( matrix < T > const & _A
                         , T            const * _b
                         )
    : A ( _A )
    , b ( _b )
    {
        b_candidate = new T [ A . m ];
    }
    ~LinearSystemEnergyL1()
    {
        delete [] b_candidate;
    }
    T operator () ( T const * x ) const
    {
        T ret = 0;
        A ( x , b_candidate );
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            ret += fabs( b_candidate[k] - b[k] );
        }
        return ret;
    }
};

template < typename T , typename functor >
class QuasiNewtonsMethodSolver
{

    void print ( std::string   str
               , T const     * a
               , int           num
               ) const
    {
        std::cout << str << "\t";
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            std::cout << a[k] << "\t";
        }
        std::cout << std::endl;
    }

    T dot ( T const * a
          , T const * b
          , int       n
          ) const
    {
        T ret = 0;
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            ret += a[k] * b[k];
        }
        return ret;
    }

    void gram ( T       * a
              , T const * b
              , T const * c
              , int       n
              , float     scal = 1
              ) const
    {
        for ( int i(0)
            , k(0)
            ; i < n
            ; ++i
            )
        {
            for ( int j(0)
                ; j < n
                ; ++j
                , ++k
                )
            {
                a[k] += scal * b[i] * c[j];
            }
        }
        // for ( int i(0)
        //     , k(0)
        //     ; i < n
        //     ; ++i
        //     )
        // {
        //     for ( int j(0)
        //         ; j < n
        //         ; ++j
        //         , ++k
        //         )
        //     {
        //         a[k] = a[i+j*n];
        //     }
        // }
    }

    void add ( T       * a
             , T const * b
             , T const * c
             , int       n
             , float     scal = 1
             ) const
    {
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            a[k] = b[k] + scal * c[k];
        }
    }

    void calculate_gradient ( functor const & f
                            , T       const * x
                            , T             * G
                            , int             num
                            , float           h
                            ) const
    {
        T * x_tmp = new T [ num ];
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            x_tmp[k] = x[k];
        }
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            x_tmp[k] += h;
            G[k] = f ( x_tmp );
            x_tmp[k] -= 2*h;
            G[k] -= f ( x_tmp );
            x_tmp[k] += h;
            G[k] /= 2*h;
        }
        delete [] x_tmp;
    }

public:

    void operator () ( functor const & f 
                     , T             * x_init  
                     , int             num
                     , T               delta = 1e-2
                     ) const
    {
        T * G_p = new T [ num ];
        T * G   = new T [ num ];
        T * X   = new T [ num ];
        T * Y   = new T [ num ];
        T * H   = new T [ num * num ];
        T * B   = new T [ num * num ];
        T * Hy  = new T [ num ];
        T * Bdx = new T [ num ];
        T * y   = new T [ num ];
        T * dx  = new T [ num ];
        T * x   = x_init;
        T scalar;
        for ( int i(0),k(0); i<num; ++i )
        for ( int j(0); j<num; ++j, ++k )
        {
            B[k] = (i==j)?1:0;
            H[k] = (i==j)?1:0;
        }
        matrix < T > H_op ( num , num , H );
        matrix < T > B_op ( num , num , B );
        H_op ( G , dx );
        calculate_gradient ( f , x , G , num , 1e-5 );
        for ( int k(0); k < num; ++k ) G_p[k] = G[k];
        for ( int iter(0)
            ; iter < 1500
            ; ++iter
            )
        {
            add ( y , G , G_p , num , -1 ); // y = G - G_prev
            H_op ( y , Hy ); 
            add ( X , dx , Hy , num , -1 ); // X = ( dx - Hy )
            // H = ( dx - Hy ) ( dx - Hy )^T / ( dx - Hy )^T y
            scalar = dot ( X , y , num );
            if ( fabs(scalar) < 1e-5 ) scalar = 1; else scalar = 1e-1 / scalar;
            gram ( H , X , X , num , scalar );
            H_op ( G , dx );
            for ( int k(0) ; k < num ; ++k ) dx[k] *= -delta;
            add ( x , x , dx , num , 1 );
            for ( int k(0); k < num; ++k ) G_p[k] = G[k];
            calculate_gradient ( f , x , G , num , 1e-5 );
        }
        std::cout << "error:" << f(x) << std::endl;
    }

};

template < typename T , typename functor >
class GradientDescentSolver
{

    void add ( T       * a
             , T const * b
             , T const * c
             , int       n
             , float     scal = 1
             ) const
    {
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            a[k] = b[k] + scal * c[k];
        }
    }

    void print ( std::string   str
               , T const     * a
               , int           num
               ) const
    {
        std::cout << str << "\t";
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            std::cout << a[k] << "\t";
        }
        std::cout << std::endl;
    }

    class FindZero
    {

        functor const * f;

        int num;

        T const * x_start;

        T const * Grad;

        T * x_candidate;

    public:

        FindZero ( functor const * _f
                 , int _num
                 , T const * _x_start
                 , T const * _Grad
                 )
        : f ( _f )
        , num ( _num )
        , x_start ( _x_start )
        , Grad ( _Grad )
        {
            x_candidate = new T [ num ];
        }

        ~FindZero ()
        {
            delete [] x_candidate;
        }

        T operator () ( T const & x ) const
        {

            for ( int k(0)
                ; k < num
                ; ++k
                )
            {
                x_candidate[k] = x_start[k] + x * Grad[k];
            }

            return (*f) ( x_candidate );

        }

    };

    BisectionSearch < T , FindZero > const bisection_search;

    void calculate_gradient ( functor const & f
                            , T       const * x
                            , T             * G
                            , int             num
                            , float           h
                            ) const
    {
        T * x_tmp = new T [ num ];
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            x_tmp[k] = x[k];
        }
        for ( int k(0)
            ; k < num
            ; ++k
            )
        {
            x_tmp[k] += h;
            G[k] = f ( x_tmp );
            x_tmp[k] -= 2*h;
            G[k] -= f ( x_tmp );
            x_tmp[k] += h;
            G[k] /= 2*h;
        }
        delete [] x_tmp;
    }

public:

    GradientDescentSolver ()
    : bisection_search ()
    {

    }

    void operator () ( functor const & f 
                     , T             * x_init  
                     , int             num
                     ) const
    {
        T delta = 1e-3 , bi_delta ;
        T * x = x_init;
        T * G = new T [ num ];
        int num_iter = 300;
        for ( int iter(0)
            ; iter < num_iter
            ; ++iter
            )
        {
            calculate_gradient ( f , x , G , num , 1e-5 );
            FindZero const find_zero ( &f , num , x , G );
            bi_delta = bisection_search ( delta , find_zero );
            add ( x , x , G , num , bi_delta );
        }
        delete [] G;
    }

};

template < typename T >
class GeneralConjugateGradientSolver
{

    T dot ( T const * a
          , T const * b
          , int       n
          ) const
    {
        T ret = 0;
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            ret += a[k] * b[k];
        }
        return ret;
    }

    void add ( T       * a
             , T const * b
             , int       n
             , float     scal = 1
             ) const
    {
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            a[k] += scal * b[k];
        }
    }

    void add_set ( T       * a
                 , T const * b
                 , T const * c
                 , int       n
                 , float     scal = 1
                 ) const
    {
        for ( int k(0)
            ; k < n
            ; ++k
            )
        {
            a[k] = b[k] + scal * c[k];
        }
    }

public:

    // solve for x in A x = b
    void operator () ( matrix < T > const & A // m x n
                     , T                  * x // n
                     , T            const * b // m
                     ) const
    {
        for ( int k(0)
            ; k < A . n
            ; ++k
            )
        {
            x[k] = 0;
        }

        T * r0 = new T [ A . m ];
        T * r  = new T [ A . m ];
        // r=b-A*x;
        A ( x , r0 );
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            r0[k] = b[k] - r0[k];
            r [k] = r0[k];
        }

        T rho , alpha , w ;

        rho = 1;

        alpha = 1;

        w = 1;

        T rho_prev , alpha_prev , w_prev ;

        rho_prev = rho ;

        alpha_prev = alpha ;

        w_prev = w ;

        T beta ;

        T * p = new T [ A . m ];

        T * v = new T [ A . m ];

        T * s = new T [ A . m ];

        T * t = new T [ A . m ];

        for ( int iter(0)
            ; iter < 10 * A . m
            ; ++iter
            )
        {
       
            // rho = (r0 , r )
            rho = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                rho += r0[k] * r[k];
            }

            if ( iter > 0 )
            {
                // beta = ( rho / rho_prev ) * ( alpha / w_prev )
                beta = ( rho / rho_prev ) * ( alpha / w_prev ) ;
                // p = r_prev - beta * ( p_prev  - w_prev * v_prev );
                for ( int k(0)
                    ; k < A . m
                    ; ++k
                    )
                {
                    p[k] = r[k] + beta * ( p[k] - w_prev * v[k] );
                }
            }
            else
            {
                for ( int k(0)
                    ; k < A . m
                    ; ++k
                    )
                {
                    p[k] = r[k];
                }
            }

            // v = A p
            A ( p , v );

            // alpha=rho/(r0 , v);
            alpha = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                alpha += r0[k] * v[k];
            }
            alpha = rho / alpha;

            // s=r_prev-alpha*v;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                s[k] = r[k] - alpha * v[k];
            }

            // if |s| is sufficiently small, set x = x_prev + alpha * p, and quit
            if ( dot ( s , s , A . m ) < 1e-20 )
            {
                add ( x , p , alpha , A . m ); 
                break;
            }

            // t = A s
            A ( s , t );

            // w = (t,s)/(t,t)
            w = dot ( t , s , A . m ) / dot ( t , t , A . m );

            // x = x_prev + alpha * p + w s
            add ( x , p , A . m , alpha );
            add ( x , s , A . m , w     );



            // r = s - w t
            add_set ( r , s , t , A.m , -w );


            w_prev = w;

            rho_prev = rho;

            alpha_prev = alpha;

        }

        delete []  r  ;
        delete []  r0 ;
        delete []  p  ;
        delete []  v  ;
        delete []  s  ;
        delete []  t  ;

    }

};

template < typename M , typename T >
class SymmetricConjugateGradientSolver
{

public:

    // solve for x in A x = b
    void operator () ( M            const & A // m x n
                     , T                  * x // n
                     , T            const * b // m
                     ) const
    {
        for ( int k(0)
            ; k < A . n
            ; ++k
            )
        {
            x[k] = 0;
        }

        T * r = new T [ A . m ];
        // r=b-A*x;
        A ( x , r );
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            r[k] = b[k] - r[k];
        }
        
        T * p = new T [ A . m ];
        // p=r;
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            p[k] = r[k];
        }

        T rsold = 0;
        // rsold=r'*r;
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            rsold += r[k] * r[k];
        }

        T rsnew , fact , alpha ;

        T * Ap = new T [ A . m ];
        for ( int iter(0)
            ; iter < A . m
            ; ++iter
            )
        {
            A ( p , Ap );
            // alpha=rsold/(p'*Ap);
            alpha = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                alpha += p[k] * Ap[k];
            }
            alpha = rsold / alpha;

            // x=x+alpha*p;
            // r=r-alpha*Ap;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                x[k] += alpha *  p[k];
                r[k] -= alpha * Ap[k];
            }

            // rsnew=r'*r;
            rsnew = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                rsnew += r[k] * r[k];
            }

            if ( rsnew < 1e-10 ) 
            {
                if ( rsnew < 0 )
                {
                    // operator not positive definitive
                    break;
                }
                break;
            }

            //p=r+(rsnew/rsold)*p;
            fact = (rsnew / rsold);
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                p[k] = r[k] + fact * p[k];
            }

            //rsold=rsnew;
            rsold = rsnew;

        }

        delete [] Ap;
        delete []  r;
        delete []  p;

    }

};

template < typename M >
class SymmetricConjugateGradientSolver < M , fftwf_complex >
{

    typedef fftwf_complex T;

public:

    // solve for x in A x = b
    void operator () ( M            const & A // m x n
                     , T                  * x // n
                     , T            const * b // m
                     , int                  num_iter = 0
                     ) const
    {
        for ( int k(0)
            ; k < A . n
            ; ++k
            )
        {
            x[k][0] = 0;
            x[k][1] = 0;
        }

        T * r = new T [ A . m ];
        // r=b-A*x;
        A ( x , r );
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            r[k][0] = b[k][0] - r[k][0];
            r[k][1] = b[k][1] - r[k][1];
        }
        
        T * p = new T [ A . m ];
        // p=r;
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            p[k][0] = r[k][0];
            p[k][1] = r[k][1];
        }

        float rsold = 0;
        // rsold=r'*r;
        for ( int k(0)
            ; k < A . m
            ; ++k
            )
        {
            rsold += r[k][0] * r[k][0] + r[k][1] * r[k][1];
        }

        float rsnew , fact , alpha ;

        T * Ap = new T [ A . m ];
        int n_iter = A . m;
        if ( num_iter > 0 )
        {
            n_iter = std::min ( A . m , num_iter );
        }
        for ( int iter(0)
            ; iter < n_iter
            ; ++iter
            )
        {
            A ( p , Ap );
            // alpha=rsold/(p'*Ap);
            alpha = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                alpha += p[k][0] * Ap[k][0] + p[k][1] * Ap[k][1];
            }
            alpha = rsold / alpha;

            // x=x+alpha*p;
            // r=r-alpha*Ap;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                x[k][0] += alpha *  p[k][0];
                x[k][1] += alpha *  p[k][1];
                r[k][0] -= alpha * Ap[k][0];
                r[k][1] -= alpha * Ap[k][1];
            }

            // rsnew=r'*r;
            rsnew = 0;
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                rsnew += r[k][0] * r[k][0] + r[k][1] * r[k][1];
            }

            if ( rsnew < 1e-10 ) 
            {
                if ( rsnew < 0 )
                {
                    // operator not positive definitive
                    break;
                }
                break;
            }

            //p=r+(rsnew/rsold)*p;
            fact = (rsnew / rsold);
            for ( int k(0)
                ; k < A . m
                ; ++k
                )
            {
                p[k][0] = r[k][0] + fact * p[k][0];
                p[k][1] = r[k][1] + fact * p[k][1];
            }

            //rsold=rsnew;
            rsold = rsnew;

        }

        delete [] Ap;
        delete []  r;
        delete []  p;

    }

};

#endif


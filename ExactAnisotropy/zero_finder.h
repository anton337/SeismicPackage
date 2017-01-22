#ifndef ZERO_FINDER_H
#define ZERO_FINDER_H

template < typename T >
class ZeroFinder
{

    typedef T value_type;

public:
    value_type operator () ( value_type const & _s
                           , value_type const & _t
                           , value_type (*apply) ( value_type const & _x 
                                                 , void * _args
                                                 )
                           , void * args
                           )
    {
        value_type r,fr;
        int n, side=0;
        /* starting values at endpoints of interval */
        value_type  s = _s;
        value_type  t = _t;
        value_type fs = (*apply)(s,args);
        value_type ft = (*apply)(t,args);
        int m = 1000000;

        for (n = 0; n < m; n++)
        {
            r = (fs*t - ft*s) / (fs - ft);
            if (fabs(t-s) < 1e-8*fabs(t+s)) 
            {
                return r;
            }
            fr = (*apply)(r,args);
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
            return -1;
        }
    }

};

#endif


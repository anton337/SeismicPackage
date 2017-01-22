#ifndef VELOCITY_MODEL_H
#define VELOCITY_MODEL_H

#include <string.h>
#include <cstdlib>

namespace raytracing 
{


    template < typename T >
    int sgn ( T const & a )
    {
        return ( a > 0.0f ) ? 1 : -1 ;
    }
    
    template < typename T > 
    T norm ( T const & x
           , T const & y
           , T const & z
           )
    {
        return sqrtf ( x*x + y*y + z*z );
    }
    
    template < typename T >
    T dot ( T const & x1
          , T const & y1
          , T const & z1
          , T const & x2
          , T const & y2
          , T const & z2
          )
    {
        return x1*x2 + y1*y2 + z1*z2;
    }
    
    template < typename T >
    void cross ( T const & x1
               , T const & y1
               , T const & z1
               , T const & x2
               , T const & y2
               , T const & z2
               , T       & x3
               , T       & y3
               , T       & z3
               )
    {
        x3 = y1*z2 - y2*z1;
        y3 = x2*z1 - x1*z2;
        z3 = x1*y2 - y1*x2;
    }
    
    template < typename T >
    T cross_norm ( T const & x1
                 , T const & y1
                 , T const & z1
                 , T const & x2
                 , T const & y2
                 , T const & z2
                 )
    {
        T x , y , z;
        cross ( x1 , y1 , z1
              , x2 , y2 , z2
              , x  , y  , z
              );
        return sgn ( dot ( x , y , z
                         , x2, y2, z2
                         )
                   )
               * norm ( x , y , z );
    }
    
    template < typename T >
    T angle ( T const & x1
            , T const & y1
            , T const & z1
            , T const & x2
            , T const & y2
            , T const & z2
            )
    {
        return atan2f ( cross_norm ( x1 , y1 , z1
                                   , x2 , y2 , z2
                                   )
                      , dot ( x1 , y1 , z1
                            , x2 , y2 , z2
                            )
                      );
    }


};

template < typename T >
class velocity_model 
{
    typedef T     value_type;
    typedef int integer_type;
    enum TYPE { ISOTROPIC
              , ANISOTROPIC
              };
    TYPE type;
    int num_x;
    int num_y;
    int num_z;
    T scale_x;
    T scale_y;
    T scale_z;
    T const * vel;
    T const * epsilon;
    T const * delta;
    T const * nux;
    T const * nuy;
    T const * nuz;
public:
    velocity_model ( int       _num_x
                   , int       _num_y
                   , int       _num_z
                   , T         _scale_x
                   , T         _scale_y
                   , T         _scale_z
                   , T const * _vel
                   , T const * _epsilon = NULL
                   , T const * _delta   = NULL
                   , T const * _nux     = NULL
                   , T const * _nuy     = NULL
                   , T const * _nuz     = NULL
                   ) 
    : num_x   ( _num_x   )
    , num_y   ( _num_y   )
    , num_z   ( _num_z   )
    , scale_x ( _scale_x )
    , scale_y ( _scale_y )
    , scale_z ( _scale_z )
    , vel     ( _vel     )
    , epsilon ( _epsilon )
    , delta   ( _delta   )
    , nux     ( _nux     )
    , nuy     ( _nuy     )
    , nuz     ( _nuz     )
    {
        if ( epsilon == NULL 
          && delta   == NULL 
          && nux     == NULL 
          && nuy     == NULL 
          && nuz     == NULL 
           )
        {
            type = ISOTROPIC;
        }
        else
        {
            type = ANISOTROPIC;
        }
    }
    inline 
    int get_index ( T const & x
                  , T const & y
                  , T const & z
                  ) const
    {
        int x_ind ( x / scale_x );
        int y_ind ( y / scale_y );
        int z_ind ( z / scale_z );
        if ( x_ind <      0 ) x_ind = 0;
        if ( y_ind <      0 ) y_ind = 0;
        if ( z_ind <      0 ) z_ind = 0;
        if ( x_ind >= num_x ) x_ind = num_x-1;
        if ( y_ind >= num_y ) y_ind = num_y-1;
        if ( z_ind >= num_z ) z_ind = num_z-1;
        //std::cout << " ind:" << x_ind << " " << y_ind << " " << z_ind << std::endl;
        //return z_ind + num_z * ( y_ind + num_y * x_ind );
        return z_ind + num_z * ( y_ind + num_y * x_ind );
    }
    inline value_type interpolate ( value_type      const * v
                                  , value_type      const & x
                                  , value_type      const & y
                                  , value_type      const & z
                                  ) const
    {
        value_type wt_x ( ( x/scale_x - (integer_type)(x/scale_x) ) );
        value_type wt_y ( ( y/scale_y - (integer_type)(y/scale_y) ) );
        value_type wt_z ( ( z/scale_z - (integer_type)(z/scale_z) ) );
        value_type c000 ( v [ get_index ( x         , y         , z         ) ] );
        value_type c001 ( v [ get_index ( x         , y         , z+scale_z ) ] );
        value_type c010 ( v [ get_index ( x         , y+scale_y , z         ) ] );
        value_type c011 ( v [ get_index ( x         , y+scale_y , z+scale_z ) ] );
        value_type c100 ( v [ get_index ( x+scale_x , y         , z         ) ] );
        value_type c101 ( v [ get_index ( x+scale_x , y         , z+scale_z ) ] );
        value_type c110 ( v [ get_index ( x+scale_x , y+scale_y , z         ) ] );
        value_type c111 ( v [ get_index ( x+scale_x , y+scale_y , z+scale_z ) ] );
        value_type c00 ( c000*(1.0f-wt_z) + c001*wt_z );
        value_type c01 ( c010*(1.0f-wt_z) + c011*wt_z );
        value_type c10 ( c100*(1.0f-wt_z) + c101*wt_z );
        value_type c11 ( c110*(1.0f-wt_z) + c111*wt_z );
        value_type c0 ( c00*(1.0f-wt_y) + c01*wt_y );
        value_type c1 ( c10*(1.0f-wt_y) + c11*wt_y );
        return c0*(1.0f-wt_x) + c1*wt_x;
    }
    inline value_type interpolate_inv ( value_type      const * v
                                      , value_type      const * eps
                                      , value_type      const * del
                                      , value_type      const * nux
                                      , value_type      const * nuy
                                      , value_type      const * nuz
                                      , value_type      const & x
                                      , value_type      const & y
                                      , value_type      const & z
                                      , value_type      const & vx
                                      , value_type      const & vy
                                      , value_type      const & vz
                                      ) const
    {
        value_type wt_x ( ( x/scale_x - (integer_type)(x/scale_x) ) );
        value_type wt_y ( ( y/scale_y - (integer_type)(y/scale_y) ) );
        value_type wt_z ( ( z/scale_z - (integer_type)(z/scale_z) ) );

        int        i000 (         get_index ( x         , y         , z         )   );
        int        i001 (         get_index ( x         , y         , z+scale_z )   );
        int        i010 (         get_index ( x         , y+scale_y , z         )   );
        int        i011 (         get_index ( x         , y+scale_y , z+scale_z )   );
        int        i100 (         get_index ( x+scale_x , y         , z         )   );
        int        i101 (         get_index ( x+scale_x , y         , z+scale_z )   );
        int        i110 (         get_index ( x+scale_x , y+scale_y , z         )   );
        int        i111 (         get_index ( x+scale_x , y+scale_y , z+scale_z )   );

        value_type s000 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i000 ] , nuy [ i000 ] , nuz [ i000 ] ) ) );
        value_type s001 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i001 ] , nuy [ i001 ] , nuz [ i001 ] ) ) );
        value_type s010 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i010 ] , nuy [ i010 ] , nuz [ i010 ] ) ) );
        value_type s011 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i011 ] , nuy [ i011 ] , nuz [ i011 ] ) ) );
        value_type s100 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i100 ] , nuy [ i100 ] , nuz [ i100 ] ) ) );
        value_type s101 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i101 ] , nuy [ i101 ] , nuz [ i101 ] ) ) );
        value_type s110 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i110 ] , nuy [ i110 ] , nuz [ i110 ] ) ) );
        value_type s111 ( sinf ( raytracing :: angle ( vx , vy , vz , nux [ i111 ] , nuy [ i111 ] , nuz [ i111 ] ) ) );

        s000 *= s000;
        s001 *= s001;
        s010 *= s010;
        s011 *= s011;
        s100 *= s100;
        s101 *= s101;
        s110 *= s110;
        s111 *= s111;

        value_type t000 ( 1 - s000 );
        value_type t001 ( 1 - s001 );
        value_type t010 ( 1 - s010 );
        value_type t011 ( 1 - s011 );
        value_type t100 ( 1 - s100 );
        value_type t101 ( 1 - s101 );
        value_type t110 ( 1 - s110 );
        value_type t111 ( 1 - s111 );

        value_type c000 ( 1 / ( v [ i000 ] * ( 1 + del [ i000 ] * s000 * t000 + eps [ i000 ] * s000 * s000 ) ) );
        value_type c001 ( 1 / ( v [ i001 ] * ( 1 + del [ i001 ] * s001 * t001 + eps [ i001 ] * s001 * s001 ) ) );
        value_type c010 ( 1 / ( v [ i010 ] * ( 1 + del [ i010 ] * s010 * t010 + eps [ i010 ] * s010 * s010 ) ) );
        value_type c011 ( 1 / ( v [ i011 ] * ( 1 + del [ i011 ] * s011 * t011 + eps [ i011 ] * s011 * s011 ) ) );
        value_type c100 ( 1 / ( v [ i100 ] * ( 1 + del [ i100 ] * s100 * t100 + eps [ i100 ] * s100 * s100 ) ) );
        value_type c101 ( 1 / ( v [ i101 ] * ( 1 + del [ i101 ] * s101 * t101 + eps [ i101 ] * s101 * s101 ) ) );
        value_type c110 ( 1 / ( v [ i110 ] * ( 1 + del [ i110 ] * s110 * t110 + eps [ i110 ] * s110 * s110 ) ) );
        value_type c111 ( 1 / ( v [ i111 ] * ( 1 + del [ i111 ] * s111 * t111 + eps [ i111 ] * s111 * s111 ) ) );

        value_type c00 ( c000*(1.0f-wt_z) + c001*wt_z );
        value_type c01 ( c010*(1.0f-wt_z) + c011*wt_z );
        value_type c10 ( c100*(1.0f-wt_z) + c101*wt_z );
        value_type c11 ( c110*(1.0f-wt_z) + c111*wt_z );
        value_type c0 ( c00*(1.0f-wt_y) + c01*wt_y );
        value_type c1 ( c10*(1.0f-wt_y) + c11*wt_y );
        return c0*(1.0f-wt_x) + c1*wt_x;
    }
    inline value_type interpolate_inv ( value_type      const * v
                                      , value_type      const & x
                                      , value_type      const & y
                                      , value_type      const & z
                                      ) const
    {
        value_type wt_x ( ( x/scale_x - (integer_type)(x/scale_x) ) );
        value_type wt_y ( ( y/scale_y - (integer_type)(y/scale_y) ) );
        value_type wt_z ( ( z/scale_z - (integer_type)(z/scale_z) ) );
        value_type c000 ( 1 / v [ get_index ( x         , y         , z         ) ] );
        value_type c001 ( 1 / v [ get_index ( x         , y         , z+scale_z ) ] );
        value_type c010 ( 1 / v [ get_index ( x         , y+scale_y , z         ) ] );
        value_type c011 ( 1 / v [ get_index ( x         , y+scale_y , z+scale_z ) ] );
        value_type c100 ( 1 / v [ get_index ( x+scale_x , y         , z         ) ] );
        value_type c101 ( 1 / v [ get_index ( x+scale_x , y         , z+scale_z ) ] );
        value_type c110 ( 1 / v [ get_index ( x+scale_x , y+scale_y , z         ) ] );
        value_type c111 ( 1 / v [ get_index ( x+scale_x , y+scale_y , z+scale_z ) ] );
        value_type c00 ( c000*(1.0f-wt_z) + c001*wt_z );
        value_type c01 ( c010*(1.0f-wt_z) + c011*wt_z );
        value_type c10 ( c100*(1.0f-wt_z) + c101*wt_z );
        value_type c11 ( c110*(1.0f-wt_z) + c111*wt_z );
        value_type c0 ( c00*(1.0f-wt_y) + c01*wt_y );
        value_type c1 ( c10*(1.0f-wt_y) + c11*wt_y );
        return c0*(1.0f-wt_x) + c1*wt_x;
    }
    inline
    T operator () ( T const & x
                  , T const & y
                  , T const & z
                  ) const
    {
        return interpolate ( vel 
                           , x 
                           , y 
                           , z 
                           );
    }
    inline
    T operator () ( T const & x
                  , T const & y
                  , T const & z
                  , T const & px
                  , T const & py
                  , T const & pz
                  ) const
    {
        switch ( type )
        {
            case ISOTROPIC :
            {
                return 1 / interpolate_inv ( vel 
                                           , x 
                                           , y 
                                           , z 
                                           );
            }
            case ANISOTROPIC :
            {
                return 1 / interpolate_inv ( vel 
                                           , epsilon
                                           , delta
                                           , nux
                                           , nuy
                                           , nuz
                                           , x 
                                           , y 
                                           , z 
                                           , px
                                           , py
                                           , pz
                                           );
            }
        }
        return 0;
    }
    inline
    void get_anisotropic_parameters_and_gradients ( int const & num_iter
                                                  , T const   & step_size

                                                  , T const   & x
                                                  , T const   & y
                                                  , T const   & z
                                                  , T const   & vx
                                                  , T const   & vy
                                                  , T const   & vz

                                                  , T const   & h

                                                  , T         & v
                                                  , T         & vel_x
                                                  , T         & vel_y
                                                  , T         & vel_z

                                                  , T         & delta_v
                                                  , T         & delta_x
                                                  , T         & delta_y
                                                  , T         & delta_z

                                                  , T         & epsilon_v
                                                  , T         & epsilon_x
                                                  , T         & epsilon_y
                                                  , T         & epsilon_z

                                                  , T         & nux_v
                                                  , T         & nux_x
                                                  , T         & nux_y
                                                  , T         & nux_z

                                                  , T         & nuy_v
                                                  , T         & nuy_x
                                                  , T         & nuy_y
                                                  , T         & nuy_z

                                                  , T         & nuz_v
                                                  , T         & nuz_x
                                                  , T         & nuz_y
                                                  , T         & nuz_z

                                                  ) const
    {
        switch ( type )
        {
            case ISOTROPIC :
            {
                v = 1 / interpolate_inv ( vel , x , y , z );
                T vel_px ( 1 / interpolate_inv ( vel , x + h , y , z ) );
                T vel_mx ( 1 / interpolate_inv ( vel , x - h , y , z ) );
                T vel_py ( 1 / interpolate_inv ( vel , x , y + h , z ) );
                T vel_my ( 1 / interpolate_inv ( vel , x , y - h , z ) );
                T vel_pz ( 1 / interpolate_inv ( vel , x , y , z + h ) );
                T vel_mz ( 1 / interpolate_inv ( vel , x , y , z - h ) );
                vel_x = (vel_px - vel_mx) / (2*h);
                vel_y = (vel_py - vel_my) / (2*h);
                vel_z = (vel_pz - vel_mz) / (2*h);

                delta_v   = 0;
                delta_x   = 0;
                delta_y   = 0;
                delta_z   = 0;

                epsilon_v = 0;
                epsilon_x = 0;
                epsilon_y = 0;
                epsilon_z = 0;

                nux_v     = 0;
                nux_x     = 0;
                nux_y     = 0;
                nux_z     = 0;

                nuy_v     = 0;
                nuy_x     = 0;
                nuy_y     = 0;
                nuy_z     = 0;

                nuz_v     = 1;
                nuz_x     = 0;
                nuz_y     = 0;
                nuz_z     = 0;
                break;
            }
            case ANISOTROPIC :
            {
                v =        1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x , y , z     , vx , vy , vz );
                T vel_px ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x + h , y , z , vx , vy , vz ) );
                T vel_mx ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x - h , y , z , vx , vy , vz ) );
                T vel_py ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x , y + h , z , vx , vy , vz ) );
                T vel_my ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x , y - h , z , vx , vy , vz ) );
                T vel_pz ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x , y , z + h , vx , vy , vz ) );
                T vel_mz ( 1 / interpolate_inv ( vel , epsilon , delta , nux , nuy , nuz , x , y , z - h , vx , vy , vz ) );
                vel_x = (vel_px - vel_mx) / (2*h);
                vel_y = (vel_py - vel_my) / (2*h);
                vel_z = (vel_pz - vel_mz) / (2*h);

                T del_px ( interpolate ( delta , x + h , y , z ) );
                T del_mx ( interpolate ( delta , x - h , y , z ) );
                T del_py ( interpolate ( delta , x , y + h , z ) );
                T del_my ( interpolate ( delta , x , y - h , z ) );
                T del_pz ( interpolate ( delta , x , y , z + h ) );
                T del_mz ( interpolate ( delta , x , y , z - h ) );
                delta_v   = interpolate ( delta , x , y , z );
                delta_x   = (del_px - del_mx) / (2*h);
                delta_y   = (del_py - del_my) / (2*h);
                delta_z   = (del_pz - del_mz) / (2*h);

                T eps_px ( interpolate ( epsilon , x + h , y , z ) );
                T eps_mx ( interpolate ( epsilon , x - h , y , z ) );
                T eps_py ( interpolate ( epsilon , x , y + h , z ) );
                T eps_my ( interpolate ( epsilon , x , y - h , z ) );
                T eps_pz ( interpolate ( epsilon , x , y , z + h ) );
                T eps_mz ( interpolate ( epsilon , x , y , z - h ) );
                epsilon_v = interpolate ( epsilon , x , y , z );
                epsilon_x = (eps_px - eps_mx) / (2*h);
                epsilon_y = (eps_py - eps_my) / (2*h);
                epsilon_z = (eps_pz - eps_mz) / (2*h);

                nux_v     = interpolate ( nux , x , y , z );
                nuy_v     = interpolate ( nuy , x , y , z );
                nuz_v     = interpolate ( nuz , x , y , z );

                T nu_r = sqrtf ( nux_v 
                               * nux_v 
                               + nuy_v 
                               * nuy_v 
                               + nuz_v 
                               * nuz_v 
                               );

                nux_v /= nu_r;
                nuy_v /= nu_r;
                nuz_v /= nu_r;

                T nux_px    = interpolate ( nux , x + h , y , z );
                T nuy_px    = interpolate ( nuy , x + h , y , z );
                T nuz_px    = interpolate ( nuz , x + h , y , z );

                nu_r = sqrtf ( nux_px 
                             * nux_px 
                             + nuy_px 
                             * nuy_px 
                             + nuz_px 
                             * nuz_px 
                             );

                nux_px /= nu_r;
                nuy_px /= nu_r;
                nuz_px /= nu_r;

                T nux_mx    = interpolate ( nux , x - h , y , z );
                T nuy_mx    = interpolate ( nuy , x - h , y , z );
                T nuz_mx    = interpolate ( nuz , x - h , y , z );

                nu_r = sqrtf ( nux_mx 
                             * nux_mx 
                             + nuy_mx 
                             * nuy_mx 
                             + nuz_mx 
                             * nuz_mx 
                             );

                nux_mx /= nu_r;
                nuy_mx /= nu_r;
                nuz_mx /= nu_r;

                T nux_py    = interpolate ( nux , x , y + h , z );
                T nuy_py    = interpolate ( nuy , x , y + h , z );
                T nuz_py    = interpolate ( nuz , x , y + h , z );

                nu_r = sqrtf ( nux_py 
                             * nux_py 
                             + nuy_py 
                             * nuy_py 
                             + nuz_py 
                             * nuz_py 
                             );

                nux_py /= nu_r;
                nuy_py /= nu_r;
                nuz_py /= nu_r;

                T nux_my    = interpolate ( nux , x , y - h , z );
                T nuy_my    = interpolate ( nuy , x , y - h , z );
                T nuz_my    = interpolate ( nuz , x , y - h , z );

                nu_r = sqrtf ( nux_my 
                             * nux_my 
                             + nuy_my 
                             * nuy_my 
                             + nuz_my 
                             * nuz_my 
                             );

                nux_my /= nu_r;
                nuy_my /= nu_r;
                nuz_my /= nu_r;

                T nux_pz    = interpolate ( nux , x , y , z + h );
                T nuy_pz    = interpolate ( nuy , x , y , z + h );
                T nuz_pz    = interpolate ( nuz , x , y , z + h );

                nu_r = sqrtf ( nux_pz 
                             * nux_pz 
                             + nuy_pz 
                             * nuy_pz 
                             + nuz_pz 
                             * nuz_pz 
                             );

                nux_pz /= nu_r;
                nuy_pz /= nu_r;
                nuz_pz /= nu_r;

                T nux_mz    = interpolate ( nux , x , y , z - h );
                T nuy_mz    = interpolate ( nuy , x , y , z - h );
                T nuz_mz    = interpolate ( nuz , x , y , z - h );

                nu_r = sqrtf ( nux_mz 
                             * nux_mz 
                             + nuy_mz 
                             * nuy_mz 
                             + nuz_mz 
                             * nuz_mz 
                             );

                nux_mz /= nu_r;
                nuy_mz /= nu_r;
                nuz_mz /= nu_r;

                nux_x     = (nux_px - nux_mx) / ( 2 * h );
                nux_y     = (nux_py - nux_my) / ( 2 * h );
                nux_z     = (nux_pz - nux_mz) / ( 2 * h );
                nuy_x     = (nuy_px - nuy_mx) / ( 2 * h );
                nuy_y     = (nuy_py - nuy_my) / ( 2 * h );
                nuy_z     = (nuy_pz - nuy_mz) / ( 2 * h );
                nuz_x     = (nuz_px - nuz_mx) / ( 2 * h );
                nuz_y     = (nuz_py - nuz_my) / ( 2 * h );
                nuz_z     = (nuz_pz - nuz_mz) / ( 2 * h );

                break;
            }
        }
    }
};

#endif


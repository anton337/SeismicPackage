#ifndef VELOCITY_MODEL_H
#define VELOCITY_MODEL_H

#include <string.h>
#include <cstdlib>

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
    T const * nu_x;
    T const * nu_y;
    T const * nu_z;
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
                   , T const * _nu_x    = NULL
                   , T const * _nu_y    = NULL
                   , T const * _nu_z    = NULL
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
    , nu_x    ( _nu_x    )
    , nu_y    ( _nu_y    )
    , nu_z    ( _nu_z    )
    {
        if ( epsilon == NULL 
          && delta   == NULL 
          && nu_x    == NULL 
          && nu_y    == NULL 
          && nu_z    == NULL 
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
                std::cout << "anisotropic not implemented yet." << std::endl;
                exit(1);
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

                                                  , T         & delta
                                                  , T         & delta_x
                                                  , T         & delta_y
                                                  , T         & delta_z

                                                  , T         & epsilon
                                                  , T         & epsilon_x
                                                  , T         & epsilon_y
                                                  , T         & epsilon_z

                                                  , T         & nux
                                                  , T         & nux_x
                                                  , T         & nux_y
                                                  , T         & nux_z

                                                  , T         & nuy
                                                  , T         & nuy_x
                                                  , T         & nuy_y
                                                  , T         & nuy_z

                                                  , T         & nuz
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
                //std::cout << " vel:" << vel_pz << " " << vel_mz << std::endl;
                vel_x = (vel_px - vel_mx) / (2*h);
                vel_y = (vel_py - vel_my) / (2*h);
                vel_z = (vel_pz - vel_mz) / (2*h);

                delta     = 0;
                delta_x   = 0;
                delta_y   = 0;
                delta_z   = 0;

                epsilon   = 0;
                epsilon_x = 0;
                epsilon_y = 0;
                epsilon_z = 0;

                nux       = 0;
                nux_x     = 0;
                nux_y     = 0;
                nux_z     = 0;

                nuy       = 0;
                nuy_x     = 0;
                nuy_y     = 0;
                nuy_z     = 0;

                nuz       = 1;
                nuz_x     = 0;
                nuz_y     = 0;
                nuz_z     = 0;
                break;
            }
            case ANISOTROPIC :
            {
                std::cout << "anisotropic not implemented yet." << std::endl;
                exit(1);
                break;
            }
        }
    }
};

#endif


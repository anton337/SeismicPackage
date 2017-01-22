#ifndef ANISOTROPIC_RAY_TRACING_H
#define ANISOTROPIC_RAY_TRACING_H

#include <iostream>
#include <math.h>
#include <vector>
#include "velocity_model.h"

namespace raytracing
{


    typedef int integer_type;


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


   
    // need : velocity
    //        delta
    //        epsilon
    //        alpha
    //        beta
    template < typename value_type 
             >
    class tti_exact_ray_step_functor_rk4
    {

        typedef velocity_model < value_type > velocity_model_type ;

    public:

        inline
        void operator ()
        ( 
          value_type    const & increment
        , value_type    const & //x0
        , value_type    const & x1
        , value_type    const & vx_in
        , value_type          & vx_out
        , value_type          & ax_out
        , value_type    const & //y0
        , value_type    const & y1
        , value_type    const & vy_in
        , value_type          & vy_out
        , value_type          & ay_out
        , value_type    const & //z0
        , value_type    const & z1
        , value_type    const & vz_in
        , value_type          & vz_out
        , value_type          & az_out
        , value_type          & delta_travel_time
        , velocity_model_type & V
        )
        {

            // Eikonal equation
            // |Grad_x T| = |Grad T| = |p| = U
            //
            // Hamiltonian
            // H = 1/2 ( |Grad T|^2 / U^2 - 1 ) =
            //
            // solve equation
            // dH = dH/dx dx + dH / dp dp = 0
            //
            // dx/ds =   dH/dp
            // dp/ds = - dH/dx
            // dT/ds = p dH/dp
            //
            // dx/ds = vel^2 p
            // dp/ds = vel^3 |Grad T|^2 Grad U
            //
            // useful properties:
            // U = 1/V
            // Grad U = -(1/V^2) Grad V
            //
            // Theoretical solution
            // value_type dx_dt  =   vel * vel * px;
            // value_type dy_dt  =   vel * vel * py;
            // value_type dz_dt  =   vel * vel * pz;
            // value_type dpx_dt = - vel * vel * slow * p_2 * vel_x;
            // value_type dpy_dt = - vel * vel * slow * p_2 * vel_y;
            // value_type dpz_dt = - vel * vel * slow * p_2 * vel_z;
            
            // value_type nux (cos(beta)*sin(alpha));
            // value_type nuy (sin(beta)*sin(alpha));
            // value_type nuz (cos(alpha));

            value_type h ( 1e-2 );
            value_type vel , slow;

            value_type vel_x , vel_y , vel_z;

            value_type epsilon , delta ;

            value_type epsilon_x , epsilon_y , epsilon_z ;

            value_type delta_x , delta_y , delta_z ;

            value_type nux , nuy , nuz;

            value_type nux_x , nuy_x , nuz_x ;
            value_type nux_y , nuy_y , nuz_y ;
            value_type nux_z , nuy_z , nuz_z ;

            V . get_anisotropic_parameters_and_gradients ( 0
                                                         , 1e+6

                                                         , x1
                                                         , y1
                                                         , z1
                                                         , vx_in
                                                         , vy_in
                                                         , vz_in

                                                         , h

                                                         , vel
                                                         , vel_x
                                                         , vel_y
                                                         , vel_z

                                                         , delta
                                                         , delta_x
                                                         , delta_y
                                                         , delta_z

                                                         , epsilon
                                                         , epsilon_x
                                                         , epsilon_y
                                                         , epsilon_z

                                                         , nux
                                                         , nux_x
                                                         , nux_y
                                                         , nux_z

                                                         , nuy
                                                         , nuy_x
                                                         , nuy_y
                                                         , nuy_z

                                                         , nuz
                                                         , nuz_x
                                                         , nuz_y
                                                         , nuz_z

                                                         );
            slow = 1.0 / vel;

            value_type nu_2    = nux*nux + nuy*nuy + nuz*nuz;

            value_type px ( vx_in );
            value_type py ( vy_in );
            value_type pz ( vz_in );

            value_type p_2    = px*px + py*py + pz*pz;

            value_type p_dot_nu ( dot ( px , py , pz , nux , nuy , nuz ) );
            value_type cos_eta_2 ( p_dot_nu * p_dot_nu / ( p_2 * nu_2 ) );
            value_type z ( cos_eta_2 );
            value_type sin_eta_2 ( 1 - cos_eta_2 );
            value_type y ( sin_eta_2 );
            value_type f_eta ( V ( x1 , y1 , z1 , px , py , pz ) );
            value_type f_eta_2 ( f_eta * f_eta );

            value_type H_0 ( 0.5 * p_2 * vel /* vel */ ); // extra vel absobed into f_eta

            value_type dz_dx  = 2 * ( p_2 * nu_2 * p_dot_nu *     dot (  px ,  py ,  pz , nux_x , nuy_x , nuz_x ) 
                                    - p_dot_nu * p_dot_nu * p_2 * dot ( nux , nuy , nuz , nux_x , nuy_x , nuz_x ) 
                                    ) / ( p_2 * p_2 * nu_2 * nu_2 );
            value_type dz_dy  = 2 * ( p_2 * nu_2 * p_dot_nu *     dot (  px ,  py ,  pz , nux_y , nuy_y , nuz_y ) 
                                    - p_dot_nu * p_dot_nu * p_2 * dot ( nux , nuy , nuz , nux_y , nuy_y , nuz_y ) 
                                    ) / ( p_2 * p_2 * nu_2 * nu_2 );
            value_type dz_dz  = 2 * ( p_2 * nu_2 * p_dot_nu *     dot (  px ,  py ,  pz , nux_z , nuy_z , nuz_z ) 
                                    - p_dot_nu * p_dot_nu * p_2 * dot ( nux , nuy , nuz , nux_z , nuy_z , nuz_z ) 
                                    ) / ( p_2 * p_2 * nu_2 * nu_2 );

            value_type df_dx  = ( 2 * f_eta *(( delta * ( y - z ) - 2 * epsilon * y ) * dz_dx + epsilon_x * y * y + delta_x * y * z + ( + 2 * epsilon * y ) * slow * vel_x ) );
            value_type df_dy  = ( 2 * f_eta *(( delta * ( y - z ) - 2 * epsilon * y ) * dz_dy + epsilon_y * y * y + delta_y * y * z + ( + 2 * epsilon * y ) * slow * vel_y ) );
            value_type df_dz  = ( 2 * f_eta *(( delta * ( y - z ) - 2 * epsilon * y ) * dz_dz + epsilon_z * y * y + delta_z * y * z + ( + 2 * epsilon * y ) * slow * vel_z ) );

            value_type dx_dt  =   /*vel * vel */ px;
            value_type dy_dt  =   /*vel * vel */ py;
            value_type dz_dt  =   /*vel * vel */ pz;
            value_type dpx_dt = - /*vel * vel */ slow * ( p_2 * vel_x ) ;
            value_type dpy_dt = - /*vel * vel */ slow * ( p_2 * vel_y ) ;
            value_type dpz_dt = - /*vel * vel */ slow * ( p_2 * vel_z ) ;

            // vel * vel is absorbed into f_eta_2

                       vx_out = (  dx_dt * f_eta_2 );
                       vy_out = (  dy_dt * f_eta_2 );
                       vz_out = (  dz_dt * f_eta_2 );
                       ax_out = ( dpx_dt * f_eta_2 - H_0 * df_dx  );
                       ay_out = ( dpy_dt * f_eta_2 - H_0 * df_dy  );
                       az_out = ( dpz_dt * f_eta_2 - H_0 * df_dz  );

            delta_travel_time = - norm ( px , py , pz ) 
                                * norm ( vx_out
                                       , vy_out
                                       , vz_out
                                       );

        }



    };
    
   
    // need : velocity
    template < typename value_type 
             >
    class isotropic_exact_ray_step_functor_rk4
    {

        typedef velocity_model < value_type > velocity_model_type ;

    public:

        inline
        void operator ()
        ( 
          value_type    const & increment
        , value_type    const & //x0
        , value_type    const & x1
        , value_type    const & vx_in
        , value_type          & vx_out
        , value_type          & ax_out
        , value_type    const & //y0
        , value_type    const & y1
        , value_type    const & vy_in
        , value_type          & vy_out
        , value_type          & ay_out
        , value_type    const & //z0
        , value_type    const & z1
        , value_type    const & vz_in
        , value_type          & vz_out
        , value_type          & az_out
        , value_type          & delta_travel_time
        , velocity_model_type & V
        )
        {

            // Eikonal equation
            // |Grad_x T| = |Grad T| = |p| = U
            //
            // Hamiltonian
            // H = 1/2 ( |Grad T|^2 / U^2 - 1 ) =
            //
            // solve equation
            // dH = dH/dx dx + dH / dp dp = 0
            //
            // dx/ds =   dH/dp
            // dp/ds = - dH/dx
            // dT/ds = p dH/dp
            //
            // dx/ds = vel^2 p
            // dp/ds = vel^3 |Grad T|^2 Grad U
            //
            // useful properties:
            // U = 1/V
            // Grad U = -(1/V^2) Grad V
            //
            // Theoretical solution
            // value_type dx_dt  =   vel * vel * px;
            // value_type dy_dt  =   vel * vel * py;
            // value_type dz_dt  =   vel * vel * pz;
            // value_type dpx_dt = - vel * vel * slow * p_2 * vel_x;
            // value_type dpy_dt = - vel * vel * slow * p_2 * vel_y;
            // value_type dpz_dt = - vel * vel * slow * p_2 * vel_z;
            
            // value_type nux (cos(beta)*sin(alpha));
            // value_type nuy (sin(beta)*sin(alpha));
            // value_type nuz (cos(alpha));

            value_type h ( 1e-2 );
            value_type denom ( 1 / (2 * h) );
            value_type vel ( V ( x1 , y1 , z1 ) );
            value_type vel_2 ( vel * vel );
            value_type vel_x ( ( V ( x1 + h , y1 , z1 ) - V ( x1 - h , y1 , z1 ) ) * denom );
            value_type vel_y ( ( V ( x1 , y1 + h , z1 ) - V ( x1 , y1 - h , z1 ) ) * denom );
            value_type vel_z ( ( V ( x1 , y1 , z1 + h ) - V ( x1 , y1 , z1 - h ) ) * denom );
            value_type slow ( 1.0 / vel );

            value_type px ( vx_in );
            value_type py ( vy_in );
            value_type pz ( vz_in );

            value_type p_2    = px*px + py*py + pz*pz;

            value_type dx_dt  =   /*vel * vel */ px;
            value_type dy_dt  =   /*vel * vel */ py;
            value_type dz_dt  =   /*vel * vel */ pz;
            value_type dpx_dt = - /*vel * vel */ slow * ( p_2 * vel_x ) ;
            value_type dpy_dt = - /*vel * vel */ slow * ( p_2 * vel_y ) ;
            value_type dpz_dt = - /*vel * vel */ slow * ( p_2 * vel_z ) ;

            // vel * vel is absorbed into f_eta_2

                       vx_out = (  dx_dt * vel_2 );
                       vy_out = (  dy_dt * vel_2 );
                       vz_out = (  dz_dt * vel_2 );
                       ax_out = ( dpx_dt * vel_2 );
                       ay_out = ( dpy_dt * vel_2 );
                       az_out = ( dpz_dt * vel_2 );

            delta_travel_time = - norm ( px , py , pz ) 
                                * norm ( vx_out
                                       , vy_out
                                       , vz_out
                                       );

        }



    };
    


 


    // fourth order Runge Kutta (RK4)
    // y_n = y_{n-1} + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    // k1 = f(t_{n-1},y_{n-1})
    // k2 = f(t_{n-1}+h/2,y_{n-1}+h*k1/2)
    // k3 = f(t_{n-1}+h/2,y_{n-1}+h*k2/2)
    // k4 = f(t_{n-1}+h,y_{n-1}+h*k3)
    template < typename value_type
             , class derivative_type
             >
    class RK4_solver
    {

        typedef velocity_model < value_type > velocity_model_type ;

        derivative_type calculate_F;

    public:


        inline void RK4 ( value_type              const & h
                        , value_type                    & x     // make local copy, so the original does not change
                        , value_type                    & vx    // make local copy, so the original does not change 
                        , value_type                    & y     // make local copy, so the original does not change 
                        , value_type                    & vy    // make local copy, so the original does not change
                        , value_type                    & z     // make local copy, so the original does not change 
                        , value_type                    & vz    // make local copy, so the original does not change
                        , value_type                    & delta_travel_time
                        , velocity_model_type           & V
                        )
        {
            value_type  k_x_1 ( vx ) ,  k_x_2 ( vx ) ,  k_x_3 ( vx ) ,  k_x_4 ( vx ) ;
            value_type  k_y_1 ( vy ) ,  k_y_2 ( vy ) ,  k_y_3 ( vy ) ,  k_y_4 ( vy ) ;
            value_type  k_z_1 ( vz ) ,  k_z_2 ( vz ) ,  k_z_3 ( vz ) ,  k_z_4 ( vz ) ;
            value_type k_vx_1 ( 0.0) , k_vx_2 ( 0.0) , k_vx_3 ( 0.0) , k_vx_4 ( 0.0) ;
            value_type k_vy_1 ( 0.0) , k_vy_2 ( 0.0) , k_vy_3 ( 0.0) , k_vy_4 ( 0.0) ;
            value_type k_vz_1 ( 0.0) , k_vz_2 ( 0.0) , k_vz_3 ( 0.0) , k_vz_4 ( 0.0) ;
            value_type  k_t_1 ( 0.0) ,  k_t_2 ( 0.0) ,  k_t_3 ( 0.0) ,  k_t_4 ( 0.0) ;
        
            // k_1 is easy, just set it directly
            calculate_F ( h
                        , x
                        , x
                        , vx   
                        , k_x_1
                        , k_vx_1
                        , y
                        , y
                        , vy   
                        , k_y_1
                        , k_vy_1
                        , z
                        , z
                        , vz   
                        , k_z_1
                        , k_vz_1
                        , k_t_1
                        , V
                        );
        
            // calculate k_2
            calculate_F ( h
                        , x
                        , x     + 0.5f*h*k_x_1
                        , vx    + 0.5f*h*k_vx_1
                        , k_x_2
                        , k_vx_2
                        , y
                        , y     + 0.5f*h*k_y_1
                        , vy    + 0.5f*h*k_vy_1
                        , k_y_2
                        , k_vy_2
                        , z
                        , z     + 0.5f*h*k_z_1
                        , vz    + 0.5f*h*k_vz_1
                        , k_z_2
                        , k_vz_2
                        , k_t_2
                        , V
                        );
        
            // calculate k_3
            calculate_F ( h
                        , x
                        , x     + 0.5f*h*k_x_2
                        , vx    + 0.5f*h*k_vx_2
                        , k_x_3
                        , k_vx_3
                        , y
                        , y     + 0.5f*h*k_y_2
                        , vy    + 0.5f*h*k_vy_2
                        , k_y_3
                        , k_vy_3
                        , z
                        , z     + 0.5f*h*k_z_2
                        , vz    + 0.5f*h*k_vz_2
                        , k_z_3
                        , k_vz_3
                        , k_t_3
                        , V
                        );
        
            // calculate k_3
            calculate_F ( h
                        , x
                        , x     + h*k_x_3
                        , vx    + h*k_vx_3
                        , k_x_4
                        , k_vx_4
                        , y
                        , y     + h*k_y_3 
                        , vy    + h*k_vy_3
                        , k_y_4
                        , k_vy_4
                        , z
                        , z     + h*k_z_3
                        , vz    + h*k_vz_3
                        , k_z_4
                        , k_vz_4
                        , k_t_4
                        , V
                        );
        
            value_type dx  ( 0.166666666f * h * (  k_x_1 + 2.0f *  k_x_2 + 2.0f *  k_x_3 + k_x_4  ) );
            value_type dy  ( 0.166666666f * h * (  k_y_1 + 2.0f *  k_y_2 + 2.0f *  k_y_3 + k_y_4  ) );
            value_type dz  ( 0.166666666f * h * (  k_z_1 + 2.0f *  k_z_2 + 2.0f *  k_z_3 + k_z_4  ) );
            value_type dvx ( 0.166666666f * h * ( k_vx_1 + 2.0f * k_vx_2 + 2.0f * k_vx_3 + k_vx_4 ) );
            value_type dvy ( 0.166666666f * h * ( k_vy_1 + 2.0f * k_vy_2 + 2.0f * k_vy_3 + k_vy_4 ) );
            value_type dvz ( 0.166666666f * h * ( k_vz_1 + 2.0f * k_vz_2 + 2.0f * k_vz_3 + k_vz_4 ) );
            x +=  dx  ;
            y +=  dy  ;
            z +=  dz  ;
            vx += dvx ;
            vy += dvy ;
            vz += dvz ;
            delta_travel_time += 0.166666666f * h * ( k_t_1 + 2.0f * k_t_2 + 2.0f * k_t_3 + k_t_4 );


        }

    };

};


#endif





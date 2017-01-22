#include <GL/glut.h>
#include <iostream>
#include <math.h>
#include "zero_finder.h"

void init(void)
{
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_BLEND);

  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
    /* aspect ratio */ 1.0,
    /* Z near */ 0.1, /* Z far */ 1000.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, 4.75,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, +2.0);
  glRotatef(90, 1.0, 0.0, 0.0);
  //glRotatef(-20, 0.0, 0.0, 1.0);
}

void idle(void)
{
  glutPostRedisplay();
}

void drawStuff(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glutSwapBuffers();
}

void display(void)
{
    drawStuff();
}

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        default: break;
    }
}

template < typename T > struct model;

template < typename T >
struct ZeroArgs
{
    T x;
    T y;
    T z;
    T vx;
    T vy;
    T vz;
    model < T > V;
    T v1;
};



template < typename T >
T find_intersection ( T const & x
                    , void * args
                    )
{
    ZeroArgs < T > * p ( (ZeroArgs<T>*)args );
    return fabs ( p->v1 
                - p->V ( p->x + x*p->vx 
                       , p->y + x*p->vy
                       , p->z + x*p->vz
                       , p->vx
                       , p->vy
                       , p->vz
                       ) 
                ) 
           > 0.0001 
           ? 
           -1 
           : 
           1
           ;
}




template < typename T >
struct ZeroVelArgs
{
    T x1;
    T y1;
    T z1;
    T vx1;
    T vy1;
    T vz1;
    T th1;
    T xc;
    T yc;
    T zc;
    T x_normal;
    T y_normal;
    T z_normal;
    T th2;
    T vx2;
    T vy2;
    T vz2;
    model < T > V;
};



template < typename T >
T find_vel_intersection ( T const & x
                        , void * args
                        )
{
    ZeroVelArgs < T > * p ( (ZeroVelArgs<T>*)args );
    T nx;
    T ny;
    T nz;
    p -> V . get_normal ( p->x_normal , p->y_normal , p->z_normal , nx , ny , nz );
    if ( dot ( p->vx1 , p->vy1 , p->vz1 , nx , ny , nz ) < 0 )
    {
        nx *= -1;
        ny *= -1;
        nz *= -1;
    }
    T th1 ( angle ( p->vx1 , p->vy1 , p->vz1 , nx , ny , nz ) );
    p -> th1 = th1;
    T X ( sin(th1) / p->V ( p->x1 , p->y1 , p->z1 , p->vx1 , p->vy1 , p->vz1 ) );
    if ( fabs ( X ) < 1e-10 )
    {
        p -> th2 = th1;
        p -> vx2 = p -> vx1;
        p -> vy2 = p -> vy1;
        p -> vz2 = p -> vz1;
        return 0;
    }
    T _v ( norm ( p->vx1 , p->vy1 , p->vz1 ) );
    T _vx ( p->vx1 / _v );
    T _vy ( p->vy1 / _v );
    T _vz ( p->vz1 / _v );
    T proj ( dot ( _vx , _vy , _vz , nx , ny , nz ) );
    T _proj_x ( proj * nx );
    T _proj_y ( proj * ny );
    T _proj_z ( proj * nz );
    T _xx ( _vx - _proj_x );
    T _xy ( _vy - _proj_y );
    T _xz ( _vz - _proj_z );
    T vx2 ( _xx + nx * x );
    T vy2 ( _xy + ny * x );
    T vz2 ( _xz + nz * x );
    T th2 ( angle ( vx2 , vy2 , vz2 , nx , ny , nz ) );
    p -> th2 = th2;
    p -> vx2 = vx2;
    p -> vy2 = vy2;
    p -> vz2 = vz2;
    return sin ( th2 ) - X * p->V ( p->xc , p->yc , p->zc , vx2 , vy2 , vz2 );
}




template < typename T >
T find_vel_intersection_critical ( T const & x
                                 , void * args
                                 )
{
    ZeroVelArgs < T > * p ( (ZeroVelArgs<T>*)args );
    T nx;
    T ny;
    T nz;
    p -> V . get_normal ( p->x_normal , p->y_normal , p->z_normal , nx , ny , nz );
    if ( dot ( p->vx1 , p->vy1 , p->vz1 , nx , ny , nz ) < 0 )
    {
        nx *= -1;
        ny *= -1;
        nz *= -1;
    }
    T th1 ( angle ( p->vx1 , p->vy1 , p->vz1 , nx , ny , nz ) );
    p -> th1 = th1;
    T X ( sin(th1) / p->V ( p->x1 , p->y1 , p->z1 , p->vx1 , p->vy1 , p->vz1 ) );
    if ( fabs ( X ) < 1e-10 )
    {
        p -> th2 = th1;
        p -> vx2 = p -> vx1;
        p -> vy2 = p -> vy1;
        p -> vz2 = p -> vz1;
        return 0;
    }
    T _v ( norm ( p->vx1 , p->vy1 , p->vz1 ) );
    T _vx ( p->vx1 / _v );
    T _vy ( p->vy1 / _v );
    T _vz ( p->vz1 / _v );
    T proj ( dot ( _vx , _vy , _vz , nx , ny , nz ) );
    T _proj_x ( proj * nx );
    T _proj_y ( proj * ny );
    T _proj_z ( proj * nz );
    T _xx ( _vx - _proj_x );
    T _xy ( _vy - _proj_y );
    T _xz ( _vz - _proj_z );
    T vx2 ( _xx - nx * x );
    T vy2 ( _xy - ny * x );
    T vz2 ( _xz - nz * x );
    T th2 ( angle ( vx2 , vy2 , vz2 , nx , ny , nz ) );
    p -> th2 = th2;
    p -> vx2 = vx2;
    p -> vy2 = vy2;
    p -> vz2 = vz2;
    return sin ( th2 ) - X * p->V ( p->x1 , p->y1 , p->z1 , vx2 , vy2 , vz2 );
}

template < typename T >
struct model
{

    inline T norm( T const & x
                     , T const & y
                     , T const & z
                     ) const
    {
        return ( sqrt ( x*x + y*y + z*z ) );
    };
    
    inline int signum ( T n ) const
    {
        return (n>=0)?1:-1;
    }
        
    inline void cross_product ( T const & ax
                              , T const & ay
                              , T const & az
                              , T       &  x
                              , T       &  y
                              , T       &  z
                              ) const
    {
        T cx (x);
        T cy (y);
        T cz (z);
        x = ay*cz - az*cy;
        y = az*cx - ax*cz;
        z = ax*cy - ay*cx;
    }
    
    
    inline void cross_product ( T const & ax
                              , T const & ay
                              , T const & az
                              , T const & bx
                              , T const & by
                              , T const & bz
                              , T       &  x
                              , T       &  y
                              , T       &  z
                              ) const
    {
        x = ay*bz - az*by;
        y = az*bx - ax*bz;
        z = ax*by - ay*bx;
    }
    
    inline T vector_cos_angle( T const & pV1x
                                 , T const & pV1y
                                 , T const & pV1z
                                 , T const & pV2x
                                 , T const & pV2y
                                 , T const & pV2z
                                 ) const
    {  /*returns the cosine of the angle between two vectors */
        return (pV1x*pV2x+pV1y*pV2y+pV1z*pV2z)
              /(norm(pV1x,pV1y,pV1z)*norm(pV2x,pV2y,pV2z));
    }
    
    inline T vector_sin_angle( T const & pV1x
                                 , T const & pV1y
                                 , T const & pV1z
                                 , T const & pV2x
                                 , T const & pV2y
                                 , T const & pV2z
                                 ) const
    {  /*returns the sine of the angle between two vectors */
        T cx,cy,cz;
        cross_product(pV1x,pV1y,pV1z,pV2x,pV2y,pV2z,cx,cy,cz);
        return norm(cx,cy,cz)/(norm(pV1x,pV1y,pV1z)*norm(pV2x,pV2y,pV2z));
    }
        
    inline T vector_angle( T const & pV1x
                             , T const & pV1y
                             , T const & pV1z
                             , T const & pV2x
                             , T const & pV2y
                             , T const & pV2z
                             ) const
    {  /* returns the measure of the angle between two vectors */
        int sdp=signum(vector_cos_angle(pV1x,pV1y,pV1z,pV2x,pV2y,pV2z));
        return ((1-sdp)*0.5f*M_PI+sdp*asin((T)vector_sin_angle(pV1x,pV1y,pV1z,pV2x,pV2y,pV2z)));
    }

    T C33 ( T vp
              , T rho
              ) const
    {
        return vp*vp*rho;
    }
    T C44 ( T vs
              , T rho
              ) const
    {
        return vs*vs*rho;
    }
    T C11 ( T epsilon
              , T c33
              ) const
    {
        return (1+2*epsilon)*c33;
    }
    T C66 ( T gamma
              , T c44
              ) const
    {
        return (1+2*gamma)*c44;
    }
    T C13 ( T vp
              , T vs
              , T rho
              , T delta
              , T c44
              ) const
    {
        return rho * sqrt ( (vp*vp-vs*vs)*(vp*vp*(1+2*delta)-vs*vs) ) - c44;
    }
    T pow2 ( T a 
               ) const
    {
        return a*a;
    }
    T M ( T c11
            , T c13
            , T c33
            , T c44
            , T th 
            ) const
    {
        return pow2((c11-c44)*pow2(sin(th)) - (c33-c44)*pow2(cos(th))) + pow2(c13+c44)*pow2(sin(2*th));
    }
    T VP ( T x
             , T y
             , T z
             ) const
    {
        return 3000;
    }
    T VS ( T x
             , T y
             , T z
             ) const
    {
        return 2500;
    }
    T get_rho ( T x
                  , T y
                  , T z
                  ) const
    {
        return 1.0 / (3000*3000);
    }
    void get_nu ( T x
                , T y
                , T z
                , T & nux
                , T & nuy
                , T & nuz
                ) const
    {
        nux = 0;
        nuy = 0;
        nuz = 1;
    }
    void get_normal ( T x
                    , T y
                    , T z
                    , T & nx
                    , T & ny
                    , T & nz
                    ) const
    {
        nx = 0;
        ny = 0;
        nz = 1;
    }
    T get_epsilon ( T x
                      , T y
                      , T z
                      ) const
    {
        return 0.05;
    }
    T get_delta   ( T x
                      , T y
                      , T z
                      ) const
    {
        return 0.03;
    }
    T get_gamma   ( T x
                      , T y
                      , T z
                      ) const
    {
        return 0;
    }
    T Vqp ( T x
              , T y
              , T z
              , T px
              , T py
              , T pz
              ) const
    {
        T vp ( VP ( x , y , z  ) );
        T vs ( VS ( x , y , z  ) );
        T rho ( get_rho ( x , y , z ) );
        T nux ( 0.0f );
        T nuy ( 0.0f );
        T nuz ( 1.0f );
        get_nu ( x , y , z , nux , nuy , nuz );
        T epsilon ( get_epsilon ( x , y , z ) );
        T delta   ( get_delta   ( x , y , z ) );
        T gamma   ( get_gamma   ( x , y , z ) );
        T th ( vector_angle ( px , py , pz , nux , nuy , nuz ) );
        T c33 ( C33 ( vp , rho ) );
        T c44 ( C44 ( vs , rho ) );
        T c11 ( C11 ( epsilon , c33 ) );
        T c66 ( C66 ( gamma , c44 ) );
        T c13 ( C13 ( vp , vs , rho , delta , c44 ) );
        T m ( M ( c11 , c13 , c33 , c44 , th ) );
        return sqrt ( (c11*pow2(sin(th)) + c33*pow2(cos(th))+c44+sqrt(m)) / (2*rho) );
    }
    T Vqs ( T x
              , T y
              , T z
              , T px
              , T py
              , T pz
              ) const
    {
        T vp ( VP ( x , y , z  ) );
        T vs ( VS ( x , y , z  ) );
        T rho ( get_rho ( x , y , z ) );
        T nux ( 0.0f );
        T nuy ( 0.0f );
        T nuz ( 1.0f );
        get_nu ( x , y , z , nux , nuy , nuz );
        T epsilon ( get_epsilon ( x , y , z ) );
        T delta   ( get_delta   ( x , y , z ) );
        T gamma   ( get_gamma   ( x , y , z ) );
        T th ( vector_angle ( px , py , pz , nux , nuy , nuz ) );
        T c33 ( C33 ( vp , rho ) );
        T c44 ( C44 ( vs , rho ) );
        T c11 ( C11 ( epsilon , c33 ) );
        T c66 ( C66 ( gamma , c44 ) );
        T c13 ( C13 ( vp , vs , rho , delta , c44 ) );
        T m ( M ( c11 , c13 , c33 , c44 , th ) );
        return sqrt ( (c11*pow2(sin(th)) + c33*pow2(cos(th))+c44-sqrt(m)) / (2*rho) );
    }
    T Vs ( T x
             , T y
             , T z
             , T px
             , T py
             , T pz
             ) const
    {
        T vp ( VP ( x , y , z  ) );
        T vs ( VS ( x , y , z  ) );
        T rho ( get_rho ( x , y , z ) );
        T nux ( 0.0f );
        T nuy ( 0.0f );
        T nuz ( 1.0f );
        get_nu ( x , y , z , nux , nuy , nuz );
        T epsilon ( get_epsilon ( x , y , z ) );
        T delta   ( get_delta   ( x , y , z ) );
        T gamma   ( get_gamma   ( x , y , z ) );
        T th ( vector_angle ( px , py , pz , nux , nuy , nuz ) );
        T c33 ( C33 ( vp , rho ) );
        T c44 ( C44 ( vs , rho ) );
        T c11 ( C11 ( epsilon , c33 ) );
        T c66 ( C66 ( gamma , c44 ) );
        T c13 ( C13 ( vp , vs , rho , delta , c44 ) );
        T m ( M ( c11 , c13 , c33 , c44 , th ) );
        return sqrt ( (c66*pow2(sin(th)) + c44*pow2(cos(th))) / (rho) );
    }
    void get_differential_update ( T const & x
                                 , T const & y
                                 , T const & z
                                 , T const & px
                                 , T const & py
                                 , T const & pz
                                 , T const & t
                                 , T       & dx
                                 , T       & dy
                                 , T       & dz
                                 , T       & dpx
                                 , T       & dpy
                                 , T       & dpz
                                 , T       & dt
                                 , T         h = 1e-3
                                 ) const
    {
        T V ( Vqp ( x , y , z , px , py , pz ) );
        T U ( 1.0 / V );
        dx = px * V * V;
        dy = py * V * V;
        dz = pz * V * V;
        T p_2 ( px*px + py*py + pz*pz );
        T dUx ( ( Vqp ( x+h , y   , z   , px , py , pz ) - Vqp ( x-h , y   , z   , px , py , pz ) ) / (2*h) );
        T dUy ( ( Vqp ( x   , y+h , z   , px , py , pz ) - Vqp ( x   , y-h , z   , px , py , pz ) ) / (2*h) );
        T dUz ( ( Vqp ( x   , y   , z+h , px , py , pz ) - Vqp ( x   , y   , z-h , px , py , pz ) ) / (2*h) );
        dpx = p_2 * V * V * V * dUx;
        dpy = p_2 * V * V * V * dUy;
        dpz = p_2 * V * V * V * dUz;
        dt = sqrt ( dx * px + dy * py + dz * pz );
    }
    void get_exact_update ( T const & x
                          , T const & y
                          , T const & z
                          , T const & vx
                          , T const & vy
                          , T const & vz
                          , T const & t
                          , T       & fx
                          , T       & fy
                          , T       & fz
                          , T       & fvx
                          , T       & fvy
                          , T       & fvz
                          , T       & ft
                          , T         ds = 1e-3
                          ) const
    {
        T v1 ( Vqp ( x , y , z , vx , vy , vz ) );
        T v ( norm ( vx , vy , vz ) );
        T _vx ( v1 * ds * vx / v );
        T _vy ( v1 * ds * vy / v );
        T _vz ( v1 * ds * vz / v );
        T _x ( x + _vx );
        T _y ( y + _vy );
        T _z ( z + _vz );
        T v2 ( Vqp (_x ,_y ,_z , vx , vy , vz ) );
        if ( fabs ( v1 - v2 ) < 0.001 )
        {
            // in the same region
            T dt = ds;
            if ( t >= dt )
            {
                x = _x;
                y = _y;
                z = _z;
                t -= dt;
            }
            else
            {
                x += v1 * t * vx / v;
                y += v1 * t * vy / v;
                z += v1 * t * vz / v;
                t = 0;
            }
        }
        else
        {
            // crossing region boundary
            // find exact crossing position
            ZeroFinder < T > zero_finder;
            ZeroArgs < T > zero_args;
            zero_args . x = x;
            zero_args . y = y;
            zero_args . z = z;
            zero_args . vx = _vx;
            zero_args . vy = _vy;
            zero_args . vz = _vz;
            zero_args . V = this;
            zero_args . v1 = v1;
            T zero ( zero_finder ( 0 , 1 , &find_intersection , &zero_args ) );
            T __x ( x + zero * _vx );
            T __y ( y + zero * _vy );
            T __z ( z + zero * _vz );
            T ___x ( x + (zero+0.001) * _vx );
            T ___y ( y + (zero+0.001) * _vy );
            T ___z ( z + (zero+0.001) * _vz );
            T x_normal ( (v1 < v2) ? x + (zero+0.001) * _vx : x + (zero-0.001) * _vx );
            T y_normal ( (v1 < v2) ? y + (zero+0.001) * _vy : y + (zero-0.001) * _vy );
            T z_normal ( (v1 < v2) ? z + (zero+0.001) * _vz : z + (zero-0.001) * _vz );
            // ray constant X = sin(th1) / v1(eta1) = sin(th2) / v2(eta2)
            ZeroVelArgs < T > zero_vel_args;
            zero_vel_args . x1 = x;
            zero_vel_args . y1 = y;
            zero_vel_args . z1 = z;
            zero_vel_args . vx1 = _vx;
            zero_vel_args . vy1 = _vy;
            zero_vel_args . vz1 = _vz;
            zero_vel_args . xc = ___x;
            zero_vel_args . yc = ___y;
            zero_vel_args . zc = ___z;
            zero_vel_args . x_normal = x_normal;
            zero_vel_args . y_normal = y_normal;
            zero_vel_args . z_normal = z_normal;
            zero_vel_args . V = this;
            T vel_zero ( zero_finder ( 0.0 , 10.0 , &find_vel_intersection , &zero_vel_args ) );
            if ( fabs ( vel_zero + 1 ) < 1e-5 )
            {
                vel_zero = ( zero_finder ( 0.0 , 10.0 , &find_vel_intersection_critical , &zero_vel_args ) );
                std::cout << vel_zero << std::endl;
            }
            x = __x;
            y = __y;
            z = __z;
            T vx2 = zero_vel_args . vx2;
            T vy2 = zero_vel_args . vy2;
            T vz2 = zero_vel_args . vz2;
            T v2 ( V ( ___x , ___y , ___z , vx2 , vy2 , vz2 ) );
            v = norm ( vx2 , vy2 , vz2 );
            vx2 /= v;
            vy2 /= v;
            vz2 /= v;
            x += (1-zero)*ds*vx2*v2;
            y += (1-zero)*ds*vy2*v2;
            z += (1-zero)*ds*vz2*v2;
            t -= ds;
            vx = vx2;
            vy = vy2;
            vz = vz2;
        }

    }
};

int main(int argc,char * argv[])
{
    float x   = 0;
    float y   = 0;
    float z   = 0;
    float px  = 0;
    float py  = 0;
    float pz  = 1;
    float t   = 0;
    float dx ;
    float dy ;
    float dz ;
    float dpx;
    float dpy;
    float dpz;
    float dt ;
    model<float> m;
    m.get_differential_update ( x
                              , y
                              , z
                              , px
                              , py
                              , pz
                              , t
                              , dx
                              , dy
                              , dz
                              , dpx
                              , dpy
                              , dpz
                              , dt
                              );
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(1000,1000);
    glutCreateWindow("machine learning");
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    init();
    glutMainLoop();

    return 0;
}


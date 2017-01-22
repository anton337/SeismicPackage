#ifndef VIEWER_H
#define VIEWER_H

#include <boost/thread.hpp>
#include <GL/glut.h>

struct vect
{
    float vx;
    float vy;
    float vz;
    vect()
    {
        vx = 0;
        vy = 0;
        vz = 0;
    }
};

extern boost::thread       * show ( int argc , char ** argv );

extern void set_timer ( float _timer );
extern void set_num_x ( int _num_x );
extern void set_num_z ( int _num_z );
extern void set_data_ptr ( float * _data );
extern void set_wave_ptr ( float * _data );
extern void set_reverse_ptr ( float * _data );
extern void set_scalar ( float _scalar );
extern void set_wave_scalar ( float _scalar );
extern void set_reverse_scalar ( float _scalar );
extern void set_data_vec_ptr ( vect * _data );

#endif


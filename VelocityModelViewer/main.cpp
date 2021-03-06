#include <GL/glut.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include <deque>
#include "sep_reader.h"
#include "sep_writer.h"

bool show_anisotropy_axis = false;

bool shift_pressed = false;

bool ctrl_pressed = false;

bool layer = true;

int type_index = 1;

int pen_radius = 20;

int num_batch;

int width  = 1000;
int height = 1000;

int num_x;
int num_t;

int delta_slice = 1;

int timer = 0;

int data_offset = 0;

float scalar = 0;

float gain = 1;

float * velocity_data     = NULL;
float * velocity_s_data   = NULL;
float * epsilon_data      = NULL;
float * delta_data        = NULL;
float * nu_x_data         = NULL;
float * nu_y_data         = NULL;
float * nu_z_data         = NULL;

float * velocity_data_2   = NULL;
float * velocity_s_data_2 = NULL;
float * epsilon_data_2    = NULL;
float * delta_data_2      = NULL;
float * nu_x_data_2       = NULL;
float * nu_y_data_2       = NULL;
float * nu_z_data_2       = NULL;

std::string horizon_in_file     ;
std::string horizon_out_file    ;

std::string script_file         ;

std::string velocity_out_file   ;
std::string velocity_s_out_file ;
std::string epsilon_out_file    ;
std::string delta_out_file      ;
std::string nu_x_out_file       ;
std::string nu_y_out_file       ;
std::string nu_z_out_file       ;

std::string velocity_in_file    ;
std::string velocity_s_in_file  ;
std::string epsilon_in_file     ;
std::string delta_in_file       ;
std::string nu_x_in_file        ;
std::string nu_y_in_file        ;
std::string nu_z_in_file        ;

bool * toggle = NULL;

int * type = NULL;

bool * visited = NULL;

int selection_index = -1;

bool move_pt = false;

bool show_multiples = false;

struct spectral_tools_window
{

};

struct hyperbolic_event
{
    float x_apex;
    float y_apex;
    float velocity;
    float coefficient;
};

float velocity = 1500;

float shift = 1;

float fact = 1000;

float epsilon = 0;

float delta = 0;

float vel_p = 0;

float vel_s = 0;

std::vector < hyperbolic_event > event(1);

struct point
{
    float x_pos;
    float y_pos;
    point ( float _x_pos
          , float _y_pos
          )
    : x_pos ( _x_pos )
    , y_pos ( _y_pos )
    {

    }
    point ()
    : x_pos ( 0 )
    , y_pos ( 0 )
    {

    }
};

std::vector < point > selection;

std::vector < point > mig_selection;

std::vector < point > input_horizon;

std::vector < point > output_horizon;

// horizon file has format:
// x0 y0
// x1 y1
// x2 y2
// ...
void write_horizon ( std::string file_name , std::vector < point > const & h )
{
    std::string line;
    std::ofstream myfile( file_name . c_str () );
    if ( myfile . is_open () )
    {
        for ( std::size_t k(0)
            ; k < h.size()
            ; ++k
            )
        {
            myfile << num_batch * h[k].x_pos / width << " " << num_t * h[k].y_pos / height << std::endl;
        }
    }
    else
    {
        std::cout << "unable to open file " << file_name << std::endl;
        exit(1);
    }
}

void print_selection ()
{
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    for ( std::size_t k(0)
        ; k < selection . size ()
        ; ++k
        )
    {
        std::cout << num_batch * selection[k].x_pos / width << " " << num_t * selection[k].y_pos / height << "     ";
    }
    std::cout << std::endl;
    input_horizon = selection;
    write_horizon ( horizon_in_file , input_horizon );
}

float x_pos;
float y_pos;

bool mute_toggle ( std::vector < point > const & sel
                 , int x_pos 
                 , int y_pos
                 )
{
    float m;
    int count = 0;
    for ( std::size_t k(0)
        ; k + 1 < sel . size ()
        ; ++k
        )
    {
        if ( x_pos > sel[k].x_pos && x_pos <= sel[k+1].x_pos )
        {
            m = (sel[k+1].y_pos - sel[k].y_pos) / (sel[k+1].x_pos - sel[k].x_pos);
            if ( y_pos < sel[k].y_pos + m * ( x_pos - sel[k].x_pos ) )
            {
                count += 1;
            }
        }
        if ( x_pos > sel[k+1].x_pos && x_pos <= sel[k].x_pos )
        {
            m = (sel[k+1].y_pos - sel[k].y_pos) / (sel[k+1].x_pos - sel[k].x_pos);
            if ( y_pos < sel[k].y_pos + m * ( x_pos - sel[k].x_pos ) )
            {
                count += 1;
            }
        }
    }
    return ( count % 2 == 0 );
}

void calculate_toggle ()
{
    for ( int k(0)
        , x(0)
        ; x < num_batch
        ; ++x
        )
    for ( int t(0)
        ; t < num_t
        ; ++t
        , ++k
        )
    {
        toggle[k] = ( mute_toggle ( selection , width * (float)x / num_batch , height * (float)t / num_t ) );
    }
}

void drawString (void * font, char const *s, float x, float y, float z)
{
     unsigned int i;
     glRasterPos3f(x, y, z);
     for (i = 0; i < strlen (s); i++)
     {
         glutBitmapCharacter (font, s[i]);
     }
}

void drawStuff(void)
{
    {
        glColor3f ( 1 , 1 , 1 );
        {
            std::stringstream str;
            str << "vel_p:" << vel_p;
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.075 );
        }
        {
            std::stringstream str;
            str << "vel_s:" << vel_s;
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.125 );
        }
        {
            std::stringstream str;
            str << "epsilon:" << epsilon;
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.175 );
        }
        {
            std::stringstream str;
            str << "delta:" << delta;
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.225 );
        }
        if ( layer )
        {
            std::stringstream str;
            str << "forward";
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.275 );
        }
        else
        {
            std::stringstream str;
            str << "reverse";
            drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), - 1 + 0.075 , 0.01 , - 1 + 0.275 );
        }
    }
    if ( !shift_pressed && !ctrl_pressed )
    {
        glColor3f ( 1 , 0 , 0 );
        glBegin(GL_LINES);
        glVertex3f( 2 * x_pos / width - 1 - 0.050 , 0 , 2 * y_pos / height - 1         );
        glVertex3f( 2 * x_pos / width - 1 + 0.050 , 0 , 2 * y_pos / height - 1         );
        glVertex3f( 2 * x_pos / width - 1         , 0 , 2 * y_pos / height - 1 - 0.050 );
        glVertex3f( 2 * x_pos / width - 1         , 0 , 2 * y_pos / height - 1 + 0.050 );
        float nux , nuz;
        if ( layer )
        {
            nux = 200*nu_x_data   [ (int)((x_pos/width)*num_batch) * num_t + (int)((y_pos/height)*num_t) ] / width  ;
            nuz = 200*nu_z_data   [ (int)((x_pos/width)*num_batch) * num_t + (int)((y_pos/height)*num_t) ] / height ;
        }
        else
        {
            nux = 200*nu_x_data_2 [ (int)((x_pos/width)*num_batch) * num_t + (int)((y_pos/height)*num_t) ] / width  ;
            nuz = 200*nu_z_data_2 [ (int)((x_pos/width)*num_batch) * num_t + (int)((y_pos/height)*num_t) ] / height ;
        }
        glVertex3f( 2 * x_pos / width - 1 + nux   , 0 , 2 * y_pos / height - 1 + nuz   );
        glVertex3f( 2 * x_pos / width - 1 - nux   , 0 , 2 * y_pos / height - 1 - nuz   );
        glEnd();
        {
            int ind ( (int)((y_pos/height)*num_t) + (int)((x_pos/width)*num_batch)*num_t );
            if ( ind >= 0 && ind < num_t*num_batch )
            {
                glColor3f ( 1 , 1 , 1 );
                {
                    std::stringstream str;
                    str << "vel_p:" << ((layer)?velocity_data[ind]:velocity_data_2[ind]);
                    drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), 2 * x_pos / width - 1 + 0.075 , 0.01 , 2 * y_pos / height - 1 + 0.075 );
                }
                {
                    std::stringstream str;
                    str << "vel_s:" << ((layer)?velocity_s_data[ind]:velocity_s_data_2[ind]);
                    drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), 2 * x_pos / width - 1 + 0.075 , 0.01 , 2 * y_pos / height - 1 + 0.125 );
                }
                {
                    std::stringstream str;
                    str << "epsilon:" << ((layer)?epsilon_data[ind]:epsilon_data_2[ind]);
                    drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), 2 * x_pos / width - 1 + 0.075 , 0.01 , 2 * y_pos / height - 1 + 0.175 );
                }
                {
                    std::stringstream str;
                    str << "delta:" << ((layer)?delta_data[ind]:delta_data_2[ind]);
                    drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), 2 * x_pos / width - 1 + 0.075 , 0.01 , 2 * y_pos / height - 1 + 0.225 );
                }
                // {
                //     std::stringstream str;
                //     str << "type:" << type[ind];
                //     drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), 2 * x_pos / width - 1 + 0.075 , 0.01 , 2 * y_pos / height - 1 + 0.275 );
                // }
            }
        }
    }
    else
    {
        glColor3f ( 1 , 0 , 0 );
        glBegin(GL_LINES);
        for ( float th(0); th < 2*M_PI; th += 0.1 )
        {
            glVertex3f( 2 * (x_pos+pen_radius*cos(th    )) / width - 1 , 0 , 2 * (y_pos+pen_radius*sin(th    )) / height - 1 );
            glVertex3f( 2 * (x_pos+pen_radius*cos(th+0.1)) / width - 1 , 0 , 2 * (y_pos+pen_radius*sin(th+0.1)) / height - 1 );
        }
        glEnd();
    }
    if ( show_multiples )
    {
    glColor3f ( 1 , 0 , 0 );
    glPointSize(3);
    glBegin(GL_POINTS);
    for ( std::size_t k(0)
        ; k < event.size()
        ; ++k
        )
    {
        float x_val;
        float y_val;
        for ( int mult(1)
            ; mult < 10
            ; ++mult
            )
        {
            for ( int x(-500)
                ; x < 500
                ; x += 5 
                )
            {
                x_val = event[k].x_apex + x;
                y_val = shift + sqrt ( (event[k].y_apex - shift)*(event[k].y_apex - shift)*mult*mult + fact*fact * (x_val-event[k].x_apex)*(x_val-event[k].x_apex) / (event[k].velocity*event[k].velocity) );
                x_val /= width;
                y_val /= height;
                glVertex3f(2*x_val-1,0,2*y_val-1);
            }
        }
    }
    glEnd();
    }
    glColor3f ( 0 , 1 , 0 );
    glBegin(GL_LINES);
    if ( move_pt )
    {
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 - 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 - 0.030 );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 - 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 + 0.030 );

        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 - 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 + 0.030 );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 + 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 + 0.030 );

        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 + 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 + 0.030 );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 + 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 - 0.030 );

        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 + 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 - 0.030 );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 - 0.030 , 0 , 2 * selection[selection_index] . y_pos / height - 1 - 0.030 );
    }
    if ( selection_index >= 0 && selection_index < (int)selection . size () )
    {
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 - 0.020 , 0 , 2 * selection[selection_index] . y_pos / height - 1         );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1 + 0.020 , 0 , 2 * selection[selection_index] . y_pos / height - 1         );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1         , 0 , 2 * selection[selection_index] . y_pos / height - 1 - 0.020 );
        glVertex3f( 2 * selection[selection_index] . x_pos / width - 1         , 0 , 2 * selection[selection_index] . y_pos / height - 1 + 0.020 );
        if ( selection_index+1 < (int)selection.size() )
        {
            glVertex3f( 2 * selection[selection_index+1] . x_pos / width - 1 - 0.020 , 0 , 2 * selection[selection_index+1] . y_pos / height - 1         );
            glVertex3f( 2 * selection[selection_index+1] . x_pos / width - 1 + 0.020 , 0 , 2 * selection[selection_index+1] . y_pos / height - 1         );
            glVertex3f( 2 * selection[selection_index+1] . x_pos / width - 1         , 0 , 2 * selection[selection_index+1] . y_pos / height - 1 - 0.020 );
            glVertex3f( 2 * selection[selection_index+1] . x_pos / width - 1         , 0 , 2 * selection[selection_index+1] . y_pos / height - 1 + 0.020 );
        }
    }
    for ( std::size_t k(0)
        ; k < selection . size ()
        ; ++k
        )
    {
        glVertex3f( 2 * selection[k] . x_pos / width - 1 - 0.020 , 0 , 2 * selection[k] . y_pos / height - 1 - 0.020 );
        glVertex3f( 2 * selection[k] . x_pos / width - 1 + 0.020 , 0 , 2 * selection[k] . y_pos / height - 1 + 0.020 );
        glVertex3f( 2 * selection[k] . x_pos / width - 1 + 0.020 , 0 , 2 * selection[k] . y_pos / height - 1 - 0.020 );
        glVertex3f( 2 * selection[k] . x_pos / width - 1 - 0.020 , 0 , 2 * selection[k] . y_pos / height - 1 + 0.020 );
        if ( k + 1 < selection . size () )
        {
            glVertex3f( 2 * selection[k  ] . x_pos / width - 1 , 0 , 2 * selection[k  ] . y_pos / height - 1 );
            glVertex3f( 2 * selection[k+1] . x_pos / width - 1 , 0 , 2 * selection[k+1] . y_pos / height - 1 );
        }
    }
    glEnd();
    glColor3f ( 0 , 0 , 1 );
    glBegin(GL_POINTS);
    for ( std::size_t k(0)
        ; k < mig_selection . size ()
        ; ++k
        )
    {
        if ( k + 1 < mig_selection . size () )
        {
            glVertex3f( 2 * mig_selection[k  ] . x_pos / width - 1 , 0 , 2 * mig_selection[k  ] . y_pos / height - 1 );
            glVertex3f( 2 * mig_selection[k+1] . x_pos / width - 1 , 0 , 2 * mig_selection[k+1] . y_pos / height - 1 );
        }
    }
    glEnd();
    if ( show_anisotropy_axis )
    {
        glColor3f(0,0,0);
        glBegin(GL_LINES);
        float nux , nuz;
        for ( int x(0)
            ; x < num_batch
            ; x += 15 
            )
        for ( int t(0)
            ; t < num_t
            ; t += 15 
            )
        {
            if ( layer )
            {
                nux = 3*nu_x_data [ x * num_t + t ]*num_batch/num_t;
                nuz = 3*nu_z_data [ x * num_t + t ];
            }
            else
            {
                nux = 3*nu_x_data_2 [ x * num_t + t ]*num_batch/num_t;
                nuz = 3*nu_z_data_2 [ x * num_t + t ];
            }
            glVertex3f 
            ( -1 + 2.0f * ((float)x+nux) / num_batch
            , 0
            , -1 + 2.0f * ((float)t+nuz) / num_t
            );
            glVertex3f 
            ( -1 + 2.0f * ((float)x-nux) / num_batch
            , 0
            , -1 + 2.0f * ((float)t-nuz) / num_t
            );
        }
        glEnd();
    }
    float val;
    float val_e;
    float val_d;
    glPointSize(2*std::max((float)height/num_t,(float)width/num_batch));
    glBegin(GL_POINTS);
    for ( int k(0)
        , x(0)
        ; x < num_batch
        ; ++x
        )
    for ( int t(0)
        ; t < num_t
        ; ++t
        , ++k
        )
    {
        if ( layer )
        {
            val = 0.1f + velocity_data[k+data_offset] * scalar * gain + velocity_s_data[k+data_offset] * scalar * gain;
            val_e = 0.5f + epsilon_data[k+data_offset] * gain;
            val_d = 0.5f + delta_data[k+data_offset] * gain;
        }
        else
        {
            val = 0.1f + velocity_data_2[k+data_offset] * scalar * gain + velocity_s_data_2[k+data_offset] * scalar * gain;
            val_e = 0.5f + epsilon_data_2[k+data_offset] * gain;
            val_d = 0.5f + delta_data_2[k+data_offset] * gain;
        }
        if ( toggle[k] )
        glColor3f
        ( val
        , val_e
        , val_d
        );
        else
        glColor3f
        ( val_e
        , val_d
        , val
        );
        glVertex3f 
        ( -1 + 2.0f * x / num_batch
        , 0
        , -1 + 2.0f * t / num_t
        );
    }
    glEnd();
    timer++;
    if ( timer > 2 )
    {
        timer = 0;
        data_offset += num_t * num_batch * delta_slice;
        if ( data_offset >= num_t * num_x )
        {
            data_offset = 0;
        }
        if ( data_offset < 0 )
        {
            data_offset = num_t * (num_x - num_batch);
        }
        delta_slice = 0;
    }
}

void display(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  drawStuff();
  glutSwapBuffers();
}

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
    /* Z near */ 1.0, /* Z far */ 1000.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, 4.75,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, +2.0);
  glRotatef(90, 1.0, 0.0, 0.0);
  //glRotatef(-20, 0.0, 0.0, 1.0);
}

struct pt
{
    int x;
    int y;
    pt ( int _x , int _y )
    {
        x = ( _x );
        y = ( _y );
    }
};

void color_region ( int k )
{
    for ( int i(0) ; i < num_t*num_batch ; ++i )
    {
        visited[i] = false;
    }
    std::deque < pt > Q;
    Q . push_back ( pt ( k/num_t , k%num_t ) );
    int my_type ( type [ k ] );
    visited         [ k ] = true    ;
    if ( layer )
    {
    velocity_data   [ k ] = vel_p   ;
    velocity_s_data [ k ] = vel_s   ;
    epsilon_data    [ k ] = epsilon ;
    delta_data      [ k ] = delta   ;
    }
    else
    {
    velocity_data_2   [ k ] = vel_p   ;
    velocity_s_data_2 [ k ] = vel_s   ;
    epsilon_data_2    [ k ] = epsilon ;
    delta_data_2      [ k ] = delta   ;
    }
    while ( ! Q . empty () )
    {
        int x = Q . front () . x;
        int y = Q . front () . y;
        Q . pop_front ();
        if ( x>0 )
        {
            int ind ( y + (x-1)*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( my_type - type [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x-1 , y ) );
                    visited         [ ind ] = true    ;
                    if ( layer )
                    {
                    velocity_data   [ ind ] = vel_p   ;
                    velocity_s_data [ ind ] = vel_s   ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                    }
                    else
                    {
                    velocity_data_2   [ ind ] = vel_p   ;
                    velocity_s_data_2 [ ind ] = vel_s   ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                    }
                }
            }
        }
        if ( y>0 )
        {
            int ind ( y-1 + x*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( my_type - type [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x , y-1 ) );
                    visited         [ ind ] = true    ;
                    if ( layer )
                    {
                    velocity_data   [ ind ] = vel_p   ;
                    velocity_s_data [ ind ] = vel_s   ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                    }
                    else
                    {
                    velocity_data_2   [ ind ] = vel_p   ;
                    velocity_s_data_2 [ ind ] = vel_s   ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                    }
                }
            }
        }
        if ( x+1<num_batch )
        {
            int ind ( y + (x+1)*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( my_type - type [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x+1 , y ) );
                    visited         [ ind ] = true    ;
                    if ( layer )
                    {
                    velocity_data   [ ind ] = vel_p   ;
                    velocity_s_data [ ind ] = vel_s   ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                    }
                    else
                    {
                    velocity_data_2   [ ind ] = vel_p   ;
                    velocity_s_data_2 [ ind ] = vel_s   ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                    }
                }
            }
        }
        if ( y+1<num_t )
        {
            int ind ( y+1 + x*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( my_type - type [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x , y+1 ) );
                    visited         [ ind ] = true    ;
                    if ( layer )
                    {
                    velocity_data   [ ind ] = vel_p   ;
                    velocity_s_data [ ind ] = vel_s   ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                    }
                    else
                    {
                    velocity_data_2   [ ind ] = vel_p   ;
                    velocity_s_data_2 [ ind ] = vel_s   ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                    }
                }
            }
        }
    }
}

void add_new_type ()
{
    for ( int k(0); k < num_batch*num_t ; ++k )
    {
        if ( !toggle[k] )
        {
            type[k] = type_index;
        }
    }
    type_index++;
}

// horizon file has format:
// x0 y0
// x1 y1
// x2 y2
// ...
void read_horizon ( std::string file_name , std::vector < point > & h )
{
    h . clear ();
    std::string line;
    std::ifstream myfile( file_name . c_str () );
    if ( myfile . is_open () )
    {
        while ( getline ( myfile , line ) )
        {
            std::stringstream ss(line);
            std::string coord;
            float x;
            float y;
            ss >> coord;
            x = atof(coord.c_str());
            ss >> coord;
            y = atof(coord.c_str());
            h.push_back(point(width*x/num_batch,height*y/num_t));
        }
        myfile . close ();
    }
    else
    {
        std::cout << "unable to open file " << file_name << std::endl;
        exit(1);
    }
}

void write_data ()
{

    // output data files

    std::vector < std::string > header_labels;
    std::vector < std::string > sort_order;
    {
        std::string out_file = velocity_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)velocity_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = velocity_s_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)velocity_s_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = epsilon_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)epsilon_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = delta_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)delta_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_x_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_x_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_y_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_y_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_z_in_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_z_data , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }



    {
        std::string out_file = velocity_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)velocity_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = velocity_s_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)velocity_s_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = epsilon_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)epsilon_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = delta_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)delta_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_x_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_x_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_y_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_y_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }
    {
        std::string out_file = nu_z_out_file;
        SEPWriter writer ( out_file . c_str () 
                         , 0 , 1 , num_t
                         , 0 , 1 , num_batch
                         , 0 , 1 , 1
                         , header_labels
                         , sort_order
                         , (out_file + std::string("@")) . c_str()
                         );
        writer . OpenDataFile ( (out_file + std::string("@")) . c_str() );
        writer . write_sepval ( (float*)nu_z_data_2 , writer . o1 , writer . o2 , writer . o3 , num_t * num_batch );
    }

    // run script

    std::cout << "begin system" << std::endl;

    std::cout << ( script_file . c_str () ) << std::endl;

    std::string cmd ( std::string ( "/bin/sh " ) + script_file );

    system ( cmd . c_str () );

    std::cout << "end system" << std::endl;

    read_horizon ( horizon_out_file , mig_selection );

}

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        case 'v': show_anisotropy_axis = !show_anisotropy_axis ; break ;
        case 'c': selection . clear () ; break ;
        case 'a': delta_slice-- ; break ;
        case 'd': delta_slice++ ; break ;
        case 'w': gain *= 1.1 ; break ;
        case 's': gain /= 1.1 ; break ;
        case 'p': print_selection () ; break ;
        case 'o': write_data () ; break ;
        case 't': velocity *= 1.01 ; break ;
        case 'g': velocity /= 1.01 ; break ;
        case 'r': shift *= 1.01 ; break ;
        case 'f': shift /= 1.01 ; break ;
        case 'y': fact *= 1.01 ; break ;
        case 'h': fact /= 1.01 ; break ;
        case 'l': layer = !layer ; break ;
        case '8': vel_p += 100 ; break ;
        case '2': vel_p -= 100 ; break ;
        case '4': vel_s -= 100 ; break ;
        case '6': vel_s += 100 ; break ;
        case '7': epsilon += 0.01 ; break ;
        case '1': epsilon -= 0.01 ; break ;
        case '9': delta += 0.01 ; break ;
        case '3': delta -= 0.01 ; break ;
        case '0': 
            if ( selection . size () == 0 ) 
            { 
                color_region ( ( (int)((y_pos/height)*num_t) + (int)((x_pos/width)*num_batch)*num_t ) ) ; 
            }
            else
            {
                std::cout << "selection in progress." << std::endl;
            }
            break ;
        case '*': 
            {
                int ind ( ( (int)((y_pos/height)*num_t) + (int)((x_pos/width)*num_batch)*num_t ) ) ; 
                if ( layer )
                {
                    vel_p   = velocity_data   [ ind ] ;
                    vel_s   = velocity_s_data [ ind ] ;
                    epsilon = epsilon_data    [ ind ] ;
                    delta   = delta_data      [ ind ] ;
                }
                else
                {
                    vel_p   = velocity_data_2   [ ind ] ;
                    vel_s   = velocity_s_data_2 [ ind ] ;
                    epsilon = epsilon_data_2    [ ind ] ;
                    delta   = delta_data_2      [ ind ] ;
                }
            }
            break;
        case '.': add_new_type () ; selection . clear () ; break ;
        case 'm': show_multiples = !show_multiples ; break ;
        case 'k':
            {
                if ( !layer )
                {
                    for ( int k(0) ; k < num_batch * num_t ; ++k )
                    {
                        velocity_data     [k] = velocity_data_2     [k];
                        velocity_s_data   [k] = velocity_s_data_2   [k];
                        epsilon_data      [k] = epsilon_data_2      [k];
                        delta_data        [k] = delta_data_2        [k];
                    }
                }
                else
                {
                    for ( int k(0) ; k < num_batch * num_t ; ++k )
                    {
                        velocity_data_2     [k] = velocity_data     [k];
                        velocity_s_data_2   [k] = velocity_s_data   [k];
                        epsilon_data_2      [k] = epsilon_data      [k];
                        delta_data_2        [k] = delta_data        [k];
                    }
                }
            }
            break;
        default : break ;
    }
    event[0].velocity = velocity;
    std::cout << "vel=" << velocity << " shift=" << shift << " fact=" << fact << std::endl;
}

bool add_new_pt = false;

int find_closest ( std::vector < point > const & sel
                 , int x_pos
                 , int y_pos
                 , float & min_dist
                 )
{
    if ( sel . size () == 0 ) { min_dist = 100000; return 0; }
    int closest_ind = -1;
    float dist = 100000;
    min_dist = 100000;
    float tmp_dist;
    float tmp_dist2;
    for ( std::size_t k(0)
        ; k < sel . size ()
        ; ++k
        )
    {
        tmp_dist = sqrtf ( (sel[k].x_pos-x_pos)*(sel[k].x_pos-x_pos) 
                         + (sel[k].y_pos-y_pos)*(sel[k].y_pos-y_pos) 
                         );
        tmp_dist2 = tmp_dist;
        if ( k > 0 && k + 1 < sel . size () )
        {
            float u_x ( sel[k+1].x_pos - sel[k].x_pos );
            float u_y ( sel[k+1].y_pos - sel[k].y_pos );
            float u_r ( sqrtf ( u_x * u_x + u_y * u_y ) );
            u_x /= u_r;
            u_y /= u_r;
            float proj ( u_x * ( x_pos - sel[k].x_pos ) + u_y * ( y_pos - sel[k].y_pos ) );
            if ( proj > 0 && proj < u_r )
            {
                tmp_dist = std :: min ( tmp_dist , sqrtf ( (x_pos-sel[k].x_pos-u_x*proj)*(x_pos-sel[k].x_pos-u_x*proj) 
                                                         + (y_pos-sel[k].y_pos-u_y*proj)*(y_pos-sel[k].y_pos-u_y*proj) 
                                                         )
                                      );
            }
        }
        if ( tmp_dist < dist )
        {
            dist = tmp_dist;
            closest_ind = k;
            if ( tmp_dist2 < min_dist )
            {
                min_dist = tmp_dist2;
            }
        }
    }
    return closest_ind;
}

void mouseClick ( int button
                , int state
                , int x
                , int y
                )
{
    if ( button == GLUT_LEFT_BUTTON )
    {
        int states = glutGetModifiers();
        if ( states == 0 )
        {
            if ( state == GLUT_DOWN )
            {
                float dist;
                selection_index = find_closest ( selection , x_pos , y_pos , dist );
                if ( dist > 20 )
                {
                    if ( selection_index > 0 && selection_index < (int)selection . size () )
                    {
                        selection_index += 1;
                    }
                    point event;
                    event . x_pos = x_pos;
                    event . y_pos = y_pos;
                    selection . insert ( selection . begin () + selection_index , event );
                }
                else
                {
                    move_pt = true;
                }
                add_new_pt = true;
            }
            else
            if ( state == GLUT_UP )
            {
                selection[selection_index].x_pos = x_pos;
                selection[selection_index].y_pos = y_pos;
                move_pt = false;
                calculate_toggle ();
                shift_pressed = false;
                ctrl_pressed = false;
            }
        }
        else
        if ( states == GLUT_ACTIVE_SHIFT )
        {
            if ( state == GLUT_DOWN )
            {
                if ( layer )
                {
                    for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
                    {
                        for ( float r(0) ; r < pen_radius ; ++r )
                        {
                            int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                            velocity_data   [ ind ] = vel_p   ;
                            velocity_s_data [ ind ] = vel_s   ;
                            epsilon_data    [ ind ] = epsilon ;
                            delta_data      [ ind ] = delta   ;
                        }
                    }
                }
                else
                {
                    for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
                    {
                        for ( float r(0) ; r < pen_radius ; ++r )
                        {
                            int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                            velocity_data_2   [ ind ] = vel_p   ;
                            velocity_s_data_2 [ ind ] = vel_s   ;
                            epsilon_data_2    [ ind ] = epsilon ;
                            delta_data_2      [ ind ] = delta   ;
                        }
                    }
                }
                shift_pressed = true;
            }
            else
            if ( state == GLUT_UP )
            {
                shift_pressed = false;
            }
        }
        else
        if ( states == GLUT_ACTIVE_CTRL )
        {
            if ( state == GLUT_DOWN )
            {
                if ( layer )
                {
                    for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
                    {
                        for ( float r(0) ; r < pen_radius ; ++r )
                        {
                            int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                            epsilon_data    [ ind ] = epsilon ;
                            delta_data      [ ind ] = delta   ;
                        }
                    }
                }
                else
                {
                    for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
                    {
                        for ( float r(0) ; r < pen_radius ; ++r )
                        {
                            int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                            epsilon_data_2    [ ind ] = epsilon ;
                            delta_data_2      [ ind ] = delta   ;
                        }
                    }
                }
                ctrl_pressed = true;
            }
            else
            if ( state == GLUT_UP )
            {
                ctrl_pressed = false;
            }
        }
    }
    if ( button == GLUT_RIGHT_BUTTON )
    {
        if ( state == GLUT_DOWN )
        {
            float dist;
            selection_index = find_closest ( selection , x_pos , y_pos , dist );
            if ( dist < 20 )
            {
                selection . erase ( selection . begin () + selection_index );
                calculate_toggle ();
            }
        }
    }
}

void mouseActiveMotion ( int x , int y )
{
    x_pos = x;
    y_pos = y;
    if ( !shift_pressed && !ctrl_pressed )
    {
        if ( add_new_pt )
        {
            selection[selection_index].x_pos = x_pos;
            selection[selection_index].y_pos = y_pos;
        }
    }
    else
    if ( shift_pressed )
    {
        if ( layer )
        {
            for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
            {
                for ( float r(0) ; r < pen_radius ; ++r )
                {
                    int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                    velocity_data   [ ind ] = vel_p   ;
                    velocity_s_data [ ind ] = vel_s   ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                }
            }
        }
        else
        {
            for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
            {
                for ( float r(0) ; r < pen_radius ; ++r )
                {
                    int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                    velocity_data_2   [ ind ] = vel_p   ;
                    velocity_s_data_2 [ ind ] = vel_s   ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                }
            }
        }
    }
    else
    if ( ctrl_pressed )
    {
        if ( layer )
        {
            for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
            {
                for ( float r(0) ; r < pen_radius ; ++r )
                {
                    int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                    epsilon_data    [ ind ] = epsilon ;
                    delta_data      [ ind ] = delta   ;
                }
            }
        }
        else
        {
            for ( float th(0) ; th < 2*M_PI ; th += 0.1 )
            {
                for ( float r(0) ; r < pen_radius ; ++r )
                {
                    int ind ( ( (int)(((y_pos+r*cos(th))/height)*num_t) + (int)(((x_pos+r*sin(th))/width)*num_batch)*num_t ) ) ;
                    epsilon_data_2    [ ind ] = epsilon ;
                    delta_data_2      [ ind ] = delta   ;
                }
            }
        }
    }
}

void mousePassiveMotion ( int x , int y )
{
    x_pos = x;
    y_pos = y;
    float dist;
    selection_index = find_closest ( selection , x_pos , y_pos , dist );
    {
        move_pt = dist < 20;
    }
    event[0].x_apex = x_pos;
    event[0].y_apex = y_pos;
    event[0].velocity = velocity;
}

void idle(void)
{
  glutPostRedisplay();
}

float find_abs_max ( float const * data , int n )
{
    float ret = 0;
    float temp;
    for ( int k(0)
        ; k < n
        ; ++k
        )
    {
        temp = fabs ( data[k] );
        if ( temp > ret )
        {
            ret = temp;
        }
    }
    return ret;
}

void flood_fill_by_param ( int k , int set_type )
{
    std::deque < pt > Q;
    Q . push_back ( pt ( k/num_t , k%num_t ) );
    float _vel_p     ( velocity_data     [k] );
    float _vel_s     ( velocity_s_data   [k] );
    float _epsilon   ( epsilon_data      [k] );
    float _delta     ( delta_data        [k] );
    while ( ! Q . empty () )
    {
        visited[k] = true;
        int x = Q . front () . x;
        int y = Q . front () . y;
        Q . pop_front ();
        if ( x>0 )
        {
            int ind ( y + (x-1)*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( _vel_p   - velocity_data     [ind] ) < 1e-5
                  && fabs ( _vel_s   - velocity_s_data   [ind] ) < 1e-5
                  && fabs ( _epsilon - epsilon_data      [ind] ) < 1e-5
                  && fabs ( _delta   - delta_data        [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x-1 , y ) );
                    visited[ind] = true;
                    type[ind] = set_type;
                }
            }
        }
        if ( y>0 )
        {
            int ind ( y-1 + x*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( _vel_p   - velocity_data     [ind] ) < 1e-5
                  && fabs ( _vel_s   - velocity_s_data   [ind] ) < 1e-5
                  && fabs ( _epsilon - epsilon_data      [ind] ) < 1e-5
                  && fabs ( _delta   - delta_data        [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x , y-1 ) );
                    visited[ind] = true;
                    type[ind] = set_type;
                }
            }
        }
        if ( x+1<num_batch )
        {
            int ind ( y + (x+1)*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( _vel_p   - velocity_data     [ind] ) < 1e-5
                  && fabs ( _vel_s   - velocity_s_data   [ind] ) < 1e-5
                  && fabs ( _epsilon - epsilon_data      [ind] ) < 1e-5
                  && fabs ( _delta   - delta_data        [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x+1 , y ) );
                    visited[ind] = true;
                    type[ind] = set_type;
                }
            }
        }
        if ( y+1<num_t )
        {
            int ind ( y+1 + x*num_t );
            if ( !visited[ind] )
            {
                if ( fabs ( _vel_p   - velocity_data     [ind] ) < 1e-5
                  && fabs ( _vel_s   - velocity_s_data   [ind] ) < 1e-5
                  && fabs ( _epsilon - epsilon_data      [ind] ) < 1e-5
                  && fabs ( _delta   - delta_data        [ind] ) < 1e-5
                   )
                {
                    Q . push_back ( pt ( x , y+1 ) );
                    visited[ind] = true;
                    type[ind] = set_type;
                }
            }
        }
    }
}

void init_types ()
{
    std::cout << "init types begin" << std::endl;
    while ( 1 )
    {
        bool all_visited = true;
        int src_index = 0;
        for ( int k(0) ; k < num_batch * num_t ; ++k )
        {
            if ( !visited[k] )
            {
                all_visited = false;
                src_index = k;
                break;
            }
        }
        if ( all_visited )
        {
            break;
        }
        flood_fill_by_param ( src_index , type_index++ );
    }
    std::cout << "init types end" << std::endl;
}

void init_data ()
{
    for ( int k(0) ; k < num_t * num_batch ; ++k )
    {
        velocity_s_data[k] = 0.5*velocity_data[k];
        velocity_s_data_2[k] = 0.5*velocity_data_2[k];
        epsilon_data[k] = 0;
        epsilon_data_2[k] = 0;
        delta_data[k] = 0;
        delta_data_2[k] = 0;
        nu_x_data[k] = 0;
        nu_x_data_2[k] = 0;
        nu_y_data[k] = 0;
        nu_y_data_2[k] = 0;
        nu_z_data[k] = 1;
        nu_z_data_2[k] = 1;
    }
}

void normalize_nu ()
{
    for ( int k(0) ; k < num_t * num_batch ; ++k )
    {
        if ( fabs(nu_x_data[k]) > 0.5 * fabs(nu_z_data[k]) ) nu_x_data[k] = 0.5 * nu_x_data[k] * fabs(nu_z_data[k]) / fabs(nu_x_data[k]);
        if ( fabs(nu_x_data_2[k]) > 0.5 * fabs(nu_z_data_2[k]) ) nu_x_data_2[k] = 0.5 * nu_x_data_2[k] * fabs(nu_z_data_2[k]) / fabs(nu_x_data_2[k]);
        float r ( sqrtf ( nu_x_data[k] * nu_x_data[k] 
                        + nu_y_data[k] * nu_y_data[k] 
                        + nu_z_data[k] * nu_z_data[k] 
                        ) 
                );
        if ( nu_z_data[k] > 0 ) r *= -1;
        float r2( sqrtf ( nu_x_data_2[k] * nu_x_data_2[k] 
                        + nu_y_data_2[k] * nu_y_data_2[k] 
                        + nu_z_data_2[k] * nu_z_data_2[k] 
                        ) 
                );
        if ( nu_z_data_2[k] > 0 ) r2 *= -1;
        nu_x_data[k] /= r;
        nu_y_data[k] /= r;
        nu_z_data[k] /= r;
        nu_x_data_2[k] /= r2;
        nu_y_data_2[k] /= r2;
        nu_z_data_2[k] /= r2;
    }
}

int main( int argc , char ** argv )
{
    if ( argc == 18 )
    {
        SEPReader   velocity_reader ( argv[1] );
        num_t = velocity_reader . n1;
        num_batch = velocity_reader . n2;
        num_x = velocity_reader . n2 
              * velocity_reader . n3
              * velocity_reader . n4
              * velocity_reader . n5
              * velocity_reader . n6
              * velocity_reader . n7
              * velocity_reader . n8
              ;
        velocity_in_file   = std::string ( argv[1 ] );
        velocity_s_in_file = std::string ( argv[2 ] );
        epsilon_in_file    = std::string ( argv[3 ] );
        delta_in_file      = std::string ( argv[4 ] );
        nu_x_in_file       = std::string ( argv[5 ] );
        nu_y_in_file       = std::string ( argv[6 ] );
        nu_z_in_file       = std::string ( argv[7 ] );
        SEPReader velocity_s_reader ( argv[2] );
        SEPReader    epsilon_reader ( argv[3] );
        SEPReader      delta_reader ( argv[4] );
        SEPReader       nu_x_reader ( argv[5] );
        SEPReader       nu_y_reader ( argv[6] );
        SEPReader       nu_z_reader ( argv[7] );
        velocity_out_file   = std::string ( argv[8 ] );
        velocity_s_out_file = std::string ( argv[9 ] );
        epsilon_out_file    = std::string ( argv[10] );
        delta_out_file      = std::string ( argv[11] );
        nu_x_out_file       = std::string ( argv[12] );
        nu_y_out_file       = std::string ( argv[13] );
        nu_z_out_file       = std::string ( argv[14] );
        SEPReader   velocity_2_reader ( argv[8 ] );
        SEPReader velocity_s_2_reader ( argv[9 ] );
        SEPReader    epsilon_2_reader ( argv[10] );
        SEPReader      delta_2_reader ( argv[11] );
        SEPReader       nu_x_2_reader ( argv[12] );
        SEPReader       nu_y_2_reader ( argv[13] );
        SEPReader       nu_z_2_reader ( argv[14] );

        horizon_in_file  = std::string ( argv[15] );
        read_horizon ( horizon_in_file , input_horizon );
        horizon_out_file = std::string ( argv[16] );
        read_horizon ( horizon_out_file , output_horizon );
        
        script_file      = std::string ( argv[17] );

        if ( num_t != velocity_s_reader . n1 && num_batch != velocity_s_reader . n2 )
        {
            std::cout << "velocity_s has inconsistent dimensions" << std::endl;
            return 0;
        }
        if ( num_t !=    epsilon_reader . n1 && num_batch !=    epsilon_reader . n2 )
        {
            std::cout << "epsilon has inconsistent dimensions" << std::endl;
            return 0;
        }
        if ( num_t !=      delta_reader . n1 && num_batch !=      delta_reader . n2 )
        {
            std::cout << "delta has inconsistent dimensions" << std::endl;
            return 0;
        }
        if ( num_t !=       nu_x_reader . n1 && num_batch !=       nu_x_reader . n2 )
        {
            std::cout << "nu_x has inconsistent dimensions" << std::endl;
            return 0;
        }
        if ( num_t !=       nu_y_reader . n1 && num_batch !=       nu_y_reader . n2 )
        {
            std::cout << "nu_y has inconsistent dimensions" << std::endl;
            return 0;
        }
        if ( num_t !=       nu_z_reader . n1 && num_batch !=       nu_z_reader . n2 )
        {
            std::cout << "nu_z has inconsistent dimensions" << std::endl;
            return 0;
        }
        std::cout << num_x << std::endl;

        velocity_data = new float [ num_x * num_t ];
        memset ( &velocity_data[0] , 0 , num_x * num_t );
        velocity_reader . read_sepval ( & velocity_data [ 0 ] , velocity_reader . o1 , velocity_reader . o2 , velocity_reader . o3 , num_x * num_t );

        velocity_s_data = new float [ num_x * num_t ];
        memset ( &velocity_s_data[0] , 0 , num_x * num_t );
        velocity_s_reader . read_sepval ( & velocity_s_data [ 0 ] , velocity_s_reader . o1 , velocity_s_reader . o2 , velocity_s_reader . o3 , num_x * num_t );

        epsilon_data = new float [ num_x * num_t ];
        memset ( &epsilon_data[0] , 0 , num_x * num_t );
        epsilon_reader . read_sepval ( & epsilon_data [ 0 ] , epsilon_reader . o1 , epsilon_reader . o2 , epsilon_reader . o3 , num_x * num_t );

        delta_data = new float [ num_x * num_t ];
        memset ( &delta_data[0] , 0 , num_x * num_t );
        delta_reader . read_sepval ( & delta_data [ 0 ] , delta_reader . o1 , delta_reader . o2 , delta_reader . o3 , num_x * num_t );

        nu_x_data = new float [ num_x * num_t ];
        memset ( &nu_x_data[0] , 0 , num_x * num_t );
        nu_x_reader . read_sepval ( & nu_x_data [ 0 ] , nu_x_reader . o1 , nu_x_reader . o2 , nu_x_reader . o3 , num_x * num_t );

        nu_y_data = new float [ num_x * num_t ];
        memset ( &nu_y_data[0] , 0 , num_x * num_t );
        nu_y_reader . read_sepval ( & nu_y_data [ 0 ] , nu_y_reader . o1 , nu_y_reader . o2 , nu_y_reader . o3 , num_x * num_t );

        nu_z_data = new float [ num_x * num_t ];
        memset ( &nu_z_data[0] , 0 , num_x * num_t );
        nu_z_reader . read_sepval ( & nu_z_data [ 0 ] , nu_z_reader . o1 , nu_z_reader . o2 , nu_z_reader . o3 , num_x * num_t );

        velocity_data_2 = new float [ num_x * num_t ];
        memset ( &velocity_data_2[0] , 0 , num_x * num_t );
        velocity_2_reader . read_sepval ( & velocity_data_2 [ 0 ] , velocity_reader . o1 , velocity_reader . o2 , velocity_reader . o3 , num_x * num_t );

        velocity_s_data_2 = new float [ num_x * num_t ];
        memset ( &velocity_s_data_2[0] , 0 , num_x * num_t );
        velocity_s_2_reader . read_sepval ( & velocity_s_data_2 [ 0 ] , velocity_s_reader . o1 , velocity_s_reader . o2 , velocity_s_reader . o3 , num_x * num_t );

        epsilon_data_2 = new float [ num_x * num_t ];
        memset ( &epsilon_data_2[0] , 0 , num_x * num_t );
        epsilon_2_reader . read_sepval ( & epsilon_data_2 [ 0 ] , epsilon_reader . o1 , epsilon_reader . o2 , epsilon_reader . o3 , num_x * num_t );

        delta_data_2 = new float [ num_x * num_t ];
        memset ( &delta_data_2[0] , 0 , num_x * num_t );
        delta_2_reader . read_sepval ( & delta_data_2 [ 0 ] , delta_reader . o1 , delta_reader . o2 , delta_reader . o3 , num_x * num_t );

        nu_x_data_2 = new float [ num_x * num_t ];
        memset ( &nu_x_data_2[0] , 0 , num_x * num_t );
        nu_x_2_reader . read_sepval ( & nu_x_data_2 [ 0 ] , nu_x_reader . o1 , nu_x_reader . o2 , nu_x_reader . o3 , num_x * num_t );

        nu_y_data_2 = new float [ num_x * num_t ];
        memset ( &nu_y_data_2[0] , 0 , num_x * num_t );
        nu_y_2_reader . read_sepval ( & nu_y_data_2 [ 0 ] , nu_y_reader . o1 , nu_y_reader . o2 , nu_y_reader . o3 , num_x * num_t );

        nu_z_data_2 = new float [ num_x * num_t ];
        memset ( &nu_z_data_2[0] , 0 , num_x * num_t );
        nu_z_2_reader . read_sepval ( & nu_z_data_2 [ 0 ] , nu_z_reader . o1 , nu_z_reader . o2 , nu_z_reader . o3 , num_x * num_t );

        toggle = new bool [ num_x * num_t ];
        for ( int k(0) ; k < num_x * num_t ; ++k ) toggle[k] = true;

        type = new int [ num_x * num_t ];
        for ( int k(0) ; k < num_x * num_t ; ++k ) type[k] = 0;

        visited = new bool [ num_x * num_t ];
        for ( int k(0) ; k < num_x * num_t ; ++k ) visited[k] = false;

        //init_types ();

        scalar = 0.45f / find_abs_max ( & velocity_data [ 0 ] , num_x * num_t );
        std::cout << scalar << std::endl;

        // init_data ();
        
        normalize_nu ();

    }
    else
    {
        std::cout << "try ./a.out sep_file.h" << std::endl;
        return 0;
    }

    std::cout << "display" << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width,height);
    glutCreateWindow(argv[1]);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseClick);
    glutMotionFunc(mouseActiveMotion);
    glutPassiveMotionFunc(mousePassiveMotion);
    glutIdleFunc(idle);
    init();
    glutMainLoop();

    return 0;
}


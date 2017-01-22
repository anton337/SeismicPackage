#include <GL/glut.h>
#include <iostream>
#include <iomanip>
#include <string.h>
#include <vector>
#include <sstream>
#include <math.h>
#include <map>
#include "sep_reader.h"

bool toggle_contours = false;

bool move_mouse = true;

int num_batch;

int width  = 1000;
int height = 1000;

int num_x;
int num_t;

int delta = 1;

int timer = 0;

int data_offset = 0;

float scalar = 0;

float gain = 1;

float * data = NULL;

bool * toggle = NULL;

std::string quote;

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

std::vector < hyperbolic_event > event(1);

struct point
{
    float x_pos;
    float y_pos;
};

std::vector < point > selection;

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
    if ( show_multiples )
    {
    glColor3f ( 1 , 0 , 0 );
    glBegin(GL_LINES);
    glVertex3f( 2 * x_pos / width - 1 - 0.050 , 0 , 2 * y_pos / height - 1         );
    glVertex3f( 2 * x_pos / width - 1 + 0.050 , 0 , 2 * y_pos / height - 1         );
    glVertex3f( 2 * x_pos / width - 1         , 0 , 2 * y_pos / height - 1 - 0.050 );
    glVertex3f( 2 * x_pos / width - 1         , 0 , 2 * y_pos / height - 1 + 0.050 );
    glEnd();
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
                {
                    y_val = shift + sqrt ( (event[k].y_apex - shift)*(event[k].y_apex - shift)*mult*mult + fact*fact * (x_val-event[k].x_apex)*(x_val-event[k].x_apex) / (event[k].velocity*event[k].velocity) );
                }
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
    if ( quote.size() > 0 )
    {
        glColor3f ( 1 , 1 , 1 );
        std::stringstream str;
        str << quote << std::endl;
        drawString(GLUT_BITMAP_HELVETICA_18, str.str().c_str(), -0.275, 0.1, -0.275);
    }
    float val;
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
        val = 0.5f + data[k+data_offset] * scalar * gain;
        if ( !toggle_contours )
        {
        if ( toggle[k] )
        glColor3f
        ( val
        , val
        , val
        );
        else
        glColor3f
        ( val
        , 0  
        , val
        );
        }
        else
        {
        if ( (int)( (val-0.5f) * 2000 ) % 20 == 0 )
        {
        glColor3f
        ( val
        , 0 
        , 0
        );
        }
        else
        {
        glColor3f
        ( val
        , val
        , val
        );
        }
        }
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
        data_offset += num_t * num_batch * delta;
        if ( data_offset >= num_t * num_x )
        {
            data_offset = 0;
        }
        if ( data_offset < 0 )
        {
            data_offset = num_t * (num_x - num_batch);
        }
        delta = 0;
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

void keyboard ( unsigned char key , int x , int y )
{
    switch ( key )
    {
        case 27 : exit(0) ; break ;
        case 'l': toggle_contours = !toggle_contours ; break ;
        case 'k': move_mouse = !move_mouse ; break ;
        case 'a': delta-- ; break ;
        case 'd': delta++ ; break ;
        case 'w': gain *= 1.1 ; break ;
        case 's': gain /= 1.1 ; break ;
        case 'p': print_selection () ; break ;
        case 't': velocity *= 1.01 ; break ;
        case 'g': velocity /= 1.01 ; break ;
        case 'r': shift *= 1.01 ; break ;
        case 'f': shift /= 1.01 ; break ;
        case 'y': fact *= 1.01 ; break ;
        case 'h': fact /= 1.01 ; break ;
        case 'm': show_multiples = !show_multiples ; break ;
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
    if ( add_new_pt )
    {
        selection[selection_index].x_pos = x_pos;
        selection[selection_index].y_pos = y_pos;
    }
}

void mousePassiveMotion ( int x , int y )
{
    if ( move_mouse )
    {
    x_pos = x;
    y_pos = y;
    }
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

int main( int argc , char ** argv )
{
    if ( argc == 2 )
    {
        SEPReader reader ( argv[1] );
        num_t = reader . n1;
        num_batch = reader . n2;
        num_x = reader . n2 
              * reader . n3
              * reader . n4
              * reader . n5
              * reader . n6
              * reader . n7
              * reader . n8
              ;
        std::cout << num_x << std::endl;
        data = new float [ num_x * num_t ];
        memset ( &data[0] , 0 , num_x * num_t );
        toggle = new bool [ num_x * num_t ];
        for ( int k(0) ; k < num_x * num_t ; ++k ) toggle[k] = true;
        reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );
        scalar = 0.45f / find_abs_max ( & data [ 0 ] , num_x * num_t );
        std::cout << scalar << std::endl;
    }
    else
    if ( argc > 2 )
    {
        SEPReader reader ( argv[1] , false );
        num_t = reader . n1;
        num_batch = reader . n2;
        num_x = reader . n2 
              * reader . n3
              * reader . n4
              * reader . n5
              * reader . n6
              * reader . n7
              * reader . n8
              ;
        std::cout << num_x << std::endl;
        data = new float [ num_x * num_t ];
        memset ( &data[0] , 0 , num_x * num_t );
        toggle = new bool [ num_x * num_t ];
        for ( int k(0) ; k < num_x * num_t ; ++k ) toggle[k] = true;
        for ( int k(2)
            ; k < argc
            ; ++k
            )
        {
            std::cout << "k=" << k-1 << "     ";
            reader . OpenDataFile ( argv[k] );
            reader . read_sepval ( & data [ (k-1) * num_t * num_batch ] , reader . o1 , reader . o2 , reader . o3 , num_batch * num_t );
        }
        scalar = 0.45f / find_abs_max ( & data [ 0 ] , num_x * num_t );
        std::cout << scalar << std::endl;
    }
    else
    {
        std::cout << "try ./a.out sep_file.h" << std::endl;
        return 0;
    }

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


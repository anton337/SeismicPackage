#include <GL/glut.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"

enum SeismicDataOrder { COMMON_SHOT = 0
                      , COMMON_CHAN = 1
                      };

SeismicDataOrder seismicDataOrder = COMMON_CHAN;

bool migrate = true;

bool show_travel_time_contours = false;

bool lock = false;

int sou_init_offset = 112;
int rec_init_offset = 8;
int t_offset = 25;

int height = 1000;
int width  = 1900;

int num_t;
int num_batch;
int num_batches;
int num_x;

int tt_num_t;
int tt_num_batch;
int tt_num_batches;
int tt_num_x;

int rec=0;
int sou=0;


int x_pos , y_pos;

float T_FACT = 1000;

float image_scale = 0.0001;

template < typename T >
inline void mem_cpy ( T * out , T * in , int num )
{
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        out[k] = in[k];
    }
}

float get_argument ( int &   arg 
                   , int     argc
                   , char ** argv
                   )
{
    ++arg;
    if ( arg > argc )
    {
        std::cout << "parameter error" << std::endl;
        exit(1);
    }
    return atof ( argv[arg] );
}

float * data;

float * chan_data;

float * image_data;

float * vel_data;

float * tt_data;

float scalar = 1;

void drawStuff(void)
{
    glPointSize(3);
    glColor3f(1,1,1);
    glBegin(GL_LINES);
    for ( int k(1)
        ; k < 10
        ; ++k
        )
    {
        glVertex3f(-1+1.0f*(1000.0*k)/9200.0,0,-1);
        glVertex3f(-1+1.0f*(1000.0*k)/9200.0,0, 1);
    }
    glEnd();
    glColor3f(1,0,1);
    glBegin(GL_LINES);
    glVertex3f(-1+2.0f*(float)x_pos/width-0.050,0,-1+2.0f*(float)y_pos/height);
    glVertex3f(-1+2.0f*(float)x_pos/width+0.050,0,-1+2.0f*(float)y_pos/height);
    glVertex3f(-1+2.0f*(float)x_pos/width,0,-1+2.0f*(float)y_pos/height-0.100);
    glVertex3f(-1+2.0f*(float)x_pos/width,0,-1+2.0f*(float)y_pos/height+0.100);
    glEnd();
    long rec_offset = (long)rec * (long)tt_num_t * (long)tt_num_batch;
    long sou_offset = (long)sou * (long)tt_num_t * (long)tt_num_batch;
    long sou_seis_offset = (long)(sou-sou_init_offset) * (long)num_t * (long)num_batch; 
    float val  = 0;
    float valr = 0;
    float valg = 0;
    float valb = 0;
    int X_pos = (2.0f*(float)x_pos/width)*tt_num_batch;
    int T_pos = ((float)y_pos/height)*tt_num_t;
    glBegin(GL_POINTS);
    for ( long k(0)
        , x(0)
        ; x < tt_num_batch
        ; ++x
        )
    {
        for ( int t(0)
            ; t < tt_num_t
            ; ++t
            , ++k
            )
        {
            if ( show_travel_time_contours )
            {
                valr = fabs(tt_data[k+rec_offset]) * scalar;
                valr = ( (int)(valr*500) % 25 == 0 ) ? 1 : valr ;
                valb = fabs(tt_data[k+sou_offset]) * scalar;
                valb = ( (int)(valb*500) % 25 == 0 ) ? 1 : valb ;
            }
            else
            {
                valr = 0;
                valb = 0;
            }
            valg = vel_data[k] * 0.0003;
            valr += image_scale*image_data[k];
            valg += image_scale*image_data[k];
            valb += image_scale*image_data[k];
            if ( pow ( (float)( X_pos - x ) / tt_num_batch , 2 ) + pow ( (float)( T_pos - t ) / tt_num_t , 2 ) < 0.00001 ) 
            {
                valr *= 0.5f;
                valg *= 0.5f;
                valb *= 0.5f;
            }
            glColor3f
            ( valr
            , valg
            , valb
            );
            glVertex3f
            ( -1 + 1.0f * x/tt_num_batch
            , 0
            , -1 + 2.0f * t/tt_num_t
            );
        }
    }
    glEnd();
    glPointSize(14);
    switch ( seismicDataOrder )
    {
        case COMMON_SHOT :
        {
            if ( sou_seis_offset >= 0 && sou_seis_offset < (long)num_t*(long)num_batch*(long)num_batches )
            {
                glBegin(GL_POINTS);
                for ( long k(0)
                    , x(0)
                    ; x < num_batch
                    ; ++x
                    )
                {
                    for ( int t(0)
                        ; t < num_t
                        ; ++t
                        , ++k
                        )
                    {
                        val = 0.5f + 0.001f*data[k+sou_seis_offset];
                        valr = val;
                        valg = val;
                        for ( int m(-40)
                            ; m <= 40
                            ; ++m
                            )
                        {
                            if ( num_x + rec_init_offset == (sou-rec-m) + 1 + x )
                            {
                                if ( T_pos >= 0 && T_pos < tt_num_t && X_pos >= 0 && X_pos < tt_num_batch )
                                {
                                    int T = T_FACT * (fabs(tt_data[rec_offset+m*tt_num_t*tt_num_batch+X_pos*tt_num_t+T_pos]) + fabs(tt_data[sou_offset+X_pos*tt_num_t+T_pos])) + t_offset;
                                    if ( fabs(T - t) < 3 ) 
                                    {
                                        valg = 1;
                                    }
                                }
                            }
                        }
                        glColor3f
                        ( valr
                        , valg
                        , val
                        );
                        glVertex3f
                        ( 1.0f * x/num_batch
                        , 0
                        , -1 + 2.0f * t/num_t
                        );
                    }
                }
                glEnd();
            }
        }
        break;
        case COMMON_CHAN :
        {
            int rec_ind = sou-rec-rec_init_offset;
            if ( rec_ind >= 0 && rec_ind < num_x )
            {
                int rec_off = (num_x-1-rec_ind) * num_batches*num_t;
                glBegin(GL_POINTS);
                for ( long x(0)
                    , k(0)
                    ; x < num_batches
                    ; ++x
                    )
                {
                    for ( int t(0)
                        ; t < num_t
                        ; ++t
                        , ++k
                        )
                    {
                        val = 0.5f + 0.001f*chan_data[rec_off + k];
                        valr = val;
                        valg = val;
                        for ( int m(0)
                            ; m < num_batches
                            ; ++m
                            )
                        {
                            if ( m+sou_init_offset >= 0 && m+sou_init_offset < tt_num_batches )
                            if ( m+sou_init_offset-rec_ind >= 0 && m+sou_init_offset-rec_ind < tt_num_batches )
                            if ( x == m+4 )
                            {
                                if ( T_pos >= 0 && T_pos < tt_num_t && X_pos >= 0 && X_pos < tt_num_batch )
                                {
                                    int T = T_FACT * (fabs(tt_data[(m+sou_init_offset)*tt_num_t*tt_num_batch+X_pos*tt_num_t+T_pos]) + fabs(tt_data[(m+sou_init_offset-rec_ind)*tt_num_t*tt_num_batch+X_pos*tt_num_t+T_pos])) + t_offset;
                                    if ( fabs(T - t) < 3 ) 
                                    {
                                        valg = 1;
                                    }
                                }
                            }
                        }
                        glColor3f
                        ( valr
                        , valg
                        , val
                        );
                        glVertex3f
                        ( 1.0f * x/num_batches
                        , 0
                        , -1 + 2.0f * t/num_t
                        );
                    }
                }
                glEnd();
            }
        }
        break;
        default:
        break;
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
        case 'w': rec++; if ( rec >= tt_num_batches ) rec = 0;                break ;
        case 's': rec--; if ( rec < 0               ) rec = tt_num_batches-1; break ;
        case 'q': rec++; if ( rec >= tt_num_batches ) rec = 0;                
                  sou++; if ( sou >= tt_num_batches ) sou = 0;                break ;
        case 'a': rec--; if ( rec < 0               ) rec = tt_num_batches-1; 
                  sou--; if ( sou < 0               ) sou = tt_num_batches-1; break ;
        case 'z': scalar *= 1.1; break;
        case 'x': scalar /= 1.1; break;
        case 'e': T_FACT *= 1.01; break;
        case 'd': T_FACT /= 1.01; break;
        case 'r': image_scale *= 1.1; break;
        case 'f': image_scale /= 1.1; break;
        case 'm': seismicDataOrder = (SeismicDataOrder)(!seismicDataOrder); break;
        case 't': sou_init_offset++; break;
        case 'g': sou_init_offset--; break;
        case 'y': rec_init_offset++; break;
        case 'h': rec_init_offset--; break;
        case 'u': t_offset++; break;
        case 'j': t_offset--; break;
        case 'p': show_travel_time_contours = !show_travel_time_contours ; break ;
        case 'l': lock = !lock ; break ;
        default : break ;
    }
    std::cout << "image_scale=" << image_scale << std::endl;
    std::cout << "scalar=" << scalar << std::endl;
    std::cout << "T_FACT=" << T_FACT << std::endl;
    std::cout << "sou_init_offset=" << sou_init_offset << std::endl;
    std::cout << "rec_init_offset=" << rec_init_offset << std::endl;
    std::cout << "t_offset=" << t_offset << std::endl;
}

void idle(void)
{
  glutPostRedisplay();
}

void mouseActiveMotion ( int x , int y )
{
    if ( !lock )
    {
        x_pos = x;
        y_pos = y;
    }
}

void mousePassiveMotion ( int x , int y )
{
    if ( !lock )
    {
        x_pos = x;
        y_pos = y;
    }
}


int main( int argc , char ** argv )
{

    if (argc == 1)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string travel_time_file_name( argv[arg] );

    arg++;

    std::string travel_time_data_file_name( argv[arg] );

    arg++;

    std::string velocity_file_name( argv[arg] );

    arg++;

    std::string velocity_data_file_name( argv[arg] );

    arg++;

    std::string seismic_file_name( argv[arg] );

    arg++;

    std::string seismic_data_file_name( argv[arg] );

    arg++;

    std::string common_chan_seismic_file_name( argv[arg] );

    arg++;

    std::string common_chan_seismic_data_file_name( argv[arg] );

    arg++;

    std::string image_file_name( argv[arg] );

    arg++;

    std::string image_data_file_name( argv[arg] );
std::cout << "p1" << std::endl;

    SEPReader reader ( seismic_file_name . c_str () , false );
    num_t = reader . n1;
    num_batch = reader . n2;
    num_x = num_batch;
    num_batches = reader . n3
                * reader . n4
                * reader . n5
                * reader . n6
                * reader . n7
                * reader . n8
                ;
    data = new float [ num_x * num_t * num_batches ];
    memset ( &data[0] , 0 , num_x * num_t * num_batches );

std::cout << "p2" << std::endl;
    reader . OpenDataFile ( (seismic_data_file_name) . c_str() );

    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t * num_batches );

    SEPReader common_chan_reader ( common_chan_seismic_file_name . c_str () , false );
    chan_data = new float [ num_x * num_t * num_batches ];
    memset ( &chan_data[0] , 0 , num_x * num_t * num_batches );
    common_chan_reader . OpenDataFile ( (common_chan_seismic_data_file_name) . c_str() );
    common_chan_reader . read_sepval ( & chan_data [ 0 ] , common_chan_reader . o1 , common_chan_reader . o2 , common_chan_reader . o3 , num_x * num_t * num_batches );

    SEPReader reader_time ( travel_time_file_name . c_str () , false );
    tt_num_t = reader_time . n1;
    tt_num_batch = reader_time . n2;
    tt_num_x = tt_num_batch;
    tt_num_batches = reader_time . n3
                * reader_time . n4
                * reader_time . n5
                * reader_time . n6
                * reader_time . n7
                * reader_time . n8
                ;
    tt_data = new float [ tt_num_x * tt_num_t * tt_num_batches ];
    memset ( &tt_data[0] , 0 , tt_num_x * tt_num_t * tt_num_batches );

std::cout << "p3" << std::endl;
    reader_time . OpenDataFile ( (travel_time_data_file_name) . c_str() );

    reader_time . read_sepval ( & tt_data [ 0 ] , reader_time . o1 , reader_time . o2 , reader_time . o3 , tt_num_x * tt_num_t * tt_num_batches );

    SEPReader vel_reader ( velocity_file_name . c_str () , false );
    vel_data = new float [ tt_num_x * tt_num_t ];
    memset ( &vel_data[0] , 0 , tt_num_x * tt_num_t );
    vel_reader . OpenDataFile ( (velocity_data_file_name) . c_str() );
    vel_reader . read_sepval ( & vel_data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , tt_num_x * tt_num_t );

std::cout << "p4" << std::endl;




    image_data = new float [ tt_num_x * tt_num_t ];
    memset ( &image_data[0] , 0 , tt_num_x * tt_num_t );

    int ind_sou , ind_rec ;
    if ( migrate )
    for ( int ind2 = 0 // sou
        ; ind2 < tt_num_batches
        ; ind2 += 1
        )
    {
        std::cout << "ind2=" << ind2 << "     num_batches=" << tt_num_batches << std::endl;
        for ( int ind1 = 0 // rec
            ; ind1 < tt_num_batches
            ; ind1 += 1
            )
        {
            for ( int x(0)
                , k(0)
                ; x < tt_num_x
                ; ++x
                )
            {
                for ( int z(0)
                    ; z < tt_num_t
                    ; ++z
                    , ++k
                    )
                {
                    {
                        ind_sou = ind2-sou_init_offset;//-115;
                        if ( ind_sou >= 0 && ind_sou < num_batches )
                        {
                            ind_rec = ind2-ind1-rec_init_offset;
                            if ( ind_rec >= 0 && ind_rec < num_x )
                            {
                                //if ( fabs ( 0.5f * ( ind_sou + ind_rec ) - x/7 ) < 140 )
                                {
                                    int _ind1 = ind1 * tt_num_x * tt_num_t;
                                    int _ind2 = ind2 * tt_num_x * tt_num_t;
                                    int T = T_FACT * (fabs(tt_data[_ind1+k]) + fabs(tt_data[_ind2+k]))+t_offset;
                                    if ( T >= 20 && T < num_t - 60 )
                                    {
                                        image_data[k] += data[ind_sou*num_batch*num_t + (num_batch-1-ind_rec)*num_t + T];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //{

    //    int num_f = tt_num_t / 2 + 1;
    //    fftwf_complex * fdata = new fftwf_complex [ tt_num_x * num_f ];

    //    float * tmp = new float [ tt_num_t ];
    //    fftwf_complex * ftmp = new fftwf_complex [ num_f ];

    //    fftwf_plan plan ( fftwf_plan_dft_r2c_1d ( tt_num_t
    //                                            , & tmp[0]
    //                                            , &ftmp[0]
    //                                            , FFTW_ESTIMATE
    //                                            )
    //                    );

    //    for ( int x(0)
    //        ; x < tt_num_x
    //        ; ++x
    //        )
    //    {
    //        mem_cpy ( & tmp[0] , & image_data[x*tt_num_t] , tt_num_t );
    //        fftwf_execute ( plan );
    //        mem_cpy < float > ( (float*)(&fdata[x*num_f]) , (float*)(&ftmp[0]) , 2*num_f );
    //    }

    //    fftwf_destroy_plan ( plan );

    //    for ( int x(0)
    //        ; x < tt_num_x
    //        ; ++x
    //        )
    //    {
    //        for ( int f(0)
    //            ; f < num_f * 0.15
    //            ; ++f
    //            )
    //        {
    //            fdata[x*num_f+f][0] = 0;
    //            fdata[x*num_f+f][1] = 0;
    //        }
    //    }

    //    fftwf_plan plan_i ( fftwf_plan_dft_c2r_1d ( tt_num_t
    //                                              , &ftmp[0]
    //                                              , & tmp[0]
    //                                              , FFTW_ESTIMATE
    //                                              )
    //                      );

    //    for ( int x(0)
    //        ; x < tt_num_x
    //        ; ++x
    //        )
    //    {
    //        mem_cpy < float > ( (float*)(&ftmp[0]) , (float*)(&fdata[x*num_f]) , 2*num_f );
    //        fftwf_execute ( plan_i );
    //        mem_cpy ( &image_data[x*tt_num_t] , & tmp[0] , tt_num_t );
    //    }

    //    fftwf_destroy_plan ( plan_i );

    //}

std::cout << "p5" << std::endl;


    //SEPWriter writer ( image_file_name . c_str () 
    //                 , reader . o1 , reader . d1 , tt_num_t
    //                 , reader . o2 , reader . d2 , tt_num_x
    //                 , reader . o3 , reader . d3 , 1
    //                 , reader . get_header_labels ()
    //                 , reader . get_sort_order ()
    //                 , (image_data_file_name) . c_str()
    //                 );

    //writer . OpenDataFile ( (image_data_file_name) . c_str() );

    //writer . write_sepval ( (float*)image_data , reader . o1 , reader . o2 , reader . o3 , tt_num_x * (tt_num_t) );

std::cout << "p6" << std::endl;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width,height);
    glutCreateWindow(argv[1]);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutIdleFunc(idle);
    glutMotionFunc(mouseActiveMotion);
    glutPassiveMotionFunc(mousePassiveMotion);
    init();
    glutMainLoop();

//    SEPWriter writer ( image_file_name . c_str () 
//                     , reader . o1 , reader . d1 , num_t
//                     , reader . o2 , reader . d2 , num_x
//                     , reader . o3 , reader . d3 , num_batches
//                     , reader . get_header_labels ()
//                     , reader . get_sort_order ()
//                     , (image_data_file_name) . c_str()
//                     );
//
//    writer . OpenDataFile ( (image_data_file_name) . c_str() );
//
//    writer . write_sepval ( (float*)image_data , reader . o1 , reader . o2 , reader . o3 , num_x * (num_t) * num_batches );


    return 0;

}


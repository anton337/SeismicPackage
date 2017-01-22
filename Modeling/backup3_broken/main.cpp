#include <iostream>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include <map>
#include "sep_writer.h"
#include "sep_reader.h"
#include "matrix.h"
#include "ray_tracing.h"

float theta = 0;

float anis_d = 0.10;
float anis_e = 0.10;

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

void construct_rotation_matrix ( float ax
                               , float ay
                               , float az
                               , float angle
                               , matrix < float > * a
                               , float * data // needs to be size 9
                               )
{
    float r = raytracing::norm ( ax , ay , az );
    if ( r > 1e-10 && fabs ( angle ) > 1e-10 )
    {
        ax /= r;
        ay /= r;
        az /= r;
        float c = cosf ( angle );
        float mc= 1-c;
        float s =-sinf ( angle );
        data[0] = c + ax*ax * mc;
        data[1] = ax * ay * mc - az * s;
        data[2] = ax * az * mc + ay * s;
        data[3] = ay * ax * mc + az * s;
        data[4] = c + ay*ay * mc;
        data[5] = ay * az * mc - ax * s;
        data[6] = az * ax * mc - ay * s;
        data[7] = az * ay * mc + ax * s;
        data[8] = c + az*az * mc;
        *a = matrix < float > ( 3 , 3 , data );
    }
    else
    {
        data[0] = 1;
        data[1] = 0;
        data[2] = 0;
        data[3] = 0;
        data[4] = 1;
        data[5] = 0;
        data[6] = 0;
        data[7] = 0;
        data[8] = 1;
        *a = matrix < float > ( 3 , 3 , data );
    }
}

void construct_bond_matrix ( matrix < float > const * a
                           , float * data   // needs to be size 36
                           , float * data_t // needs to be size 36
                           , matrix < float > * m
                           , matrix < float > * m_t
                           )
{

// M = [ a11^2      a12^2       a13^2       2 a12 a13            2 a13 a11            2 a11 a12         ]
//     [ a21^2      a22^2       a23^2       2 a22 a23            2 a23 a21            2 a21 a22         ]
//     [ a31^2      a32^2       a33^2       2 a32 a33            2 a33 a31            2 a31 a32         ]
//     [ a21 a31    a22 a32     a23 a33     a22 a33 + a23 a32    a21 a33 + a23 a31    a22 a31 + a21 a32 ]
//     [ a31 a11    a32 a12     a33 a13     a12 a33 + a13 a32    a13 a31 + a11 a33    a11 a32 + a12 a31 ]
//     [ a11 a21    a12 a22     a13 a23     a12 a23 + a13 a22    a13 a21 + a11 a23    a11 a22 + a12 a21 ]

    data[0]  = a -> array [0] * a -> array [0];
    data[1]  = a -> array [1] * a -> array [1];
    data[2]  = a -> array [2] * a -> array [2];

    data[3]  = 2 * a -> array [1] * a -> array [2];
    data[4]  = 2 * a -> array [2] * a -> array [0];
    data[5]  = 2 * a -> array [0] * a -> array [1];

    data[6]  = a -> array [3] * a -> array [3];
    data[7]  = a -> array [4] * a -> array [4];
    data[8]  = a -> array [5] * a -> array [5];

    data[9]  = 2 * a -> array [4] * a -> array [5];
    data[10] = 2 * a -> array [5] * a -> array [3];
    data[11] = 2 * a -> array [3] * a -> array [4];

    data[12] = a -> array [6] * a -> array [6];
    data[13] = a -> array [7] * a -> array [7];
    data[14] = a -> array [8] * a -> array [8];

    data[15] = 2 * a -> array [7] * a -> array [8];
    data[16] = 2 * a -> array [8] * a -> array [6];
    data[17] = 2 * a -> array [6] * a -> array [7];

    data[18] = a -> array [3] * a -> array [6];
    data[19] = a -> array [4] * a -> array [7];
    data[20] = a -> array [5] * a -> array [8];

    data[21] = a -> array [4] * a -> array [8] + a -> array [5] * a -> array [7];
    data[22] = a -> array [3] * a -> array [8] + a -> array [5] * a -> array [6];
    data[23] = a -> array [4] * a -> array [6] + a -> array [3] * a -> array [7];

    data[24] = a -> array [6] * a -> array [0];
    data[25] = a -> array [7] * a -> array [1];
    data[26] = a -> array [8] * a -> array [2];

    data[27] = a -> array [1] * a -> array [8] + a -> array [2] * a -> array [7];
    data[28] = a -> array [2] * a -> array [6] + a -> array [0] * a -> array [8];
    data[29] = a -> array [0] * a -> array [7] + a -> array [1] * a -> array [6];

    data[30] = a -> array [0] * a -> array [3];
    data[31] = a -> array [1] * a -> array [4];
    data[32] = a -> array [2] * a -> array [5];

    data[33] = a -> array [1] * a -> array [5] + a -> array [2] * a -> array [4];
    data[34] = a -> array [2] * a -> array [3] + a -> array [0] * a -> array [5];
    data[35] = a -> array [0] * a -> array [4] + a -> array [1] * a -> array [3];

    *m = matrix < float > ( 6 , 6 , data );

    for ( int x(0) ,k(0) ; x < 6 ; ++x )
    for ( int y(0) ; y < 6 ; ++y , ++k )
    data_t[k] = data[x+y*6];

    *m_t = matrix < float > ( 6 , 6 , data_t );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// C_TTI = M C_VTI M^T
//
// a - rotation matrix
//
//     [ a11 a12 a13 ]
// a = [ a21 a22 a23 ]
//     [ a31 a32 a33 ]
//
// rotation matrix around arbitrary axis
//
//     [ c + ax^2 (1-c)         ax ay (1-c) - az s          ax az (1-c) + ay s ]
// a = [ ay ax (1-c) + az s     c + ay^2 (1-c)              ay az (1-c) - ax s ]
//     [ az ax (1-c) - ay s     az ay (1-c) + ax s          c + az^2 (1-c)     ]
//
// M - bond matrix
//
// M = [ a11^2      a12^2       a13^2       2 a12 a13            2 a13 a11            2 a11 a12         ]
//     [ a21^2      a22^2       a23^2       2 a22 a23            2 a23 a21            2 a21 a22         ]
//     [ a31^2      a32^2       a33^2       2 a32 a33            2 a33 a31            2 a31 a32         ]
//     [ a21 a31    a22 a32     a23 a33     a22 a33 + a23 a32    a21 a33 + a23 a31    a22 a31 + a21 a32 ]
//     [ a31 a11    a32 a12     a33 a13     a12 a33 + a13 a32    a13 a31 + a11 a33    a11 a32 + a12 a31 ]
//     [ a11 a21    a12 a22     a13 a23     a12 a23 + a13 a22    a13 a21 + a11 a23    a11 a22 + a12 a21 ]
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void convert_Thomson_parameters_to_moduli_tti ( int           size
                                              , float       * A11
                                              , float       * A12
                                              , float       * A13
                                              , float       * A21
                                              , float       * A22
                                              , float       * A23
                                              , float       * A31
                                              , float       * A32
                                              , float       * A33
                                              , float       * Q11
                                              , float       * Q12
                                              , float       * Q13
                                              , float       * Q21
                                              , float       * Q22
                                              , float       * Q23
                                              , float       * Q31
                                              , float       * Q32
                                              , float       * Q33
                                              , float       * C11
                                              , float       * C12
                                              , float       * C13
                                              , float       * C21
                                              , float       * C22
                                              , float       * C23
                                              , float       * C31
                                              , float       * C32
                                              , float       * C33
                                              , float       * G11
                                              , float       * G12
                                              , float       * G13
                                              , float       * G21
                                              , float       * G22
                                              , float       * G23
                                              , float       * G31
                                              , float       * G32
                                              , float       * G33
                                              , float       * nux
                                              , float       * nuy
                                              , float       * nuz
                                              , float const * rho      
                                              , float const * vp
                                              , float const * delta    = NULL
                                              , float const * epsilon  = NULL
                                              , float const * vs       = NULL
                                              , float const * gamma    = NULL
                                              )
{
    float c33;
    float c44;
    float c11;
    float c66;
    float c13;
    float array [36];
    float stiffness [36];
    float stiffness_t [36];
    float stiffness_x [36];
    float ax = 0;
    float ay = 1;
    float az = 0;
    float angle;
    matrix < float > * rotation = new matrix < float > ();
    float data   [36];
    float data_t [36];
    matrix < float > * bond   = new matrix < float > ();
    matrix < float > * bond_t = new matrix < float > ();
    float rotation_data[9];
    for ( int k(0)
        ; k < size
        ; ++k
        )
    {
        c33 = vp[k] * vp[k] * rho[k];
        c44 = (vs)?vs[k] * vs[k] * rho[k]:c33;
        c11 = ( 1 + 2 * ((epsilon)?epsilon[k]:0) ) * c33;
        c66 = ( 1 + 2 * ((gamma  )?gamma  [k]:0) ) * c44;
        // c13 = sqrtf( 2 * ((delta)?delta[k]:0) * c33 * ( c33 - c44 ) + ( c33 - c44 ) * ( c33 - c44 ) ) - c44;
        c13 = rho[k] * sqrtf( ( vp[k]*vp[k] - vs[k]*vs[k] ) * ( vp[k]*vp[k]*(1+2*((delta)?delta[k]:0)) - vs[k]*vs[k] ) ) - c44;
        array[0] = c11; array[1] = (c11 - 2*c66); array[2] = c13; array[3] = 0; array[4] = 0; array[5] = 0;
        array[6] = (c11 - 2*c66); array[7] = c11; array[8] = c13; array[9] = 0; array[10] = 0; array[11] = 0;
        array[12] = c13; array[13] = c13; array[14] = c33; array[15] = 0; array[16] = 0; array[17] = 0;
        array[18] = 0; array[19] = 0; array[20] = 0; array[21] = c44; array[22] = 0; array[23] = 0;
        array[24] = 0; array[25] = 0; array[26] = 0; array[27] = 0; array[28] = c44; array[29] = 0;
        array[30] = 0; array[31] = 0; array[32] = 0; array[33] = 0; array[34] = 0; array[35] = c66;
        if ( k == 0 )
        for ( int x(0) , k(0); x<6; ++x )
        {
            for ( int y(0); y<6; ++y , ++k )
            {
                std::cout << array[k] << " \t";
            }
            std::cout << std::endl;
        }
        angle = raytracing :: angle ( nux[k] , nuy[k] , nuz[k] , 0.0f , 0.0f , 1.0f );
        construct_rotation_matrix ( ax
                                  , ay
                                  , az
                                  , angle
                                  , rotation
                                  , rotation_data
                                  );
        construct_bond_matrix ( rotation
                              , data   // needs to be size 36
                              , data_t // needs to be size 36
                              , bond
                              , bond_t
                              );
        bond   -> operator () ( &array[0 ] , &stiffness[0 ] );
        bond   -> operator () ( &array[6 ] , &stiffness[6 ] );
        bond   -> operator () ( &array[12] , &stiffness[12] );
        bond   -> operator () ( &array[18] , &stiffness[18] );
        bond   -> operator () ( &array[24] , &stiffness[24] );
        bond   -> operator () ( &array[30] , &stiffness[30] );
        for ( int x(0) ,k(0) ; x < 6 ; ++x )
        for ( int y(0) ; y < 6 ; ++y , ++k )
        stiffness_x[k] = stiffness[x+y*6];
        *bond_t = matrix < float > ( 6 , 6 , stiffness_x );
        bond_t     -> operator () ( &data[0 ] , &stiffness_t[0 ] );
        bond_t     -> operator () ( &data[6 ] , &stiffness_t[6 ] );
        bond_t     -> operator () ( &data[12] , &stiffness_t[12] );
        bond_t     -> operator () ( &data[18] , &stiffness_t[18] );
        bond_t     -> operator () ( &data[24] , &stiffness_t[24] );
        bond_t     -> operator () ( &data[30] , &stiffness_t[30] );
        // 11 -> 0
        // 12 -> 1
        // 13 -> 2
        // 14 -> 3
        // 15 -> 4
        // 16 -> 5
        // 21 -> 6
        // 22 -> 7
        // 23 -> 8
        // 24 -> 9
        // 25 -> 10
        // 26 -> 11
        // 31 -> 12
        // 32 -> 13
        // 33 -> 14
        // 34 -> 15
        // 35 -> 16
        // 36 -> 17
        // 41 -> 18
        // 42 -> 19
        // 43 -> 20
        // 44 -> 21
        // 45 -> 22
        // 46 -> 23
        // 51 -> 24
        // 52 -> 25
        // 53 -> 26
        // 54 -> 27
        // 55 -> 28
        // 56 -> 29
        // 61 -> 30
        // 62 -> 31
        // 63 -> 32
        // 64 -> 33
        // 65 -> 34
        // 66 -> 35
        // A = [ c11 c16 c15 ]
        //     [ c16 c66 c56 ]
        //     [ c15 c56 c55 ]
        //     c55 = c44 for VTI
        if(A11)A11[k] = stiffness_t[0] ; if(A21)A21[k] = stiffness_t[5] ; if(A31)A31[k] = stiffness_t[4] ;
        if(A12)A12[k] = stiffness_t[30]; if(A22)A22[k] = stiffness_t[35]; if(A32)A32[k] = stiffness_t[29];
        if(A13)A13[k] = stiffness_t[24]; if(A23)A23[k] = stiffness_t[34]; if(A33)A33[k] = stiffness_t[28];
        // Q = [ c55 c45 c35 ]
        //     [ c45 c44 c34 ]
        //     [ c35 c34 c33 ]
        //     c55 = c44 for VTI
        if(Q11)Q11[k] = stiffness_t[28]; if(Q21)Q21[k] = stiffness_t[22]; if(Q31)Q31[k] = stiffness_t[16];
        if(Q12)Q12[k] = stiffness_t[27]; if(Q22)Q22[k] = stiffness_t[21]; if(Q32)Q32[k] = stiffness_t[15];
        if(Q13)Q13[k] = stiffness_t[26]; if(Q23)Q23[k] = stiffness_t[34]; if(Q33)Q33[k] = stiffness_t[14];
        // C = [ c15 c14 c13 ]
        //     [ c56 c46 c36 ]
        //     [ c55 c45 c35 ]
        if(C11)C11[k] = stiffness_t[4] ; if(C21)C21[k] = stiffness_t[3] ; if(C31)C31[k] = stiffness_t[2] ;
        if(C12)C12[k] = stiffness_t[29]; if(C22)C22[k] = stiffness_t[23]; if(C32)C32[k] = stiffness_t[17];
        if(C13)C13[k] = stiffness_t[28]; if(C23)C23[k] = stiffness_t[22]; if(C33)C33[k] = stiffness_t[16];
        // G = [ c15 c56 c55 ]
        //     [ c14 c46 c45 ]
        //     [ c13 c36 c35 ]
        if(G11)G11[k] = stiffness_t[4] ; if(G21)G21[k] = stiffness_t[29]; if(G31)G31[k] = stiffness_t[28];
        if(G12)G12[k] = stiffness_t[3] ; if(G22)G22[k] = stiffness_t[23]; if(G32)G32[k] = stiffness_t[22];
        if(G13)G13[k] = stiffness_t[2] ; if(G23)G23[k] = stiffness_t[17]; if(G33)G33[k] = stiffness_t[16];
        if ( k == 0 )
        {
            {
                std::cout << G11[k] << " \t" << G12[k] << " \t" << G13[k] << " \t" << std::endl;
                std::cout << G21[k] << " \t" << G22[k] << " \t" << G23[k] << " \t" << std::endl;
                std::cout << G31[k] << " \t" << G32[k] << " \t" << G33[k] << " \t" << std::endl;
            }
        }
    }
}

void convert_Thomson_parameters_to_moduli ( int           size
                                          , float       * A11
                                          , float       * A12
                                          , float       * A13
                                          , float       * A21
                                          , float       * A22
                                          , float       * A23
                                          , float       * A31
                                          , float       * A32
                                          , float       * A33
                                          , float       * Q11
                                          , float       * Q12
                                          , float       * Q13
                                          , float       * Q21
                                          , float       * Q22
                                          , float       * Q23
                                          , float       * Q31
                                          , float       * Q32
                                          , float       * Q33
                                          , float       * C11
                                          , float       * C12
                                          , float       * C13
                                          , float       * C21
                                          , float       * C22
                                          , float       * C23
                                          , float       * C31
                                          , float       * C32
                                          , float       * C33
                                          , float       * G11
                                          , float       * G12
                                          , float       * G13
                                          , float       * G21
                                          , float       * G22
                                          , float       * G23
                                          , float       * G31
                                          , float       * G32
                                          , float       * G33
                                          , float const * rho      
                                          , float const * vp
                                          , float const * delta    = NULL
                                          , float const * epsilon  = NULL
                                          , float const * vs       = NULL
                                          , float const * gamma    = NULL
                                          )
{
    float c33;
    float c44;
    float c11;
    float c66;
    float c13;
    for ( int k(0)
        ; k < size
        ; ++k
        )
    {
        c33 = vp[k] * vp[k] * rho[k];
        c44 = (vs)?vs[k] * vs[k] * rho[k]:c33;
        c11 = ( 1 + 2 * ((epsilon)?epsilon[k]:0) ) * c33;
        c66 = ( 1 + 2 * ((delta  )?delta  [k]:0) ) * c44;
        c13 = sqrtf( 2 * ((delta)?delta[k]:0) * c33 * ( c33 - c44 ) + ( c33 - c44 ) * ( c33 - c44 ) ) - c44;
        // A = [ c11 c16 c15 ]
        //     [ c16 c66 c56 ]
        //     [ c15 c56 c55 ]
        //     c55 = c44 for VTI
        if(A11)A11[k] = c11 ; if(A12)A12[k] = 0   ; if(A13)A13[k] = 0   ;
        if(A21)A21[k] = 0   ; if(A22)A22[k] = c66 ; if(A23)A23[k] = 0   ;
        if(A31)A31[k] = 0   ; if(A32)A32[k] = 0   ; if(A33)A33[k] = c44 ;
        // Q = [ c55 c45 c35 ]
        //     [ c45 c44 c34 ]
        //     [ c35 c34 c33 ]
        //     c55 = c44 for VTI
        if(Q11)Q11[k] = c44 ; if(Q12)Q12[k] = 0   ; if(Q13)Q13[k] = 0   ;
        if(Q21)Q21[k] = 0   ; if(Q22)Q22[k] = c44 ; if(Q23)Q23[k] = 0   ;
        if(Q31)Q31[k] = 0   ; if(Q32)Q32[k] = 0   ; if(Q33)Q33[k] = c33 ;
        // C = [ c15 c14 c13 ]
        //     [ c56 c46 c36 ]
        //     [ c55 c45 c35 ]
        if(C11)C11[k] = 0   ; if(C12)C12[k] = 0      ; if(C13)C13[k] = c13 ;
        if(C21)C21[k] = 0   ; if(C22)C22[k] = 0      ; if(C23)C23[k] = 0   ;
        if(C31)C31[k] = c44 ; if(C32)C32[k] = 0      ; if(C33)C33[k] = 0   ;
        // G = [ c15 c56 c55 ]
        //     [ c14 c46 c45 ]
        //     [ c13 c36 c35 ]
        if(G11)G11[k] = 0   ; if(G12)G12[k] = 0      ; if(G13)G13[k] = c44 ;
        if(G21)G21[k] = 0   ; if(G22)G22[k] = 0      ; if(G23)G23[k] = 0   ;
        if(G31)G31[k] = c13 ; if(G32)G32[k] = 0      ; if(G33)G33[k] = 0   ;
    }
}

void propagate_anisotropic ( float         h 
                           , int           num_x 
                           , int           num_t 
                           , float       * p2x
                           , float       * p2y
                           , float       * p2z
                           , float       * p1x
                           , float       * p1y
                           , float       * p1z
                           , float       * cx
                           , float       * cy
                           , float       * cz
                           , float       * vel 
                           , float       * rho
                           , float const * A11
                           , float const * A12
                           , float const * A13
                           , float const * A21
                           , float const * A22
                           , float const * A23
                           , float const * A31
                           , float const * A32
                           , float const * A33
                           , float const * Q11
                           , float const * Q12
                           , float const * Q13
                           , float const * Q21
                           , float const * Q22
                           , float const * Q23
                           , float const * Q31
                           , float const * Q32
                           , float const * Q33
                           , float const * C11
                           , float const * C12
                           , float const * C13
                           , float const * C21
                           , float const * C22
                           , float const * C23
                           , float const * C31
                           , float const * C32
                           , float const * C33
                           , float const * G11
                           , float const * G12
                           , float const * G13
                           , float const * G21
                           , float const * G22
                           , float const * G23
                           , float const * G31
                           , float const * G32
                           , float const * G33
                           , bool absorbing_boundary_condition = true
                           )
{
    float fact ( h * h );
    if ( cx )
    {
        for ( int t(1)
            ; t+1 < num_t
            ; ++t
            )
        {
            for ( int x(1)
                ; x+1 < num_x
                ; ++x
                )
            {
                {
                    //
                    // A = [ (c11)  c16   c15  ]
                    //     [  c16  (c66)  c56  ]
                    //     [  c15   c56  (c55) ]
                    //     c55 = c44 for VTI
                    //
                    // Q = [ (c55)  c45   c35  ]
                    //     [  c45  (c44)  c34  ]
                    //     [  c35   c34  (c33) ]
                    //     c55 = c44 for VTI
                    //
                    // C = [  c15  c14 (c13) ]
                    //     [  c56  c46  c36  ]
                    //     [ (c55) c45  c35  ]
                    //
                    // G = [  c15  c56 (c55) ]
                    //     [  c14  c46  c45  ]
                    //     [ (c13) c36  c35  ]
                    //
                    cx[t+x*num_t] = 0.5 * ( fact / rho[t+x*num_t] ) * ( ( ( (A11&&p1x)?( ( A11[t+(x+1)*num_t] + A11[t+x*num_t] ) * ( p1x[t    +(x+1)*num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A12&&p1y)?( ( A12[t+(x+1)*num_t] + A12[t+x*num_t] ) * ( p1y[t    +(x+1)*num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A13&&p1z)?( ( A13[t+(x+1)*num_t] + A13[t+x*num_t] ) * ( p1z[t    +(x+1)*num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        ) 
                                                                      - ( ( (A11&&p1x)?( ( A11[t+(x-1)*num_t] + A11[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A12&&p1y)?( ( A12[t+(x-1)*num_t] + A12[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A13&&p1z)?( ( A13[t+(x-1)*num_t] + A13[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[t    +(x-1)*num_t] ) ): 0 )
                                                                        ) 
                                                                      + ( ( (Q11&&p1x)?( ( Q11[t+(x+1)*num_t] + Q11[t+x*num_t] ) * ( p1x[(t+1)+x    *num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q12&&p1y)?( ( Q12[t+(x+1)*num_t] + Q12[t+x*num_t] ) * ( p1y[(t+1)+x    *num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q13&&p1z)?( ( Q13[t+(x+1)*num_t] + Q13[t+x*num_t] ) * ( p1z[(t+1)+x    *num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        ) 
                                                                      - ( ( (Q11&&p1x)?( ( Q11[t+(x-1)*num_t] + Q11[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q12&&p1y)?( ( Q12[t+(x-1)*num_t] + Q12[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q13&&p1z)?( ( Q13[t+(x-1)*num_t] + Q13[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[(t-1)+x    *num_t] ) ): 0 )
                                                                        )
                                                                      + 0.5 * ( ( ( (C11&&p1x)? ( C11[t  +(x+1)*num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t-1+(x+1)*num_t] ) - C11[t  +(x-1)*num_t] * ( p1x[t+1+(x-1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (C12&&p1y)? ( C12[t  +(x+1)*num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t-1+(x+1)*num_t] ) - C12[t  +(x-1)*num_t] * ( p1y[t+1+(x-1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (C13&&p1z)? ( C13[t  +(x+1)*num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t-1+(x+1)*num_t] ) - C13[t  +(x-1)*num_t] * ( p1z[t+1+(x-1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 )
                                                                                ) 
                                                                              + ( ( (G11&&p1x)? ( G11[t+1+x    *num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t+1+(x-1)*num_t] ) - G11[t-1+x    *num_t] * ( p1x[t-1+(x+1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (G12&&p1y)? ( G12[t+1+x    *num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t+1+(x-1)*num_t] ) - G12[t-1+x    *num_t] * ( p1y[t-1+(x+1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (G13&&p1z)? ( G13[t+1+x    *num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t+1+(x-1)*num_t] ) - G13[t-1+x    *num_t] * ( p1z[t-1+(x+1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 )
                                                                                )
                                                                              )
                                                                      )
                                  + 2 * p1x[t+x*num_t]
                                  -     p2x[t+x*num_t]
                                  ;
                }
            }
        }
        if ( absorbing_boundary_condition )
        {
            // absorbing boundary t = 0
            for ( int x(0)
                , t(0)
                ; x < num_x
                ; ++x
                )
            {
                cx[t+x*num_t] = ( 2 * p1x[t+x*num_t] 
                                    - p2x[t+x*num_t] 
                                    - cx [t+2+x*num_t] 
                                + 2 * p1x[t+2+x*num_t] 
                                    - p2x[t+2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cx[t+2+x*num_t] 
                                                          - p2x[t+2+x*num_t] 
                                                          + p2x[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1x[t+2+x*num_t] 
                                                                                - 2 * p1x[t+1+x*num_t] 
                                                                                    + p1x[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary t = num_t-1
            for ( int x(0)
                , t(num_t-1)
                ; x < num_x
                ; ++x
                )
            {
                cx[t+x*num_t] = ( 2 * p1x[t+x*num_t] 
                                    - p2x[t+x*num_t] 
                                    - cx [t-2+x*num_t] 
                                + 2 * p1x[t-2+x*num_t] 
                                    - p2x[t-2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cx[t-2+x*num_t] 
                                                          - p2x[t-2+x*num_t] 
                                                          + p2x[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1x[t-2+x*num_t] 
                                                                                - 2 * p1x[t-1+x*num_t] 
                                                                                    + p1x[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = 0
            for ( int x(0)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cx[t+x*num_t] = ( 2 * p1x[t+x*num_t] 
                                    - p2x[t+x*num_t] 
                                    - cx [t+(x+2)*num_t] 
                                + 2 * p1x[t+(x+2)*num_t] 
                                    - p2x[t+(x+2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cx[t+(x+2)*num_t] 
                                                          - p2x[t+(x+2)*num_t] 
                                                          + p2x[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1x[t+(x+2)*num_t] 
                                                                                - 2 * p1x[t+(x+1)*num_t] 
                                                                                    + p1x[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = num_x-1
            for ( int x(num_x-1)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cx[t+x*num_t] = ( 2 * p1x[t+x*num_t] 
                                    - p2x[t+x*num_t] 
                                    - cx [t+(x-2)*num_t] 
                                + 2 * p1x[t+(x-2)*num_t] 
                                    - p2x[t+(x-2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cx[t+(x-2)*num_t] 
                                                          - p2x[t+(x-2)*num_t] 
                                                          + p2x[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1x[t+(x-2)*num_t] 
                                                                                - 2 * p1x[t+(x-1)*num_t] 
                                                                                    + p1x[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p2x[k] = p1x[k];
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p1x[k] = cx[k];
        }
    }
    if ( cy )
    {
        for ( int t(1)
            ; t+1 < num_t
            ; ++t
            )
        {
            for ( int x(1)
                ; x+1 < num_x
                ; ++x
                )
            {
                {
                    cy[t+x*num_t] = 0.5 * ( fact / rho[t+x*num_t] ) * ( ( ( (A21&&p1x)?( ( A21[t+(x+1)*num_t] + A21[t+x*num_t] ) * ( p1x[t    +(x+1)*num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A22&&p1y)?( ( A22[t+(x+1)*num_t] + A22[t+x*num_t] ) * ( p1y[t    +(x+1)*num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A23&&p1z)?( ( A23[t+(x+1)*num_t] + A23[t+x*num_t] ) * ( p1z[t    +(x+1)*num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        )
                                                                      - ( ( (A21&&p1x)?( ( A21[t+(x-1)*num_t] + A21[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A22&&p1y)?( ( A22[t+(x-1)*num_t] + A22[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A23&&p1z)?( ( A23[t+(x-1)*num_t] + A23[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[t    +(x-1)*num_t] ) ): 0 )
                                                                        )
                                                                      + ( ( (Q21&&p1x)?( ( Q21[t+(x+1)*num_t] + Q21[t+x*num_t] ) * ( p1x[(t+1)+x    *num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q22&&p1y)?( ( Q22[t+(x+1)*num_t] + Q22[t+x*num_t] ) * ( p1y[(t+1)+x    *num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q23&&p1z)?( ( Q23[t+(x+1)*num_t] + Q23[t+x*num_t] ) * ( p1z[(t+1)+x    *num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        )
                                                                      - ( ( (Q21&&p1x)?( ( Q21[t+(x-1)*num_t] + Q21[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q22&&p1y)?( ( Q22[t+(x-1)*num_t] + Q22[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q23&&p1z)?( ( Q23[t+(x-1)*num_t] + Q23[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[(t-1)+x    *num_t] ) ): 0 )
                                                                        )
                                                                      + 0.5 * ( ( ( (C21&&p1x)? ( C21[t  +(x+1)*num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t-1+(x+1)*num_t] ) - C21[t  +(x-1)*num_t] * ( p1x[t+1+(x-1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                + ( (C22&&p1y)? ( C22[t  +(x+1)*num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t-1+(x+1)*num_t] ) - C22[t  +(x-1)*num_t] * ( p1y[t+1+(x-1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                + ( (C23&&p1z)? ( C23[t  +(x+1)*num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t-1+(x+1)*num_t] ) - C23[t  +(x-1)*num_t] * ( p1z[t+1+(x-1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                )
                                                                              + ( ( (G21&&p1x)? ( G21[t+1+x    *num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t+1+(x-1)*num_t] ) - G21[t-1+x    *num_t] * ( p1x[t-1+(x+1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                + ( (G22&&p1y)? ( G22[t+1+x    *num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t+1+(x-1)*num_t] ) - G22[t-1+x    *num_t] * ( p1y[t-1+(x+1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                + ( (G23&&p1z)? ( G23[t+1+x    *num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t+1+(x-1)*num_t] ) - G23[t-1+x    *num_t] * ( p1z[t-1+(x+1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 ) 
                                                                                )
                                                                              )
                                                                      )
                                  + 2 * p1y[t+x*num_t]
                                  -     p2y[t+x*num_t]
                                  ;
                }
            }
        }
        if ( absorbing_boundary_condition )
        {
            // absorbing boundary t = 0
            for ( int x(0)
                , t(0)
                ; x < num_x
                ; ++x
                )
            {
                cy[t+x*num_t] = ( 2 * p1y[t+x*num_t] 
                                    - p2y[t+x*num_t] 
                                    - cy [t+2+x*num_t] 
                                + 2 * p1y[t+2+x*num_t] 
                                    - p2y[t+2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cy[t+2+x*num_t] 
                                                          - p2y[t+2+x*num_t] 
                                                          + p2y[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1y[t+2+x*num_t] 
                                                                                - 2 * p1y[t+1+x*num_t] 
                                                                                    + p1y[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary t = num_t-1
            for ( int x(0)
                , t(num_t-1)
                ; x < num_x
                ; ++x
                )
            {
                cy[t+x*num_t] = ( 2 * p1y[t+x*num_t] 
                                    - p2y[t+x*num_t] 
                                    - cy [t-2+x*num_t] 
                                + 2 * p1y[t-2+x*num_t] 
                                    - p2y[t-2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cy[t-2+x*num_t] 
                                                          - p2y[t-2+x*num_t] 
                                                          + p2y[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1y[t-2+x*num_t] 
                                                                                - 2 * p1y[t-1+x*num_t] 
                                                                                    + p1y[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = 0
            for ( int x(0)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cy[t+x*num_t] = ( 2 * p1y[t+x*num_t] 
                                    - p2y[t+x*num_t] 
                                    - cy [t+(x+2)*num_t] 
                                + 2 * p1y[t+(x+2)*num_t] 
                                    - p2y[t+(x+2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cy[t+(x+2)*num_t] 
                                                          - p2y[t+(x+2)*num_t] 
                                                          + p2y[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1y[t+(x+2)*num_t] 
                                                                                - 2 * p1y[t+(x+1)*num_t] 
                                                                                    + p1y[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = num_x-1
            for ( int x(num_x-1)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cy[t+x*num_t] = ( 2 * p1y[t+x*num_t] 
                                    - p2y[t+x*num_t] 
                                    - cy [t+(x-2)*num_t] 
                                + 2 * p1y[t+(x-2)*num_t] 
                                    - p2y[t+(x-2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cy[t+(x-2)*num_t] 
                                                          - p2y[t+(x-2)*num_t] 
                                                          + p2y[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1y[t+(x-2)*num_t] 
                                                                                - 2 * p1y[t+(x-1)*num_t] 
                                                                                    + p1y[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p2y[k] = p1y[k];
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p1y[k] = cy[k];
        }
    }
    if ( cz )
    {
        for ( int t(1)
            ; t+1 < num_t
            ; ++t
            )
        {
            for ( int x(1)
                ; x+1 < num_x
                ; ++x
                )
            {
                {
                    cz[t+x*num_t] = 0.5 * ( fact / rho[t+x*num_t] ) * ( ( ( (A31&&p1x)?( ( A31[t+(x+1)*num_t] + A31[t+x*num_t] ) * ( p1x[t    +(x+1)*num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A32&&p1y)?( ( A32[t+(x+1)*num_t] + A32[t+x*num_t] ) * ( p1y[t    +(x+1)*num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (A33&&p1z)?( ( A33[t+(x+1)*num_t] + A33[t+x*num_t] ) * ( p1z[t    +(x+1)*num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        ) 
                                                                      - ( ( (A31&&p1x)?( ( A31[t+(x-1)*num_t] + A31[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A32&&p1y)?( ( A32[t+(x-1)*num_t] + A32[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[t    +(x-1)*num_t] ) ): 0 )
                                                                        + ( (A33&&p1z)?( ( A33[t+(x-1)*num_t] + A33[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[t    +(x-1)*num_t] ) ): 0 )
                                                                        ) 
                                                                      + ( ( (Q31&&p1x)?( ( Q31[t+(x+1)*num_t] + Q31[t+x*num_t] ) * ( p1x[(t+1)+x    *num_t] - p1x[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q32&&p1y)?( ( Q32[t+(x+1)*num_t] + Q32[t+x*num_t] ) * ( p1y[(t+1)+x    *num_t] - p1y[t    +x    *num_t] ) ): 0 )
                                                                        + ( (Q33&&p1z)?( ( Q33[t+(x+1)*num_t] + Q33[t+x*num_t] ) * ( p1z[(t+1)+x    *num_t] - p1z[t    +x    *num_t] ) ): 0 )
                                                                        ) 
                                                                      - ( ( (Q31&&p1x)?( ( Q31[t+(x-1)*num_t] + Q31[t+x*num_t] ) * ( p1x[t    +x    *num_t] - p1x[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q32&&p1y)?( ( Q32[t+(x-1)*num_t] + Q32[t+x*num_t] ) * ( p1y[t    +x    *num_t] - p1y[(t-1)+x    *num_t] ) ): 0 )
                                                                        + ( (Q33&&p1z)?( ( Q33[t+(x-1)*num_t] + Q33[t+x*num_t] ) * ( p1z[t    +x    *num_t] - p1z[(t-1)+x    *num_t] ) ): 0 )
                                                                        )
                                                                      + 0.5 * ( ( ( (C31&&p1x)? ( C31[t  +(x+1)*num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t-1+(x+1)*num_t] ) - C31[t  +(x-1)*num_t] * ( p1x[t+1+(x-1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (C32&&p1y)? ( C32[t  +(x+1)*num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t-1+(x+1)*num_t] ) - C32[t  +(x-1)*num_t] * ( p1y[t+1+(x-1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (C33&&p1z)? ( C33[t  +(x+1)*num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t-1+(x+1)*num_t] ) - C33[t  +(x-1)*num_t] * ( p1z[t+1+(x-1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 )
                                                                                ) 
                                                                              + ( ( (G31&&p1x)? ( G31[t+1+x    *num_t] * ( p1x[t+1+(x+1)*num_t] - p1x[t+1+(x-1)*num_t] ) - G31[t-1+x    *num_t] * ( p1x[t-1+(x+1)*num_t] - p1x[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (G32&&p1y)? ( G32[t+1+x    *num_t] * ( p1y[t+1+(x+1)*num_t] - p1y[t+1+(x-1)*num_t] ) - G32[t-1+x    *num_t] * ( p1y[t-1+(x+1)*num_t] - p1y[t-1+(x-1)*num_t] ) ) :0 )
                                                                                + ( (G33&&p1z)? ( G33[t+1+x    *num_t] * ( p1z[t+1+(x+1)*num_t] - p1z[t+1+(x-1)*num_t] ) - G33[t-1+x    *num_t] * ( p1z[t-1+(x+1)*num_t] - p1z[t-1+(x-1)*num_t] ) ) :0 )
                                                                                )
                                                                              )
                                                                      )
                                  + 2 * p1z[t+x*num_t]
                                  -     p2z[t+x*num_t]
                                  ;
                }
            }
        }
        if ( absorbing_boundary_condition )
        {
            // absorbing boundary t = 0
            for ( int x(0)
                , t(0)
                ; x < num_x
                ; ++x
                )
            {
                cz[t+x*num_t] = ( 2 * p1z[t+x*num_t] 
                                    - p2z[t+x*num_t] 
                                    - cz [t+2+x*num_t] 
                                + 2 * p1z[t+2+x*num_t] 
                                    - p2z[t+2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cz[t+2+x*num_t] 
                                                          - p2z[t+2+x*num_t] 
                                                          + p2z[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1z[t+2+x*num_t] 
                                                                                - 2 * p1z[t+1+x*num_t] 
                                                                                    + p1z[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary t = num_t-1
            for ( int x(0)
                , t(num_t-1)
                ; x < num_x
                ; ++x
                )
            {
                cz[t+x*num_t] = ( 2 * p1z[t+x*num_t] 
                                    - p2z[t+x*num_t] 
                                    - cz [t-2+x*num_t] 
                                + 2 * p1z[t-2+x*num_t] 
                                    - p2z[t-2+x*num_t] 
                                    + vel[t+x*num_t] * h * ( cz[t-2+x*num_t] 
                                                          - p2z[t-2+x*num_t] 
                                                          + p2z[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1z[t-2+x*num_t] 
                                                                                - 2 * p1z[t-1+x*num_t] 
                                                                                    + p1z[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = 0
            for ( int x(0)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cz[t+x*num_t] = ( 2 * p1z[t+x*num_t] 
                                    - p2z[t+x*num_t] 
                                    - cz [t+(x+2)*num_t] 
                                + 2 * p1z[t+(x+2)*num_t] 
                                    - p2z[t+(x+2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cz[t+(x+2)*num_t] 
                                                          - p2z[t+(x+2)*num_t] 
                                                          + p2z[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1z[t+(x+2)*num_t] 
                                                                                - 2 * p1z[t+(x+1)*num_t] 
                                                                                    + p1z[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
            // absorbing boundary x = num_x-1
            for ( int x(num_x-1)
                , t(0)
                ; t < num_t
                ; ++t
                )
            {
                cz[t+x*num_t] = ( 2 * p1z[t+x*num_t] 
                                    - p2z[t+x*num_t] 
                                    - cz [t+(x-2)*num_t] 
                                + 2 * p1z[t+(x-2)*num_t] 
                                    - p2z[t+(x-2)*num_t] 
                                    + vel[t+x*num_t] * h * ( cz[t+(x-2)*num_t] 
                                                          - p2z[t+(x-2)*num_t] 
                                                          + p2z[t+x*num_t] 
                                                          ) 
                                    - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1z[t+(x-2)*num_t] 
                                                                                - 2 * p1z[t+(x-1)*num_t] 
                                                                                    + p1z[t+x*num_t] 
                                                                                    ) 
                                ) / ( 1 + vel[t+x*num_t] * h );
            }
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p2z[k] = p1z[k];
        }
        for ( int k(0)
            ; k < num_t * num_x
            ; ++k
            )
        {
            p1z[k] = cz[k];
        }
    }
}

void propagate ( float   h 
               , int     num_x 
               , int     num_t 
               , float * p2 
               , float * p1 
               , float * c 
               , float * vel 
               , bool absorbing_boundary_condition = true
               )
{
    float fact ( h * h );
    for ( int t(1)
        ; t+1 < num_t
        ; ++t
        )
    {
        for ( int x(1)
            ; x+1 < num_x
            ; ++x
            )
        {
            if ( fabs ( p2[t+x*num_t] ) > 1e-40 
              || fabs ( p1[t+1+x*num_t] ) > 1e-40 
              || fabs ( p1[t-1+x*num_t] ) > 1e-40 
              || fabs ( p1[t+(x+1)*num_t] ) > 1e-40 
              || fabs ( p1[t+(x-1)*num_t] ) > 1e-40 
               )
            {
                c[t+x*num_t] = fact * vel[t+x*num_t] * vel[t+x*num_t] * ( -4*p1[t+x*num_t] 
                                                                        +    p1[t+(x+1)*num_t]
                                                                        +    p1[t+(x-1)*num_t]
                                                                        +    p1[(t+1)+x*num_t]
                                                                        +    p1[(t-1)+x*num_t]
                                                                        )
                             + 2 * p1[t+x*num_t]
                             -     p2[t+x*num_t]
                             ;
            }
        }
    }
    if ( absorbing_boundary_condition )
    {
        // absorbing boundary t = 0
        for ( int x(0)
            , t(0)
            ; x < num_x
            ; ++x
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+2+x*num_t] 
                           + 2 * p1[t+2+x*num_t] 
                               - p2[t+2+x*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+2+x*num_t] 
                                                      - p2[t+2+x*num_t] 
                                                      + p2[t+x*num_t] 
                                                      ) 
                                - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+2+x*num_t] 
                                                                            - 2 * p1[t+1+x*num_t] 
                                                                                + p1[t+x*num_t] 
                                                                                ) 
                            ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary t = num_t-1
        for ( int x(0)
            , t(num_t-1)
            ; x < num_x
            ; ++x
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t-2+x*num_t] 
                           + 2 * p1[t-2+x*num_t] 
                               - p2[t-2+x*num_t] 
                               + vel[t+x*num_t] * h * ( c[t-2+x*num_t] 
                                                     - p2[t-2+x*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t-2+x*num_t] 
                                                                           - 2 * p1[t-1+x*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary x = 0
        for ( int x(0)
            , t(0)
            ; t < num_t
            ; ++t
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+(x+2)*num_t] 
                           + 2 * p1[t+(x+2)*num_t] 
                               - p2[t+(x+2)*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+(x+2)*num_t] 
                                                     - p2[t+(x+2)*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+(x+2)*num_t] 
                                                                           - 2 * p1[t+(x+1)*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
        }
        // absorbing boundary x = num_x-1
        for ( int x(num_x-1)
            , t(0)
            ; t < num_t
            ; ++t
            )
        {
            c[t+x*num_t] = ( 2 * p1[t+x*num_t] 
                               - p2[t+x*num_t] 
                               - c [t+(x-2)*num_t] 
                           + 2 * p1[t+(x-2)*num_t] 
                               - p2[t+(x-2)*num_t] 
                               + vel[t+x*num_t] * h * ( c[t+(x-2)*num_t] 
                                                     - p2[t+(x-2)*num_t] 
                                                     + p2[t+x*num_t] 
                                                     ) 
                               - 2 * vel[t+x*num_t] * vel[t+x*num_t] * h * h * ( p1[t+(x-2)*num_t] 
                                                                           - 2 * p1[t+(x-1)*num_t] 
                                                                               + p1[t+x*num_t] 
                                                                               ) 
                           ) / ( 1 + vel[t+x*num_t] * h );
        }
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p2[k] = p1[k];
    }
    for ( int k(0)
        ; k < num_t * num_x
        ; ++k
        )
    {
        p1[k] = c[k];
    }
}

inline float pow_2 ( float const & x )
{
    return x*x;
}

float find_max ( float const * x , int num )
{
    float m ( 0 ) , t;
    for ( int k(0)
        ; k < num
        ; ++k
        )
    {
        t = fabs(x[k]);
        if ( t > m )
        {
            m = t;
        }
    }
    return m;
}

struct ray
{
    float  x ;
    float  y ;
    float  z ;
    float vx ;
    float vy ;
    float vz ;
    float  t ;
    float th_init;
    ray(){}
    ray ( float  _x 
        , float  _y 
        , float  _z 
        , float _vx 
        , float _vy 
        , float _vz 
        )
    :  x (  _x )
    ,  y (  _y )
    ,  z (  _z )
    , vx ( _vx )
    , vy ( _vy )
    , vz ( _vz )
    ,  t (   0 )
    {

    }
};

int main( int argc , char ** argv )
{

    if (argc < 3)
    {
        std::cout << "wrong number of arguments" << std::endl;
        exit(1);
    }

    int arg = 0;

    arg++;

    std::string file_name( argv[arg] );

    arg++;

    std::string output_file_name( argv[arg] );

    int num_iter = get_argument ( arg , argc , argv );

    float x_frac = get_argument ( arg , argc , argv );

    float t_frac = get_argument ( arg , argc , argv );

    SEPReader reader ( file_name . c_str () );
    int num_t = reader . n1;
    int num_x = reader . n2 
              * reader . n3
              * reader . n4
              * reader . n5
              * reader . n6
              * reader . n7
              * reader . n8
              ;
    std::cout << num_t << std::endl;
    std::cout << num_x << std::endl;
    float * data = new float [ num_x * num_t ];
    memset ( &data[0] , 0 , num_x * num_t );
    reader . read_sepval ( & data [ 0 ] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    int num_x2 = num_x;
    int num_t2 = num_t;
    float * vel = new float [ num_x2 * num_t2 ];
    float * rho = new float [ num_x2 * num_t2 ];
    float * delta  = new float [ num_x2 * num_t2 ];
    float * epsilon= new float [ num_x2 * num_t2 ];
    float * nux    = new float [ num_x2 * num_t2 ];
    float * nuy    = new float [ num_x2 * num_t2 ];
    float * nuz    = new float [ num_x2 * num_t2 ];
    float * vs     = new float [ num_x2 * num_t2 ];

    for ( int x(0)
        , T
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t2
            ; ++t
            , ++k
            )
        {
            vs [k]     = 1500;
            vel[k]     = 3000;//data[t+x*num_t2];
            rho[k]     = 1/(vel[k]*vel[k]);
            delta[k]   = anis_d;
            epsilon[k] = anis_e;
            nux[k]     = sin(theta*M_PI/180);
            nuy[k]     = 0;
            nuz[k]     = cos(theta*M_PI/180);
        }
    }

    float scale_x ( 1 );
    float scale_y ( 1 );
    float scale_z ( 1 );

    velocity_model < float > model ( num_x 
                                   , 1 
                                   , num_t 
                                   , scale_x
                                   , scale_y 
                                   , scale_z 
                                   , vel
                                   , epsilon
                                   , delta
                                   , nux
                                   , nuy
                                   , nuz
                                   );



    raytracing :: RK4_solver < float 
                             // , raytracing :: isotropic_exact_ray_step_functor_rk4 < float >
                             , raytracing :: tti_exact_ray_step_functor_rk4 < float >
                             > solver;

    float * p2x= new float [ num_x * num_t2 ]; memset ( &p2x[0] , 0 , num_x2 * num_t2 );
    float * p2y= new float [ num_x * num_t2 ]; memset ( &p2y[0] , 0 , num_x2 * num_t2 );
    float * p2z= new float [ num_x * num_t2 ]; memset ( &p2z[0] , 0 , num_x2 * num_t2 );
    float * p2d= new float [ num_x * num_t2 ]; memset ( &p2d[0] , 0 , num_x2 * num_t2 );
    float * p1x= new float [ num_x * num_t2 ]; memset ( &p1x[0] , 0 , num_x2 * num_t2 );
    float * p1y= new float [ num_x * num_t2 ]; memset ( &p1y[0] , 0 , num_x2 * num_t2 );
    float * p1z= new float [ num_x * num_t2 ]; memset ( &p1z[0] , 0 , num_x2 * num_t2 );
    float * p1d= new float [ num_x * num_t2 ]; memset ( &p1d[0] , 0 , num_x2 * num_t2 );
    float * cx = new float [ num_x * num_t2 ]; memset ( &cx [0] , 0 , num_x2 * num_t2 );
    float * cy = new float [ num_x * num_t2 ]; memset ( &cy [0] , 0 , num_x2 * num_t2 );
    float * cz = new float [ num_x * num_t2 ]; memset ( &cz [0] , 0 , num_x2 * num_t2 );
    float * cd = new float [ num_x * num_t2 ]; memset ( &cd [0] , 0 , num_x2 * num_t2 );
    float * A11= new float [ num_x * num_t2 ]; memset ( &A11[0] , 0 , num_x2 * num_t2 );
    float * A12= new float [ num_x * num_t2 ]; memset ( &A12[0] , 0 , num_x2 * num_t2 );
    float * A13= new float [ num_x * num_t2 ]; memset ( &A13[0] , 0 , num_x2 * num_t2 );
    float * A21= new float [ num_x * num_t2 ]; memset ( &A21[0] , 0 , num_x2 * num_t2 );
    float * A22= new float [ num_x * num_t2 ]; memset ( &A22[0] , 0 , num_x2 * num_t2 );
    float * A23= new float [ num_x * num_t2 ]; memset ( &A23[0] , 0 , num_x2 * num_t2 );
    float * A31= new float [ num_x * num_t2 ]; memset ( &A31[0] , 0 , num_x2 * num_t2 );
    float * A32= new float [ num_x * num_t2 ]; memset ( &A32[0] , 0 , num_x2 * num_t2 );
    float * A33= new float [ num_x * num_t2 ]; memset ( &A33[0] , 0 , num_x2 * num_t2 );
    float * Q11= new float [ num_x * num_t2 ]; memset ( &Q11[0] , 0 , num_x2 * num_t2 );
    float * Q12= new float [ num_x * num_t2 ]; memset ( &Q12[0] , 0 , num_x2 * num_t2 );
    float * Q13= new float [ num_x * num_t2 ]; memset ( &Q13[0] , 0 , num_x2 * num_t2 );
    float * Q21= new float [ num_x * num_t2 ]; memset ( &Q21[0] , 0 , num_x2 * num_t2 );
    float * Q22= new float [ num_x * num_t2 ]; memset ( &Q22[0] , 0 , num_x2 * num_t2 );
    float * Q23= new float [ num_x * num_t2 ]; memset ( &Q23[0] , 0 , num_x2 * num_t2 );
    float * Q31= new float [ num_x * num_t2 ]; memset ( &Q31[0] , 0 , num_x2 * num_t2 );
    float * Q32= new float [ num_x * num_t2 ]; memset ( &Q32[0] , 0 , num_x2 * num_t2 );
    float * Q33= new float [ num_x * num_t2 ]; memset ( &Q33[0] , 0 , num_x2 * num_t2 );
    float * C11= new float [ num_x * num_t2 ]; memset ( &C11[0] , 0 , num_x2 * num_t2 );
    float * C12= new float [ num_x * num_t2 ]; memset ( &C12[0] , 0 , num_x2 * num_t2 );
    float * C13= new float [ num_x * num_t2 ]; memset ( &C13[0] , 0 , num_x2 * num_t2 );
    float * C21= new float [ num_x * num_t2 ]; memset ( &C21[0] , 0 , num_x2 * num_t2 );
    float * C22= new float [ num_x * num_t2 ]; memset ( &C22[0] , 0 , num_x2 * num_t2 );
    float * C23= new float [ num_x * num_t2 ]; memset ( &C23[0] , 0 , num_x2 * num_t2 );
    float * C31= new float [ num_x * num_t2 ]; memset ( &C31[0] , 0 , num_x2 * num_t2 );
    float * C32= new float [ num_x * num_t2 ]; memset ( &C32[0] , 0 , num_x2 * num_t2 );
    float * C33= new float [ num_x * num_t2 ]; memset ( &C33[0] , 0 , num_x2 * num_t2 );
    float * G11= new float [ num_x * num_t2 ]; memset ( &G11[0] , 0 , num_x2 * num_t2 );
    float * G12= new float [ num_x * num_t2 ]; memset ( &G12[0] , 0 , num_x2 * num_t2 );
    float * G13= new float [ num_x * num_t2 ]; memset ( &G13[0] , 0 , num_x2 * num_t2 );
    float * G21= new float [ num_x * num_t2 ]; memset ( &G21[0] , 0 , num_x2 * num_t2 );
    float * G22= new float [ num_x * num_t2 ]; memset ( &G22[0] , 0 , num_x2 * num_t2 );
    float * G23= new float [ num_x * num_t2 ]; memset ( &G23[0] , 0 , num_x2 * num_t2 );
    float * G31= new float [ num_x * num_t2 ]; memset ( &G31[0] , 0 , num_x2 * num_t2 );
    float * G32= new float [ num_x * num_t2 ]; memset ( &G32[0] , 0 , num_x2 * num_t2 );
    float * G33= new float [ num_x * num_t2 ]; memset ( &G33[0] , 0 , num_x2 * num_t2 );

    float * vp     = vel;
    float * gamma  = NULL;

    std::cout << "p1" << std::endl;
    convert_Thomson_parameters_to_moduli_tti ( num_x2 * num_t2
                                             , A11
                                             , A12
                                             , A13
                                             , A21
                                             , A22
                                             , A23
                                             , A31
                                             , A32
                                             , A33
                                             , Q11
                                             , Q12
                                             , Q13
                                             , Q21
                                             , Q22
                                             , Q23
                                             , Q31
                                             , Q32
                                             , Q33
                                             , C11
                                             , C12
                                             , C13
                                             , C21
                                             , C22
                                             , C23
                                             , C31
                                             , C32
                                             , C33
                                             , G11
                                             , G12
                                             , G13
                                             , G21
                                             , G22
                                             , G23
                                             , G31
                                             , G32
                                             , G33
                                             , nux
                                             , nuy
                                             , nuz
                                             , rho      
                                             , vp
                                             , delta   
                                             , epsilon 
                                             , vs      
                                             , gamma   
                                             );
    std::cout << "p2" << std::endl;


    float width = 0.5;

    for ( int x(0)
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t2
            ; ++t
            , ++k
            )
        {
            float disp (sqrtf(pow_2(x-num_x*x_frac)+pow_2(t-num_t*t_frac)));
            p2x[ k ] = ((x>num_x*x_frac)?1:-1)*0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            p2y[ k ] = 0;
            p2z[ k ] = ((t>num_t*t_frac)?1:-1)*0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            // p2x[ k ] = ((x>num_x*x_frac)?1:-1)*0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            // p2y[ k ] = 0*0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            // p2z[ k ] = ((t>num_t*t_frac)?1:-1)*0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            p2d[ k ] = 0.05*exp(-disp*disp/(128.0*width))*sin(2*M_PI*disp/(16.0*sqrt(width)));
            p1x[ k ] = 0;
            p1y[ k ] = 0;
            p1z[ k ] = 0;
            p1d[ k ] = 0;
            cx [ k ] = 0;
            cy [ k ] = 0;
            cz [ k ] = 0;
            cd [ k ] = 0;
        }
    }
    std::cout << "p3" << std::endl;

    float h = 0.00020;

    for ( int k(0) ; k < num_iter ; ++k )
    {
        propagate_anisotropic ( h 
                              , num_x2 
                              , num_t2
                              , p2x
                              , p2y
                              , p2z
                              , p1x
                              , p1y
                              , p1z
                              , cx
                              , cy
                              , cz
                              , vel 
                              , rho
                              , A11
                              , A12
                              , A13
                              , A21
                              , A22
                              , A23
                              , A31
                              , A32
                              , A33
                              , Q11
                              , Q12
                              , Q13
                              , Q21
                              , Q22
                              , Q23
                              , Q31
                              , Q32
                              , Q33
                              , C11
                              , C12
                              , C13
                              , C21
                              , C22
                              , C23
                              , C31
                              , C32
                              , C33
                              , G11
                              , G12
                              , G13
                              , G21
                              , G22
                              , G23
                              , G31
                              , G32
                              , G33
                              , false
                              );
        propagate ( h 
                  , num_x2 
                  , num_t2
                  , p2d
                  , p1d
                  , cd
                  , vel 
                  );
    }

    std::vector < ray > rays;

    for ( float th(-M_PI); th <= M_PI; th += M_PI/1000 )
    {
        ray r ( scale_x*num_x*x_frac 
              , scale_y*0 
              , scale_z*num_t*t_frac
              , sin(th)
              , 0 
              , cos(th) 
              );
        r . th_init = th;
        float slow ( 1 / model ( r . x 
                               , r . y
                               , r . z
                               , r . vx
                               , r . vy
                               , r . vz
                               ) 
                   );
        r . vx *= slow;
        r . vy *= slow;
        r . vz *= slow;
        rays . push_back ( r );
    }

    for ( int k(0) ; k + 5 < num_iter ; ++k )
    {
        for ( int i(0) ; i < rays . size () ; ++i )
        {
            ray & r ( rays [ i ] );
            solver . RK4 ( h
                         , r . x 
                         , r . vx
                         , r . y 
                         , r . vy
                         , r . z 
                         , r . vz
                         , r . t
                         , model
                         );
            if ( i % 100 == 0 )
            {
                int z_ind ( (int)(r.z/scale_z) );
                int x_ind ( (int)(r.x/scale_x) );
                if ( z_ind >= 0 && z_ind < num_t && x_ind >= 0 && x_ind < num_x )
                data[ z_ind + num_t * x_ind ] = 0;
            }
        }
    }

    std::cout << "p4" << std::endl;

    float inv_max_w = 1 / ( find_max ( cx , num_x*num_t2 ) + find_max ( cy , num_x*num_t2 ) + find_max ( cz , num_x*num_t2 ) ) ;
    float max_v = find_max ( data   , num_x*num_t  );

    for ( int i(0) ; i < rays . size () ; ++i )
    {
        ray & r ( rays [ i ] );
        int z_ind ( (int)(r.z/scale_z) );
        int x_ind ( (int)(r.x/scale_x) );
        if ( z_ind >= 0 && z_ind < num_t && x_ind >= 0 && x_ind < num_x )
        data[ z_ind + num_t * x_ind ] = 0;
    }

    for ( int x(0)
        , k(0)
        ; x < num_x
        ; ++x
        )
    {
        for ( int t(0)
            ; t < num_t
            ; ++t
            , ++k
            )
        {
            data[k] -= 0.5 * (cx[t+x*num_t2] + cy[t+x*num_t2] + cz[t+x*num_t2]) * max_v * inv_max_w;
        }
    }

    SEPWriter writer ( output_file_name . c_str () 
                     , reader . o1 , reader . d1 , num_t
                     , reader . o2 , reader . d2 , num_x
                     , reader . o3 , reader . d3 , 1
                     , reader . get_header_labels ()
                     , reader . get_sort_order ()
                     , (output_file_name + std::string("@")) . c_str()
                     );

    writer . OpenDataFile ( (output_file_name + std::string("@")) . c_str() );

    writer . write_sepval ( (float*)&data[0] , reader . o1 , reader . o2 , reader . o3 , num_x * num_t );

    return 0;

}


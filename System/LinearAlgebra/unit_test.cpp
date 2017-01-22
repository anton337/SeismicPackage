#include <math.h>
#include <iostream>
#include <fftw3.h>
#include "solvers.h"

int main()
{
    SymmetricConjugateGradientSolver < matrix < float > , float > solver_sym;
    GeneralConjugateGradientSolver < float > solver_asym;
    float * data = new float [9];
    data[0] = 4;
    data[1] = 3;
    data[2] = 2;
    data[3] = 3;
    data[4] = 4;
    data[5] = 3;
    data[6] = 2;
    data[7] = 3;
    data[8] = 4;
    matrix < float > A ( 3 , 3 , data );
    float * b = new float [3];
    b[0] = 1;
    b[1] = 2;
    b[2] = 3;
    float * x = new float [3];
    float * x_correct = new float [3];
    float * x1 = new float [3];
    solver_sym ( A , x , b );
    solver_sym ( A , x_correct , b );
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    data[2] += 0.5;
    solver_asym ( A , x , b );
    solver_asym ( A , x1 , b );
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    GradientDescentSolver < float , LinearSystemEnergyL2 < float > > solver_grad;
    data[2] -= 0.5;
    LinearSystemEnergyL2 < float > energy ( A , b );
    solver_grad ( energy , x , 3 );
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    QuasiNewtonsMethodSolver < float , LinearSystemEnergyL2 < float > > solver_quas;
    data[2] += 0.5;
    x[0] = 0;
    x[1] = 0;
    x[2] = 0;
    solver_quas ( energy , x , 3 );
    std::cout << x[0] << " " << x[1] << " " << x[2] << std::endl;
    std::cout << x_correct[0] << " " << x_correct[1] << " " << x_correct[2] << std::endl;
    std::cout << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;

    fftwf_complex * data_c = new fftwf_complex [ 3 ];
    
    data_c[0][0] = 1;
    data_c[0][1] = 0;
    data_c[1][0] = 2;
    data_c[1][1] =-3;
    data_c[2][0] = 4;
    data_c[2][1] =-5;

    fftwf_complex * b_c = new fftwf_complex [ 3 ];
    fftwf_complex * x_c = new fftwf_complex [ 3 ];

    b_c[0][0] = 1;
    b_c[0][1] = 2;
    b_c[1][0] = 3;
    b_c[1][1] = 4;
    b_c[2][0] = 5;
    b_c[2][1] = 6;

    SymmetricConjugateGradientSolver < toeplitz_matrix < fftwf_complex > , fftwf_complex > solver_sym_c;

    toeplitz_matrix < fftwf_complex > A_c ( 3 , 3 , data_c );

    solver_sym_c ( A_c , x_c , b_c );

    std::cout << x_c[0][0] << " " << x_c[0][1] << "       " 
              << x_c[1][0] << " " << x_c[1][1] << "       "
              << x_c[2][0] << " " << x_c[2][1] << std::endl;

    return 0;
}


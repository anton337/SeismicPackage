#include <iostream>
#include <fftw3.h>
#include <math.h>

void hilbert ( int n , double * original_input , double * hilbert_output )
{
    fftw_complex * in  = new fftw_complex[n];
    fftw_complex * out = new fftw_complex[n];
    fftw_complex * hlb  = new fftw_complex[n];

    for ( int i=0; i < n; i++ )
    {
        in[i][0] = original_input[i];
        in[i][1] = 0;
    }

    fftw_plan fwd_plan = fftw_plan_dft_1d ( n , in , out , FFTW_FORWARD , FFTW_ESTIMATE );
    fftw_execute ( fwd_plan );
    double tmp;
    for ( int f(0) ; f < n ; f++ )
    {
        if ( f < n/2 )
        {
            out[f][1] *= -1;
        }
        else
        if ( f == n/2 )
        {
            out[f][0] = 0;
            out[f][1] = 0;
        }
        else
        if ( f > n/2 )
        {
            out[f][0] *= -1;
        }
        else
        if ( f == 0 )
        {
            out[f][0] = 0;
            out[f][1] = 0;
        }
        tmp = out[f][1];
        out[f][1] = -out[f][0];
        out[f][0] = -tmp;
    }
    fftw_plan inv_plan = fftw_plan_dft_1d ( n , out , hlb , FFTW_BACKWARD , FFTW_ESTIMATE );
    fftw_execute ( inv_plan );
    for ( int i(0) ; i < n ; i++ )
    {
        hlb[i][0] /= n;
        hlb[i][1] /= n;
        hilbert_output[i] = hlb[i][0];
    }

}

double angle(double x1,double y1,double x2,double y2)
{
    double phase1 = atan2(x1,y1);
    double phase2 = atan2(x2,y2);
    return atan(tan(phase1 - phase2));
}

double envelope(double x,double y)
{
    return sqrt(x*x+y*y);
}

double frequency(int n,double * x,double * y,double * freq)
{
    double factor ( 1.0 / (2*M_PI*1) );
    for ( int i(0) ; i+1 < n ; i++ )
    {
        freq[i] = atan2(x[i]*y[i+1] - x[i+1]*y[i],x[i]*x[i+1] + y[i]*y[i+1]) * factor;
    }
}

int main()
{
    std::cout << "Hello World - Hilbert Transform" << std::endl;
    int n = 1024;
    double * in1 = new double[n];
    double * in2 = new double[n];
    double * freq = new double[n];
    double f = 1.0 / 128;
    for ( int i(0) ; i < n ; i++ )
    {
        in1[i] = cos(2*M_PI*i*f);
        //in2[i] = cos(2*M_PI*i*0.1+0.3-0.6*(double)i/(double)n);
        //in2[i] = cos(2*M_PI*i*0.01+0.3);
        in2[i] = cos(2*M_PI*i*f+5*(f*2*M_PI));
        //std::cout << in1[i] << "\t" << in2[i] << std::endl;
    }
    double * hlb1 = new double[n];
    double * hlb2 = new double[n];
    double * ang = new double[n];
    hilbert(n,in1,hlb1);
    hilbert(n,in2,hlb2);
    frequency(n,in1,hlb1,freq);
    for ( int i(0); i < n; i++ )
    {
        ang[i] = angle(in1[i],hlb1[i],in2[i],hlb2[i]);
        //std::cout << ang[i] << "\t" << freq[i] << std::endl;
        std::cout << ((freq[i]!=0)?(ang[i] / (2*M_PI*freq[i])):0) << std::endl;
    }
    return 0;
}


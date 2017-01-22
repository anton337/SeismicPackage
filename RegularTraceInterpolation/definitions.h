/**********************************************************************************//**
 *
 *       Filename:  definitions.h
 *
 *    Description:  variables for module.
 *
 *        Version:  1.0
 *        Created:  04/04/2012 8:00:00 AM
 *       Revision:  $id$
 *       Compiler:  g++
 *
 *         Author:  Norman Bramblett (), norman.bramblett@iongeo.com
 *                  Jochen Heyd      (), jochen.heyd     @iongeo.com
 *        Company:  ION
 *
 **************************************************************************************/
#ifndef  DEFINITIONS_H_
#define  DEFINITIONS_H_

#define module_type                 FKinterp2D

#include <map>

#include <boost/shared_ptr.hpp>

#include "thetis/ThetisModuleTypes.h"
#include "utils/FFTWInterface.h"
#include "utils/FFTUtils.hpp"

namespace utils  = tpcommon::utils;
namespace linalg = tpcommon::linalg;
namespace stlint = tpcommon::stlint;
namespace thetis = tpcommon::thetis;

/***********************************************************************************//**
 *     \struct:  Tile
 *      \brief:  Defines tile location, and other attributes.
 *               Also, contains the tile data
 **************************************************************************************/
struct Tile
{

    typedef tpcommon::linalg::alg_vector < float > vector_type;

    vector_type   m_array;
    int           m_x_index;
    int           m_t_index;
    int           m_x_width;
    int           m_t_width;

    Tile ( int p_x_index
         , int p_t_index
         , int p_x_width
         , int p_t_width
         )
    : m_x_index ( p_x_index )
    , m_t_index ( p_t_index )
    , m_x_width ( p_x_width )
    , m_t_width ( p_t_width )
    {
        m_array . resize (m_x_width*m_t_width);
    }

    Tile ()
    {

    }

};

/***********************************************************************************//**
 *      \class:  FKInterpolator
 *      \brief:  Output / Client class for NetworkInterface
 **************************************************************************************/
class FKInterpolator
{

  
    /*-----------------------EXCEPTIONS-------------------------------------------*/

public:
    /***********************************************************************************//**
     *       \class  FKInterpolator
     *       \brief  Parent exception class
     ***************************************************************************************/
    class FKInterpolator_error : public std::runtime_error
    {

    public:

        explicit
        FKInterpolator_error ( std::string const & p_arg )
        : std::runtime_error( p_arg ){}
    };

        
    /*-----------------------TYPEDEFS-------------------------------------------------*/

public:

    typedef thetis::ThetisModuleTypes::record_type                  record_type         ;
    typedef thetis::ThetisModuleTypes::record_field_type            record_field_type   ;
    typedef thetis::ThetisModuleTypes::float_array_type             float_array_type    ;

    typedef utils::FFTWInterface< float, std::complex< float > >    fft_interface_type  ;
    typedef fft_interface_type::time_vector     ::iterator          time_iterator       ;
    typedef fft_interface_type::frequency_vector::iterator          frequency_iterator  ;
    typedef tpcommon::linalg::alg_vector < float >                  vector_type         ;
    typedef tpcommon::linalg::alg_vector < float >::iterator        vector_iterator     ;

    /*-----------------------LIFECYCLE------------------------------------------------*/

public:

    /***********************************************************************************//**
     *      \class:  FKInterpolator
     *         \fn:  Constructor
     *      \brief:  Constructs interpolator including FFTInterfaces etc.
     **************************************************************************************/
    FKInterpolator( std::size_t p_L
                  , std::size_t p_trace_count
                  , std::size_t p_sample_count
                  )
    : m_L               ( p_L            )
    , m_trace_count     ( p_trace_count  )
    , m_sample_count    ( p_sample_count )
    {
        // Live traces
        m_live_trace_count = m_trace_count;
        {
            std::size_t u_groups ( m_live_trace_count / m_L );
            m_live_trace_count -= u_groups * ( m_L - 1 );
        }

        m_full_sample_count = m_L * m_sample_count;

        // FFT padding for efficiency
        m_padded_live_trace_count  = 
            utils::FFTUtils::compute_optimal_size< float, std::complex< float > >( m_live_trace_count              );
        m_padded_sample_count      = 
            utils::FFTUtils::compute_optimal_size< float, std::complex< float > >( m_sample_count                  );
        m_padded_trace_count       = 
            utils::FFTUtils::compute_optimal_size< float, std::complex< float > >( m_L * m_padded_live_trace_count );
        m_padded_full_sample_count = 
            utils::FFTUtils::compute_optimal_size< float, std::complex< float > >( m_L * m_padded_sample_count     );

        // Sizes in FK domain
        m_k_count  =   m_padded_live_trace_count             ;
        m_f_count  = ( m_padded_sample_count       >> 1 ) + 1;
        m_Lk_count =   m_padded_trace_count                  ;
        m_Lf_count = ( m_padded_full_sample_count  >> 1 ) + 1;

        // Set sizes for the FFTs
        linalg::alg_vector< int > v_input_size     ( 2 );
        linalg::alg_vector< int > v_input_repl_size( 2 );
        linalg::alg_vector< int > v_output_size    ( 2 );
        v_input_size      [ 0 ] = m_padded_live_trace_count;
        v_input_size      [ 1 ] = m_padded_sample_count;
        v_input_repl_size [ 0 ] = m_padded_trace_count;
        v_input_repl_size [ 1 ] = m_padded_full_sample_count;
        v_output_size     [ 0 ] = m_padded_trace_count;
        v_output_size     [ 1 ] = m_padded_sample_count;

        // Set up FFTInterfaces
        m_input_interface     .get_exclusive_mutex().lock();
        m_input_repl_interface.get_exclusive_mutex().lock();
        m_output_interface    .get_exclusive_mutex().lock();

        m_input_interface     .set_fft_size( v_input_size      );
        m_input_repl_interface.set_fft_size( v_input_repl_size );
        m_output_interface    .set_fft_size( v_output_size     );

        m_input_interface     .create_plans();
        m_input_repl_interface.create_plans();
        m_output_interface    .create_plans();

        // Allocate temporary arrays
        m_temp_1.resize( m_Lk_count * m_Lf_count );
        m_temp_2.resize( m_Lk_count * m_Lf_count );
        m_H     .resize( m_Lk_count *  m_f_count );
    }

    /***********************************************************************************//**
     *      \class:  FKInterpolator
     *         \fn:  Destructor
     *      \brief:  Mainly unlocks FFT interfaces
     **************************************************************************************/
    ~FKInterpolator()
    {
        m_input_interface     .get_exclusive_mutex().unlock();
        m_input_repl_interface.get_exclusive_mutex().unlock();
        m_output_interface    .get_exclusive_mutex().unlock();
    }

    /*-----------------------METHODS--------------------------------------------------*/




    /***********************************************************************************//**
     *      \class:  FKInterpolator
     *         \fn:  interpolate
     *      \brief:  Interpolates the traces
     **************************************************************************************/
    void interpolate( vector_type         & p_input_output
                    )
    {

        if ( p_input_output . size () != m_trace_count * m_sample_count )
        {
            throw FKInterpolator_error ( "Input vector should be of size trace_count * sample_count." );
        }

        // Zero everything
        m_input_interface     .zero_time_buffer     ();
        m_input_repl_interface.zero_time_buffer     ();
        m_output_interface    .zero_frequency_buffer();

        tpcommon::tpmemset( &( m_temp_2[ 0 ] ), 0, sizeof( float ) * m_Lk_count * m_Lf_count );

        m_output_interface.zero_frequency_buffer();

        // Copy input data into FFT interfaces
        {
            time_iterator   input_time_iter     ( m_input_interface     .  get_time_buffer_begin() );
            time_iterator   input_repl_time_iter( m_input_repl_interface.  get_time_buffer_begin() );
            vector_iterator input_begin_iter    ( p_input_output        .                  begin() );
            vector_iterator input_end_iter      ( p_input_output        . begin() + m_sample_count );
            for( std::size_t x( 0 )
               ; x < m_trace_count
               ; ++x
               )
            {
                {
                    if ( x % m_L == 0 )
                    {
                        // Dump trace data into FFTWInterfaces
                        std::copy( input_begin_iter, input_end_iter, input_time_iter      );
                        std::copy( input_begin_iter, input_end_iter, input_repl_time_iter );

                        // Increment FFT interface iters by padded trace lengths
                        input_time_iter      += m_padded_sample_count;
                        input_repl_time_iter += m_padded_full_sample_count;
                    }
                    input_begin_iter     += m_sample_count;
                    if ( x+1 >= m_trace_count )
                    {
                        input_end_iter        = p_input_output . end();
                    }
                    else
                    {
                        input_end_iter       += m_sample_count;
                    }    
                }
            }
        }

        // Transform input
        m_input_interface     .transform();
        m_input_repl_interface.transform();

        // Generate the larger complex via L copies of the transformed input
        {
            fft_interface_type::frequency_vector::iterator  input_freq_begin(  m_input_interface.get_frequency_buffer_begin() );
            fft_interface_type::frequency_vector::iterator output_freq_begin( m_output_interface.get_frequency_buffer_begin() );
            fft_interface_type::frequency_vector::iterator  input_freq_iter;
            fft_interface_type::frequency_vector::iterator output_freq_iter;
            for( std::size_t k( 0 ); k < m_k_count; ++k ) {

                input_freq_iter = input_freq_begin + ( k * m_f_count );

                for( std::size_t l( 0 ); l < m_L; ++l ) {

                    output_freq_iter = output_freq_begin + ( l * m_k_count + k ) * m_f_count;

                    std::copy( input_freq_iter, input_freq_iter + m_f_count, output_freq_iter );
                }
            }
        }

        // Construct filter: First temporary
        {
            float           const   f_thres( 1.0e-7 );

            stlint::abs   < float > f_abs;
            stlint::sqrt  < float > f_sqrt;
            stlint::square< float > f_square;

            frequency_iterator src_iter( m_input_repl_interface.get_frequency_buffer_begin() );
            frequency_iterator end_iter( m_input_repl_interface.get_frequency_buffer_end  () );
                 time_iterator     iter( m_temp_1.begin()                                    );
            float              f_re, f_im, f_temp;
            for ( ; src_iter != end_iter; ++iter, ++src_iter ) {

                f_re = f_abs( src_iter->real() );
                f_im = f_abs( src_iter->imag() );

                       if ( f_re < f_thres ) {

                    *iter  = f_im;

                } else if ( f_im < f_thres ) {

                    *iter  = f_re;

                } else if ( f_re > f_im    ) {

                    f_temp = f_im / f_re;
                    *iter  = f_re * f_sqrt( 1.0f + f_square( f_temp ) );

                } else {

                    f_temp = f_re / f_im;
                    *iter  = f_im * f_sqrt( 1.0f + f_square( f_temp ) );
                }
            }
        }

        // Construct filter: Second temporary (apparently some wrap-around sum)
        {
            float const f_L_inv( 1.0f / static_cast< float >( m_L ) );

            std::size_t ind_0;
            std::size_t ind_1;
            for ( std::size_t k( 0 ); k < m_Lk_count; ++k ) {
                ind_0 = k*m_Lf_count;
                for ( std::size_t l( 0 ); l < m_L       ; ++l ) {
                    if ( k + l*m_k_count < m_Lk_count ) {
                        ind_1 = ( k + l*m_k_count              ) * m_Lf_count;
                    }
                    else
                    {
                        ind_1 = ( k + l*m_k_count - m_Lk_count ) * m_Lf_count;
                    }
                    for ( std::size_t f( 0 ); f < m_Lf_count; ++f ) {
                        m_temp_2[ ind_0 + f ] += m_temp_1[ ind_1 + f ];
                    }
                }

                for ( std::size_t f( 0 ); f < m_Lf_count; ++f ) {
                    // Now scale
                    m_temp_2[ ind_0 + f ] *= f_L_inv;
                }
            }
        }

        // Find the some sort of "mean" of v_temp_2
        float f_mean( 0.0f );
        {
            stlint::abs< float        > f_abs;
            stlint::max< float, float > f_max;

            for ( std::size_t i( 0 ); i < m_temp_2.size(); ++i ) {

                f_mean = f_max( f_mean, f_abs( m_temp_2[ i ] ) );
            }

            // "Post-processing"
            f_mean *= 0.01f;
            f_mean  = f_max( f_mean, 0.01f );
            f_mean /= m_Lk_count * m_Lf_count;
        }

        // Construct filter: The actual thing (also accumulates a max used for scaling)
        float f_scale( 0.0 );
        {
            stlint::abs< float        > f_abs;
            stlint::max< float, float > f_max;

            float       f_H_max( 0.01f );

            std::size_t t_idx, h_idx;
            float       f_temp_2;
            for ( std::size_t k( 0 ); k < m_Lk_count; ++k ) {

                t_idx = k * m_Lf_count;
                h_idx = k * m_f_count;

                for ( std::size_t f( 0 ); f < m_f_count; ++f, ++t_idx, ++h_idx ) {

                    f_temp_2 = m_temp_2[ t_idx ];

                    // Protect against small numbers
                    if ( f_abs( f_temp_2 ) < f_mean ) {

                        if ( f_temp_2 >= 0.0f ) {
                            f_temp_2 =  f_mean;
                        } else {
                            f_temp_2 = -f_mean;
                        }
                    }

                    // Here's the filter
                    m_H[ h_idx ] = m_temp_1[ t_idx ] / f_temp_2;

                    // And its max
                    f_H_max      = f_max( f_H_max, f_abs( m_H[ h_idx ] ) );
                }
            }

            // Calc scale factor
            f_scale = static_cast< float >( m_L ) / ( f_H_max * m_H.size() );
        }

        // Finally, apply the filter
        {
            frequency_iterator begin_iter( m_output_interface.get_frequency_buffer_begin() );
            frequency_iterator       iter;
            std::size_t                                     h_idx;
            for ( std::size_t k( 0 ); k < m_Lk_count; ++k ) {

                iter  = begin_iter + k * m_f_count;
                h_idx =              k * m_f_count;

                for ( std::size_t f( 0 ); f < m_f_count; ++f, ++h_idx, ++iter ) {

                    *iter *= f_scale * m_H[ h_idx ];
                }
            }
        }

        // Invert output
        m_output_interface.invert();

        // Store interpolated traces
        {
            {
                time_iterator iter( m_output_interface.get_time_buffer_begin() );
                vector_iterator output_iter( p_input_output . begin () );
                for( std::size_t x( 0 )
                   ; x < m_trace_count 
                   ; ++x
                   )
                {
                    {
                        std::copy ( iter 
                                  , iter + m_sample_count
                                  , output_iter
                                  );
                    }
                    iter += m_padded_sample_count;
                    output_iter += m_sample_count;
                }
            }
        }

        // SAB
        {
            stlint::sqrt  < float > f_sqrt;
            float sample_count;
            vector_type rms(m_trace_count);
            for ( std::size_t s(0)
                , block_size (3)
                ; s < m_trace_count
                ; s += block_size
                )
            {
                for ( std::size_t x(0)
                    , k(s)
                    ; x < m_trace_count
                    ; ++x
                    )
                {
                    rms[x] = 0;
                    k = x*(m_sample_count) + s;
                    sample_count = 0;
                    for ( std::size_t t(s)
                        ; t < m_sample_count && t < s + block_size
                        ; ++t
                        , ++k
                        )
                    {
                        rms[x] += p_input_output[k]*(p_input_output[k]);
                        sample_count += 1;
                    }
                    rms[x] = f_sqrt(rms[x]/sample_count);
                }
                float factor;
                float scale;
                for ( std::size_t x(0)
                    , k(s)
                    ; x < m_trace_count 
                    ; ++x
                    )
                {
                    if ( x % m_L != 0 )
                    {
                        factor = static_cast<float>((x%m_L)*rms[m_L*(x/m_L+1)] + (m_L - static_cast<int>(x%m_L))*rms[m_L*(x/m_L)])/m_L;
                        if ( x + m_L >= m_trace_count+1 )
                        {
                            factor = static_cast<float>(rms[m_L*(x/m_L)]);
                        }
                        scale = (rms[x] > 1e-20) ? factor / rms[x] : 0.0;
                        k = x*(m_sample_count) + s;
                        for ( std::size_t t(s)
                            ; t < m_sample_count && t < s + block_size
                            ; ++t
                            , ++k
                            )
                        {
                            p_input_output[k] *= scale;
                        }
                    }
                    else
                    {
                        k += m_sample_count;
                    }
                }
            }
        }


    }

    /*-----------------------PRIVATE DATA---------------------------------------------*/

private:

    std::size_t                     m_L;

    std::size_t                     m_trace_count;
    std::size_t                     m_live_trace_count;
    std::size_t                     m_sample_count;
    std::size_t                     m_full_sample_count;

    std::size_t                     m_padded_trace_count;
    std::size_t                     m_padded_live_trace_count;
    std::size_t                     m_padded_sample_count;
    std::size_t                     m_padded_full_sample_count;

    std::size_t                     m_k_count;
    std::size_t                     m_f_count;
    std::size_t                     m_Lk_count;
    std::size_t                     m_Lf_count;

    fft_interface_type              m_input_interface;
    fft_interface_type              m_input_repl_interface;
    fft_interface_type              m_output_interface;

    fft_interface_type::time_vector m_temp_1;
    fft_interface_type::time_vector m_temp_2;
    fft_interface_type::time_vector m_H;

}; // -----  end of class FKInterpolator  -----

struct module_structures
{
    struct global_resources
    {
        std::string  s_pad_flag;  
        std::string  s_interp_flag;
        std::string  s_time_series_field;

        std::size_t  L;
    };

    struct worker_resources
    {
        typedef std::pair< std::size_t, std::size_t >   interpolation_type;
        typedef boost::shared_ptr< FKInterpolator >     interpolator_pointer;
        typedef std::map< interpolation_type
                        , interpolator_pointer
                        >                               interpolator_map_type;

        interpolator_map_type                           interpolator_map;

        boost::shared_ptr< FKInterpolator > interpolator;
    };
};

#endif   // ----- #ifndef DEFINITIONS_H_  -----

/** \file ADIDerotator.hpp
 * \author Jared R. Males
 * \brief Defines a generic ADI derotator class.
 * \ingroup hc_imaging_files
 *
 */

#include <mx/ioutils/fits/fitsHeader.hpp>

#ifndef ADIDerotator_hpp
#define ADIDerotator_hpp

#include <mx/math/geo.hpp>

namespace mx
{

namespace improc
{

/// A generic ADI derotator class.
/** This class is used to calculate the derotation angle for angular differential imaging.
 *
 * \ingroup hc_imaging
 *
 */
template <typename _realT, class verboseT=mx::verbose::vvv>
struct ADIDerotator
{
    typedef _realT realT;

    typedef mx::fits::fitsHeader<verboseT> fitsHeaderT;
    typedef mx::fits::fitsHeaderCard<verboseT> fitsHeaderCardT;

    /// Vector of keywords to extract from the fits headers
    std::vector<std::string> m_keywords;

    /// Vector(s) to hold the keyword values
    std::vector<realT> m_angles;

    std::string m_angleKeyword; ///< The keyword for the angle attribute.  Do not set this directly.

    /// Set the angle keyword
    /** Populates the kewords vector appropriately.
     */
    void angleKeyword( const std::string &akw /**< [in] The angle keyword */ )
    {
        m_angleKeyword = akw;
        m_keywords = { akw };
    }

    realT m_angleScale{ 0 };    ///< The scale to multiply the angle by
    realT m_angleConstant{ 0 }; ///< The constant to add to the scaled-angle.

    ADIDerotator()
    {
    }

    /// To allow ADIobservation to check for errors.
    bool isSetup()
    {
        if( ( m_angleKeyword == "" || m_keywords.size() == 0 ) || ( m_angleScale == 0 && m_angleConstant == 0 ) )
        {
            return false;
        }

        return true;
    }

    /// Method called by ADIobservation to get keyword-values
    /**
      * \returns mx::error_t::noerror on success, and \p bad will be empty.
      * \returns mx::error_t::error if any header fails conversion for \ref m_angleKeyword. \p bad will
      *          then contain the indices of \p heads for which conversion failed.
      */
    mx::error_t
    extractKeywords( std::vector<fitsHeaderT> &heads, /**< [in] The headers from the images being reduced.*/
                     std::vector<size_t> &bad /**< [in] indices of any headers which fail to convert */ )
    {
        m_angles.clear();
        return fits::headersToValues<realT>( m_angles, bad, heads, m_angleKeyword );
    }

    /// Calculate the derotation angle for a given image number
    /**
     * \returns the angle in radians by which to de-rotate the image c.c.w.
     */
    realT derotAngle( size_t imno /**< [in] the image number */ ) const
    {
        realT derot = m_angleScale * m_angles[imno] + m_angleConstant;
        derot = math::angleMod<math::degreesT<realT>>( derot );
        return math::dtor( derot );
    }
};

///@}

extern template struct ADIDerotator<float, mx::verbose::o>;
extern template struct ADIDerotator<float, mx::verbose::v>;
extern template struct ADIDerotator<float, mx::verbose::vv>;
extern template struct ADIDerotator<float, mx::verbose::vvv>;

extern template struct ADIDerotator<double, mx::verbose::o>;
extern template struct ADIDerotator<double, mx::verbose::v>;
extern template struct ADIDerotator<double, mx::verbose::vv>;
extern template struct ADIDerotator<double, mx::verbose::vvv>;

} // namespace improc
} // namespace mx

#endif // ADIDerotator_hpp

/*! 
 *   \file pstruct.h
 *   \brief Definition of the polar representation structure
 *
 *   Structure holding complex numbers in polar coordinates, i.e. m.e^(j.a) 
 *   where j=sqrt(-1)
*/ 
#ifndef EGS_Polar_H_Seen
#define EGS_Polar_H_Seen

#include <complex>
#include "math.h"
using namespace std;

/** 
 * A structure holding complex number in polar representation
 */
namespace egsim {

    typedef struct {
      double m ; /**< Magnitude */
      double a ; /**< Angle     */
    } pstruct ;

    pstruct c2p(const complex<double> a) ; /**< Convert complex to polar */
    complex<double> p2c(const pstruct a) ; /**< Convert polar to complex */

}

#endif

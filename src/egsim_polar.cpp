#include "egsim_pstruct.hpp"

egsim::pstruct egsim::c2p(const complex<double> a) { 

  egsim::pstruct b ;

  b.m = abs(a);
  b.a = arg(a);

  return ( b ) ;

}

complex<double> egsim::p2c(const egsim::pstruct a) { 
  return ( polar(a.m,a.a) ) ;
}

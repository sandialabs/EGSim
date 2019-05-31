#include <iostream>
#include <fstream>
#include "egsim.hpp"

using namespace std;

void egsim::buildAdmMat(Array1D<bus> &bdata, Array1D<line> &ldata, Array2D< complex<double> > &Y) {

  int nbus = (int) bdata.XSize() ;
  int nlin = (int) ldata.XSize() ;

  /* build y matrix */
  Y.Resize(nbus,nbus,complex<double>(0.0,0.0)) ;

  for ( int k = 0; k<nlin; k++ ) {
    if ( ldata(k).get_state() == 1 ) {
      int i = ldata(k).get_frbus();
      int j = ldata(k).get_tobus();
      int u = ldata(k).get_state();
      double chrg = ((double) u)*ldata(k).get_lchrg() ;
      complex<double> impd = ((double) u)*ldata(k).get_impd()  ;
      complex<double> ladm = ((double) u)*1.0/impd ;
      double tr  = ldata(k).get_tr();
      double psa = ldata(k).get_psa();
      complex<double> tap = polar(tr,psa);
      double tap2=tr*tr;
      Y(i,j) -= ladm/conj(tap) ;
      Y(j,i) -= ladm/tap ;
      Y(i,i) += ladm/tap2+complex<double>(0.0,0.5*chrg);
      Y(j,j) += ladm     +complex<double>(0.0,0.5*chrg);
    } /* done if line is on */
  } /* done loop over all lines */
  
  return ;

} /* done building the admittance matrix */

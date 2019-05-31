#include "egsim.hpp"

namespace egsim {

void calc_Gy(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array2D<double> &J) {

  calc_Gy_line(bdata, Yadm, y, J) ;
  calc_Gy_PQ  (bdata, y, J) ;
  calc_Gy_PV  (bdata, J) ;
  calc_Gy_SW  (bdata, J) ;

  return ;

}
/*
   ____                 _ _            
  / ___|_   _          | (_)_ __   ___ 
 | |  _| | | |  _____  | | | '_ \ / _ \
 | |_| | |_| | |_____| | | | | | |  __/
  \____|\__, |         |_|_|_| |_|\___|
        |___/                          
 */
void calc_Gy_line(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array2D<double> &J) {

  int nbus = (int) bdata.XSize() ;

  Array2D< complex<double> > dS(nbus,nbus,complex<double> (0.0,0.0));
  Array2D< complex<double> > dR(nbus,nbus,complex<double> (0.0,0.0));

  Array1D< complex<double> > Vc(nbus,complex<double> (0.0,0.0));
  for ( int i = 0; i < nbus; i++ ) Vc(i) = polar(y(i+nbus),y(i)) ;

  /* Compute dR - A in the text */
  for ( int i = 0; i < nbus; i++ ) {
      complex<double> csum = complex<double> (0.0,0.0);
      for ( int j = 0; j < nbus ; j++ ) csum += Yadm(i,j)*Vc(j);
      dR(i,i) = conj(Vc(i)) * csum ;
  }

  for ( int i = 0; i < nbus; i++ ) {
      for ( int j = 0; j < nbus ; j++ ) 
          dR(i,j) -= conj(Vc(i)) * Yadm(i,j) * Vc(j) ;
  }

  /* Compute dS - B in the text */
  for ( int i = 0; i < nbus; i++ ) {
    for ( int j = 0; j < nbus ; j++ ) 
      dS(i,j) = Vc(i)*conj(Yadm(i,j)*polar(1.0,y(j)));
  }
  for ( int i = 0; i < nbus; i++ ) {
    complex<double> csum = complex<double> (0.0,0.0);
    for ( int j = 0; j < nbus ; j++ ) 
      csum += Yadm(i,j)*Vc(j) ;
    dS(i,i) += conj(csum) * polar(1.0,y(i)) ;
  }

  /* 1st quarter */
  for ( int i = 0; i < nbus; i++ )
      for ( int j = 0; j < nbus ; j++ ) 
          J(i,j) = dR(i,j).imag() ;

  /* 2nd quarter */
  for ( int i = 0; i < nbus; i++ )
      for ( int j = 0; j < nbus ; j++ ) 
          J(i,j+nbus) = dS(i,j).real() ;

  /* 3rd quarter */
  for ( int i = 0; i < nbus; i++ )
      for ( int j = 0; j < nbus ; j++ ) 
          J(i+nbus,j) = dR(i,j).real() ;

  /* 4th quarter */
  for ( int i = 0; i < nbus; i++ )
      for ( int j = 0; j < nbus ; j++ ) 
          J(i+nbus,j+nbus) = dS(i,j).imag() ;

  return ;

}
/*
    ____                 ____   ___  
   / ___|_   _          |  _ \ / _ \ 
  | |  _| | | |  _____  | |_) | | | |
  | |_| | |_| | |_____| |  __/| |_| |
   \____|\__, |         |_|    \__\_\
         |___/                       
 */
void calc_Gy_PQ(Array1D< egsim::bus > &bdata, Array1D<double> &y, Array2D<double> &J) {

  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus ; i++ ) {

      if ( bdata(i).get_type() != PQld ) continue ;

      if ( bdata(i).get_shAdm() ) {
        double vmn2 = pow(bdata(i).get_Vmn(),2);
        J(i,     nbus+i) += 2.0 * bdata(i).get_pload() * y(nbus+i) / vmn2 ;
        J(nbus+i,nbus+i) += 2.0 * bdata(i).get_qload() * y(nbus+i) / vmn2 ;
      }

  }

  return ;
}

/*
    ____                 ______     __
   / ___|_   _          |  _ \ \   / /
  | |  _| | | |  _____  | |_) \ \ / / 
  | |_| | |_| | |_____| |  __/ \ V /  
   \____|\__, |         |_|     \_/   
         |___/                        
 */
void calc_Gy_PV(Array1D< egsim::bus > &bdata, Array2D<double> &J) {

  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus; i++ ) {

      if ( bdata(i).get_type() != PVgn ) continue ;

      for ( int j = 0; j < 2*nbus ; j++ ) {
          J(i+nbus,j) = 0.0;
          J(j,i+nbus) = 0.0;
      }
      J(i+nbus,i+nbus) = 1.0 ;

  }

  return ;

}

/*
  ____                 ______        __
 / ___|_   _          / ___\ \      / /
| |  _| | | |  _____  \___ \\ \ /\ / / 
| |_| | |_| | |_____|  ___) |\ V  V /  
 \____|\__, |         |____/  \_/\_/   
       |___/                           
 */
void calc_Gy_SW(Array1D< egsim::bus > &bdata, Array2D<double> &J) {

  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus; i++ ) {

      if ( bdata(i).get_type() != SW ) continue ;

      for ( int j = 0; j < 2*nbus ; j++ ) {
          J(i,j) = 0.0;
          J(j,i) = 0.0;
      }
      J(i,i) = 1.0 ;

  }

  for ( int i = 0; i < nbus; i++ ) {

      if ( bdata(i).get_type() != SW ) continue ;

      for ( int j = 0; j < 2*nbus ; j++ ) {
          J(i+nbus,j) = 0.0;
          J(j,i+nbus) = 0.0;
      }
      J(i+nbus,i+nbus) = 1.0 ;

  }

  return ;
}

}

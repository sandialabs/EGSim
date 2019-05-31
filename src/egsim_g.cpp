#include "egsim.hpp"

//#define DEBUG

namespace egsim {

void calc_g(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, 
            Array1D<double> &y, Array1D<double> &g) {

  calc_g_line(bdata, Yadm, y, g); /* build rhs - line data             */
  calc_g_PQ  (bdata,       y, g); /* build rhs - PQ load bus data      */
  calc_g_PV  (bdata,          g); /* build rhs - PV generator bus data */
  calc_g_SW  (bdata,          g); /* build rhs - SW bus data           */

  return ;

}
/*
                  _ _            
   __ _          | (_)_ __   ___ 
  / _` |  _____  | | | '_ \ / _ \
 | (_| | |_____| | | | | | |  __/
  \__, |         |_|_|_| |_|\___|
  |___/                          
*/
void calc_g_line(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, 
                 Array1D<double> &y, Array1D<double> &g) {

  int nbus = (int) bdata.XSize() ;

  /* build rhs - line data */
  Array1D< complex<double> > Vc(nbus);
  for ( int i = 0; i < nbus ; i++ ) 
    Vc(i) = polar(y(i+nbus),y(i));
  Array1D< complex<double> > VYV(nbus);
  for ( int i = 0; i < nbus ; i++ ) {
    VYV(i) = complex<double> (0.0,0.0);
    for ( int j = 0; j < nbus ; j++ ) 
      VYV(i) += Yadm(i,j) * Vc(j);
    VYV(i) = Vc(i) * conj( VYV(i) ) ;
  }
  for ( int i = 0; i < nbus ; i++ ) {
    g(i     ) = VYV(i).real();
    g(i+nbus) = VYV(i).imag();
  }

#ifdef DEBUG
    /* Display g based on line data */
  cout<<"--------------g-line--------------------"<<endl;
  for ( int i = 0; i < nbus ; i++ ) 
      cout << i+1 << " " << g(i) << g(i+nbus) << endl;
  cout<<"----------------------------------------"<<endl;
#endif

  return ;

}
/*
                    ____   ___  
     __ _          |  _ \ / _ \ 
    / _` |  _____  | |_) | | | |
   | (_| | |_____| |  __/| |_| |
    \__, |         |_|    \__\_\
    |___/                       
 */
void calc_g_PQ(Array1D< egsim::bus > &bdata, Array1D<double> &y, Array1D<double> &g) {

  /* build rhs - contributions from PQ loads */
  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus ; i++ ) {
      if ( bdata(i).get_type() != PQld ) continue ;

      g(i     ) += bdata(i).get_pload();
      g(i+nbus) += bdata(i).get_qload();

      if ( bdata(i).get_shAdm() ) {
          g(i     ) += bdata(i).get_pload()*pow(y(nbus+i),2)/pow(bdata(i).get_Vmn(),2)
	                    -bdata(i).get_pload();
          g(i+nbus) += bdata(i).get_qload()*pow(y(nbus+i),2)/pow(bdata(i).get_Vmn(),2)
                      -bdata(i).get_qload();
      }

  }

#ifdef DEBUG
    /* Display g after line and PQ buses */
  cout<<"--------------g-linePQ------------------"<<endl;
  for ( int i = 0; i < nbus ; i++ ) 
      cout << i+1 << " " << g(i) << g(i+nbus) << endl;
  cout<<"----------------------------------------"<<endl;
#endif

  return ;

}
/*
                   ______     __
    __ _          |  _ \ \   / /
   / _` |  _____  | |_) \ \ / / 
  | (_| | |_____| |  __/ \ V /  
   \__, |         |_|     \_/   
   |___/                        
*/
void calc_g_PV(Array1D< egsim::bus > &bdata, Array1D<double> &g) {

  /* build rhs - contributions from PV generators */
  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus ; i++ ) {

      if ( bdata(i).get_type() != PVgn ) continue ;

      g(i     ) += bdata(i).get_pload()-bdata(i).get_pgen();
      g(i+nbus)  = 0.0 ;

  }

  return ;
}
/*
                  ______        __
   __ _          / ___\ \      / /
  / _` |  _____  \___ \\ \ /\ / / 
 | (_| | |_____|  ___) |\ V  V /  
  \__, |         |____/  \_/\_/   
  |___/                           
*/
void calc_g_SW(Array1D< egsim::bus > &bdata, Array1D<double> &g) {

  /* build rhs - zero-out the slack contribution */
  int nbus = (int) bdata.XSize() ;

  for ( int i = 0; i < nbus ; i++ ) {

      if ( bdata(i).get_type() != SW ) continue ;

      g(i     ) = 0.0 ;
      g(i+nbus) = 0.0 ;

  }

  return ;
}

}
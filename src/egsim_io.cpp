#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace std;

#include "egsim.hpp"

void egsim::saveAdmMat(Array2D< complex<double> > &Y,string admout) {

  ofstream myfile, myfileR,myfileI;

  myfile.open  ("admmat.dat");
  myfileR.open ("admmatR.dat");
  myfileI.open ("admmatI.dat");
  if ( admout == string("full") ) {
    /* Output all entries */
    for ( int i = 0; i < (int) Y.XSize(); i++ )  {
      for ( int j = 0; j < (int) Y.YSize(); j++ ) {
        myfile<<Y(i,j)<<" ";
        myfileR<<Y(i,j).real()<<" ";
        myfileI<<Y(i,j).imag()<<" ";
      }
      myfile<<endl;
      myfileR<<endl;
      myfileI<<endl;
    }

  }
  else if ( admout == string("sparse") ) {
    /* Output non-zero entries */
    for ( int j = 0; j < (int) Y.YSize(); j++ )
      for ( int i = 0; i < (int) Y.XSize(); i++ )  {
        if (abs(Y(i,j))>0.001)
          myfile<<i+1<<","<<j+1<<": "<<Y(i,j)<<endl;
    }
  }
  else {
    std::cout<<"egsim::saveAdmMat() ERROR: Unknown output format"<<std::endl;
  }
  myfile.close();
  myfileR.close();
  myfileI.close();


  return ;

}

void egsim::saveStaticSol(Array1D< egsim::bus > &bdata) {

  int nbus = (int) bdata.XSize() ;

  ofstream myfile;
  myfile.open ("statsol.dat");
  for (  int i = 0; i<nbus; i++ )
    myfile << setw(3) << i+1 
           << std::fixed << setprecision(4) << setw(9) << bdata(i).get_vm() << " " << setprecision(4) << setw(9) << bdata(i).get_va()*180/PI << endl ;
  myfile.close();
  
  return ;

} /* done computing line powers */


void egsim::getLinePow(Array1D< egsim::bus > &bdata, Array1D< egsim::line > &ldata) {

  int nbus = (int) bdata.XSize() ;
  int nlin = (int) ldata.XSize() ;

  ofstream myfile;
  myfile.open ("linepow.dat");
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
          double tap2 = tr * tr;

          complex<double> Vi=bdata(i).get_vri();
          complex<double> Vj=bdata(j).get_vri();

          complex<double> MWs = Vi*conj(Vi*(ladm+complex<double>(0.0,0.5*chrg))/tap2 
                               -Vj*ladm/conj(tap));
          complex<double> MWr = Vj*conj(Vj*(ladm+complex<double>(0.0,0.5*chrg))
                               -Vi*ladm/tap);

          myfile<<setw(3)<<k+1<<"  "<<setw(3)<<i+1<<"  "<<setw(3)<<j+1<<"  "<<fixed
                << std::fixed << setprecision(4) << setw(9) 
                << MWs.real()<<"  "<<MWs.imag()<<"  "
                << MWr.real()<<"  "<<MWr.imag()<<endl;

      } /* done if line is on */
  } /* done loop over all lines */
  myfile.close();
  
  return ;

} /* done computing line powers */


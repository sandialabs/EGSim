#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "math.h"
#include <string.h>
#include <vector>
#include <complex>

#include "egsim.hpp"

using namespace std;

#define LOCSMALL 1.e-20         
#define EPS      1.0e-9

void getPowers(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Y) ;

void usage() {

  cout<<"Command-line options:"<<endl;
  cout<<"  -h: this message"<<endl;
  cout<<"  -d <gridfile>: grid file"<<endl;
  cout<<"  -e <eps>: solution tolerance"<<endl;
  cout<<"  -a <full/sparse>: save admittance matrix (full or sparse format)"<<endl;
  cout<<"  -p: save line powers"<<endl;
  // cout<<"  -s <solm>:     method for steady state solution"<<endl;
  // cout<<"                 -> GaussSeidel"<<endl;
  // cout<<"                 -> Newton"<<endl;
  // return ;
}

int main(int argc, char *argv[])
{
  Array1D< egsim::bus  >  bdata;
  Array1D< egsim::line >  ldata;

  /* Default settings */
  double eps         = EPS  ;
  bool   saveAdm     = false ;
  bool   saveLinePow = false ;
  string confile("sixbus.dat") ;
  string solm("Newton") ;
  string admout;

  /* Get command-line arguments */
  int argi = 1 ;
  while ( argi < argc )
  {
    if ( strcmp(argv[argi], "-h") == 0 ) {
      usage() ;
      return ( 0 ) ;
    }
    else if ( strcmp(argv[argi], "-d") == 0 ) {
      argi++;
      confile.assign(argv[argi++]);
    }
    // else if ( strcmp(argv[argi], "-s") == 0 )
    // {
    //   argi++;
    //   solm.assign(argv[argi++]);
    // }
    else if ( strcmp(argv[argi], "-e") == 0 ) {
      argi++;
      eps = atof(argv[argi++]);
    }
    else if ( strcmp(argv[argi], "-a") == 0 ) {
      argi++;
      saveAdm = true;
      admout.assign(argv[argi++]);
    }
    else if ( strcmp(argv[argi], "-p") == 0 ) {
      argi++;
      saveLinePow = true;
    }
    else {
      std::cout<<"egsim::main() ERROR: unknown flag "<<argv[argi]<<std::endl;
      exit(1);
    }
    /* done with if over command line input */
  } 

  /* Echo arguments */
  cout<<"1. Reading from :             "<<confile<<endl<<flush ;
  //cout<<"2. Solution method :          "<<solm   <<endl<<flush ;
  cout<<"2. Tolerance :                "<<eps    <<endl<<flush ;
  cout<<"3. Save admittance matrix ? : "<<saveAdm<<endl<<flush ;
  cout<<"---------------------------------------"<<endl<<flush ;

  /* Input data and echo to file */
  egsim::datainputIEEE( confile.c_str(), bdata, ldata) ;

  //setSchInj( bdata ) ;

#ifdef DEBUG
  for ( int i = 0; i<bdata.XSize(); i++ )
    cout<<"Bus "<<i+1<<" "<<bdata(i).get_vri()<<endl;
#endif

  /* Build Admittance Matrix */
  Array2D< complex<double> > Y ;
  egsim::buildAdmMat(bdata, ldata, Y) ;
  if ( saveAdm ) egsim::saveAdmMat(Y,admout);

  /* Newton solution solution */
  egsim::newton_solution(bdata, Y, eps) ;
 
  /* Get powers at nodes */
  //getPowers(bdata, Y) ;

  /* Output solution */
  egsim::saveStaticSol(bdata) ;
  if ( saveLinePow ) getLinePow(bdata, ldata);

  return ( 0 ) ;

} /* end of main */



void getPowers(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Y) {
 
  int nbus = bdata.XSize() ;

  Array1D< complex<double> > VYV(nbus);
  for ( int i = 0; i < nbus ; i++ ) {

      VYV(i) = complex<double> (0.0,0.0);
      for ( int j = 0; j < nbus ; j++ ) 
          VYV(i) += Y(i,j)*bdata(j).get_vri();
      VYV(i) = (bdata(i).get_vri())*conj(VYV(i)) ;
  }

  for ( int i = 0; i < nbus ; i++ ) {

      if ( bdata(i).get_type() != PQld ) 
          bdata(i).set_pgen(VYV(i).real());
      else
          bdata(i).set_pgen(VYV(i).real()+bdata(i).get_pload());

      if ( bdata(i).get_type() != PQld )
          bdata(i).set_qgen(VYV(i).imag());
      else
          bdata(i).set_qgen(VYV(i).imag()+bdata(i).get_qload());

  }  

  return ;

}


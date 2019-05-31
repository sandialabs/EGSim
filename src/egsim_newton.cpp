#include "egsim.hpp"

#define NITERMAX 10
//#define VERBOSE

void egsim::newton_solution(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, double eps) {


  int nbus = (int) bdata.XSize() ;

  /* Set dgesv parameters */
  /* lapack variables */
  int info ;            // return flag
  int ny   = 2 * nbus;  // system size
  int nrhs = 1;         // no. of unknown sets
  Array1D<int> ipiv(ny) ;

  Array1D<double> y(ny,0.0);
  Array1D<double> g(ny,0.0);
  Array2D<double> J(ny,ny,0.0);

  /* Initial condition */
  for ( int i = 0; i < nbus ; i++ ) y(i     ) = bdata(i).get_va() ;
  for ( int i = 0; i < nbus ; i++ ) y(i+nbus) = bdata(i).get_vm() ;

  /* Iteration loop */
  cout << "Beginning Newton iteration"<<endl;
  for ( int k = 0 ; k < NITERMAX ; k++ ) {

      calc_g (bdata, Yadm, y, g) ; /* Compute rhs      */
      calc_Gy(bdata, Yadm, y, J) ; /* Compute Jacobian */
  
#ifdef VERBOSE
      /* Display g and Jacobian matrix */
      for ( int i = 0; i < 2*nbus ; i++ ) 
          cout<<i+1<<" "<<g(i)<<endl;
        cout<<endl;
      for ( int j = 0; j < 2*nbus ; j++ ) 
          for ( int i = 0; i < 2*nbus ; i++ ) 
              if ( abs(J(i,j)) > 1.e-10 )
                  cout<<"("<<i+1<<","<<j+1<<"): "<<J(i,j)<<endl;

#endif

      dgesv_( &ny, &nrhs, J.GetArrayPointer(), &ny, ipiv.GetArrayPointer(), 
              g.GetArrayPointer(), &ny, &info ) ;

      for ( int i = 0; i < ny ; i++ ) y(i) -= g(i) ;

#ifdef VERBOSE
      /* Display increment and current solution */
      cout << "Current increment" << endl;
      for ( int i = 0; i < ny ; i++ ) 
          cout<<"("<<i+1<<"): "<<g(i+nbus)<<" "<<g(i)<<endl;
      cout<<endl;
      cout << "Current solution" << endl;
      for ( int i = 0; i < nbus ; i++ ) 
          cout<<"("<<i+1<<"): "<<y(i+nbus)<<" "<<y(i)<<endl;
#endif

      /* Compute absolute error */
      double dmax = 0.0 ;
      for ( int i = 0; i < ny ; i++ )
          if ( dmax < fabs(g(i)) ) dmax = fabs(g(i)) ;
      cout << " Iteration : "<<k+1<<", error: "<<dmax<<endl;
      if ( dmax < eps ) break ;

  }
  cout<<"---------------------------------------"<<endl<<flush ;

  /* Store voltages in the bus containers */
  for ( int i = 0; i < nbus ; i++ ) {

      bdata(i).set_va( y(i)      );
      bdata(i).set_vm( y(i+nbus) );
      bdata(i).resetV() ;

  }

  return ;

}

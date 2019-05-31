#ifndef EGS_H_SEEN
#define EGS_H_SEEN

#include "copyright.hpp"
#include "array1D.hpp"
#include "array2D.hpp"
#include "egsim_pstruct.hpp"
#include "egsim_bus.hpp"
#include "egsim_line.hpp"

#define PI 3.1415926535897932384

namespace egsim {

	/**
	* \brief Parse transmission model in IEEE format
	*/
	void datainputIEEE(const char *filename, Array1D<bus> &bdata, Array1D<line> &ldata) ;

	/**
	* \brief Compute admittance matrix, full storage
	*/
	void buildAdmMat(Array1D<bus> &bdata, Array1D<line> &ldata, Array2D< complex<double> > &Y) ;

	/**
	* \brief Output admittance matrix to file, 'full'/'sparse' formats
	*/
	void saveAdmMat(Array2D< complex<double> > &Y, string admout) ;

	/**
	* \brief Save static solution (bus voltages)
	*/
	void saveStaticSol(Array1D< bus > &bdata) ;

	/**
	* \brief Compute and save line powers 
	*/
	void getLinePow(Array1D< bus > &bdata, Array1D< line > &ldata) ;

	/**
	* \brief Solution to the non-linear system via Newton iterations
	*/
	void newton_solution(Array1D< bus > &bdata, Array2D< complex<double> > &Y, double eps);

	/**
	* \brief Overall rhs
	*/
	void calc_g     (Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array1D<double> &g) ;
	/**
	* \brief Contributions to rhs from transmission lines
	*/
	void calc_g_line(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array1D<double> &g) ;
	/**
	* \brief Contributions to rhs from PQ loads
	*/
	void calc_g_PQ  (Array1D< egsim::bus > &bdata, Array1D<double> &y, Array1D<double> &g) ;
	/**
	* \brief Contributions to rhs from PV generators
	*/
	void calc_g_PV  (Array1D< egsim::bus > &bdata, Array1D<double> &g) ;
	/**
	* \brief Contributions to rhs from swing node
	*/
	void calc_g_SW  (Array1D< egsim::bus > &bdata, Array1D<double> &g) ;

	/**
	* \brief Overall Jacobian
	*/
	void calc_Gy     (Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array2D<double> &J) ;
	/**
	* \brief Contributions to Jacobian from transmission lines
	*/
	void calc_Gy_line(Array1D< egsim::bus > &bdata, Array2D< complex<double> > &Yadm, Array1D<double> &y, Array2D<double> &J) ;
	/**
	* \brief Contributions to Jacobian from PQ loads
	*/
	void calc_Gy_PQ  (Array1D< egsim::bus > &bdata, Array1D<double> &y, Array2D<double> &J) ;
	/**
	* \brief Contributions to Jacobian from PV generators
	*/
	void calc_Gy_PV  (Array1D< egsim::bus > &bdata, Array2D<double> &J) ;
	/**
	* \brief Contributions to Jacobian from swing node
	*/
	void calc_Gy_SW  (Array1D< egsim::bus > &bdata, Array2D<double> &J) ;


}

/* LAPACK libraries */
extern "C" {
  void dgesv_( int *n, int *nrhs, double *a, int *lda, int *ipiv,
               double *b, int *ldb, int *info ) ;
  void dgeev_ ( char *, char *, int*, double *, int *,
                double *, double *, double *, int *, double *,
                int *, double *, int *, int * );
  void dgemm_ (char *, char *,int *,int *,int *,double *, double *,
               int *,double *,int *,double *,double *,int *) ;
}

#endif



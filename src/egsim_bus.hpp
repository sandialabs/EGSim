/*! \file bus.h
    \brief Header file for bus class
*/ 
#ifndef BUS_H_SEEN
#define BUS_H_SEEN

#include <iostream>
#include <complex>
#include <assert.h>
#include "egsim_pstruct.hpp"

using namespace std ;

namespace egsim {

#define SW   1
#define PVgn 2
#define PQld 3

class bus {

 public:
  bus()  {

    bno_ = -1 ;
    pqG_ = complex<double> (0.0,0.0) ;
    pqL_ = complex<double> (0.0,0.0) ;
    GBs_ = complex<double> (0.0,0.0) ;
    bt_  = -1 ;
    qMx_ = -1.e99 ;
    qMn_ = -1.e99 ;

    pqNetInj_ = complex<double> (0.0,0.0) ;
    pqSchInj_ = complex<double> (0.0,0.0) ;
    netPow_   = complex<double> (0.0,0.0) ;

    lVar_ = rVar_ = hVar_ = false ;

  } 

  ~bus() {} 

  /* Setup methods */
  void setup(const int bno, const egsim::pstruct vma, 
	     const complex<double> pqgen, const complex<double> pqload, 
	     const complex<double> GBshunt, 
	     const int btype, const double qmax, const double qmin)
  {
    bno_ = bno    ;
    vMA_ = vma    ;
    pqG_ = pqgen  ;
    pqL_ = pqload ;
    GBs_ = GBshunt;
    bt_  = btype  ;
    qMx_ = qmax   ;
    qMn_ = qmin   ;

    shunt_ = false ;
    vRI_  = polar(vMA_.m,vMA_.a) ;
    lVar_ = hVar_ = false ;
    rVar_ = true ;

    return ;
  }

  /**
   * \brief Set complex voltage value to the one stored in polar coordinates
   */
  void resetV() { vRI_ = p2c(vMA_); return ; }
  /**
   * \brief Set voltage magnitude
   */
  void set_vm( const double vmag ) {  vMA_.m = vmag ; return ; }
  /**
   * \brief Set voltage angle
   */
  void set_va( const double vang ) {  vMA_.a = vang ; return ; }
  /**
   * \brief Set complex voltage, both polar and complex representations
   */
  void set_vma ( const egsim::pstruct vma ) { vMA_  = vma ; vRI_ = egsim::p2c(vMA_); return ; }
  /**
   * \brief Set complex voltage, both complex and polar representations
   */
  void set_vri ( const complex<double> vri ) { vRI_  = vri ; vMA_ = egsim::c2p(vRI_); return ; }
  /**
   * \brief Set active power for a generator bus
   */
  void set_pgen (const double pgen ) { pqG_.real(pgen) ; return ; }
  /**
   * \brief Set reactive power for a generator bus
   */
  void set_qgen (const double qgen ) { pqG_.imag(qgen) ; return ; }
  /**
   * \brief Set active power for a load bus
   */
  void set_pload(const double pload) { pqL_.real(pload) ; return ; }
  /**
   * \brief Set reactive power for a load bus
   */
  void set_qload(const double qload) { pqL_.imag(qload) ; return ; }
  /**
   * \brief Set apparent power for a generator bus
   */
  void set_gen (const complex<double> pqG ) { pqG_ = pqG ; return ; }
  /**
   * \brief Set apparent power for a load bus
   */
  void set_load(const complex<double> pqL ) { pqL_ = pqL ; return ; }
  /**
   * \brief Set max/min reactive powers
   */
  void set_qmxmn (const double qmax, const double qmin) { 
    qMx_ = qmax   ;
    qMn_ = qmin   ;
    return ; 
  }
  /**
   * \brief Set reference active powers
   */
  void set_pg(const double pg ) { Pg0_ = pg ; return ; }
  /**
   * \brief Set shunt status
   */
  void set_shAdm(const bool s ) {  shunt_ = s ; return ; }
  /**
   * \brief Retrieve shunt status
   */
  bool get_shAdm() {  return(shunt_) ; }
  /**
   * \brief Set max voltage magnitude
   */
  void set_Vmx(const double v ) {  vMx_ = v ; return ; }
  /**
   * \brief Set min voltage magnitude
   */
  void set_Vmn(const double v ) {  vMn_ = v ; return ; }
  void set_z  (const double z ) {  z_   = z ; return ; }
  /**
   * \brief Retrieve max voltage magnitude
   */
  double get_Vmx() {  return (vMx_) ; }
  /**
   * \brief Retrieve min voltage magnitude
   */
  double get_Vmn() {  return (vMn_) ; }

  /**
   * \brief Set net apparent power injected at the bus
   */
  void set_netInj(const complex<double> pqNetInj ) { pqNetInj_ = pqNetInj ; return ; }
  /**
   * \brief Set apparent power scheduled to be injected at the bus
   */
  void set_schInj(const complex<double> pqSchInj ) { pqSchInj_ = pqSchInj ; return ; }
  /**
   * \brief Set net apparent power at the bus
   */
  void set_netPow(const complex<double> netPow   ) { netPow_   = netPow   ; return ; }

  void set_pnetInj(const double pNetInj ) { pqNetInj_.real(pNetInj) ; return ; }
  void set_pschInj(const double pSchInj ) { pqSchInj_.real(pSchInj) ; return ; }
  void set_pnetPow(const double pnetPow ) { netPow_.real(pnetPow)   ; return ; }

  void set_qnetInj(const double qNetInj ) { pqNetInj_.imag(qNetInj) ; return ; }
  void set_qschInj(const double qSchInj ) { pqSchInj_.imag(qSchInj) ; return ; }
  void set_qnetPow(const double qnetPow ) { netPow_.imag(qnetPow)   ; return ; }

  void set_lVar(const bool vlim) { lVar_ = vlim ; if ( vlim ) { rVar_ = hVar_ = false; } return ; }
  void set_rVar(const bool vlim) { rVar_ = vlim ; if ( vlim ) { lVar_ = hVar_ = false; } return ; }
  void set_hVar(const bool vlim) { hVar_ = vlim ; if ( vlim ) { lVar_ = rVar_ = false; } return ; }

  /**
   * \brief Return bus type: slack (1)/PV generator (2)/PQ load (3)
   */
  int get_type( ) { return bt_ ; }
  /**
   * \brief Return voltage magnitude
   */
  double get_vm( ) { return vMA_.m ; }
  /**
   * \brief Return voltage angle
   */
  double get_va( ) { return vMA_.a ; }
  /**
   * \brief Return voltage (polar representation)
   */
  egsim::pstruct get_vma ( ) { return vMA_ ; }
  /**
   * \brief Return voltage (complex representation)
   */
  complex<double> get_vri ( ) { return vRI_ ; }

  double get_pgen ( ) { return pqG_.real() ; }
  double get_qgen ( ) { return pqG_.imag() ; }
  double get_pload( ) { return pqL_.real() ; }
  double get_qload( ) { return pqL_.imag() ; }

  complex<double> get_gen ()  { return pqG_ ; }
  complex<double> get_load()  { return pqL_ ; }
  complex<double> get_shunt() { return GBs_ ; }

  /**
   * \brief Return maximum reactive power
   */
  double get_qmx () { return qMx_ ;} 
  /**
   * \brief Return minimum reactive power
   */
  double get_qmn () { return qMn_ ;} 

  complex<double> get_netInj( ) { return pqNetInj_ ; }
  complex<double> get_schInj( ) { return pqSchInj_ ; }
  complex<double> get_netPow( ) { return netPow_   ; }

  double get_pnetInj( ) { return pqNetInj_.real(); }
  double get_pschInj( ) { return pqSchInj_.real(); }
  double get_pnetPow( ) { return netPow_.real()  ; }
        
  double get_qnetInj( ) { return pqNetInj_.imag() ; }
  double get_qschInj( ) { return pqSchInj_.imag() ; }
  double get_qnetPow( ) { return netPow_.imag()   ; }

  bool get_lVar() { return lVar_ ; }
  bool get_rVar() { return rVar_ ; }
  bool get_hVar() { return hVar_ ; }

  /**
   * \brief Set connectivity flag
   */
  void set_u (int u) { u_ = u ; return ;}
  /**
   * \brief Return connectivity flag
   */
  int  get_u (     ) { return ( u_ );   }

  private:

  int bno_ ; /**< bus number */

  double pbase_ ; /**< base active power */
  double vbase_ ; /**< base voltage */

  egsim::pstruct  vMA_ ; /**< voltage, polar coordinates  */ 
  complex<double> vRI_ ; /**< voltage, complex representation */

  complex<double> pqG_ ; /**< active and reactive generator power         */
  complex<double> pqL_ ; /**< active and reactive load power              */
  complex<double> GBs_ ; /**< shunt impedance (resistance and reactance)  */ 

  complex<double> pqNetInj_, pqSchInj_, netPow_ ;

  int    bt_  ; /**< 1-swing bus, 2-PV bus (generator), 3-PQ (load) */ 
  double qMx_ ; /**< Max reactive power  */
  double qMn_ ; /**< Min reactive power  */
  double vMx_ ; /**< Max voltage */
  double vMn_ ; /**< Min voltage */

  double Pg0_ ;
  double gamma_ ;
  int z_ ;
  int u_ ; /**< Connectivity flag: 1-connected, 0-not connected */ 

  bool lVar_, rVar_, hVar_ ;

  bool shunt_ ;

  int An_ ; /**< area number   */
  int Rn_ ; /**< region number */

  double srtng_ ; /**< power rating   */
  double vrtng_ ; /**< voltage rating */

} ;

}
#endif

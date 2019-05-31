/*! \file line.h
    \brief Header file for line class
*/ 
#ifndef EGS_LINE_H_SEEN
#define EGS_LINE_H_SEEN

#include <complex>
#include "egsim_bus.hpp"

namespace egsim {

class line {

 public:

  line() { state_ = 1; }
  ~line() {}

  void setup(const int frombus,const int tobus, const complex<double> impd,
             const double chrg, const double tap, const double phase_shift_ang)
  {
    fB_  = frombus ;
    tB_  = tobus   ;
    impd_= impd    ;
    lC_  = chrg    ;
    tR_  = tap     ;
    psa_ = phase_shift_ang ;
    if ( tR_ < 1.0e-6 ) {
      tR_ = 1.0 ;
      psa_ = 0.0;
    }
    return ;
  }

  void setup_ln(const int frombus,const int tobus,const double srtng,const double vrtng,
                const double frtng,const double llen, const complex<double> impd,
                const double chrg,const double iMax,const double pMax,const double pqMax)
  {
    fB_  = frombus ;
    tB_  = tobus   ;
    impd_= impd    ;
    lC_  = chrg    ;

    srtng_ = srtng ;
    vrtng_ = vrtng ;
    frtng_ = frtng ;

    llen_  = llen ;

    iMax_  = iMax  ;
    pMax_  = pMax  ;
    pqMax_ = pqMax ;

    return ;
  }

  /**
   * \brief Switches line state to off (0)
   */
  void set_off() { state_ = 0 ; return ; } 
  /**
   * \brief Switches line state to on (1)
   */
  void set_on()  { state_ = 1 ; return ; } 

  /* get methods */
  /**
   * \brief Returns state of the line: 1 (on) or 0 (off)
   */
  int get_state() { return state_ ; }
  /**
   * \brief Returns 'from' bus id
   */
  int get_frbus() { return fB_    ; }
  /**
   * \brief Returns 'to' bus id
   */
  int get_tobus() { return tB_    ; }
  /**
   * \brief Returns line impedance
   */
  complex<double> get_impd() { return impd_ ; }
  /**
   * \brief Returns line charging
   */
  double get_lchrg() { return lC_ ; }
  /**
   * \brief Returns tap ratio
   */
  double get_tr() { return tR_ ; }
  /**
   * \brief Returns phase shift angle
   */
  double get_psa() { return psa_ ; }

 private:

  int    state_ ; /**< line status - on (1)/off (0)  */
  int    fB_    ; /**< from bus id's */
  int    tB_    ; /**< to bus id's */
  
  complex<double> impd_ ; /**< line impedance=resistance + j reactance */

  double lC_  ; /**< line charging */ 
  double tR_  ; /**< tap ratio     */ 
  double psa_ ; /**< phase shift angle */
 
  double srtng_ ; /**< power rating     */
  double vrtng_ ; /**< voltage rating   */
  double frtng_ ; /**< frequency rating */

  double llen_ ; /**< line length */

  double iMax_  ; /**< current limit        */
  double pMax_  ; /**< active power limit   */
  double pqMax_ ; /**< apparent power limit */

} ;

}
#endif

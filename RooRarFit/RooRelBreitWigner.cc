/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: RooRelBreitWigner.cc,v 1.10 2014/09/24 17:09:29 fwilson Exp $
 * Authors:                                                                  
 *                                             
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This is an implementation of Relativistic Breit Wigner for Spin 0,1,2 particles
// with Blatt-Weisskopf form factors and barrier functions for RooRarFit. 
// The meson radius is set to 3.1 GeV^-1 which is similar to the EvtGen number 
// but may not be appropriate for every meson. See Phys Rev D 72, 052002 (2005).
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is an implementation of Relativistic Breit Wigner for Spin 0,1,2 particles
// with Blatt-Weisskopf form factors and barrier functions for RooRarFit. 
// The meson radius is set to 3.1 GeV^-1 which is similar to the EvtGen number 
// but may not be appropriate for every meson. See Phys Rev D 72, 052002 (2005).
// END_HTML
//

#include "rarVersion.hh"

#include "Riostream.h"

#include "RooRelBreitWigner.hh"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "TMath.h"

#include <complex>
typedef std::complex<Double_t> dcmplx;

ClassImp(RooRelBreitWigner)

//------------------------------------------------------------------
RooRelBreitWigner::RooRelBreitWigner(const char *name, const char *title,
				     RooAbsReal& _x, RooAbsReal& _mean,
				     RooAbsReal& _width, RooAbsReal& _radius, 
				     RooAbsReal& _mass_a, RooAbsReal& _mass_b, 
				     RooAbsReal& _spin) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  mean("mean","Mean",this,_mean),
  width("width","Width",this,_width),
  radius("radius","Radius",this,_radius),  // meson radius e.g. 3.1 GeV^-1
  mass_a("mass_a","Mass of first daughter",this,_mass_a), // e.g. pion from K*
  mass_b("mass_b","Mass of second daugher",this,_mass_b), // e.g. kaon from K*
  spin("spin","Spin",this,_spin)
{
}

//------------------------------------------------------------------
RooRelBreitWigner::RooRelBreitWigner(const RooRelBreitWigner& other, 
				     const char* name) : 
  RooAbsPdf(other,name), 
  x("x",this,other.x), 
  mean("mean",this,other.mean),
  width("width",this,other.width),
  radius("radius",this,other.radius),
  mass_a("mass_a",this,other.mass_a),
  mass_b("mass_b",this,other.mass_b),
  spin("spin",this,other.spin)
{
}

//------------------------------------------------------------------
Double_t RooRelBreitWigner::evaluate() const
{
  Double_t temp = mean*getWidth();
  // numerator
  dcmplx T(sqrt(temp),0.0);
  // denominator
  dcmplx denom(mean*mean-x*x,-1*temp);
  T = T /denom;     // Transition probability
  return(std::abs(T) * std::abs(T)); // Intensity (arbitrary scale)
}

//------------------------------------------------------------------
Double_t RooRelBreitWigner::getWidth() const
{
  Double_t q  = getQ(x);
  Double_t q0 = getQ(mean);
  Double_t result(0.0);

  if (q>0 && q0>0 && x>0 && mean>0) {
    result = width * getQterm(q,q0) * (mean/x) 
      * (getFF(q) / getFF(q0));
  }
  return (result);
}

//------------------------------------------------------------------
Double_t RooRelBreitWigner::getQterm(Double_t q, Double_t q0) const
{
  Double_t result = pow(q/q0,2*spin+1);
  return (result); // result = (q/q0)^(2*spin+1);
}

//------------------------------------------------------------------
Double_t RooRelBreitWigner::getFF(Double_t q) const
{
  
  Double_t z = q * radius;
  Double_t result(1.0);

  Double_t zsq = z*z;
  if (spin==1) {
    result = 1.0/(1+zsq);
  } else if (spin==2) {
    result = 1.0/(9 + 3*zsq + zsq*zsq);
  }

  return (result); //square of the Blatt-Weisskopf form factor
}

//-----------------------------------------------------
Double_t RooRelBreitWigner::getQ(Double_t mass) const
{

  if (mass < (mass_b+mass_a)) {return(0);}

  const Double_t mDaugSumSq  = (mass_b+mass_a)*(mass_b+mass_a);
  const Double_t mDaugDiffSq = (mass_b-mass_a)*(mass_b-mass_a);

  Double_t q = sqrt((mass*mass-mDaugSumSq)*(mass*mass-mDaugDiffSq))/(2*mass);
  return(q);
}

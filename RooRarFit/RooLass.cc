/******************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: RooLass.cc,v 1.6 2014/09/14 17:33:47 fwilson Exp $
 * Authors:                                                       *
 * Fergus Wilson, RAL, fwilson@slac.stanford.edu                  *
 *                                                                *
 * Copyright (C) 2005-2012, RAL                                   *
 *                                                                *
 ******************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This is an implentation of the LASS function for RooRarFit
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is an implentation of the LASS function for RooRarFit
// END_HTML
//

#include "Riostream.h"

#include "RooLass.hh"
#include "RooAbsReal.h"
#include "RooRealVar.h"

#include "TMath.h"

#include <complex>

typedef std::complex<Double_t> dcmplx;

using std::cout;
using std::cerr;
using std::endl;

ClassImp(RooLass)

RooLass::RooLass(const char * name, const char * title, 
		 RooAbsReal & fitVariable, RooAbsReal & theMean,
		 RooAbsReal & theWidth, RooAbsReal & effectiveRange,
		 RooAbsReal & scatteringLength, RooAbsReal & turnOffValue)
  : RooAbsPdf(name,title),
    x("x","Dependent",this,fitVariable),
    mean("mean","Mean",this,theMean),
    width("width","Width",this,theWidth),
    effRange("effRange","Effective Range",this,effectiveRange),
    scatLen("scatLen","Scattering Length",this,scatteringLength),
    turnOffVal("turnOffVal","Turn Off Value",this,turnOffValue)
{
}

RooLass::RooLass(const RooLass & other, const char * name)
  : RooAbsPdf(other,name),
    x("x",this,other.x),
    mean("mean",this,other.mean),
    width("width",this,other.width),
    effRange("effRange",this,other.effRange),
    scatLen("scatLen",this,other.scatLen),
    turnOffVal("turnOffVal",this,other.turnOffVal)
{
}

//------------------------------------------------------------
Double_t RooLass::getQ(Double_t mass) const
{

  const Double_t m_Kaon = 0.493677;
  const Double_t m_Pion = 0.13957018;

  if (mass < (m_Kaon+m_Pion)) {return(0);}

  const Double_t mDaugSumSq  = (m_Kaon+m_Pion)*(m_Kaon+m_Pion);
  const Double_t mDaugDiffSq = (m_Kaon-m_Pion)*(m_Kaon-m_Pion);

  Double_t q  = sqrt((mass*mass-mDaugSumSq)*(mass*mass-mDaugDiffSq))/(2*mass);
  return(q);
}

//------------------------------------------------------------
Double_t RooLass::evaluate() const
{
  //return (kmatrix());
  return (smatrix());
}

//----------------------------------------------------
Double_t RooLass::kmatrix() const
{
  Double_t mass = x;

  Double_t q  = getQ(mass);
  if (q==0) {return(0);}

  Double_t q0 = getQ(mean);
  if (q0==0) {return(0);}

  Double_t rho  = 2*q/mass;
  Double_t rho0 = 2*q0/mean;

  Double_t g0sqr = mean*mass*width/rho0;
  // K matrix for BW resonance
  Double_t Khat = g0sqr/(mean*mean - mass*mass);

  if (mass <= turnOffVal) {
    // K matrix for S-wave background + BW resonance
    Khat += ((scatLen*mass)/(2+scatLen*effRange*q*q));
  }

  dcmplx T(Khat,0.0);
  dcmplx denom(1.0,-Khat*rho);
  T = T / denom; // Transition probability

  return (std::abs(T) * std::abs(T)); // Amplitude (arbitrary scale);
}

//----------------------------------------------------
Double_t RooLass::smatrix() const
{
  Double_t mass = x;

  const Double_t q  = getQ(mass);
  if (q==0) {return(0);}

  const Double_t qR = getQ(mean);
  if (qR==0) {return(0);}

  // Breit-Wigner resonace
  dcmplx T2(width*mean*mean/qR,0);
  dcmplx bw_denom(mean*mean - mass*mass, -mean*width*q*mean/(mass*qR));
  
  T2 = T2 / bw_denom;

  dcmplx T(T2);

  if (mass < turnOffVal) {
    // effective range parameterization
    Double_t cot_db = 1.0/ (scatLen*q) + (effRange*q)*0.5;
    Double_t db = TMath::ATan(1.0/cot_db);

    // relative phase exp^(2 i db)
    dcmplx phase(cos(2*db),sin(2*db));

    // S-wave
    dcmplx T1(mass, 0);
    dcmplx s_denom(q*cot_db, -q);

    T = (T1/s_denom) + (phase * T); 
  }

  return (std::abs(T) * std::abs(T)); // Amplitude (arbitrary scale);
}

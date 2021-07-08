/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: RooThreshold.cc,v 1.2 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 * 
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/

// -- CLASS DESCRIPTION [PDF] --
// This is an implementation of a generic threshold function
//
//
//////////////////////////////////////////////////////
//
// BEGIN_HTML
// This is an implementation of a generic threshold function
//
// END_HTML
//
// Begin_Latex Amplitude : A(m) = ({m-m_{0})^{p} exp(P00 x + P01 x^{2} + ...) End_Latex

#include "Riostream.h"
#include "TMath.h"

#include "RooAbsReal.h"
#include "RooRealVar.h"

//#include "rarVersion.hh"
#include "RooThreshold.hh"

using namespace std;

ClassImp(RooThreshold)

//------------------------------------------------------------------
RooThreshold::RooThreshold(const char *name, const char *title,
			   RooAbsReal& _x, RooAbsReal& _m0,
			   RooAbsReal& _power, const RooArgList& _coeffs) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  m0("mean","Mean",this,_m0),
  power("power","power",this,_power),
  coeffs("coeffs","coeffs",this)
{
  // copy over coefficeints
  TIterator* coefIter = _coeffs.createIterator() ;
  RooAbsArg* coef(0) ;
  while ((coef = (RooAbsArg*)coefIter->Next())) {
    if (!dynamic_cast<RooAbsReal*>(coef)) {
      cout << "RooTheshold::ctor(" << GetName() << ") ERROR: coefficient " << coef->GetName() 
	   << " is not of type RooAbsReal" << endl ;
      assert(0) ;
    }
    coeffs.add(*coef) ;
  }
  delete coefIter ;
}

//------------------------------------------------------------------
RooThreshold::RooThreshold(const RooThreshold& other, 
		     const char* name) : 
  RooAbsPdf(other,name), 
  x("x",this,other.x), 
  m0("m0",this,other.m0),
  power("power",this,other.power),
  coeffs("coeffs",this, other.coeffs)
{
}

//------------------------------------------------------------------
Double_t RooThreshold::evaluate() const
{
  // calculate Threshold

  if (x<m0) {return(0);} // below threshold
  Double_t result(0);

  // loop over coefficient
  TIterator* iter = coeffs.createIterator() ;
  RooRealVar* coef(0) ;
  Int_t n(0);
  while ((coef = (RooRealVar*) iter->Next())) {
    n++;
    result += coef->getVal() * TMath::Power(x,n);
  }
  delete iter ;

  result = TMath::Exp(result) * TMath::Power(x-m0,power);

  //  cout << result << " " << x << " " << power << " " << m0 << endl;
  return(result);
}


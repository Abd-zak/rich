/*****************************************************************************
* Project: BaBar detector at the SLAC PEP-II B-factory
* Package: RooRarFit
 *    File: $Id: RooOsipDisc.rdl,v 1.5 2014/09/14 17:33:47 fwilson Exp $
 * Authors:                                                                  *
 *   WV, Wouter Verkerke, UC Santa Barbara, verkerke@slac.stanford.edu       *
 *   DK, David Kirkby,    UC Irvine,         dkirkby@uci.edu                 *
 *                                                                           *
 * Copyright (C) 2005-2012, Regents of the University of California          *
 *                    and Stanford University. All rights reserved.    *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/
#ifndef ROO_OSIPDISC	
#define ROO_OSIPDISC

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooOsipDisc : public RooAbsPdf {
public:
  RooOsipDisc(const char *name, const char *title,
	      RooAbsReal& _x, RooAbsReal& _r, RooAbsReal& _b, RooAbsReal& _a, RooAbsReal& _h, RooAbsReal& _s );
  RooOsipDisc(const RooOsipDisc& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooOsipDisc(*this,newname); }
  inline virtual ~RooOsipDisc() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

  Int_t getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK=kTRUE) const;
  void generateEvent(Int_t code);

protected:

  RooRealProxy x ;
  RooRealProxy r ;
  RooRealProxy bb ;
  RooRealProxy a ;
  RooRealProxy h ;
  RooRealProxy s ;

  Double_t evaluate() const ;

private:

  ClassDef(RooOsipDisc,0) 
};

#endif

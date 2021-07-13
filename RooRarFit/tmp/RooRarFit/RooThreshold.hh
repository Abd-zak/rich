/*****************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: RooRarFit
 *    File: $Id: RooThreshold.rdl,v 1.2 2014/09/14 17:33:48 fwilson Exp $
 * Authors: Fergus Wilson
 * History:
 *
 * Copyright (C) 2005-2012, RAL
 *****************************************************************************/
#ifndef ROO_THRESHOLD
#define ROO_THRESHOLD

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooListProxy.h"

class RooRealVar;
class RooArgList;

class RooThreshold : public RooAbsPdf {

public:
  RooThreshold(const char *name, const char *title,
	       RooAbsReal& _x, RooAbsReal& _m0, 
	       RooAbsReal& _power, const RooArgList& _coeffs);

  RooThreshold(const RooThreshold& other, 
		    const char* name=0) ;
  virtual TObject* clone(const char* newname) const { 
    return new RooThreshold(*this,newname); 
  }
  inline virtual ~RooThreshold() { }

protected:

  RooRealProxy x ;     // 
  RooRealProxy m0 ;  // threshold value
  RooRealProxy power ;  // power
  RooListProxy coeffs ; // 
 
  Double_t evaluate() const ;

private:
 
  ClassDef(RooThreshold,0) // Threshold PDF 
};

#endif

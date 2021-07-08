// The function is a fourthorder Polynomial for 0<x<r and a Gaussian for r<x<1
// The restrictions are that: 1.The derivative be 0 at x=0, 2. That the function be continuous at r, 3. That the first derivative be continuous at r.
// by Osipenkov, Ilya.

#include "rarVersion.hh"

#include "Riostream.h"
#include "TMath.h"

#include "RooOsipDisc.hh"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooRandom.h"
#include "RooMath.h"

ClassImp(RooOsipDisc)

RooOsipDisc::RooOsipDisc(const char *name, const char *title,
			 RooAbsReal& _x, RooAbsReal& _r,
			 RooAbsReal& _b,RooAbsReal& _a,RooAbsReal& _h,RooAbsReal& _s) :
  RooAbsPdf(name,title),
  x("x","Dependent",this,_x),
  r("r","r",this,_r),
  bb("b","b",this,_b),
  a("a","a",this,_a),
  h("h","h",this,_h),
  s("s","s",this,_s)
{
}

RooOsipDisc::RooOsipDisc(const RooOsipDisc& other, const char* name) : 
  RooAbsPdf(other,name), x("x",this,other.x), r("r",this,other.r),
  bb("b",this,other.bb),a("a",this,other.a),h("h",this,other.h),s("s",this,other.s)
{
}


Double_t RooOsipDisc::evaluate() const
{
  Double_t t = -(2*a*r+3*bb*r*r+4*s*r*r*r)*h*h/(1+a*r*r+bb*r*r*r+s*r*r*r*r);
  Double_t d = (1+a*r*r+bb*r*r*r+s*r*r*r*r)*exp(t*t/(2*h*h));
  if (x<r) {
    return 1+a*x*x+bb*x*x*x+s*x*x*x*x;
  }
  if ((x>r)||(x==r)) {
  return d*exp(-0.5*(x-r+t)*(x-r+t)/(h*h));
  }
  return 0;
}

Int_t RooOsipDisc::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,x)) return 1 ;
  return 0 ;
}

Double_t RooOsipDisc::analyticalIntegral(Int_t code, const char* /* rangeName */) const 
{
  assert(code==1) ;

  const Double_t rootPiBy2 = sqrt(TMath::PiOver2());
  Double_t t = -(2*a*r+3*bb*r*r+4*s*r*r*r)*h*h/(1+a*r*r+bb*r*r*r+s*r*r*r*r);
  Double_t d = (1+a*r*r+bb*r*r*r+s*r*r*r*r)*exp(t*t/(2*h*h));
  Double_t xscale = TMath::Sqrt2() * h;

  return  r+a*r*r*r/3+bb*r*r*r*r/4+s*r*r*r*r*r/5+d*rootPiBy2*h*(TMath::Erf((1-r+t)/xscale)-TMath::Erf(t/xscale)); 
}

Int_t RooOsipDisc::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t /*staticInitOK*/) const
{
  if (matchArgs(directVars,generateVars,x)) return 1 ;  
  return 0 ;
}

void RooOsipDisc::generateEvent(Int_t code)
{
  assert(code==1) ;
  Double_t xgen ;
  while(1) {    
    xgen = 0.5; //RooRandom::randomGenerator()->OsipDisc(r,bb,a,h,s);
    if (xgen<x.max() && xgen>x.min()) {
      x = xgen ;
      break;
    }
  }
  return;
}



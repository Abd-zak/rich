// --------------------------------------------------------------
//  A parameterisation of LASS line shape
//
//  See BAD 1282 pages 100-101 for K matrix version
//  See BAD 1342 pages 9-10 for S matrix version
//
//  Effective Range   = 3.32 +/- 0.34 (Gev/c)^-1
//  Scattering length = 2.07 +/- 0.10 (Gev/c)^-1
//  turnOffVal = 1.6 GeV/c^2
//
// --------------------------------------------------------------

#ifndef ROO_KMATRIXLASS
#define ROO_KMATRIXLASS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooLass : public RooAbsPdf {
  
public:

  RooLass(const char * name, const char * title, 
	RooAbsReal & fitVariable, RooAbsReal & theMean,
    	RooAbsReal & theWidth, RooAbsReal & effectiveRange,
	RooAbsReal & scatteringLength, RooAbsReal & turnOffValue);

  RooLass(const RooLass & other, const char * name = 0);
  virtual TObject * clone(const char * newname) const {return new RooLass(*this,newname);}
  inline virtual ~RooLass() {}
  
  //Int_t getAnalyticalIntegral(RooArgSet & allVars, RooArgSet & analVars, const char * rangeName = 0) const;
  //Double_t analyticalIntegral(Int_t code, const char * rangeName = 0) const;
  
protected:
  
  RooRealProxy x;
  RooRealProxy mean;
  RooRealProxy width;
  TString      spec;
  RooRealProxy effRange;
  RooRealProxy scatLen;
  RooRealProxy turnOffVal;
  
  Double_t evaluate() const;
  
private:
  Double_t getQ(Double_t mass) const;
  Double_t kmatrix() const; // K-matrix version
  Double_t smatrix() const; // S-matrox version from BAD 1341 p9
  
  ClassDef(RooLass,0)
    
};

#endif

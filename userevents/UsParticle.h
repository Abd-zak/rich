#ifndef UsParticle_h
#define UsParticle_h
/*!
  \class V0
  \brief Reconstructed V0

  \author Yann.Bedfer@cern.ch
*/

#include "TLorentzVector.h"

class UsParticle
{
 public:
  UsParticle(const TLorentzVector &lv0,
	 int ivs, int iET1, int iET2) {
    lv = TLorentzVector(lv0);
    iv = ivs; i1 = iET1; i2 = iET2;
    type = 0;
  }
  UsParticle(const TVector3 &v0,int ivs, int iET1, int iET2,
	     double M2) {  // PDG mass squared
    double E = sqrt(M2+v0.Mag2()); lv = TLorentzVector(v0,E);
    iv = ivs; i1 = iET1; i2 = iET2;
    type = 0;
  }
  TLorentzVector lv;
  int iv, i1, i2, type;
};
#endif

// $Id: MCLambda.h,v 1.3 2020/07/30 10:55:51 ybedfer Exp $

// Class for housing the Lambda acceptance and efficiency determinations

// $Log: MCLambda.h,v $
// Revision 1.3  2020/07/30 10:55:51  ybedfer
// No change.
//
// Revision 1.2  2012/11/05 19:58:13  ybedfer
// Upgrade c++ source.
//
// Revision 1.1  2009/07/03 14:05:11  ybedfer
// Initial revision
//

#ifndef MCLambda_h
#define MCLambda_h

#include "TH1.h"
#include "PaEvent.h"
#include "MCInfo.h"

class MCLambda {
public:
  MCLambda(double yLow, double yUp, double pTCut); // Constructor. Book histos.
  static MCLambda *Ptr() { return address; }  // Return pointer to MCLambda
  void Update();  // Import the MCInfo flags of the current PaEvent's MC tracks and check (and store) MC Q2 and y are w/in range
  // Fill: Retuned value = bit pattern: LS 2 bits = 2^iLaL, if argument (p,pi) corresponds to a genuine MC Lambda w/ C = "iLaL"
  int Fill(PaEvent &e, // Only used to document the error messages
	   PaTrack &tp, PaTrack &tpi, // p and pi reco'd tracks
	   int iLaL,      // Lambda or anti-Lambda
	   bool winMCut);  // Reco passed mass cut
private:
  static MCLambda *address;        // Address of this->object
  MCInfo *MC;
  int *MCTypes; // Array of MCInfo flags of current PaEvent
  int MCDIS;    // Bit pattern: 0x1 = MC Q2>1, 0x2 = y w/in Range 
  //! Kinematical and event cuts
  double yLowCut, yUpCut;     // Event cuts
  double V0pTCut, V0cthCut;   // Kin cuts

  // Histograms for Lambda (Sigma*) acceptance, vs. cos(theta*) and vs ppi
  // Ar = Accepted (Q2,y) and reconstructed
  // r  = reconstructed outside Acceptance
  // Rr = Reconstructible and reconstructed
  TH1D *hc_ArL[2][2], *hc_ArLk[2][2], *hc_ArLi[2][2], *hc_ArLki[2][2];
  TH1D *hc_rL[2][2],  *hc_rLk[2][2],  *hc_rLi[2][2],  *hc_rLki[2][2];
  TH1D *hc_RrL[2][2], *hc_RrLk[2][2], *hc_RrLi[2][2], *hc_RrLki[2][2];
  TH1D *hc_ArS[2], *hc_RrSki[2];
  TH1D *hc_rS[2];
  
};
#endif

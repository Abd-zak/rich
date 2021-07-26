// $Id: CSEventData.h,v 1.5 2021/04/05 12:42:42 ybedfer Exp $

// Based on Luigi's DISEventData.h


#ifndef CSEventData_h
#define CSEventData_h 1

#include <math.h>
#include "Rtypes.h"
#include <vector>

struct CSEventData
{
  CSEventData(): runNo(-1),spillNo(-1),evtNo(-1),
    trigMask(0),primaryPat(0),spillOKPat(0),
    Xp(0),Yp(0),Zp(0),    // pVertex
    E0(0),Q2(0),y(0),
    index(0),prodIndex(0),piThr(0),
    dEK(100), // A large number, that won't pass any exclusivity (what about incl. phi's?)
    nOuts(0), nKs(0), nTrksRIt(0), nTrksRIb(0),
    KThr(0), pThr(0)
  {}
  virtual ~CSEventData();

  Int_t runNo,spillNo,evtNo;
  UInt_t trigMask;
  // ***** "PRIMARY" FLAG:
  // 0x1  : There is p(rimary)Vertex.
  // 0x2  : If muon beam, pV's BMS w/in range (=> implicitly, genuine BMS exists)
  // 0x4  : If muon beam, pV's backpropagation to BMS OK
  // 0x8  : If muon beam, pV has scattered mu.
  // 0x10 : If muon beam, pV's scattered mu abides by PaVertex::isMuPrim
  // 0x20 : If muon beam, pV's scattered mu has chi2/NDF w/in cut
  // 0x40 : W/in target zone
  // 0x80 : W/in target (w/, possibly, a tighter requirement than 0x10).
  // 0x100: If !hadronPhysics, fulfills y requirements
  // 0x200: If !hadronPhysics, fulfills Q2 cut 
  UShort_t primaryPat;
  // ***** BAD SPILL FLAG:
  // 0x1: Hodo trigger stability
  // 0x2: Calo
  // 0x4: RICH
  UShort_t spillOKPat;
  Float_t  Xp, Yp, Zp; // pVertex
  Float_t  E0, Q2, y;
  Float_t  index, prodIndex; // Index-1: @ PHAST time, @ mass prod. time
  Float_t  piThr; // Cherenkov threshold for pi. KThr, pTHr: cf. "Expand" infra.

  // ***** EXCLUSIVITY phi
  Float_t  dEK; // Exclusivity (EMiss) evaluated w/ K mass

  // ***** ESTIMATORS of EVENT MULTIPLICITY (e.g. in RICH)
  Int_t    nOuts; // # of charged tracks in BPV
  Int_t    nKs;   // # of charged K's in BPV
  Int_t    nTrksRIt, nTrksRIb;  // #tracks in top/bottom of RICH

  void Reset()
  {
    runNo = 0; spillNo = 0; evtNo = 0; 
    trigMask = 0; primaryPat = 0; spillOKPat = 0;
    Xp = 0; Yp = 0; Zp = 0; // pVertex
    E0 = 0; Q2 = 0; y = 0;
    index = 0; prodIndex = 0; piThr = 0;
    dEK = 100;  // A large number, that won't pass any exclusivity cut (what about incl. phi's?)
    nOuts = 0; nKs = 0; nTrksRIt=nTrksRIb = 0;
    KThr=pThr = 0;
  }

  void Expand()
  {
    KThr = piThr*3.53712/* M_K/M_pi */, pThr = piThr*6.72258/* M_p/M_pi */;
  }

  //float P; int q; // Added at some point. No longer remember why...
  float KThr, pThr;

  ClassDef(CSEventData,3); // Must be the last item before the closing '};'
};



struct CSHadronData {

CSHadronData(): Px(0), Py(0), Pz(0), qP(0),XX0(0),ZFirst(0),ZLast(0),chi2(0),ECAL(0),phiR(0),thC(0),
    MCpid(-1),
    hasR(false),XR(0),YR(0),tgXR(0),tgYR(0),
    P(0),q(0),thR(0)
  {
    for (int i = 0; i<6; ++i) LH[i] = 0;
    for (int i = 0; i<3; ++i) dLHdI[i] = 0;
  }
  virtual ~CSHadronData();

  double   Px, Py, Pz; // Momentum @ pVertex.
  Float_t  qP;  // Momentum at, or close to, RICH
  Float_t  XX0;
  Float_t  ZFirst, ZLast, chi2;
  Float_t  ECAL;
  Float_t  phiR; // Azimuth at, or close to, RICH

  Float_t  LH[6];    // pi,K,p,e,mu,back. Note: small LHs (<1e-6) are set =0.
  Float_t  dLHdI[3]; // dLH/dI: pi,K,p
  Float_t  thC;      // Cherenkov angle (we chose max. LH angle)

  Short_t  MCpid;

  Bool_t   hasR;        // Has RICH: i.e. a RICH block w/ a finite background LH
  Float_t  XR, YR;
  Float_t  tgXR, tgYR;

  void Reset()
  {
    Px=Py=Pz = 0;
    qP = 0;
    XX0 = 0;
    ZFirst = 0; ZLast = 0; chi2 = 0;
    ECAL = 0; phiR = 0;
    for (int i = 0; i<6; ++i) LH[i] = 0;
    for (int i = 0; i<3; ++i) dLHdI[i] = 0;
    thC = 0;
    MCpid = -1;
    hasR = false; XR=YR = 0; tgXR=tgYR = 0;
    P = 0; q = 0;    
    thR = 0;
  }

  void Expand()
  {
    P = fabs(qP); q = qP >0 ? 1 : -1;
    thR = acos(1/sqrt(1.+ tgXR*tgXR + tgYR*tgYR));
  }

  float P; int q;
  float thR;

  ClassDef(CSHadronData,3); // Must be the last item before the closing '};'
};


struct CSResonanceData{

CSResonanceData(): phiPat(0),K0Pat(0),LambdaPat(0),
    h1(-1),h2(-1),m(0),
    Xs(0),Ys(0),Zs(0),    // sVertex
    D(0),dD(0),cth(0),pT(0),alpha(0),chi2(0)
  {}
  virtual ~CSResonanceData();
  // ***** EXCLUSIVE phi SELECTION
  // 0x01: 2 ``hadron'' tracks in pV and no more.
  //      (Note: The above flags excl. phi. But incl. phi is the 0x400 infra.
  //       The remaining being in the quasi-excl. phi w/ "(phiPat&0x401)==0".)
  // 0x02: h+h-				     
  // 0x04: pT > cut
  // 0x08: dEK < cut
  // 0x10: h1 w/in RICH acceptance and  piTHr<P1<Pmx
  // 0x20: h2 w/in RICH acceptance and  piTHr<P2<Pmx
  // 0x40: KID of h1 (above or below KThr; piTHr*f<P1<Pmx)
  // 0x80: KID of h2
  // 0x100: No detached track compatible w/ vertex
  // 0x200: No, in time, activity in ECALS
  // 0x400: Incl. phi (Note: 0x300 then meaningless.)
  // 0x800: rho
  UShort_t phiPat;
  // ***** V0 SELECTION
  // 0x01: Vertex distance D/dD
  // 0x02: Collinearity
  // 0x04: pT > cut
  // 0x08: Chi2 of Vertex
  // 0x10: h1 w/in RICH acceptance & piTHr[*f]<P1<Pmx (Lambda's f,Pmx != K0's)
  // 0x20: h2 w/in RICH acceptance & piTHr[*f]<P2<Pmx
  // 0x40: PID of h1
  // 0x80: PID of h2
  // 0x100: Tighter selection: D>>dD, Zs>ZT => disregard upstream of target
  // 0x200: eVeto of h1 if pion
  // 0x400: eVeto of h2 if pion
  // 0x800: Histo range OK for \bar{Lambda}
  UShort_t K0Pat, LambdaPat;

  Short_t  h1, h2; // Indices in CSHadronData
  Float_t  m;
  Float_t  Xs, Ys, Zs; // sVertex
  Float_t  D, dD; // Vertex distance and uncertainty
  Float_t  cth;   // Cosine of collinearity
  Float_t  pT;    // pT w.r.t. to the pair system
  Float_t  alpha; // Armenteros alpha
  Float_t  chi2;  // Vertex chi2 (phi: pVertex, V0's: sVertex)

  void Reset()
  {
    phiPat = 0; K0Pat = 0; LambdaPat = 0;
    h1 = -1; h2 = -1; m = 0; 
    Xs = 0; Ys = 0; Zs = 0;     // sVertex
    D = 0; dD = 0; cth = 0; pT = 0; alpha = 0; chi2 = 0;
  }

  ClassDef(CSResonanceData,3); // Must be the last item before the closing '};'
};

// The following in order to avoid: 
// "Error in <TTree::Branch>: The pointer specified for Hadrons is not of a class or type known to ROOT"

/*
  ClassDef(std::vector\<CSHadronData\>,1)
*/

/*
  ClassDef(std::vector\<CSResonanceData\>,1)
*/

#endif
/*
  gSystem->Load("libPhysics.so");
  gSystem->Load("libGraf.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("lib/libPhast.so");
  gSystem->Load("user/libUser.so");

  TCut pWinR = "Hs[h1].XR*XR+YR*YR>4 && abs(Hs[h1].XR)<161 && abs(Hs[h1].YR)<121 && Hs[h1].ZFirst<750 && Hs[h1].ZLast>360"; // ZRICH ~= 750, ZSM1 ~= 360

  //static double LHCut = 1.01, LHBckCut = 1.20;
  //static double subKThrLHpiVeto = 0.92;
  //static double thrMargin = 1.05; // Safety margin: require P > Thr*margin
  TCut KpID = "Hs[h1].LH[1]>1.01*Hs[h1].LH[0] && Hs[h1].LH[1]>1.01*Hs[h1].LH[2] && Hs[h1].LH[1]>1.2*Hs[h1].LH[5] && piThr*3.53712<Hs[h1].qP && Hs[h1].qP<60"&&pWinR;
  // Note: for subpThr PID, we allow for RICH being absent (i.e. either no RICH
  // block or background LH =0) to mean lack of response and hence that particle
  // is below threshold. This must be resticted to particles:
  // - w/in the acceptance of the RICH,
  // - w/in reach of the RICH sofware, i.e. ZFirst<< and ZLast>> and P<<.
  // Note that in earlier software versions, processing of events w/o photon,
  // not even background ones, was given up, yielding tracks w/o RICH block.
  //static double pIDPmax = 80;
  //static double pIDPmin     = 1.05 (* times pThr)
  //static double subpThrPmin = 1.1  (* times piThr)
  TCut KpSubID = "(!Hs[h1].hasR && piThr<Hs[h1].qP || Hs[h1].LH[0]<0.92*Hs[h1].LH[5] && Hs[h1].LH[3]<0.92*Hs[h1].LH[5] && piThr*1.05<Hs[h1].qP) && Hs[h1].qP<piThr*3.53712"&&pWinR;

  TCut pipID = "Hs[h1].LH[0]>1.01*Hs[h1].LH[1] && Hs[h1].LH[0]>1.01*Hs[h1].LH[2] && Hs[h1].LH[0]>1.2*Hs[h1].LH[5] && piThr<Hs[h1].qP && Hs[h1].qP<60"&&pWinR;
  TCut peVeto = "!Hs[h1].hasR || Hs[h1].LH[3]==0 || Hs[h1].LH[3]<1.5*Hs[h1].LH[0]"&&pWinR;

  // Note: pID includes pLH>eLH, which may not be relevant at higher P
  TCut ppID = "Hs[h1].LH[2]>1.01*Hs[h1].LH[0] && Hs[h1].LH[2]>1.01*Hs[h1].LH[1] && Hs[h1].LH[2]>1.01*Hs[h1].LH[3] && Hs[h1].LH[2]>1.2*Hs[h1].LH[5] && piThr*6.72258*1.05<Hs[h1].qP && Hs[h1].qP<80"&&pWinR;
  TCut ppSubID = "(!Hs[h1].hasR || Hs[h1].LH[0]<1.4*Hs[h1].LH[5] && (Hs[h1].qP<piThr*3.53712 || Hs[h1].LH[1]<1.5*Hs[h1].LH[5]) && Hs[h1].LH[3]<1.5*Hs[h1].LH[5]) && piThr*1.1<Hs[h1].qP && Hs[h1].qP<piThr*6.72258*1.05"&&pWinR;

  TCut mWinR = "Hs[h2].XR*XR+YR*YR>4 && abs(Hs[h2].XR)<161 && abs(Hs[h2].YR)<121 && Hs[h2].ZFirst<750 && Hs[h2].ZLast>360";
  TCut KmID = "Hs[h2].LH[1]>1.01*Hs[h2].LH[0] && Hs[h2].LH[1]>1.01*Hs[h2].LH[2] && Hs[h2].LH[1]>1.2*Hs[h2].LH[5] && -60<Hs[h2].qP && Hs[h2].qP<-piThr*3.53712"&&mWinR;
  TCut KmSubID = "(!Hs[h2].hasR && Hs[h2].qP<-piThr || Hs[h2].LH[0]<0.92*Hs[h2].LH[5] && Hs[h2].LH[3]<0.92*Hs[h2].LH[5] && Hs[h2].qP<-piThr*1.05) && -piThr*3.53712<Hs[h2].qP"&&mWinR;

  TCut pimID = "Hs[h2].LH[0]>1.01*Hs[h2].LH[1] && Hs[h2].LH[0]>1.01*Hs[h2].LH[2] && Hs[h2].LH[0]>1.2*Hs[h2].LH[5] && -60<Hs[h2].qP && Hs[h2].qP<-piThr"&&mWinR;
  TCut meVeto = "!Hs[h2].hasR || Hs[h2].LH[3]==0 || Hs[h2].LH[3]<=1.5*Hs[h2].LH[0]"&&mWinR;

  TCut pmID = "Hs[h2].LH[2]>1.01*Hs[h2].LH[0] && Hs[h2].LH[2]>1.01*Hs[h2].LH[1] && Hs[h2].LH[2]>1.01*Hs[h2].LH[3] && Hs[h2].LH[2]>1.2*Hs[h2].LH[5] && -80<Hs[h2].qP && Hs[h2].qP<-piThr*6.72258*1.05"&&mWinR;
  TCut pmSubID = "(!Hs[h2].hasR || Hs[h2].LH[0]<1.4*Hs[h2].LH[5] && (-piThr*3.53712<Hs[h2].qP || Hs[h2].LH[1]<1.5*Hs[h2].LH[5]) && Hs[h2].LH[3]<1.5*Hs[h2].LH[5]) && -piThr*6.72258*1.05<Hs[h2].qP && Hs[h2].qP<-piThr*1.1"&&mWinR;

*/

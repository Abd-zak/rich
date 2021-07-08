// $Id: MCInfo.h,v 1.16 2020/07/30 10:46:55 ybedfer Exp $


#ifndef MCInfo_h
#define MCInfo_h
/*!
  \class MCInfo
  \brief Sort out usefull from MC objects:<br>
  - Reconstructibility.
  - Kinematics.

  \author Yann.Bedfer@cern.ch
*/

#include <iostream>
#include <cmath>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "Phast.h"

class MCInfo {
 public:
  //! Constructor. Book histograms.
  MCInfo(int year);

  //! Returns pointer to MCInfo 
  static MCInfo *Ptr()
  {
    return address;
  }

  void HitsPerZone(PaEvent &e, int imct, int imct2);

  //! Book ``denominator'' histograms for tracking efficiency + some initialisations
  void BookTrackEfficiency(int nQ2, double *Q2bins, int nxB, double *xBbins,
			   int ny, double yMin, double yMax);

  //! Process MC raw info of argument event and fill reconstructibility
  // histograms.
  void SortOut(PaEvent &e,
	       bool looseBremsstrahlCut = false,
	       bool piK2muEval = true); // Evaluate reco'ibility on mu continuation of pi or K


  //! Book ``denominator'' histograms for D0 efficiency + some initialisations
  void BookD0Efficiency(double cLowS, double cUpS, double zS,
			double cLow0, double cUp0, double z0,
			double yLow, double yUp, double R2,
			bool richAngleCut, double richCut);
  //! Process MC raw info of argument event for D0
  void SortD0Out(PaEvent &e, bool richAcceptance);
  //! Process MC LUND info of argument event for D0
  void SortD0Out_LUND(PaEvent &e, bool richAcceptance);
  //! Subroutines of SortD0Out[_LUND]
  void D0EffHistosAndFlags(PaEvent &e, int C);
  void DSEffHistosAndFlags(PaEvent &e, int C);

  //! Book ``denominator'' histograms for Lambda acceptance + some initialisations
  void BookLambdaAcceptance(double yLow, double yUp, double R2,
			    double pTCut);
  //! Process MC raw info of argument event for Lambda
  void SortLambdaOut(PaEvent &e, bool richAcceptance);

  //! Book gamma histograms
  void Bookgamma();
  //! gamma reconstructibility
  void SortgammaOut(PaEvent &e);
  //! gamma efficiency and reco accuracy
  void Fillgamma(PaEvent &e, int imctp, int imcte, double pg);

  //! Book 3/5pi histograms
  void BookNpi(int Npis);
  //! 3/5pi analysis
  void SortNpiOut(PaEvent &e);

  //! Kinematical and event cuts
  double cLowCutS, cUpCutS, zCutS;   // Kinematical cut for the D* selection
  double cLowCut0, cUpCut0, zCut0;   // Kinematical cut for the D0 selection
  double yLowCut, yUpCut, R2Cut;     // Event cut
  bool RICHAngleCut; double RICHCut;
  double V0pTCut, V0cthCut;

  // ***** scattered mu = iMCT_mu. If none, = -1 and kinematics undef'd *****
  int iMCT_mu; double Pmup;
  //! Kinematics
  double Q2, yB, xB, nu, sqrts, XpV, YpV, ZpV;
  TLorentzVector lvq;
  TVector3 betacm, vqUnit;

  // ***** (a)D0 *****
  int iMCT_D0[2]; //! Track index of (a)D0: obsolete because not always available, e.g. when (a)D0 is made to decay in the MC generator (as opposed to COMGeant)
  unsigned int D0Type[2]; //! Type flag of the, possibly not available, (a)D0 MC track
  int iMCD0_K[2], iMCD0_pi[2]; //! MC decay track indices for (a)D0 -> Kpi[+pi0]
  double PD0[2], PD0_K[2], PD0_pi[2];
  double zD0[2], cthstarD0[2];
  double D0Decay[2]; // 0x1: K-pi+, 0x2: K-pi+pi0, 0x4: K-rho+(pi+pi0), 0x8: K*-(K-pi0)pi+, 0x10: K*0(K-pi+)pi0, 0x100: K*(892) (Note: K*0(892)->K-pi+ is outside our Kpi range), 0x200: K*(1430), 0x400: K*(1680), ..., 0x1000: K-pi+pi+pi-?
  // ***** D*+/- *****
  int iMCT_DS[2], iMCDS_D0[2]; //! Track indices for D*+/- and D0(<-D*): obsolete because not always available
  unsigned int DSType[2]; //! Type flag of the, possibly not available, D* MC track
  int iMCDS_piS[2]; //! MC track index for soft pi (piS) from D*->D0pi
  double PDS[2];

  // ***** (anti-)Lambda
  int iMCT_L[2], iMCL_p[2], iMCL_pi[2];
  double PL[2], PL_p[2], PL_pi[2];
  double xFL[2], cthstarL[2];
  // ***** Sigma*/Xi/Omega -> (anti)Lambda + pi
  int iMCT_S[2], iMCS_L[2], iMCS_pi[2];
  double PS[2];

  // ***** gamma *****
  vector<int> iMCT_gamma, iMCgamma_e,  iMCgamma_p, gammaHasAvatar;

  // ***** 3/5pi *****
  int MCNpis; // 3 or 5
  int iMCNpi_pims[3], iMCNpi_pips[2];
  double PMCNpi_pims[3], PMCNpi_pips[2]; // Momentum @ 1st encountered detector
  double MMCNpi, TMCNpi; // Mass, t of the Npi system
  double mpippimMCNpi[2];

  //! Accessors
  int *GetTypes() const { return Types; };
  int *GetRecos() const { return Recos; };
  int **GetNHits() { return NHits; };

  //#define DEBUG_RECOIBILITY
#ifdef DEBUG_RECOIBILITY
  int ribilityFilled;
#endif	

 private:
  static MCInfo *address;        //! (-) address of this object

  //! Enable track efficiency ``denominator'' histograms
  bool TrackHistos;

  //! Histograms for track reconsructibility
  TH1D *hf_pt;                      // Generated: primaries
  TH1D *hf_at;                      // Acceptance
  TH1D *hf_Rt, *hf_lRt, *hf_sRt, *hf_aRt;// Tracking (all,LAS,SAS,LSAS)...
  TH1D *hf_fRt;                            // ...and FF
  TH1D *hf_jRt;                            // ...and MWA
  TH2D *hf_jID;                            // ...MWA as a f(ID)
  TH1D *hf_pVt, *hf_pVtP, *hf_pVaP; // pVertexing: all/primary, vs. ThetaV...
  TH1D *hf_pVtF;                           // ...and primary FF
  // Scattered muon
  TH2D *hf_RQ2, *hf_Ry, *hf_RxB;   // 2nd criterion
  TH2D *hf_RyvsQ2[4];
  // 2D
  TH2D *haa_pVtP, *hpa_pVtP;
  // Reinteraction
  TH2D *hi_pID, *hi_iID; // Rate vs. ID: p = all primaries, i = reinteracting
  TH2D *hi_pId, *hi_iId; // Rate in downstream (of center) part of target
  TH2D *hi_IID, *hi_IId; // Excluding soft reinteraction
  TH1D *hi_yID, *hi_ypT; // Yield of reinteraction vs. ID, vs. pT

  //! Histograms for D0->Kpi reconstructibility, vs. pD0, pK, cos(theta*) and zD0 (w/ or w/o cut on theta*, w/ all kin cuts)
  TH1D *hp_aD0[4],*hp_RD0[4], *hp_aD0c[4],*hp_RD0c[4], *hp_aD0k[4],*hp_RD0k[4];
  TH1D *hp_aD0ki[4],*hp_RD0ki[4], *hp_aD0kI[4],*hp_RD0kI[4];
  // ***** D0->Kpi ACCEPTANCE w/ NOT ALL KIN CUTS
  TH1D *hp_aD0cz[4];
  //! Histograms for D*->D0pi reconstructibility, vs. pD0, pK, cos(theta*) and zD0
  TH1D *hp_aDS[4], *hp_RDSk[4], *hp_RDSki[4], *hp_RDSkI[4];
  // ***** KINEMATICAL CUT *****
  TH1D *hD_cs; 

  //! Histograms for Lambda (Sigma*) acceptance, vs. cos(theta*) and vs ppi
  TH1D *hc_aL[2][2], *hc_aLk[2][2], *hc_aLi[2][2], *hc_aLki[2][2];
  TH1D *hc_RL[2][2], *hc_RLk[2][2], *hc_RLi[2][2], *hc_RLki[2][2];
  TH1D *hc_aS[2], *hc_RSki[2];
  TH1D *hD_ReLambda;

  //! Histograms for gamma
  TH1D *hp_ag, *hp_Ag, *hv_Ag, *hp_ACg, *hv_ACg, *hp_Irgee, *hp_Frgee;
  TH2D *hp_Agee, *hp_Rgee, *hv_Rgee, *hp_Rrgee;
  TProfile *pe_gdpvsp[3], *pe_gdpvsZ[3];

  //! Histograms for 3/5pi
  TH2D *hf_Dalitz; TH1D *hf_tNpi, *hf_mNpi, *hf_pNpi;
  //! Zone boundaries: IDs of LAST detector (0|target|1|SM1|2|SM2|3|muWall2|4)
  int iZoneLasts[5];
  //! Sub-zone of MuWallA: IDs of FIRST and LAST MAs, plane pattern
  int FirstMuWallA, LastMuWallA; unsigned int PlPatMuWallA;

  //! Sorted Track info: type of tracks, #hits per zone

  //  0x0001:   Primary particle
  //  0x0002:   gamma decay
  //        :   pVertexible, i.e. either primary or secondary w/ mother vertex
  //            close to primary.
  //#define PSEUDO_pV_VIA_VERTEX
#ifdef PSEUDO_pV_VIA_VERTEX
  //  0x0004:   Close in term of vertex resolution => dZ<1.0~=2sigma as measured
  //            CG,v7 D0s and dR<.025
  //  0x0008:   dZ<2.0 and dR<.05
#else
  //  0x0004:   Close in term of pseudo-CDA => dR<.02 (as measure on D*)
  //  0x0008:   dR<.04
#endif

  //  0x0010:    Fringe Field reconstructibility
  //  0x0020:    LAS Reconstructibility criterion
  //  0x0040:    SAS Reconstructibility criterion
  //  0x0080;    Any Reconstructibility
  //  0x0100:    Fringe Field Reconstructed

  //  0x0200:    LAS Reconstructed
  //  0x0400:    SAS Reconstrctible & Reconstructed
  //  0x0800:    pVertexible & pVertexed
  //  0x1000:    Reco'd

  //  0x4000:  Pile-up

  //  0X00010000: Scattered muon
  //  0x00020000: Scattered muon reco'ibility
  //  0x00040000: MuID'ibility in MWB
  //  0x00080000: MuID'ibility in MWA

  //  0x00100000: (a)D0->Kpi | (a)Lambda->ppi
  //  0x00200000: (a)D0               w/in cth* cut (D0)
  //  0x00400000: (a)D0 | a(Lambda)   w/in kin cuts: cth*,z cuts (D0); pT (L)
  //  0x00800000: (a)D0               w/in loose kin cth*,z cuts (for D* sel'n)
  //  0x01000000: K <- (a)D0              or   p  <- (a)Lambda
  //  0x02000000: pi <- (a)D0                  pi <- (a)Lambda
  //  0x04000000: pi <- D*+/-                  pi <- Sigma*/Xi/Omega
  //  0x10000000: RICH ID'able
  //  0x20000000: with RICH Cut (angle or radius)
  //  0x40000000; K RICH ID'd
  //  0x80000000; pi RICH misID'd as K
  int *Types, *Recos, *NHits[5];

  int Year;
  
  // RICH efficiencies
  int nRICHaBins, nRICHpBins;
  TMatrixD eRICHKm, eRICHKp, eRICHpim, eRICHpip;
  double *RICHaBins, *RICHpBins;
  double RICHpLow;
};

// Alex 2006 efficiencies. Retrieved from sophie's UserEvent22.cc on 2010/01.
const int nAL06aBins = 8, nAL06pBins = 9;
const double eAL06Kp[nAL06pBins*nAL06aBins] = {
  0,     0,     0,     0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,     0,     0,     0,
  0.860, 0.669, 0.779, 0.971, 0.965, 0.968, 0.968, 0.968,
  0.629, 0.466, 0.832, 0.986, 0.979, 0.977, 0.982, 0.955,
  0.385, 0.384, 0.875, 0.984, 0.982, 0.977, 0.980, 0.977,
  0.251, 0.331, 0.873, 0.986, 0.972, 0.977, 0.976, 0.976,
  0.100, 0.265, 0.828, 0.970, 0.972, 0.972, 0.972, 0.972,
  0.023, 0.187, 0.784, 0.940, 0.931, 0.928, 0.928, 0.928 };

const double eAL06Km[nAL06pBins*nAL06aBins] = {
  0,     0,     0,     0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,     0,     0,     0,
  0,     0,     0,     0,     0,     0,     0,     0,
  0.815, 0.581, 0.801, 0.967, 0.962, 0.965, 0.969, 0.956,
  0.532, 0.407, 0.877, 0.983, 0.971, 0.969, 0.981, 0.960,
  0.332, 0.349, 0.915, 0.983, 0.977, 0.979, 0.990, 0.981,
  0.184, 0.315, 0.916, 0.983, 0.981, 0.981, 0.981, 0.981,
  0.075, 0.279, 0.887, 0.976, 0.968, 0.921, 0.965, 0.965,
  0.025, 0.223, 0.866, 0.938, 0.894, 0.892, 0.892, 0.892 };

const double eAL06pip[nAL06pBins*nAL06aBins] = {
  0.838, 0.839, 0.858, 0.852, 0.867, 0.935, 0.944, 0.931,
  0.862, 0.860, 0.869, 0.869, 0.933, 0.949, 0.939, 0.945,
  0.920, 0.914, 0.846, 0.951, 0.967, 0.964, 0.950, 0.950,
  0.937, 0.799, 0.874, 0.983, 0.982, 0.976, 0.980, 0.973,
  0.717, 0.658, 0.912, 0.983, 0.980, 0.974, 0.978, 0.976,
  0.509, 0.617, 0.936, 0.979, 0.976, 0.970, 0.972, 0.964,
  0.379, 0.591, 0.939, 0.970, 0.967, 0.962, 0.945, 0.941,
  0.253, 0.538, 0.863, 0.853, 0.877, 0.864, 0.875, 0.868,
  0.081, 0.463, 0.664, 0.530, 0.603, 0.600, 0.597, 0.597 };

const double eAL06pim[nAL06pBins*nAL06aBins] = {
  0.809, 0.837, 0.872, 0.850, 0.866, 0.939, 0.948, 0.940,
  0.814, 0.854, 0.864, 0.865, 0.941, 0.954, 0.950, 0.947,
  0.903, 0.892, 0.834, 0.958, 0.970, 0.967, 0.960, 0.953,
  0.897, 0.763, 0.881, 0.985, 0.982, 0.978, 0.984, 0.978,
  0.658, 0.645, 0.923, 0.985, 0.978, 0.977, 0.984, 0.978,
  0.493, 0.632, 0.952, 0.982, 0.977, 0.977, 0.987, 0.967,
  0.347, 0.621, 0.956, 0.978, 0.974, 0.971, 0.983, 0.982,
  0.241, 0.567, 0.884, 0.871, 0.891, 0.892, 0.847, 0.847,
  0.150, 0.519, 0.664, 0.554, 0.630, 0.630, 0.625, 0.625 };

const double AL06aBins[nAL06aBins] = {
  0.005, 0.010, 0.020, 0.030, 0.040, 0.060, 0.090, 0.120 };
const double AL06pBins[nAL06pBins] = {
  5.0, 7.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0 };

#endif
// $Log: MCInfo.h,v $
// Revision 1.16  2020/07/30 10:46:55  ybedfer
// No change.
//
// Revision 1.14  2020/07/06 22:03:15  ybedfer
// - Merge v1.14 from CERN and 1.13 from CCIN2P3.
// - Bug fix in reco'ibility criterion for D*: piS has to be (L|S)AS reco'ible
//  (or if FF-reco'ible if FF not rejected !(pi_FF_REJECTION&0x2).
//
// Revision 1.14  2020/07/06 21:45:21  ybedfer
// No change.
//
// Revision 1.13  2020/05/26 00:22:12  ybedfer
// Improved C++ programming.
//
// Revision 1.13  2016/12/06 00:32:51  ybedfer
// - SortOut:
//  - Add argument "piK2muEval", in order to be able to disable that feature (
//   i.e. evaluate reco'ibility of piK tracks + their possible mu extension).
//  - assert ->abort.
// - SortD0Out: Bug fix: Reco'ibility of piS(<-D*) corresponds to 0x60 (and no
//  longer to 0x80 (since that meaning, 0x80 <-> (L|S)AS reco, is no longer
//  valid).
// - SortLambdaOut: LorentzVector "lv(1|2)" of decays no longer references (
//  note that I don't remember why...).
// - "Get(Types|Recos)" now const.
//
// Revision 1.12  2015/12/31 20:36:04  ybedfer
// - Changing the meaning of "0x80": used to be "0x2|0x4" (not very useful),
//  now = "resonance reco'ibility".
// - Reinteraction: fill histos/flags, in module "flagReinteraction".
// - Member arrays of PaMCtrack attributes:
//   - "Recos", declared next to "Types". And init'd = 0.
//   - "nHits" -> "NHits".
// - MWA:
//   - New "MuWallA"-type variables (including HG02).
//   - New histos "_jR" = Reco'ible&Reco'd MWA (also "hf_jID").
//   - MWA(|B) reco'ibility flag (0x00040000|0x00080000) extended to all-out (
//    not restricted to mu') case.
// - "BookTrackEfficiency" now also as a f(y).
// - More modularity: "HitsPerZone", "iPiKDecay".
// - "MC_GEOM_ACCEPTANCE": revisited, no loger default.
// - cpp:
//   - New "MC_DEBUG_pimu".
//   - "MC_DEBUG_ASSOCIATION" cleaned.
// - Bremsstrahlung: somewhat (how much?...) revisited.
//
// Revision 1.10  2012/05/27 01:08:14  ybedfer
//  - "3pi" -> "Npi", valid for 5 as well as 3 pions: "MCNpis = 3|5".
//    - BookNpi(int Npis)", "SortNpiOut".
//  - "hf_pVaP", analogue "hf_pVtp" for angle instead of momentum.
//  - Axis title for histo "hf_".
//
// Revision 1.9  2012/05/26 23:46:32  ybedfer
//  - "MC_GEOM_ACCEPTANCE":
//    - cpp directive.
//    - When defined, reco'ibility = primaries w/in acceptance.
//    - In particular, accept when no hit whatsoever or if decay vertex < DC01,
//     < ST03, < GM09. But require particle w/in ST03, w/in PA05.
//  - Get more from "PaSetup":
//    - "zSM1", "mTarget" (target magnet).
//    - Cancel cpp "REQUIRE_WIN_SMC"; use instead "mTarget".
//    - Set z/xFI06 depending upon "mTarget".
//  - Reco'ible tracks as a f(theta).
//  - 180 -> 200 GeV.
//  - 3pi: accept mu as beam.
//
// Revision 1.8  2011/05/07 16:12:04  ybedfer
//  - RICH efficiency matrices.
//  - New methods "Book3pi" and "Sort3piOut".
//
// Revision 1.7  2010/02/01 01:17:13  ybedfer
//  Fill c(theta*) and zD0 distributions, in addition to pD0 and pK.
//
// Revision 1.6  2009/12/08 11:21:50  ybedfer
//  - "MCInfo::SortD0Out_LUND" when generator info is available.
//  - New MCInfo members "D0Type[2]" and "DSType[2]":
//    - No longer have to check that D0="MC->iMCT_D0[i]" is defined (i.e.
//     >=0) prior to using it as the index into the "MCTypes" array in
//     order to check the reco'ibility of D0.
//     (Btw, is it useful to check D0 independently of K and pi? Cf. the
//     block bookkeeping reco'iblities in array "Rs[2][3]".)
//  - "iMCT_D0" and "iMCT_DS" obsolete.
//  - New members "D0Type[2]" and "DSType[2]": Used to store the type (in
//   the "MCInfo::Types" sense) of (a)D0 and D*+/-, now that we allow their
//   PaMCtracks not to exist (which is the case when they are decayed in
//   the physics generator, as opposed to COMGeant).
//  - New members "D0Decay".
//  - Fix: the retrieving the PaMCtrack indices of D0 and piS from D* decay
//   was wrongly coded (but the bug was harmless, due to the fact that the
//   D0 track comes always first n the list of event's PaMCtracks).
//  - New methods "D0/DSEffHistosAndFlags".
//    - Histos, reco'ibility flags ("Types" and "D0/DSType").
//    - Also momenta of D0, D* and K,pi decays. For D0 and D*, the value is
//     computed from the decays, even if their very PaMCtracks exist.
//    - Also cos(theta*).
//    - Check mass of D*, as reco'd from its decay PaMCtracks.
//  - New "MCInfo::SortD0Out_LUND" method:
//    - It parallels "SortD0Out", except that it retrieves the info from the
//     LUJET blocks.
//    - It includes a better monitoring of D0->Kpipi0, w/ provision for
//     2-body decay into an intermediate resonance: Krho, K*pi and K*0pi0.
//      - All K* excitation are considered.
//      - The result of the monitoring is stored in new "D0Decay" member.
//
// Revision 1.5  2009/12/03 00:49:23  ybedfer
// No change w.r.t. v1.4.
//
// Revision 1.4  2009/07/03 14:02:54  ybedfer
//  - Include "Lambda" analysis.
//  - More data members for event kinematics.
//
// Revision 1.3  2009/06/07 11:18:39  ybedfer
//  Adding new member "Pmup".
//
// Revision 1.2  2008/04/16 18:09:31  ybedfer
//  - Redefinition of "Types":
//    - 0x00002000: mu+/-
//    - 0x00080000: MuID'ibility in MWA
//
// Revision 1.1  2004/02/12 17:07:43  ybedfer
// Initial revision
//

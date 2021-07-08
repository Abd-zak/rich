// $Id: MCInfo.cc,v 1.16 2020/07/30 10:46:34 ybedfer Exp $

#include <cassert>
#include "MCInfo.h"
#include "Masses.h"
#include "TLorentzVector.h"
#include "PaAlgo.h"

// ***** cpp macros
// - "MC_HADRON_SETUP": to not require a muon to be the scattered particle
//#define MC_HADRON_SETUP

void setAxisTitle(const char*);

static const int nHitMapWs = 13; // Cf. "CSTRACK_MAPSIZE" in CORAL's "CsTrack".

MCInfo* MCInfo::address = 0; // Init static pointer

MCInfo::MCInfo(int year)
{
  if (address) {
    printf("** MCInfo::MCInfo:\a MCInfo already instantiated\n"); assert(false);
  }
  address = this;

  // *************** VARIOUS OTHER INITIALISATIONS **********
  Year = year;
  Types=Recos = 0;// ********** INIT "Types" and "NHits" BIT-PATTERNS **********
  for (int i = 0; i<5; i++) NHits[i] = 0;

  // ********** ZONE BOUNDARIES **********

  // They will serve to speed up the computation of the #hits per zone
  // (stored in "NHits[]")
  const PaSetup &setup = PaSetup::Ref();
  PaMagInfo *m = setup.PtrMagField()->getMagInfo();
  int nMags = setup.PtrMagField()->getNumOfMags(), iSM2 = nMags-1;
  int iSM1 = iSM2-1; double zSM1 = m[iSM1].zcm/10, zSM2 = m[iSM2].zcm/10;
  printf("\n * MCInfo: SM1 %.0f cm SM2 %.0f cm\n",zSM1,zSM2);
  const double zmuW2 = 3400;
  for (int i = 0; i<5; i++) iZoneLasts[i] = -1;
  FirstMuWallA=LastMuWallA = -1; PlPatMuWallA = 0;
  float zTarget = setup.TargetCenterZ();
  int iw, idet; for (iw = 0, idet = -1; iw<nHitMapWs; iw++) {
    for (int ib = 0; ib<32; ib++) {
      if (++idet>=setup.NDetectors()) continue;
      const PaDetect &d = setup.Detector(idet); double z = d.Z();
      if (d.Name().find("MA")==0 || d.Name().find("HG02")==0) {
	if (FirstMuWallA<0) FirstMuWallA = idet; LastMuWallA = idet;
	PlPatMuWallA |= 1<<idet-FirstMuWallA;
      }
      if      (z<zTarget) iZoneLasts[0] = idet;
      else if (z<zSM1)    iZoneLasts[1] = idet;
      else if (z<zSM2)    iZoneLasts[2] = idet;
      else if (z<zmuW2)   iZoneLasts[3] = idet;
      else                iZoneLasts[4] = idet;
    }
  }
  for (int i = 0; i<5; i++) {
    if (iZoneLasts[i]==-1) {
      printf("\n** MCInfo:\a No iZoneLasts for zone %d\n",i); exit(1);
    }
    for (int j = i+1; j<5; j++) if (iZoneLasts[j]<iZoneLasts[i]) {
      printf("\n** MCInfo:\a Inconsistent iZoneLasts [%d]=%d < [%d]=%d\n",
	     j,iZoneLasts[j],i,iZoneLasts[i]); exit(1);
    }
  }
  if (FirstMuWallA==-1) {
    printf("\n** MCInfo:\a MuWallA sub-zone undefined\n"); exit(1);
  }

  TrackHistos = false;  // ***** DISABLE TRACK EFFICIENCY HISTOGRAMS *****

  //           ********** RICH EFFICIENCIES **********
  // Fill member data w/:
  //  - Either of the efficiencies available on the market, depending upon
  //    - argument "year",
  //    - cpp macro "MC_RICH_PERFS".
  //  - Or, default, ideal RICH.
  TMatrixD *eRICHs[] = { &eRICHKp, &eRICHKm, &eRICHpip, &eRICHpim};
  if (!year) {                  // ***** Default: IDEAL RICH *****
    nRICHaBins=nRICHpBins = 1;
    for  (int Kpipm = 0; Kpipm<4; Kpipm++) {
      TMatrixD *eRICH = eRICHs[Kpipm];
      eRICH->ResizeTo(1,1); (*eRICH)(0,0) = 1;
    }
    RICHpBins = new double[1]; RICHaBins = new double[1];
    RICHaBins[0]=RICHpBins[0] = 0;
    RICHpLow = 0; // Means binning is expressed as a f(P-Threshold)
  }
  else {
    int error = 0; string *RICHPerfs = 0;
#define MC_RICH_PERFS 1 // 1=ALEX
#if MC_RICH_PERFS == 1
    // Cf. "~korzenev/public/RICH_prob/MakePurBin_Prob_2006_orig_2_P.txt"
    // and "~compass/publication/results/2008/september_deltaq/SIDIS_release_note.pdf"
    RICHPerfs = new string("Alex");
    if (year==2006) {
      error = 0;
      nRICHaBins = nAL06aBins; nRICHpBins = nAL06pBins;
      const double *eAL06s[] = { eAL06Kp, eAL06Km, eAL06pip, eAL06pim};
      for  (int Kpipm = 0; Kpipm<4; Kpipm++) {
	TMatrixD *eRICH = eRICHs[Kpipm]; const double *eAL06 = eAL06s[Kpipm];
	eRICH->ResizeTo(nRICHpBins,nRICHaBins);
	for (int ip = 0; ip<nRICHpBins; ip++)
	  for (int ia = 0; ia<nRICHaBins; ia++)
	    (*eRICH)(ip,ia) = eAL06[ip*nRICHaBins+ia];
	RICHpBins = new double[nRICHpBins];
	for (int ip = 0; ip<nRICHpBins; ip++) RICHpBins[ip] = AL06pBins[ip];
	RICHaBins = new double[nRICHaBins];
	for (int ia = 0; ia<nRICHaBins; ia++) RICHaBins[ia] = AL06aBins[ia];
      }
      RICHpLow = 3;
    }
#endif
    if (error==0) {
      if (!RICHPerfs) {
 printf("** MCInfo:\a RICH_PERFS inconsistency\n"); exit(1);
      }
      printf("** MCInfo:\a RICH_PERFS for year %d taken from \"%s\"\n",
	year,RICHPerfs->c_str());
      if (RICHpLow) printf("P    : %.3f",RICHpLow);
      else          printf("P-Thr: 0.00");
      for (int ip = 0; ip<nRICHpBins; ip++) printf(" < %.3f",RICHpBins[ip]);
      printf("\nA    : 0.000");
      for (int ip = 0; ip<nRICHaBins; ip++) printf(" < %.3f",RICHaBins[ip]);
      printf("\n");
      for (int Kpipm = 0; Kpipm<4; Kpipm++) {
	printf("%s%c(p,A)\n",Kpipm/2?"pi":"K",Kpipm%2?'-':'+');
	for (int ip = 0; ip<nRICHpBins; ip++) {
	  for (int ia = 0; ia<nRICHaBins; ia++) {
	    switch (Kpipm) {
	    case  0: printf(" %5.3f",eRICHKp [ip][ia]); break;
	    case  1: printf(" %5.3f",eRICHKm [ip][ia]); break;
	    case  2: printf(" %5.3f",eRICHpip[ip][ia]); break;
	    default: printf(" %5.3f",eRICHpim[ip][ia]);
	    }
	  }
	  printf("\n");
	}
      }
    }
    else if (RICHPerfs) {
 printf("** MCInfo:\a No \"%s\" RICH_PERFS available for year %d\n",
	RICHPerfs->c_str(),year); exit(1);
    }
    else {
 printf("** MCInfo:\a No RICH_PERFS available for year %d\n",year); exit(1);
    }
  }

  // ***** Reinteraction *****
  // Rate of reinteraction: p = all primaries, i/I = all/hard reinteractions
  // Independently downstream-most part of target
  double Pbins[] = {0,3.5,10,40}; int nPbins = sizeof(Pbins)/sizeof(double)-1;
  hi_pID = new TH2D("hi_pID","Prod. per ID vs. P",29,0.5,29.5,nPbins,Pbins);
  hi_iID = new TH2D("hi_iID","ReInter. per ID vs. P",29,0.5,29.5,nPbins,Pbins);
  hi_pId = new TH2D("hi_pId","Prod. per ID vs. P (Zv>>)",29,0.5,29.5,nPbins,Pbins);
  hi_iId = new TH2D("hi_iId","ReInter. per ID vs. P (Zv>>)",29,0.5,29.5,nPbins,Pbins);
  hi_IID = new TH2D("hi_IID","Hard ReInter. per ID vs. P",29,0.5,29.5,nPbins,Pbins);
  hi_IId = new TH2D("hi_IId","Hard ReInter. per ID vs. P (Zv>>)",29,0.5,29.5,nPbins,Pbins);
  // Yield of reinteraction: vs. ID and vs. pT
  hi_yID = new TH1D("hi_yID","ReInter. yield vs. ID",29,0.5,29.5);
  hi_ypT = new TH1D("hi_ypT","ReInter. yield vs. pT(#gamma*)", 80,0,4);
}

void MCInfo::BookTrackEfficiency(int nQ2, double *Q2bins,
				 int nxB, double *xBbins,
				 int ny, double yMin, double yMax)
{
  // ********** ``DENOMINATOR'' HISTOGRAMS FOR TRACKING EFFICIENCY **********
  TrackHistos = true;
  // Acceptance
  hf_pt  = new TH1D("hf_pt",  "Primary Tracks",                   400,0,200);
  hf_at  = new TH1D("hf_at",  "All Tracks",                       400,0,200);
  // Tracking
  hf_Rt  = new TH1D("hf_Rt",  "Reconstructible tracks",           400,0,200);
  hf_lRt = new TH1D("hf_lRt", "Reconstructible LAS tracks",       400,0,200);
  hf_sRt = new TH1D("hf_sRt", "Reconstructible SAS tracks",       400,0,200);
  hf_aRt = new TH1D("hf_aRt", "Reconstructible L+SAS tracks",     400,0,200);
  hf_fRt = new TH1D("hf_fRt", "Reconstructible fringe-field",     400,0,200);
  hf_jRt = new TH1D("hf_jRt", "Reco'ible muWallA tracks",         400,0,200);
  hf_jID = new TH2D("hf_jID", "Reco'ible muWallA tracks ID=?,#mu',#mu^{#pm},#pi^{#pm},K^{#pm},#gamma",400,0,200,9,-0.5,8.5);
  // Vertexing
  hf_pVtP= new TH1D("hf_pVtP","Reco'ible primary tracks",         400,0,200);
  hf_pVaP= new TH1D("hf_pVaP","Reco'ible primary tracks vs. #Theta_{V}",
		    200,0,200);
  hf_pVt = new TH1D("hf_pVt", "Vertexible tracks",                400,0,200);
  hf_pVtF= new TH1D("hf_pVtF","Vertexible fringe-field",          400,0,200);
  // ***** Scattered muon *****
  // (as a function of reo'ibility criterion: 0: none, 1: loose, 2: strict)
  hf_RQ2  = new TH2D("hf_RQ2","Reco'ible Q2",nQ2,Q2bins,8,-.5,7.5);
  hf_Ry =   new TH2D("hf_Ry", "Reco'ible yB",ny,yMin,yMax, 8,-.5,7.5);
  hf_RxB  = new TH2D("hf_RxB","Reco'ible xB",nxB,xBbins,8,-.5,7.5);

  // 2D y vs. Q2
  char h2N[] = "hf_RyvsQ2A"; const char labels[] = "aABL";
  const char *tags[] = {"all","#muID:B","#muID:A&!B","#muID:A"};
  char h2T[] = "Reco'ible y vs. Q2 - #muID:B&!A";
  double *q2bins = new double[nQ2/2+1];
  for (int bin = 0; bin<=nQ2; bin += 2) q2bins[bin/2] = Q2bins[bin];
  for (int muIDType = 0; muIDType<4; muIDType++) {
    sprintf(h2N,"hf_RyvsQ2%c",labels[muIDType]);
    sprintf(h2T,"Reco'ible y vs. Q2 - %s",tags[muIDType]);
    hf_RyvsQ2[muIDType] = new TH2D(h2N,h2T,nQ2/2,q2bins,ny/4,yMin,yMax);
  }
  delete q2bins;
  // 2D
  haa_pVtP = new TH2D("haa_pVtP","Reco'ible pTracks #ThetaY vs. #ThetaX @ Vtx",
		      60,-300,300,60,-300,300); 
  hpa_pVtP = new TH2D("hpa_pVtP","Reco'ible pTracks #Theta vs. P @ Vtx",
		      100,0,200,30,0,300); 
  setAxisTitle("P (GeV/c)");
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~        sortOut          ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int iPiKDecay(PaEvent &e, PaMCtrack &trk);
void flagReinteraction(PaEvent &e,
		       PaMCvertex vtx,        // Vertex of origin
		       int imct, int imother, // Particles's index and mother's
		       int *Types);
void MCInfo::SortOut(PaEvent &e,
		     bool looseBremsstrahlCut,
		     bool piK2muEval) // Evaluate reco'ibility on mu continuation of pi or K
{
  int runNum = e.RunNum(), evNum = (int)e.UniqueEvNum();
  // ********** PROCESS MC RAW INFO OF ARGUMENT EVENT **********

  // => Fill array of bit patterns "MCInfo::Types", flagging tracks
  //   for various reconstructibility criteria
  // => Fill arrays of MC hits per zone
  // => Fill, or rather init =0, array  "MCInfo::Recos".
  // => Fill reconstructibility histograms
  // => Fill reinteraction histo.

  // Reconstructibility:
  // - Standard Reconstructibility is requiring, in this "MCInfo", that track
  //  maintains its integrity (w/in a given reconstruction zone depending
  //  upon whether one considers LAS, SAS, etc...).
  // - This can be overridden (Caveat: work in progress...)
  // I) Arg. "looseBremsstrahlCut": Some loss, e.g. via bremsstrahlung, is still
  //   allowed. Which one can loosen here. The idea is to provide from a more
  //   detailed analysis of what's going on in the reco of e^+/e^-'s from gamma
  //   conversion, cf. "SortgammaOut".
  // II) cpp #define: "MC_GEOM_ACCEPTANCE", which only checks that particles are
  //   w/in the geometrical aperture of the detection (basically DC01).
  //   This allows to decompose the acceptance in three factors:
  //  1) Pure geometry (which can be calculated w/o any knowledge of COMPASS for
  //    any particular physics channel): N = w/in aperture, D = all tracks.
  //  2) Acceptance depending on COMPASS innards: including decay (which rate
  //    depends upon COMPASS compactness) and re-interaction (which rate depends
  //    upon COMPASS material): N = ideal reco, D = w/in aperture.
  //  3) Reconstruction efficiency, i.e. deficiency of the reconstruction
  //    software: N = reco, D = ideal reco.
  //   Note that this concerns primary tracks. We still apply requirements on
  //  "zFirst" to be upstream of some detector.
  //#define MC_GEOM_ACCEPTANCE

  // ***** GEOMETRY INFO for ACCEPTANCE-BASED RECONSTRUCTIBILITY  *****
  // - For when requiring w/in acceptance or no decay/interaction w/in abscissae
  //  range.
  // - "zTarget" is used to decide whether we're in 2006 (=> acceptance is w/in
  //  1st DC) or 2004 (=> w/in MMs and, upon define, also w/in SMC solenoid, if
  //  vertex of origin is upstream of magnet's exit).
  //  N.B.: - Would be better to look at how many DCs upstream of SM1.
  //        - Hadron setup not really considered. 2004 run will behave as
  //         a 2002:2004 muon setup (because target was @ a <0 abscissa).
  const PaSetup &setup = PaSetup::Ref();
  // ***** ABSCISSAE of MAGNETS
  const PaMagInfo *m = setup.PtrMagField()->getMagInfo();
  int nMags = setup.PtrMagField()->getNumOfMags(), iSM2 = nMags-1;
  int iSM1 = iSM2-1; double zSM1 = m[iSM1].zcm/10;
  int mTarget = 0; // Target magnet. Default = no magnet. 1=OD, 2=SMC.
  if (m[0].fsc) // Caveat: this does not provide for a magnet existing, but not powered
    mTarget = m[0].flg1/2==1 /* I.e. SMC: 2=soleno, 3=dipole */ ? 2 : 1;
  static bool targetMagnetWarning = true; if (targetMagnetWarning) {
    printf(" * MCInfo: Assuming %s target magnet\n",mTarget?(mTarget==2?"SMC":"OD"):"NO");
    targetMagnetWarning = false;
  }
  double zAcc, hwAcc;
  if (mTarget==2) { // 2002:2004 setup: MM01
    hwAcc = 20; zAcc = 142.0;
  }
  else {            // 2006 on: DC00
    hwAcc = 60; zAcc = 160.0;
  }
  double hwSMC = 15, zSMCExit = 118.4;

  const double zMM02V = 190; // Just upstream of MM02V. Assumed to define the...
  // ...ENTRANCE to LAS, i.e. the farthest abscissa upstream of which any given
  // track must originate in order to be reco'd in LAS.
  // DC01: 1/2 aperture is equated to #channels/2 * pitch. We ignore stagering:
  //      not straightforward whether it must be added or subtracted...
  const double hwDC01 = 0.7*(176/2);
  const double zDC01 = 260.;  // Approx.
  const double hwxST = 150.0, hwyST = 120.0, zST03U1 = 540.;  // Approx.
  const double hwFI05 = 4.0, zFI05 =  585.0;
#ifndef MC_GEOM_ACCEPTANCE
  // Caveat: parameters infra valid only for the muons (including DVCS 2009 = d)
  const double hwFI06 = 5.0, zFI06m = 1500.0, xFI06m = 3.0;
  const double               zFI06d = 1800.0, xFI06d = 1.0 /* Sic! cf. detectors.79433.mu+.dat */;
  double zFI06, xFI06; if (mTarget==0) { zFI06 = zFI06d; xFI06 = xFI06d; }
  else                                 { zFI06 = zFI06m; xFI06 = xFI06m; }
#endif
  const double zPS01 = 950., zGM09 = 2085.;
  const double zPA05 = 2075., hwxPA = 90., hwyPA = 60.;
  // We shall consider also tracks through ECAL, HCAL and muF central holes
  // Need their position and central hole: take it from
  // "~compass/delivery/comgeant/0-0-7.00/data/geom >geom_general_070.ffr"
  // MU1F:
  // position : 1431.5-1491.5 in the HALL  'BOX '  3   30.   70.   35. 
  const double zMF1 = 1461.5, hwxMF1 = 70, hwyMF1 = 35;   // cm, half sizes
  // HC1:
  // equal to 1195.-1340. in the HALL      'BOX '  3   10.   62.   31.
  const double zHC1 = 1267.0, hwxHC1 = 62, hwyHC1 = 31.;  // cm, half sizes

  static PaTPar heli2;

  // ********** ALLOCATE BIT PATTERNS... **********
  int nTrks = e.NMCtrack();
  //! Sorted Track info: type of tracks, #hits per zone

  //  0x0001:   Primary particle
  //  0x0002:   !e^+/- but from prim.gamma || (semi-)leptonic decay of neutrals (out of neutral, if gamma)
  //        :   pVertexible, i.e. either primary or secondary w/ mother vertex
  //            close to primary.
#ifdef PSEUDO_pV_VIA_VERTEX
  //  0x0004:   Close in term of vertex resolution => dZ<1.0~=2sigma as measured
  //            CG,v7 D0s and dR<.025
  //  0x0008:   dZ<2.0 and dR<.05
#else
  //  0x0004:   Close in term of pseudo-CDA => dR<.02 (as measure on D*)
  //  0x0008:   dR<.04
#endif

  //  0x00000010:    Fringe Field reconstructibility
  //  0x00000020:    LAS Reconstructibility criterion
  //  0x00000040:    SAS Reconstructibility criterion
  //  0x0080;        Reconstructibility of resonances (D0, Lambda, Sigma*,...)

  //  0x0100;    Reco'ible FF and only that(?)
  //  0x0800     Vertexible and vertexed in 'UserEvent4" (to be removed: bad use of the MCType concept...)

  //  0x00001000: pi+,K+,K0,p,n,L,S+/- and CC. Has reinteracted.

  //  0x00002000: Is a muon +/-

  //  0x00004000: Pile-up (i.e. any primary but 1st one, cf. "PaEventImportMC")

  //  0x00008000: Product of reinteraction

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
  //  0x10000000: RICH'able, i.e. w/in RICH acceptance
  //  0x20000000: RICH'able & with RICH Cut (angle or radius)
  //  0x40000000; K RICH ID'd
  //  0x80000000; pi RICH ID'd (D0) or misID'd (Lambda)
  // MC::Recos
  //  0x0010:    Fringe Field Reconstructed (and not more)
  //  0x0020:    LAS Reconstructible & Reconstructed
  //  0x0040:    SAS Reconstructible & Reconstructed
  //  0x0800:
  //  0x1000:    Reco'd
  //  0x00080000:MWA-ible and MWA reco'd
  if (Types) delete Types; Types = new int[nTrks];
  if (Recos) delete Recos; Recos = new int[nTrks];
  for (int i = 0; i<5; i++) {
    if (NHits[i]) delete [] NHits[i];
    NHits[i] = new int[nTrks];
  }
  for (int imct = 0; imct<e.NMCtrack(); imct++) {
    Types[imct] = 0; Recos[imct] = 0; //...INIT
  }

  int imcv, imcvP;
  for (imcv = 0, imcvP = -1; imcv<e.NMCvertex(); imcv++) {
    // ********** LOOP ON MCvertex's **********
    const PaMCvertex &vtx = e.vMCvertex()[imcv];
    if (vtx.IsPrimary()) {      // ***** FIND PRIMARY, SET PRIMARY'S XYZpV *****
      if (imcvP>=0) {
	printf("\n **MCInfo:\a Run %d Evt %d: 2 Primaries %d,%d => Exit...\n",
	       runNum,evNum,imcvP,imcv); abort();
      }
      XpV = vtx.Pos(0); YpV = vtx.Pos(1); ZpV = vtx.Pos(2); imcvP = imcv;
    }
  }
  if (imcvP<0) {
 printf("\n ** MCInfo:\a Run %d Evt %d: No Primary Vertex => Exit...\n",
	runNum,evNum); assert(false);
  }


#ifndef MC_HADRON_SETUP
  int pid1st = -1;
#endif
  int imct; for (imct = 0, iMCT_mu = -1; imct<e.NMCtrack(); imct++) {

    // *************** LOOP ON MCtrack's, Excluding incident mu ***************

    PaMCtrack  &trk = e.vMCtrack()[imct];

    //#define MC_DEBUG_pimu
#ifdef MC_DEBUG_pimu
    int pid = trk.Pid();
#endif


    // ********** RESET Track <-> MCtrack X-REFERENCES **********
    // - There can be several reco'd tracks associated to a MC track. Let's
    //  determine the best one. And assign it to "PaMCtrack::iBestTrkRef".
    // - Update the reverse association (track -> MC track) accordingly by
    //  retaining one at most "PaTrack::iMCtrack" (per MC track).
    set<int>::iterator it;
    int iBest; static double chi2Best; static int nHitsBest;
    //#define MC_DEBUG_ASSOCIATION
#ifdef MC_DEBUG_ASSOCIATION
    vector<int> iTrWMoms;
#endif
    for (it = trk.sTrkRef().begin(), iBest = -1; it!=trk.sTrkRef().end();
	 it++) {
      const PaTrack &t = e.vTrack(*it); const PaTPar &h = t.vTPar()[0];
      if (h.HasMom()) {          // ...REJECT ALL MOMENTUMLESS...
#ifdef MC_DEBUG_ASSOCIATION
	iTrWMoms.push_back(*it);
#endif
	double chi2 = t.Chi2tot()/(t.NHits()-5); int nHits = t.NHits();
	if (iBest==-1 ||
	    nHits>nHitsBest &&   // ...RETAIN THAT W/ LARGEST #Hits
	    chi2<chi2Best+5) {   // ...AND STILL REASONABLE chi2

	  iBest = *it; chi2Best = chi2; nHitsBest = nHits;
	}
      }
    }
#ifdef MC_DEBUG_ASSOCIATION
    if ((int)iTrWMoms.size()>1) {
      printf("\n");
      for (int i = 0; i<(int)iTrWMoms.size(); i++) {
	int iTr = iTrWMoms[i]; 
	const PaTrack &t = e.vTrack(iTr); const PaTPar &h = t.vTPar()[0];
	double chi2 = t.Chi2tot()/(t.NHits()-5); int nHits = t.NHits();
	printf("PaTrack %3d: %d Hits %.1f %.2f   %3d  %d %.2f\n",
	       iTr,t.NHits(),h(0),chi2,iBest,nHitsBest,chi2Best);
      }
    }
#endif
    // ***** SET Track <-> MCtrack X-REFERENCES *****
    trk.SetBestTrkRef(iBest);                           // MCtrack -> BEST Track
    for (it = trk.sTrkRef().begin(); it!=trk.sTrkRef().end(); it++) {
      PaTrack &t = const_cast<PaTrack&>(e.vTrack(*it));
      if (*it!=iBest) t.SetMCtrack(-1);                 // Track   -> MCtrack
    }
    if (imct==0) continue;

    // ********** FIRST: EXCLUDE UNINTERESTING MCTrack's:... **********
    // ...pileups, neutrals, e^+/- of the delta ray type which are
    // lumped by COMGeant into chimeric, tracks. In order to discard the
    // latter one discards all e^+/- except when primary or gamma decay or
    // what is tentatively determined to be (semi-)letonic decay of hadrons (to
    // this end one restricts this to the decay of neutrals so as not to be
    // confused w/ delta rays and, if neutral = gamma, one requires gamma to
    // originate from a neutral also).

    int ivtxm = trk.iVertexOfOrigin();
    const PaMCvertex &vtx = e.vMCvertex()[ivtxm];
    if (vtx.IsPileup()) {                                  // ***** NOT A PILEUP
      Types[imct] |= 0x00004000; continue;
    }

    //            ********** REINTERACTION **********
    // - 0x00008000: particle coming from reinteraction, direct- or indirect-ly.
    // - 0x00001000: particle having reinteracted
    // (Notes:
    //  - We do this before excluding neutrals so that the inheritance mechanism
    //   below also addresses neutrals' daughters.
    //  - We're only interested in harmful reinteractions, that prevent the
    //   affected track from being reco'd => Let's cut on Z of the track's end
    //   vertex @ 550 cm, i.e. something preserving the track's integrity down
    //   to the last station of straws (which, btw, does not take properly
    //   into account the special case of the VSAT tracks; too bad...).)
    if (!vtx.IsPrimary() && vtx.iTrackOfOrigin()!=-1 &&
	vtx.Pos(2)<550) {
      int imother = vtx.iTrackOfOrigin();
      if (Types[imother]&0x00009000)
	// I) Inheritance (0x00008000): flag reinteraction offsprings
	// II) Flag reinteraction sisters (0x00001000): it speeds up processing
	Types[imct] |= 0x00008000;
      else
	flagReinteraction(e,vtx,imct,imother,Types);
    }
    if (Types[imct]&0x00008000) { // Fill histos of reinteraction yield
      hi_yID->Fill(trk.Pid()); hi_ypT->Fill(trk.LzVec().Perp(lvq.Vect()));
    }

    //#define MC_DEBUG_NEUTRAL_RECOIBILITY
#ifndef MC_DEBUG_NEUTRAL_RECOIBILITY
    if (!trk.Q()) {                              // ***** NOT A NEUTRAL PARTICLE
      if (imct==1) {
	// Case where 1st particle is not a muon
	// This might happen. For it's not a feature of pythia that
	// mu is systematically first. It should have been corrected by Vadim
	// Nevertheless...checking (and also so as to cope w/ files produced
	// prior to Vadim's fix).
	pid1st = trk.Pid();
#  ifdef VERBOSE
 printf("\n * MCInfo: Run %d Evt %d: 1st track = %s. Looking for mu in 2nd..\n",
	   runNum,evNum,trk.Name().c_str());
#  endif
      }
      continue;
    }
#else
    int idebug = 0;
    if (idebug || (int)e.UniqueEvNum()==1048598 && imct==39) {
      int imctm = vtx.iTrackOfOrigin();
      const PaMCtrack *mother = imctm==-1 ? 0 : &e.vMCtrack()[imctm];
      printf("MC%d %s <- %d (%.2f %.2f %.2f) <- Trk %d %s\n",
	     imct,trk.Name().c_str(),
	     trk.iVertexOfOrigin(),vtx.Pos(0),vtx.Pos(1),vtx.Pos(2),
	     imctm,imctm==-1 ? "incident" : mother->Name().c_str());
    }
#endif

    // ********** #HITS per ZONE **********
    HitsPerZone(e,imct,imct /* Do reset */);
    int imct2; if (piK2muEval) {
      // PION IN-FLIGHT DECAY: 
      imct2 = iPiKDecay(e,trk); if (imct2>=0) HitsPerZone(e,imct,imct2);
    }
    else
      imct2 = -1;

    TVector3 v3MC = trk.Mom3();
    double pMC = v3MC.Mag(), thvMC = v3MC.Theta(), phMC = v3MC.Phi();
    double thxMC = thvMC*cos(phMC), thyMC = thvMC*sin(phMC); 
    if (vtx.iTrackOfOrigin()==-1) {           // ***** PRIMARY PARTICLE MCType
      Types[imct] |= 0x0001; if (TrackHistos) hf_pt->Fill(pMC);
    }
#ifdef PSEUDO_pV_VIA_VERTEX
    {                                         // ***** CLOSE to PRIMARY MCType
      double xV = vtx.Pos(0), yV = vtx.Pos(1), zV = vtx.Pos(2);
      // Close in term of vertex resolution => dZ<1.0~=2sigma as measured
      // CG,v7 D0s and dR<.025
      // (Nota: there's a second chance to get those 0x0004 and 0x0008 bits,
      // even for primary particles, see infra.)
      if ((xV-XpV)*(xV-XpV)+(yV-YpV)*(yV-YpV)<.025*.025 &&
	  (zV-ZpV)*(zV-ZpV)<1.0*1.0) Types[imct] |= 0x0004;
      if ((xV-XpV)*(xV-XpV)+(yV-YpV)*(yV-YpV)<.05*.05 &&
	  (zV-ZpV)*(zV-ZpV)<2.0*2.0) Types[imct] |= 0x0008;
    }
#endif
    if (trk.Pid()== 2 || trk.Pid()==3) {      // ***** INTERESTING e^+/-? *****

      if (vtx.iTrackOfOrigin()!=-1) {  // Name() is G3partName[trk.pid]
	const PaMCtrack &mother = e.vMCtrack()[vtx.iTrackOfOrigin()];
	if (mother.Q()==0) {             // Require e^+/- decays from neutral...
	  if (mother.Pid()==1) {         // If gamma...
	    PaMCvertex &grandVtx = e.vMCvertex()[mother.iVertexOfOrigin()];
	    if (grandVtx.iTrackOfOrigin()==-1) // ...Require it to be primary...
	      Types[imct] |= 0x0002;
	    else {
	      const PaMCtrack &grandMother =
		e.vMCtrack()[grandVtx.iTrackOfOrigin()];
	      if (grandMother.Q()==0 &&        // ...OR to decay from neutral...
		  grandMother.Name().find("Gamma"))  // ...HADRON
		Types[imct] |= 0x0002;
	    }
	  }
	  else Types[imct] |= 0x0002;
	}
      }
    }
    else {
      Types[imct] |= 0x0002;                  // ***** !e+/- but prim.... MCType
      if (trk.Name()=="Muon +  " || trk.Name()=="Muon -  ")
	Types[imct] |= 0x00002000; // ***** mu+/-
    }

#ifdef PSEUDO_pV_VIA_VERTEX
    if (!(Types[imct]&0x000f)) continue;
#else
    if (!(Types[imct]&0x0003)) continue;
#endif

    // **********************************************************************
    // ***** => ALL INTERESTING PARTICLES = PRIMARY, pV-gamma decay, !e^+/-
    // **********************************************************************
    if (TrackHistos) hf_at->Fill(pMC);

#ifdef MC_HADRON_SETUP
    bool scatteredParticle = imct == 1;
#else
    bool scatteredParticle = false;
    if (imct==1 || pid1st>=0 && imct==2) {
      if (trk.Pid()!=5 && trk.Pid()!=6) { // mu+=5, mu-=6. Will have to make this an argument of "SortOut"
	if (imct==1) {
#  ifdef VERBOSE
 printf("\n* MCInfo: Run %d Evt %d: 1st track = %s. Looking for mu' in 2nd..\n",
	  runNum,evNum,trk.Name().c_str());
#  endif
	  pid1st = trk.Pid();
	}
	else {
 printf("\n * MCInfo: Run %d Evt %d: 1st track %s 2nd one %s => No mu'\n",
	  runNum,evNum,G3partName[pid1st],trk.Name().c_str());
	}
      }
      else scatteredParticle = true;
    }
#endif
    if (scatteredParticle) {
      iMCT_mu = imct; Pmup = 1/trk.Pinv();

      // ************************************************************
      //          ********** SCATTERED MUON (or otherwise) **********
      // ************************************************************

      const TLorentzVector k = e.vMCtrack()[0].LzVec();
      const TLorentzVector ki(-k.Vect(),k.E()); // Rotate beam track by 2 M_PI
      const TLorentzVector ks = trk.LzVec();
      lvq = ki-ks;
      const TLorentzVector lvp(0,0,0,M_p);
      double pq = lvp.Dot(lvq), pk = lvp.Dot(ki);
      Q2 = -lvq.M2();  xB = Q2/(2*pq); yB = pq/pk; nu = yB*ki.E();
      TLorentzVector lvcm = lvq+lvp; sqrts = lvcm.M();
      betacm = -lvcm.BoostVector(); vqUnit = lvq.Vect().Unit();

      int smu_ible = 0;
      Types[imct] |= 0x00010000; // I.e. flag that it _is_ scattered muon
      // ***** smu RECONSTRUCTIBILITY: #HITS POST... and post MUwALL *****
      if (NHits[2][imct]>=4) { Types[imct] |= 0x00020000; smu_ible |= 1;} // 0x4
      if (NHits[3][imct]>=4) { Types[imct] |= 0x00040000; smu_ible |= 2;} // MWB
      if (NHits[4][imct]>=9) { Types[imct] |= 0x00080000; smu_ible |= 4;} // MWA
      if (TrackHistos) {
	hf_RQ2->Fill(Q2,smu_ible); hf_Ry->Fill(yB,smu_ible);
	hf_RxB->Fill(xB,smu_ible);
	hf_RyvsQ2[0]->Fill(Q2,yB);
	if      (smu_ible&2) hf_RyvsQ2[1]->Fill(Q2,yB); 
	else if (smu_ible&4) hf_RyvsQ2[3]->Fill(Q2,yB); // Exclusive MWB
	if      (smu_ible&4) hf_RyvsQ2[2]->Fill(Q2,yB); // All MWB
      }
    } // End scattered mu
    else {    // ***** MuWALL RECO'IBILITY: BASED #HITS IN MuWALL DETECTORS
      if (NHits[3][imct]>=4) Types[imct] |= 0x00040000; // MWB
      if (NHits[4][imct]>=9) Types[imct] |= 0x00080000; // MWA
    }

    //                 ***** NO HIT WHATSOEVER! *****
    // I guess this corresponds to cases when the particle has reinteracted!?
    // E.g. in "~valexakh/test/pythia.160.full_p-up_56.2002.03.outpipe.fz.1"
    // Evt #2 MCTrack #130 Positron  1.78 GeV/c
    if (NHits[0][imct]==0 && NHits[1][imct]==0 &&
	NHits[2][imct]==0 && NHits[3][imct]==0) {
#ifndef MC_GEOM_ACCEPTANCE
      continue;
#else
      if (!(Types[imct]&0x0001)) continue;
#endif
    }

    // **********************************************************************
    // ********** 1ST and LAST POINTS (i.e. origin and decay vtces) *********
    // **********************************************************************
    double zFirst = vtx.Pos(2);
    // Determine decay vertex = "zLast" (certified by "endID!=-1")
    int endID = -1; static double zLast;            // ***** DECAY VERTEX? *****
    double pi = pMC, pip = pi; const vector<int> &imcvs = trk.vMCvertex();
    int iv, bremsstrahl; for (iv=bremsstrahl = 0; iv<(int)imcvs.size(); iv++) {
      const PaMCvertex &decayV = e.vMCvertex()[imcvs[iv]];
      double zDecay = decayV.Pos(2);
      const vector<int> &decayTrks = decayV.vMCtrack();
      int nDecays = (int)decayTrks.size();
      if (endID!=-1 && zDecay<zLast) {
	printf("** MCInfo:\a Run %d Evt %d MC %d decay#%d: %s",
	       runNum,evNum,imct,iv,trk.Name().c_str());
	for (int jv = 0; jv<(int)imcvs.size(); jv++) {
	  const PaMCvertex &dV = e.vMCvertex()[imcvs[jv]];
	  double zD = dV.Pos(2); printf(" -> %.3f",zD);
	  const vector<int> &dTrks = dV.vMCtrack();
	  for (int jt = 0; jt<(int)dTrks.size(); jt++)
	    printf(" %s",e.vMCtrack()[dTrks[jt]].Name().c_str());
	}
	printf("\n");
      }
      else zLast = zDecay;
      endID = 1; // Flag certifies that "zLast" has been set
      for (int idcy = 0; idcy<nDecays; idcy++) {
	int imctd = decayTrks[idcy]; if (imctd==imct) {
	  // I understand this case could correspond to some soft reinteraction.
	  // And, anyway, the integrity of the track is maintained and no
	  // continuation, w/ a different "imct", is created whose reco'ibility
	  // could be tracked => Could simply update current momentum "pip"
	  // (to be later compared to initial one "pi"). BUT let's assumed this
	  // is not possible
	  printf("** MCInfo:\a Run %d Evt %d MC %d -> %d: %s ->:",
		 runNum,evNum,imct,imct,trk.Name().c_str());
	  for (int j = 0; j<nDecays; j++)
	    printf(" %s",e.vMCtrack()[decayTrks[j]].Name().c_str());
	  printf("\n"); assert(false);
	}
	//else if (pid<=3 && // i.e. e^+ or e^/- (given neutral are disregarded)
	//  	   nDecays==1 && decay.Pid()==pid)
	//  pip = decay.Mom3().Mag();
	// This is another kind of process one can track. Not considered yet.
	//else... In all other cases (reinteraction) one assumes that the
	//       process does not affect the reco'ibility...
      }
      if (nDecays==1) {
	const PaMCtrack &decay = e.vMCtrack()[decayTrks[0]];
	if (decay.Pid()==1) { // Bremsstrahlung
	  // - The distortion this process introduces can be tracked fairly well
	  // - Its probability is quite high.
	  // Therefore it's interesting to try to determine whether this
	  // affects the reco'ibility and remove corresponding tracks.
	  pip -= decay.Mom3().Mag(); bremsstrahl = 1;
	}
      }
      if (!bremsstrahl) continue;
      //#define MC_DEBUG_gamma
#ifdef MC_DEBUG_gamma
      int jdebug = 0;
      if (jdebug)
	printf("%d %s %.2f  -> %.2f (%.1f%%) @ %.f\n",
	       imct,trk.Name().c_str(),pi,pip,pip/pi*100,zDecay);
#endif
      if (looseBremsstrahlCut) {  // ***** LOOSE ``BREMSSTRAHLUNG CUT''... *****
	if (zDecay<zMM02V) { // Bremsstrahlung happens before ENTRANCE to LAS...
	  pi = pip; endID = -1;     // ...reset (initial momentum, flag) 
	  // ...The idea is to let "SortOut" determine the reco'ibility of the
	  // the resulting track piece. Leaving the case of the orginal track
	  // to be decided later, cf. "SortgammaOut"
	}
	if ((zDecay<zSM1 && // Decay upstream of SM1 => momentum is bound to...
	     pip<pi*.33) ||   // ...be biased upon reco => have a stricter cut,
	     pip<pi*.20)    // Else: looser cut... leaving the case of SAS...
	  break;              // ...reco'ibility somewhat wronged
      }
      else {                                              // ***** ...ELSE *****
	if ((zDecay<zSM1 && // Decay upstream of SM1 => momentum is bound to...
	     pip<pi*.90) ||   // ...be biased upon reco => have a strict cut,
	     pip<pi*.75)    // Else: looser cut... leaving the case of SAS...
	  break;              // ...reco'ibility somewhat wronged
      }
      bremsstrahl = 0;  // Reset "bremsstrahl" flag in view of next "iv" loop 
    }  // End loop on daughter vertices

    if (!bremsstrahl && // NO bremsstrahlung OR it has NOT done too much damage
	endID!=-1) {
      // After the last accident (i.e. intervening vertex) the track may
      // in fact have a continuation: determine this on view of its hits list.
      if (zLast<zDC01 &&
	  NHits[0][imct]>=6) zLast = zDC01+.1; // 6 allows to get scifi-tracks
      if (zLast<zST03U1 &&
	  NHits[1][imct]>=8) zLast = zST03U1+.1; // This leaves the case of scifi-tracks unaddressed
      if (zLast<zGM09 &&
	  NHits[2][imct]>=8) zLast = zGM09+.1; // This leaves the case of scifi-tracks unaddressed
    }

    if (trk.Mom3()[2]<0) continue;                 // ***** NO BACKWARD POINTING

    PaTPar helix = trk.ParInVtx();
    if (fabs(helix(5))>1/(.1+.4))   // ***** TOO LOW MOMENTUM Runge-Kutta) *****
      // (cf. "../lib/PaTParExtrap.cc") + margin to keep at a low rate
      // messages "PaAlgo::RKutta ==> Do not getting closer..."
      // + some real margin because we posit that reco'ibility requires it
      continue;
    if (fabs(helix(3))>.5 || fabs(helix(4))>.5)         // ***** TOO LARGE ANGLE
      continue; // +/- 500 mrd: in order to accomodate OD magnet acceptance
			       
    // ***** pVERTEXIBILITY: CDA < .2 (Value derived from RMS of plot of
    // pseudo-CDA (=d(Track's Intercept @ pV - pV)) for MC D*
    helix.Extrapolate(ZpV,heli2,false);
    if (sqrt((heli2(1)-XpV)*(heli2(1)-XpV)+(heli2(2)-YpV)*(heli2(2)-YpV))<.02)
      Types[imct] |= 0x0004;
    else if (Types[imct]&0x0001) {
      printf("\n* MCInfo: Run %d Evt %d: MCTrack %d primary and !vertexible\n",
	     runNum,evNum,imct);
    }
    if (sqrt((heli2(1)-XpV)*(heli2(1)-XpV)+(heli2(2)-YpV)*(heli2(2)-YpV))<.04)
      Types[imct] |= 0x0008;

    // **********************************************************************
    // ***** LAS RECO'IBILITY:
    // - W/IN ACCEPTANCE MM01&DC01&(ST03||FI06)
    // - zFirst < MM02 
    // - #ifndef MC_GEOM_ACC: zLast > ST03 or FI06 for VSATrack's
    // ***** FF RECO'IBILE: zLast<zST03U1
    // **********************************************************************

    if (zFirst>zMM02V)                   // ***** ORIGIN VERTEX < MM02V... *****
      goto SAS_reco;                             // ...i.e. ENTRANCE to LAS

#ifndef MC_GEOM_ACCEPTANCE
    if (endID!=-1 && zLast<zDC01)                  // ***** DECAY VERTEX  > DC01
      // This rejects pions decaying into muons and yield continuous, reco'ible
      // pi-/mu-on tracks to which one can assign a piID in their entirety.
      // Yet, one can argue that those decaying upstream of SM1 will be
      // difficult to reco.
      continue;
#endif

    if (mTarget==2 &&                             // ***** Case of SMC magnet...
	zFirst<zSMCExit) { // Do this only if Origin vertex upstream of its exit
      helix.Extrapolate(zSMCExit,heli2,false);
      if (sqrt(heli2(1)*heli2(1)+heli2(2)*heli2(2))>hwSMC)
	goto SAS_reco;                                      // ***** ...W/IN SMC
    }
    helix.Extrapolate(zAcc,heli2,false);
    if (fabs(heli2(1))>hwAcc || fabs(heli2(2))>hwAcc)
      goto SAS_reco;                                          // ***** W/IN MM01
    heli2.Extrapolate(zDC01,helix,false);
    if (fabs(helix(1))>hwDC01 || fabs(helix(2))>hwDC01)
      goto SAS_reco;                                          // ***** W/IN DC01

    if (endID!=-1 && zLast<zST03U1) { // ********** RECO'IBLE FF... **********
      Types[imct] |= 0x0100;          // ...DECAY VERTEX  < ST03
      if (TrackHistos) {
	hf_fRt->Fill(pMC);
	//                               ********** FF pVERTEXIBILITY **********
	if (Types[imct]&0x0004)
	  hf_pVtF->Fill(pMC); // ***** All Tracks *****
      }
#ifndef MC_GEOM_ACCEPTANCE
      continue;                // ***** ifndef MC_GEOM_ACC: DECAY VERTEX  > ST03
#endif
    }
    helix.Extrapolate(zST03U1,heli2,false);
    if (fabs(heli2(1))>hwxST || fabs(heli2(2))>hwyST)
      goto SAS_reco;                                          // ***** W/IN ST03
#ifndef MC_GEOM_ACCEPTANCE
    if (fabs(heli2(1))<hwFI05 && fabs(heli2(2))<hwFI05) {
      //                     ***** ifndef MC_GEOM_ACC: W/IN VSAT => zLast > FI06
      heli2.Extrapolate(zFI06,helix,false);
      if (fabs(helix(1)-xFI06)<hwFI06 && fabs(helix(2))<hwFI06)
	if (endID!=-1 && zLast<zFI06) goto SAS_reco;
    }
#endif
    // **********************************************************************
    //          ********** RECO'IBLE LAS: FLAGS, HISTOS **********
    // **********************************************************************
    Types[imct] |= 0x00000020;
    if (TrackHistos) {
      hf_Rt->Fill(pMC); hf_lRt->Fill(pMC); // reco'ible in i) any ii) LAS
      //                                    ********** pVERTEXIBILITY **********
      if (Types[imct]&0x0001) {
	hf_pVtP->Fill(pMC); hf_pVaP->Fill(thvMC*1000);// ***** Primaries: P/theta
	haa_pVtP->Fill(thxMC,thyMC); hpa_pVtP->Fill(pMC,thvMC);
      }
      if (Types[imct]&0x0004)
	hf_pVt ->Fill(pMC);                                 // ***** All Tracks
      if (Types[imct]&0x00080000) {
	// Fill denominator for all-out (not restricted to mu) MuWallA efficiency
	hf_jRt->Fill(pMC);
	int pid = trk.Pid(), origin;
	if (pid==5 || pid==6) { // mu+/-
	  if (imct==1) origin = 1;                     // 1:   scattered mu
	  else origin = pid-3;                         // 2/3: mu+/-
	}
	else if (imct2>0) {// "imct2=0" although programmatically allowed is excluded: would mean the particle is a decay of the incident track
	  // Expected 8|9(pi+/-)|11|12(K+/-)
	  if      (pid==8 || pid==9)   origin = pid-4; // 4/5: pi+/-
	  else if (pid==11 || pid==12) origin = pid-5; // 6/7: K+/-
	  else origin = -1;                            // -1: inconsistency
	}
	else { // Could be a gamma conversion
	  origin = 0;                                  // 0: unknown
	  if (vtx.iTrackOfOrigin()!=-1) {
	    const PaMCtrack &mother = e.vMCtrack()[vtx.iTrackOfOrigin()];
	    if (mother.Pid()==1) origin = 8;           // 8: gamma
	  }
	}
	hf_jID->Fill(pMC,origin);
      }
    }

    // **********************************************************************
    // ***** SAS RECO'IBILITY:
    // - W/IN ACCEPTANCE MF1 && HCAL
    // - zFirst < PS01 or FI05 for VSAT
    // - P > 3 GeV
    // - #ifndef MC_GEOM_ACC:
    // **********************************************************************
  SAS_reco:
#ifndef MC_GEOM_ACCEPTANCE
    //                       *****  ifndef MC_GEOM_ACC: NO HIT DOWNSTREAM of SM2
    // This may correspond to cases when the particle has reinteracted!?
    // and which cannot be flagged as such because there's no end vertex
    // E.g. in "~valexakh/test/pythia.160.full_p-up_56.2002.03.outpipe.fz.1"
    // Evt #7 MCTrack #7 Pion +  1.52 GeV/c
    if (NHits[2][imct]==0 && NHits[3][imct]==0) continue;
#endif
    // ***** RECO'IBILITY SAS: W/IN ACCEPTANCE MF1 && HCAL &&
    //   zPS01 or FI05 for VSATrack's < ABSCISSA < ST05(2003) or PA03
    if (fabs(helix(5))>1/(3.))           // ***** TOO LOW MOMENTUM for SM2 *****
      // (cf. Alex in tracking compass note, where > 4 GeV/c)
      continue;

    if (zPS01<zFirst) continue;            // ***** IN ANY CASE: zFirst < PS01...
#ifndef MC_GEOM_ACCEPTANCE
    if (endID!=-1 && zLast<zGM09)  // ***** ifndef MC_GEOM_ACC: ... GM09 < zLast
      continue;
#endif
    helix.Extrapolate(zFI05,heli2,false);
    if (fabs(heli2(1))<hwFI05 && fabs(heli2(2))<hwFI05
	&& zFirst>zFI05)                  // ***** W/IN VSAT => REQUIRE z < FI05
      continue;
    helix.Extrapolate(zPS01,heli2,false);
    if (fabs(heli2(1))>hwxST || fabs(heli2(2))>hwyST) // ***** IN ANY CASE...
      continue;                 // ***** LOOSE ENTRANCE ACCEPTANCE CUT (@ ~PS01)
    heli2.Extrapolate(zHC1,helix,false);
    if (fabs(helix(1))>hwxHC1 || fabs(helix(2))>hwyHC1)
      continue;                                 // ***** W/IN HCAL1 CENTRAL HOLE
    helix.Extrapolate(zMF1,heli2,false);
    if (fabs(heli2(1))>hwxMF1 || fabs(heli2(2))>hwyMF1)
      continue;                                  // ***** W/IN muF1 CENTRAL HOLE
    heli2.Extrapolate(zPA05,helix,false);
    if (fabs(helix(1))>hwxPA || fabs(helix(2))>hwyPA)
      continue;                                               // ***** W/IN PA05

    Types[imct] |= 0x00000040;  // ********** RECO'IBLE SAS **********
    if (TrackHistos) {                          // Reco'ible in L+SAS...
      if (Types[imct]&0x00000020) hf_aRt->Fill(pMC);// ...both
      else                        hf_Rt->Fill(pMC); // ...either one (Fill it once!)
      hf_sRt->Fill(pMC);                        // Reco'ible in SAS
    }
  }
  // ***** SUMMARY OF REINTERACTION
  // - Histogram the ID of all primary tracks and of reInteracting ones
  // - Set MCTrack -> Track to assign to the reinteracting, mother, MCTrack
  //                       the reco daughter 
  // - Set   Track -> MCTrack to assign to the reco daughter the mother MCTrack
  //   In the process the reco'ibility of the daughter MCTrack is transfered to
  //  its mother.
  for (imct = 1; imct<e.NMCtrack(); imct++) {
    PaMCtrack &trk = e.vMCtrack()[imct]; int &Type = Types[imct];
    int ivtxm = trk.iVertexOfOrigin();
    const PaMCvertex &vtx = e.vMCvertex()[ivtxm];
    if (!vtx.IsPrimary()) continue;
    int pid = trk.Pid();
    double Px = trk.P(0), Py = trk.P(1), Pz = trk.P(2), P = sqrt(Px*Px+Py*Py+Pz*Pz);
    hi_pID->Fill(pid,P); if (Type&0x00001000) hi_iID->Fill(pid,P);
    if (vtx.Pos(2)>setup.TargetCenterZ()) { // Independently, most downstream pVs
      hi_pId->Fill(pid,P); if (Type&0x00001000) hi_iId->Fill(pid,P);
    }
    if (Type&0x00001000) { // Reevaluate interacting mother...
      const PaMCvertex &dVtx = e.vMCvertex(trk.vMCvertex().front());
      const vector<int> &idaughters = dVtx.vMCtrack();
      int ndaughters = idaughters.size();
      double E = trk.E();
      for (int i = 0; i<ndaughters; i++) {
	int idaughter = idaughters[i];
	const PaMCtrack &daughter = e.vMCtrack()[idaughter];
	if (daughter.Pid()!=pid) continue;
	double Ep = daughter.E(), dE = E-Ep; if (dE<E*.05) {
	  //#define MC_DEBUG_ReInter
#ifdef MC_DEBUG_ReInter
	  printf(" * MCInfo::SortOut: Evt#%d MC %d(%d,%.1fGeV) -> %d(%.1fGeV)\n",
		 evNum,imct,pid,E,idaughter,Ep);
#endif
	  // ...if soft reinteraction (dE<3%) => undo the re-interaction flag
	  Type ^= 0x00001000;
	  // ...and try and update MC<->reco X-references
	  if (trk.iBestTrkRef()<0) {
	    int iref = daughter.iBestTrkRef(); if (iref>0) {
	      trk.SetBestTrkRef(iref);                  // MCtrack -> BEST Track
	      PaTrack &t = const_cast<PaTrack&>(e.vTrack(iref));
	      t.SetMCtrack(imct);                       // Track   -> MCtrack
	      // Transfer reco'ibility of daughter to mother
	      // Transfer also the compatibility w/ vertex
	      Type |= (Types[idaughter]&(0x0060|0x0004));
	    }
	  }
	  else
	    // - The reinteracing track is reco'd before the, late, interaction.
	    // - Since we build, cf. supra, BestTrkRef requiring reco'd momentum,
	    //  among other things, the reco'd track is already a pretty good
	    //  reco of the MC.
	    // => Do no update
	    printf("** MCInfo::SortOut: Evt#%d Reinteracting MC %d(%d,%.1fGeV)"
		   " -> %d(%.1fGeV) already reco'd\n",
		   evNum,imct,pid,E,idaughter,Ep);
	}
      }
      if (Type&0x00001000) { // If hard reinteraction => histogram.
	hi_IID->Fill(pid,P);
	if (vtx.Pos(2)>setup.TargetCenterZ()) hi_IId->Fill(pid,P);
      }
    }
  }
  // ***** WARNING: IF NO MUON FOUND
  if (iMCT_mu==-1)
 printf("\n** MCInfo:\a Run %d Evt %d: mu not found\n",runNum,evNum);
}
void MCInfo::HitsPerZone(PaEvent &e, int imct, int imct2)
{
  // - Standard is "imct2==imct".
  // - "imct2!=imct" correspond to the pi(K)->mu case.
  //   - pi = "imct",
  //   - mu = "imct2",
  //   - "NHits[*][imct]" is updated (and not reset) w/ the hits of ("imct2").
  PaMCtrack  &trk = e.vMCtrack()[imct2];
  int iw, idet, izone;
  if (imct2==imct) for (izone = 0; izone<5; izone++) NHits[izone][imct] = 0;
  for (iw = 0, idet = 0, izone = -1; iw<nHitMapWs; iw++) {
    for (int ib = 0; ib<32; ib++, idet++) {
      if (idet>iZoneLasts[izone+1]) {
	izone++;
	// Hence "izone=0" means target<->SM1, cf. def. of "iZoneLasts" 
      }
      if (izone==4) break;
      if (izone>=0 && (1<<ib&trk.HitMapArray()[iw])) {
	NHits[izone][imct]++;
	if (izone==1 && FirstMuWallA<=idet && idet<=LastMuWallA &&
	    (1<<idet-FirstMuWallA&PlPatMuWallA))
	  NHits[4][imct]++;
      }
    }
    if (izone==4) break;
  }
}
int iPiKDecay(PaEvent &e, PaMCtrack &trk)
{
  // Criteria:
  // - pi+/-
  // - First decay vertex is pi+/- -> mu+/- + one other MC track
  // - mu's MC track does not differ too much from pi's in momentum.
  // Return the "imct" of the daughter muon MC track.
  int C; if (trk.Pid()==8 || trk.Pid()==11) C = 0;  // pi+ or K+
  else if   (trk.Pid()==9 || trk.Pid()==12) C = 1;  // pi- or K-
  else return -1;
  const vector<int> &imcvs = trk.vMCvertex(); if (imcvs.empty()) return -1;
  const PaMCvertex &decayV = e.vMCvertex()[imcvs[0]];
  const vector<int> &decayTrks = decayV.vMCtrack();
  int nDecays = decayTrks.size();
  if (nDecays<=2) { // A priori, we would like "nDecays==2", since we're after
    // pi(K) -> mu+nu. But it looks like nu are bypassed
    for (int idcy = 0; idcy<nDecays; idcy++) {
      int imct = decayTrks[idcy];
      const PaMCtrack &decay = e.vMCtrack()[imct];
      if (decay.Pid()==5+C/* mu+ or mu- */) {
	TLorentzVector pi_K = trk.LzVec(), mu = decay.LzVec();
	mu.Boost(-pi_K.BoostVector()); double angle = mu.Angle(pi_K.Vect());
	if (angle<90) return imct;
	else          return -1;
      }
    }
  }
  return -1;    
}
void flagReinteraction(PaEvent &e,
		       PaMCvertex vtx,        // Vertex of origin
		       int imct, int imother, // Particles's index and mother's
		       int *Types)
{
  // ********** REINTERACTION **********
  // - Assigning 0x00008000 to MCType (and 0x00001000 to mother) 
  // - Simplifications:
  //   - Only pi/K+/-, K0S, proton, neutron, Lambda and Sigma+/- reinteractions.
  //    And offsprings, cf. step preceding the call to "flagReinteraction".
  //    (This in order to exclude em showers. And muons being left aside because
  //    weakly (re)interacting.)
  //   - Reject decays by checking whether mother and sisters are consistent w/
  //    a decay branch w/ a significant branching ratio:
  //    - pi+/- -> mu+/- (+ neutrino)
  //    - K+/-:
  //      - 1->1: e/mu+/- (+ neutrino)
  //      - 1->2: pi0 + e/mu+/- (+neutrino)
  //              pi+/- + pi0
  //      - 1->3: pi+/- + 2pi
  //    - Lambda: p+pi, n+pi0
  //      Sigma+: p+p0, n+pi-
  //      Sigma-: n+pi-, 
  //   - These pi/K, etc.. are required to come from primary or from the decay of
  //    main resonances (rho/phi/omega/K*)
  //   - MC particles are assumed to be handed over in lineage order (which is
  //    reasonable).
  const PaMCtrack &mother = e.vMCtrack()[imother];
  int motherId = mother.Pid();
  bool motherpi =    motherId==8  || motherId==9;  // pi+/-
  bool motherK =     motherId==11 || motherId==12; // K+/-
  bool motherN =     13<=motherId && motherId<=15; // n,p,A-p
  bool motherK0 =    motherId==16;                 // K0S
  bool motherHyper = motherId==18 || motherId==19; // Lambda/Sigma+
  motherHyper |=     motherId==26 || motherId==27; // C.C.
  motherHyper |=     motherId==21 || motherId==29; // Sigma-
  if (motherpi || motherK || motherN || motherK0 || motherHyper ||
      motherId==25 /* A-neutron */) {
    const vector<int> &isisters = vtx.vMCtrack(); int nsisters = isisters.size();
    if (!nsisters) {
      printf("\n **MCInfo:\a Run %d Evt %d: Part %d(%s) inconsistent decay\n",
	     e.RunNum(),(int)e.UniqueEvNum(),imother,mother.Name().c_str());
      assert(false);
    }
    int id1 = e.vMCtrack()[isisters[0]].Pid(), decay = 0;
    if (motherK) {                     // Special case of K+/-
      if      (nsisters==1) {     // Reject lepto decays: e/mu+/- (+ neutrino)
	if (id1==motherId-8 || id1==motherId-6) decay = 1;
      }
      else if (nsisters==2) {     // Reject pi+/-pi0 decays
	int id2 = e.vMCtrack()[isisters[1]].Pid();
	if      (id1==motherId-3 && id2==7 ||
	    id2==motherId-3 && id1==7) decay = 2;
	//                           Reject semilepto: pi0 + e/mu+/- (+ neutrino)
	else if (id1==7 && (id2==motherId-8 || id2==motherId-6) ||
	    id2==7 && (id1==motherId-8 || id1==motherId-6))
	  decay = 4;
      }
      else if (nsisters==3) {       // Reject 3pi decays
	int id2 = e.vMCtrack()[isisters[1]].Pid(),
	  id3 = e.vMCtrack()[isisters[2]].Pid();
	int sameSign = motherId-3, cc = 17-sameSign;
	if ((id1==sameSign && (id2==sameSign && id3==cc ||
			       id3==sameSign && id2==cc || id2==7 && id3==7)) ||
	    (id2==sameSign && (id3==sameSign && id1==cc ||
			       id1==sameSign && id3==cc || id3==7 && id1==7)) ||
	    (id3==sameSign && (id1==sameSign && id2==cc ||
			       id2==sameSign && id1==cc || id1==7 && id2==7)))
	  decay = 3;
      }
    }
    else if (motherK0) {               // Special case of K0S
       if (nsisters==2) {
	int id2 = e.vMCtrack()[isisters[1]].Pid();
	if (id1==8 && id2==9 ||
	    id2==8 && id1==9 || id1==7 && id2==7) decay = 5;
       }
    }
    else if (motherpi) {               // Special case of pi+/-
      if (nsisters==1) {            // Reject lepto decays: e|mu+/- + neutrino
	if (id1==motherId-3 || id1==motherId-5) decay = 10;
      }      
    }
    else if (motherId==13) {               // Special case of n
      if (nsisters==2) {            // Reject p + e-
	int id2 = e.vMCtrack()[isisters[1]].Pid();
	if (id1==14 && id2==3 ||
	    id2==14 && id1==3) decay = 11;
      }      
    }
    else if (motherId==25) {               // Special case anti-n
      if (nsisters==2) {            // Reject p- + e+
	int id2 = e.vMCtrack()[isisters[1]].Pid();
	if (id1==15 && id2==2 ||
	    id2==15 && id1==2) decay = 11;
      }      
    }
    else if (motherId==18 || motherId==26) {// Special case of Lambda
      if (nsisters==2) {
	int id2 = e.vMCtrack()[isisters[1]].Pid(), cc = motherId==18?0:1;
	if (id1==14+cc && id2==9-cc || id1==13+cc*12 && id2==7 ||
	    id2==14+cc && id1==9-cc || id2==13+cc*12 && id1==7) decay = 20;
      }
    }
    else if (motherId==19 || motherId==27) {// Special case of Sigma+
      if (nsisters==2) {
	int id2 = e.vMCtrack()[isisters[1]].Pid(), cc = motherId==19?0:1;
	if (id1==14+cc && id2==7 || id1==13+cc*12 && id2==8+cc ||
	    id2==14+cc && id1==7 || id2==13+cc*12 && id1==8+cc) decay = 21;
      }
    }
    else if (motherId==21 || motherId==29) {// Special case of Sigma-
      if (nsisters==2) {
	int id2 = e.vMCtrack()[isisters[1]].Pid(), cc = motherId==21?0:1;
	if (id1==13+cc*12 && id2==9-cc ||
	    id2==13+cc*12 && id2==9-cc) decay = 22;
      }
    }
    if (!decay) {
      const PaMCvertex &grandVtx = e.vMCvertex()[mother.iVertexOfOrigin()];
      if (grandVtx.IsPrimary()) {
	Types[imct] |= 0x00008000;    // Current particle flagged: product of reinteraction
	Types[imother] |= 0x00001000; // Mother flagged: has reinteracted
      }
      else {
	// Not clear what we have to do. Typically we get here in case of
	// grandma = K0S or Lambda and pi or pi daughter reinteracts.
	// Do we want to include them in the reinteraction? Keeping in mind
	// that we are after
	// - primary hadron loss via reinteraction,
	// - undue association to the pV of products of reinteractions,
	// and that we want to histogram this (cf. "hi_" histos).
	// For the time being, let's put a cut of the decay dist (of K0S,
	// Lambda,... or whatever) and primary vertex.
	int igrandma = grandVtx.iTrackOfOrigin();
	if (igrandma!=-1) {
	  const PaMCtrack &grandma = e.vMCtrack()[igrandma];
	  //int grandmaId = grandma.Pid();
	  PaMCvertex &gdgdVtx = e.vMCvertex()[grandma.iVertexOfOrigin()];
	  if (gdgdVtx.IsPrimary()) {
	    double X = gdgdVtx. Pos(0), Y = gdgdVtx. Pos(1), Z = gdgdVtx. Pos(2);
	    double x = grandVtx.Pos(0), y = grandVtx.Pos(1), z = grandVtx.Pos(2);
	    double d = sqrt((X-x)*(X-x)+(Y-y)*(Y-y)+(Z-z)*(Z-z));
	    if (d<1) {
	      Types[imct] |= 0x00008000;    // Current particle flagged: product of reinteraction
	      Types[imother] |= 0x00001000; // Mother flagged: has reinteracted
	    }
	  }
	}
      }
    } // End is not K decay
  } // End piKpn mother
}
 
void MCInfo::BookD0Efficiency(double cLowS, double cUpS, double zS, // D* select
			      double cLow0, double cUp0, double z0, // D0 select
			      double yLow, double yUp,
			      double R2, // If 2016/17, apply PaAlgo::InTarget in addition
			      bool richAngleCut, double richCut)
{
  // ********** BOOK ``DENOMINATOR'' HISTOGRAMS **********
  // + VARIOUS INITIALISATIONS

  cLowCutS = cLowS; cUpCutS = cUpS; zCutS = zS;  // ***** KINEMATICAL CUTS *****
  cLowCut0 = cLow0; cUpCut0 = cUp0; zCut0 = z0;
  yLowCut = yLow; yUpCut = yUp; R2Cut = R2;
  RICHAngleCut = richAngleCut; RICHCut = richCut;
#ifdef MC_DEBUG_RECOIBILITY
  ribilityFilled = 0;
#endif	

  // ***** D0 -> Kpi: ACCEPTANCE, RECONSTRUCTIBILITY *****
  char title[] =
    "Reco'ible D^{o} vs. c#theta*  -0.85<c#theta*>0.85,z>0.25,0.35<y<0.85,RICH        ";
  char RRange[] = "A>20,";
  if (RICHAngleCut) sprintf(RRange,"A>%.0f,",RICHCut*1000);
  else              sprintf(RRange,"R>%.0f,",RICHCut);
  for (int iDKcz = 0; iDKcz<4; iDKcz++) {  // Vs. pD, pK, cos(theta*) and zD
    // Successively all D0s
    //              reco'ible
    //              reco'ible RICH'able
    // For all D0s
    //     D0s w/in cth* cut
    //     D0s w/in cth* cut and z cut
    //     D0s w/in all kin cuts
    char hName[] = "hP_RD0ki"; char Pp, vs[] = "c#theta*"; double xMin, xMax;
    switch (iDKcz) {
    case 0: Pp = 'P'; sprintf(vs,"pD"); xMin = 0; xMax = 100; break;
    case 1: Pp = 'p'; sprintf(vs,"pK"); xMin = 0; xMax = 100; break;
    case 2: Pp = 'c'; sprintf(vs,"c#theta*"); xMin = -1; xMax = 1; break;
    default:
    case 3: Pp = 'z'; sprintf(vs,"zD"); xMin = 0; xMax = 1;
    }
    sprintf(hName,"h%c_aD0",Pp);
    sprintf(title,"All D^{o} vs. %s",vs);
    hp_aD0[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    sprintf(hName,"h%c_RD0",Pp);
    sprintf(title,"Reco'ible D^{o} vs. %s",vs);
    hp_RD0[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    // KIN CUTS: theta*
    sprintf(hName,"h%c_aD0c",Pp);
    sprintf(title,"All D^{o} vs. %s - %4.2f<c#theta*<%3.2f",       vs,cLowCut0,cUpCut0);
    hp_aD0c[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    sprintf(hName,"h%c_RD0c",Pp);
    sprintf(title,"Reco'ible D^{o} vs. %s - %4.2f<c#theta*<%3.2f", vs,cLowCut0,cUpCut0);
    hp_RD0c[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    // KIN CUTS: theta* + z
    sprintf(hName,"h%c_aD0cz",Pp);
    sprintf(title,"All D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f",
	    vs,cLowCut0,cUpCut0,zCut0);
    hp_aD0cz[iDKcz] = new TH1D(hName,title,20,xMin,xMax);
    // All KIN CUTS: theta*,z + y
    sprintf(hName,"h%c_aD0k",Pp);
    sprintf(title,
	    "All D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut);
    hp_aD0k[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    sprintf(hName,"h%c_RD0k",Pp);
    sprintf(title,"Reco'ible D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut);
    hp_RD0k[iDKcz]  = new TH1D(hName,title,20,xMin,xMax);
    // All KIN CUTS, AND RICH'able
    sprintf(hName,"h%c_aD0ki",Pp);
    sprintf(title,"All D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,RICH",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut);
    hp_aD0ki[iDKcz] = new TH1D(hName,title,20,xMin,xMax);
    sprintf(hName,"h%c_RD0ki",Pp);
    sprintf(title,"Reco'ible D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,RICH",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut);
    hp_RD0ki[iDKcz] = new TH1D(hName,title,20,xMin,xMax);
    // All KIN CUTS, AND RICH'able and w/in angle or radius cut
    sprintf(hName,"h%c_aD0kI",Pp);
    sprintf(title,"All D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,%s",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut,RRange);
    hp_aD0kI[iDKcz] = new TH1D(hName,title,20,xMin,xMax);
    sprintf(hName,"h%c_RD0kI",Pp);
    sprintf(title,"Reco'ible D^{o} vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,%s",
	    vs,cLowCut0,cUpCut0,zCut0,yLowCut,yUpCut,RRange);
    hp_RD0kI[iDKcz] = new TH1D(hName,title,20,xMin,xMax);
  }

  // ***** D* -> D0pi: ACCEPTANCE, RECONSTRUCTIBILITY *****
  char titl2[] =
    "Reco'ible D^{o} vs. c#theta* - -0.85<c#theta*<.085,z>0.25,0.35<y<0.85,RICH       ";
  // Reco'ible D^{o} vs. c#theta* - -0.85<c#theta*<0.90,z>0.25,0.35<y<0.85,R>12,
 for (int iDKcz = 0; iDKcz<4; iDKcz++) {  // Vs. pD, pK, cos(theta*) and zD
    // Successively all D0s
    //              reco'ible
    //              reco'ible RICH'able
    // For all D0s
    //     D0s w/in all kin cuts
   char hName[] = "hP_RDSki"; char Pp, vs[] = "c#theta*"; double xMin, xMax;
    switch (iDKcz) {
    case 0: Pp = 'P'; sprintf(vs,"pD"); xMin = 0; xMax = 100; break;
    case 1: Pp = 'p'; sprintf(vs,"pK"); xMin = 0; xMax = 100; break;
    case 2: Pp = 'c'; sprintf(vs,"c#theta*"); xMin = -1; xMax = 1; break;
    default:
    case 3: Pp = 'z'; sprintf(vs,"zD"); xMin = 0; xMax = 1; 
    }
    sprintf(hName,"h%c_aDS",Pp);
    sprintf(titl2,"All D* vs. %s",vs);
    hp_aDS[iDKcz]   = new TH1D(hName,titl2,20,xMin,xMax);
    sprintf(hName,"h%c_RDSk",Pp);
    sprintf(titl2,"Reco'ible D* vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f",
	    vs,cLowCutS,cUpCutS,zCutS,yLowCut,yUpCut);
    hp_RDSk[iDKcz]  = new TH1D(hName,titl2,20,xMin,xMax);
    sprintf(hName,"h%c_RDSki",Pp);
    sprintf(titl2,"Reco'ible D* vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,RICH",
	    vs,cLowCutS,cUpCutS,zCutS,yLowCut,yUpCut);
    hp_RDSki[iDKcz] = new TH1D(hName,titl2,20,xMin,xMax);
    sprintf(hName,"h%c_RDSkI",Pp);
    sprintf(titl2,"Reco'ible D* vs. %s - %4.2f<c#theta*<%3.2f,z>%3.2f,%3.2f<y<%3.2f,%s",
	    vs,cLowCutS,cUpCutS,zCutS,yLowCut,yUpCut,RRange);
    hp_RDSkI[iDKcz] = new TH1D(hName,titl2,20,xMin,xMax);
  }

  // ***** KINEMATICAL CUT *****
  hD_cs    = new TH1D("hD_cs"      ,"c#theta*  -  K#pi=D^{o}",
		      200,-1,1);
}
// ***************************************************************************
// ******************************   SortD0Out   ******************************
// ***************************************************************************
void MCInfo::SortD0Out(PaEvent &e, bool richAcceptance)
{
  // ********** PROCESS MC RAW INFO FOR D0 & D* RECONSTRUCTIBILITIES **********

  // D0 -> Kpi  DS -> Kpipi

  // Expecting "SortOut" to have been already executed


  for (int i = 0; i<2; i++) {
    iMCT_D0[i]=iMCD0_K[i]=iMCD0_pi[i] = -1;
    iMCT_DS[i]=iMCDS_D0[i]=iMCDS_piS[i] = -1;
    D0Decay[i] = 0; D0Type[i]=DSType[i] = 0;
  }

  int runNum = e.RunNum(), evNum = (int)e.UniqueEvNum();

  if (richAcceptance) {

    // *************** RICH ACCEPTANCE && TABULATED PID for ***************
    // ***************       PIONS and KAONS                ***************

    // Will have to extrapolate to RICH and determine whether track w/in its
    // overall acceptance and w/in its special angle or radius cut
    static PaTPar heli2;
    // ***** RICH ACCEPATNCE: GEOMETRY INFO *****
    const double zSM1 = 350;
    const double zRICH = 615.6, rRICH = 5, r2RICH = rRICH*rRICH,
      XRICH = 161, YRICH = 121;
    // ***** TABULATED EFFICIENCY (AND misID EFFICIENCY) *****
    // - RICH perfs can be taken from MC simulation of the RICH response, or
    //  derived from a table of the real perfs as a function of phase space (as
    //  measured e.g. by the K0/phi method).
#define RICH_EFFICIENCY_VS_PS_2 1
#if   RICH_EFFICIENCY_VS_PS_2 == 1
    // Ideal
    const double eRICHps[9] = {1,1,1,1,1,1,1,1,1};
    const double mRICHps[9] = {0,0,0,0,0,0,0,0,0};
#elif RICH_EFFICIENCY_VS_PS_2 == 2
    // 2sigmas veto 3sigmas
    const double eRICHps[9] = {.243,.285,.078,
			       .358,.384,.113,
			       .339,.609,.190};  // for >0
    // ``Efficiency'' of misID
    const double mRICHps[9] = {.046,.081,.076,
			       .039,.042,.030,
			       .016,.022,.026};
#endif
#if RICH_EFFICIENCY_VS_PS_2 < 3
    static double pRICHps[9];  // Purity
    static bool first = true; if (first) {
      first = false;
      for (int i = 0; i<9; i++) {
	// Determine purity from efficiency and misID, assuming K/pi = 1/10
	double e = eRICHps[i]; pRICHps[i] = e/(e+10*mRICHps[i]);
      }
    }
#endif

    int imct; for (imct = 1; imct<e.NMCtrack(); imct++) {

      // ********** LOOP ON MC TRACKS, Excluding incident muon **********

      int &Type = Types[imct];
      if (!(Type&0x0003)) continue;      // ***** PRIMARY, pV-gamma decay, !e^+-
      if (NHits[0][imct]==0 && NHits[1][imct]==0 &&// ***** NO HIT WHATSOEVER...
	  NHits[2][imct]==0 && NHits[3][imct]==0) continue;     // ...=> GIVE-UP

      const PaMCtrack  &trk = e.vMCtrack()[imct];
      bool ispi = trk.Pid()==8  || trk.Pid()==9;
      bool isK  = trk.Pid()==11 || trk.Pid()==12;
      if (!ispi && !isK) continue;                         // ***** pi's AND K's

      const PaMCvertex &vtx = e.vMCvertex()[trk.iVertexOfOrigin()];
      double pMC = trk.Mom3().Mag();
      const double KThr = 10;
      if (pMC<=KThr) continue;  // Have to do something for pi's

      // ********** 1ST and LAST POINTS (i.e. origin and decay vtces) **********
      double zFirst = vtx.Pos(2);
      int imcv, endID; static double zLast;     // ***** DECAY VERTEX? *****
      int nvs = (int)trk.vMCvertex().size();
      for (imcv = 0, endID = -1; imcv<nvs; imcv++) {
	PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex()[imcv]];
	if (endID== -1 || decayV.Pos(2)>zLast) {
	  endID = imcv; zLast = decayV.Pos(2);
	}
      }
      if (zSM1<zFirst ||                              // ***** W/IN RICH SPAN...
	  (endID!= -1 && zLast<zRICH &&
	  // ...Yet pions allowed to decay upstream, because their decay into
	  // muons yields a continuous, reco'ible pi-/mu-on track, which RICH
	  // signature is almost identical to that of a pure pion track.
	   !ispi))
	continue;

      //         ********** EXTRAPOLATE to RICH **********
      if (trk.Mom3()[2]<0) continue;  // No backward pointing
      PaTPar helix = trk.ParInVtx();
      if (fabs(helix(5))>1/(.1+.05))  // Too low a momentum for Runge-Kutta
	continue;
      if (fabs(helix(3))>.5 || fabs(helix(4))>.5) // Too large angle
	continue;
      helix.Extrapolate(zRICH,heli2,false);

      if (heli2(1)*heli2(1)+heli2(2)*heli2(2)>r2RICH &&
	  fabs(heli2(1))<XRICH && fabs(heli2(2))<YRICH) {
	//                                         
	Type |= 0x10000000; // ***** RICH'able =  W/IN RICH RANGE and ACCEPTANCE
	if (heli2(1)*heli2(1)+heli2(2)*heli2(2)>RICHCut)
	  Type |= 0x20000000;                           // ***** W/IN INNER RICH
      }
      else continue;
      // Note: we do give up at this point. It may not be fully correct, in so
      // far as tabulated RICH perfs may have been determined w/o a RICH
      // acceptance cut and hence may include, although not in a fully
      // satisfactory way (they're derived from K0 and and phi decays which have
      // a definite phase space distribution), the RICH inefficiency due to
      // tracks outside acceptance.

      //          ********** PHASE SPACE BINNING **********
#if RICH_EFFICIENCY_VS_PS_2 < 3
      double A = acos(1/sqrt(1+heli2(3)*heli2(3)+heli2(4)*heli2(4)));
      int ipa = 0;
      if      (pMC-KThr>20) ipa += 2;
      else if (pMC-KThr>10) ipa += 1;
      if      (A>.024) ipa += 6;
      else if (A>.012) ipa += 3;

      if      (trk.Pid()==8 || trk.Pid()==9) {
	if (drand48()<=mRICHps[ipa]) Type |= 0x80000000;
      }
      else if (trk.Pid()==11 || trk.Pid()==12) {
	if (drand48()<=eRICHps[ipa]) Type |= 0x40000000;
      }
#else          // ***** RICH EFFICIENCY from PERFS TABULATED as a f(PHASE SPACE)
      double A = acos(1/sqrt(1+helix(3)*helix(3)+helix(4)*helix(4)));
      int ia, ja; for (ja = 0, ia = -1; ja<nRICHaBins; ja++)
		    if (A<RICHaBins[ja]) { ia = ja; break; }
      if (ia<0) continue;//ia = nRICHaBins-1;
      double pTable; // Either P or P-Thr, depending upon "RICHpLow=0" or not.
      if (RICHpLow) pTable = pMC;
      else          pTable = ispi ? pMC-piThr : pMC-KThr;
      if (pMC>RICHpLow) { // If P not lower than lowest tabulated bin
	int ip, jp; for (jp = 0, ip = -1; jp<nRICHpBins; jp++)
		      if (pTable<RICHpBins[jp]) { ip = jp; break; }
	if (ip<0) continue;//ip = nRICHpBins-1;

	if      (trk.Pid()==11) {
	  if (drand48()<eRICHKp[ip][ia]) Type |= 0x40000000;
	}
	else if (trk.Pid()==12) {
	  if (drand48()<eRICHKm[ip][ia]) Type |= 0x40000000;
	}
	else if (trk.Pid()==8) {
	  if (drand48()<eRICHpip[ip][ia]) Type |= 0x80000000;
	}
	else if (trk.Pid()==9) {
	  if (drand48()<eRICHpim[ip][ia]) Type |= 0x80000000;
	}
      }
#endif
    }  // End loop on MC tracks

  }  // ********** END DETERMINE "richAcceptance" **********

  for (int iter = 0; iter<2; iter++) {

    // ********** 1ST ITER: D0->K-pi+ **********
    // ********** 2ND ITER: D0*->D0pi+ (using results of 1st iter) **********

    int imct = 0, nD0s[2] = {0,0};
    imct++;                                       // ***** EXCLUDE INCIDENT MUON
    for (;  imct<e.NMCtrack(); imct++) {

      int &Type = Types[imct];
      if (Type&0x00004000) continue;                       // ***** NOT A PILEUP
      const PaMCtrack  &trk = e.vMCtrack()[imct];

      if (iter==0 && (trk.Pid()==37 || trk.Pid()==38)) {        // ***** D0, aD0
	//         ******************************************
	//         ********** 1ST ITER: D0/anti-D0 **********
	//         ******************************************
	int C = trk.Pid()==37 ? 1 : 0; // D0: C = 1! So that C = 0 means K+
	//#define NO_aD0
#ifdef NO_aD0
	if (C==0) continue;
#endif
	if (nD0s[C])
 printf("\n** MCInfo:\a Run#%d Evt %d: >1 %s (MCs %d,%d)\n",
	  runNum,evNum,C?"D0":"aD0",iMCT_D0[C],imct);
	nD0s[C]++;
	iMCT_D0[C] = imct;
	if (trk.vMCvertex().size()>1) {
 printf("** MCInfo:\a Run#%d Evt %d MC %d: D0 w/ >1 decay vertices\n",
	  runNum,evNum,imct); assert(false);
	}
	const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
	const vector<int> &decayTrks = decayV.vMCtrack();
	//#define MCInfo_ALLOW_Kpipi0
#ifdef MCInfo_ALLOW_Kpipi0
	int imct1, imct2, branch;
	if      (decayTrks.size()==2) {
	  imct1 = decayTrks[0]; imct2 = decayTrks[1];
	  branch = 0x1; // Kpi not yet sure, but non Kpi 2-body later rejected
	}
	else if (decayTrks.size()==3) {          // ***** 3 DAUGHTERS, 1 NEUTRAL
	  int idecay, nChargedDecays, imct3, *imcts[3] = {&imct1,&imct2,&imct3};
	  for (idecay=nChargedDecays = 0; idecay<3; idecay++) {
	    int imcti = decayTrks[idecay];
	    if (e.vMCtrack()[imcti].Q()) *imcts[nChargedDecays++] = imcti;
	  }
	  if      (nChargedDecays==3) {
 printf("** MCInfo:\a Run %d Evt %d MC %d: D0 w/ 3 decays\n",
	    runNum,evNum,imct); continue;
	  }
	  else if (nChargedDecays<2) continue;
	  branch = 0x2;
	}
	else continue;
#else
	if (decayTrks.size()!=2) continue;                  // ***** 2 DAUGHTERS
	int imct1 = decayTrks[0], imct2 = decayTrks[1];
#endif
	const PaMCtrack &decay1 = e.vMCtrack()[imct1];
	const PaMCtrack &decay2 = e.vMCtrack()[imct2];
	if      (decay1.Pid()==9-C && decay2.Pid()==11+C) {      // ***** 1: piK
	  iMCD0_pi[C] = imct1; iMCD0_K[C] = imct2;
	}
	else if (decay2.Pid()==9-C && decay1.Pid()==11+C) {      // ***** 2: Kpi
	  iMCD0_pi[C] = imct2; iMCD0_K[C] = imct1;
	}
	else continue;
#ifdef MCInfo_ALLOW_Kpipi0
	D0Decay[C] = branch;
#else
	D0Decay[C] = 0x1;
#endif
	D0EffHistosAndFlags(e,C);
      }  // End D0 block
      else if (iter==1 && (trk.Pid()==55 || trk.Pid()==56)) {     // ***** D*+/-
	//         ******************************************
	//         ********** 2ND ITER: D*->D0pi **********
	//         ******************************************
	int C = trk.Pid()==55 ? 1 : 0;
	iMCT_DS[C] = imct;
	if (trk.vMCvertex().size()>1) {
 printf("** MCInfo:\a Run %d Evt %d MC %d: D* w/ >1 decay vertices\n",
	  runNum,evNum,imct); assert(false);
	}
	const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
	const vector<int> &decayTrks = decayV.vMCtrack();
	if (decayTrks.size()!=2) {                          // ***** 2 DAUGHTERS
 printf("\n** MCInfo:\a Run#%d Evt %d MC %d: D* w/ %d daughters\n",
	  runNum,evNum,imct,(int)decayTrks.size()); continue;
	}
	int imct1 = decayTrks[0], imct2 = decayTrks[1];
	const PaMCtrack &decay1 = e.vMCtrack()[imct1];
	const PaMCtrack &decay2 = e.vMCtrack()[imct2];
	if      (decay1.Pid()==38-C && decay2.Pid()==9-C) {    // D0pi
	    iMCDS_D0[C] = imct1; iMCDS_piS[C] = imct2;
	}
	else if (decay1.Pid()==9-C  && decay2.Pid()==38-C) {   // piD0
	    iMCDS_D0[C] = imct2; iMCDS_piS[C] = imct1;
	}
	else continue;
	DSEffHistosAndFlags(e,C);
      }  // End D0,D* blocks 
    }  // End loop on tracks
  }  // End loop on iterations 1 = D0->Kpi, 2 = D+*->D0pi
}
// ***************************************************************************
// ***************************   SortD0Out_LUND    ***************************
// ***************************************************************************
void MCInfo::SortD0Out_LUND(PaEvent &e, bool richAcceptance)
{
  // ********** PROCESS LUND INFO FOR D0 & D* RECONSTRUCTIBILITIES **********

  // D0 -> Kpi  DS -> Kpipi

  // Expecting "SortOut" to have been already executed


  for (int i = 0; i<2; i++) {
    iMCT_D0[i]=iMCD0_K[i]=iMCD0_pi[i] = -1;
    iMCT_DS[i]=iMCDS_D0[i]=iMCDS_piS[i] = -1;
    D0Decay[i] = 0; D0Type[i]=DSType[i] = 0;
  }

  int runNum = e.RunNum(), evNum = (int)e.UniqueEvNum();

  if (richAcceptance) {

    // *************** RICH ACCEPTANCE && TABULATED PID for ***************
    // ***************       PIONS and KAONS                ***************

    // Will have to extrapolate to RICH and determine whether track w/in its
    // overall acceptance and w/in its special angle or radius cut
    static PaTPar heli2;
    // ***** RICH ACCEPATNCE: GEOMETRY INFO *****
    const double zSM1 = 350;
    const double zRICH = 615.6, rRICH = 5, r2RICH = rRICH*rRICH,
      XRICH = 161, YRICH = 121;
    // ***** TABULATED EFFICIENCY (AND misID EFFICIENCY) *****
    // - RICH perfs can be taken from MC simulation of the RICH response, or
    //  derived from a table of the real perfs as a function of phase space (as
    //  measured e.g. by the K0/phi method).
#if   RICH_EFFICIENCY_VS_PS_2 == 1
    // Ideal
    const double eRICHps[9] = {1,1,1,1,1,1,1,1,1};
    const double mRICHps[9] = {0,0,0,0,0,0,0,0,0};
#elif RICH_EFFICIENCY_VS_PS_2 == 2
    // 2sigmas veto 3sigmas
    const double eRICHps[9] = {.243,.285,.078,
			       .358,.384,.113,
			       .339,.609,.190};  // for >0
    // ``Efficiency'' of misID
    const double mRICHps[9] = {.046,.081,.076,
			       .039,.042,.030,
			       .016,.022,.026};
#endif
#if RICH_EFFICIENCY_VS_PS_2 < 3
    static double pRICHps[9];  // Purity
    static bool first = true; if (first) {
      first = false;
      for (int i = 0; i<9; i++) {
	// Determine purity from efficiency and misID, assuming K/pi = 1/10
	double e = eRICHps[i]; pRICHps[i] = e/(e+10*mRICHps[i]);
      }
    }
#endif

    int imct; for (imct = 1; imct<e.NMCtrack(); imct++) {

      // ********** LOOP ON MC TRACKS, Excluding incident muon **********

      int &Type = Types[imct];
      if (!(Type&0x0003)) continue;      // ***** PRIMARY, pV-gamma decay, !e^+-
      if (NHits[0][imct]==0 && NHits[1][imct]==0 &&// ***** NO HIT WHATSOEVER...
	  NHits[2][imct]==0 && NHits[3][imct]==0) continue;     // ...=> GIVE-UP

      const PaMCtrack  &trk = e.vMCtrack()[imct];
      bool ispi = trk.Pid()==8  || trk.Pid()==9;
      bool isK  = trk.Pid()==11 || trk.Pid()==12;
      if (!ispi && !isK) continue;                         // ***** pi's AND K's

      const PaMCvertex &vtx = e.vMCvertex()[trk.iVertexOfOrigin()];
      double pMC = trk.Mom3().Mag();
      const double KThr = 10;
      if (pMC<=KThr) continue;  // Have to do something for pi's

      // ********** 1ST and LAST POINTS (i.e. origin and decay vtces) **********
      double zFirst = vtx.Pos(2);
      int imcv, endID; static double zLast;     // ***** DECAY VERTEX? *****
      int nvs = (int)trk.vMCvertex().size();
      for (imcv = 0, endID = -1; imcv<nvs; imcv++) {
	PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex()[imcv]];
	if (endID== -1 || decayV.Pos(2)>zLast) {
	  endID = imcv; zLast = decayV.Pos(2);
	}
      }
      if (zSM1<zFirst ||                              // ***** W/IN RICH SPAN...
	  (endID!= -1 && zLast<zRICH &&
	  // ...Yet pions allowed to decay upstream, because their decay into
	  // muons yields a continuous, reco'ible pi-/mu-on track, which RICH
	  // signature is almost identical to that of a pure pion track.
	   !ispi))
	continue;

      //         ********** EXTRAPOLATE to RICH **********
      if (trk.Mom3()[2]<0) continue;  // No backward pointing
      PaTPar helix = trk.ParInVtx();
      if (fabs(helix(5))>1/(.1+.05))  // Too low a momentum for Runge-Kutta
	continue;
      if (fabs(helix(3))>.5 || fabs(helix(4))>.5) // Too large angle
	continue;
      helix.Extrapolate(zRICH,heli2,false);

      if (heli2(1)*heli2(1)+heli2(2)*heli2(2)>r2RICH &&
	  fabs(heli2(1))<XRICH && fabs(heli2(2))<YRICH) {
	//                                         
	Type |= 0x10000000; // ***** RICH'able =  W/IN RICH RANGE and ACCEPTANCE
	if (heli2(1)*heli2(1)+heli2(2)*heli2(2)>RICHCut)
	  Type |= 0x20000000;                           // ***** W/IN INNER RICH
      }
      else continue;
      // Note: we do give up at this point. It may not be fully correct, in so
      // far as tabulated RICH perfs may have been determined w/o a RICH
      // acceptance cut and hence may include, although not in a fully
      // satisfactory way (they're derived from K0 and and phi decays which have
      // a definite phase space distribution), the RICH inefficiency due to
      // tracks outside acceptance.

      //          ********** PHASE SPACE BINNING **********
#if RICH_EFFICIENCY_VS_PS_2 < 3
      double A = acos(1/sqrt(1+heli2(3)*heli2(3)+heli2(4)*heli2(4)));
      int ipa = 0;
      if      (pMC-KThr>20) ipa += 2;
      else if (pMC-KThr>10) ipa += 1;
      if      (A>.024) ipa += 6;
      else if (A>.012) ipa += 3;

      if      (trk.Pid()==8 || trk.Pid()==9) {
	if (drand48()<=mRICHps[ipa]) Type |= 0x80000000;
      }
      else if (trk.Pid()==11 || trk.Pid()==12) {
	if (drand48()<=eRICHps[ipa]) Type |= 0x40000000;
      }
#else          // ***** RICH EFFICIENCY from PERFS TABULATED as a f(PHASE SPACE)
      double A = acos(1/sqrt(1+helix(3)*helix(3)+helix(4)*helix(4)));
      int ia, ja; for (ja = 0, ia = -1; ja<nRICHaBins; ja++)
		    if (A<RICHaBins[ja]) { ia = ja; break; }
      if (ia<0) ia = nRICHaBins-1;
      double pTable; // Either P or P-Thr, depending upon "RICHpLow=0" or not.
      if (RICHpLow) pTable = pMC;
      else          pTable = ispi ? pMC-piThr : pMC-KThr;
      if (pMC>RICHpLow) { // If P not lower than lowest tabulated bin
	int ip, jp; for (jp = 0, ip = -1; jp<nRICHpBins; jp++)
		      if (pTable<RICHpBins[jp]) { ip = jp; break; }
	if (ip<0) ip = nRICHpBins-1;

	if      (trk.Pid()==11) {
	  if (drand48()<eRICHKp[ip][ia])  Type |= 0x40000000;
	}
	else if (trk.Pid()==12) {
	  if (drand48()<eRICHKm[ip][ia])  Type |= 0x40000000;
	}
	else if (trk.Pid()==8) {
	  if (drand48()<eRICHpip[ip][ia]) Type |= 0x80000000;
	}
	else if (trk.Pid()==9) {
	  if (drand48()<eRICHpim[ip][ia]) Type |= 0x80000000;
	}
      }
#endif
    }  // End loop on MC tracks

  }  // ********** END DETERMINE "richAcceptance" **********

  const vector<PaMCgen> &vMCgen = e.vMCgen(); if (vMCgen.empty()) return;
  const int k1DS = 413, k1D0 = 421, k1K = 321, k1pi = 211;
  //#define MCInfo_ALLOW_Kpipi0 1
#ifdef MCInfo_ALLOW_Kpipi0          // ********** ALLOW D0 -> Kpi + X **********
  const int k1rho = 213, k1pi0 = 111;
  const int qqKS  = 320; // K*+(892)=323, K*+(1430)=10321, K*+(1680)=30323
  const int qqKS0 = 310; // K*0(892)=313, K*0(1430)=10311, K*0(1680)=30313
  const int k892 = 3, k1430 = 10001, k1680 = 30003;  
#endif
  int igenD0s[2] = {0,0}; for (int iter = 0; iter<2; iter++) {

    // ********** 1ST ITER: D0->K-pi+ **********
    // ********** 2ND ITER: D0*->D0pi+ (using results of 1st iter) **********

    for (int igen = 0; igen<(int)vMCgen.size(); igen++) {

      // ********** LOOP OVER LUJET INFO. BLOCKS **********

      if (vMCgen[igen].vInt().back()!=102) continue;//   ***** It's LUJET *****
      LUJET lj = vMCgen[igen].LundJet();

      if (iter==0 && abs(lj.k[1])==k1D0) {
	//         ******************************************
	//         ********** 1ST ITER: D0/anti-D0 **********
	//         ******************************************
	int C = lj.k[1]==k1D0 ? 1 : 0, pm = 1-2*C /* K polarity */, mp = -pm;
	if (igenD0s[C])
 printf("\n** MCInfo:\a Run#%d Evt %d: >1 %s (LJs %d,%d)\n",
	  runNum,evNum,C?"D0":"aD0",igenD0s[C],igen);
	igenD0s[C] = igen;
	int nDecays = lj.k[4]-lj.k[3]+1, idecay1, idecay2;
#ifdef MCInfo_ALLOW_Kpipi0          // ********** ALLOW D0 -> Kpi + X **********
	int ipi0 = 0, branch = 0; // 0x1: Krho, 0x2:iKS0pi, 0x4: iKSpi; 
	// We evaluate 5 possibilites (times charge conjugation):
	//  2 daughters: D0 -> K-pi+, K-rho+, K*0pi0, K*pi+
	//  3 daughters: D0 -> K-pi+pi0
	// Cases w/ 4 particles in the final state (K-pi+p+pi-) cannot be
	// integrated in the present scheme (MC class w/ solely "iMCD0_K" and
	// "iMCD0_pi"), because of the ambiguity on which pi+ to retain.
	// The ouputs of the evaluation are vMCgen indices "idecay1/2", which
	// identify K and pi indifferently.
	if (nDecays==2) {                                // ***** 2 DAUGHTERS...
	  idecay1 = lj.k[3]; idecay2 = lj.k[4];
	  LUJET lj1 = vMCgen[idecay1].LundJet();
	  LUJET lj2 = vMCgen[idecay2].LundJet();
	  int k11 = lj1.k[1], k12 = lj2.k[1];
	  int qq1 = abs(k11)%1000/10*10, qq2 = abs(k12)%1000/10*10;
	  if      (k11==pm*k1K && k12==mp*k1rho ||             // *****...K-rho+
		   k12==pm*k1K && k11==mp*k1rho) {
	    int irho; if (k11==pm*k1K) { irho = idecay2; }
	    else {                       irho = idecay1; idecay1 = idecay2; }
	    LUJET ljrho = vMCgen[irho].LundJet();
	    if (ljrho.k[4]-ljrho.k[3]+1!=2) continue;
	    LUJET ljrho1 = vMCgen[ljrho.k[3]].LundJet();
	    LUJET ljrho2 = vMCgen[ljrho.k[4]].LundJet();
	    if      (ljrho1.k[1]==mp*k1pi && ljrho2.k[1]==k1pi0) {
	      idecay2 = ljrho.k[3]; ipi0 = ljrho.k[4];
	    }
	    else if (ljrho2.k[1]==mp*k1pi && ljrho1.k[1]==k1pi0) {
	      idecay2 = ljrho.k[4]; ipi0 = ljrho.k[3];
	    }
	    else continue; // Non pipi0 decay of rho
	    branch = 0x6;  // 0x2=Kpipi0 + 0x4=Krho
	  }
	  else if (k11*pm>0 && qq1==qqKS && k12==mp*k1pi ||    // *****...K*-pi+
		   k12*pm>0 && qq2==qqKS && k11==mp*k1pi) {
	    int iKS; if (qq2==qqKS) { iKS = idecay2; }
	    else {                    iKS = idecay1; idecay1 = idecay2; }
	    LUJET ljKS = vMCgen[iKS].LundJet();	    
	    if (ljKS.k[4]-ljKS.k[3]+1!=2) continue;
	    LUJET ljKS1 = vMCgen[ljKS.k[3]].LundJet();
	    LUJET ljKS2 = vMCgen[ljKS.k[4]].LundJet();
	    if      (ljKS1.k[1]==pm*k1K && ljKS2.k[1]==k1pi0) {
	      idecay2 = ljKS.k[3]; ipi0 = ljKS.k[4];
	    }
	    else if (ljKS2.k[1]==pm*k1K && ljKS1.k[1]==k1pi0) {
	      idecay2 = ljKS.k[4]; ipi0 = ljKS.k[3];
	    }
	    else continue; // Non Kpi0 decay of K*
	    branch = 0xa;  // 0x2=Kpipi0 + 0x8=K*pi
	    int k1001 = pm*ljKS.k[1]-qqKS;
	    if      (k1001==k892)  branch |= 0x100;
	    else if (k1001==k1430) branch |= 0x200;
	    if      (k1001==k1680) branch |= 0x400;
	  }
	  else if (k11*pm>0 && qq1==qqKS0 && k12==k1pi0 ||     // *****...K*0pi0
		   k12*pm>0 && qq2==qqKS0 && k11==k1pi0) {
	    int iKS0; if (qq1==qqKS0) { iKS0 = idecay1; ipi0 = idecay2; }
	    else {                      iKS0 = idecay2; ipi0 = idecay1; }
	    LUJET ljKS0 = vMCgen[iKS0].LundJet();
	    if (ljKS0.k[4]-ljKS0.k[3]+1!=2) continue;
	    LUJET ljKS01 = vMCgen[ljKS0.k[3]].LundJet();
	    LUJET ljKS02 = vMCgen[ljKS0.k[4]].LundJet();
	    if      (ljKS01.k[1]==pm*k1K && ljKS02.k[1]==mp*k1pi ||
		     ljKS02.k[1]==pm*k1K && ljKS01.k[1]==mp*k1pi) {
	      idecay1 = ljKS0.k[4]; idecay2 = ljKS0.k[3];
	    }
	    else continue; // Non Kpi decay of K*0
	    branch = 0x12; // 0x2=Kpipi0 + 0x10=K*0pi0
	    int k1001 = pm*ljKS0.k[1]-qqKS0;
	    if      (k1001==k892)  branch |= 0x100;
	    else if (k1001==k1430) branch |= 0x200;
	    if      (k1001==k1680) branch |= 0x400;
	  }
	  else if (k11==pm*k1K && k12==mp*k1pi ||               // *****...K-pi+
		   k12==pm*k1K && k11==mp*k1pi) {
	    idecay1 = lj.k[3]; idecay2 = lj.k[4];
	    branch = 0x1;
	  }
	  //#  define MCInfo_ALLOW_KK
#  ifdef MCInfo_ALLOW_KK          // ********** ALLOW D0 -> K+K- **********
	  else if (k11==pm*k1K && k12==mp*k1K ||                 // *****...K-K+
		   k12==pm*k1K && k11==mp*k1K) {
	    idecay1 = lj.k[3]; idecay2 = lj.k[4];
	    branch = 0x1000;
	  }
#  endif
	  else continue;
	}
	else if (nDecays==3) {                           // ***** 3 DAUGHTERS...
	  int idecay3, *idecays[3] = {&idecay1,&idecay2,&idecay3};
	  int idecay, nChargedDecays, semilepton;
	  for (idecay = lj.k[3], nChargedDecays=semilepton = 0; idecay<=lj.k[4];
	       idecay++) {
	    int k1 = vMCgen[idecay].LundJet().k[1];
	    if (abs(k1)<100) { semilepton = true; break; }
	    if ((k1/10%10+k1/100%10)%2) // q and \bar{q} of unlike parity
	      *idecays[nChargedDecays++] = idecay;
	    else if (k1==111) ipi0 = idecay;
	  }
	  if (semilepton) continue;
	  else if (nChargedDecays==3) {
 printf("** MCInfo:\a Run %d Evt %d LU %d: D0 w/ 3 charged decays: %d,%d,%d\n",
	    runNum,evNum,igen,idecay1,idecay2,idecay3); continue;
	  }
	  else if (nChargedDecays<2) continue;
	  if (!ipi0) continue;                              // ***** ...w/ pi0
	  branch = 0x2;
	}
	else continue; // >3 decays
#else                              // ********** ALLOW ONLY D0 -> Kpi **********
	if (nDecays==2) {
	  idecay1 = lj.k[3]; idecay2 = lj.k[4];
	}
	else continue;
#endif	  
	LUJET lj1 = vMCgen[idecay1].LundJet();
	LUJET lj2 = vMCgen[idecay2].LundJet();
	int k11 = lj1.k[1], k12 = lj2.k[1];
	int imct1 = lj1.lu2kine, imct2 = lj2.lu2kine;
	if (imct1<0 || imct2<0) { // No associated MC tracks...
	  if ((k11==pm*k1K && k12==mp*k1pi) ||// ...should not happen for K and pi
	      (k12==pm*k1K && k11==mp*k1pi))
 printf("\n** MCInfo:\a Run#%d Evt %d LJ %d: D0->Kpi[X] w/ no K|pi MC track\n",
	    runNum,evNum,igen); continue;
	}
	if      (k11==pm*k1K && k12==mp*k1pi) {                 // ***** -> K pi
	  iMCD0_K[C] = imct1; iMCD0_pi[C] = imct2;
	}
	else if (k12==pm*k1K && k11==mp*k1pi) {                 // ***** -> K pi
	  iMCD0_K[C] = imct2; iMCD0_pi[C] = imct1;
	}
#  ifdef MCInfo_ALLOW_KK
	else if (k11==pm*k1K && k12==mp*k1K ||                  // ***** -> K+K-
		 k12==pm*k1K && k11==mp*k1K) {
	  iMCD0_K[C] = imct1; iMCD0_pi[C] = imct2;
	}
#  endif
	else continue;
#ifdef MCInfo_ALLOW_Kpipi0
	D0Decay[C] = branch;
#else
	D0Decay[C] = 0x1;
#endif
	D0EffHistosAndFlags(e,C);
      }  // End of D0 block
      else if (iter==1 && abs(lj.k[1])==k1DS) {     // ***** D*+/-
	//         ******************************************
	//         ********** 2ND ITER: D*->D0pi **********
	//         ******************************************
	int C = lj.k[1]==k1DS ? 1 : 0, mp = 2*C-1 /* piS polarity */;
	int nDecays = lj.k[4]-lj.k[3]+1, idecay1, idecay2;
	if (nDecays!=2) {
 printf("\n** MCInfo:\a Run#%d Evt %d LUJET %d: D* w/ %d daughters\n",
  	  runNum,evNum,igen,nDecays); continue;          // ***** 2 DAUGHTERS...
	}
	idecay1 = lj.k[3]; idecay2 = lj.k[4];
	LUJET lj1 = vMCgen[idecay1].LundJet();
	LUJET lj2 = vMCgen[idecay2].LundJet();
	int k11 = lj1.k[1], k12 = lj2.k[1], idecayD0 = 0;
	if      (k11==mp*k1D0 && k12==mp*k1pi) {                 // *****...D0pi
	  idecayD0 = idecay1; iMCDS_piS[C] = lj2.lu2kine;
	}
	else if (k12==mp*k1D0 && k11==mp*k1pi) {                // ***** ...piD0
	  idecayD0 = idecay2; iMCDS_piS[C] = lj1.lu2kine;
	}
	else continue;
	if (iMCDS_piS[C]<0) {// No associated MC track: shouldn't happen for pi
 printf("\n** MCInfo:\a Run#%d Evt %d LJ %d: D*->D0pi w/ no pi MC track\n",
  	    runNum,evNum,igen); continue;
	}
	if (idecayD0!=igenD0s[C]) {
 printf("\n** MCInfo:\a Run#%d Evt %d: >1 %s: LJ#%d (<-D*%c) and LJ#%d\n",
	  runNum,evNum,C?"D0":"aD0",igen,C?'+':'-',igenD0s[C]); continue;
	}
	DSEffHistosAndFlags(e,C);
      }  // End D0/D* blocks
    }  // End loop on LUJET info
  }
}

// ***************************************************************************
// *************************    D0EffHistAndFlags    *************************
// ***************************************************************************
void MCInfo::D0EffHistosAndFlags(PaEvent &e, int C)
{
  // - Histos, reco'ibility flags ("Types" and "D0/DSType").
  // - Also momenta of D0, D* and K,pi decays. For D0 and D*, the value is
  //  computed from the decays, even if their very PaMCtracks exist.
  // - Compute also cos(theta*).

  int runNum = e.RunNum(), evNum = (int)e.UniqueEvNum();

  //       ********** COMPUTE D0 KINEMATICAL VARIABLES **********

  const PaMCtrack &K  = e.vMCtrack()[iMCD0_K[C]];
  const PaMCtrack &pi = e.vMCtrack()[iMCD0_pi[C]];
  const TLorentzVector lvK = K.LzVec(), lvpi = pi.LzVec();
  TLorentzVector lvD0 = lvpi; lvD0 += lvK; TLorentzVector lvK2 = lvK;
  if (D0Decay[C]==0x1 &&          // If it's a (a)D0->Kpi...
      fabs(lvD0.M()-M_D0)>.001) { // ...X-check D0 mass
 printf("\n** MCInfo:\a Run %d Evt %d MC %d+%d: %s mass = %.3f GeV\n",
    runNum,evNum,iMCD0_K[C],iMCD0_pi[C],C?"D0":"aD0",lvD0.M()); assert(false);
  }
  double pK = lvK.P(), ppi = lvpi.P(), pD0 = lvD0.P();
  PD0[C] = pD0; PD0_K[C] = pK; PD0_pi[C] = ppi;
  lvK2.Boost(-lvD0.BoostVector()); // Boost K to D0 frame
  double cths = cos(lvK2.Angle(lvD0.Vect())); cthstarD0[C] = cths;
  double zD = lvD0.E()/nu; zD0[C] = zD;
  hD_cs->Fill(cths);

  //       ********** SET "Types" FLAGS AND FILL HISTOS **********

  Types[iMCD0_K[C]] |= 0x01000000; Types[iMCD0_pi[C]] |= 0x02000000; 
  hp_aD0[0]->Fill(pD0);  hp_aD0[1]->Fill(pK);
  hp_aD0[2]->Fill(cths); hp_aD0[3]->Fill(zD);
  unsigned int &Type = D0Type[C]; Type |= 0x00100000; // Flag (a)D0 -> Kpi
  int KType = Types[iMCD0_K[C]];
  if (cLowCutS<cths && cths<cUpCutS) {
    //    ***** D0 W/IN LOOSE (for eventual D* TAGGED SELECTION) cth* CUT *****
    if (lvD0.E()/nu>zCutS) Type |= 0x00800000; // Flag w/in D* kin cuts
    if (cLowCut0<cths && cths<cUpCut0) {
      Type |= 0x00200000;                              // ***** D0 W/IN cth* CUT
      hp_aD0c[0]->Fill(pD0);  hp_aD0c[1]->Fill(pK);
      hp_aD0c[2]->Fill(cths); hp_aD0c[3]->Fill(zD);
      if (zD>zCut0) {
	Type |= 0x00400000;                         // ***** D0 W/IN cth*,z CUTS
	hp_aD0cz[0]->Fill(pD0);  hp_aD0cz[1]->Fill(pK);
	hp_aD0cz[2]->Fill(cths); hp_aD0cz[3]->Fill(zD);
	bool winTarget = XpV*XpV+YpV*YpV<R2Cut;
	if (Year==2012 || Year==2016 || Year==2017) {
	  winTarget |= PaAlgo::InTarget(XpV,YpV,ZpV,'O',PaEvent::Ptr()->RunNum());
	}
	if (yLowCut<yB && yB<yUpCut &&                 // ***** EVENT CUTS *****
	    winTarget) {
	  hp_aD0k[0]->Fill(pD0);  hp_aD0k[1]->Fill(pK);
	  hp_aD0k[2]->Fill(cths); hp_aD0k[3]->Fill(zD);
	  if (KType&0x10000000) {
	    hp_aD0ki[0]->Fill(pD0);  hp_aD0ki[1]->Fill(pK);
	    hp_aD0ki[2]->Fill(cths); hp_aD0ki[3]->Fill(zD);
	    if (KType&0x20000000) {
	      hp_aD0kI[0]->Fill(pD0);  hp_aD0kI[1]->Fill(pK);
	      hp_aD0kI[2]->Fill(cths); hp_aD0kI[3]->Fill(zD);
	    }
	  }
	}
      }
    }
  }
  if ((Types[iMCD0_K[C]]&0x00000060) && (Types[iMCD0_pi[C]]&0x00000060)) {
    // Reconstructibility histos = denominator of efficiencies
    // Notes:
    //  - Reconstructibility includes w/in target fiducial volume, cf. supra
    //  - The numerators (of efficiencies), cf. "UserEvent5" have to be
    //   reco'ibility * reco'd and not reco'ibility alone.
    Type |= 0x0080;                    // ***** FLAG D0 AS RECO'IBLE *****
#ifdef MC_DEBUG_RECOIBILITY
    ribilityFilled |= 1<<C;
#endif
    hp_RD0[0]->Fill(pD0);  hp_RD0[1]->Fill(pK);
    hp_RD0[2]->Fill(cths); hp_RD0[3]->Fill(zD);
    if (Type&0x200000) {
      hp_RD0c[0]->Fill(pD0);  hp_RD0c[1]->Fill(pK);
      hp_RD0c[2]->Fill(cths); hp_RD0c[3]->Fill(zD);
      if ((Type&0x00400000) /* w/in z cut */ && yLowCut<yB && yB<yUpCut) {
	hp_RD0k[0]->Fill(pD0);  hp_RD0k[1]->Fill(pK);
	hp_RD0k[2]->Fill(cths); hp_RD0k[3]->Fill(zD);
	if (KType&0x10000000) {    // ***** K is RICH'able *****
	  hp_RD0ki[0]->Fill(pD0);  hp_RD0ki[1]->Fill(pK);
	  hp_RD0ki[2]->Fill(cths); hp_RD0ki[3]->Fill(zD);
	  if (KType&0x20000000) {  // ***** K w/in RICH Cut *****
	    hp_RD0kI[0]->Fill(pD0);  hp_RD0kI[1]->Fill(pK);
	    hp_RD0kI[2]->Fill(cths); hp_RD0kI[3]->Fill(zD);
	  }
	}
      }
    }
  }
}
// ***************************************************************************
// *************************    DSEffHistAndFlags    *************************
// ***************************************************************************
void MCInfo::DSEffHistosAndFlags(PaEvent &e, int C)
{
  if (!(D0Type[C]&0x00100000)) return;
  int runNum = e.RunNum(), evNum = (int)e.UniqueEvNum();

  //       ********** COMPUTE D* KINEMATICAL VARIABLES **********

  const PaMCtrack &K   = e.vMCtrack()[iMCD0_K[C]];
  const PaMCtrack &pi  = e.vMCtrack()[iMCD0_pi[C]];
  const PaMCtrack &piS = e.vMCtrack()[iMCDS_piS[C]];
  const TLorentzVector lvK = K.LzVec(), lvpi = pi.LzVec(), lvpiS = piS.LzVec();
  TLorentzVector lvDS = lvpi; lvDS += lvK; lvDS += lvpiS;
  if (D0Decay[C]==0x1 &&          // If it's a (a)D0->Kpi...
      fabs(lvDS.M()-M_DS)>.001) { // X-check D* mass
 printf("\n** MCInfo:\a Run %d Evt %d MC %d%d+%d: D*%c Mass = %.3f GeV\n",
    runNum,evNum,iMCD0_K[C],iMCD0_pi[C],iMCDS_piS[C],C?'+':'-',lvDS.M());
    assert(false);
  }
  PDS[C] = lvDS.P(); double pD0 = PD0[C], pK = PD0_K[C];
  double cths = cthstarD0[C], zD = zD0[C];

  //       ***** SET "Type" FLAGS AND FILL HISTOS for D* *****

  // If daughter D0 (be it #1 or #2) decays -> Kpi
  Types[iMCDS_piS[C]] |= 0x04000000;
  hp_aDS[0]->Fill(pD0);  hp_aDS[1]->Fill(pK);
  hp_aDS[2]->Fill(cths); hp_aDS[3]->Fill(zD);
#define pi_FF_REJECTION 3
#if ! defined pi_FF_REJECTION || ! ( pi_FF_REJECTION & 2 ) 
  // If one allows fringe-field soft pion, this particle may have
  // 0x0100. While D0 can only be 0x0080 or nul.
  bool recoibleDS = (D0Type[C]&0x0080) && (Types[iMCDS_piS[C]]&0x0160);
#else
  bool recoibleDS = (D0Type[C]&0x0080) && (Types[iMCDS_piS[C]]&0x0060);
#endif
  if (recoibleDS) {
    DSType[C] |= 0x0080;                                  // ***** FLAG D* *****
    if ((D0Type[C]&0x00800000)
	// If daughter D0 (be it #1 or #2) w/in loose kin cuts
	&& yLowCut<yB && yB<yUpCut) {       // And event is w/in y cut
      hp_RDSk[0]->Fill(pD0);  hp_RDSk[1]->Fill(pK);
      hp_RDSk[2]->Fill(cths); hp_RDSk[3]->Fill(zD);
      int KType = Types[iMCD0_K[C]];
      if (KType&0x10000000) {                      // ***** K is RICH'able *****
	hp_RDSki[0]->Fill(pD0);  hp_RDSki[1]->Fill(pK);
	hp_RDSki[2]->Fill(cths); hp_RDSki[3]->Fill(zD);
	if (KType&0x20000000) {                   // ***** K w/in RICH Cut *****
	  hp_RDSkI[0]->Fill(pD0);  hp_RDSkI[1]->Fill(pK);
	  hp_RDSkI[2]->Fill(cths); hp_RDSkI[3]->Fill(zD);
	}
      }
    }
  }
}

// ***************************************************************************
// *************************   BookLambdaAcceptance  *************************
// ***************************************************************************
void MCInfo::BookLambdaAcceptance(double yLow, double yUp,
				  double R2, // If 2016/17, apply PaAlgo::InTarget in addition
				  double pTCut)
{
  //   ********** BOOK ``DENOMINATOR'' HISTOGRAMS **********
  // + VARIOUS INITIALISATIONS
  // Note: Would be good to create these histo in a "MC/Lambda" TDIrectory

  yLowCut = yLow; yUpCut = yUp;                         // ***** INITIALIZATIONS
  R2Cut = R2;
  V0pTCut = pTCut; V0cthCut = 0;

  // ***** Lambda -> ppi: ACCEPTANCE, RECONSTRUCTIBILITY *****
  char title[] =
    "Reco'ible #bar{#Lambda}* vs. c#theta* - pT>25MeV,VV.V0>0.99990,RICH      ";
  //"Reco'ible #bar{#Lambda}* vs. c#theta* (Q2>1,.35<y<.85)      ";
  for (int iLaL = 0; iLaL<2; iLaL++) {  // Loop on Lambda/anti-Lambda

    // Successively: all Ls
    //               - w/in kin cuts
    //               - RICH'able w/in RICH acceptance cuts
    //               - w/in kin and RICH'able
    //               Reco'ible
    //               - w/in...
    char hName[] = "hc_aALki"; char LaL[] = "#bar{#Lambda}", LAL[] = "AL";
    if (iLaL) { sprintf(LaL,"bar{#Lambda}"); sprintf(LAL,"AL"); }
    else      { sprintf(LaL,"#Lambda");      sprintf(LAL,"L"); }
    const char *vars[] = {"c#theta*","Pp"}, cps[] = "cp"; 
    for (int icp = 0; icp<2; icp++) {
      const char *var = vars[icp], cp = cps[icp];
      int nBins; double xMn, xMx;
      if (icp==0) { nBins = 16; xMn = -1; xMx = 1; }
      else        { nBins = 40; xMn =  0; xMx = 80; }
      // "All" histos: hc_a..
      sprintf(hName,"h%c_a%s",cp,LAL);
      sprintf(title,"All %s vs. %s (Q2>1,%.2f<y<%.2f)",LaL,var,yLow,yUp);
      hc_aL[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_a%sk",cp,LAL);
      sprintf(title,"All %s vs. %s - pt>%.0fMeV,V0.VV>%.5f",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_aLk[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_a%si",cp,LAL);
      sprintf(title,"All %s vs. %s - RICH",LaL,var);
      hc_aLi[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_a%ski",cp,LAL);
      sprintf(title,"All %s vs. %s - pt>%.0fMeV,V0.VV>%.5f,RICH",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_aLki[icp][iLaL] = new TH1D(hName,title,nBins,xMn,xMx);
      // "Reconstructible" histos: hc_R..
      sprintf(hName,"h%c_R%s",cp,LAL);
      sprintf(title,"Reco'ible %s vs. %s (Q2>1,%.2f<y<%.2f)",LaL,var,yLow,yUp);
      hc_RL[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_R%sk",cp,LAL);
      sprintf(title,"Reco'ible %s vs. %s - pt>%.0fMeV,V0.VV>%.5f",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_RLk[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_R%si",cp,LAL);
      sprintf(title,"Reco'ible %s vs. %s - RICH",LaL,var);
      hc_RLi[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_R%ski",cp,LAL);
      sprintf(title,"Reco'ible %s vs. %s - pt>%.0fMeV,V0.VV>%.5f,RICH",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_RLki[icp][iLaL] = new TH1D(hName,title,nBins,xMn,xMx);
    }

    // ***** Sigma* -> Lpi: ACCEPTANCE, RECONSTRUCTIBILITY *****
    // Successively all Sigma*
    //              Reco'ible w/in kin cuts and RICH'able
    if (iLaL) { sprintf(LAL,"AS"); }
    else      { sprintf(LAL,"S"); }
    sprintf(hName,"hc_a%s",LAL);
    sprintf(title,"All %s* vs. c#theta* (Q2>1,%.2f<y<%.2f)",LaL,yLow,yUp);
    hc_aS[iLaL]   = new TH1D(hName,title,16,-1,1);
    sprintf(hName,"hc_R%ski",LAL);
    sprintf(title,"Reco'ible %s* vs. c#theta* - pt>%.0fMeV,V0.VV>%.5f,RICH",
	    LaL,V0pTCut*1000,V0cthCut);
    hc_RSki[iLaL] = new TH1D(hName,title,16,-1,1);
  }

  hD_ReLambda =
    new TH1D("hD_ReLambda","#Lambda(#bar{#Lambda}) from ReInteractions: D (cm)",
	     100,0,1000);
}
// ***************************************************************************
// ****************************** SortLambdaOut ******************************
// ***************************************************************************
void MCInfo::SortLambdaOut(PaEvent &e, bool richAcceptance)
{
  // ********** PROCESS MC RAW INFO FOR (anti-)Lambda ACCEPTANCE **********

  // L (resp. anti-L) -> p + pi^- (resp. bar{p} + pi^+)

  // Expecting "SortOut" to have been already executed


  for (int i = 0; i<2; i++) {
    iMCT_L[i]=iMCL_p[i]=iMCL_pi[i] = -1;
    iMCT_S[i]=iMCS_L[i]=iMCS_pi[i] = -1;
  }

  if (richAcceptance) {

    // *************** RICH ACCEPTANCE && TABULATED PID for ***************
    // ***************       PIONS and KAONS                ***************

    // Extrapolate to RICH and determine whether track w/in its overall
    // acceptance and with its special angle or radius cut
    static PaTPar heli2;
    // ***** RICH ACCEPATNCE: GEOMETRY INFO *****
    const double zSM1 = 350;
    const double zRICH = 750, rRICH = 5, r2RICH = rRICH*rRICH,
      XRICH = 161, YRICH = 121;
    // ***** TABULATED EFFICIENCY (AND misID EFFICIENCY) *****
    // ***** TABULATED EFFICIENCY (AND misID EFFICIENCY) *****
    // - RICH perfs can be taken from MC simulation of the RICH response, or
    //  derived from a table of the real perfs as a function of phase space (as
    //  measured e.g. by the K0/phi method).
#define RICH_EFFICIENCY_VS_PS 1
#if   RICH_EFFICIENCY_VS_PS == 1
    // Ideal
    const double eRICHps[9] = {1,1,1,1,1,1,1,1,1};
#elif RICH_EFFICIENCY_VS_PS == 2
    // 2sigmas veto 3sigmas
    const double eRICHps[9] = {.243,.285,.078,
			       .358,.384,.113,
			       .339,.609,.190};  // for >0
#endif
    // ``Efficiency'' of misID
    const double mRICHps[9] = {.046,.081,.076,
			       .039,.042,.030,
			       .016,.022,.026};
    static double pRICHps[9];  // Purity
    static bool first = true;
    if (first) {
      first = false;
      for (int i = 0; i<9; i++) {
	// Determine "purity" from "efficiency" and "misID"
	// assuming K/pi = 1/10
	double e = eRICHps[i];
	pRICHps[i] = e/(e+10*mRICHps[i]);
      }
    }

    int imct; for (imct = 1; imct<e.NMCtrack(); imct++) {

      // ********** LOOP ON MC TRACKS, Excluding incident muon **********

      int &Type = Types[imct];
      if (!(Type&0x0003)) continue;// ***** OK= primaries, pV-gamma decay, !e^+-
      if (NHits[0][imct]==0 && NHits[1][imct]==0 &&  // ***** NO HIT WHATSOEVER!
	  NHits[2][imct]==0 && NHits[3][imct]==0) continue;         // => BYPASS

      const PaMCtrack  &trk = e.vMCtrack()[imct];
      bool ise  = trk.Pid()==2  || trk.Pid()==3;
      bool ispi = trk.Pid()==8  || trk.Pid()==9;
      bool isK  = trk.Pid()==11 || trk.Pid()==12;
      bool isp  = trk.Pid()==14 || trk.Pid()==15;
      if (!(ise || ispi || isK || isp)) continue;         // ***** e, pi, K or p

      const PaMCvertex &vtx = e.vMCvertex()[trk.iVertexOfOrigin()];

      // ***** CUT on P > THRESHOLD.
      double pMC = trk.Mom3().Mag();
      //printf("evt#%d: %d pid %d %s %6.1f %f  %f\n",(int)e.UniqueEvNum(),imct,trk.Pid(),trk.Name().c_str(),pMC,vtx.Pos(2),zSM1);
      // - What threshold? I.e. Which of UV or visible? What actual index(ices)?
      // - For the time being: let's have n(UV) = 1.001522
      // - Will have to re-evaluate the matter later on....
      const double piThr   = 2.5; //  KThr = 8.9, pThr = 17.0;
      // - W/ n(visible) = 1.001367
      //const double piThr = 2.67, KThr = 9.4, pThr = 17.9;
      double thr = 0; // Case of e^+/-
      if (ispi || isp || isK) thr = piThr; //if (isK) thr = KThr; if (isp) thr = pThr;
      if (pMC<=thr) continue;

      // ********** 1ST and LAST POINTS (i.e. origin and decay vtces) **********
      double zFirst = vtx.Pos(2);
      int imcv, endID; static double zLast;     // ***** DECAY VERTEX? *****
      int nvs = (int)trk.vMCvertex().size();
      for (imcv = 0, endID = -1; imcv<nvs; imcv++) {
	PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex()[imcv]];
	if (endID== -1 || decayV.Pos(2)>zLast) {
	  endID = imcv; zLast = decayV.Pos(2);
	}
      }
      if (zSM1<zFirst ||
	  (endID!= -1 && zLast<zRICH)) continue;   // ***** W/IN RICH SPAN *****

      if (trk.Mom3()[2]<0) continue;         // ***** NO BACKWARD POINTING *****

      PaTPar helix = trk.ParInVtx();
      if (fabs(helix(5))>1/(.1+.05))// ***** TOO LOW MOMENTUM Runge-Kutta) *****
	continue;
      if (fabs(helix(3))>.5 || fabs(helix(4))>.5) // ***** TOO LARGE ANGLE *****
	continue;
			       
      // ********** EXTRAPOLATE to RICH **********
      helix.Extrapolate(zRICH,heli2,false);

      //#define DEBUG_Lambda_RICH
#ifdef DEBUG_Lambda_RICH
      if (isp) printf("evt#%d: %d pid %d %s %6.1f %f  %f %f %s\n",(int)e.UniqueEvNum(),imct,trk.Pid(),trk.Name().c_str(),pMC,sqrt(heli2(1)*heli2(1)+heli2(2)*heli2(2)),heli2(1),heli2(2),!(heli2(1)*heli2(1)+heli2(2)*heli2(2)<r2RICH || fabs(heli2(1))>XRICH || fabs(heli2(2))>YRICH)?"-> OK":"");
#endif
      if (heli2(1)*heli2(1)+heli2(2)*heli2(2)<r2RICH ||
	  fabs(heli2(1))>XRICH || fabs(heli2(2))>YRICH)
	continue;                            // ***** W/IN RICH ACCEPTANCE *****
      Type |= 0x10000000;  // RICH'able

      // ********** PHASE SPACE BIN **********
      double A = acos(1/sqrt(1+heli2(3)*heli2(3)+heli2(4)*heli2(4)));
      int ipa = 0;
      if      (pMC-thr>20) ipa += 2;
      else if (pMC-thr>10) ipa += 1;
      if      (A>.024) ipa += 6;
      else if (A>.012) ipa += 3;

      if      (ispi || isK) {
	if (drand48()<=mRICHps[ipa]) Type |= 0x80000000;
      }
      else if (isp) {
	if (drand48()<=eRICHps[ipa]) Type |= 0x40000000;
      }

    }  // End loop on MC tracks

  }  // ********** END DETERMINE "richAcceptance" **********

  if (e.vMCvertex().empty()) return;
  const PaMCvertex &vDIS = e.vMCvertex()[0];

  for (int iter = 0; iter<2; iter++) {

    // ********** 1ST ITER: L  -> p pi **********
    // ********** 2ND ITER: S*/Xi -> L pi (using results of 1st iter) **********

    int imct = 0, nLs[2] = {0,0};
    imct++;                                       // ***** EXCLUDE INCIDENT MUON
    for (;  imct<e.NMCtrack(); imct++) {

      int &Type = Types[imct];
      if (Type&0x00004000) continue;                       // ***** NOT A PILEUP
      const PaMCtrack  &trk = e.vMCtrack()[imct];

      //printf("Evt %d imct %d pid %d \"%s\"\n",(int)e.UniqueEvNum(),imct,trk.Pid(),trk.Name().c_str());
      if (iter==0 && (trk.Pid()==18 || trk.Pid()==26)) {         // ***** L, a-L
	int C = trk.Pid()==18? 0 : 1;
	//#define NO_aL
#ifdef NO_aL
	if (C==1) continue;
#endif
	//            ***** (a)L->ppi RECONSTRUCTIBILITY
	// ***** Lambda TO ARISE (directly or via heavy hyperon) FROM DIS *****
	// Hence, discrading those Lambda's that come from the re-interaction
	// of primaries w/ hte materila of the spectrometer, which may not
	// be negligible given that we have strangeness in the event.
	int ivtxm = trk.iVertexOfOrigin(), imctm;
	const PaMCvertex &vtxm = e.vMCvertex()[ivtxm];
	if ((imctm = vtxm.iTrackOfOrigin())==-1) { // No mother track..
	  // ...means primary DIS, since pile-up tracks are discarded supra. Let
	  // us check, to be on the safe side.
	  if (ivtxm!=0) {
  printf("\n** MCInfo:\a Run#%d Evt %d: Primary !pileup vertex w/ index(=%d) != 0\n",
	    e.RunNum(),(int)e.UniqueEvNum(),ivtxm); exit(1);
	  }
	}
	else {// ...Else could be a re-interaction or a heavy hyperon decay
	  // We single out most heavy hyperon by asking that mother vertex close
	  // to DIS vertex (case of Sigma0, Sigma*).
	  // N.B.: Is this useful? Are Sigma0 and Sigma* in vMCtrack? Usually
	  // not. But I guess it does happen if Sigma0/* are, e.g., forbidden to
	  // decay in LEPTO.
	  double dx = vtxm.Pos(0)-vDIS.Pos(0), dy = vtxm.Pos(1)-vDIS.Pos(1);
	  double dz = vtxm.Pos(2)-vDIS.Pos(2), d = sqrt(dx*dx+dy*dy+dz*dz);
	  if (d>.01 /* not too strict to make for numerical precision */) {
	    const PaMCtrack &mother = e.vMCtrack()[imctm];
	    int mid = mother.Pid(), ivtxgm = mother.iVertexOfOrigin();
	    if (!((22<=mid && mid<=24) || (30<=mid && mid<=32)) ||
		ivtxgm!=0) {      // Single out primary Xi & Omega..
	      hD_ReLambda->Fill(d);
	      //#define DEBUG_ReLambda
#ifdef DEBUG_ReLambda
  printf("\n** MCInfo:\a Run#%d Evt %d: ",e.RunNum(),(int)e.UniqueEvNum());
	      if (imctm==iMCT_L[C])  // Primary Lambda has re-interacted
		printf("Primary \"%s\" re-nteracting\n",trk.Name().c_str());
	      else {
		const PaMCvertex &vtxgm = e.vMCvertex()[ivtxgm];
		dx = vtxm.Pos(0)-vtxgm.Pos(0); dy = vtxm.Pos(1)-vtxgm.Pos(1);
		dz = vtxm.Pos(2)-vtxgm.Pos(2); d = sqrt(dx*dx+dy*dy+dz*dz);
		if (d<.01) {
		  const PaMCtrack &gm = e.vMCtrack()[vtxgm.iTrackOfOrigin()];
  printf("\"%s\" re-interacting -> \"%s\" -> \"%s\"\n",
	 gm.Name().c_str(),mother.Name().c_str(),trk.Name().c_str());
		}
		else 
  printf("\"%s\" re-interacting -> \"%s\"\n",
	 mother.Name().c_str(),trk.Name().c_str());
  printf(" DIS vertex (%.2f,%.2f,%.2f) -> Re-interaction (%.2f,%.2f,%.2f)\n",
	 vDIS.Pos(0),vDIS.Pos(1),vDIS.Pos(2),
	 vtxm.Pos(0),vtxm.Pos(1),vtxm.Pos(2));
	      }
#endif
	      continue;
	    }
	  }
	}
	if (nLs[C]) { // Lambda of same C-polarity already encountered...
  printf ("\n** MCInfo:\a Run#%d Evt %d: >1%sLambda (MCs %d,%d), 2nd ignored\n",
	  e.RunNum(),(int)e.UniqueEvNum(),C?" a-":" ",iMCT_L[C],imct);
          continue;  // ...too complex a case => let's ignore the 2nd Lambda
	}
	nLs[C]++; iMCT_L[C] = imct;
	if (trk.vMCvertex().size()>1) {
  printf("** MCInfo:\a Run %d Evt %d MC %d: Lambda w/ >1 decay vertices\n",
	  e.RunNum(),(int)e.UniqueEvNum(),imct); exit(1);
	}
	else if (trk.vMCvertex().size()==0) continue;
	const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
	const vector<int> &decayTrks = decayV.vMCtrack();
	if (decayTrks.size()!=2) continue;                  // ***** 2 DAUGHTERS
	int imct1 = decayTrks[0], imct2 = decayTrks[1];
	const PaMCtrack &decay1 = e.vMCtrack()[imct1];
	const PaMCtrack &decay2 = e.vMCtrack()[imct2];
	const TLorentzVector *lvp = NULL, *lvpi = NULL; double pp = 0, ppi = 0;
	const TLorentzVector lv1 = decay1.LzVec(), lv2 = decay2.LzVec();
	if (decay1.Pid()==14+C && decay2.Pid()==9-C) {           // ***** 1: ppi
	  lvp = &lv1;               lvpi = &lv2;
	  pp = decay1.Mom3().Mag(); ppi = decay2.Mom3().Mag();
	  iMCL_p[C] = imct1;        iMCL_pi[C] = imct2;
	}
	if (decay2.Pid()==14+C && decay1.Pid()==9-C) {           // ***** 2: pip
	  lvp = &lv2;               lvpi = &lv1;
	  pp = decay2.Mom3().Mag(); ppi = decay1.Mom3().Mag();
	  iMCL_p[C] = imct2;        iMCL_pi[C] = imct1;
	}
	if (!lvpi) continue;
	Types[iMCL_p[C]] |= 0x01000000; Types[iMCL_pi[C]] |= 0x02000000;
	PL[C] = trk.Mom3().Mag(); PL_p[C] = pp; PL_pi[C] = ppi;
	TLorentzVector lvL = *lvpi; lvL += *lvp; TLorentzVector lvp2 = *lvp;
	double pT = lvpi->Perp(lvL.Vect());
	if (fabs(lvL.M()-M_Lam)>.001) {   // X-check mass
 printf("** MCInfo:\a Run %d Evt %d MC %d: L Mass = %f\n",
	  e.RunNum(),(int)e.UniqueEvNum(),imct,lvL.M()); exit(1);
	}
	lvp2.Boost(-lvL.BoostVector()); // Boost p to L frame
	if (iMCT_mu!=-1) {
	  double cthstar = cos(lvp2.Angle(lvq.Vect())); cthstarL[C] = cthstar;
	  lvL.Boost(betacm);
	  // Maximal momentum is obtained w/ gamm*N -> Lambda K
	  // sqrt(s) = EL+EK = sqrt(mL^2+p^2)+sqrt(mK^2+p^2)
	  // A-2p^2 = 2sqrt(mL^2+p^2)(mK^2+p^2), A=s-mL^2-mK^2
	  // P^2 = A^2-4mK^2mL^2 / 4 s
	  double A = sqrts*sqrts-M2_Lam-M2_K, B = A*A-4*M2_Lam*M2_K;
	  if (B<=0) {
 printf("** MCInfo:\a y(=%.2f=>sqrt(s)=%.1f) too small for Lambda+K\n",
	    yB,sqrts);
            xFL[C] = -1;
	  }
	  else {
	    double pZMx = sqrt(B)/2/sqrts; xFL[C] = lvL.Vect().Dot(vqUnit)/pZMx;
	  }
	  Type |= 0x00100000; // Flag (a)Lambda -> ppi
	  bool winTarget = XpV*XpV+YpV*YpV<R2Cut;
	  if (Year==2012 || Year==2016 || Year==2017) {
	    winTarget |= PaAlgo::InTarget(XpV,YpV,ZpV,'O',PaEvent::Ptr()->RunNum());
	  }
	  if (Q2>1 && yLowCut<yB && yB<yUpCut && winTarget) {
	    // Set "Types" flag and fill histos
	    hc_aL[0][C]->Fill(cthstar); hc_aL[1][C]->Fill(pp);
	    int pType = Types[iMCL_p[C]];
#ifdef DEBUG_Lambda_RICH
	    static int nOKs = 0, nKOs = 0;
	    if (pType&0x10000000) nOKs++; else nKOs++;
	    printf("Evt#%d  %d %d  %d %6.1f %s\n",(int)e.UniqueEvNum(),nOKs,nKOs,iMCL_p[C],e.vMCtrack()[iMCL_p[C]].Mom3().Mag(),(pType&0x10000000)?"":"-> KO");
#endif
	    if (pT>V0pTCut) { // Have to add something about collinearity cut
	      Type |= 0x00400000; // Flag (a)Lambda w/in kin cuts
	      hc_aLk[0][C]->Fill(cthstar); hc_aLk[1][C]->Fill(pp);
	      if (pType&0x10000000) {
		hc_aLki[0][C]->Fill(cthstar); hc_aLki[1][C]->Fill(pp);
	      }
	    }
	    if (pType&0x10000000) {
	      hc_aLi[0][C]->Fill(cthstar); hc_aLi[1][C]->Fill(pp);
	    }
	    if ((Types[imct1]&0x00000060) && (Types[imct2]&0x00000060)) {
	      // Reconstructibility histos (denominator of efficiencies)
	      Type |= 0x0080;                 // ***** FLAG L AS RECO'IBLE *****
	      hc_RL[0][C]->Fill(cthstar); hc_RL[1][C]->Fill(pp);
	      if (Type&0x00400000) { // If also w/in kin cuts
		hc_RLk[0][C]->Fill(cthstar); hc_RLk[1][C]->Fill(pp);
		if (pType&0x10000000) {
		  hc_RLki[0][C]->Fill(cthstar); hc_RLki[1][C]->Fill(pp);
		}
	      }
	      if (pType&0x10000000) {            // ***** K is RICH'able *****
		hc_RLi[0][C]->Fill(cthstar); hc_RLi[1][C]->Fill(pp);
	      }
	    }
	  }
	}
	else {
 printf("** MCInfo:\a Run %d Evt %d: No scattered mu => no Lambda cos(theta*)\n",e.RunNum(),(int)e.UniqueEvNum());
           cthstarL[C] = 2;
	}
      }  // End (a)L block
      else if (iter==1 && (trk.Pid()==23 || trk.Pid()==31 ||  // ***** Xi-/+
			   trk.Pid()==24 || trk.Pid()==32)) { // ***** Omega-/+
	// ***** Sigma+/-*(Xi-/Omega-) -> Lpi+/-(pi-) RECONSTRUCTIBILITY

	int C = (trk.Pid()==23 || trk.Pid()==24) ? 0 : 1;
	iMCT_S[C] = imct;
	if (trk.vMCvertex().size()>1) {
 printf("** MCInfo:\a Run %d Evt %d MC %d: Xi/Omega w/ >1 decay vertices\n",
	e.RunNum(),(int)e.UniqueEvNum(),imct); continue;
	}
	else if (trk.vMCvertex().size()==0) continue;
	const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
	const vector<int> &decayTrks = decayV.vMCtrack();
	if (decayTrks.size()!=2) continue;                  // ***** 2 DAUGHTERS
	int imct1 = decayTrks[0], imct2 = decayTrks[1];
	const PaMCtrack &decay1 = e.vMCtrack()[imct1];
	const PaMCtrack &decay2 = e.vMCtrack()[imct2];
	if ((decay1.Pid()==18+8*C && (decay2.Pid()==8 || decay2.Pid()==9)) ||
	    (decay2.Pid()==18+8*C && (decay1.Pid()==8 || decay1.Pid()==9))) {

	  if (decay1.Pid()==18+8*C) {
	    iMCS_L[C] = imct1; iMCS_pi[C] = imct2;
	  }
	  else {
	    iMCS_L[C] = imct2; iMCS_pi[C] = imct1;
	  }
	  // ***** SET "Type" FLAG AND FILL HISTOS for Sigma* *****

	  if (Types[iMCS_L[C]]&0x00100000) {
	    // If daughter L (be it #1 or #2) decays -> Lpi
	    Types[iMCS_pi[C]] |= 0x04000000; // Flag pi as <- Sigma*
	    PS[C] = trk.Mom3().Mag();
	    double cthstar = cthstarL[C];
	    hc_aS[C]->Fill(cthstar);
	    bool recoibleS = (Types[iMCS_L[C]]&0x0080) &&
	      /* */          (Types[iMCS_pi[C]]&0x00000060);
	    if (recoibleS) {
	      Type |= 0x0080;                     // ***** FLAG Sigma* Reco'ible
	      int pType = Types[iMCL_p[C]];
	      if (pType&0x10000000 &&                 // ***** p is RICH'able...
		  (Types[iMCT_L[C]]&0x00400000))// ***** ...and Lambda Reco'ible
		hc_RSki[C]->Fill(cthstar);
	    }
	  }
	}  // End D* -> D0piS
      }  // End D0,D* blocks 
    }  // End loop on tracks
  }  // End loop on iterations 1 = D0->Kpi, 2 = D+*->D0pi
}
// ***************************************************************************
// ******************************   Bookgamma   ******************************
// ***************************************************************************
void MCInfo::Bookgamma()
{
  float zTarget = PaSetup::Ref().TargetCenterZ();

  hp_ag   = new TH1D("hp_ag",   "All #gamma vs. P",     100,0,100);
  hp_Ag   = new TH1D("hp_Ag",   "Accepted #gamma vs. P",100,0,100);
  hv_Ag   = new TH1D("hv_Ag",   "Acc'ed #gamma vs. Zconversion",1000,-200,800);
  hp_ACg  = new TH1D("hp_ACg",  "Acc'ed converted #gamma vs. P",100,0,100);
  hv_ACg  = new TH1D("hv_ACg",  "Acc'ed converted #gamma vs. Z",1000,-200,800);
  hp_Agee = new TH2D("hp_Agee", "#gamma#rightarrowe^{+}+e^{-} vs. P",
		     100,0,100,3,-.5,2.5);
  hp_Rgee = new TH2D("hp_Rgee", "Reco'ible #gamma#rightarrowe^{+}e^{-} vs. P",
		     100,0,100,3,-.5,2.5);
  hv_Rgee = new TH2D("hv_Rgee", "Reco'ible #gamma#rightarrowe^{+}e^{-} vs. Z",
		     1000,-200,800,3,-.5,2.5);

  hp_Rrgee= new TH2D("hp_Rrgee","Reco'ible&reco'd #gamma vs. P", 100,0,100,
		     3,-.5,2.5);
  hp_Irgee= new TH1D("hp_Irgee","!Reco'ible&reco'd #gamma vs. P",100,0,100);
  hp_Frgee= new TH1D("hp_Frgee","Fake&reco'd #gamma vs. P",      100,0,100);
  pe_gdpvsp[0] = new TProfile("pe_gdpvsp0","#gamma #Deltap/p vs. p",
			      100,0,100,"S");
  pe_gdpvsp[1] = new TProfile("pe_gdpvsp1","#gamma #Deltap/p vs. p",
			      100,0,100,"S");
  pe_gdpvsp[2] = new TProfile("pe_gdpvsp2","#gamma #Deltap/p vs. p",
			      100,0,100,"S");
  pe_gdpvsZ[0] = new TProfile("pe_gdpvsZ0","#gamma #Deltap/p vs. Zconversion",
			      40,zTarget-65,zTarget+135,"S");
  pe_gdpvsZ[1] = new TProfile("pe_gdpvsZ1","#gamma #Deltap/p vs. Zconversion",
			      40,zTarget-65,zTarget+135,"S");
  pe_gdpvsZ[2] = new TProfile("pe_gdpvsZ2","#gamma #Deltap/p vs. Zconversion",
			      40,zTarget-65,zTarget+135,"S");
}

void MCInfo::SortgammaOut(PaEvent &e)
{
#define REQUIRE_g_WIN_MM01
#ifdef REQUIRE_g_WIN_MM01
  const double zAcc = 142.0, hwAcc = 20;
#endif
  //#define REQUIRE_g_WIN_DC01  // This is assumed to correspond approx. to the acceptance of an ECAL
#ifdef REQUIRE_g_WIN_DC01
  const double zAcc = 260, hwAcc = 60;
#endif
  const double zMM02V = 190; // Just upstream of
  const double zSM1 = 350;
  const double zST03U1 = 540;  // Approx.

  iMCT_gamma.clear(); iMCgamma_e.clear(); iMCgamma_p.clear();
  gammaHasAvatar.clear();
  int imct; for (imct = 1; imct<e.NMCtrack(); imct++) {

    // ********** LOOP ON MC TRACKS, Excluding incident muon **********

    const PaMCtrack  &trk = e.vMCtrack()[imct];
    if (trk.Q()) continue;                         // ***** REQUIRE A NEUTRAL...
    if (trk.Pid()!=1) continue;                                      // ...gamma

    const PaMCvertex &vtx = e.vMCvertex()[trk.iVertexOfOrigin()];
    if (vtx.IsPileup()) continue;                          // ***** NOT A PILEUP

    //#define MC_DEBUG_gamma
#ifdef MC_DEBUG_gamma
    static int idebug = 0;
    if (idebug) {
      int itrk = vtx.iTrackOfOrigin();
      const PaMCtrack *mother = itrk==-1 ? 0 : &e.vMCtrack()[itrk];
      printf("(%.2f %.2f %.2f) <- Trk %d %s\n",
	     vtx.Pos(0),vtx.Pos(1),vtx.Pos(2),itrk,
	     itrk==-1 ? "incident" : mother->Name().c_str());
    }
#endif

    if (vtx.iTrackOfOrigin()!=-1) {
      const PaMCtrack &mother = e.vMCtrack()[vtx.iTrackOfOrigin()];
      if (mother.Pid()==2 || mother.Pid()==3) // ! an e^+ nor an e^-
	continue;                                          // ***** NOT A SHOWER
    }

    double pg = trk.Mom3().Mag(); hp_ag->Fill(pg);      // ***** ALL GOOD gammas

    //                                                          ***** ACCEPTANCE
    if (trk.Mom3()[2]<0) continue; // No backward pointing
    const PaTPar &helix = trk.ParInVtx();
    double z0 = helix(0), x0 = helix(1), y0 = helix(2), dx = helix(3), dy = helix(4);
    if (fabs(x0+dx*(zAcc-z0))>hwAcc || fabs(y0+dy*(zAcc-z0))>hwAcc)
      continue;
    hp_Ag->Fill(pg);

    if (trk.vMCvertex().size()>0) {                     // ***** DECAY -> e^+e^-
      const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
      double zDecay = decayV.Pos(2); hv_Ag->Fill(zDecay);
      const vector<int> &decayTrks = decayV.vMCtrack();
      if (decayTrks.size()!=2) {
	if (decayTrks.size()>2) {
	  printf("** MCInfo:\a Run %d Evt %d MC %d: gamma w/ %d decay tracks:",
		 e.RunNum(),(int)e.UniqueEvNum(),imct,(int)decayTrks.size());
	  for (int i = 0; i<(int)decayTrks.size(); i++)
	    printf(" %s",e.vMCtrack()[decayTrks[i]].Name().c_str());
	  printf("\n"); exit(1);
	}
	continue;
      }
      int imct1 = decayTrks[0], imct2 = decayTrks[1];
      const PaMCtrack *decays[2];
      decays[0] = &e.vMCtrack()[imct1]; decays[1] = &e.vMCtrack()[imct2];
      if (decays[0]->Pid()*decays[1]->Pid()!=6) {
	printf("** MCInfo:\a Run %d Evt %d MC %d: gamma decays = %s %s\n",
	       e.RunNum(),(int)e.UniqueEvNum(),imct,
	       decays[0]->Name().c_str(),decays[1]->Name().c_str());
	exit(1);
      }
      hv_ACg->Fill(zDecay); if (zDecay>zMM02V) continue;
      hp_ACg->Fill(pg);                     // ***** DECAY UPSTREAM of LAS *****
      int i, decaysOK; for (i = 0, decaysOK = 2; i<2; i++) {
	// Determine whether showering develops early enough that reco is not
	// possible, and what kind of reco it allows:
	//  - either reco of the original track that would have kept, to some
	//   extent, its integrity throughout:         "decaysOK=2"
	//  - or the reco of some avatar of the track: "decaysOK=0,1"
	//  - The 2 e^+ and e^- are successively considered
	//  - Onset of showering is caracterised by:
	//     - either a  5% loss in momentum in LAS upper ([MM02V:SM1]) arm
	//     - or     a 15% loss before ST03U,
	//  - Lesser bremsstrahlung being disregarded.
	//  - Higher loss (80% upper or 66%) => reco altogether impossible.
	const PaMCtrack *decay = decays[i];
	int pid = decay->Pid();
	double pei = decay->Mom3().Mag();  // Initial momentum
	double pep = pei, pec = pei;         // Updated mpmenta
	const vector<int> &imcvs = decay->vMCvertex();
	for (int iv = 0; iv<(int)imcvs.size(); iv++) {
	  const PaMCvertex &daughterV = e.vMCvertex()[imcvs[iv]];
	  double zDaughter = daughterV.Pos(2);
	  const vector<int> &daughterTrks = daughterV.vMCtrack();
	  for (int idght = 0; idght<(int)daughterTrks.size(); idght++) {
	    const PaMCtrack &daughter = e.vMCtrack()[daughterTrks[idght]];
	    if (daughter.Pid()==pid) {
	      pep = daughter.Mom3().Mag(); break;
	    }
	    if (daughter.Pid()!=1) {
	      printf("** MCInfo:\a Run %d Evt %d MC %d: %s ->:",
		 e.RunNum(),(int)e.UniqueEvNum(),imct,decay->Name().c_str());
	      for (int j = 0; j<(int)daughterTrks.size(); j++)
		printf(" %s",e.vMCtrack()[daughterTrks[j]].Name().c_str());
	      printf("\n"); decaysOK = -1; break;
	    }
	    pep -= daughter.Mom3().Mag();
	  }
	  if (zDaughter<zMM02V) {    // Bremsstrahlung before LAS entrance...
	    // ...In that case, reco'ibility poses no problem
	    pec = pep;                       // ...update initial momentum.
	  }
	  else if (zDaughter<zSM1) { // Bremsstrahlung w/in LAS upper arm...
	    if (pep<pec*.33) {             // ...loss > maximum maximorum=66%...
	      decaysOK = -1; break;            // ...=> no reco'ibility
	    }
	    if (pep<pec*.95) {             // ...loos > 5% but still reasonable
	      if (pep<pec*.66) decaysOK = 0;
	      else             decaysOK = 1;
	    }
	  }
	  else if (zDaughter>zST03U1) break;
	  else {                     // Bremsstrahl. w/in LAS downstream arm...
	    if (pep<pec*.20) {             // ...loss > maximum maximorum=50%...
	      // (This is assumed to imply too skewed a trajectory to be reco'd)
	      decaysOK = -1; break;            // ...=> no reco'ibility
	    }
	    if (pep<pec*.85) {             // ...loos > 15% but still reasonable
	      if (pep<pec*.50) decaysOK = 0;
	      else             decaysOK = 1;
	    }
	  }
	} // End loop over decay vertices for a given e^+ or e^- decay
#ifdef MC_DEBUG_gamma
	static int idebug = 0;
	int jdebug = idebug && 55<zDecay && zDecay<80;
	if (jdebug) {
	  if (i==0) printf("(%.2f,%.2f,%.2f) %d %s %.2f %.2f %.2f",
			   decayV.Pos(0),decayV.Pos(1),decayV.Pos(2),
			   imct1,decay->Name().c_str(),pei,pec,pep);
	  else      printf("%d %s %.2f %.2f %.2f\n",
			   imct2,decay->Name().c_str(),pei,pec,pep);
	}
#endif
	if (pec<.95*pei && decaysOK>=1) {// This corresponds to the case when
	  // bremsstrahlung loss occuring before entrance to LAS, summed up,
	  // goes beyond limit.   // => simply update "decaysOK"
	  if      (pep<pei*.66) decaysOK = 0;
	  else if (pep<pei*.95) decaysOK = 1;
	}
	if (decaysOK==-1) /* No wasting CPU on any further loop=> */ break;
      } // End loop over e^+/e^-
      if (decaysOK==-1) continue; int gRible = 2-decaysOK;
      hp_Agee->Fill(pg,gRible);      // ***** NO SHOWERING OR LATE ENOUGH in LAS
      int &Type = Types[imct]; Type = Types[imct1]&Types[imct2]&0xe0;
      if (!Type) continue;                  // ***** e^+e^ DECAYS RECONSTUCTIBLE
      hp_Rgee->Fill(pg,gRible); hv_Rgee->Fill(zDecay,gRible);
      iMCT_gamma.push_back(imct); gammaHasAvatar.push_back(gRible);
      if (decays[0]->Pid()==3) { // 1st is electron
	iMCgamma_e.push_back(imct1); iMCgamma_p.push_back(imct2);
      }
      else {
	iMCgamma_e.push_back(imct2); iMCgamma_p.push_back(imct1);
      }
    }
  }
}
void MCInfo::Fillgamma(PaEvent &e, int imctp, int imcte, double pg)
{
  int i, ig; for (i = 0, ig = -1; i<(int)iMCT_gamma.size(); i++) {
    if (iMCgamma_p[i]==imctp && iMCgamma_e[i]==imcte) { ig = i; break; }
  }
  if (ig<0) {
#ifdef MC_DEBUG_gamma
      int idebug = 0;
#endif
    int imctg, isg; for (i=isg = 0, imctg = -1; i<2; i++) {
      int imct = i?imcte:imctp; const PaMCtrack &trk = e.vMCtrack()[imct];
      int iv = trk.iVertexOfOrigin();
      const PaMCvertex &vtx = e.vMCvertex()[iv];
      int imctm = vtx.iTrackOfOrigin();
      const PaMCtrack *mother = imctm==-1 ? 0 : &e.vMCtrack()[imctm];
#ifdef MC_DEBUG_gamma
      if (idebug) {
	printf("%s MC%d <- %d (%.2f %.2f %.2f) <- Trk %d %s\n",
	       trk.Name().c_str(),imct,iv,vtx.Pos(0),vtx.Pos(1),vtx.Pos(2),
	       imctm,imctm==-1 ? "incident" : mother->Name().c_str());
      }
#endif
      if (i==0) imctg = imctm;
      else if (imctm!=imctg) imctg = -1;
      if (imctg>=0 && mother->Pid()==1) isg = 1;
      else                              isg = 0;
    }
#ifdef MC_DEBUG_gamma
    if (idebug) {
      if (imctg!=-1) {
	const PaMCtrack *mother = &e.vMCtrack()[imctg];
	int iv = mother->iVertexOfOrigin();
	const PaMCvertex &vtx = e.vMCvertex()[iv];
	int imct = vtx.iTrackOfOrigin();
	const PaMCtrack *gdmother = imct==-1 ? 0 : &e.vMCtrack()[imct];
	printf("g  : MC%d <- %d (%.2f %.2f %.2f) <- Trk %d %s\n",
	       imctg,iv,vtx.Pos(0),vtx.Pos(1),vtx.Pos(2),imct,
	       imct==-1 ? "incident" : gdmother->Name().c_str());
      }
    }
#endif	
    hp_Irgee->Fill(pg); if (!isg) hp_Frgee->Fill(pg);
    return;
  }
  const PaMCtrack  &trk = e.vMCtrack()[iMCT_gamma[ig]];
  double pmc = trk.Mom3().Mag();
  int hasAvatar = gammaHasAvatar[ig]; // One of e^+ or e^- is only partly (through one of its avatar) reco'ible 
  hp_Rrgee->Fill(pmc,hasAvatar); pe_gdpvsp[hasAvatar]->Fill(pmc,(pg-pmc)/pmc);
  const PaMCvertex &decayV = e.vMCvertex()[trk.vMCvertex().front()];
  pe_gdpvsZ[hasAvatar]->Fill(decayV.Pos(2),(pg-pmc)/pmc);
}

// ***************************************************************************
// ******************************    BookNpi    ******************************
// ***************************************************************************
void MCInfo::BookNpi(int Npis) {
  if (Npis!=3 && Npis!=5) {
    printf("** MCInfo::BookNpi: Bad argument =%d, neither = 3, nor = 5\n",
	   Npis); assert(false);
  }
  MCNpis = Npis;
  char hT[] = "3#pi mass  ";
  if (Npis==3)
    hf_Dalitz = new TH2D("hf_Dalitz","3#pi Dalitz",500,0,.2,50,0,.2);
  sprintf(hT,"%d#pi -t",Npis);   hf_tNpi = new TH1D("hf_tNpi",hT,200,0,1);
  sprintf(hT,"%d#pi mass",Npis); hf_mNpi = new TH1D("hf_mNpi",hT,100,0.5,4.5);
  sprintf(hT,"%d#pi Pv",Npis);   hf_pNpi = new TH1D("hf_pNpi",hT,200,0,200);
}
// ***************************************************************************
// ******************************   SortNpiOut  ******************************
// ***************************************************************************
void MCInfo::SortNpiOut(PaEvent &e) {
  int imct, im, ip, npips, npims; if (MCNpis==3) { npips = 1; npims = 2; }
  else {                                           npips = 2; npims = 3; }

  //    ***** INITIALIZE MCInfo 3/5pi DATA *****
  // (Note: PNpi = 0: means no P assigned, even if corresponding particle found,
  // iMCNpi != -1, due to the fact that no MChits (out of those output to mDST)
  //  could be found.)
  for (im = 0; im<npims; im++) { iMCNpi_pims[im] = -1; PMCNpi_pims[im] = 0; }
  for (ip = 0; ip<npips; ip++) { iMCNpi_pips[ip] = -1; PMCNpi_pips[ip] = 0; }

  TLorentzVector lvpims[3], lvpips[2];
  int run = e.RunNum(), event = (int)e.UniqueEvNum();

  const PaMCtrack &beam = e.vMCtrack()[0];
  if (beam.Pid()!=9 && beam.Pid()!=44 && beam.Pid()!=5) {
 printf("** SortNpi:\a run #%d, evt #%d: beam ID = %d != pi-(9), Prot. sp.(44) or mu-(5)\n",
	run,event,beam.Pid()); assert(false);
  }
  
  for (imct = 1, im=ip = 0; imct<e.NMCtrack(); imct++) {
    int &Type = Types[imct]; if (!(Type&0x001)) break; 

    //         ***** LOOP ON PRIMARY TRACKS *****

    const PaMCtrack  &trk = e.vMCtrack()[imct]; int id = trk.Pid();
    double P = 0;
    const set<int> &mcHits = trk.sMChitRef(); if (!mcHits.empty()) {
      const PaMChit &h = e.MChits()[*mcHits.begin()];
      double px = h.Px(), py = h.Py(), pz = h.Pz(); P = sqrt(px*px+py*py+pz*pz);
    }
    if      (id==8) {                                              // ***** PI+
      if (ip==npips) {
	printf("** SortNpi:\a run #%d, evt #%d: # of pi+ > %d\n",
	       run,event,ip); assert(false);
      }
      iMCNpi_pips[ip] = imct; PMCNpi_pips[ip] = P; lvpips[ip++] = trk.LzVec();
    }
    else if (id==9) {                                              // ***** PI-
      if (im==npims) {
	printf("** SortNpi:\a run #%d, evt #%d: # of pi- > %d\n",
	       run,event,ip); assert(false);
      }
      iMCNpi_pims[im] = imct; PMCNpi_pims[im] = P; lvpims[im++] = trk.LzVec();
    }
    else if (id!=14) { // Recoil p    
      printf("** SortNpi:\a run #%d, evt #%d: primary ! a pi+/- (id=%d)\n",
	     run,event,id); assert(false);
    }
  }
  if (im!=npims) {
    printf("** SortNpi:\a run #%d, evt #%d of pi- < %d\n",
	   run,event,npims); assert(false);
  }
  if (ip!=npips) {
    printf("** SortNpi:\a run #%d, evt #%d of pi+ < %d\n",
	   run,event,npips); assert(false);
  }

  TLorentzVector lvNpi = lvpips[0]; hf_pNpi->Fill(lvNpi.P());
  if (MCNpis==3) {
    TLorentzVector lvpippims[2];
    for (im = 0; im<2; im++ ) {
      TLorentzVector &lvpim = lvpims[im]; hf_pNpi->Fill(lvpim.P());
      lvpippims[im] = lvNpi+lvpim; lvNpi += lvpim;
    }
    mpippimMCNpi[0] = lvpippims[0].M2(); mpippimMCNpi[1] = lvpippims[1].M2(); 
    hf_Dalitz->Fill(mpippimMCNpi[0],mpippimMCNpi[1]);
    hf_Dalitz->Fill(mpippimMCNpi[1],mpippimMCNpi[0]);
  }
  else {
    for (im = 0; im<3; im++ ) {
      TLorentzVector &lvpim = lvpims[im]; hf_pNpi->Fill(lvpim.P());
      lvNpi += lvpim;
    }
    TLorentzVector &lvpip = lvpips[1]; hf_pNpi->Fill(lvpip.P());
    lvNpi += lvpips[1];
  }
  MMCNpi = lvNpi.M(); hf_mNpi->Fill(MMCNpi);

  // Beam may be a special particle (non interacting, designed to behave well
  // when propagated back in time) 
  // => compute its energy w/ correct pion mass assignment.
  TVector3 vbeam = beam.Mom3(); double pbeam = vbeam.Mag();
  double ebeam = sqrt(pbeam*pbeam+M2_pi);
  TLorentzVector lvt = TLorentzVector(-vbeam,ebeam); // Beam goes backward in time
  lvt -= lvNpi; // t = p1 - p3
  TMCNpi = lvt.M2(); hf_tNpi->Fill(-TMCNpi);

}
#include "TObject.h"
#include "TKey.h"
#include "TList.h"
#include "TDirectory.h"
#include "TAxis.h"
void setAxisTitle(const char *title)
{
  TList *l =  gDirectory->GetListOfKeys();
  TKey *key = (TKey*)l->First(), *oldkey = 0;  do {
    //   ***** SEARCH DIRECTORY CONTENTS FOR HISTOS "hf_*t" *****
    const char *name = key->GetName();
    // Keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),name)) continue;
    oldkey = key;
    if (!strncmp(name,"hf_",3) && strstr(name,"t")) {
      TH1D *hf = (TH1D*)gDirectory->Get(name);
      hf->GetXaxis()->SetTitle(title);
    }
  } while ((key = (TKey*)l->After(key))); // End loop on gDirectory contents
}

// $Log: MCInfo.cc,v $
// Revision 1.16  2020/07/30 10:46:34  ybedfer
// - Scattered mu:
// + Can be mu-.
// + "MC_HADRON"_SETUP": to not require a muon to be the scattered particle.
//  Note: better be retrieved from DataTakingDB and then as an arg.
// - "flagReinteraction": mother idetification: bug fixes.
// - "D0EffHistosAndFlags": Bug fix.
//
// Revision 1.14  2020/07/06 22:03:15  ybedfer
// - Merge v1.14 from CERN and 1.13 from CCIN2P3.
// - Bug fix in reco'ibility criterion for D*: piS has to be (L|S)AS reco'ible
//  (or if FF-reco'ible if FF not rejected !(pi_FF_REJECTION&0x2).
//
// Revision 1.14  2020/07/06 21:44:30  ybedfer
// - Set track <-> MCtrack X-reference also for incident particle.
// - Allow decays to not be ordered by increasing abscissa.
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
// Revision 1.13  2020/05/26 00:22:12  ybedfer
// Improved C++ programming.
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
//  - 180 -> 200 GeV.// 
//  - 3pi: accept mu as beam.
//
// Revision 1.8  2011/05/07 16:12:25  ybedfer
//  - Ctor
//   - Accepts new argument "int year".
//   - Initializes RICH efficiency matrices. (Yet "RICH_EFFICIENCY_VS_PS_2"
//    in "SortD0Out".
//  - 0x80000000; pi RICH misID'd as K -> ID'd (D0) or misID'd (Lambda).
//  - New histo "hf_pt"  = "Primary Tracks".
//  - "SortD0Out[_LUND]":
//    - Pions allowed to decay upstream of RICH.
//    - Reject tracks pointing backward.
//    - Upon "RICH_EFFICIENCY_VS_PS_2>=3", use tables initialized in ctor.
//  - New methods "Book3pi" and "Sort3piOut".
//
// Revision 1.7  2010/02/01 01:17:13  ybedfer
//  Fill c(theta*) and zD0 distributions, in addition to pD0 and pK.
//
// Revision 1.6  2009/12/08 11:21:50  ybedfer
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
// Revision 1.5  2009/12/03 00:49:56  ybedfer
//  - Bug fixes:
//    - Size of "nHitMapWs": 10 -> 13.
//    - "NHits" array: enforce size of 5 (in init.
//    - While setting "Type" flag for D*: set "iMCDS_piS" and not
//     "iMCS_pi" (bug which must have been introduced when implementing
//     MCInfo for the Lambda case).
//  - Abscissae of SM1, SM2 retrieved from PaSetup.
//  - delete -> delete[].
//  - "MC_HADRON_SETUP":
//    - MC index of scattered particle is = 1..
//
// Revision 1.4  2009/07/03 14:02:44  ybedfer
//  - Include "Lambda" analysis.
//  - More data members for event kinematics.
//  - RICH acceptance no longer X/Y symmetric.
//  - D0 analysis:
//    - Bug fix: type "iMCDS_pi" and not "iMCDS_piS" set if
//    - Some simplified programming in the D0 analysis.
//  - gamma analysis: obvious simplification.
//  - Notice that the routine is not able to cope w/ <0 muons!
//
// Revision 1.3  2009/06/07 11:18:39  ybedfer
//  Adding new member "Pmup".
//
// Revision 1.2  2008/04/16 18:10:29  ybedfer
//  - Redefinition of "Types":
//    - 0x00002000: mu+/-
//    - 0x00080000: MuID'ibility in MWA
//  - Smu reconstructibility: A 0x7 bit pattern (w/ MAs included) instead
//   of a 3-option choice. => Binning of kinematics histos.
//
// Revision 1.1  2004/02/12 17:07:43  ybedfer
// Initial revision
//

// $Id: UserEvent103.h,v 1.13 2021/04/05 13:27:15 ybedfer Exp $

// ******************** HISTOGRAMS for UserEvent103 ********************


// ***** The histogram booking is done according the following cpp defintions
//#define U3_FILL_Lambda
//#define U3_FILL_RICH_PERFS
//#define U3_GENERALPURPOSE_HISTOS
//#define U3_K0_HIGHLEVEL 	// High level = incidence related histos
//#define U3_ALLOUT_KINE   	// Enable all-out kinematics histos

//#define U3_PLOT_CHI2
//#define U3_PLOT_dTheta

//#define U3_DEBUG_W38_3
#ifdef U3_DEBUG_W38_3
#  define U3_ALLOUT_KINE
static TTree *tree;
#endif

//   ************************* GENERAL PURPOSE *************************
#ifdef U3_GENERALPURPOSE_HISTOS
static TH1D *h_triggers, *h_spills, *h_primary;
// Special histo of #delta(1/P_BMS) to monitor whether the cut on cov(1/P,1/P)
// used to reject beams w/o BMS is appropriate, for, over the course of the
// successive mass productions, it may have happened that the value assigned to
// cov(1/P,1/P) by default of BMS was badly set.
static TH1D *h_pBMS, *h_dqPBMS;
#endif
#ifdef U3_FILL_RICH_PERFS
//    ************************* RICH *************************
static TH2D *hR_dvsr[3];  // ***** dTheta vs. run *****
static TH1D *hR_dKPD[2][4];                                 // ***** Resolutions
static TProfile *hR_LHK[2], *hR_LHKpi[2], *hR_LHSpi[2];     // ***** LH Profiles
//    ****************************** phi ******************************
const int npB = 4; int pBins[npB+1] = {0,10,20,35,50}; // (p-KThr) bins
static TProfile *hR_LHKvsP[2][npB];
static TH2D *hR_dKvst[2][2];              // ***** Resolution vs. theta x angle
static TH2D *hR_nphvsp[2], *hR_chivsp[2];         // ***** #photons, chi^2 vs. p
static TH1D *hR_phiP[4], *hR_phiPK[4], *hR_phiAK, *hR_phiRK;
// ***** phi for RICH Eff *****
#  ifdef U3_PLOT_CHI2        // ********** Either RICH CHI2 PID... **********
static TH2D *hC_phi2, *hC_phiM2, *hC_phiMA2;     // ***** Chi2 Eff vs. Selection
#  else                      // ********** ...or  RICH LH PID **********
static TH2D *hL_phi2, *hL_phiM2, *hL_phiMA2;     // ***** Chi2 Eff vs. Selection
#  endif
static TH2D *hL_phiM[2];                              // ***** Effs w/in p Range
static TH2D *hL_phiPD[2][4], *hL_phiH[2][4];// ***** f(pseudoPD), f(Time in Day)
static TH2D *hL_phiS[2];                               // ***** Effs sub p Range
const int nphiBins = 4;                             // ***** As a f(Phase Space)
const int nphB = nphiBins*npB+1; // 1 ``sum'' bin  + nphBins * npB momentum bins
static TH2D *hL_phioP[2][nphB], *hL_phioPT[2][nphB], *hL_phioPV[2][nphB];
static TH1D *hR_phioP[2][nphB], *hR_phioA[2][nphB];
#  ifdef U3_PLOT_CHI2        // ********** RICH CHI2 PID **********
static TH2D *hC_phiM[2];                              // ***** Effs w/in p Range
static TH2D *hC_phiPD[2][4], *hC_phiH[2][4];// ***** f(pseudoPD), f(Time in Day)
static TH2D *hC_phiS[2];                               // ***** Effs sub p Range
static TH2D *hC_phioP[2][nphB], *hC_phioPT[2][nphB], *hC_phioPV[2][nphB];
#  endif
#  ifdef U3_PLOT_dTheta      // ********** RICH dTheta PID **********
static TH2D *hT_phiM[2];                              // ***** Effs w/in p Range
static TH2D *hT_phiPD[2][4];                                // ***** f(pseudoPD)
static TH2D *hT_phioP[2][nphB], *hT_phioPT[2][nphB], *hT_phioPV[2][nphB];
#  endif
//       ****************************** K0 ******************************
static TH1D *hR_dpiPD[2][16];                               // ***** Resolutions
static TProfile *hR_LHpi[2];                        // ***** Likelihoods Profile
static TH2D *hR_dpivst[2][2][2];     // ***** Resolution vs. theta x +/- x angle
//       *************** K0 for RICH Eff and Purity ***************
#  ifdef U3_PLOT_CHI2
static TH2D *hC_K02, *hC_K0M2, *hC_K0MA2;        // ***** Chi2 Eff vs. Selection
#  else
static TH2D *hL_K02, *hL_K0M2, *hL_K0MA2;          // ***** LH Eff vs. Selection
#  endif
static TH2D *hL_K0M[2], *hM_K0M[2];             // ***** Effs/misID w/in p Range
static TH2D *hM_K0S[2];                               // ***** misID sub p Range
static TH2D *hL_K0PD[2][4], *hM_K0PD[2][4];         // ***** (mis)ID f(pseudoPD)
static TH2D *hM_K0H[2][4];                               // misID f(Time in Day)
const int nK0Bins = 6;                              // ***** As a f(Phase Space)
const int nK0B = nK0Bins*npB+1;  // 1 ``sum'' bin  + nphBins * npB momentum bins
static TH2D *hL_K0oP[2][nK0B], *hL_K0oPV[2][nK0B], *hL_K0oPT[2][nK0B];
// Same sequence of eff vs. PS, but for K misID this time, different binning
static TH2D *hM_K0oP[2][nphB], *hM_K0oPV[2][nphB], *hM_K0oPT[2][nphB];
static TH1D *hR_K0oP[2][nphB], *hR_K0oA[2][nphB];
#  ifdef U3_PLOT_CHI2
static TH2D *hC_K0M[2], *hD_K0M[2];             // ***** Effs/misID w/in p Range
static TH2D *hD_K0S[2];                               // ***** misID sub p Range
static TH2D *hC_K0PD[2][4], *hD_K0PD[2][4];         // ***** (mis)ID f(pseudoPD)
static TH2D *hD_K0H[2][4];                               // misID f(Time in Day)
static TH2D *hC_K0oP[2][nK0B], *hC_K0oPV[2][nK0B], *hC_K0oPT[2][nK0B];
// Same sequence of eff vs. PS, but for K misID this time, different binning
static TH2D *hD_K0oP[2][nphB], *hD_K0oPV[2][nphB], *hD_K0oPT[2][nphB];
#  endif
#  ifdef U3_PLOT_dTheta
static TH2D *hT_K0PD[2][4], *hU_K0PD[2][4], *hT_K0M[2], *hU_K0M[2];
static TH2D *hT_K0oP[2][nK0B], *hT_K0oPV[2][nK0B], *hT_K0oPT[2][nK0B];
static TH2D *hU_K0oP[2][nphB], *hU_K0oPV[2][nphB], *hU_K0oPT[2][nphB];
#  endif
// ********** Lambda **********
static TH1D *hR_dpPD[4];                                    // ***** Resolutions
static TProfile *hR_LHp, *hR_LHppi;                // ***** Likelihoods Profiles
// ***** Lambda for RICH Eff *****
#  ifdef U3_PLOT_CHI2
static TH2D *hC_L, *hC_LM, *hC_LMA;              // ***** Chi2 Eff vs. Selection
#  else                     // Note: For LH Eff: consider also anti-Lambda 
static TH2D *hL_L[2], *hL_LM[2], *hL_LMA[2];       // ***** LH Eff vs. Selection
#  endif
#  ifdef U3_PLOT_dTheta
static TH2D *hT_LM;
#  endif
#endif

// ********** KINEMATICS **********
static TH1D *hQ2phi, *hyBphi;               // Kine histos for excl. phi
static TH2D *hQ2K0,  *hyBK0;                // Kine histos for K0
#ifdef U3_ALLOUT_KINE
static TH2D *hX_T, *hY_T, *hXp_T, *hYp_T, *hqoP_T; // mu' track pars vs. trig
static TH2D *hPar_Ts[5];
static TH2D *hQ2_T, *hxB_T, *hyB_T;
static TH1D *hQ2,   *hxB,   *hyB;
static TH1D *hQ2h1, *hxBh1, *hyBh1;
static TH1D *hQ2h2, *hxBh2, *hyBh2;

// *************** VERTICES ***************
static TH1D *hv_pZ, *hv_pZT, *hv_sZ;// Z of pV (w/in Target), of 2ndary vertices
static TH2D *hv_pXY, *hv_pXYT;      // X,Y of pV (w/in Target)
static TH1D *hv_sZV0;               // Z of sVertex conditioned by V0
static TH1D *hv_pZh1, *hv_pZh2;     // Z of pV for pT>.7, SpT2>2.5
static TH1D *hv_pZQ2,  *hv_pZDI;    // Z of pV for Q2>.1, Q2>1
static TH1D *hv_pZQ2h, *hv_pZDIh;   // Z of pV w/ 1h for Q2>.1, Q2>1 (Golden)
static TH1D *hv_pZGMP, *hv_pZGDC;   // Golden w/ only MPs, w/ DC01 
static TH2D *hv_XYQ2,  *hv_XYDI;    // XY of pV for Q2>.1, Q2>1
static TH2D *hv_XYQ2h, *hv_XYDIh;   // XY of pV w/ 1h for Q2>.1, Q2>1
static TH2D *hv_FIQ2s[3],  *hv_FIDIs[3];        // Same for vertices in FI03s
static TH2D *hv_FIQ2hs[3], *hv_FIDIhs[3];
static TH2D *hv_MPQ2s[2],  *hv_MPDIs[2];        // Same for vertices in MP01UVs
static TH2D *hv_MPQ2hs[2], *hv_MPDIhs[2];
static TH2D *hv_KdRvsR;             // Beam dRdZ vs. R @ pVertex
#endif
//       ************************* phi *************************
// ********** phi Selection **********
static TH1D *hk_dE, *hk_dEK, *hk_dEKphi;          // Exclusivity
static TH1D *hm_phiAll, *hm_phi, *hm_phiID;     // phi: w/o|w pT cut, w/o|w ID.
static TH1D *hm_rhoAll, *hm_rho;                // rho: w/o|w pT cut
static TH1D *hk_dEW, *hm_phiIDW;                // Same sign pairs
static TH1D *hk_dEKcphi, *hm_phiIDc;            // ELoss correction
static TH2D *hm_phiVsID;                        // phi w/ K+/-ID 
static TH1D *hk_pTphi;
static TH1D *hm_phill, *hm_phisl, *hm_phiss;    // Excl. K+K- in LAS, S+LAS, SAS
static TH1D *hm_phiLL, *hm_phiHH;  // phi for Low (resp. H) K momenta
static TH2D *hk_Xphi, *hk_Yphi;                 // @ RICH
static TH2D *hk_PThphi;
static TH2D *hR_thCPphi, *hR_ThCPphi;           // ThC vs. P (uncorr'd/cor'd)
static TH2D *hR_THCPphi;                        // ThC-ThCK
static TH2D *hm_iphibb, *hm_iphiab, *hm_iphiaa; // Incl. phi sup-/sub-KThr, w/ K0
static TH1D *hm_Iphibb, *hm_Iphiab, *hm_Iphiaa; // Incl. phi sup-/sub-KThr
// *************** K0 phi ***************
static TH1D *hm_K0phi, *hm_K0phiC, *hm_K0phiA, *hm_K0phiAC;
static TH1D *hS_K0phi, *hS_K0phiC;
static TH1D *hs_K0phipi, *hs_K0phipiC;
// ***** KINEMATICAL CUT *****
static TH1D *hk_tK0phi, *hk_tKSpi[2];


//    ************************* V0 Selection *************************

// ***** LENGTH, ANGLE of VERTEX LINE, V CHI2, Pt, TIME, Theta *****
static TH1D *hs_dK0, *hs_dL, *hs_aK0, *hs_aL, *hs_vK0, *hs_vL, *hs_pK0, *hs_pL, *hs_tK0, *hs_tL;
static TH1D *hn_dK0, *hn_dL, *hn_aK0, *hn_aL, *hn_vK0, *hn_vL, *hn_pK0, *hn_pL, *hn_tK0, *hn_tL;
#ifdef U3_K0_HIGHLEVEL
static TH1D *hs_TK0, *hn_TK0;
#endif

//            ******************** K0 ********************
// ********** K0 MASS CONDITIONED by Z of 2NDARY VERTEX **********
static TH1D *hm_K0a, *hm_K0, *hm_K0t, *hm_K0ll, *hm_K0sl, *hm_K0ss;
// *************** K0 MASS conditioned by VERTEX LINE ***************
static TH1D *hm_K0c11, *hm_K0c12, *hm_K0c21, *hm_K0c22;
static TH1D *hm_K0cll, *hm_K0csl, *hm_K0css;  // In LAS^2, LAS*(L)SAS, (L)SAS^2
static TH1D *hm_K0eV;               // eVeto
static TH2D *hm_K0VsID;             // Vs. PID

#ifdef U3_K0_HIGHLEVEL
static TH1D *hm_K0cSSp, *hm_K0cSSm; // |Ppi|>40GeV
static TH1D *hm_K0cSSP, *hm_K0cSSM; // |Ppi|>50GeV
static TH2D *hm_K0cllT, *hm_K0cssT;           // In LAS^2/(L)SAS^2 vs. Theta
static TH2D *hm_K0cllx, *hm_K0cssx;           // In LAS^2/(L)SAS^2 vs. Theta x
static TH2D *hm_K0cslx, *hm_K0clsx;           // In LAS^2/(L)SAS^2 vs. Theta x
static TH2D *hm_K0clly, *hm_K0cssy;           // In LAS^2/(L)SAS^2 vs. Theta y
static TH2D *hm_K0cllX, *hm_K0cssX;           // In LAS^2/(L)SAS^2 vs. Theta X
static TH2D *hm_K0clla, *hm_K0cssa;           // In LAS^2/(L)SAS^2 vs. alpha
static TH2D *hm_K0cllb, *hm_K0cssb;           // In LAS^2/(L)SAS^2 vs. beta
static TH1D *hm_K0cnn;                        // No LAS: tracking downstream SM1
static TH1D *hm_K0cNNp, *hm_K0cNNm;           // No LAS, w/ |Ppi|>40GeV
static TH2D *hm_K0cllP, *hm_K0cssP;
static TH1D *hm_K0cssS, *hm_K0cssL;
#endif

static TH1D *hm_K0C22;                 // Requiring mu/mu'
static TH1D *hm_K0i, *hm_K0o;
static TH2D *hk_XK0, *hk_YK0;          // @ RICH
static TH2D *hk_PK0, *hk_PThK0;
static TH2D *hR_thCPK0, *hR_ThCPK0;    // ThC vs. P (w/o, w/ index correction)
static TH2D *hR_THCPK0;                // ThC-ThCpi
// *************** K0 pi ***************
static TH1D *hm_K0pi, *hm_K0piC;
// *************** K0 K ***************
static TH1D *hm_K0K[2], *hm_K0KC[2];
// *************** K0 p ***************
static TH1D *hm_K0p[2], *hm_K0pC[2], *hm_K0pCM[2], *hc_missM2;
// *************** K*pi ***************
static TH1D *hm_KSpi[2], *hm_KSpiC[2], *hm_KSpit[2], *hm_KSpiT[2];
static TH1D *hS_KSpi[4], *hS_KSpiC[4], *hS_KSpit[4], *hS_KSpiT[4];
static TH1D *hm_KSpipi[4], *hm_KSpipiC[4];

#ifdef U3_FILL_Lambda
//           ******************** Lambda ********************
static TH1D *hk_PpL[4], *hk_PpiL[4];// Decay momenta, [a]Lambda, w/ p,piID
static TH2D *hk_XL[2], *hk_YL[2];   // p,pi @ RICH
static TH2D *hk_PThL[2];
static TH2D *hm_LVsID[2];           // Mass vs. PID
static TH2D *hR_thCPL, *hR_ThCPL;   // ThC vs. P (w/o, w/ index correction)
static TH2D *hR_THCPL;              // ThC-ThCp
static TH1D *hm_Lpi[4];             // Sigma* -> [a]Lambda pi+/-
static TH1D *hm_LpiID[4];           // Sigma* -> [a]Lambda pi+/-, w/ piID
static TH1D *hm_LID[2];             // [a]Lambda for Sigma* search
static TH2D *hm_LcthS[5], *hm_aLcthS[5];  // ctheta* vs. mass vs. as a f(xF) 
static TH1D *hQ2_L[2], *hxB_L[2], *hyB_L[2], *hxF_L[2];
static TH1D *hv_pZ_L[2];
// ***** CASCADES
static TH2D *hm_LpiCa[2];           // Ksi   -> [a]Lambda pi-(pi+), w/ piID
static TH2D *hm_LpiCa3[2], *hm_LpiCa6[2];        // D and cth cascade cuts
static TH2D *hm_LpiCaD[2], *hm_LpiCaC[2];        // D/dD and cth cascade cuts
static TH2D *hm_LKCa[2];            // Omega -> [a]Lambda K- (K+), w/ KID
static TH2D *hm_LKCa3[2],  *hm_LKCa6[2];         // D and cth cascade cuts
static TH2D *hm_LKCaD[2],  *hm_LKCaC[2];         // D/dD and cth cascade cuts
static TH1D *hs_d3hXi, *hs_D3hXi, *hn_d3hXi, *hn_D3hXi;
static TH1D *hs_dpcXi, *hs_DpcXi, *hn_dpcXi, *hn_DpcXi;
static TH1D *hs_cpcXi, *hs_CpcXi, *hn_cpcXi, *hn_CpcXi;
static TH1D *hs_dcsXi, *hs_DcsXi, *hn_dcsXi, *hn_DcsXi;
static TH1D *hs_cpsLa, *hs_CpsLa, *hn_cpsLa, *hn_CpsLa;
#endif

//#define UserEvent_NTuple 3
#ifdef UserEvent_NTuple
  static TNtuple *ntDth;
#endif

// $Log: UserEvent103.h,v $
// Revision 1.13  2021/04/05 13:27:15  ybedfer
// No change.
//
// Revision 1.12  2021/03/31 17:34:59  ybedfer
// - "fillLambdaX" selection modified.
// - "book(K0|Exclusive)": more correct titles.
// - V0 cuts:
// + "bookLambdaMisc": array of "cth" and "pT" cuts.
// + V0DsdD[2] now used, w/ tightened value.
// - theta vs. P:
// + Tighter "cth/pT" cuts to "hR_thC(PK0|L)" histos and the like,...
// + No LH/LHBck cut.
// + "hR_THCP" of diff. w.r.t. theory for mean index.
// - RICH perf. for Lambda: No more impact on the following.
//
// Revision 1.11  2020/10/14 14:26:31  ybedfer
// - "U3_uDST_PRESELECTION":
// + Loose "(E|I)phi_pTCUT".
// + Looser "V0pTCut".
// - Exclusivity:
// + "fillExclusive" sets "KEMiss" status flag.
// + Select/Histo rho's, w/ its piID.
// - Incl. phi:
// + Require w/in target zone.
// + Better enforce non-exclusivity: even if bad BMS (based on "KEMiss" flag).
// - Cascades:
// + "fillLCascades": Reshuffled, w/ extra "ppiIDs" arg.
// + "subpiThrLHeVeto": loose eVeto below pTHr.
// + Dedicated "bookLCascades".
// - "K0LambdaphiXi_SELECTION" -> "K0LambdaphiXi_WRITEBACK".
// - Tree:
// + Cancel "nK0s". "nSVtx" -> "nOuts".
// + "fPa2fHs": Prevent "CSHadronData" from entering twice, when in two resonances.
// - "bookEphi" renamed "bookExclusive".
// - Binning Z of vertex: "getFIMMAbscissae", "getTargetBinning".
// - "WinTargetZone": extra X,Y arg.
//
// Revision 1.10  2020/07/30 10:20:58  ybedfer
// - "fillIphi": Bug fixes:
// + pT cut is minimum pT cut.
// + Update "phiPat" prior to filling TTree.
// - "U3_TRIGGER": built-in setting. (Note: better be retrieved from DataTakingDB.)
// - Upgrade to CSEventData,v1.4.
// + "copyCSHadronData":
// ++ Add "Px,Py,Pz".
// ++ New arg. list: "P" (cancelled: no re-scaling) -> "iV".
// + "copyCSResonanceData": New "alpha" arg.
// - "MuPrim_MASK": Bug fix: cpp-based setting had not been taken into account.
// - "aRs":
// + Skip tracks starting dowstream of RICH.
// + Take advantage of available smooth helices.
// + =-2 <-> downstream of RICH. Double check ">=-1".
// - RICH indices:
// + Fundamental one is UV = APV.
// + Histo thC vs. P w/o and w/ index evolution.
// + Histo'ing: tighter safety margin above threshold.
// - "allTrigs": shorter -> 0xffff.
// - Target cuts:
// + Bug fix: if "TARGET_CUT" =5, looser primary requirement.
// + In "primaryPat", account for "DataTakingDB::WinTarget" being unsigned int.
// + Kinematics: histo "Yp" vs. "Xp", w/o and w/ full target cuts.
// - "fillEphi_RICHPerf": "pRange": bug fix?
// - Rescaling: simplify further w.r.t. v1.9.
// - "setPeriodDep": New arg. list: "year" -> "e".
// - Update a few global cpp options, but not "uDST_PRESELECTION".
// - "muonBeam" flag, replacing "dpBMS".
// - Globals: add "Run,EvNum,Year".
// - Some, but no systematics, sprintf -> snprintf for histo titles.
//
// Revision 1.8  2018/04/26 08:34:21  ybedfer
// - Inclusive phi: not sure what's the situation:
// + Requirements = "(primaryPat&0x50)==0x50" and w/o "nTrks>3" condition:
//  (afair) after Quiela had complained that she did not see the expected
//  amount of inclusive phi's in her tree.
// - Tree:
// + Align w/ "CSEventData.h,v1.3".
// - cpp #idef:
// + "U3_TARGET_y_CUTS" => "targetOption".
// + "CHI2_CUT" -> "CHI2_CUT_OPT".
// + "BMS_CHI2_CUT" -> "BMS_LH_CUT".
// - "DataTakingDB":
// + "WinTargetZone".
// + "HadronPhysics" -> "(Hadron|Primakoff_mu|DY)".
// + "PBMS" -> "PBeam", "Y(Low|Up)Cut" -> "y(Low|Up)Cut".
// + "dpBMS" is now checked to ascertain beam is muon.
// - TLorentzVector's "lvk(i|s)" in "parseBestPV/fill(K0|Lambda)X".
// - "fillExclusive":
// + "parseVertex" -> "checkExclusive"
// + "fillK0X": more parentheses.
// - Kinematics: vertex distributions:
// + New binning. More histos (Q2>.1 and Q2>1, w/ and w/o hadron required).
// + "U3_FIMP_VTCS"
// - "getRICHMultiplicity".
// - "getTrigType".
//
// Revision 1.7  2013/05/28 17:50:46  ybedfer
// - Introduce TTree output"U3_OUTPUT_TREE".
//   - Count #sVertices, #K0s, #Ks ("nChargedKsInPV").
// - Overall:
//   - Fringe field rejection: "zone==0x1" -> "zone&0x2" => accounts for BoT.
//   - Variable "thrMargin" replacing constants 1.05(KThr), 1.1(pThr).
//   - All resonances R = K0|phi|L: "hm_RVsID", "hk_PThR", "hk_XR" histos.
//     - Z of V0 vertices: to be able to spot fake V0's from re-interaction.
// - phi:
//   - More requirements on event's pVertex: "primaryPat&0x10" -> "0x32".
//   - Relaxed pT cut (now passed by arg. to "fillExclusive"): 50 -> 20 MeV.
//   - "U3_QUASIEXCL_PHI": allow one extra track for exclusivity selection.
//   - New subroutines: "parseVertex,detachedHadrons,nInTimeNeutrals".
//   - Dedicated "fillIphi", conditioned by strict "primaryPat"
//    => Inclusive phi's: now w/o coinc. w/ K0; saved to uDST.
// - V0:
//   - "U3_V0_CUT: V0DsdD/V0cthCut[3]":
//     - Simplified and relaxed, on the ground that now, w/ the TTree output,
//       we can refined the cut later on...
//     - Tighter option for "U3_uDST_ANALYSIS" (w/ output to TTree).
//   - Vertex chi2 -> chi2/NDF => Since cut value hasn't been updated, it is
//    now exceedingly loose.
//   - Early requirement on mass to be w/in histo range, and, for K0, on
//    D/dD and theta_c:
//     - Should have no impact.
//     - Previously these cuts had been applied to the uDST selection.
//   - PID: fully reshuffled.
//   - "piS_e_REJECTION"
//     - Evaluate PID in any case, even if "piS_e_REJECTION" is not defined.
//     - Cancel "piS_e_ALL_REJECTION".
//     - "eLH>LHeVeto*LHBckCut" dropped. "eLH>piLHeVeto*piLH" only applied.
//   - In/outside target (definition of "cleanV0", histo): setup dependent.
//   - "U3_K0/L_HIGHLEVEL" (for dedicated mass shifts studies)
//     - Disabled by default.
//     - Now a dedicated "fillK0HighLevel" routine.
//   - Cut on "primaryPat" moved upstream. (This should have no implication
//    and is mentioned for the record.)
//   - "U3_ERASE_K0s":
//     - Looser cut on s_dist,s_theta: >=2 -> >=1.
//     - More consistent set of histos disabled.
//   - Histo of Z abscissa of sV: to visualize fake V0 vs. re-interaction.
// - Lambda:
//   - First block w/ loosest pT cut: should not affect uDST selection.
//   - Maybe some more changes.
// - RICH
//   - "RICH_TUBE_RADIUS": default now a more reasonable "4" (instead of "1").
//   - "PS_USE_RADIUS" not default: do not remember why.
//   - Index of Cherenkov angle estimator: now a variable (<-"DataTakingDB").
//   - Threshold cuts now made on P from tracking, in principle, more
//    appropriate than on P from vertex.
//   - RICH_perfs: Extrapolation to RICH entrance is one one in "transport".
//   - "transport":
//     - Require track ZFirst to be upstream of RICH.
//     - (Not clear) Bug fix: cutting now of RRICH and not RRICH^2.
// - BMS:
//   - "primaryPat" 0x2 bit (pV beam w/in range): reject makeshift BMS.
//    (Corresponding cut on beam P uncertainty is set in "DataTakingDB".)
//   - Histogram P and dP (upon "U3_GENERALPURPOSE_HISTOS").
// - Bad spills under UserEvent control: "U3_BAD_SPILLS".
//   - Bad spill files are to found in "~/phast.utils/badSpills/<year>".
//   - W/ syntax like "BadSpill_06W32_2.dat" or "BadSpill_11W30_1.RICH.dat".
//   - Remember previous spill#: "static int Spill".
// - Flowchart:
//   - Arrays of track/vertex attributes: instantiation moved downstream,
//    but still done in any case.
//   - All cpp definitions move to the .cc. The inclusion of the header
//    ("UserEvent3.h") being delayed until all relevant definitions are
//    done.
//   - Always require pV => no longer condition V0 w/ pV by "primarypat".
//   - "parseBestPV" no returns also beam energy "E0".
//   - All extrapolations done in "transport" (stored in "PaTrack::vAux").
// - CPP macros:
//   - "(E|I)_phi_pTCUT" -> "(E|I)phi_pTCUT".
//   - Remove some GLOBAL definitions: "U3_V0_WO_PID, CHI2_REJECTION,
//    U3_DEBUG_e_REJECTION, U3_Lambda_DEBUG,...".
//   - Introduce some consistency checks (mutual exclusion) in GLOBAL CPP'S.
//   - Check (partially) that only one GLOBAL cpp definition is enabled.
// - Bug fixes:
//   - BMS P spread was wrongly retrieved from "DataTakingDB".
//   - Assignments of "p|K|pThr" were comma-separated.
//   - "trigType" now made depending upon data taking year (but this remains
//    to be moved to "DataTakingDB").
//   - "mom1/2" is re-scaled (for when non-Marcin re-scaling enabled).
//
// Revision 1.6  2012/11/05 21:16:35  ybedfer
// - In general, try to systemtically store the outcome of cuts into flags (
//  to be later copied to output TTree (not yet done)) rather than exit right
//  away. Example of such flag: "primaryPat".
// - Introduce (and select by default) Marcin's re-scaling of SM2.
// - Interfacing w/ "DataTakingDB", from which per year/period defaults are
//  taken.
// - "InitPIDs" replaced by simpler "GetPIDSs", which no longer bother about
//  mu-ID.
// - "SelectBeam" is dropped: now resort to more standard methods to deal
//  w/ bad beam rejection.
// - pID specifics ("pIDPmax, pIDPmin, subpThrPmin"): simplified.
// - More modularity: added "parseBestPV, bookKinematics, fillAllOutKine,
//  bookV0Selection, fillK0woPV, get/fillK0Incidence, fillK0X, bookLambdaMisc,
//  fillLambdaX, fillLCascades, get/fillK/L0GM0709".
// - Histos now created in a hierarchy of sub-TDirectory's.
// - Excl. phi: now also plotting/selecting same sign pairs.
//
// Revision 1.5  2009/07/03 13:49:08  ybedfer
//  - "fillExclusive" in an independent routine and simplified: no independent
//   histo of "dE" upon mu'ID type.
//
// Revision 1.4  2009/06/07 00:51:56  ybedfer
//  Revison corresponding to "UserEvent3.cc,v1.4".
//
// Revision 1.1  2008/04/17 20:28:58  ybedfer
// Initial revision
//
// Revision 1.2  2004/05/20 23:09:25  ybedfer
// ...
//

// $Id: UserEvent103.cc,v 1.13 2021/04/05 13:26:00 ybedfer Exp $

// K0, Lambda, phi, Xi (Lambda cascades)
//   I) Selection ("K0LambdaphiXi_SELECTION") of these channels for writing
//     back microDSTs (uDSTs).
//  II) Histogramming documenting the selection: they may include some cuts
//     that do not affect the output uDST.
//      The optional ("U3_OUTPUT_TREE") output TTree doesn't include these cuts
//     either, but keeps track of them via a series of per resonance flags.
// III) Application to the evaluation of RICH Efficiency and Purity
//  IV) Search for K0+X and Lambda+X resonances

// (I) The selection (cf. "uDSTSelection" in the source code) is based on:
// - Mass range cuts (cf. "dM_K0|Lambda|phi", and there are also cuts for
//  cascades Xsi and Omega.)
// - pT cut.
// - K0/Lambda: V0 cuts, viz.: VV distance and collinearity. Cf.
//             "s_dist>=1 && s_ctheta>=1". 
//              sVertex chi2 cut, but it is exceedingly loose.
// - Excl. phi: mu+mu'+2*h events, w/ cut on inelasticity, cf. "|dEK|<dECut".
//             (There is provision for an extra, slow h, upon option.)
//   Incl. phi: one of the decay K's is KID'd (either below or above threshold).
// - Cascades : No Lambda collinearity, mass cut on Lambda, cut on CDA of
//             bachelor, Dist and coll. pV/cascade, Dist Cascade/Lambda.
// - In any case: exclude fringe-field (but keep them in the inventory of
//  tracks on which exclusive events are singled out for the phi selection). 
// The selection process is documented in a number of histos, showing:
// - the selected events in their mass distributions,
// - in addition, those events that are retained by looser and tighter cuts,
// - various kinematical distributions.
//  The cpp macro "U3_ERASE_K0s" removes from the events the decay pions of
// K0's selected w/ looser cut. It may be helpful when used in conjunction w/
// the D0(->Kpi)-dedicated "UserEvent5", to reduce the combinatorics.

// (III) The default evaluation is that based on likelihoods.
//      Methods relying on Chi2 or Cherenkov angle can also be evaluated, cf.
//     "U3_USE_CHCHI2", "U3_PLOT_dTheta". Note that they are now obsolete.

// Bad Spills: they can be read in and saved into output TTree, rather than
// supplied to the PHAST command line. This provides increased versatility.

// ***** WORKFLOW
// - EVENT SELECTION SWITCH
// - EXCLUSIVE PHI
// - V0 SELECTION
//   - FIRST PARTICLE = POSITIVE

// NOTABENE: CPP DEFINITIONS
//  The source code comes in several flavours, depending upon cpp macro
// definitions ("#define"). See "UserEvent103 VERSIONS" infra.
//  For convenience, some GLOBAL CPP DEFINITIONS recapitulate the settings
// recommended for some particular application (e.g. "U3_uDST_PRESELECTION").
//  They are listed at the top (of the sequence of cpp definitions). But are 
// specified at the bottom (so that they overwrite all intervening definitions).
// Beware still that some of these definitions must be made consistent between
// this "UserEvent103.cc" file and the header files included therein.

// Dependencies:
// - UserEvent103.h
// - Masses.h, GetPIDs.cc, RICH.cc|h, MCInfo.cc|h, DataTakingDB.cc|h,
//  HistoLambda.cc, MCLambda.cc|h, UsParticle.h

#include <math.h>
#include <iostream>
#include <cmath>

#include "TVectorD.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"

#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "PaPid.h"
#include "PaUtils.h"
#include "G3part.h"

#include "RICH.h"
#include "MCInfo.h"
#include "MCLambda.h"
#include "Masses.h"
#include "UsParticle.h"
#include "DataTakingDB.h"

// *************************************************************************
// *************************  UserEvent103 VERSIONS  *************************
// *************************************************************************

//       ********** LIST OF *GLOBAL* CPP DEFINITIONS **********
#define 	U3_uDST_PRESELECTION	// Default GOBAL CPP: all purpose K0Lambdaphi preselection
//#define 	U3_uDST_ANALYSIS 	// Re-analyzing preselected uDSTs w/ tighter cuts
//#define 	U3_Lambda_uDST_PRESEL	// Dedicated selection: aimed at Lambda physics
//#define 	U3_Lambda_ANALYSIS
//#define 	U3_L_MC_ANALYSIS
//#define 	U3_ERASE_K0s	// Single out K0's and erase their pion decays: to be used in conjunction w/ another UserEvent. 

//     ********** CPP OPTIONS INDEPENDENT of GLOBALS **********
#define BESTpV_OPTION 1 	// =1: CORAL if newerCoral, =2: CORAL in any case, else: PHAST
#define MuPrim_OPTION 1 		// =1: PHAST, =2: coral's mu'-ID
// One can overwrite the default PaVertex::iMuPrim set of arguments:
// - Default = 0xd, w/ bit 0x1 possibly if prescribed by "DataTakingDB".
// - Here we do not aim at DIS => No "reject2MuEvents", nor "checkCanBeMuon".
#define MuPrim_MASK 8
// Momentum ReScaling
// =1: Use Marcin's "../contrib/RescaleMom.cc" in 2006
// !!NOTA BENE!!: DO NOT enable Marcin on microDSTs !!
#define SCALE_MOMENTA
//#U3_SPECIAL_SCHEME 2008	// Special features can be enabled this way (see source code)
//#define U3_Lambda_LOOSE 	// For debugging
//#define U3_USE_CHCHI2		// Use Cherenkov chi2 instead of LH
//#define U3_L_HIGHLEVEL  	// Mass vs. incidence, cf. "./HistoLambda.cc"
//#define U3_BAD_SPILLS   	// Bad spills under the control of UserEvent103

// ##### HISTOGRAMMING
#define U3_GENERALPURPOSE_HISTOS
#define U3_ALLOUT_KINE		// Extra histos, showing input data kinematics
#define U3_FILL_RICH_PERFS
#define U3_FILL_Lambda
#define U3_K0_HIGHLEVEL 	// High level = incidence related histos

//      ********** CPP OPTIONS THAT MAY BE RESET BY GLOBALS **********
// - IMPACT uDST SELECTION
// Write-back can be disabled on phast command line. Yet, in case UE103 is
// executed in parallel w/ other UE's, it may be useful to disable it from w/in.
#define K0LambdaphiXi_WRITEBACK // Enable write back of selected events
// Trigger: TODO: Get it from DataTakingDB.
#define U3_TRIGGER 0x21e	// BUILT-IN! Corresponding to 2016 
// Chi2 cut: it's taken from DataTakingDB
#define CHI2_CUT_OPT 1		// =1: Loose cut, =2: Period dependent cut
// Target cuts:
// =0: No cuts.
// =1: Standard: PaAlgo::InTarget
// =2: Home made
// =3: DISTarget if relevant (i.e. MC w/ SVN# available), else as for =1
// =4: Between bT and spectro, or (MP01U+MP01X)/2 in 2016/17
// =5: As =4 but then use =1 to document TTree
#define U3_TARGET_CUT 5
// y cut:
//  Since we are not here interested in DIS (where e.g., y=1 region can be
// interesting), let's have the basic DataTakingDB cut: something to guarantee
// a reliable reconstruction.
#define U3_y_CUT
// Q2 cut:
// =1: No cut.
// =2: Q2>2.10-2: well defined pVertex => meaningful D/dD cut.
// =3: Q2>.9: preselection of DIS.
// =4: Q2>1: DIS.
#define U3_Q2_CUT 2
// BMS related selection (rejection of events w/ bestpV w/ beam not OK)
// I) BMS-beamTelescope association: exists (BMS momentum not a fixed default
//   value), w/in range
// (Note that this is required, whether cpp defined or not, for phi selection.)
#define BMS_P_CUT	// The relevance (muon beam) and range is retrieved from "DataTakingDB"
// II) Backpropagation
// =1: No cut
// =2: PaParticle::BackPropLH or Chi2CutFlag as prescribed by "DataTakingDB"
#define BMS_LH_CUT 1	//
// Muons other than THE scattered one are most of the time pions, which have
// undergone a mu+nu decay, but which trajectories are reco'd as a single,
// continuous, track up- and down-stream of the decay.
// => One may then not want to DISCARD_MUONS
//#define DISCARD_MUONS
// pT cut for phi: pT is w.r.t. phi line of flight; cuts away e+e-
// It's set low: CSEvtTree analysis may set it tighter
#define Ephi_pTCUT 5  // pT ('hadron' w.r.t VM) cut for excl. phi uDST selection
#define Iphi_pTCUT 5  // pT CUT for incl. phi selection


// - DO NOT IMPACT uDST SELECTION
// (Note: RICH cuts do not impact selection, which is purely kinematics)
// - piS_e_REJECTION: Reject e polluting pi: by a severe cut on e/pi LH (which
//  hence affects Soft particles). This impacts on:
//   - K0 sample (actual one w/ mass cut about K0 peak, as opposed to pipi
//    sample w/in K0 histogramming range) for kinematical histogramming (Q2, y,
//    decay momenta), and for K0+X resonances.
//   - In Lambda case:
//     - "FillLambda" (cf. "./FillLambda.cc"), "U3_L_HIGHLEVEL" and histo of
//      kinematics of (p,pi) system,
//     - Lambda sample (actual one w/ mass cut about Lambda peak).
//     - Search for cascades.
#define piS_e_REJECTION
#define REJECT_SubKThr_e
#define Lambda_p_ID		// ID of p(\bar{p}) in [a]Lambda histogramming: e,pi,K veto below p-Threshold, pID above.
// RICH CENTRAL TUBE
// - According to our COMGeant setup, 5 cm is the design value (at least until
//  2012).
//   4 cm is some loose setting meant to explore the impact of explictitly
//  rejecting tracks impinging the RICH in its central tube (though it's, in
//  principle implicit in the RICH likelihood).
// - RICH_TUBE_RADIUS is also used to evaluate the RICH multiplicity (which is
//  output to CSEventData TTree).
#define RICH_TUBE_RADIUS 4

//#define RICH_ANGLE_CUT
//#define PS_USE_RADIUS

// REJECTION of UNRELIABLE TRACKS
// (Used in determining some auxiliary quantities:
// - # of detached hadrons, cf. "detachedHadrons".
// - RICH multiplicity, cf. "getRICHMultiplicity".
// - ...)
#define TRACK_NDF_MIN 4  // NdF = 4, means 8 or 9 hits (whether there's momentum or not)

//  **********************************************************************
//  ***** OVERWRITE some of the above for some SPECIFIC APPLICATIONS *****
//  **********************************************************************

#ifdef U3_uDST_PRESELECTION	// Default GOBAL CPP: all purpose K0Lambdaphi
// - CPP options influencing the uDST selection are set at their loosest.
// - But some CPP options, impacting only histogramming, are also set so.
// - No Tree output: the cuts, being too loose, we would  get too large a tree.
//  (This illustrates the diff. between tree and uDST. The latter serves as
//  a reservoir of all interesting events, that one would possibly want to
//  revisit after a first analysis is performed. (E.g. check the validity of
//  a RICH parameterization obtained on events w/ a pVertex in target region,
//  and hence expected to have been nicely reco'd, on more exotic events.) The
//  former is tailored so as to ease the analysis.)
// The following are set at their loosest by default. Yet, to be on the safe side...
#  undef  CHI2_CUT_OPT
#  define CHI2_CUT_OPT 1	// Chi2: Loosest cut
#  undef  U3_TARGET_CUT
#  define U3_TARGET_CUT 5       // Zone (w/ possibly ReZone), std to docu. TTree
#  undef  U3_y_CUT
#  define U3_y_CUT              // It's anyway always loose
#  undef  U3_Q2_CUT
#  define U3_Q2_CUT 2		// Well defined pVertex => meaningful D/dD cut
#  undef RICH_TUBE_RADIUS
// RICH_TUBE_RADIUS serves among other purposes to evaluate RICH multiplicity
// => Even though the incidence of each "hadron" track on RICH is output to tree
//   and can then be used to tailor RICH acceptance, one needs here a reasonable
//   value.
#  define RICH_TUBE_RADIUS 4	// Loosest cut
#  undef  DISCARD_MUONS
#  undef  Ephi_pTCUT
#  define Ephi_pTCUT 5		// Loose pT CUT for excl. phi uDST selection
#  undef  Iphi_pTCUT
#  define Iphi_pTCUT 5          // Loose pT CUT for incl. phi uDST selection
// The following are reset to their loosest...
// ...so as to have the widest possible selection
#  undef  BMS_P_CUT
#  undef  BMS_LH_CUT
// ...so as to get ab initio histos
#  undef  piS_e_REJECTION	// e-Rejection is still applied for RICH perfs.
#  undef  Lambda_p_ID
#  undef  U3_V0_CUT
// Allow one extra track for exclusivity selection
#  define U3_QUASIEXCL_PHI	// Enabled exceptionally: to explore the case
// Not Tree output
#  undef  U3_OUTPUT_TREE
#  define U3_OUTPUT_TREE	// Enable TTree (at least temporarily, for debugging)
// The following two sets of histos are only for dedicated mass shifts studies
#  undef  U3_K0_HIGHLEVEL
#  undef  U3_L_HIGHLEVEL
// Fill Lambda's: to have a look
#  undef  U3_FILL_Lambda
#  define U3_FILL_Lambda
// Only one GLOBAL CPP DEFINITION is allowed: check (at least partially) it's the case
#  ifdef U3_uDST_ANALYSIS
#    error: "Inconsistency! U3_uDST_ANALYSIS while U3_uDST_PRESELECTION also defined"
#  endif
#endif

#ifdef U3_uDST_ANALYSIS		// Re-analyzing preselected uDSTs w/ tighter cuts
// - CPP options influencing the event selection as well as those impacting
//  histogramming are tightened.
// - Adapted to e.g., RICH performances analysis: in particular the CSEventData
//  TTree is booked and filled. Athough enabling the CSEventData TTree at the
//  "U3_uDST_PRESELECTION" seems feasible: tree size keeping being manageable.
// - May not be adapted to physics analyses (K0 multiplictiy, Lambda physics)
#  undef CHI2_CUT_OPT
#  define CHI2_CUT_OPT 2	// Tighter cut
#  undef U3_TARGET_CUT
#  define U3_TARGET_CUT 1	// Standard cut
#  undef  U3_y_CUT
#  define U3_y_CUT              // It's anyway always loose
#  undef U3_Q2_CUT
#  define U3_Q2_CUT 2		// Well defined pVertex => meaningful D/dD cut
#  undef RICH_TUBE_RADIUS
// RICH_TUBE_RADIUS serves among other purposes to evaluate RICH multiplicity
// => Even though the incidence of each "hadron" track on RICH is output to tree
//   and can then be used to tailor RICH acceptance, one needs here a reasonable
//   value.
#  define RICH_TUBE_RADIUS 5	// Design value
// Let's not enable DISCARD_MUONS: it's thought to be:
// I) Inefficient, because few muons are produced, aprt from the scattered one.
// II) Harmful: it in fact erases pions decaying to muons and froming w/ these
//    a one and single track.  
#  undef  DISCARD_MUONS
#  undef  Ephi_pTCUT
#  define Ephi_pTCUT 25		// Tighter pT CUT for excl. phi uDST selection
#  undef  Iphi_pTCUT
#  define Iphi_pTCUT 25         // Tighter pT CUT for incl. phi uDST selection
// - The following two are not expected to affect K0 and Lambda reconstructions.
// The 1st is applied to phi reconstruction in any case. The 2nd one affects
//  directly exclusivity (via muon kine) but only marginally, the phi selection.
// - Yet, they don't affect much statistics either => Enable them
//#  undef BMS_P_CUT		// 
//#  undef BMS_LH_CUT		// 
#  undef  piS_e_REJECTION	// e-Rejection would still be applied for RICH perfs.
#  define piS_e_REJECTION	// Its anyway a well justified cut
#  undef  Lambda_p_ID
// The following is not in "U3_uDST_PRESELECTION"
#  define U3_V0_CUT 1		// Slightly tighter cuts than U3_uDST_PRESELECTION
// Extra track for exclusivity selection: preliminary evaluation were negative
#  undef  U3_QUASIEXCL_PHI	// To be enabled exceptionally: only to explore the case
#  define U3_OUTPUT_TREE	// Do enable TTree output this time
#  undef U3_ALLOUT_KINE		// These can in any case be filled from TTree
// The following two sets of histos are only for dedicated mass shifts studies.
#  undef  U3_K0_HIGHLEVEL
#  undef  U3_L_HIGHLEVEL
// No Lambda's here: reserved for a dedicated analysis
#  undef  U3_FILL_Lambda
// Only one GLOBAL CPP DEFINITION is allowed: check (at least partially) it's the case
#  ifdef U3_uDST_PRESELECTION
#    error: "Inconsistency! U3_uDST_PRESELECTION while U3_uDST_ANALYSIS also defined"
#  endif
#endif

#ifdef U3_Lambda_uDST_PRESEL	// Preselection aimed at physics. Not sure it's still useful...
// - Tighter selection than U3_uDST_PRESELECTION, else what's the use?
#  undef  CHI2_CUT_OPT
#  define CHI2_CUT_OPT 2	// Tighter cut
#  undef  U3_TARGET_CUT
#  define U3_TARGET_CUT 1       // Standard (Lambda on FI in 2016/17 makes no sense)
#  undef  U3_y_CUT
#  define U3_y_CUT              // It's anyway always loose
#  undef  U3_Q2_CUT
#  define U3_Q2_CUT 2		// Well defined pVertex => meaningful D/dD cut
#  undef RICH_TUBE_RADIUS
#  define RICH_TUBE_RADIUS 5	// Design value
#  undef  DISCARD_MUONS
// The following are out of Lambda scope
#  undef  Ephi_pTCUT
#  undef  Iphi_pTCUT
// Explore the following makes no sense for Lambda
//#  undef  BMS_P_CUT
//#  undef  BMS_LH_CUT
#  undef  piS_e_REJECTION
#  define piS_e_REJECTION	// Well justified. Not sure it's used for Lambda
#  undef  U3_FILL_Lambda
#  define U3_FILL_Lambda
//
#  undef  Lambda_p_ID		// Do undefine, in order to get ab initio histos
#  ifndef U3_FILL_Lambda
#    #error: "Inconsistency! U3_Lambda_uDST_PRESEL while U3_FILL_Lambda not defined"
#  endif
// Only one GLOBAL CPP DEFINITION is allowed: check (at least partially) it's the case
#  ifdef U3_uDST_PRESELECTION
#    error: "Inconsistency! U3_Lambda_uDST_PRESEL while U3_uDST_PRESELECTION also defined"
#  endif
#endif

#ifdef U3_Lambda_ANALYSIS	// For re-analyzing "Lambda" uDSTs
// - Standard target cuts (PaAlgo::InTarget)
#  define CHI2_CUT_OPT 2	// Tighter cut, as prescribed by "DataTakingDB"
#  define U3_Q2_CUT 2		// So as to have well defined pVertex
#  define U3_TARGET_CUT 1	// 1 = Standard target cuts = PaAlgo::Intarget
#  if U3_SPECIAL_SCHEME == 2008
#    undef piS_e_REJECTION
#  endif
#  define Lambda_p_ID		// ID of p(\bar{p} in Lambda(\bar{Lambda}) selection: e,pi,K veto below p-Threshold, pID above.
#  define U3_V0_CUT 1		// Slightly tighter cuts than U3_uDST_PRESELECTION
#  ifndef U3_FILL_Lambda
#    #error: "Inconsistency! U3_Lambda_uDST_PRESEL while U3_FILL_Lambda not defined"
#  endif
// The following two sets of histos are only for dedicated mass shifts studies.
#  undef  U3_K0_HIGHLEVEL
#  undef  U3_L_HIGHLEVEL
// Only one GLOBAL CPP DEFINITION is allowed: check (at least partially) it's the case
#  if defined U3_uDST_PRESELECTION || defined U3_L_MC_ANALYSIS
#    error: "Inconsistency! U3_uDST_PRESELECTION || U3_L_MC_ANALYSIS while U3_Lambda_ANALYSIS also defined"
#  endif
#endif

#ifdef U3_L_MC_ANALYSIS 	//
#  define U3_TARGET_CUT 1	// Retain target and (if !hadron setup) y cuts
#  define RICH_TUBE_RADIUS 5
#  undef  REJECT_SubKThr_e	// No PID in MC, yet
#  define Lambda_p_ID		// ID of p(\bar{p} in Lambda(\bar{Lambda}) selection: e,pi,K veto below p-Threshold, pID above.
#  define U3_IDEAL_RICH		// !!!!!! BEWARE: SPECIAL IDEAL RICH SETTING !!!
#  define U3_V0_CUT 1		// Take advantage of benefits brought by RICH to loosen VV dist and pT cut
#  ifndef U3_FILL_Lambda
#    #error: "Inconsistency! U3_Lambda_uDST_PRESEL while U3_FILL_Lambda not defined"
#  endif
// The following two sets of histos are not relevant to MC
#  undef  U3_K0_HIGHLEVEL
#  undef  U3_L_HIGHLEVEL
// Only one GLOBAL CPP DEFINITION is allowed: check (at least partially) it's the case
#  if defined U3_uDST_PRESELECTION || defined U3_Lambda_ANALYSIS
#    error: "Inconsistency! U3_uDST_PRESELECTION || U3_Lambda_ANALYSIS while U3_L_MC_ANALYSIS also defined"
#  endif
#endif

#ifdef U3_ERASE_K0s	        // Single out K0 and erase their pion decays: to be used in conjunction w/ another UserEvent....
// ...=> Care to have tight enough V0 cuts (else K0 purity becomes bad)
#  undef  U3_V0_CUT
#  define U3_V0_CUT 1
// ...=> Cancel "uDSTSelection"
#  undef K0LambdaphiXi_WRITEBACK
// ...=> Do undefine the booking of histos w/ either too high a level or too
//      general a purpose.
#  undef  U3_FILL_Lambda
#  undef  U3_FILL_RICH_PERFS
#  undef  U3_GENERALPURPOSE_HISTOS
#  undef  U3_K0_HIGHLEVEL
#  undef  U3_L_HIGHLEVEL
#  define U3_Q2_CUT 1		// Q2: no cut
#  ifdef U3_uDST_PRESELECTION
#    error: "Inconsistency! U3_uDST_PRESELECTION while U3_ERASE_K0s also defined"
#  endif
#endif


//    ***** CHECK/ENSURE CONSISTENCY among CPP DEFINITIONS *****
//    ***** #INCLUDE WHAT NEED BE                          *****
#ifndef U3_FILL_Lambda
#  ifdef U3_FILL_Lambda
#    warning U3_L_HIGHLEVEL defined whereas U3_FILL_Lambda not defined
#  endif
#  undef U3_L_HIGHLEVEL
#endif
// ***** BEST pV *****
#ifndef BESTpV_OPTION
#  warning BESTpV_OPTION was undefined: set =1
#  define BESTpV_OPTION 1	// Option "Best pV from CORAL if useful" is default
#endif
// ***** OUTPUT TREE *****
#ifdef U3_OUTPUT_TREE
#  include "TTree.h"
#  include "CSEventData.h"
#endif

//  **********************************************************************
//  *****            INCLUDE "UserEvent103.h" ONLY NOW...
//      ...so that histos are booked according to cpp definitions
//  **********************************************************************
#include "UserEvent103.h"  // Histograms

// *************************************************************************
// ************************ END of BLOCK of CPP VERSIONS *******************
// *************************************************************************


// ********** GLOBALS: PARAMETERS and CUTS **********
// ***** GLOBALS which setting depends upon CPP DEFINITIONS
static int Run = 0, EvNum;
static int Year;
static double chi2Cut = 10; // Will be updated in any case: see CHI2_CUT_OPT
static bool   hadronPhysics = false, muonBeam = false, getBestPVfromCORAL = true;
static double Q2Cut = 0;
static double yLowCut = 0, yUpCut = 2;
static double pBMS, dpBMS = 0;
static double dqP2BMSCut = 0; // Cut on cov(1/P,1/P) to reject beams w/o BMS
static unsigned short spillOKPat = 0;
#ifdef U3_BAD_SPILLS
static set<unsigned> badSpills, badSpillsCalo, badSpillsRICH;
#endif

// ***** Other GLOBALS
static double D0piLow  = -.0028, D0piUp  = .0032;
static double zCutS = .20, cLowCutS = -.85, cUpCutS = .85;
static double zCut0 = .25, cLowCut0 = -.50, cUpCut0 = .50;
// RICH GLOBALS
static RICH *rich;
static double piThr, KThr, pThr;
// FI/MM ASCISSAE
static double ZFI02X, ZMM02X;
#define U3_FIMP_VTCS
#ifdef U3_FIMP_VTCS
// Z abscissae of FI03 and MP01UV planes, to histogram the vertices
// reco'd thanks to the ReZone option of CORAL.
static double ZFI03s[3], ZMP01UVs[2];
static double corrFI03 = 55.55-(52+52.5+53.6)/3;
static double corrMP01UV = 141.4-(140.932+142.6055)/2;
static int mpWs[3], dcWs[3]; static unsigned int mpBs[3], dcBs[3];
#endif

// **************************************************************************
// ****************************** RICH ID CUTS ******************************
// **************************************************************************
//           ********** CUT ON MOMENTUM **********
#if defined U3_FILL_RICH_PERFS || defined U3_FILL_Lambda
static double PRICHCut = 45; // 45 is typicall of 2002-4. May be updated infra.
#endif

#if defined U3_FILL_Lambda || defined Lambda_p_ID || defined U3_FILL_RICH_PERFS
// pID
static double pIDPmax = 80;
static double pIDPmin     = 1.05 /* times pThr*/;
static double subpThrPmin = 1.1  /* times piThr */;
#endif
// pID
// (Note: over much of the P range, there is no overlap between pi (or K) and
// p => No need to have such a tight veto as in the KID case. But this no longer
// holds true at higher P, where the cut setting infra may become too loose...)
static double subpThrLHpiVeto = 1.40, subpThrLHVeto = 1.50;

#ifdef U3_USE_CHCHI2
// Cherenkov chi2 (used instead of LH)
static double CHChi2Cut = 0.95, CHChi2Mx = 5, subKThrChi2piVeto = 4;
#endif

//             ********** CUT on LH **********
static double LHCut = 1.01, LHBckCut = 1.20;
static double subKThrLHpiVeto = 0.92; // Also veto 2nd largest, even if not pi
static double subKThrPmin = 1.5;
static double thrMargin = 1.05; // Safety margin: require P > Thr*margin

// ***** ELECTRON REJECTION
//  Note: 3 different stories:
// i)   eVeto for pi selection "piS_e_REJECTION"
// ii)  eVeto for pi selection below piTHr
// iii) SubThr K selection.
//  In (i and (ii), we arbitrate the efficiency/purity tradeoff in favour of
// efficiency (of pion selection).
//  In (iii), we are below threshold and arbitrate in favour of the purity (of
// e.g. the K selection).
static double piLHeVeto = 1.5;
static double subpiThrLHeVeto = 1.5;
#ifdef REJECT_SubKThr_e
static double subKThrLHeVeto = 0.92;
#endif

//           ********** CUT ON dTheta **********
// (Note: This is obsolete. Only kept just in case...)
static double nSigmas = 3, nVeto = 2;  // # of sigmas used for pID and Veto cuts

//      ********** RICH ACCEPTANCE **********
#ifdef RICH_TUBE_RADIUS
const double rRICH = RICH_TUBE_RADIUS;
#else
const double rRICH = 0;
#endif
  // 2008/8/6 Giulia Pesaro <gpesaro@cern.ch>:
  // z_rich= 615.6
  //  TMath::Abs(X)<161 && TMath::Abs(Y)<121
  // ((pow(X,2)+pow(Y,2))>25)
const double ZRICH = 615.6, r2RICH = rRICH*rRICH, XRICH = 161, YRICH = 121;
static int thCIndex;
// Following two affect only RICHPerfs?
static double rRICHCut = 12;  // ***** Radius cut *****
#ifdef RICH_ANGLE_CUT
static double aRICHCut = .025;  // ***** Angle cut *****
#endif

//      ********** SUBTHRESHOLD **********
// (Note: In earlier RICH code, events w/o any Cherenkov photons are deprived
// of RICH block. Let's consider these just as good a KID as standard piVeto.)
#define K_BELOW_THR 2 	// >1: No RICH block, w/in RICH & range => KID

//   ***** ARRAYS of RICH INFO for ALL TRACKS *****
static int *PIDs;
static double *richLHs[5], *richChi2s[5];
static double *richDths; // dtheta/dn/d1/beta * DeltaN
#ifdef U3_FILL_RICH_PERFS
#  ifdef PS_USE_RADIUS
static double RRBin = 6;
#  else
static double ARBin = .012;
#  endif
#endif
static int *tZones, *aRs;

static double ZSM1, ZSM2;// Z absissae of magnets (retrieved from PaSetup)...
static double ZTarget;   // ...and target

#ifdef U3_OUTPUT_TREE
static CSEventData *fCSEvt;
static vector<CSHadronData> fHadrons;
static map<int,int> fPa2fH;  // X-reference PaParticle::MyIndex -> fHadrons
static vector<CSResonanceData> fResonances;
static TTree *fCSEvtTree;
static PaPid *fPid;
#endif

// ********** INTERFACES **********
#ifdef SCALE_MOMENTA
void RescaleMom(PaEvent & e, bool faster_mode);
#endif
void GetPIDs(PaEvent &e, RICH *rich,
	     double runIndex, double prodIndex, double deltaIndex, int thCIndex,
	     int *PIDs, double **richLHs, double **richChi2s,
	     int *tZones,
	     double *richDths);
void parseBestPV(PaEvent &e, int ipV,
		 unsigned short &primaryPat,
		 double &Xp, double &Yp, double &Zp, TMatrixD &CovP,
		 TLorentzVector &lvp, double &E0,
		 int &imuS, const PaParticle *&muS,
		 TLorentzVector &lvq,
		 double &yB, double &Q2, double &xB);
void bookKinematics(double phiMassCut, double K0MassCut,
		    const double *dECuts, int Ephi_pTCut,
		    const double *V0DsdD, const double *V0cthCut, const double *V0pTCut);
int GetTrigType(unsigned int evTrig);
#ifdef U3_ALLOUT_KINE
void fillAllOutKine(const PaEvent &e, int trigType,
		    int ipV, const PaVertex &pV, unsigned short primaryPat,
		    int imup, const PaParticle *mup,
		    double Q2, double xB, double yB, const TLorentzVector &lvq);
#endif
void fillExclusive(PaEvent &e, int ipV, int imuS,
		   const TLorentzVector &lvp, const TLorentzVector &lvq,
		   const double *dECuts, int pTCut,
		   double dM_rho, double dM_phi, double phiMassCut,
		   int &KEMiss, bool &uDSTSelection, bool &isphi);
void bookExclusive(const double *dECuts, int pTCut, double dM_rho,
		   double dM_phi, double phiMassCut);
void bookEphi_RICHPerf(double dM_phi, double phiMassCut);
void fillEphi_RICHPerf(PaEvent &e, int iET1, int iET2,
		       const TLorentzVector &lvK1, const TLorentzVector &lvK2,
		       double m_KK, bool isphi);
void fillIphi(PaEvent &e, int ipV, int imuS, double dM_phi,
	      bool &uDSTSelection);
void bookV0Selection();
void bookK0(const double *V0DsdD, const double *V0cthCut, const double *V0pTCut,
	    double dM_K0, double K0MassCut);
void fillK0woPV(double m_pipi, double Xs, double Ys, double Zs,
		double pT, double zL1, double zL2, const TLorentzVector &lvpipi);
#ifdef U3_K0_HIGHLEVEL
void fillK0HighLevel(const PaTrack &trk1, const PaTrack &trk2,
		     double m_pipi, const TLorentzVector &lvpipi,
		     double Theta, double Thetx, double Thety, double ThetX,
		     double alpha, double beta);
void getK0Incidence(const PaTrack &trk1, const PaTrack &trk2,
		    int zones1, int zones2,
		    const TVector3 &v13, const TVector3 &v23,
		    const PaTPar &hi1, const PaTPar &hi2,
		    double &Theta, double &Thetx, double &Thety, double &ThetX,
		    double &alpha, double &beta);
void fillK0IncidenceLAS(double m_pipi, const TLorentzVector &lvpipi,
			double Theta, double Thetx, double Thety, double ThetX,
			double alpha, double beta);
void fillK0IncidenceSAS(double m_pipi, const TLorentzVector &lvpipi,
			const PaTPar &hv1, const PaTPar &hv2,
			double Theta, double Thetx, double Thety, double ThetX,
			double alpha, double beta);
#endif
void bookK0X(double dM_phi, double K0MassCut, double KSMassCut);
void fillK0X(PaEvent &e, int ipV, int imuS,
	     const TLorentzVector &lvq,
	     int iET1, int iET2, const TLorentzVector &lvK0,
	     double KSMassCut, double dM_phi, double phiMassCut);
void bookK0_RICHPerf(double dM_K0, double K0MassCut);
void bookK0_RICHPur(double dM_K0);
void fillK0_RICHPerf(PaEvent &e, int iET1, int iET2,
		     const TLorentzVector &lvpi1, const TLorentzVector &lvpi2,
		     double m_pipi, bool isK0, bool cleanV0, int Run);
#ifdef U3_FILL_Lambda
void BookLambda(const double *V0DsdD, const double *V0cthCut, double dM_Lambda,
		double pTCut, int PIDScheme, int histoLevel = 0);
void FillLambda(double m_ppi, double m_pip,
		int s_dist, int s_angle,  // Vertex line cut indices
		bool hasmumuS,            // p-V has mu/mu'
		const PaTrack &trk1,      // >0 track
		const PaTrack &trk2,      // <0 track
		bool downstreamOfTarget,  // s-Vertex downstream of target
		int pID);                 // 2-bit pattern, w/ LSB = >0 particle
void bookLambdaMisc(const double *V0DsdD, const double *V0cthCut, const double *V0pTCut,
		    double dM_Lambda, double LMassCut, double *xFbins,
		    int PIDScheme);
void printLambdaCuts(double DsdDCut, double cthCut, double pTCut, double vChi2Cut);
#  ifdef U3_L_HIGHLEVEL
void FillLambda(double m_ppi, int histoLevel);
#  endif
void fillLambdaX(PaEvent &e, int ipV, int imuS,
		 const TLorentzVector &lvq,
		 int iET1, int iET2, int iLaL, const TLorentzVector *lvL);
void fillLCascades(PaEvent &e, int ipV, TMatrixD &CovP, int imuS,
		   int isV, TMatrixD &CovS, int iET1, int iET2, int iLaL,
		   const TLorentzVector *lvL, unsigned short ppiID,
		   bool &uDSTSelection);
//void FillLambdaMC(double m, int pID);
void bookLambda_RICHPerf(double dM_Lambda, double LMassCut);
void fillLambda_RICHPerf(PaEvent &e, int iET1, int iET2, int iETp,
			 const TLorentzVector &lvp1,
			 double m_ppi, bool &isLambda, int Run);
#endif
void transport(PaEvent &e,
	       PaTrack &trk,    // 
	       int &aR,         // RICH acceptance flag
	       int zones);
#if defined U3_K0_HIGHLEVEL || defined U3_L_HIGHLEVEL
void getGM0709WB(int *gm0709Ws, unsigned int *gm0709Bs);
#  ifdef U3_K0_HIGHLEVEL
void fillK0GM0709(double m_pipi, const PaTrack &trk1, const PaTrack &trk2,
		  const int *gm0709Ws, const unsigned int *gm0709Bs);
#  endif
#  ifdef U3_L_HIGHLEVEL
void fillLGM0709(double m_ppi, const PaTrack &trk1, const PaTrack &trk2,
		 const int *gm0709Ws, const unsigned int *gm0709Bs);
#  endif
#endif
#ifdef U3_OUTPUT_TREE
void copyCSEventHeader(PaEvent &e,
		       unsigned short primaryPat, unsigned short spillOKPat,
		       const PaVertex &pV,
		       double E0, double Q2, double xB, double yB,
		       double index, double prodIndex);
void copyCSHadronData(PaEvent &e, const PaParticle *pa, int iV);
void copyCSResonanceData(const PaVertex &v,
			 double m, double D, double dD, double cth,
			 double pT, double alpha);
int getHadronIndex(const PaParticle *pa);
int nChargedKsInPV(PaEvent &e, const PaVertex &pV);
#endif
void getFIMMAbscissae();
void getTargetBinning(int &nZbins, double &ZMn, double &ZMx);

//#define DEBUG_K0
#ifdef DEBUG_K0
static FILE *dbK0_o = 0;
// DEBUG FLAG:
// - ==-3: No input
// - ==-2: Init
// - ==-1: No K0 expected in current event
// - >=0: Bit pattern housing the amount of achievement:
//   0x1, 0x2: >0 and <0 track, 0x3
static int dbK0Event = -2;
static int dbK0p1p2 = 0; // 0x3 if current 2nd vertex contain expected p1 p2
static float dbK0p[2], dbK0a[2]; // Characteristics of K0 decays
static float dbK0v1, dbK0v2; // Value outside cut
void dbK0GetEvent(const PaEvent &e);
#endif

void setPeriodDep(PaEvent &e)
{
  // *************** PERIOD DEPENDENT PARAMETERS and SETTINGS ***************
  // Settings are retrieved from DataTakingDB
#ifdef U3_Lambda_LOOSE
  int RICHOption = 2;
#else
  int RICHOption = 1;
#endif
#ifdef U3_TARGET_CUT
#  if U3_TARGET_CUT == 5
  int targetOption = 1;
#  else
  int targetOption = U3_TARGET_CUT;
#  endif
#else
  int targetOption = 0;
#endif
  DataTakingDB *DTDB = DataTakingDB::Ptr();
  if (!DTDB) DTDB = new DataTakingDB();
  DTDB->Init(e,BESTpV_OPTION,CHI2_CUT_OPT,RICHOption,targetOption);
  Year = DTDB->Year;
  chi2Cut = DTDB->Chi2Cut;
  // "hadronPhysics": can be muon beam: Primakoff
  hadronPhysics =
    DTDB->PhysicsType->find("Hadron")==0 ||
    DTDB->PhysicsType->find("Primakoff_mu")==0 ||
    DTDB->PhysicsType->find("DY")==0;
  // Muon beam
  muonBeam = !hadronPhysics || DTDB->PhysicsType->find("Primakoff_mu")==0;
  pBMS = DTDB->PBeam; dpBMS = DTDB->DPBMS; dqP2BMSCut = DTDB->DqP2BMSCut;
  yLowCut = DTDB->yLowCut; yUpCut = DTDB->yUpCut;
  getBestPVfromCORAL = DTDB->GetBestPVfromCORAL;    // Best p(rimary)Vertex
#ifdef U3_BAD_SPILLS
  DTDB->GetBadSpills(badSpills,badSpillsCalo,badSpillsRICH);
#endif
  // RICH settings
  LHCut = DTDB->LHCut; LHBckCut = DTDB->LHBckCut;
#if defined U3_FILL_RICH_PERFS || defined U3_FILL_Lambda
  PRICHCut = DTDB->PRICHCut;
#endif
  pIDPmax = DTDB->pIDPmax; pIDPmin = DTDB->pIDPmin; subpThrPmin = DTDB->SubpThrPmin;

  return;
}

// ************************************************************
// ************************* UserEvent ************************
// ************************************************************

void UserEvent103(PaEvent &e)
{
  int i, j, k; // Generic loop indices

  // **********************************************************************
  //           ******************** CUTS ********************
  // **********************************************************************

  // ********** INELASTICITY **********
  const double dECuts[2] = {2.5,2.0}; // w/o, w/ ELoss correction; the latter not used, in fact

  //          ********** V0 SELECTION **********
  // VV dist[0]: Event selection, [1]: Histogramming, [2]: Clean signal, Cascades
  // Angle  [0]: Event selection, [1]: Histogramming, [2]: ??
  // pT Cut [0] intended for K0, [1] for Lambda
#ifndef U3_V0_CUT  // To be used for uDST selection
  const double V0DsdD[3] = {1,2,5}, V0cthCut[3] = {.9998, .9999, .99995 };
  const double V0pTCut[2] = {.010,.020 }; 
#else
  // Slightly tighter cuts, e.g. for TTree output
#  if   U3_V0_CUT == 1
  const double V0DsdD[3] = {2,3,5}, V0cthCut[3] = {.9998, .9999, .99995 };
  const double V0pTCut[2] = {.020,.025 }; 
#  endif
#endif

  //           ***** MASS SELECTIONS *****
  // ***** Cuts for uDST selection = Crude mass cuts
  const double dM_phi = 0.100, dM_K0 = 0.075; // dM_Lambda: cf. infra
  const double dM_rho = 0.450;
  // ***** Cuts for signal selection = Tight mass cuts 
  const double K0MassSigma = .0076;  // Measured om hm_K0c22(DsdD=5), no scaling
  const double K0MassCut = 2*K0MassSigma;
#ifdef U3_FILL_Lambda // Special mass (and pT) cuts for Lambda
  const double dM_Lambda = 0.050;
  const double LMassSigma = .0026;   // Measured om hm_Lc22(DsdD=5), no scaling
  const double LMassCut = 2*LMassSigma;
#endif
  const double phiMassSigma = .0036; // Measured om hm_phi(1:1)(dE<2.5), no scaling
  const double phiMassCut = 2*phiMassSigma;
  const double KSMassCut = .050;


  // **********************************************************************
  // ********************    RUN (=> period), SPILL    ********************
  // **********************************************************************

  static int Spill = -1;
  int prvRun = Run; Run = e.RunNum();
  int prvSpill = Spill; Spill = e.SpillNum();
  EvNum = (int)e.UniqueEvNum();

  //                          ******************** NEW RUN ********************
  if (Run!=prvRun) // Mainly: Call DataTakingDB and retrieve its settings
    setPeriodDep(e);

#ifdef U3_BAD_SPILLS
  if (Spill!=prvSpill) {  // ******************** NEW SPILL ********************
    unsigned rs = (Run<<12)+Spill;
    // There is room for THREE bad spill lists. Yet, most of the time only one
    // is available. Then, by default, the corresponding OK bits are set.
    set<unsigned> *bs[] = {&badSpills, &badSpillsCalo, &badSpillsRICH};
    int ibs; for (ibs = 0, spillOKPat = 0x7; ibs<3; ibs++) { 
      set<unsigned> *b = bs[ibs];
      if (b->find(rs)!=b->end()) spillOKPat ^= 1<<ibs;
    }
  }
#endif

  // ***** MOMENTUM RESCALING ALTERNATIVES *****
#ifdef SCALE_MOMENTA
  if (Year==2006) RescaleMom(e,true);
#endif

  // **********************************************************************
  //        ******************** INITIALIZATIONS ********************
  // **********************************************************************
#ifdef U3_FILL_Lambda
  static double xFbins[5] = {.05,.17,.24,.34,.50};
  static MCLambda *LambdaMC = 0;
#endif
  static bool first(true); if (first) {
    first=false;
    
    // ***** INSTANTIATE RICH **********
    if ((rich = RICH::Ptr())==0) rich = new RICH(nSigmas,nVeto);

    // *************** INSTANTIATE "MCInfo" "MCLambda" ***************
#ifdef U3_FILL_Lambda
    if (e.IsMC()) {
      if (MCInfo::Ptr()==0) new MCInfo(Year);
      double RTCut = DataTakingDB::Ptr()->RTCut, R2TCut = RTCut*RTCut;
      MCInfo::Ptr()->BookLambdaAcceptance(yLowCut,yUpCut,R2TCut,V0pTCut[1]);
      LambdaMC = new MCLambda(yLowCut,yUpCut,V0pTCut[1]);
    }
#endif

    //   ************************* TREE  *************************
#ifdef U3_OUTPUT_TREE
    gDirectory->cd("/");
    fCSEvt = new CSEventData;
    //fHadronsPtr = &fHadrons; fResonancesPtr = &fResonances;
    fCSEvtTree = new TTree("CSEvtTree","CS event");
    fCSEvtTree->Branch("CSEvt","CSEventData",&fCSEvt);
    //fCSEvtTree->Branch("Hadrons","std::vector<CSHadronData>",&fHadronsPtr);
    fCSEvtTree->Branch("Hs",&fHadrons);
    //fCSEvtTree->Branch("Resonances","std::vector<CSResonanceData>",&fResonancesPtr);
    fCSEvtTree->Branch("Rs",&fResonances);
    fCSEvtTree->SetMaxTreeSize(1000000000);
    fPid = new PaPid;
#endif

    // **********************************************************************
    // ************************* HISTOGRAMS BOOKING *************************
    // **********************************************************************

    //            *************** GENERAL PURPOSE ***************
#ifdef U3_GENERALPURPOSE_HISTOS
    h_triggers = new TH1D("h_triggers","Trigger Mask", 512,-.5,511.5);
    h_spills   = new TH1D("h_spills","Spills",250,-.5,249.5);
    h_primary  = new TH1D("h_primary","pV parsing",1024,-0.5,1023.5);
    h_pBMS     = new TH1D("h_pBMS","|P| BMS (GeV)",
			  12*int(dpBMS),pBMS-1.5*dpBMS,pBMS+1.5*dpBMS);
    h_dqPBMS   = new TH1D("h_dqPBMS","#deltaq/P BMS",100,0,2.5e-4);
#endif
    //             *************** GENERAL RICH ***************
#ifdef U3_FILL_RICH_PERFS
    int strtRun = DataTakingDB::Ptr()->StrtRun;
    int stopRun = DataTakingDB::Ptr()->StopRun;
    int nRuns = stopRun-strtRun+1;
    for (i = 0; i<3; i++) {         // ********** dtheta vs. RUN # **********
      // For K+/-, pi+ and pi-. This is in view of adjusting index(run#).
      // (Since it needs "nRuns,strtRun,stopRun" => not in RICH specific
      // booking subroutines.)
      char hN[] = "hR_dvsr0";
      char hT[] = "RICH #pi#pm#delta#theta  -  K0#pm15Mev  ";
      sprintf(hN,"hR_dvsr%d",i);
      switch (i) {
      case 0: sprintf(hT,"RICH #pi+#delta#theta  -  K0#pm%.0fMev",
		      K0MassCut*1000);
      hR_dvsr[0] = new TH2D(hN,hT,nRuns,strtRun,stopRun+1,100,-10,10); break;
      case 1: sprintf(hT,"RICH #pi-#delta#theta  -  K0#pm%.0fMev",
		      K0MassCut*1000);
      hR_dvsr[1] = new TH2D(hN,hT,nRuns,strtRun,stopRun+1,100,-10,10); break;
      case 2: sprintf(hT,"RICH K#pm#delta#theta  -  #phi#pm%.0fMev",
		      phiMassCut*1000);
      hR_dvsr[2] = new TH2D(hN,hT,nRuns,strtRun,stopRun+1,100,-10,10); break;
      }
    }
#endif

    //                   ********** KINEMATICS **********
    bookKinematics(phiMassCut,K0MassCut,
		   dECuts,Ephi_pTCUT,V0DsdD,V0cthCut,V0pTCut);

    //             *************** EXCLUSIVE phi ***************
    bookExclusive(dECuts,Ephi_pTCUT,dM_rho,dM_phi,phiMassCut);
#ifdef U3_FILL_RICH_PERFS
    // ***** RICH PERFORMANCES from EXCLUSIVE phi
    bookEphi_RICHPerf(dM_phi,phiMassCut);
#endif

    // ***** V0 SELECTION HISTOS: VV DISTANCE, COLLINEARITY, V-CHI2, pT *****
    bookV0Selection();

    //                 *************** K0 ***************
    bookK0(V0DsdD,V0cthCut,V0pTCut,dM_K0,K0MassCut);
    // ***** K0X  resonances
    bookK0X(dM_phi,K0MassCut,KSMassCut);
#ifdef U3_FILL_RICH_PERFS
    // ***** RICH PERFORMANCES from K0
    bookK0_RICHPerf(dM_K0,K0MassCut); bookK0_RICHPur(dM_K0);
#endif

    //               *************** LAMBDA ***************
#ifdef U3_FILL_Lambda
    int PIDScheme = 0;
#  ifdef Lambda_p_ID
    PIDScheme |= 1;
#  endif
#  ifdef piS_e_REJECTION
    PIDScheme |= 2;
#  endif
#  ifdef U3_L_HIGHLEVEL
    int histoLambdaLevel = 0x1f;
#  else
    int histoLambdaLevel = 0x1f;
#  endif
    if (e.IsMC()) histoLambdaLevel |= 0x100;
    BookLambda(V0DsdD,V0cthCut,dM_Lambda,V0pTCut[1],PIDScheme,histoLambdaLevel);
#  ifdef U3_FILL_RICH_PERFS
    bookLambda_RICHPerf(dM_Lambda,LMassCut);
#  endif
    bookLambdaMisc(V0DsdD,V0cthCut,V0pTCut,dM_Lambda,LMassCut,xFbins,
		   PIDScheme);
#endif

    //       ********** END OF HISTOGRAM BOOKING **********

    // **********************************************************************
    //   ********** DETERMINE SOME GLOBAL VARIABLES (from cpp, PaSetup) *****
    // **********************************************************************
    // Q2Cut
#if   U3_Q2_CUT == 1
    Q2Cut = 0;
#elif U3_Q2_CUT == 2
    Q2Cut = 2e-2;
#elif U3_Q2_CUT == 3
    Q2Cut = .9;
#else
    Q2Cut = 1.;
#endif
    
    const PaSetup &setup = PaSetup::Ref();

    ZTarget = setup.TargetCenterZ();

    // ***** ABSCISSAE of MAGNETS
    PaMagInfo *m = setup.PtrMagField()->getMagInfo();
    int nms = setup.PtrMagField()->getNumOfMags(), iSM2 = nms-1, iSM1 = iSM2-1;
    ZSM1 = m[iSM1].zcm/10; ZSM2 = m[iSM2].zcm/10;

  }

  // **********************************************************************
  // ******************** END OF INITALISATION ********************
  // **********************************************************************

  //            ********** RETRIEVE RICH INDICES(run) **********
  static double runIndex, prodIndex, runIndVS, prodIndVS;
  if (Run!=prvRun) {
    rich->UpdateIndices(Run,e.IsMC(),runIndVS,prodIndVS,runIndex,prodIndex);
    //#define U3_CORR_INDEX
#ifdef U3_CORR_INDEX
    double corr = +40e-6; runIndex *= 1+corr; runIndVS *= 1+corr;
    rich->SetNRICH(runIndVS,runIndex); rich->UpdateThresholds();
#endif
    printf(" * U3: Indices for RICH PID (run/prod[/mean]): %.6f/%.6f/%.6f %.6f/%.6f\n",
	   runIndex,prodIndex,DataTakingDB::Ptr()->MeanIndex,
	   runIndVS,prodIndVS);

  }
  //                                                      ***** UPDATE THRSHOLDS
  piThr = rich->piThr; KThr = rich->KThr; pThr = rich->pThr;


  //       ******************** DEBUGGING ********************
  static int iEvt = 0;
  //#define DEBUG_U3 68600
#ifdef DEBUG_U3 
#  if DEBUG_U3 > 10000
  if (iEvt<DEBUG_U3) { // If it's a large #, let's delay the printing out 
    iEvt++; return;
  }
  printf("%d  Run %d Evt %d\n",iEvt,Run,EvNum);
#  else
  if (!(iEvt%DEBUG_U3))
    printf("Evt %d  Run %d Evt %d\n",iEvt,Run,EvNum);
#  endif
#endif
  iEvt++;


  //   ********** MCINFO PROCESSES MC DATA TO SORT OUT USEFUL INFO **********
#ifdef U3_FILL_Lambda
  if (e.IsMC()) {
    MCInfo::Ptr()->SortOut(e);
#  ifdef Lambda_p_ID
    MCInfo::Ptr()->SortLambdaOut(e,true);
#  else
    MCInfo::Ptr()->SortLambdaOut(e,false);
#  endif
  }
#endif

  //         ******************** TRIGGER ********************
  unsigned int allTrigs = 0xffff, trigMask = e.TrigMask();
  unsigned int evTrig = trigMask&allTrigs;
#ifdef U3_GENERALPURPOSE_HISTOS
  h_triggers->Fill(evTrig);
  if (Spill!=prvSpill && h_spills->GetBinContent(Spill+1)==0)
    h_spills->Fill(Spill);
#endif
  if (!e.IsMC()) {
    if (!(evTrig&U3_TRIGGER)) return;
  }
#ifdef U3_ALLOUT_KINE
  int trigType = GetTrigType(evTrig);
#endif

  // ******************** EVENT SELECTION SWITCH ********************
  bool uDSTSelection = false;

  // **************************************************************************
  //   ************************ PRIMARY VERTEX ************************
  // **************************************************************************

  // ***** Select best pV
  // ***** Bad beam rejection
  // ***** (upon cpp options) Apply cuts on:
  //  - target,
  //  - y,
  //  - Q2.

  // ***** "PRIMARY" FLAG: "primaryPat"
  // 0x1  : There is p(rimary)Vertex.
  // 0x2  : If muon beam, pV's BMS w/in range (and hence genuine BMS exists)
  // 0x4  : If muon beam, pV's backpropagation to BMS OK
  // 0x8  : If muon beam, pV muon.
  // 0x10 : If muon beam, pV's scattered mu abides by PaVertex::isMuPrim
  // 0x20 : If muon beam, pV's scattered mu has chi2/NDF w/in cut
  // 0x40 : W/in target zone
  // 0x80 : W/in target (w/, possibly, a tighter requirement than 0x40).
  // 0x100: If muonBeam, fulfills y requirements
  // 0x200: If muonBeam, fulfills Q2 cut 
  unsigned short primaryPat = 0;
  unsigned short primaryRequirements = 0x1;
#ifdef BMS_P_CUT
  if (dpBMS) primaryRequirements |= 0x2; // I.e. if width of BMS cut specified, not the case for hadron data, cf. DataTakingDB
#endif
#ifdef BMS_LH_CUT
  if (dpBMS) primaryRequirements |= 0x4; // If BMS P cut (cf. supra), makes sense to also require its reliability
#endif
#ifdef U3_TARGET_CUT
  primaryRequirements |= 0x40;
#  if U3_TARGET_CUT != 4 && U3_TARGET_CUT != 5
  primaryRequirements |= 0x80;
#  endif
#endif
  if (muonBeam) {
#ifdef U3_y_CUT 
    primaryRequirements |= 0x138; // Require mu' and y
#endif
#if defined U3_Q2_CUT && U3_Q2_CUT > 1
    primaryRequirements |= 0x238; // Require mu' and Q2
#endif
  }
  //                                                          ***** BEST pVERTEX
  // Two alternatives: CORAL'S or PHAST'S
  // - The latter considers (criterion #1) the #tracks in pV, then discrimates
  //  among pVs w/ same #tracks based on (criterion #2) chi2.
  // - The former can be better when it comes to events w/ only neutrals in
  //  addition to the scattered particle, which don't differ in criter. #1 from
  //  non-interaction vertices (i.e. fake vertices built on a non-interacting
  //  beam track and its continuation into the spectrometer): it takes into
  //  account the momentum transfered to the scattered particle.
  // - BUT it's only so in later versions of CORAL.
  int ipV; if (getBestPVfromCORAL) {
    ipV = e.iBestCoralPrimaryVertex();
    if (ipV>=0) {// Let's do some consistency check.
      if (!e.vVertex(ipV).IsPrimary()) {
	printf("\n** UB:\a Event #%d Vertex #%d is BPV while not primary!!\n",
	       EvNum,ipV); ipV = -1;
      }
    }
  }
  else ipV = e.iBestPrimaryVertex();

  // ******************** ARRAYS of TRACK/VERTEX ATTRIBUTES ********************
  // (Note: the present block, as well as the following one, must be executed 
  // before "parseBestPV", because the latter implies executing "transport" on
  // the scattered muon, which in turn implies, for CPU optimization reasons,
  // a number of checks, involving RICH and tZone info.)
  int nTrks = e.vTrack().size();
  PIDs  = new int[nTrks];
  richDths = new double[nTrks];
  for (i = 0; i<5; i++) richLHs[i]   = new double[nTrks];
  for (i = 0; i<5; i++) richChi2s[i] = new double[nTrks];
  tZones = new int[nTrks]; aRs = new int[nTrks];
  memset((void*)tZones,0,nTrks*sizeof(int));
  memset((void*)aRs,   0,nTrks*sizeof(int));
  memset((void*)PIDs,  0,nTrks*sizeof(int));

  //          ******************** PREPARE for pID ********************
  // "deltaIndex": Diff between "runIndex" and "meanIndex"(averaged over period)
  // (Note that it is then expected there be only one period per PHAST job)):
  // will serve to build a corrected Theta, corresponding to the mean index, to
  // be histogrammed in 2-D Theta vs. p 
  double deltaIndex = DataTakingDB::Ptr()->MeanIndex-runIndex;
  thCIndex = DataTakingDB::Ptr()->ThCIndex;
#ifdef U3_USE_CHCHI2
  GetPIDs(e,rich,runIndex,prodIndex,deltaIndex,thCIndex,
	  PIDs,richLHs,richChi2s,tZones,richDths);
#else
  GetPIDs(e,rich,runIndex,prodIndex,deltaIndex,thCIndex,
	  PIDs,richLHs,0,        tZones,richDths);
#endif
  //                                                         ***** PARSE BEST pV
  static double  Xp,  Yp,  Zp;  // Primary coordinates
  TMatrixD CovP(3,3);
  int imuS = -1; const PaParticle *muS = 0;
  TLorentzVector lvq, lvp(0,0,0,M_p);
  double yB = 0,  Q2 = 0, xB = 0, E0 = 0;
  if (ipV>=0)
    parseBestPV(e,ipV,primaryPat,
		Xp,Yp,Zp,CovP,lvp,E0,imuS,muS,lvq,yB,Q2,xB);

  //                                 ********** CHECK REQUIREMENTS on pVERTEX...
  h_primary->Fill(primaryPat);
  if ((primaryPat&primaryRequirements)!=primaryRequirements) {
    delete PIDs; delete tZones; delete aRs; delete richDths;
    for (i = 0; i<5; i++) { delete richLHs[i]; delete richChi2s[i]; }
    //                               ****************************************
    return;                       // ********** ... EXIT if NOT FULFILLED
    //                               ****************************************
  }

  const PaVertex &pV = e.vVertex(ipV); int nTrksPV = pV.NOutParticles();

#ifdef U3_OUTPUT_TREE
  copyCSEventHeader(e,primaryPat,spillOKPat,pV,E0,Q2,xB,yB,runIndex,prodIndex);
  fHadrons.clear(); fPa2fH.clear();
  fResonances.clear();
#endif
#ifdef U3_ALLOUT_KINE
  fillAllOutKine(e,trigType,ipV,pV,primaryPat,imuS,muS,Q2,xB,yB,lvq);
#endif

  // *****************************************************************
  // ************************* EXCLUSIVE PHI *************************
  // *****************************************************************
#ifdef U3_QUASIEXCL_PHI	// May be enabled exceptionally: to explore the case.
  bool candidatExcl = nTrksPV==3 || nTrksPV==4;
#else
  bool candidatExcl = nTrksPV==3;
#endif
  int KEMiss = 0; // K-eval'd EMiss: -1 | 0(=undefined) | 1 | 2 | 3
  if ((primaryPat&0x31)==0x31 && // Need a clean evaluation of ``exclusivity''
      // => require mu (w/in range => genuine) and mu' (ID'd, good chi2)...
      //  Yet, we also want to single out potentially exclusive events, that
      // cannot be fuly confirmed for lack of BMS or BMS out of range, so that
      // we can discard them from the inclusive search below, based on their
      // "KEMiss==3" flag. "fillExclusive" does remove them from its exclusive
      // selection, but after only "KEMiss" is set.
      //  Also let's not require "backpropagation to BMS OK" since we want to
      // keep the possibility to study its impact (e.g. what about the rate
      // of <0 missing energy events when it's dis-/en-abled?).
      //  The other pimaryPat bits can also be useful to clean our exclusive
      // sample: they are not applied either, but "primaryPat" is anyway saved
      // in the output TTree.
      candidatExcl) {   // 4 tracks = beam + scattered mu,X+,X-
    bool isphi;
    fillExclusive(e,ipV,imuS,lvp,lvq,dECuts,Ephi_pTCUT,dM_rho,dM_phi,
		  phiMassCut,KEMiss,uDSTSelection,isphi);
    if (isphi) { // Kinematics of Excl. phi
      hQ2phi->Fill(Q2); hyBphi->Fill(yB);
    }
  }
  // *****************************************************************
  // ************************* INCLUSIVE PHI *************************
  // *****************************************************************
  if ((primaryPat&0x71)==0x71 && // I.e. 0x30 (as supra) and pVertex w/in target zone...
      // ...These conditions are there to ensure a clean event reconstruction
      // for this channel that is known to be difficult
      //  => We require same conditions as above in the excl. case and in
      //    addition add target zone.
      //     Note that anyway the selection criteria required above in the
      //    exclsuive case need be fulfilled for "KEMiss" to be meaningfull.
      //  To prevent double counting w/ the exclusive case
      //nTrksPV>3 && // Note: this somehow turned out not to be needed...
      //  For further avoiding double counting, we have to consider the case of
      // quasi-exclusive phi, obtained supra when "nTrksPV==4" while enabling
      // "U3_QUASIEXCL_PHI": this is not done fully rigorously, since in
      // quasi-exclusive, we only allow the pair of ''hadrons'' w/ highest P to
      // build a phi, whereas here we would also consider the other two pairs.)
      KEMiss!=1 && KEMiss!=2)
    fillIphi(e,ipV,imuS,dM_phi,uDSTSelection);

#ifdef U3_FILL_Lambda
  // Now that cuts on muon kine are passed: update internal state of "MCLambda"
  if (e.IsMC()) LambdaMC->Update();
#endif

  //           ********** gamm*N(=p) CENTER of MASS **********
  static double sqrts; TVector3 betacm, vqUnit;
  if (primaryPat&0x10) {
    TLorentzVector lvcm = lvq+lvp; sqrts = lvcm.M();
    betacm = -lvcm.BoostVector(); vqUnit = lvq.Vect().Unit();
    //#define U3_DEBUG_CM_BOOST
#ifdef U3_DEBUG_CM_BOOST
    lvcm.Boost(betacm); lvcm.Vect().Print("\nCM");
#endif
  }


  // *****************************************************************
  // *************************  V0 SELECTION *************************
  // *****************************************************************

#ifdef DEBUG_K0
  dbK0GetEvent(e);  // Is K0 expected in current event?...
#endif
  // ***** V0 SELECTION PATTERN
  unsigned short DdDOK =  0x01;      // 0x1:  Vertex distance
  unsigned short collOK = 0x02;      // 0x2:  Collinearity
  unsigned short pTOK =   0x04;      // 0x4:  pT > cut
  unsigned short chi2OK = 0x08;      // 0x8:  Chi2 of Vertex
  //                                    0x30: h1/2 w/in RICH acc. and P range
  //                                    0xc0: PID of h1/2
  unsigned short tightOK = 0x100;    // 0x100: Tighter D/dD, Zs>ZTarget
  //                                    0x600: eVeto of h1/2

  int isV; for (isV = 0; isV<e.NVertex(); isV++) {
    const PaVertex& sV = e.vVertex(isV);

    //                                     ******************** LOOP ON VERTICES

    if (sV.IsPrimary()) continue;                      // ***** BYPASS PRIMARIES
    if (sV.NOutParticles()!=2) {           // ***** BYPASS BADLY FORMED 2NDARIES
      printf("\n** U3:\a PaVertex Inconsistency: NOutParticles = %d\n",
	     sV.NOutParticles());
      continue;
    }

    double Xs = sV.Pos(0), Ys = sV.Pos(1), Zs = sV.Pos(2);

    //                                        ***** GET TRACKS, HELICES @ VERTEX
    int iEP1 = sV.iOutParticle(0), iEP2 = sV.iOutParticle(1);
    if (ipV>=0 && (iEP1==imuS || iEP2==imuS))
      continue;                                 // ***** (if exists) EXCLUDE mu'
    const PaParticle *pa1 = &e.vParticle(iEP1), *pa2 = &e.vParticle(iEP2);
    if (pa1->NFitPar()==0 || pa2->NFitPar()== 0 ||
	pa1->Q()==pa2->Q()) {
      printf("\n** U3:\a h+h- Inconsistency: Evt %d#%d Vtx %d = (%d,%d)   %d %d\n",
	     Run,EvNum,isV,
	     pa1->NFitPar(),pa2->NFitPar(),pa1->Q(),pa2->Q());
      continue;
    }
    if (pa1->Q()<pa2->Q()) swap(pa1,pa2);      // ***** FIRST PARTICLE = POSITIVE

    int iET1 = pa1->iTrack(), iET2 = pa2->iTrack();
    PaTrack &trk1 = e.vTrack()[iET1], &trk2 = e.vTrack()[iET2];
    double zL1 = trk1.ZLast(), zL2 = trk2.ZLast();

#ifdef DEBUG_K0
    if (dbK0Event>=0) {
      const PaTPar *h; int i12; 
      for (i12 = 0, h = &trk1.vTPar()[0], dbK0p1p2 = 0; i12<2; i12++) {
	if (fabs(h->Mom()-dbK0p[i12])<.1*dbK0p[i12] &&
	    fabs(h->Theta()-dbK0a[i12])<.1*fabs(dbK0a[i12])) dbK0p1p2 |= 1<<i12;
	h = &trk2.vTPar()[0];
      }
      dbK0Event |= dbK0p1p2;
    }
#endif

    int zones1 = tZones[iET1], zones2 = tZones[iET2];
    if (!(zones1&0x2) || !(zones2&0x2))            // ***** EXCLUDE FRINGE FIELD
      continue; // (Note: Excluding also possible 0x11 tracks...)
#ifdef DISCARD_MUONS                           // ***** (OPTIONAL) EXCLUDE MUONS
    if (trk1.XX0()>15 || trk2.XX0()>15) continue;
#endif
    if (trk1.Chi2tot()/trk1.Ndf()>chi2Cut ||            // ***** CUT ON CHI2/NDF
	trk2.Chi2tot()/trk2.Ndf()>chi2Cut) continue;

    //                                          ***** SOME MORE TRACK ATTRIBUTES
    float time1 = trk1.MeanTime(), time2 = trk2.MeanTime();
    // Determine acceptance by RICH and exit angle for SM1 and incidence/exit
    // angle for SM2.
    // (Note that if Marcin's rescaling is in effect (cf. the definition of
    // "SCALE_MOMENTA" supra), the rescaling is already done. However, this
    // is limited to an 0.8% effect in SM2 alone. It is then not expected to
    // affect the extrapolation to the RICH, given the small SM1 bending
    // strength. For what concerns the extrapolation to downstream of SM2, which
    // is also performed by "transport", for track bridged over SM2, it's done
    // relative to the end helix, and hence insensitive to any magnetic effect.
    // This, provided that the end helix is available, though: which is the case
    // in the latest mass productions.)
    if (aRs[iET1]==0) /* If not yet done */ transport(e,trk1,aRs[iET1],zones1);
    if (aRs[iET2]==0)                       transport(e,trk2,aRs[iET2],zones2);
    if (aRs[iET1]<-1 || aRs[iET2]<-1)
      continue; // Skip if either track starting downstream of RICH
    const PaTPar &hv1 = pa1->ParInVtx(isV), &hv2 = pa2->ParInVtx(isV);
    TVector3 v13 = hv1.Mom3(), v23 = hv2.Mom3();
    // For P @ RICH, best would be a smoothing point in zone 0x2 (note that we
    // do it in "copyCSHadronData"). Next best...
    const PaTPar &hi1 = trk1.vTPar()[0], &hi2 = trk2.vTPar()[0];
    double mom1 = fabs(1/hi1(5)), mom2 = fabs(1/hi2(5));

    //    *************** PI+PI- MASS ASSUMPTION ***************
    double E1 = sqrt(v13.Mag2()+M2_pi), E2 = sqrt(v23.Mag2()+M2_pi);
    TLorentzVector lvpi1(v13,E1), lvpi2(v23,E2), lvpipi = lvpi1+lvpi2;
    double pT = lvpi1.Perp(lvpipi.Vect()), m_pipi = lvpipi.M();
    double PL1 = sqrt(v13.Mag2()-pT*pT), PL2 = sqrt(v23.Mag2()-pT*pT);
    double alpha = (PL1-PL2)/(PL1+PL2);
    
    fillK0woPV(m_pipi,Xs,Ys,Zs,pT,zL1,zL2,lvpipi);  // ***** K0 W/O pV CONDITION

    // ******************** REQUIRE PRIMARY => VERTEX LINE ********************

    TMatrixD CovS(3,3);
    for (i=k = 0; i<3; i++) for (j = 0; j<=i; j++, k++) {
	CovS(i,j) = sV.Cov(k); if (i!=j) CovS(j,i) = sV.Cov(k);
      }
    TVector3  vv3(Xs-Xp,Ys-Yp,Zs-Zp);	// Vertex Line
    double vvv0 = vv3*lvpipi.Vect(), dist = vv3.Mag();
    double ctheta = vvv0/dist/lvpipi.Vect().Mag();
    dist = ctheta>0 ? dist : -dist;

    double ddist = 1; // dctheta = 1;
    {
      TVectorD vv(0,2,Xs-Xp,Ys-Yp,Zs-Zp,"END");
      TVectorD tmps = vv, tmpp = vv;   tmps *= CovS; tmpp *= CovP;
      ddist = tmps*vv + tmpp*vv; ddist = sqrt(ddist)/fabs(dist);
      /*
	TVectorD ddx = vv; ddx *= -vvv0/lvpipi.Mag2();
	TVectorD tmp0(0,2,lvpipi(0),lvpipi(1),lvpipi(2),"END"); ddx += tmp0;
	tmps = ddx; tmpp = ddx;         tmps *= CovS; tmpp *= CovP;
	dctheta = tmps*ddx+tmpp*ddx;
	if (0) {
	  TVectorD ddv = vv; tmp0 *= -vvv0/lvpipi.Mag2(); ddv += tmp0;
	  TVectorD tmp1 = ddv, tmp2 = ddv; tmp1 *= CovS[1]; tmp2 *= CovS[2];
	  dctheta += tmp1*ddv+tmp2*ddv;
	  TVectorD tmp12 = ddv;            tmp12 *= CovS[5];
	  dctheta += 2*(tmp12*ddv);
	  tmp1 = ddv; tmp2 = ddv;         tmp1 *= CovS[3]; tmp2 *= CovS[4];
	  dctheta += 2*(tmp1*ddx+tmp2*ddx);
	  dctheta = sqrt(dctheta)/lvpipi.Mag()/dist;
	}
      */
    }

    //   ********** V0 SELECTION CRITERIA ********** 
    //                                              ***** INDICES of FULFILLMENT
    double chi2_s = sV.Chi2()/sV.Ndf();
    int s_dist = 0, s_ctheta = 0, s_pT = 0, s_chi2 = 0;
    //int s_t = 0; // Time selection: not used yet
    if      (ctheta>V0cthCut[2])   s_ctheta = 3;
    else if (ctheta>V0cthCut[1])   s_ctheta = 2;
    else if (ctheta>V0cthCut[0])   s_ctheta = 1;
    if      (pT>V0pTCut[1])        s_pT = 2;
    else if (pT>V0pTCut[0])        s_pT = 1;
    if      (chi2_s<10) s_chi2 = 2; // Exceedingly loose, since chi2 is per NDF
    else if (chi2_s<15) s_chi2 = 1;
    if      (dist/ddist>V0DsdD[2]) s_dist = 3;
    else if (dist/ddist>V0DsdD[1]) s_dist = 2;
    else if (dist/ddist>V0DsdD[0]) s_dist = 1;
    //                     ***** DISTRIBUTION HISTOS for K0 SIGNAL and SIDEBANDS
    if      (fabs(m_pipi-M_K0)<.008) {
      if (s_ctheta && s_pT     && s_chi2)   hs_dK0->Fill(dist/ddist);
      if (s_pT     && s_chi2   && s_dist)   hs_aK0->Fill(ctheta);
      if (s_chi2   && s_dist   && s_ctheta) hs_pK0->Fill(pT);
      if (s_dist   && s_ctheta && s_pT)     hs_vK0->Fill(chi2_s);
      if (s_ctheta && s_pT     && s_chi2  && s_dist) {
	hs_tK0->Fill(time1); hs_tK0->Fill(time2);
      }
    }
    else if (fabs(m_pipi-(M_K0-5.25*.008))<.004 ||
	     fabs(m_pipi-(M_K0+5.25*.008))<.004) {
      if (s_ctheta && s_pT     && s_chi2)   hn_dK0->Fill(dist/ddist);
      if (s_pT     && s_chi2   && s_dist)   hn_aK0->Fill(ctheta);
      if (s_chi2   && s_dist   && s_ctheta) hn_pK0->Fill(pT);
      if (s_dist   && s_ctheta && s_pT)     hn_vK0->Fill(chi2_s);
      if (s_ctheta && s_pT     && s_chi2  && s_dist) {
	hn_tK0->Fill(time1); hn_tK0->Fill(time2);
      }
    }

#ifdef DEBUG_K0
    if (dbK0p1p2==0x3) dbK0Event |= 0x4;
#endif
#ifdef U3_K0_HIGHLEVEL  // K0 vs. Incidence and related
    double Theta, Thetx, Thety, ThetX, alpha, beta;
    getK0Incidence(trk1,trk2,zones1,zones2,v13,v23,hi1,hi2,
		   Theta,Thetx,Thety,ThetX,alpha,beta);
    if      (fabs(m_pipi-M_K0)<.008)             hs_TK0->Fill(Theta);
    else if (fabs(m_pipi-(M_K0-5.25*.008))<.004 ||
	     fabs(m_pipi-(M_K0+5.25*.008))<.004) hn_TK0->Fill(Theta);
#endif

#ifdef U3_ALLOUT_KINE
    if (s_pT>=1 && s_chi2>=2 && s_dist>=3 && s_ctheta>=2) {
      //                          ***** VERY STRICT V0 SELECTION => Z OF SVERTEX
      // (This is meant to evidence possible flocking of V0's at some particular
      // Zs corresponding to actual primary- or re-interactions.)
      hv_sZV0->Fill(Zs);
    }
#endif

    // ****************************************************************
    // ****************************** K0 ******************************
    // ****************************************************************

    //                                      ********** LOOSE V0 SELECTION for K0
    if (s_pT>=1 && s_chi2>=2 &&      // ***** CUT on PT(loose), V CHI2 and...
	fabs(m_pipi-M_K0)<dM_K0 &&   // ***** ...K0 HISTO MASS RANGE
	s_dist>=1 && s_ctheta>=1) {  // ***** ...(CRUDE) VERTEX LINE
      unsigned short K0Pat = pTOK|chi2OK;

      uDSTSelection = true;                     // ********** uDST SELECTION: K0

      //                                                      ***** EVALUATE PID
      unsigned short id = 0; // 0x9: piThr<P<Pmax, 0x12: !e+e-, 0x24: actual piID
      unsigned short idpipOrpim = 0, isNotee = 0;
      int pm, iET; double mom; for (pm = 0, iET = iET1, mom = mom1; pm<2;
				    pm++) {
	if (aRs[iET]>0) {
	  bool winRange = piThr<mom && mom<PRICHCut;
	  if (winRange) {
	    K0Pat |= (0x1<<pm)<<4; id |= 0x1<<(pm*3);
	  }
	  if (PIDs[iET]&0x10) {
	    double eLH = richLHs[3][iET], piLH = richLHs[0][iET];
	    double KLH = richLHs[1][iET], pLH =  richLHs[2][iET];
	    double hadronLH = KLH>pLH ? KLH : pLH;
	    if (eLH==0 || // eLH=piLH=0: eVeto even above piThr!
		eLH<=piLHeVeto*piLH) {
	      // eVeto: don't require piThr<P<Pmax. It's:
	      // - Unnecessary: "piLHeVeto" is expected to be large => P>Pmax,
	      //  where eLH=piLH, is anyway rejected,
	      // - Harmful: P<piThr can be efficiently vetoed.
	      // (Note that piID (i.e. id&0x4) here does not place any condition
	      // on pi/eLH, which let well ID'd e+/e- pass through. One could
	      // consider a cut on pi/eLH in the region below asymptotic pions,
	      // to correct for this. But we want here to arbitrate the purity
	      // vs. efficiency trade-off in favour of pion, efficiency, since
	      // pion yield is large compared to e+e- and reasonable purity is
	      // kind of automatically granted. Anyway, a consequence is that an
	      // OR of pionID (id&0x4) and eVeto yields more than eVeto alone.)
	      id |= 0x2<<(pm*3); isNotee |= 0x1<<pm;
	    }
	    if (winRange && piLH>LHCut*hadronLH && piLH>LHBckCut) {
	      id |= 0x4<<(pm*3); idpipOrpim |= 0x1<<pm;
	    }
	  }
	  else {
	    // Absence of RICH block or LH =0: eVeto, even above piThr, cf. comment in Lambda PID
	    id |= 0x2<<(pm*3); isNotee |= 0x1<<pm;
	  }
	}
	iET = iET2; mom = mom2;
      }
      K0Pat |= idpipOrpim<<6; K0Pat |= isNotee<<9;

      if (s_dist>=1 && s_ctheta>=1)
	hm_K0c11->Fill(m_pipi);              // ********** ALL K0s HISTOGRAMMING
      if (s_dist>=1 && s_ctheta>=2)
	hm_K0c12->Fill(m_pipi);
      if (s_dist>=2 && s_ctheta>=1)
	hm_K0c21->Fill(m_pipi);

#ifdef U3_ERASE_K0s
      // Let's erase as many pions as possible from the pVertex (in order to
      // reduce its combinatorics) => Loose cuts (particularly important being
      // a loose setting for D/dD, which allows to catch pions attached to the
      // pVertex despite the fact they come from a V0 decay).
      if (fabs(m_pipi-M_K0)<K0MassCut) {
	PaParticle &pi1 = const_cast<PaParticle&>(e.vParticle(iEP1));
	PaParticle &pi2 = const_cast<PaParticle&>(e.vParticle(iEP2));
	pi1.SetPID(0); pi2.SetPID(0); // Tag decay pions
      }
#endif

      //          ********** STRICT V0 SELECTION for K0: STRICT VERTEX LINE CUTS
      if (s_dist>=2 && s_ctheta>=2) {
	K0Pat |= DdDOK; K0Pat |= collOK;

	//                                 ********** STRICT V0 K0-HISTOGRAMMING
	hm_K0c22->Fill(m_pipi);

	if (primaryPat&0x10)                                  // ***** pV HAS MU'
	  hm_K0C22->Fill(m_pipi);
	//                         ********* LAS and SAS K0s **********
	if      (zL1<ZSM2 && zL2<ZSM2) hm_K0cll->Fill(m_pipi);    // ***** LAS^2
	else if (zL1<ZSM2 || zL2<ZSM2) hm_K0csl->Fill(m_pipi);  // ***** LASxSAS
	else                           hm_K0css->Fill(m_pipi);    // ***** SAS^2

#ifdef U3_K0_HIGHLEVEL  // K0 vs. Incidence and related
	fillK0HighLevel(trk1,trk2,m_pipi,lvpipi,
			Theta,Thetx,Thety,ThetX,alpha,beta);
#endif
	//                                               ***** IN/OUTSIDE TARGET
	// I.e. single out those V0's that cannot be confused w/ interactions
	// w/in the physics target.
	// (Note: At pesent, there's info on the target extension only if
	// DataTakingDB was invoked w/ targetOption =1.)
	if (DataTakingDB::Ptr()->TargetOption==1) {
	  if (Zs<DataTakingDB::Ptr()->ZD_2) hm_K0i->Fill(m_pipi);
	  else                              hm_K0o->Fill(m_pipi);
	}
	hm_K0VsID->Fill(m_pipi,id);                   // ***** pipi MASS vs. PID

	/*
	// Examples of use of hm_K0VsID in ROOT:
	TH2D *H = (TH2D*)gDirectory->Get("hm_K0VsID");
	int nBinsX = H->GetNbinsX(); double xMn = H->GetBinLowEdge(1);
	double xMx = H->GetBinLowEdge(nBinsX)+H->GetBinWidth(nBinsX);
	unsigned int id, jd; int bin, strt, pm, ok;
	TH1D *Hi; char HiN[] = "h64_64";
	// "idOK", i.e. piID or subpThr elseVeto = id&0x4||id0x2
	TH1D *HOK = new TH1D("hm_K0OK","K0 id&0x4||id&0x2",nBinsX,xMn,xMx);

	for(bin=1,strt=-1;bin<=64;bin++){for(pm=0,ok=1,id=bin-1;pm<2;pm++){id>>=3*pm;jd=id&0x7;if(!((jd&0x4)||(jd&0x2)))ok=0;};if(ok){printf("  %d 0x%x(0x%x,0x%x)",bin,bin-1,(bin-1)&0x7,(bin-1)>>3);if (strt<0)strt=bin;}else if (strt>=0){sprintf(HiN,"h%d_%d",strt,bin-1);if((hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();hi=H->ProjectionX(HiN,strt,bin-1);HOK->Add(hi);printf("  Add [%d,%d]",strt,bin-1);strt=-1;}};if (ok) {sprintf(HiN,"h%d_%d",strt,64);if((hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();hi=H->ProjectionX(HiN,strt,64);HOK->Add(hi);printf("  Add [%d,%d]\n",strt,64);} else printf("\n");
	// In order to condition only p (or pi), simply change the bounds of the 2nd for-loop.
	// In order to have an OR instead of an AND, initialize "id=0x3" and reset it w/ "id^=1<<pm".
	*/

	bool cleanV0; if (Zs>ZTarget && s_dist>=3) {     // ***** TAG CLEAN V0'S
	  cleanV0 = true; K0Pat |= tightOK;
	}
	else cleanV0 = false;


#ifdef DEBUG_K0
	if (dbK0p1p2==0x3) { dbK0Event |= 0x8; dbK0v1 = m_pipi; }
#endif

	bool isK0 = fabs(m_pipi-M_K0)<K0MassCut;               // ***** TAG K0'S

	if (isK0&& (idpipOrpim || isNotee)) {
	  //                         ***** DECAY P and (X,Y), THETA vs. P @ RICH
	  // (This is to visualize the K0 kinematical domain for RICH perfs:
	  // there, when targeting pi+/-, we allow ourselves to ID pi-/+.)
	  const PaTrack *trk; double mom; const TLorentzVector *lvpi;
	  for (pm = 0, iET = iET1, trk = &trk1, mom = mom1, lvpi = &lvpi1; pm<2;
	       pm++) {	      
	    if ((idpipOrpim&1<<(1-pm)) || (isNotee&1<<(1-pm))) {
	      int mp = 1-2*pm; hk_PK0->Fill(mom,mp);
	      const vector<Float_t> &aux = trk->vAux();
	      if (aRs[iET]<-1) { // Ensure then "aRs" are meaningful
		printf("** U3:\a Evt %d,%d Track %d(<-K0): Inconsistency: aRs = %d\n",
		       e.RunNum(),(int)e.UniqueEvNum(),iET,aRs[iET]);
		abort();
	      }
	      float tgx = aux[2], tgy = aux[3];
	      double thRICH = acos(1/sqrt(1.+ tgx*tgx + tgy*tgy));
	      hk_PThK0->Fill(mom,thRICH);
	      hk_XK0->Fill(aux[0],mp); hk_YK0->Fill(aux[1],mp);
	      if (s_dist>=3 && s_ctheta>=3 && s_pT>=2 &&
		  // Exclude pion decaying to muon, possibly upstream of RICH
		  trk->XX0()<15 &&
		  PIDs[iET]&0x10) {
		double thC = trk->RichInf(thCIndex);
		// W/o and w/ index evolution.
		hR_thCPK0->Fill(mp*mom,thC);
		double beta = lvpi->Beta(), thCorr = thC+richDths[iET]/beta;
		hR_ThCPK0->Fill(mp*mom,thCorr);
		double idx = DataTakingDB::Ptr()->MeanIndex;
		double thCpi = 1000*acos(1/beta/idx);
		hR_THCPK0->Fill(mp*mom,thCorr-thCpi);
	      }
	    }
	    iET = iET2; trk = &trk2; mom = mom2; lvpi = &lvpi2;
	  }
#ifdef U3_ALLOUT_KINE
	  if (primaryPat&0x10) {      // pV HAS MU', pi+|pi-ID
	    hQ2K0->Fill(Q2,trigType); hyBK0->Fill(yB,trigType);
	  }
#endif
	}

	// ************************************************************
	// ********** RICH PERFORMANCES (and related) from K0 *********
	// ************************************************************
#ifdef U3_FILL_RICH_PERFS
#  ifdef DEBUG_K0
	if (dbK0p1p2==0x3 && cleanV0) {
	  dbK0Event |= 0x10; dbK0v1 = Zs; dbK0v2 = dist/ddist;
	}
#  endif
	fillK0_RICHPerf(e,iET1,iET2,lvpi1,lvpi2,m_pipi,isK0,cleanV0,Run);
#endif

	if (!isNotee) { // h1|h2 ID'd as e+|e-
#ifdef piS_e_REJECTION
	  isK0 = false;          // ***** OPTIONAL piS_e_REJECTION => "isK0" = 0
#endif
	}
	else hm_K0eV->Fill(m_pipi); 

	if (isK0) {
	  // ************************************************************
	  // ********************    K0 SELECTION    ********************
	  // ***************    - MASS CUT           ********************
	  // ***************    - (OPTIONALLY) eVETO ********************
	  // ************************************************************
	  double EK0 = sqrt(M2_K0+lvpipi.Vect().Mag2());  // Assuming M_K0
	  TLorentzVector lvK0(lvpipi.Vect(),EK0);

	  fillK0X(e,ipV,imuS,lvq,                                // ***** K0 + X
		  iET1,iET2,lvK0,KSMassCut,dM_phi,phiMassCut);

#ifdef U3_ALLOUT_KINE
	  // ***** SOME MORE HISTOS OF K0 CHARACTERISTICS *****
	  // ...   UNRELATED to RICH PERFS
	  const PaParticle &beam = e.vParticle(pV.InParticle());
	  const PaTPar &hv = beam.ParInVtx(ipV);
	  double R = sqrt(Xp*Xp+Yp*Yp), dRdZ = (Xp*hv(3)+Yp*hv(4))/R;
	  hv_KdRvsR->Fill(R,dRdZ);
#endif
	}  // End Strict K0 selection
      }  // End strict V0 selection for K0

#ifdef U3_OUTPUT_TREE
      //                                                  ***** FILL OUTPUT TREE
      copyCSResonanceData(sV,m_pipi,dist,ddist,ctheta,pT,alpha);
      CSResonanceData &r = fResonances.back(); r.K0Pat = K0Pat;
      const PaParticle *pa; for (pm = 0, pa = pa1; pm<2; pm++) {
	int fHInd = getHadronIndex(pa);
	if (pm==0) r.h1 = fHInd; //adrons.size();
	else       r.h2 = fHInd; //adrons.size();
	if (fHInd==(int)fHadrons.size()) copyCSHadronData(e,pa,isV);
	pa = pa2;
      }
#endif
    } // End loose V0 selection for K0
#ifdef U3_FILL_Lambda
    // ****************************************************************
    // **************************** Lambda ****************************
    // ****************************************************************

    //      *************** P,PI PI,P MASS ASSUMTIONS ***************
    E1 = sqrt(v13.Mag2()+M2_p); E2 = sqrt(v23.Mag2()+M2_p);
    TLorentzVector lvp1(v13,E1), lvp2(v23,E2);
    TLorentzVector lvppi = lvp1 + lvpi2, lvpip = lvpi1 + lvp2;
    double m_ppi = lvppi.M(), m_pip = lvpip.M();

    int winLRange = 0;
    if      (fabs(m_ppi-M_Lam)<dM_Lambda) winLRange |= 0x1;
    else if (fabs(m_pip-M_Lam)<dM_Lambda) winLRange |= 0x2;

    //                               ********** PRELIMINARY SELECTION for Lambda
    if (s_pT>=1 && s_chi2>=1 &&        // ***** LOOSE CUT on pT, V CHI2 and...
	winLRange) {                   // ***** ...[a]Lambda HISTO MASS RANGE


      int iLaL;
      //                                                      ***** EVALUATE PID
      unsigned short ppiIDs[2] = {0,0}; // 0x9: piTHr<P<Pmax, 0x12: subThr, 0x24: actual (p|pi)ID
      unsigned short idpOrms[2] = {0,0};// ID above/below p threshold, 0x1:p(lus), 0x2: m(minus)
      unsigned short RICHOKs[2] = {0,0}, isNotees[2] = {0,0};
      int pm, iET; double mom; for (pm = 0, iET = iET1, mom = mom1; pm<2;
				    pm++) { 
	if (aRs[iET]==0) {
	  printf("** U3:\a Evt %d,%d Track %d(<-Lambda) RICH not yet checked\n\n",
		 e.RunNum(),(int)e.UniqueEvNum(),iET); assert(false);
	}
	for (iLaL = 0; iLaL<2; iLaL++) {
	  if (aRs[iET]<0) continue;
	  unsigned short idpm = 0;
	  if (iLaL==pm) {                   // Evaluate pID
	    if (mom<piThr*subpThrPmin || pIDPmax<mom) continue;
	    RICHOKs[iLaL] |= 0x1<<pm;
#  ifdef U3_IDEAL_RICH
	    idpm |= 0x6;
#  endif
	    idpm |= 0x1;
	    if (PIDs[iET]&0x10) {// => LHs are availble
	      double piLH = richLHs[0][iET], KLH = richLHs[1][iET];
	      double pLH  = richLHs[2][iET], eLH = richLHs[3][iET];
	      double elseLH = eLH>piLH?eLH:piLH; if (KLH>elseLH) elseLH = KLH;
	      if (pThr*pIDPmin<mom) {
		// "pIDPmin" extra tuning knob that we do not have for the KID
		if (pLH>LHCut*elseLH && pLH>LHBckCut) idpm |= 0x4;
	      }
	      else {
		bool veto = eLH>subpThrLHVeto || piLH>subpThrLHpiVeto;
		if (KThr<mom) {// The momentum cut here need not be precise: it only cuts out that relatively small fraction of background that the K's constitute => therefore no need to be efficient.
		  veto |= KLH>subpThrLHVeto;
		}
		if (!veto) idpm |= 0x2;
	      }
	    }
	    else if (mom<pThr*pIDPmin)
	      // - Absence of RICH block, i.e. "!(PIDs[iET]&0x8)" means, in
	      // older RICH package, no photon => particle is below threshold.
	      // - In addition, one may surmise that backgr. LH =0, i.e.
	      //  "!(PIDs[iET]&0x10)", also means particle is below threshold.
	      //  To be conditioned yet by w/in RICH acceptance, granted (
	      //  that's still depending upon acceptance description, e.g.
	      //  "rRICH") if "aR>0", and w/in the reach of RICH software,
	      //  assumed to imply ZFirst<ZRICH, also granted by "aR>0".
	      idpm |= 0x2;
	    if ((idpm&0x4) ||     // Actual pID
		(idpm&0x2))       // SubpThreshold else Veto
	      idpOrms[iLaL] |= 0x1<<pm;
	  }
	  else {                            // Evaluate piID
	    bool winRange = piThr<mom && mom<PRICHCut;
	    if (winRange) { RICHOKs[iLaL] |= 0x1<<pm; idpm |= 0x1; }
	    if (PIDs[iET]&0x10) {
	      double eLH = richLHs[3][iET], piLH = richLHs[0][iET];
	      double KLH = richLHs[1][iET], pLH =  richLHs[2][iET];
	      double hadronLH = KLH>pLH ? KLH : pLH;
	      if (eLH==0 || // eLH=piLH=0: eVeto even above piThr!
		  eLH<piLHeVeto*piLH) {
		idpm |= 0x2; isNotees[iLaL] |= 0x1<<pm;
	      }
	      if (winRange && piLH>LHCut*hadronLH && piLH>LHBckCut) idpm |= 0x4;
	    }
	    else // Absence of RICH block or LH =0: eVeto even above piThr, cf. supra
	      idpm |= 0x2;
	    if (idpm&0x4 || // Actual piID
		idpm&0x2)   // eVeto
	      idpOrms[iLaL] |= 1<<pm;
	  }
	  ppiIDs[iLaL] |= idpm<<(3*(pm*(1-iLaL)+(1-pm)*iLaL));
	} // End loop over Lambda/antiLambda
	iET = iET2; mom = mom2;
      } // End loop on >0, <0 decay particles
#  ifdef Lambda_p_ID  // ***** pID (pID or e,pi,KVETO BELOW THRESHOLD)... *****
      int papID; for (iLaL = 0, papID = 0; iLaL<2; iLaL++)
		   if (idpOrms[iLaL]&0x1<<iLaL) papID |= 1<<iLaL;
#  else
      int papID = 0x3;
#  endif
#  ifdef piS_e_REJECTION       // ***** piID: e-REJECTION *****
      for (iLaL = 0; iLaL<2; iLaL++)
	if (!(isNotees[iLaL]&0x1<<(1-iLaL)))
	  papID &= 1<<(1-iLaL); // I.e. retain complementary of current "iLaL"
#  endif

      // ***** DISTRIBUTION HISTOS for Lambda (!aLambda) SIGNAL and SIDEBANDS
      // (Note: Contrary to what's done in the V0 case, where the  distributions
      // are histo'd prior to any V0 cut, here it's done w/in the preliminary
      // "LOOSE V0 SELECTION for Lambda" and pID: else the Lambda is too much
      // diluted.)
      if      ((papID&0x1) && fabs(m_ppi-M_Lam)<.0025) {
	if (s_ctheta && s_pT     && s_chi2)   hs_dL->Fill(dist/ddist);
	if (s_pT     && s_chi2   && s_dist)   hs_aL->Fill(ctheta);
	if (s_chi2   && s_dist   && s_ctheta) hs_pL->Fill(pT);
	if (s_dist   && s_ctheta && s_pT)     hs_vL->Fill(chi2_s);
	if (s_ctheta && s_pT     && s_chi2  && s_dist) {
	  hs_tL->Fill(time1); hs_tL->Fill(time2);
	}
      }
      else if ((papID&0x1) && (fabs(m_ppi-(M_Lam-5.25*.0025))<.00125 ||
			       fabs(m_ppi-(M_Lam+5.25*.0025))<.00125)) {
	if (s_ctheta && s_pT     && s_chi2)   hn_dL->Fill(dist/ddist);
	if (s_pT     && s_chi2   && s_dist)   hn_aL->Fill(ctheta);
	if (s_chi2   && s_dist   && s_ctheta) hn_pL->Fill(pT);
	if (s_dist   && s_ctheta && s_pT)     hn_vL->Fill(chi2_s);
	if (s_ctheta && s_pT     && s_chi2  && s_dist) {
	  hn_tL->Fill(time1); hn_tL->Fill(time2);
	}
      }

      unsigned short LambdaPat = 0;

      //                                ********** LOOSE V0 SELECTION for Lambda
      if (s_pT>=2 && s_chi2>=2 &&       // ***** STRICT CUT on pT, V CHI2 and...
	  s_dist>=1 && s_ctheta>=1) {   // ***** ...(CRUDE) VERTEX LINE
	LambdaPat = pTOK|chi2OK;

	uDSTSelection = true;               // ********** uDST SELECTION: Lambda

	static bool firstLambda = true; if (firstLambda) {
	  firstLambda = false; // ********** PRINT a SUMMARY of CUTS **********
	  printLambdaCuts(V0DsdD[1],V0cthCut[1],V0pTCut[1],10.);
	}
	FillLambda(m_ppi,m_pip,                   // ***** Lambda HISTOs *****
		   s_dist,s_ctheta,primaryPat&0x10,trk1,trk2,Zs>ZTarget+65,papID);

#  ifdef U3_L_HIGHLEVEL      // Lambdas w/ either p or pi in pV
	if (s_dist>=2 && s_ctheta>=2 && (papID&0x1)) {
	  int ppiInpV = 0;
	  for (int iv1 = 0; iv1<pa1->NVertex(); iv1++)
	    if (pa1->iVertex(iv1)==ipV) { ppiInpV |= 0x1; break; }
	  for (int iv2 = 0; iv2<pa2->NVertex(); iv2++)
	    if (pa2->iVertex(iv2)==ipV) { ppiInpV |= 0x2; break; }
	  if (ppiInpV) FillLambda(m_ppi,0x2);
	}
	//                      Lambdas w/ p outside RICH, outside K0
	if (s_dist>=2 && // (Would be better if!) Tight V0 vertex dist cut
	    s_ctheta>=2 && !(papID&0x1) && aRs[iET1]<0) {
	  FillLambda(m_ppi,0x4);
	  if (fabs(m_pipi-M_K0)>.018) // 3 sigma cut, assuming sigma = 6 MeV
	    FillLambda(m_ppi,0x8);
	}
	//                      SAT ot LAT at the exit of SM2?
	if (s_dist>=2 && s_ctheta>=2 && (papID&0x1) &&
	    (zones1&0x4) && (zones2&0x4)) {
	  fillLGM0709(m_ppi,trk1,trk2,gm0709Ws,gm0709Bs);
	}
#  endif
      }

      //                               ********** STRICT V0 SELECTION for Lambda
      if (s_pT>=2 && s_chi2>=2 &&       // ***** STRICT CUT on PT, V CHI2 and...
	  s_dist>=2 && s_ctheta>=2) {   // ***** ...(STRICT) VERTEX LINE
	// Cf. "hm_[a]Lc22" in "FillLambda"
	LambdaPat |= DdDOK; LambdaPat |= collOK;

	for (iLaL = 0; iLaL<2; iLaL++) {         // ***** ppi (pip) MASS vs. PID
	  if (!(winLRange&1<<iLaL)) continue;
	  hm_LVsID[iLaL]->Fill(iLaL==0 ? m_ppi : m_pip,ppiIDs[iLaL]);
	}
	/*
	// Examples of use of hm_LVsID in ROOT:
	TH2D *H = (TH2D*)gDirectory->Get("hm_LVsID"); TAxis *ax = H->GetXaxis();
	int nBinsX = H->GetNbinsX(); double xMn = ax->GetBinLowEdge(1);
	double xMx = ax->GetBinLowEdge(nBinsX)+ax->GetBinWidth(nBinsX);
	unsigned int id, jd; int bin, strt, ppi, ok;
	TH1D *Hi; char HiN[] = "h64_64";
	// "idOK", i.e. pID or elseVeto subpThr = id&0x4||id0x2
	TH1D *HOK = new TH1D("hm_LOK","L id&0x4||id&0x2",nBinsX,xMn,xMx);

	for(bin=1,strt=-1;bin<=64;bin++){for(ppi=0,ok=1,id=bin-1;ppi<1;ppi++){id>>=3*ppi;jd=id&0x7;if(!((jd&0x4)||(jd&0x2)))ok=0;};if(ok){printf("  %d 0x%x(0x%x,0x%x)",bin,bin-1,(bin-1)&0x7,(bin-1)>>3);if (strt<0)strt=bin;}else if (strt>=0){sprintf(HiN,"h%d_%d",strt,bin-1);if((Hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();Hi=H->ProjectionX(HiN,strt,bin-1);HOK->Add(Hi);printf("  Add [%d,%d]",strt,bin-1);strt=-1;}};if (ok) {sprintf(HiN,"h%d_%d",strt,64);if((Hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();Hi=H->ProjectionX(HiN,strt,64);HOK->Add(Hi);printf("  Add [%d,%d]\n",strt,64);} else printf("\n");
	// In order to condition both p+ and pi-, simply change the bounds of the 2nd for-loop.
	// In order to have an OR instead of an AND, initialize "id=0x3" and reset it w/ "id^=1<<ppi".
	*/

	TLorentzVector *lvLambda = 0,    // ***** VECTORS w/ *EXACT* Lambda MASS
	  *lvaLambda = 0;

	//                        ***** HISTO (p,pi) SYSTEM KINEMATICS AS A F(m)
	if (winLRange&0x1) {
	  double EL = sqrt(M2_Lam+lvppi.Vect().Mag2());
	  lvLambda = new TLorentzVector(lvppi.Vect(),EL);
	}
	if (winLRange&0x2) {
	  double EL = sqrt(M2_Lam+lvpip.Vect().Mag2());
	  lvaLambda = new TLorentzVector(lvpip.Vect(),EL);
	}
	double xF[2]; for (iLaL = 0; iLaL<2; iLaL++) {      // ***** FEYNMAN's X
	  if (!(winLRange&0x1<<iLaL)) continue;
	  TLorentzVector lvL = iLaL==0 ? *lvLambda : *lvaLambda;
	  lvL.Boost(betacm);
	  //xF[iLaL] = 2*lvL.Vect().Dot(vqUnit)/sqrts; // Approximation
	  // Maximal momentum is obtained w/ gamm*N -> Lambda K
	  // sqrt(s) = EL+EK = sqrt(mL^2+p^2)+sqrt(mK^2+p^2)
	  // A-2p^2 = 2sqrt(mL^2+p^2)(mK^2+p^2), A=s-mL^2-mK^2
	  // P^2 = A^2-4mK^2mL^2 / 4 s
	  double A = sqrts*sqrts-M2_Lam-M2_K, B = A*A-4*M2_Lam*M2_K;
	  if (B<=0) {
	    //#define U3_DEBUG_LambdaK
#ifdef U3_DEBUG_LambdaK
	    printf(" * U3:\a y(=%.2f=>sqrt(s)=%.1f) too small for Lambda+K\n",
		   yB,sqrts);
#endif
	    xF[iLaL] = -1;
	  }
	  else {
	    double pZMx = sqrt(B)/2/sqrts;
	    xF[iLaL] = lvL.Vect().Dot(vqUnit)/pZMx;
	  }
	}

	double cthstar;                  // ***** cos(theta*) MASS DISTRIBUTIONS
	if (winLRange&0x1) {
	  TLorentzVector lvp1_L = lvp1, lvq_L = lvq;
	  lvp1_L.Boost(-lvLambda->BoostVector());   // Boost p to L frame
	  lvq_L.Boost(-lvLambda->BoostVector());    // Boost gamma* to L frame
	  cthstar = cos(lvp1_L.Angle(lvq_L.Vect()));
	  if (papID&0x1) {
	    hm_LcthS[0]->Fill(m_ppi,cthstar);
	    int ixF, bin; for (bin = 1, ixF = 0; bin<5; bin++) {
	      if (xFbins[bin-1]<xF[0] && xF[0]<xFbins[bin]) {
		ixF = bin; break;
	      }
	    }
	    if (ixF) hm_LcthS[ixF]->Fill(m_ppi,cthstar);
	  }
	}
	if (winLRange&0x2) {
	  TLorentzVector lvp2_aL = lvp2, lvq_aL = lvq;
	  lvp2_aL.Boost(-lvaLambda->BoostVector());// Boost bar{p} to aL frame
	  lvq_aL.Boost(-lvaLambda->BoostVector()); // Boost gamma* to aL frame
	  cthstar = cos(lvp2_aL.Angle(lvq_aL.Vect()));
	  if (papID&0x2) {
	    hm_aLcthS[0]->Fill(m_pip,-cthstar);
	    int ixF, bin; for (bin = 1, ixF = 0; bin<5; bin++) {
	      if (xFbins[bin-1]<xF[1] && xF[1]<xFbins[bin]) {
		ixF = bin; break;
	      }
	    }
	    if (ixF) hm_aLcthS[ixF]->Fill(m_pip,cthstar);
	  }
	}

	//   ********** CANDIDATE Lambda SELECTION **********
	int isLambda = winLRange;
	if (fabs(m_ppi-M_Lam)>LMassCut) isLambda &= 0x2;
	if (fabs(m_pip-M_Lam)>LMassCut) isLambda &= 0x1;
	if (isLambda) {
	  if (isLambda&0x1) {               // ***** L DECAY MOMENTA w/o RICH ID
	    hk_PpL[0]->Fill(mom1); hk_PpiL[0]->Fill(mom2);
	  }
	  if (isLambda&0x2) {
	    hk_PpL[1]->Fill(mom2); hk_PpiL[1]->Fill(mom1);
	  }
	  //                                     ***** (X,Y), THETA vs. P @ RICH
	  // (This is to visualize the Lambda kinematical domain for RICH perfs:
	  // there, when targeting +/-, we allow ourselves to ID -/+.)
	  for (iLaL = 0; iLaL<2; iLaL++) {
	    if (!(isLambda&0x1<<iLaL)) continue;
	    const PaTrack *trk; double mom; const TLorentzVector *lvp;
	    for (pm = 0, iET = iET1, trk = &trk1, mom = mom1, lvp = &lvp1; pm<2;
		 pm++) {	      
	      if (idpOrms[iLaL]&1<<(1-pm)) {
		int ppi = pm==iLaL ? 0 : 1, mp = 1-2*pm;
		const vector<Float_t> &aux = trk->vAux();
		if (aRs[iET]<-1) { // Ensure then "aRs" are meaningful
		  printf("** U3:\a Evt %d,%d Track %d(<Lambda): Inconsistency: aRs = %d\n",
			 e.RunNum(),(int)e.UniqueEvNum(),iET,aRs[iET]);
		  abort();
		}
		float tgx = aux[2], tgy = aux[3];
		double thRICH = acos(1/sqrt(1.+ tgx*tgx + tgy*tgy));
		hk_PThL[ppi]->Fill(mom,thRICH);
		hk_XL[ppi]->Fill(aux[0],mp); hk_YL[ppi]->Fill(aux[1],mp);
		if (ppi==0 &&
		    s_dist>=3 && s_ctheta >=3 &&
		    PIDs[iET]&0x10) {
		  double thC = trk->RichInf(thCIndex);
		  // W/o and w/ index evolution.
		  hR_thCPL->Fill(mp*mom,thC);
		  double beta = lvp->Beta(), thCorr = thC+richDths[iET]/beta;
		  hR_ThCPL->Fill(mp*mom,thCorr);
		  double idx = DataTakingDB::Ptr()->MeanIndex;
		  double thCp = 1000*acos(1/beta/idx);
		  hR_THCPL->Fill(mp*mom,thCorr-thCp);
		}
	      }
	      iET = iET2; trk = &trk2; mom = mom2; lvp = &lvp2;
	    }
	  }
	} // End of "isLambda"
	  

	if (s_dist>=3) {                                // ***** VERY CLEAN V0'S
	  LambdaPat |= tightOK;
	  // ************************************************************
	  // ********** RICH EFFICIENCY (and related) from Lambda *******
	  // ************************************************************
	  // - No impact on the rest of the processing.
#  ifdef U3_FILL_RICH_PERFS
	  if (winLRange&0x1) {
	    bool isL = true;
	    fillLambda_RICHPerf(e,iET1,iET2,iET1,lvp1,m_ppi,isL,Run);
	    //if (!isL) isLambda &= 0x2;   // No to have any impact
	  }
	  if (winLRange&0x2) {
	    bool isaL = true;
	    fillLambda_RICHPerf(e,iET1,iET2,iET2,lvp1,m_pip,isaL,Run);
	    //if (!isaL) isLambda &= 0x1;  // No to have any impact
	  }
	}
#  endif

	if (e.IsMC()) {
	  int iLaLMC = 0;
	  if      (papID&0x1)
	    iLaLMC |= LambdaMC->Fill(e,trk1,trk2,0,isLambda&0x1);
	  else if (papID&0x2)
	    iLaLMC |= LambdaMC->Fill(e,trk2,trk1,1,isLambda&0x2);
	  //if      (iLaLMC&0x1) FillLambdaMC(m_ppi,0x1);
	  //else if (iLaLMC&0x2) FillLambdaMC(m_pip,0x2);
	}

	unsigned int isIDdLambda;
	for (iLaL = 0, isIDdLambda = 0; iLaL<2; iLaL++) {
	  if (!(1<<iLaL&isLambda)) continue;
	  if ((idpOrms[iLaL]&0x1<<iLaL) && (isNotees[iLaL]&0x1<<(1-iLaL)))
	    isIDdLambda |= 1<<iLaL;
	}
	for (iLaL = 0; iLaL<2; iLaL++) {
	  if (!(1<<iLaL&isIDdLambda)) continue;
	  // ************************************************************
	  // ********************  Lambda SELECTION  ********************
	  // - MASS CUT
	  // - (OPTIONALLY) pID
	  // ************************************************************
	  TLorentzVector *lvL = iLaL==0 ? lvLambda : lvaLambda;
	  hm_LID[iLaL]->Fill(iLaL==0 ? m_ppi : m_pip);
	  if      (iLaL==0) { // ***** DECAY MOMENTA and KINEMATICS w/ RICH ID
	    hk_PpL[2]->Fill(mom1); hk_PpiL[2]->Fill(mom2);
	  }
	  else {
	    hk_PpL[3]->Fill(mom2); hk_PpiL[3]->Fill(mom1);
	  }
	  //                     ***** EVENT KINEMATICS 
	  if (primaryPat&0x10) { // Is this useful??
	    hQ2_L[iLaL]->Fill(Q2); hxB_L[iLaL]->Fill(xB);
	    hyB_L[iLaL]->Fill(yB);
	    hxF_L[iLaL]->Fill(xF[iLaL]); // Requires mu', because of sqrt(s)
	    hv_pZ_L[iLaL]->Fill(Zp);
	  }
	  //                     ***** Lambda+X: Sigma* (and Xi)
	  fillLambdaX(e,ipV,imuS,lvq,iET1,iET2,iLaL,lvL);
	}
	if (lvLambda) delete lvLambda; if (lvaLambda) delete lvaLambda;
      } // End STRICT Vertex line cuts

#  ifdef U3_OUTPUT_TREE
      //                                     ********** LOOSE V0 SELECTION AGAIN
      //                                                  ***** FILL OUTPUT TREE
      if (s_pT>=2 && s_chi2>=2 && s_dist>=1 && s_ctheta>=1) {
	for (iLaL = 0; iLaL<2; iLaL++) if (winLRange&1<<iLaL) {
	    double m = iLaL==0 ? m_ppi : m_pip;
	    copyCSResonanceData(sV,m,dist,ddist,ctheta,pT,alpha);
	    CSResonanceData &r = fResonances.back();
	    r.LambdaPat =
	      LambdaPat|RICHOKs[iLaL]<<4|idpOrms[iLaL]<<6|isNotees[iLaL]<<9;
	    if (iLaL) r.LambdaPat |= 0x800;
	    const PaParticle *pa; for (pm = 0, pa = pa1; pm<2; pm++) {
	      int fHInd = getHadronIndex(pa);
	      if (pm==0) r.h1 = fHInd; //adrons.size();
	      else       r.h2 = fHInd; //adrons.size();
	      if (fHInd==(int)fHadrons.size()) copyCSHadronData(e,pa,isV);
	      pa = pa2;
	    }
	  }
      }
#  endif

	//     *************************************************************
	//        ********** SEARCH FOR CASCADES -> Lambda pi/K **********
	//     *************************************************************

      if (s_dist>=3 && winLRange) {

	//   ********** STRICT DISTANCE, NO ANGLE CUT **********

	unsigned int isLambda = 0x3;
	if (fabs(m_ppi-M_Lam)>LMassCut) isLambda &= 0x2;
	if (fabs(m_pip-M_Lam)>LMassCut) isLambda &= 0x1;
	for (iLaL = 0; iLaL<2; iLaL++) if (1<<iLaL&isLambda) {
	    TLorentzVector *lvL = 0; // ***** VECTORS w/ *EXCACT* Lambda MASS
	    if (iLaL==0) {
	      double EL = sqrt(M2_Lam+lvppi.Vect().Mag2());
	      lvL = new TLorentzVector(lvppi.Vect(),EL);
	    }
	    else {
	      double EL = sqrt(M2_Lam+lvpip.Vect().Mag2());
	      lvL = new TLorentzVector(lvpip.Vect(),EL);
	    }
	    fillLCascades(e,ipV,CovP,imuS,isV,CovS,iET1,iET2,iLaL,lvL,ppiIDs[iLaL],
			  uDSTSelection);
	  }
      }  // End strict dist, loose angle
    }  // End of preliminary selection for Lambda
#endif
  } // End of loop over vertices for V0's

#ifdef DEBUG_K0
  if (dbK0Event>=0 && dbK0Event!=0xff) {
    // K0 was expected in current event but no 2 reco'd tracks found to match
    // it decay characteristics (else "dbK0Event" would have been set |= 0x3)
    fprintf(dbK0_o,"%5d %3d %5d %10Ld 0x%02x",
	    Run,e.SpillNum(),e.EvInSpill(),e.UniqueEvNum(),dbK0Event);
    if      (dbK0Event==0x3f || dbK0Event==0xf)
      fprintf(dbK0_o,"  %6.3f",dbK0v1);
    else if (dbK0Event==0x1f) fprintf(dbK0_o,"  %6.3f %6.3f",dbK0v1,dbK0v2);
    fprintf(dbK0_o,"  %ld\n",e.UniqueEvNum());
  }
#endif

  //
  // Event selection 
  //
  // if output stream is specified, current event will be saves
  // if TagToSave() function had been called at least once.

  if (uDSTSelection) {
#ifdef K0LambdaphiXi_WRITEBACK
    e.TagToSave();
#endif
#ifdef U3_OUTPUT_TREE
    // Get # of charged K+ (w/in RICH ID capabilities)
    fCSEvt->nKs = nChargedKsInPV(e,pV);
    fCSEvtTree->Fill();
#endif
  }

  delete PIDs; delete tZones; delete aRs; delete richDths;
  for (i = 0; i<5; i++) { delete richLHs[i]; delete richChi2s[i]; }

}
//*********************************************************************
// *************************     parseBestPV    ***********************
//*********************************************************************
void parseBestPV(PaEvent &e, int ipV,
		 unsigned short &primaryPat,
		 double &Xp, double &Yp, double &Zp, TMatrixD &CovP,
		 TLorentzVector &lvp, double &E0,
		 int &imuS, const PaParticle *&muS,
		 TLorentzVector &lvq,
		 double &yB, double &Q2, double &xB)
{
  DataTakingDB *DTDB = DataTakingDB::Ptr();

  // ********** RETRIEVE PRESELECTED BEST PRIMARY VERTEX **********
  const PaVertex pV = e.vVertex(ipV);
  primaryPat |= 0x1;

  // pV coordinates and cov. matrix
  Xp = pV.Pos(0); Yp = pV.Pos(1); Zp = pV.Pos(2);
  int i, j, k; for (i=k = 0; i<3; i++) for (j = 0; j<=i; j++, k++) {
      CovP(i,j) = pV.Cov(k); if (i!=j) CovP(j,i) = pV.Cov(k);
    }

  const PaParticle &muI = e.vParticle(pV.InParticle());
  if (muI.NFitPar()==0) {
    printf("\n** U3:\a MuI Inconsistency: Run %d Evt %d nPar %d Q %d\n\n",
	   e.RunNum(),(int)e.UniqueEvNum(),muI.NFitPar(),muI.Q());
    assert(false);
  }

  if (muonBeam) {                                        // ***** BEAM SELECTION
    TVector3 v3I = muI.ParInVtx(ipV).Mom3();
    double pI = v3I.Mag(), pI2 = pI*pI;
    const PaTPar &hiI = e.vTrack(muI.iTrack()).vTPar()[0];
    double dqP2BMS = hiI(5,5);
    if (dqP2BMS<dqP2BMSCut) { // Genuine BMS
#ifdef U3_GENERALPURPOSE_HISTOS
      h_pBMS->Fill(pI);
#endif
      if (fabs(pI-pBMS)<dpBMS) primaryPat |= 0x2;
    }
#ifdef U3_GENERALPURPOSE_HISTOS
    h_dqPBMS->Fill(sqrt(dqP2BMS));
#endif
    if (!DTDB->BadBackProp(muI)) primaryPat |= 0x4;
    E0 = sqrt(pI2+M2_mu); TLorentzVector lvki(v3I,E0);

    int nTrksPV = pV.NOutParticles();                  // ***** GET SCATTERED MU
    for (int ip = 0; ip<nTrksPV; ip++) { 
      int iEP = pV.iOutParticle(ip); const PaParticle &pa = e.vParticle(iEP);
      int iET = pa.iTrack(); const PaTrack &trk = e.vTrack()[iET];
      if (trk.XX0()>15) { primaryPat |= 0x8; break; }
    }

#if MuPrim_OPTION == 1
    // ***** MU' OPTIONS: Update default (see infra) w/ MuPrim_MASK (if defined)
    // Setting of U3 default mask
    static bool checkYokeSM2; // Yoke is under control in latest mass productions
    static bool reject2MuEvents = false; // Imho, selecting highest P is better strategy
    static bool checkCanBeMuon = true, checkHodos = true; // Konrad's defaults
    static bool resetMuSMask = true; if (resetMuSMask) {
      resetMuSMask = false;
      checkYokeSM2 = DTDB->CheckYokeSM2;
      bool *muSMask[4] =
	{ &checkYokeSM2, &reject2MuEvents, &checkCanBeMuon, &checkHodos};
      int ibit; unsigned int defaultMask, mask;
      for (ibit = 0, defaultMask=mask = 0; ibit<4; ibit++) {
	unsigned int bit = 1<<ibit;
	if (*muSMask[ibit]) defaultMask |= bit;
#  ifdef MuPrim_MASK
	if (bit&MuPrim_MASK) *muSMask[ibit] = true;
	else                 *muSMask[ibit] = false;
#  endif
	if (*muSMask[ibit]) mask |= bit;
      }
      printf("\n * U3:\a Mask in Konrad's mu'-ID = 0x%x",mask);
      if (mask!=defaultMask)
	printf(" (instead of default 0x%x)\n",defaultMask);
      else
	printf("\n");
    }
    imuS = pV.iMuPrim(checkYokeSM2,reject2MuEvents,checkCanBeMuon,checkHodos);
#elif MuPrim_OPTION == 2                   // ***** ALTERNATIVE MU'-ID = CORAL'S
    imuS = pV.iMuPrimCoral();
#endif
    if (imuS>=0) {
      muS = &(e.vParticle(imuS));

      //      ********** RETRIEVE SCATTERED MU => KINEMATICS **********
      if (muS->NFitPar()==0) {
	printf("\n** U3:\a Run %d Evt %d Leading particle w/o P\n",
	       e.RunNum(),(int)e.UniqueEvNum());
	assert(false);
      }
      primaryPat |= 0x10;
      int iETs = muS->iTrack(); PaTrack &trks = e.vTrack()[iETs];
      if (trks.Chi2tot()/trks.Ndf()<chi2Cut)                   // ***** CHI2 CUT
	primaryPat |= 0x20;
      // Prior to, possibly re-scale momentum, determine exit angle for SM1 and
      // incidence/exit angle for SM2. May need to extrapolate, if so take
      // advantage to set also acceptance by RICH.
      int zoness = tZones[iETs];
      if (aRs[iETs]==0) // If not yet done
	transport(e,trks,aRs[iETs],zoness);
      TVector3 vS3 = muS->ParInVtx(ipV).Mom3();
      double Es = sqrt(vS3.Mag2()+M2_mu); lvq = lvki-TLorentzVector(vS3,Es);
      double pq = lvp.Dot(lvq), pk = lvp.Dot(lvki);

      //                                                   ***** MUON KINEMATICS
      Q2 = -lvq.M2(); xB = Q2/(2*pq); yB = pq/pk;
      if (Q2>Q2Cut) primaryPat |= 0x200;                         // ***** Q2 CUT
      if (yLowCut<yB && yB<yUpCut) primaryPat |= 0x100;          // ****** y CUT
    }
  }

  if (DTDB->WinTargetZone(Xp,Yp,Zp)) {       	      // ***** TARGET ZONE CUT
    primaryPat |= 0x40;
    if (DTDB->WinTarget(Xp,Yp,Zp)==0x3)               // ***** TARGET PROPER CUT
      primaryPat |= 0x80;
  }
}
//*********************************************************************
// *************************   bookKinematics   ***********************
//*********************************************************************
void bookKinematics(double phiMassCut, double K0MassCut,
		    const double *dECuts, int Ephi_pTCut,
		    const double *V0DsdD, const double *V0cthCut, const double *V0pTCut)
{
  DataTakingDB *DTDB = DataTakingDB::Ptr();

  char hT[] =
    "Q2 - K0#pm10Mev,D>2.0#deltaD,cos#theta>0.99990,pT>100MeV,#pi+|#pi-ID|eVeto  ";
  //"#SigmapT2>2.5 - #chi3<10,0.15<y<0.85,Target#4,Q2>0.02;ZpV (cm)  ";
  //"Q2 - #Lambda#pm20Mev - #chi3<10,0.15<y<0.85,Target#4,Q2>0.02  ";
  //"#mu' q/P  ";
  //"Q2 - #phi#pm10Mev,|dE|<10.2GeV,pT>100MeV,K+|K-ID  ";
  size_t len = strlen(hT);
  char cCCut[] = "#chi3<10"; size_t size = strlen(cCCut)+1;
  if (snprintf(cCCut,size,"#chi2<%.0f",chi2Cut)>=(int)size) {
    printf("\n** bookKinematics:\a Error writing string \"%s...\"\n\n",cCCut);
    assert(false);
  }
  string select(cCCut);
#ifdef U3_TARGET_CUT
  char cTCut[] = ",Target#4"; sprintf(cTCut,",Target#%d",U3_TARGET_CUT);
  select += string(cTCut);
#endif
#ifdef U3_y_CUT
  char cYCut[] = ",0.15<y<0.85"; size = strlen(cYCut)+1;
  if (snprintf(cYCut,size,",%4.2f<y<%4.2f",yLowCut,yUpCut)>=(int)size) {
    printf("\n** bookKinematics:\a Error writing string \"%s...\"\n\n",cYCut);
    assert(false);
  }
  select += string(cYCut);
#endif
#if U3_Q2_CUT > 1
  char cQ2Cut[] = ",Q2>0.02"; size = strlen(cQ2Cut)+1;
  if (snprintf(cQ2Cut,size,",Q2>%4.2f",Q2Cut)>=(int)size) {
    printf("\n** bookKinematics:\a Error writing string \"%s...\"\n\n",cQ2Cut);
    assert(false);
    select += string(cQ2Cut);
  }
#endif  
  int pTCut = int(V0pTCut[0]*1000+.5);

#ifdef U3_ALLOUT_KINE
  gDirectory->cd("/");
  TDirectory *dKine; if (!(dKine = (TDirectory*)gDirectory->Get("Kine")))
    gDirectory->mkdir("Kine","Kinematics of input events");
  if (!(dKine = (TDirectory*)gDirectory->Get("Kine"))) {
    printf("\n** U3:\a No creating subdir \"Kine\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dKine->cd();

  // ***** VERTICES
  const PaSetup &setup = PaSetup::Ref();
  
  getFIMMAbscissae();
  double ZMn, ZMx; int nZbins; getTargetBinning(nZbins,ZMn,ZMx);
  snprintf(hT,len,"%s;ZpV (cm)",select.c_str());
  hv_pZ =   new TH1D("hv_pZ",  hT,nZbins,ZMn,ZMx);
  snprintf(hT,len,"%s;XpV (cm);YpV (cm)",select.c_str());
  hv_pXY =  new TH2D("hv_pXY", hT,100,-2.5,2.5,100,-2.5,2.5);
  snprintf(hT,len,"%s,w/in Target;ZpV (cm)",select.c_str());
  hv_pZT =  new TH1D("hv_pZT", hT,nZbins,ZMn,ZMx);
  snprintf(hT,len,"%s,w/inTarget;XpV (cm);YpV (cm)",select.c_str());
  hv_pXYT = new TH2D("hv_pXYT",hT,100,-2.5,2.5,100,-2.5,2.5);
  sprintf(hT,"%s;ZsV (cm)",select.c_str());
  hv_sZ  = new TH1D("hv_sZ",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"Zs V0 - D>%.0f#deltaD,cos#theta>%.4f,pT>%dMeV",
	  V0DsdD[2],V0cthCut[1],pTCut);
  hv_sZV0  = new TH1D("hv_sZV0",hT,nZbins,ZMn,ZMx);
  hv_KdRvsR = new TH2D("hv_KdRvsR","Beam: dR/dZ vs. R @ pV - K0",
		       50,0,5,100,-.005,.005);
#endif

  if (hadronPhysics) return;

  // ***** Q2 and x BINNING
  double sq16_10; sq16_10 = sqrt(sqrt(sqrt(sqrt(10.))));
  double Q2bin; int bin; double Q2bins[105];
  for (bin = 104, Q2bin = 100*sq16_10*sq16_10; bin>=0; bin--) {
    Q2bins[bin] = Q2bin; Q2bin /= sq16_10;
  }
  double xBbin; double xBbins[105];
  for (bin = 104, xBbin = 1*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10;
       bin>=0; bin--) {
    xBbins[bin] = xBbin; xBbin /= sq16_10;
  }
  
#ifdef U3_ALLOUT_KINE

  // ***** pVERTEX w/ pT, Q2
  sprintf(hT,"pT>.7,%s;ZpV (cm)",   select.c_str());
  hv_pZh1 =  new TH1D("hv_pZh1",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"#SigmapT2>2.5%s;ZpV (cm)",select.c_str());
  hv_pZh2 =  new TH1D("hv_pZh2",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"Q2>.1,%s;ZpV (cm)",   select.c_str());
  hv_pZQ2 =  new TH1D("hv_pZQ2",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"Q2>1,%s;ZpV (cm)",    select.c_str());
  hv_pZDI =  new TH1D("hv_pZDI",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"Q2>1,1h,%s", select.c_str());
  hv_pZQ2h = new TH1D("hv_pZQ2h;ZpV (cm)",hT,nZbins,ZMn,ZMx);
  sprintf(hT,"Q2>1,1h,%s;ZpV (cm)", select.c_str());
  hv_pZDIh = new TH1D("hv_pZDIh",hT,nZbins,ZMn,ZMx);

  double ZTMn = DTDB->ZTMn, ZTMx = DTDB->ZTMx;
  snprintf(hT,len,"%.0f<Z<%.0f,Q2>.1,%s;XpV (cm); YpV (cm)",   ZTMn,ZTMx,select.c_str());
  hv_XYQ2 =   new TH2D("hv_XYQ2", hT,120,-3,3,120,-3,3);
  snprintf(hT,len,"%.0f<Z<%.0f,Q2>1,%s;XpV (cm); YpV (cm)",    ZTMn,ZTMx,select.c_str());
  hv_XYDI =   new TH2D("hv_XYDI", hT,120,-3,3,120,-3,3);
  snprintf(hT,len,"%.0f<Z<%.0f,Q2>.1,1h,%s;XpV (cm); YpV (cm)",ZTMn,ZTMx,select.c_str());
  hv_XYQ2h =  new TH2D("hv_XYQ2h",hT,120,-3,3,120,-3,3);
  snprintf(hT,len,"%.0f<Z<%.0f,Q2>1,1h,%s;XpV (cm); YpV (cm)", ZTMn,ZTMx,select.c_str());
  hv_XYDIh =  new TH2D("hv_XYDIh",hT,120,-3,3,120,-3,3);

#ifdef U3_FIMP_VTCS
  if (Year==2016 || Year==2017) {
    // (X,Y) for slices in Z corresponding to FI03 and MP01UV planes
    char hN[] = "hv_FIXQ2h", coords[] = "XYUV"; double ZLw, ZUp;
    for (int co = 0; co<3; co++) {
      ZLw = ZFI03s[co]; ZUp = ZLw+.25; ZLw -= .25; char CO = coords[co];
      sprintf(hN,"hv_FI%cQ2",CO);
      snprintf(hT,len,"FI03%c - %.2f<Z<%.2f,Q2>.1,%s;XpV (cm); YpV (cm)",   CO,ZLw,ZUp,select.c_str());
      hv_FIQ2s[co] =  new TH2D(hN, hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_FI%cDI",CO);
      snprintf(hT,len,"FI03%c - %.2f<Z<%.2f,Q2>1,%s;XpV (cm); YpV (cm)",    CO,ZLw,ZUp,select.c_str());
      hv_FIDIs[co] =  new TH2D(hN,hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_FI%cQ2h",CO);
      snprintf(hT,len,"FI03%c - %.2f<Z<%.2f,Q2>.1,1h,%s;XpV (cm); YpV (cm)",CO,ZLw,ZUp,select.c_str());
      hv_FIQ2hs[co] = new TH2D(hN,hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_FI%cDIh",CO);
      snprintf(hT,len,"FI03%c - %.2f<Z<%.2f,Q2>1,1h,%s;XpV (cm); YpV (cm)", CO,ZLw,ZUp,select.c_str());
      hv_FIDIhs[co] = new TH2D(hN,hT,240,-3,3,240,-3,3);
    }
    for (int co = 0; co<2; co++) {
      ZLw = ZMP01UVs[co]; ZUp = ZLw+1.5; ZLw -= 1.5; char CO = coords[co+2];
      sprintf(hN,"hv_MP%cQ2",CO);
      snprintf(hT,len,"MP01%c - %.2f<Z<%.2f,Q2>.1,%s;XpV (cm); YpV (cm)",   CO,ZLw,ZUp,select.c_str());
      hv_MPQ2s[co] =  new TH2D(hN, hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_MP%cDI",CO);
      snprintf(hT,len,"MP01%c - %.2f<Z<%.2f,Q2>1,%s;XpV (cm); YpV (cm)",    CO,ZLw,ZUp,select.c_str());
      hv_MPDIs[co] =  new TH2D(hN,hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_MP%cQ2h",CO);
      snprintf(hT,len,"MP01%c - %.2f<Z<%.2f,Q2>.1,1h,%s;XpV (cm); YpV (cm)",CO,ZLw,ZUp,select.c_str());
      hv_MPQ2hs[co] = new TH2D(hN,hT,240,-3,3,240,-3,3);
      sprintf(hN,"hv_MP%cDIh",CO);
      snprintf(hT,len,"MP01%c - %.2f<Z<%.2f,Q2>1,1h,%s;XpV (cm); YpV (cm)", CO,ZLw,ZUp,select.c_str());
      hv_MPDIhs[co] = new TH2D(hN,hT,240,-3,3,240,-3,3);
    }
    // Debugging the Z abscissae of zone 0x1: we histogram independently the
    // material of FI03 reconstructed w/ different patterns of detectors: MPs
    // only w/ only little contribution from DCs, and vice versa.
    // For this, we use the Golden channel G = DIh
    int iw, jw, idet; for (jw = 0; jw<2; jw++) { mpWs[jw]=dcWs[jw] = -1; }
    int mpNWs, dcNWs; for (iw=idet=mpNWs=dcNWs = 0; iw<HIT_MAP_SIZE; iw++)
      for (int ib = 0; ib<32; ib++, idet++) {
	if (idet>=setup.NDetectors()) continue;
	const PaDetect &d = setup.Detector(idet); const string &name = d.Name();
	if      (name.find("MP")==0) {
	  int kw; for (jw = 0, kw = -1; jw<mpNWs; jw++) {
	    if (mpWs[jw]==iw) { kw = jw; break; }
	  }
	  if (kw<0) {
	    if (mpNWs>2) {
	      printf("** U3_FIMP_VTCS: Too many Words needed to describe MPs\n");
	      abort();
	    }
	    kw = mpNWs++; mpWs[kw] = iw; mpBs[kw] = 0;
	  }
	  mpBs[kw] |= 0x1<<ib;
	}
	else if (name.find("DC00")==0 || name.find("DC01")==0) {
	  int kw; for (jw = 0, kw = -1; jw<dcNWs; jw++) {
	    if (dcWs[jw]==iw) { kw = jw; break; }
	  }
	  if (kw<0) {
	    if (dcNWs>2) {
	      printf("** U3_FIDC_VTCS: Too many Words needed to describe DCs\n");
	      abort();
	    }
	    kw = dcNWs++; dcWs[kw] = iw; dcBs[kw] = 0;
	  }
	  dcBs[kw] |= 0x1<<ib;
	}
      }
    printf(   " * U3_FIMP_VTCS: MP bits:");
    for (int jw = 0; jw<3; jw++)
      if (mpWs[jw]>=0) printf(" %d=0x%x",mpWs[jw],mpBs[jw]);
    printf("\n * U3_FIMP_VTCS: DC bits:");
    for (int jw = 0; jw<3; jw++)
      if (dcWs[jw]>=0) printf(" %d=0x%x",dcWs[jw],dcBs[jw]);
    printf("\n");
    snprintf(hT,len,"W/ MPs - Q2>1,1h,%s", select.c_str());
    hv_pZGMP = new TH1D("hv_pZGMP",hT,nZbins,ZMn,ZMx);
    snprintf(hT,len,"W/ DCs - Q2>1,1h,%s", select.c_str());
    hv_pZGDC = new TH1D("hv_pZGDC",hT,nZbins,ZMn,ZMx);
  }
#endif

  // Track parameters for mu', vs. trigger
  snprintf(hT,len,"#mu' X - %s",  select.c_str());
  hX_T  = new TH2D("hX_T",  hT,200,-2.5,2.5,6,-.5,5.5);
  snprintf(hT,len,"#mu' Y - %s",  select.c_str());
  hY_T  = new TH2D("hY_T",  hT,200,-2.5,2.5,6,-.5,5.5);
  snprintf(hT,len,"#mu' X' - %s", select.c_str());
  hXp_T = new TH2D("hXp_T", hT,200,-.25,.25,6,-.5,5.5);
  snprintf(hT,len,"#mu' Y' - %s", select.c_str());
  hYp_T = new TH2D("hYp_T", hT,200,-.25,.25,6,-.5,5.5);
  snprintf(hT,len,"#mu' q/P - %s",select.c_str());
  hqoP_T= new TH2D("hqoP_T",hT,200, .0 ,.10,6,-.5,5.5);
  hPar_Ts[0] = hX_T;  hPar_Ts[1] = hY_T;
  hPar_Ts[2] = hXp_T; hPar_Ts[3] = hYp_T;
  hPar_Ts[4] = hqoP_T;

  snprintf(hT,len,"%s;Q2",select.c_str());
  hQ2    = new TH1D("hQ2",  hT,104,Q2bins);
  snprintf(hT,len,"pT>.7,%s;Q2",select.c_str());
  hQ2h1  = new TH1D("hQ2h1",hT,104,Q2bins);
  snprintf(hT,len,"#SigmapT2>2.5,%s;Q2",select.c_str());
  hQ2h2  = new TH1D("hQ2h2",hT,104,Q2bins);
  snprintf(hT,len,"%s;Q2;Trigger",select.c_str());
  hQ2_T  = new TH2D("hQ2_T",hT,104,Q2bins,   6,-.5,5.5);
  snprintf(hT,len,"%s;x",select.c_str());
  hxB    = new TH1D("hxB",  hT,104,xBbins);
  snprintf(hT,len,"%s;y",select.c_str());
  hyB    = new TH1D("hyB",  hT,110,-.05,1.05);
  snprintf(hT,len,"pT>.7,%s;x",select.c_str());
  hxBh1  = new TH1D("hxBh1",hT,104,xBbins);
  snprintf(hT,len,"pT>.7,%s;y",select.c_str());
  hyBh1  = new TH1D("hyBh1",hT,110,-.05,1.05);
  snprintf(hT,len,"#SigmapT2>2.5,%s;x",select.c_str());
  hxBh2  = new TH1D("hxBh2",hT,104,xBbins);
  snprintf(hT,len,"y  - #SigmapT2>2.5%s",select.c_str());
  hyBh2  = new TH1D("hyBh2",hT,110,-.05,1.05);
  snprintf(hT,len,"%s;x;Trigger",select.c_str());
  hxB_T  = new TH2D("hxB_T",hT,104,xBbins,   6,-.5,5.5);
  snprintf(hT,len,"%s;y;Trigger",select.c_str());
  hyB_T  = new TH2D("hyB_T",hT,110,-.05,1.05,6,-.5,5.5);
#endif

  //    ********** KINEMATICS FOR K0 OR phi EVENTS **********
  gDirectory->cd("/");
  TDirectory *dEphi; if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi")))
    gDirectory->mkdir("Ephi","Exclusive phi Histos");
  if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi"))) {
    printf("\n** U3:\a No creating subdir \"Ephi\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dEphi->cd();
  sprintf(hT,"Q2 - #phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,K+|K-ID",
	  phiMassCut*1000,dECuts[0],Ephi_pTCut);
  hQ2phi = new TH1D("hQ2phi",hT,104,Q2bins);
  sprintf(hT,"y  - #phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,K+|K-ID",
	  phiMassCut*1000,dECuts[0],Ephi_pTCut);
  hyBphi = new TH1D("hyBphi",hT,110,-.05,1.05);
  gDirectory->cd("/");
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
    gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();
  sprintf(hT,"Q2 - K0#pm%.0fMeV,D>%.0f#deltaD,cos#theta>%.4f,pT>%dMeV,#pi+|#pi-ID|eVeto",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hQ2K0  = new TH2D("hQ2K0", hT,104,Q2bins,   6,-.5,5.5);
  sprintf(hT,"y  - K0#pm%.0fMeV,D>%.0f#deltaD,cos#theta>%.4f,pT>%dMeV,#pi+|#pi-ID|eVeto",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hyBK0  = new TH2D("hyBK0", hT,110,-.05,1.05,6,-.5,5.5);
  // Lambda kinematics: cf. "bookLambdaMisc".
}
#ifdef U3_ALLOUT_KINE
//*********************************************************************
// *************************   fillAllOutKine   ***********************
//*********************************************************************
void fillAllOutKine(const PaEvent &e, int trigType,
		    int ipV, const PaVertex &pV, unsigned short primaryPat,
		    int imuS, const PaParticle *muS,
		    double Q2, double xB, double yB, const TLorentzVector &lvq)
{
  // ***** BASIC KINEMATICS HISTOS *****

  double Xp = pV.Pos(0), Yp = pV.Pos(1), Zp = pV.Pos(2);
  hv_pZ ->Fill(Zp); hv_pXY->Fill(Xp,Yp); // W/in at least target zone
  if (primaryPat&0x80) {                 // W/in target
    hv_pZT ->Fill(Zp); hv_pXYT->Fill(Xp,Yp);
  }
  // ***** Z of SECONDARY VERTICES
  for (int iv = 0; iv<e.NVertex(); iv++) { // Loop on all vertices
    const PaVertex &v = e.vVertex(iv); if (!v.IsPrimary())
      hv_sZ->Fill(v.Pos(2));    // ***** HISTOGRAM 2NDARY VERTICES *****
  }

  if (!(primaryPat&0x30)) return;                // ***** pV HAS, GOOD CHI2, MU'

  int nTrksPV = pV.NOutParticles();
  int nHighPTs, ip; float pT1, pT2;
  for (ip=nHighPTs = 0, pT1=pT2 = 0; ip<nTrksPV; ip++) { 

    // ********** SEARCH for HIGH pT (pT1/2>.7 SpT2>2.5) **********

    int iEP = pV.iOutParticle(ip);
    if (iEP==imuS) continue;                                // ***** EXCLUDE mu'
    const PaParticle &pa = e.vParticle(iEP);
    int iET = pa.iTrack(); const PaTrack &trk = e.vTrack()[iET];
    int zones = tZones[iET]; if (!(zones&0x2)) continue; // ***** EXCLUDE FRINGE
    if (trk.Chi2tot()/trk.Ndf()>chi2Cut)
      continue;                                         // ***** CUT ON CHI2/NDF
    const PaTPar &hv = pa.ParInVtx(ipV); TVector3 v3 = hv.Mom3();
    float pT = v3.Perp(lvq.Vect());
    if (pT<0) {
      printf("\n** U3:\a Run %d Evt %d Vtx %d Out %d: pT(=%f)<0!\n\n",
	     e.RunNum(),(int)e.UniqueEvNum(),ipV,iEP,pT);
      assert(false);
    }
    if (pT>.7) {
      nHighPTs++; if (pT>pT2) {
	if (pT>pT1) { pT2 = pT1; pT1 = pT; }
	else                     pT2 = pT;
      }
    }
  }
  int highPT = 0; if (nHighPTs>=1) {
    highPT = 0x1; if (nHighPTs>=2 && pT1*pT1+pT2*pT2>2.5) highPT = 0x3;
  }

  hQ2->Fill(Q2); hxB->Fill(xB); hyB->Fill(yB);
  if (Q2>.1) {
    double Xp = pV.Pos(0), Yp = pV.Pos(1);
    bool winTZ = DataTakingDB::Ptr()->WinTargetZone(Xp,Yp,Zp);
    hv_pZQ2->Fill(Zp); if (winTZ) hv_XYQ2->Fill(Xp,Yp); 
    if (nTrksPV>1) { hv_pZQ2h->Fill(Zp); if (winTZ) hv_XYQ2h->Fill(Xp,Yp); }
    if (Q2>1) {
      hv_pZDI->Fill(Zp); if (winTZ) hv_XYDI->Fill(Xp,Yp); 
      if (nTrksPV>1) {
	hv_pZDIh->Fill(Zp); if (winTZ) hv_XYDIh->Fill(Xp,Yp);
#ifdef U3_FIMP_VTCS
	int ip, nMPs, nDCs; for (ip = 0, nMPs=nDCs = 0; ip<nTrksPV; ip++) { 
	  int iEP = pV.iOutParticle(ip); const PaParticle &pa = e.vParticle(iEP);
	  int iET = pa.iTrack(); const PaTrack &trk = e.vTrack()[iET];
	  const unsigned int *found = trk.FoundHitsBitMap();
	  for (int jw = 0; jw<3; jw++) {
	    int iw = mpWs[jw]; unsigned int bs = mpBs[jw]&found[iw];
	    for (unsigned int b = 0x1; b<=0x10000000; b <<= 1) if (b&bs) nMPs++;
	  }
	  for (int jw = 0; jw<3; jw++) {
	    int iw = dcWs[jw]; unsigned int bs = dcBs[jw]&found[iw];
	    for (unsigned int b = 0x1; b<=0x10000000; b <<= 1) if (b&bs) nDCs++;
	  }
	}
	if      (nMPs && nDCs<3) hv_pZGMP->Fill(Zp);
	else if (nDCs && nMPs<3) hv_pZGDC->Fill(Zp);
#endif
      }
    }
#ifdef U3_FIMP_VTCS
    int co, cp; for (co = 0, cp = -1; co<3; co++) {
      double Z = ZFI03s[co]; if (Z-.25<Zp && Zp<Z+.25) { cp = co; break; }
    }
    if (cp>=0) {
      hv_FIQ2s[cp]->Fill(Xp,Yp); if (nTrksPV>1) hv_FIQ2hs[cp]->Fill(Xp,Yp);
      if (Q2>1) {
	hv_FIDIs[cp]->Fill(Xp,Yp); if (nTrksPV>1) hv_FIDIhs[cp]->Fill(Xp,Yp);
      }
      for (co = 0, cp = -1; co<2; co++) {
	double Z = ZMP01UVs[co]; if (Z-1.5<Zp && Zp<Z+1.5) { cp = co; break; }
      }
      if (cp>=0) {
	hv_MPQ2s[cp]->Fill(Xp,Yp); if (nTrksPV>1) hv_MPQ2hs[cp]->Fill(Xp,Yp);
	if (Q2>1) {
	  hv_MPDIs[cp]->Fill(Xp,Yp); if (nTrksPV>1) hv_MPDIhs[cp]->Fill(Xp,Yp);
	}
      }
    }
#endif
  }
  if (highPT) {
    hQ2h1->Fill(Q2); hxBh1->Fill(xB); hyBh1->Fill(yB);
    hv_pZh1->Fill(Zp);
    if (highPT&0x2) {
      hQ2h2->Fill(Q2); hxBh2->Fill(xB); hyBh2->Fill(yB);
      hv_pZh2->Fill(Zp);
    }
  }
  for (int i = 0; i<5; i++) hPar_Ts[i]->Fill(muS->ParInVtx(ipV)(i+1),trigType);
  hQ2_T->Fill(Q2,trigType); hxB_T->Fill(xB,trigType); hyB_T->Fill(yB,trigType);

}
#endif
//*********************************************************************
// *************************   bookExclusive    ***********************
//*********************************************************************
void bookExclusive(const double *dECuts, int pTCut, double dM_rho,
		   double dM_phi, double phiMassCut)
{
  gDirectory->cd("/");
  TDirectory *dEphi; if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi")))
    gDirectory->mkdir("Ephi","Exclusive phi Histos");
  if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi"))) {
    printf("\n** U3:\a No creating subdir \"Ephi\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dEphi->cd();

  // ******************** BOOKING for EXCLUSIVE phi ********************

  char hT[] =
    "#phi#pm%15Mev,|dE|<2.5GeV,pT>100MeV,K#mpID;P#pm (GeV);#ThetaC+d#ThetaC-#ThetaCK(P) (mrd)   ";
  //"dE - h+h-=K+K-,|dEK|<2.5GeV,pT>100MeV,#phi#pm15Mev,K+&K-ID;dE (GeV)   ";
  //"PK#pm vs. L|SAS - |dE|<2.5GeV,pT>10MeV,#phi#pm20Mev  ";
  //"|dE|<2.5GeV,pT>10MeV,#phi#pm20Mev;M#pi+#pi-"
  size_t sT = strlen(hT);
  //                                                    ***** EXCLUSIVITY HISTOS
  snprintf(hT,sT,"dE - pT>%dMeV",pTCut);
  hk_dE =   new TH1D("hk_dE",  hT,100,-5,10);
  snprintf(hT,sT,"dEK - h+h-=K+K-,pT>%dMeV",pTCut);
  hk_dEK =  new TH1D("hk_dEK", hT,100,-5,10);
  snprintf(hT,sT,"dE - h#pmh#pm,pT>%dMeV",pTCut);
  hk_dEW =  new TH1D("hk_dEW", hT,100,-5,10);

  //                                          ***** KINEMATICS OF phi CANDIDATES
  snprintf(hT,sT,"pT(h#pm) - |dEK|<%3.1fGeV,#phi#pm%.0fMev,p>%.2f#piThr,K+|K-ID",
	   dECuts[0],phiMassCut*1000,thrMargin);
  hk_pTphi =   new TH1D("hk_pTphi", hT,40,0,4);
  snprintf(hT,sT,"dEK - h+h-=K+K-,pT>%dMeV,#phi#pm%.0fMeV,p>%.2f#piThr,K+|K-ID",
	   pTCut,phiMassCut*1000,thrMargin);
  hk_dEKphi =  new TH1D("hk_dEKphi", hT,100,-5,10);
  snprintf(hT,sT,"dEKc - h+h-=K+K-,pT>%dMeV,#phi#pm%.0fMeV,p>%.2f#piThr,K+|K-ID",
	   pTCut,phiMassCut*1000,thrMargin);
  hk_dEKcphi = new TH1D("hk_dEKcphi",hT,100,-5,10);
  //                                             ***** (X,Y), theta vs. P @ RICH
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K#pmID;P (GeV);#thetaRICH (rd)",
	   phiMassCut*1000,dECuts[0],pTCut,thrMargin);
  hk_PThphi = new TH2D("hk_PThphi",hT,100,0,100,80,0,0.8);
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K#pmID;P (GeV);XRICH (cm)",
	   phiMassCut*1000,dECuts[0],pTCut,thrMargin);
  hk_Xphi = new TH2D("hk_Xphi",hT,160,-160,160,2,-2,2);
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K#pmID;P (GeV);YRICH (cm)",
	   phiMassCut*1000,dECuts[0],pTCut,thrMargin);
  hk_Yphi = new TH2D("hk_Yphi",hT,160,-160,160,2,-2,2);
  //                                                 ***** THETA CHERENKOV vs. P
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,K#mpID;P#pm (GeV);#ThetaC (mrd)",
	   phiMassCut*1000,dECuts[0],pTCut);
  hR_thCPphi  = new TH2D("hR_thCPphi",hT,320,-80,80,240,0,60);
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,K#mpID;P#pm (GeV);#ThetaC+d#ThetaC (mrd)",
	   phiMassCut*1000,dECuts[0],pTCut);
  hR_ThCPphi  = new TH2D("hR_ThCPphi",hT,320,-80,80,240,0,60);
  snprintf(hT,sT,"#phi#pm%.0fMev,|dE|<%3.1fGeV,pT>%dMeV,K#mpID;P#pm (GeV);#ThetaC+d#ThetaC-#ThetaCK(P) (mrd)",
	   phiMassCut*1000,dECuts[0],pTCut);
  hR_THCPphi  = new TH2D("hR_THCPphi",hT,320,-80,80,128,-4,4);

  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,No pT cut",dECuts[0]);
  hm_phiAll = new TH1D("hm_phiAll",hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV",dECuts[0],pTCut);
  hm_phi =    new TH1D("hm_phi",   hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,%s",dECuts[0],pTCut,"LAS");
  hm_phill  = new TH1D("hm_phill", hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,%s",dECuts[0],pTCut,"L+SAS");
  hm_phisl  = new TH1D("hm_phisl", hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,%s",dECuts[0],pTCut,"SAS");
  hm_phiss  = new TH1D("hm_phiss", hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,%s",dECuts[0],pTCut,"p<10GeV");
  hm_phiLL  = new TH1D("hm_phiLL", hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,%s",dECuts[0],pTCut,"p>10GeV");
  hm_phiHH  = new TH1D("hm_phiHH", hT,200,2*M_K-.02,M_phi+dM_phi);

  snprintf(hT,sT,"|dE|<%3.1fGeV;M#pi+#pi-",dECuts[0]);
  hm_rhoAll = new TH1D("hm_rhoAll",hT,200,2*M_pi-.02,M_rho+dM_rho);
  snprintf(hT,sT,"|dE|<%3.1fGeV,pT>%dMeV;M#pi+#pi-",dECuts[0],pTCut);
  hm_rho =    new TH1D("hm_rho",   hT,200,2*M_pi-.02,M_rho+dM_rho);

  char cID[] =
    "K+K- - KID:0x9=1.05#piThr<P<100GeV 0x12=P<KThr,Else<0.92Back 0x24=K>1.20Back,1.01Hadron";
  size_t size = strlen(cID)+1;
  if (snprintf(cID,size,
    "K+K- - KID:0x9=%.2f#piThr<P<%.0fGeV 0x12=P<KThr,Else<%.2fBack 0x24=K>%.2fBack,%.2fHadron",
	       thrMargin,PRICHCut,subKThrLHpiVeto,LHBckCut,LHCut)>=(int)size) {
    printf("\n** bookEphi:\a Error writing string \"%s...\"\n\n",cID);
    assert(false);
  }
  hm_phiVsID = new TH2D("hm_phiVsID",cID,200,2*M_K-.02,M_phi+dM_phi,64,-0.5,63.5);
  snprintf(hT,sT,"K+K- - |dE|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K+|K-ID",
	  dECuts[0],pTCut,thrMargin);
  hm_phiID =   new TH1D("hm_phiID", hT,200,2*M_K-.02,M_phi+dM_phi);
  snprintf(hT,sT,"K#pmK#pm - |dE|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K+|K-ID",
	  dECuts[0],pTCut,thrMargin);
  hm_phiIDW =  new TH1D("hm_phiIDW",hT,200,2*M_K-.02,M_phi+dM_phi);

  // Diff in ELoss between pi and K corrected for:
  snprintf(hT,sT,"K+K- - |dEc|<%3.1fGeV,pT>%dMeV,p>%.2f#piThr,K+|K-ID,C",
	  dECuts[0],pTCut,thrMargin);
  hm_phiIDc =  new TH1D("hm_phiIDc",hT,200,2*M_K-.02,M_phi+dM_phi);

}
// **********************************************************************
// *************************   fillExclusive   **************************
// **********************************************************************
void checkExclusive(PaEvent &e, int ipV, int imuS,
		    const PaParticle *&pa1, const PaParticle *&pa2, const PaParticle *&pa3,
		    unsigned short &exclPhiPat);
bool detachedHadrons(PaEvent &e, const PaVertex &pV);
int  nInTimeNeutrals(PaEvent &e, const PaVertex &pV);
void fillExclusive(PaEvent &e, int ipV, int imuS,
		   const TLorentzVector &lvp, const TLorentzVector &lvq,
		   const double *dECuts, int pTCut /* Caveat: in MeV */,
		   double dM_rho, double dM_phi, double phiMassCut,
		   int &KEMiss, bool &uDSTSelection, bool &isphi)
{
  // Init before any early exit
  isphi = false; KEMiss = 0;
  
  // phi uDST Selection:
  // - On input: expecting pV, w/ mu, mu' and 2 h-tracks
  //   - Various cuts expected to be fulfilled as cpp defines.
  //   - Particular case of (incoming) mu: its BMS P should be genuine (as
  //    opposed to fixed value assigned by default) and w/in range (and this
  //    independent of "BMS_P_CUT"), since unprecise BMS impacts directly (and
  //    adversely) on the exclusivity condition.
  // - Chi2 cuts on h-tracks.
  // - Require h-tracks upstream of RICH.
  // - Exclusivity assuming K mass |dEK| < dECuts[0]
  // - If "pTCut" finite (ifdef Ephi_pTCUT in calling program), pT of h-track
  //  w.r.t. h-system require > "pTCut".
  // phi Selection for histogramming, in addition:
  // - Opposite sign h-tracks.
  // Special case of U3_QUASIEXCL_PHI: allow for one extra slow particle. This
  // is to evaluate whether one gets a better phi signal compared to the
  // strictly exclusive case.
  // Returned values:
  // - KEMiss: =-1: dEK<-dECuts[0],
  //           =0:  Undef,
  //           =1: |dEK|<dECuts[0],
  //           =2: No BMS and dEK+dEK/dqP*dqP2BMS<dECuts[0], i.e. may be excl.
  //           =3: dEK>dECuts[0]
  // - uDSTSelection: Updated if excl. phi selected
  // - isphi: |m_KK-M_phi|<phiMassCut

  // ***** ALGORITHM:
  // - CHECK #TRACKS, ORDER THEM
  // - pT CUT
  // - REQUIRE MISSING ENERGY ~= NUCLEON MASS
  //   1) ASSUMING pi+pi-  => PLOT "dE" DISTRIBUTION
  //   2) ASSUMING K+K-    => GO ON TO EXCLUSIVE phi SELECTION
  // - Set CSEvt dEK.
  
  // ***** EXCLUSIVITY phi PATTERN
  unsigned short twoHs =  0x01;      // 0x1:  2 hadron tracks in pV and no more
  unsigned short pmHs =   0x02;      // 0x2:  h+h-
  unsigned short pTOK =   0x04;      // 0x4:  pT > cut
  unsigned short dEKOK =  0x08;      // 0x8:  dEK < cut
  //                                    0x30: h1/2 w/in RICH acc. and P range
  unsigned short hasBMS = 0x40;      // 0x40: BMS exists
  unsigned short BMSOK =  0x80;      // 0x80: BMS w/in range
  unsigned short detachedOK = 0x100; // 0x1000: No detached track compatible w/ vertex
  unsigned short neutralsOK = 0x200; // 0x2000: No, in time, activity in ECALS
  unsigned short exclPhiPat = 0, exclPhiRequired;
  exclPhiRequired = twoHs|pmHs; // Initialized. Later to be updated.
  unsigned short RICHOK = 0;

  //                                       ***** CHECK #TRACKS, ORDER THEM *****
  const PaVertex &pV = e.vVertex(ipV);
  if (pV.NOutParticles()>4) { // Double-check that calling program does enforce this condition
    printf("** U3::fillExclusive: Evt %d#%d Vtx %d w/ #OutParticles(=%d) > 4\n",
	   Run,EvNum,ipV,pV.NOutParticles()); abort();
  }
  const PaParticle *pa1 = 0, *pa2 = 0, *pa3 = 0; int pm;
  checkExclusive(e,ipV,imuS,pa1,pa2,pa3,exclPhiPat); // Assign twoHs=0x1, inter alia

  //                                                  ***** REJECTION CUTS *****
  int iET1 = pa1->iTrack(), iET2 = pa2->iTrack();
  PaTrack &trk1 = e.vTrack()[iET1], &trk2 = e.vTrack()[iET2];
  if (trk1.Chi2tot()/trk1.Ndf()>chi2Cut ||                     // ***** CHI2 CUT
      trk2.Chi2tot()/trk2.Ndf()>chi2Cut) return;
  int zones1 = tZones[iET1], zones2 = tZones[iET2];
  if (!(zones1&0x2) || !(zones2&0x2))              // ***** EXCLUDE FRINGE FIELD
    return; // (Note: Excluding also possible 0x11 tracks...)
  // Prior to, possibly re-scale momentum, determine acceptance by RICH and
  // exit angle for SM1 and incidence/exit angle for SM2.
  // (Cf. also comment supra, in V0 block, concenring Marcin's rescaling.)
  if (aRs[iET1]==0)/* Not yet done? */ transport(e,trk1,aRs[iET1],zones1);
  if (aRs[iET2]==0)                    transport(e,trk2,aRs[iET2],zones2);
  if (aRs[iET1]<-1 || aRs[iET2]<-1)
    return; // Skip if either track starting downstream of RICH
  const PaTPar &hv1 = pa1->ParInVtx(ipV), &hv2 = pa2->ParInVtx(ipV);
  TVector3 v13 = hv1.Mom3(), v23 = hv2.Mom3();
  // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
  const PaTPar &hi1 = trk1.vTPar()[0], &hi2 = trk2.vTPar()[0];
  double mom1 = fabs(1/hi1(5)), mom2 = fabs(1/hi2(5));
  double E1 = sqrt(v13.Mag2()+M2_pi),  E2 = sqrt(v23.Mag2()+M2_pi);
  TLorentzVector lv1 = TLorentzVector(v13,E1);  // >0 hadron @ vertex
  TLorentzVector lv2 = TLorentzVector(v23,E2);  // <0 hadron @ vertex
  TLorentzVector lv0 = lv1 + lv2;               // V0
  TLorentzVector lvPR =  lvp + lvq - lv0;   // Recoil proton assuming pi+pi-
  double dE = ( lvPR.M2() - M_p*M_p ) /2/M_p;   // Inelasticity
  double pT = lv1.Perp(lv0.Vect());
  if (pTCut) {                                                   // ***** pT CUT
    if (pT>pTCut*.001) exclPhiPat |= pTOK;
    exclPhiRequired |= pTOK;
  }
  double PL1 = sqrt(v13.Mag2()-pT*pT), PL2 = sqrt(v23.Mag2()-pT*pT);
  double alpha = (PL1-PL2)/(PL1+PL2);

  //        ********** K MASS ASSIGNMENT **********
  E1 = sqrt(v13.Mag2()+M2_K); E2 = sqrt(v23.Mag2()+M2_K);
  TLorentzVector lvK1 = TLorentzVector(v13,E1); // >0 hadron @ vertex
  TLorentzVector lvK2 = TLorentzVector(v23,E2); // <0 hadron @ vertex
  TLorentzVector lvKK = lvK1 + lvK2;            // V0
  //            ***** BMS: EVALUATE SOME LOW LIMIT OF EMiss
  // - BMS is necessary to get a meaningful selection of exclusive.
  // - Yet we also want here to single out potentially exclusive events lacking
  //  BMS so that they are later dircarded from the inclusive phi search. 
  const PaParticle &muI = e.vParticle(pV.InParticle());
  const PaTPar &hvI = muI.ParInVtx(ipV); double dqP2BMS = hvI(5,5);
  TVector3 v3I = hvI.Mom3(); double pI = v3I.Mag();
  if (dqP2BMS<dqP2BMSCut) {
    exclPhiPat |= hasBMS;
    if (fabs(pI-pBMS)<dpBMS) exclPhiPat |= BMSOK;
  }
  if (!(exclPhiPat&hasBMS) && 
      (exclPhiPat&twoHs)) { // Set aside the case w/ extra, 3rd, particle in pV
    double E0 = sqrt(pI*pI+M2_mu); TLorentzVector lvki(v3I,E0);
    double pJ = pI-dpBMS; v3I *= pJ/pI; E0 = sqrt(pJ*pJ+M2_mu);
    TLorentzVector lvqj(v3I,E0); lvqj += lvq-lvki;
    TLorentzVector lvKPR = lvp + lvqj - lvKK;  // Recoil proton
    double dEK = ( lvKPR.M2() - M_p*M_p )/2/M_p;  // Inelasticity
    if (dEK<dECuts[0]) KEMiss = 2;
    return;
  }
  if (!(exclPhiPat&hasBMS)) // Don't store dEK in CSEvt: not reliable...
    return;                 // ...exit rightaway
  TLorentzVector lvKPR = lvp + lvq - lvKK;  // Recoil proton assuming K+K-
  if (pa3) {
    const PaTPar &hv3 = pa3->ParInVtx(ipV); TVector3 v33 = hv3.Mom3();
    double E3 = sqrt(v33.Mag2()+M2_e); // Assuming 3rd particle is e+|e-
    lvKPR -= TLorentzVector(v33,E3);
  }
  double dEK = ( lvKPR.M2() - M_p*M_p )/2/M_p;  // Inelasticity
  hk_dEK->Fill(dEK);
 
  if      (dEK<-dECuts[0]) KEMiss = -1;
  else if (dEK> dECuts[0]) KEMiss = 3;
  else                     KEMiss = 1;
  if (dEK>10 || dE>10) return;  // ***** CRUDE dE CUT for uDST SELECTION *****

  // Now the, slightly, CPU extensive infra can take place
  //                                         ***** FURTHER EXCLUSIVITY ESTIMATES
  if (!detachedHadrons(e,pV)) exclPhiPat |= detachedOK;
  if (!nInTimeNeutrals(e,pV)) exclPhiPat |= neutralsOK;
  //                                                                   ***** PID
  unsigned short id = 0;     // 0x11: piTHr<P<PMax, 0x22: subKthr, 0x44: KID
  unsigned short idKpm = 0;  // ID below/above K threshold, 0x1: +, 0x2: -
  unsigned short idpipm = 0; // piID and !eID, 0x1: +, 0x2: -
  int iET; double mom; for (pm = 0, iET = iET1, mom = mom1; pm<2; pm++) { 
    if (aRs[iET]>0 && piThr<mom && mom<PRICHCut) {
      RICHOK |= 0x1<<pm;
      unsigned short idpm =
	thrMargin*piThr<mom ? 0x1 :// Above pi threshold w/ safety margin...
	// ...in order to avoid the region where pions generate few photons
	0;
      if (PIDs[iET]&0x10) {
	double piLH = richLHs[0][iET], KLH = richLHs[1][iET];
	double pLH  = richLHs[2][iET], eLH = richLHs[3][iET];
	double hadronLH = pLH>piLH ? pLH : piLH;
	double elseLH = eLH>hadronLH ? eLH : hadronLH;
	if (thrMargin*piThr<mom && mom<KThr &&
	    elseLH<subKThrLHpiVeto)                         idpm |= 0x2;
	if (KThr<mom && KLH>LHBckCut && KLH>hadronLH*LHCut) idpm |= 0x4;
	if (!(idpm&0x6) &&
	    eLH<piLHeVeto*piLH && piLH>pLH*LHCut)           idpipm |= 1<<pm;
      }
      else if (mom<KThr) // Which results in piTHr<P<KThr, w/o margin, this time
	// Absence of RICH block or LH =0: cf. comment in PID for Lambda
	/* */                                               idpm |= 0x2;
      if ((idpm&0x4) ||     // Actual KID
	  (idpm&0x2))       // SubKThreshold piVeto
	idKpm |= 1<<pm;
      id |= idpm<<(3*pm);
    }
    iET = iET2; mom = mom2;
  }

  double m_KK = lvKK.M();
  // ********** RHO: EXCLUSIVITY, HISTO dE and m_pipi, FILL TTree **********
  unsigned short exclPhiOK = exclPhiPat&exclPhiRequired;
  if ((exclPhiOK|pmHs)==exclPhiRequired) {
    bool isExcl = dE<dECuts[0]; double m_pipi = lv0.M();
    bool winMRange = fabs(m_KK-M_phi)>dM_phi && m_pipi<M_rho+dM_rho;
    bool rhoInvMass = // Here same sign pairs are included
      fabs(m_KK-M_phi)<dM_rho*2/3;
    if ((exclPhiPat&pmHs) && isExcl) hm_rhoAll->Fill(m_pipi);
    if ((exclPhiPat&pTOK) && isExcl) {
      // Note: We have in mind to evaluate RICH perf on these exclusive rho's.
      // => We need:
      // - Clean pi+pi-'s, since rho too wide for Gaussian fit selection.
      // - No RICH bias.
      // - RICH acceptance and P range for both +/-: future RICH perf. for one
      //  and clean present selection for the other.
      if ((id&0x11)==0x11) { // Both w/in RICH acceptance and P range
	if (idpipm) {        // One at least piID (piID proper || !eID)
	  if (exclPhiPat&pmHs) {
	    if (rhoInvMass)    hk_dE->Fill(dE);
	    if (isExcl)        hm_rho->Fill(m_pipi);
	  }
	  else if (rhoInvMass) hk_dEW->Fill(dE);
	}
	if (winMRange &&
	    isExcl) { // As opposed to phi case where exploring dEK in incl. vs. excl.
	  uDSTSelection = true;
#ifdef U3_OUTPUT_TREE
	  //                                              ***** FILL OUTPUT TREE
	  fCSEvt->dEK = dE; // Nota Bene: Misleading leaf name!...
	  copyCSResonanceData(pV,m_pipi,0,0,0,pT,alpha);
	  CSResonanceData &r = fResonances.back();
	  r.phiPat = exclPhiPat|RICHOK<<4|idpipm<<6;
	  r.phiPat |= 0x800;
	  const PaParticle *pa; for (pm = 0, pa = pa1; pm<2; pm++) {
	    int fHInd = getHadronIndex(pa);
	    if (pm==0) r.h1 = fHInd; //adrons.size();
	    else       r.h2 = fHInd; //adrons.size();
	    if (fHInd==(int)fHadrons.size()) copyCSHadronData(e,pa,ipV);
	    pa = pa2;
	  }
#endif
	}
      }
    }
  }
  if (dEK>10) return; // Remove dE<10 && dE>10
    
  //             ********** PHI **********
  isphi = (exclPhiPat&pmHs) && fabs(m_KK-M_phi)<phiMassCut;

  if (KEMiss==1) exclPhiPat |= dEKOK;              // ********** EXCLUSIVITY CUT
  exclPhiRequired |= dEKOK;

  // ***** HISTOGRAM DOCUMENTING THE IMPACT of SELECTION on phi SIGNAL *****
  exclPhiOK = exclPhiPat&exclPhiRequired;
  if ((exclPhiOK|dEKOK)==exclPhiRequired &&
      isphi && idKpm) hk_dEKphi->Fill(dEK);

  //                                                         ***** KK MASS PLOTS
  bool winMRange = fabs(m_KK-M_phi)<dM_phi;
  if ((exclPhiOK|pTOK)==(exclPhiRequired|pTOK) &&
      winMRange) {
    // This is to document the impact of the pT cut: let's disregard it
    hm_phiAll->Fill(m_KK);
    if (isphi && idKpm) {
      hk_pTphi->Fill(lv1.Perp(lvq.Vect()));
      hk_pTphi->Fill(lv2.Perp(lvq.Vect()));
    }
    if (exclPhiPat&pTOK) {
      hm_phi->Fill(m_KK); if (idKpm) hm_phiID->Fill(m_KK);
      if (isphi && idKpm) {
	//                                       ***** (X,Y), THETA vs. P @ RICH
	// (This is to visualize the phi kinematical domain for RICH perfs:
	// there, when targeting K+/-, we allow ourselves to ID K-/+.)
	const PaTrack *trk; double mom; const TLorentzVector *lvK;
	for (pm = 0, iET = iET1, trk = &trk1, mom = mom1, lvK = &lvK1; pm<2;
	     pm++) {
	  if (idKpm&1<<(1-pm)) {
	    const vector<Float_t> &aux = trk->vAux();
	    if (aRs[iET]<-1) { // Ensure then "aRs" are meaningful
	      printf("** U3:\a Evt %d,%d Track %d(<-phi): Inconsistency: aRs = %d\n",
		     e.RunNum(),(int)e.UniqueEvNum(),iET,aRs[iET]);
	      abort();
	    }
	    float tgx = aux[2], tgy = aux[3];
	    double thRICH = acos(1/sqrt(1.+ tgx*tgx + tgy*tgy));
	    hk_PThphi->Fill(mom,thRICH);
	    int mp = 1-2*pm; hk_Xphi->Fill(aux[0],mp); hk_Yphi->Fill(aux[1],mp);
	    //                                       ***** THETA CHERENKOV vs. P
	    if (PIDs[iET]&0x10) {
	      double thC = trk->RichInf(thCIndex);
	      // W/o and w/ index evolution.
	      hR_thCPphi->Fill(mp*mom,thC);
	      // Use reference index to have all events plotted in a same histo
	      // => Correct theta for the diff between current and reference one
	      double beta = lvK->Beta(), thCorr = thC+richDths[iET]/beta;
	      hR_ThCPphi->Fill(mp*mom,thCorr);
	      double idx = DataTakingDB::Ptr()->MeanIndex;
	      double thCK = 1000*acos(1/beta/idx);
	      hR_THCPphi->Fill(mp*mom,thCorr-thCK);
	    }
	  }
	  iET = iET2; trk = &trk2; mom = mom2; lvK = &lvK2;
	}
      }
      if (!(zones1&0x4)) {                               // ***** KK vs. LAS/SAS
	if (!(zones2&0x4)) hm_phill->Fill(m_KK);
	else               hm_phisl->Fill(m_KK);
      }
      else {
	if (!(zones2&0x4)) hm_phisl->Fill(m_KK);
	else               hm_phiss->Fill(m_KK);
      }
      if      (mom1<10 && mom2<10) hm_phiLL->Fill(m_KK);
      else if (mom1>10 && mom2>10) hm_phiHH->Fill(m_KK);

      hm_phiVsID->Fill(m_KK,id);                        // ***** KK MASS vs. PID
      /*
	// Examples of use of hm_phiVsID in ROOT:
	TH2D *H = (TH2D*)gDirectory->Get("hm_phiVsID");
	int nBinsX = H->GetNbinsX(); double xMn = H->GetBinLowEdge(1);
	double xMx = H->GetBinLowEdge(nBinsX)+H->GetBinWidth(nBinsX);
	unsigned int id, jd; int bin, strt, pm, ok;
	TH1D *Hi; char HiN[] = "h64_64";
	// "idOK", i.e. KID or elseVeto subKThr = id&0x4||id0x2
	TH1D *HOK = new TH1D("hm_phiOK","phi id&0x4||id&0x2",nBinsX,xMn,xMx);
	for(bin=1,strt=-1;bin<=64;bin++){for(pm=0,ok=1,id=bin-1;pm<2;pm++){id>>=3*pm;jd=id&0x7;if(!((jd&0x4)||(jd&0x2)))ok=0;};if(ok){printf("  %d 0x%x(0x%x,0x%x)",bin,bin-1,(bin-1)&0x7,(bin-1)>>3);if (strt<0)strt=bin;}else if (strt>=0){sprintf(HiN,"h%d_%d",strt,bin-1);if((Hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();Hi=H->ProjectionX(HiN,strt,bin-1);HOK->Add(Hi);printf("  Add [%d,%d]",strt,bin-1);strt=-1;}};if (ok) {sprintf(HiN,"h%d_%d",strt,64);if((Hi=(TH1D*)gDirectory->Get(HiN)))Hi->Delete();Hi=H->ProjectionX(HiN,strt,64);HOK->Add(Hi);printf("  Add [%d,%d]\n",strt,64);} else printf("\n");
	// In order to condition only K+ (or K-), simply change the bounds of the 2nd for-loop.
	// In order to have an OR instead of an AND, initialize "id=0x3" and reset it w/ "id^=1<<pm".
	*/

#ifdef U3_FILL_RICH_PERFS
      // ************************************************************
      // ********** RICH EFFICIENCY (and related) from phi **********
      // ************************************************************
      fillEphi_RICHPerf(e,iET1,iET2,lvK1,lvK2,m_KK,isphi);
#endif
    } // End cut on pT

  } // End cut on twoHs,h+h-,dEK
  else if ((exclPhiOK|pmHs)==exclPhiRequired &&
	   !(exclPhiPat&pmHs) &&
	   idKpm)
    hm_phiIDW->Fill(m_KK);

  if (winMRange &&                                       // ***** uDST SELECTION
      // Allow more than 2 hadrons and wrong sign pairs and extended dEK
      (exclPhiOK|twoHs|pmHs|dEKOK)==exclPhiRequired) {
    uDSTSelection = true;                // ********** uDST SELECTION: Excl. phi
#ifdef U3_OUTPUT_TREE
    //                                                    ***** FILL OUTPUT TREE
    fCSEvt->dEK = dEK;
    copyCSResonanceData(pV,m_KK,0,0,0,pT,alpha);
    CSResonanceData &r = fResonances.back();
    r.phiPat = exclPhiPat|RICHOK<<4|idKpm<<6;
    const PaParticle *pa; for (pm = 0, pa = pa1; pm<2; pm++) {
      int fHInd = getHadronIndex(pa);
      if (pm==0) r.h1 = fHInd; //adrons.size();
      else       r.h2 = fHInd; //adrons.size();
      if (fHInd==(int)fHadrons.size()) copyCSHadronData(e,pa,ipV);
      pa = pa2;
    }
#endif
  }

  isphi &=    // Condition "isphi" by strict selection, cf. "hm_phiID"
    exclPhiRequired && idKpm;

  if ((exclPhiOK|dEKOK)!=exclPhiRequired) return;

  // ***** E(pi+pi-) Loss CORRECTION *****
  const double phScale = 1.5, Zp = pV.Pos(2); PaTPar Hout;
  hi1.Extrapolate(Zp,Hout,true);
  v13 *= 1+phScale*(Hout.Mom()-mom1)/mom1;
  E1 = sqrt(v13.Mag2()+M2_K); lvK1 = TLorentzVector(v13,E1);
  hi2.Extrapolate(Zp,Hout,true);
  v23 *= 1+phScale*(Hout.Mom()-mom2)/mom2;
  E2 = sqrt(v23.Mag2()+M2_K); lvK2 = TLorentzVector(v23,E2);

  lvKK = lvK1 + lvK2;              // V0
  lvPR = lvp + lvq - lvKK; // Recoil proton
  double dEKc = ( lvPR.M2() - M_p*M_p )/2/M_p; // Inelasticity
  double m_KKc = lvKK.M(); bool isphic = fabs(m_KKc-M_phi)<phiMassCut;

  if (isphic && idKpm) hk_dEKcphi->Fill(dEKc);

  if (fabs(dEK)<dECuts[0] && // Same cut as supra: want to compare same phi evts
      idKpm && pTOK)
    hm_phiIDc->Fill(m_KKc);
}
//*********************************************************************
// *************************  bookV0Selection *************************
//*********************************************************************
void bookV0Selection()
{
  // ***** V0 SELECTION HISTOS: VV DISTANCE, COLLINEARITY, V-CHI2, pT *****
  // ***** K0
  gDirectory->cd("/");
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
		     gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();
  // SIGNAL
  hs_dK0  = new TH1D("hs_dK0","K0  D/#deltaD",    100,0,100);
  hs_aK0  = new TH1D("hs_aK0","K0  cos#theta",    100,.9998,1);
  hs_vK0  = new TH1D("hs_vK0","K0  #chi^{2}",     100,0.,10.);
  hs_pK0  = new TH1D("hs_pK0","K0  p_{T}",        100,0.,.25);
  hs_tK0  = new TH1D("hs_tK0","K0  Track Time",   100,0.,25.);
  // BACKGROUND
  hn_dK0  = new TH1D("hn_dK0","!K0  D/#deltaD",   100,0,100);
  hn_aK0  = new TH1D("hn_aK0","!K0  cos#theta",   100,.9998,1);
  hn_vK0  = new TH1D("hn_vK0","!K0  #chi^{2}",    100,0.,10.);
  hn_pK0  = new TH1D("hn_pK0","!K0  p_{T}",       100,0.,.25);
  hn_tK0  = new TH1D("hn_tK0","!K0  Track Time",  100,0.,25.);
#ifdef U3_K0_HIGHLEVEL  // K0 vs. Incidence and related
  hs_TK0  = new TH1D("hs_TK0","K0  #Theta",       80,0.,.4);
  hn_TK0  = new TH1D("hn_TK0","!K0  #Theta",      80,0,.4);
#endif
  gDirectory->cd("/");
  // ***** Lambda
  TDirectory *dLambda; if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda")))
			 gDirectory->mkdir("Lambda","Lambda Histos");
  if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda"))) {
    printf("\n** U3:\a No creating subdir \"Lambda\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dLambda->cd();
  // SIGNAL
  hs_dL   = new TH1D("hs_dL", "#Lambda  D/#deltaD",  100,0,100);
  hs_aL   = new TH1D("hs_aL", "#Lambda  cos#theta",  100,.9998,1);
  hs_vL   = new TH1D("hs_vL", "#Lambda  #chi^{2}",   100,0.,10.);
  hs_pL   = new TH1D("hs_pL", "#Lambda  p_{T}",      100,0.,.15);
  hs_tL   = new TH1D("hs_tL", "#Lambda  Track Time", 100,0.,25.);
  // BACKGROUND
  hn_dL   = new TH1D("hn_dL" ,"!#Lambda  D/#deltaD", 100,0,100);
  hn_aL   = new TH1D("hn_aL", "!#Lambda  cos#theta", 100,.9998,1);
  hn_vL   = new TH1D("hn_vL" ,"!#Lambda  #chi^{2}",  100,0.,10.);
  hn_pL   = new TH1D("hn_pL" ,"!#Lambda  p_{T}",     100,0.,.15);
  hn_tL   = new TH1D("hn_tL" ,"#!Lambda  Track Time",100,0.,25.);
  gDirectory->cd("/");
  // hs|n_aK0|L: The X axis lebelling, w/ labels such as "0.99995", is bound
  // to get crowded => let's make it somewhat sparser.
  hs_aK0->SetNdivisions(505); hn_aK0->SetNdivisions(505);
  hs_aL->SetNdivisions(505);  hn_aL->SetNdivisions(505);
}
// **********************************************************************
// ****************************** bookK0 ********************************
// **********************************************************************
void bookK0(const double *V0DsdD, const double *V0cthCut, const double *V0pTCut,
	    double dM_K0, double K0MassCut)
{
  gDirectory->cd("/");
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
    gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();

  // *************** K0 MASS CONDITIONED by Z of 2NDARY VERTEX ***************
  hm_K0a = new TH1D("hm_K0a", "#pi+#pi- - Z>35cm",               // Z>35
		    150,M_K0-dM_K0,M_K0+dM_K0);
  hm_K0  = new TH1D("hm_K0",  "#pi+#pi- - Z>50cm",               // Z>50
		    150,M_K0-dM_K0,M_K0+dM_K0);
  hm_K0t = new TH1D("hm_K0t", "#pi+#pi- - Z>50cm,TargetPointing",// Z>50+TargetPointing
		    150,M_K0-dM_K0,M_K0+dM_K0);
  hm_K0ll= new TH1D("hm_K0ll","#pi+#pi- - Z>50cm - LAS",         // LAS
		    150,M_K0-dM_K0,M_K0+dM_K0);
  hm_K0sl= new TH1D("hm_K0sl","#pi+#pi- - Z>50cm - L+SAS",
		    150,M_K0-dM_K0,M_K0+dM_K0);
  hm_K0ss= new TH1D("hm_K0ss","#pi+#pi- - Z>50cm - SAS",         // SAS
		    150,M_K0-dM_K0,M_K0+dM_K0);

  // *************** V0 MASS conditioned by VERTEX LINE ***************
  char V0Title[] = 
    "#pi+#pi- - D>2.0#deltaD,cos#theta>0.99990 - LAS#times(L)SAS vs. #Theta_{X}  ";
  //"#pi+#pi- - D>2.0#deltaD,cos#theta>0.99990 - #mu/#mu'  ";
  char K0Format[] = "#pi+#pi- - D>%.0f#deltaD,cos#theta>%.5f%s";
  //                             ***** CUT on D/dD and Angle: INCREASINGLY TIGHT
  sprintf(V0Title,K0Format,V0DsdD[0],V0cthCut[0]," ");
  hm_K0c11 = new TH1D("hm_K0c11",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[0],V0cthCut[1]," ");
  hm_K0c12 = new TH1D("hm_K0c12",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[0]," ");
  hm_K0c21 = new TH1D("hm_K0c21",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," ");
  hm_K0c22 = new TH1D("hm_K0c22",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - #mu/#mu'");
  hm_K0C22 = new TH1D("hm_K0C22",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  //                                                      ***** IN/OUT of TARGET
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - OUT");
  hm_K0o   = new TH1D("hm_K0o",  V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - IN");
  hm_K0i   = new TH1D("hm_K0i",  V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  //                                     *****  K0 as a function of spectrometer
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2}");
  hm_K0css = new TH1D("hm_K0css",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS#times(L)SAS");
  hm_K0csl = new TH1D("hm_K0csl",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2}");
  hm_K0cll = new TH1D("hm_K0cll",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  //                                                            ***** K0 vs. PID
   char cID[] =
     "#pi+#pi- - #piID:0x9=#piThr<P<60GeV 0x12=e<0.92#pi 0x24=#pi>1.20Back,1.01Hadron";
  size_t size = strlen(cID)+1;
  if (snprintf(cID,size,
     "#pi+#pi- - piID:0x9=#piThr<P<%.0fGeV 0x12=e<%.2f#pi 0x24=#pi>%.2fBack,%.2fHadron",
	       PRICHCut,piLHeVeto,LHBckCut,LHCut)>=(int)size) {
    printf("\n** bookK0:\a Error writing string \"%s...\"\n\n",cID);
    assert(false);
  }
  hm_K0VsID = new TH2D("hm_K0VsID",cID,150,M_K0-dM_K0,M_K0+dM_K0,64,-0.5,63.5);

#ifdef U3_K0_HIGHLEVEL
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - SAS^{2}");
  hm_K0cnn = new TH1D("hm_K0cnn",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - #pi>50GeV");
  hm_K0cSSp = new TH1D("hm_K0cSSp",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - #pi<-50GeV");
  hm_K0cSSm = new TH1D("hm_K0cSSm",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - #pi>60GeV");
  hm_K0cSSP = new TH1D("hm_K0cSSP",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - #pi<-60GeV");
  hm_K0cSSM = new TH1D("hm_K0cSSM",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - SAS^{2},#pi>50GeV");
  hm_K0cNNp = new TH1D("hm_K0cNNp",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - SAS^{2},#pi<-50GeV");
  hm_K0cNNm = new TH1D("hm_K0cNNm",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  //                                     ***** K0 as a function of opening angle
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. #Theta");
  hm_K0cllT = new TH2D("hm_K0cllT",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,3,-.5,2.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. #Theta");
  hm_K0cssT = new TH2D("hm_K0cssT",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,3,-.5,2.5);
  //                                ***** K0 as a function of opening angle in X
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. #Theta_{x}");
  hm_K0cllx = new TH2D("hm_K0cllx",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. #Theta_{x}");
  hm_K0cssx = new TH2D("hm_K0cssx",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS#timesLAS vs. #Theta_{x}");
  hm_K0cslx = new TH2D("hm_K0cslx",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS#times(L)SAS vs. #Theta_{x}");
  hm_K0clsx = new TH2D("hm_K0clsx",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  //                                ***** K0 as a function of opening angle in Y
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. #Theta_{y}");
  hm_K0clly = new TH2D("hm_K0clly",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. #Theta_{y}");
  hm_K0cssy = new TH2D("hm_K0cssy",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  //                                   ***** K0 as a function of incidence angle
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. incidence");
  hm_K0clla = new TH2D("hm_K0clla",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. incidence");
  hm_K0cssa = new TH2D("hm_K0cssa",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  //                                   ***** K0 as a function of exit opening
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. #Theta_{X}");
  hm_K0cllX = new TH2D("hm_K0cllX",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,6,-.5,5.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. #Theta_{X}");
  hm_K0cssX = new TH2D("hm_K0cssX",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,6,-.5,5.5);
  //                                   ***** K0 as a function of exit opening
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. mean exit");
  hm_K0cllb = new TH2D("hm_K0cllb",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,6,-.5,5.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. mean exit");
  hm_K0cssb = new TH2D("hm_K0cssb",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,6,-.5,5.5);
  //                                   ***** K0 as a function of K0 P
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - LAS^{2} vs. P");
  hm_K0cllP = new TH2D("hm_K0cllP",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,5,-.5,4.5);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - (L)SAS^{2} vs. P");
  hm_K0cssP = new TH2D("hm_K0cssP",V0Title,150,M_K0-dM_K0,M_K0+dM_K0,7,-.5,6.5);

  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - SAS/SAT");
  hm_K0cssS  = new TH1D("hm_K0cssS",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - SAS/LAT");
  hm_K0cssL  = new TH1D("hm_K0cssL",V0Title,150,M_K0-dM_K0,M_K0+dM_K0);
#endif

  //                                       ***** eVeto on pi <- K0 for K0X search
  sprintf(V0Title,K0Format,V0DsdD[1],V0cthCut[1]," - eVeto");
  hm_K0eV  = new TH1D("hm_K0eV", V0Title,150,M_K0-dM_K0,M_K0+dM_K0);

  //                ********** p+- (<-K0) momentum **********
  int pTCut = int(V0pTCut[0]*1000+.5);
  char K0PTitle[] = "K0#pm20MeV,D>2.0#deltaD,c#theta>0.99995,pT>100MeV,%X0<15,#pi#mpID|eVeto;P#pm (GeV);#ThetaC+d#ThetaC-#ThetaC#pi (mrd)  ";
  size_t sT = strlen(K0PTitle)+1;
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.4f,pT>%dMeV,#pi#pmID|eVeto;P (GeV);-/+",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hk_PK0 = new TH2D("hk_PK0",K0PTitle,100,0,100,2,-2,2);
  //            ********** (X,Y), THETA vs. P @ RICH **********
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.4f,pT>%dMeV,#pi#pmID|eVeto;P (GeV);XRICH (cm)",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hk_XK0 = new TH2D("hk_XK0",K0PTitle,160,-160,160,2,-2,2);
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.4f,pT>%dMeV,#pi#pmID|eVeto;P (GeV); YRICH (cm)",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hk_YK0 = new TH2D("hk_YK0",K0PTitle,160,-160,160,2,-2,2);
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.4f,pT>%dMeV,#pi#mpID|eVeto;P (GeV);#thetaRICH (rd)",
	  K0MassCut*1000,V0DsdD[1],V0cthCut[1],pTCut);
  hk_PThK0 = new TH2D("hk_PThK0",K0PTitle,100,0,100,80,0,0.8);
  //                                                 ***** THETA CHERENKOV vs. P
  // N.B.: Stricter cthCut and pTCut
  int pTStrictCut = int(V0pTCut[1]*1000+.5);
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.5f,pT>%dMeV,%%X0<15,#pi#mpID|eVeto;P (GeV);#ThetaC (mrd)",
	   K0MassCut*1000,V0DsdD[2],V0cthCut[2],pTStrictCut);
  hR_thCPK0   = new TH2D("hR_thCPK0" ,K0PTitle,320,-80,80,240,0,60);
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.5f,pT>%dMeV,%%X0<15,#pi#mpID|eVeto;P (GeV);#ThetaC+d#ThetaC (mrd)",
	   K0MassCut*1000,V0DsdD[2],V0cthCut[2],pTStrictCut);
  hR_ThCPK0   = new TH2D("hR_ThCPK0" ,K0PTitle,320,-80,80,240,0,60);
  snprintf(K0PTitle,sT,"K0#pm%.0fMeV,D>%.0f#deltaD,c#theta>%.5f,pT>%dMeV,%%X0<15,#pi#mpID|eVeto;P#pm (GeV);#ThetaC+d#ThetaC-#ThetaC#pi(P) (mrd)",
	   K0MassCut*1000,V0DsdD[2],V0cthCut[2],pTStrictCut);
  hR_THCPK0   = new TH2D("hR_THCPK0" ,K0PTitle,320,-80,80,128,-4,4);
}
//*******************************************************************
// *************************   fillK0woPV   *************************
//*******************************************************************
void fillK0woPV(double m_pipi, double Xs, double Ys, double Zs,
		double pT, double zL1, double zL2, const TLorentzVector &lvpipi)
{
  // ***** K0 W/O pV CONDITION
  if (pT>.02) {
    if (Zs>35) {                                          // Cut on Zs
      hm_K0a->Fill(m_pipi);
      if (Zs>50) {
	hm_K0->Fill(m_pipi);
	if (zL1<ZSM2) {                                   // In LAS, SAS...
	  if (zL2<ZSM2) hm_K0ll->Fill(m_pipi);
	  else          hm_K0sl->Fill(m_pipi);
	}
	else {
	  if (zL2<ZSM2) hm_K0sl->Fill(m_pipi);
	  else          hm_K0ss->Fill(m_pipi);
	}
	if (fabs(Xs-lvpipi.Vect()[0]/lvpipi.Vect()[2]*Zs)<5 &&
	    fabs(Ys-lvpipi.Vect()[1]/lvpipi.Vect()[2]*Zs)<5)
	  hm_K0t ->Fill(m_pipi);                          // Target Pointing
      }
    }
  }
}
#ifdef U3_FILL_RICH_PERFS
// *********************************************************************
// ************************* bookEphi_RICHPerf *************************
// *********************************************************************
void bookEphi_RICHPerf(double dM_phi, double phiMassCut)
{
  // ******************** EXCLUSIVE phi * RICH PERFS ********************
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dEphi; if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi")))
    gDirectory->mkdir("Ephi","Exclusive phi Histos");
  if (!(dEphi = (TDirectory*)gDirectory->Get("Ephi"))) {
    printf("\n** U3:\a No creating subdir \"Ephi\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dEphi->cd();
  TDirectory *dRICH; if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs")))
    gDirectory->mkdir("RICHPerfs","RICH Performances from Excl. phi decay");
  if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs"))) {
    printf("\n** U3:\a No creating subdir \"RICHPerfs\" in TDirectory \"%s/Ephi\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dRICH->cd();

#  ifdef U3_USE_CHCHI2
  char cChi2Cut[] = "K#chi^{2}<5.0,<0.95#pi,p  ";
  sprintf(cChi2Cut,"K#chi^{2}<%.1f,<#pi,p",CHChi2Mx); string sChi2KCut(cChi2Cut);
  // Cuts w/ pi veto
  sprintf(cChi2Cut,"K#chi^{2}<%.1f,<%.2f#pi,p",CHChi2Mx,CHChi2Cut);
  string sChi2Cut(cChi2Cut);
  sprintf(cChi2Cut,"K#chi^{2}<%.2f#pi,p",CHChi2Cut);  string sChi2VCut(cChi2Cut);
  sprintf(cChi2Cut,"#pi#chi^{2}>%.1f",subKThrChi2piVeto);
  string sChi2SCut(cChi2Cut);
#  endif
  char cLHCut[] = "KLH>1.005bck,>0.999#pi,p  ";
  int precLHBC = LHBckCut-1 ? int(-log10(fabs(LHBckCut-1))+.999) : 0;
  int precLHC  = LHCut-1    ? int(-log10(fabs(LHCut   -1))+.999) : 0;
  sprintf(cLHCut,"KLH>%#.*fbck,>#pi,p",precLHBC,LHBckCut);
  string sLHKCut(cLHCut);
  sprintf(cLHCut,"KLH>%#.*fbck,>%#.*f#pi,p",precLHBC,LHBckCut,precLHC,LHCut);
  string sLHCut(cLHCut);
  sprintf(cLHCut,"KLH>%#.*f#pi,p",precLHC,LHCut);
  string sLHVCut(cLHCut);
  int precLHSC = subKThrLHpiVeto-1 ? int(-log10(fabs(subKThrLHpiVeto-1))+.999)
    : 0;
  sprintf(cLHCut,"#piLH<%#.*fbck",precLHSC,subKThrLHpiVeto);
  string sLHSCut(cLHCut);

#  ifdef U3_PLOT_dTheta
  char cdTKCut[] = "#delta_{K}#theta<2.5#sigma  ";
  sprintf(cdTKCut,"#delta_{K}#theta<%.1f#sigma",nSigmas);
  char cpiVCut[] = "#delta_{#pi}#theta>2.5#sigma  ";
  sprintf(cpiVCut,"#delta_{#pi}#theta>%.1f#sigma",nVeto);
  string sdTKCut(cdTKCut), spiVCut(cpiVCut);
  string sdTCut(cdTKCut+string("#times")+cpiVCut);
#  endif
#  ifdef RICH_ANGLE_CUT
  char cRCut[] = ",A>25"; sprintf(cRCut,",R>%.0f",1000*aRICHCut);
#  else    // ***** Radius cut *****
  char cRCut[] = ",R>12"; sprintf(cRCut,",R>%.0f",rRICHCut);
#  endif
  string sRCut(cRCut);

  char cPRange[] = ",KThr<P<100GeV"; sprintf(cPRange,",KThr<P<%.0fGeV",PRICHCut);
  string KK("K+K- - "), sPRange(cPRange), sKID("K#pm:");
  string sPSub(",#piThr*1.5<P<KThr");

  //   ********** COMBINED +/- HISTOS FOR CORRELATED EFFICIENCIES **********
#  ifdef U3_PLOT_CHi2
  //             *************** Chi2 ***************
  // Raw efficiency
  hC_phi2   = new TH2D("hC_phi",  string(KK+sKID+sChi2Cut).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
  // Momenta over threshold
  hC_phiM2  = new TH2D("hC_phiM", string(KK+sKID+sChi2Cut+sPRange).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
  // Momenta over threshold AND radius(angle) < limit
  hC_phiMA2 = new TH2D("hC_phiMA",string(KK+sKID+sChi2Cut+sPRange+sRCut).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
#  else
  //              *************** LH ***************
  // Raw efficiency
  hL_phi2   = new TH2D("hL_phi",  string(KK+sKID+sLHCut).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
  // Momenta over threshold
  hL_phiM2  = new TH2D("hL_phiM", string(KK+sKID+sLHCut+sPRange).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
  // Momenta over threshold AND radius(angle) < limit
  hL_phiMA2 = new TH2D("hL_phiMA",string(KK+sKID+sLHCut+sPRange+sRCut).c_str(),
		      200,2*M_K-.02,M_phi+dM_phi,4,-.5,3.5);
#  endif

  // *************** Chi2, dTheta, LikeliHood ***************
  // (Successively K+ and K-)
  char cKpm[] = "K+:", hName[] = "hC_phiMp";
  for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
    sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);
    // ***** KID Eff FOR MOMENTA OVER THRESHOLD *****
#  ifdef U3_PLOT_CHi2
    // chi2
    sprintf(hName,"hC_phiM%c",pm?'m':'p');
    hC_phiM[pm] = new TH2D(hName,string(KK+spm+sChi2Cut+sPRange).c_str(),
		 	   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
    // SubThreshold chi2
    sprintf(hName,"hC_phiS%c",pm?'m':'p');
    hC_phiS[pm] = new TH2D(hName,string(KK+spm+sChi2SCut+sPSub).c_str(),
		 	   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
    // dTheta
    sprintf(hName,"hT_phiM%c",pm?'m':'p');
    hT_phiM[pm] = new TH2D(hName,string(KK+spm+sdTCut+sPRange).c_str(),
			   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
    // Likelihood
    sprintf(hName,"hL_phiM%c",pm?'m':'p');
    hL_phiM[pm] = new TH2D(hName,string(KK+spm+sLHCut+sPRange).c_str(),
			   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
    // SubThreshold Likelihood
    sprintf(hName,"hL_phiS%c",pm?'m':'p');
    hL_phiS[pm] = new TH2D(hName,string(KK+spm+sLHSCut+sPSub).c_str(),
			   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
  }

  for (int ipa = 0; ipa<nphB; ipa++) {

    // *************** KID EFF FOR P > THRES AS A F(PHASE SPACE) ***************

    char hN[] = "hC_phioPTp10";
    int PMn, PMx; char cABin[] = ",12<A<24";
    if (ipa==0) { PMn = 0;           PMx = 60; }
    else {
      //PMn = (ipa-1)%3*10; PMx = PMn+10;
      int ip = (ipa-1)%npB; PMn = pBins[ip]; PMx = pBins[ip+1];
    }
#  ifdef PS_USE_RADIUS
    char V = 'R'; double RBin = RRBin;
#  else
    char V = 'A'; double RBin = ARBin*1000;
#  endif
    if      (ipa==0)               sprintf(cABin,"%s","\0");
    else if (0<ipa && ipa<=npB)    sprintf(cABin,",%c<%.0f",V,RBin); 
    else if (npB*(nphiBins-1)<ipa) sprintf(cABin,",%.0f<%c",(ipa-1)/npB*RBin,V);
    else                           sprintf(cABin,",%.0f<%c<%.0f",
				       (ipa-1)/npB*RBin,V,((ipa-1)/npB+1)*RBin);
    char cPBin[] = ",15<P-KThr<60GeV ";
    if      (ipa==0)             sprintf(cPBin,",KThr<P<45GeV");
    else if ((ipa-1)%npB==npB-1) sprintf(cPBin,",KThr+%d<P<60GeV",PMn);
    else                         sprintf(cPBin,",%d<P-KThr<%dGeV",PMn,PMx);
    string sABin(cABin), sPBin(cPBin), sPABin(sPBin+sABin);

    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_USE_CHCHI2
      // ********** USING CHI2 CUTS **********
      sprintf(hN,"hC_phioP%c%d",pm?'m':'p',ipa);
      hC_phioP[pm][ipa]  = new TH2D(hN,string(KK+spm+sChi2KCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hC_phioPT%c%d",pm?'m':'p',ipa);
      hC_phioPT[pm][ipa] = new TH2D(hN,string(KK+spm+sChi2Cut +sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // pi enhanced rejection (pi Veto)
      sprintf(hN,"hC_phioPV%c%d",pm?'m':'p',ipa);
      hC_phioPV[pm][ipa] = new TH2D(hN,string(KK+spm+sChi2VCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      // ********** USING dTheta CUTS **********
      sprintf(hN,"hT_phioP%c%d",pm?'m':'p',ipa);
      hT_phioP[pm][ipa]  = new TH2D(hN,string(KK+spm+sdTKCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hT_phioPT%c%d",pm?'m':'p',ipa);
      hT_phioPT[pm][ipa] = new TH2D(hN,string(KK+spm+sdTCut +sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // pi Veto Only
      sprintf(hN,"hT_phioPV%c%d",pm?'m':'p',ipa);
      hT_phioPV[pm][ipa] = new TH2D(hN,string(KK+spm+spiVCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
      // Radius Distribution
      sprintf(hN,"hR_phioP%c%d",pm?'m':'p',ipa);
      hR_phioP[pm][ipa]  = new TH1D(hN,string(KK+spm+sPABin).c_str(),
				    80,0,80);
      // Momentum Distribution
      sprintf(hN,"hR_phioA%c%d",pm?'m':'p',ipa);
      hR_phioA[pm][ipa]  = new TH1D(hN,string(KK+spm+sPABin).c_str(),
				    80,0,40);
      // ********** USING LH CUTS **********
      sprintf(hN,"hL_phioP%c%d",pm?'m':'p',ipa);
      hL_phioP[pm][ipa]  = new TH2D(hN,string(KK+spm+sLHKCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hL_phioPT%c%d",pm?'m':'p',ipa);
      hL_phioPT[pm][ipa] = new TH2D(hN,string(KK+spm+sLHCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
      // pi enhanced rejection (pi Veto)
      sprintf(hN,"hL_phioPV%c%d",pm?'m':'p',ipa);
      hL_phioPV[pm][ipa] = new TH2D(hN,string(KK+spm+sLHVCut+sPABin).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
    }
  }
  for (int ipd = 0; ipd<4; ipd++) {

    // *************** KID EFF FOR P > THRES AS A F(pseudoPD) ***************

    char hN[] = "hT_phiPDp0", cPD[] = ",PD0";
    sprintf(cPD,",PD%d",ipd); string sPD(cPD);
    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_PLOT_CHi2
      sprintf(hN,"hC_phiPD%c%d",pm?'m':'p',ipd);
      hC_phiPD[pm][ipd] = new TH2D(hN,string(KK+spm+sChi2Cut+sPD).c_str(),
				   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      sprintf(hN,"hT_phiPD%c%d",pm?'m':'p',ipd);
      hT_phiPD[pm][ipd] = new TH2D(hN,string(KK+spm+sdTCut+sPD).c_str(),
				   200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
      sprintf(hN,"hL_phiPD%c%d",pm?'m':'p',ipd);
      hL_phiPD[pm][ipd]  = new TH2D(hN,string(KK+spm+sLHCut+sPD).c_str(),
				    200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
    }
  }
  for (int ih = 0; ih<4; ih++) {

    // *************** KID EFF FOR P > THRES AS A F(Time in Day) ***************

    char hN[] = "hT_phiHp0";
    char h1[] = "03:09", h2[] = "09:15", h3[] = "15:21", h4[] = "21:03";
    char *hours[4] = { h1, h2, h3, h4}; string sH(hours[ih]);
    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_PLOT_CHi2
      sprintf(hN,"hC_phiH%c%d",pm?'m':'p',ih);
      hC_phiH[pm][ih] = new TH2D(hN,string(KK+spm+sChi2Cut+sH).c_str(),
				 200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
#  endif
      sprintf(hN,"hL_phiH%c%d",pm?'m':'p',ih);
      hL_phiH[pm][ih] = new TH2D(hN,string(KK+spm+sLHCut+sH).c_str(),
				 200,2*M_K-.02,M_phi+dM_phi,2,-.5,1.5);
    }
  }

  for (int pm = 0; pm<2; pm++) {    // ********** LOOP ON K+ K- **********

    //     *************** CERENKOV ANGLE RESOLUTIONS... ***************
    //     ***** (N.B.: No ID of counterpart) *****

    char hName[] = "hR_dKvstp1";
    char hTitle[] = "RICH K+#delta#theta vs. #theta - #phi#pm20Mev,A<100  ";
    //char hTit[] = "#pi+ (X,Y) @ Z=750 - #phi#pm20Mev  ";  // XY
    char cpm = pm?'-':'+';
    for (int ipd = 0; ipd<4; ipd++) {          // ***** ...as a f(pseudoPD, +/-)
      sprintf(hName,"hR_dKPD%c%x",pm?'m':'p',ipd);
      sprintf(hTitle,"RICH K%c#delta#theta - #phi#pm%.0fMev,PD%d",
	      cpm,phiMassCut*1000,ipd);
      hR_dKPD[pm][ipd] = new TH1D(hName,hTitle,100,-10,10);
    }
    for (int ia = 0; ia<2; ia++) {  // ***** ...vs. Theta as f(+/-, SM1/2)
      sprintf(hName,"hR_dKvst%c%d",pm?'m':'p',ia);
      sprintf(hTitle,"RICH K%c#delta#theta vs. #theta - #phi#pm%.0fMev,%s",
	    cpm,phiMassCut*1000,ia?"Y>10":"A<10");
      hR_dKvst[pm][ia] = new TH2D(hName,hTitle,6,0,60,120,-6,6);
    }
  }
  for (int pm = 0; pm<2; pm++) {

    // *************** LIKELIHOOD PROFILES ***************

    // w/ RICH and P range
    char hName[] = "hR_LHKvsPp0";
    char hTitle[] = "K/#pi+LH - 15<P-KThrK<60GeV  ";
    //char hTitle[] = "#pi+LH - 1.5#piThr<P<KThr"
    sprintf(hName,"hR_LHK%c",pm?'m':'p');
    sprintf(hTitle,"K%cLH - KThr<P<%.0fGeV",pm?'-':'+',PRICHCut);
    hR_LHK[pm]   = new TProfile(hName,hTitle,60,2*M_K-.02,M_phi+dM_phi,"S");
    for (int ip = 0; ip<npB; ip++) {   // As a f(P)
      int PMn = pBins[ip], PMx = pBins[ip+1]; char cPBin[] = "15<P-KThr<60GeV";
      if   (ip==npB-1) sprintf(cPBin,"KThr+%d<P<60GeV",PMn);
      else             sprintf(cPBin,"%d<P-KThr<%dGeV",PMn,PMx);
      sprintf(hName,"hR_LHKvsP%c%d",pm?'m':'p',ip);
      sprintf(hTitle,"K%cLH - %s",pm?'-':'+',cPBin);
      hR_LHKvsP[pm][ip] =
	new TProfile(hName,hTitle,60,2*M_K-.02,M_phi+dM_phi,"S");
    }
    sprintf(hName,"hR_LHKpi%c",pm?'m':'p');
    sprintf(hTitle,"K/#pi%cLH - KThr<P<%.0fGeV",pm?'-':'+',PRICHCut);
    hR_LHKpi[pm] = new TProfile(hName,hTitle,60,2*M_K-.02,M_phi+dM_phi,"S");
    sprintf(hName,"hR_LHSpi%c",pm?'m':'p');
    sprintf(hTitle,"#pi%cLH - 1.5#piThr<P<KThr",pm?'-':'+');
    hR_LHSpi[pm] = new TProfile(hName,hTitle,60,2*M_K-.02,M_phi+dM_phi,"S");
  }

  for (int ipd = 0; ipd<4; ipd++) {

    // ********** K+- (<-phi) momentum as a f(pseudoPD) **********

    char hName[] = "hR_phiPK0";
    char hTitle[] = "PK#pm - #phi#pm20Mev,K#mp ID  ";
    sprintf(hName,"hR_phiP%d",ipd);
    sprintf(hTitle,"PK#pm - #phi#pm%.0fMev,w/in RICH",phiMassCut*1000);
    hR_phiP[ipd]  = new TH1D(hName,hTitle,200,-100,100);
    sprintf(hName,"hR_phiPK%d",ipd);
    sprintf(hTitle,"PK#pm - #phi#pm%.0fMev,K#mp ID",phiMassCut*1000);
    hR_phiPK[ipd] = new TH1D(hName,hTitle,200,-100,100);
    if (ipd==0) {

	// ********** Angular/Radius distribution **********

	sprintf(hTitle,"AK#pm - #phi#pm%.0fMev,K#mp ID",phiMassCut*1000);
	hR_phiAK = new TH1D("hR_phiAK",hTitle,150,0,150);
	sprintf(hTitle,"RK#pm - #phi#pm%.0fMev,K#mp ID",phiMassCut*1000);
	hR_phiRK = new TH1D("hR_phiRK",hTitle,150,0,150);
      }
    }

}
// *********************************************************************
// ************************* fillEphi_RICHPerf *************************
// *********************************************************************
void fillEphi_RICHPerf(PaEvent &e, int iET1, int iET2,
		       const TLorentzVector &lvK1, const TLorentzVector &lvK2,
		       double m_KK, bool isphi)
{
  int Run = e.RunNum();

  // ********** DETERMINE TRACKS' ACCEPTANCE by RICH **********

  const PaTrack &trk1 = e.vTrack()[iET1], &trk2 = e.vTrack()[iET2];
  // Start hlx
  const PaTPar &h1 = trk1.vTPar()[0], &h2 = trk2.vTPar()[0];
  double XR1, YR1, RR1, XR2, YR2, RR2;
  double *XRs[2] = {&XR1,&XR2}, *YRs[2] = {&YR1,&YR2}, *RRs[2] = {&RR1,&RR2};
  double AR1, AR2, *ARs[2] = {&AR1,&AR2};
  int acc = 0, pRange = 0, pDomain = 0; // Acceptance and p Range (+/-) Flags...
  int pm, iET; const PaTrack *trk; for (pm = 0, iET = iET1, trk = &trk1; pm<2;
					pm++) {
    if (aRs[iET]==0) {
 printf("** U3/Ephi RICHPerf:\a Evt %d,%d Track %d RICH not yet checked\n\n",
	e.RunNum(),(int)e.UniqueEvNum(),iET); assert(false);
    }
    if (aRs[iET]>0) acc |= 0x1<<pm;
    const vector<Float_t> &aux = trk->vAux();
    if (aRs[iET]<-1) { // Ensure then "aRs" are meaningful
      printf("** U3:\a Evt %d,%d Track %d(<-Ephi): Inconsistency: aRs = %d\n",
	     e.RunNum(),(int)e.UniqueEvNum(),iET,aRs[iET]);
      abort();
    }
    double XR = aux[0], YR = aux[1], RR = sqrt(XR*XR+YR*YR);
    *XRs[pm] = XR; *YRs[pm] = YR; *RRs[pm] = RR;
#  ifndef PS_USE_RADIUS
    float tgx = aux[2], tgy = aux[3];
    *ARs[pm] = acos(1/sqrt(1.+ tgx*tgx + tgy*tgy));
#  endif
    iET = iET2; trk = &trk2;
  }
  int subDomain = 0, subRange = 0;
  // More precision: which part of acceptance is hit...
  int ipd1  = 0, ipd2 = 0;  // .. PD = pseudo Photon Detector
  if (XR1>0) ipd1 |= 0x1; if (YR1>0) ipd1 |= 0x2;
  if (XR2>0) ipd2 |= 0x1; if (YR2>0) ipd2 |= 0x2;
  double mom1 = h1.Mom(), mom2 = h2.Mom();
  if (isphi) {                            // ***** KK == phi... *****
    //                    ...HISTOGRAM CHARACTERISTICS OF K's FROM phi
    hR_phiP[ipd1]->Fill(1/h1(5));                  // Momenta
    hR_phiP[ipd2]->Fill(1/h2(5));
  }

  if (acc) {
    // *********************************************************************
    // ******************** AT LEAST 1 HADRON W/IN RICH ********************
    // *********************************************************************

    if (piThr<mom1 && mom1<KThr) subDomain |= 0x1;
    if (piThr<mom2 && mom2<KThr) subDomain |= 0x2;
    // ********** SET A NUMBER OF RICH-RELATED VARIABLES **********
    int thid = 0;                // "th" = theta availability
    if (PIDs[iET1]&0x8) thid |= 0x1; if (PIDs[iET2]&0x8) thid |= 0x2;
    double thRich1 = 0, thRich2 = 0, thMom1 = 0, thMom2 = 0;
#  ifdef U3_PLOT_dTheta    // deltaTheta KID: Retain only 2 LS bits
    unsigned int dthok1 = rich->KID(trk1,thRich1,thMom1)&0x3;
    unsigned int dthok2 = rich->KID(trk2,thRich2,thMom2)&0x3;
#  else
    if (PIDs[iET1]&0x8) thRich1 = trk1.RichInf(RICH_ESTIMATOR);
    if (PIDs[iET2]&0x8) thRich2 = trk2.RichInf(RICH_ESTIMATOR);
    if (KThr+.1<mom1)
      thMom1 = acos(sqrt(M2_K/mom1/mom1+1)/rich->CurrentIndAPV)*1000;
    if (KThr+.1<mom2)
      thMom2 = acos(sqrt(M2_K/mom2/mom2+1)/rich->CurrentIndAPV)*1000;
#  endif
#  ifdef U3_USE_CHCHI2 
    // ***** KID: chi2 method *****
    int chiok1 = 0, chiok2 = 0, chiSubok1 = 0, chiSubok2 = 0;
    if (PIDs[iET1]&0x10) {
      double pichi2 = richChi2s[0][iET1], Kchi2  = richChi2s[1][iET1];
      double pchi2  = richChi2s[2][iET1];
      if (Kchi2<          pichi2 && Kchi2<          pchi2 &&
	  Kchi2<CHChi2Mx)                                  chiok1 |= 0x1;
      if (Kchi2<CHChi2Cut*pichi2 && Kchi2<CHChi2Cut*pchi2) chiok1 |= 0x2;
      if ((subDomain&0x1) && pichi2>subKThrChi2piVeto) {
#    ifdef REJECT_SubKThr_e
	bool eid = richChi2s[3][iET1]<subKThrChi2piVeto;
#    else
	bool eid = false;
#    endif
	if (!eid) chiSubok1 = 1;
      }
    }
    if (PIDs[iET2]&0x10) {
      double pichi2 = richChi2s[0][iET2], Kchi2  = richChi2s[1][iET2];
      double pchi2  = richChi2s[2][iET2];
      if (Kchi2<          pichi2 && Kchi2<          pchi2 &&
	  Kchi2<CHChi2Mx)                                  chiok2 |= 0x1;
      if (Kchi2<CHChi2Cut*pichi2 && Kchi2<CHChi2Cut*pchi2) chiok2 |= 0x2;
      if ((subDomain&0x2) && pichi2>subKThrChi2piVeto) {
#    ifdef REJECT_SubKThr_e
	bool eid = richChi2s[3][iET2]<subKThrChi2piVeto;
#    else
	bool eid = false;
#    endif
	if (!eid) chiSubok2 = 1;
      }
    }
#  endif
    // ***** KID: LH method *****
    int lhok1 = 0, lhok2 = 0, lhSubok1 = 0, lhSubok2 = 0;
    if (PIDs[iET1]&0x10) {
      double piLH = richLHs[0][iET1], KLH  = richLHs[1][iET1];
      double pLH  = richLHs[2][iET1];
      if (KLH>      piLH && KLH>      pLH && KLH>LHBckCut) lhok1 |= 0x1;
      if (KLH>LHCut*piLH && KLH>LHCut*pLH)                 lhok1 |= 0x2;
      if ((subDomain&0x1) && piLH<subKThrLHpiVeto)  {
#  ifdef REJECT_SubKThr_e
	bool eid = richLHs[3][iET1]>subKThrLHeVeto;
#  else
	bool eid = false;
#  endif
	if (!eid)                                          lhSubok1 = 1;
      }
    }
#  if K_BELOW_THR > 1
    else if (acc&0x1) {
      // SubThreshold KID: no tighter requirement on momentum than >piThr here,
      // but later on the "lhSubKok1" will be enabled only in a resticted range.
      if (subDomain&0x1)                                   lhSubok1 = 1;
      //else... in fact, no real need to ID pion below threshold any further
      //       than what's achieved by eVeto.
    }
#  endif
    if (PIDs[iET2]&0x10) {
      double piLH = richLHs[0][iET2], KLH  = richLHs[1][iET2];
      double pLH  = richLHs[2][iET2];
      if (KLH>      piLH && KLH>      pLH && KLH>LHBckCut) lhok2 |= 0x1;
      if (KLH>LHCut*piLH && KLH>LHCut*pLH)                 lhok2 |= 0x2;
      if ((subDomain&0x2) && piLH<subKThrLHpiVeto) {
#  ifdef REJECT_SubKThr_e
	bool eid = richLHs[3][iET2]>subKThrLHeVeto;
#  else
	bool eid = false;
#  endif
	if (!eid)                                          lhSubok2 = 1;
      }
    }	      
#  if K_BELOW_THR > 1
    else if (acc&0x2) {
      // SubThreshold KID: no tighter requirement on momentum than >piThr here
      if (subDomain&0x2)                                   lhSubok2 = 1;
      //else... cf. supra.
    }
#  endif
#  ifdef UserEvent_NTuple
    double phD1 = 0, phD2 = 0;                       // Photon density
    if (PIDs[iET1]&0x8) phD1 = trk1.RichInf(9)/thRich1;
    if (PIDs[iET2]&0x8) phD2 = trk2.RichInf(9)/thRich2;
#  endif

    //            ***** COMBINED ID FLAGS *****
#  ifdef U3_USE_CHCHI2
    int chiok_pm = 0;
    if (chiok1==0x3) chiok_pm |= 0x1; if (chiok2==0x3) chiok_pm |= 0x2;
#  else
    int lhok_pm = 0;
    if (lhok1==0x3)  lhok_pm  |= 0x1; if (lhok2==0x3)  lhok_pm  |= 0x2;
#  endif

    if (isphi) {                               // ********** KK = phi **********
      // ***** Histogram characteristics of Ks from ID'ed phi *****
      if (lhok2==0x3 && (acc&0x2)) {                     // K-ID: plot K+
	hR_phiPK[ipd1]->Fill(1/h1(5)); hR_phiAK->Fill(AR1*1000);
	hR_phiRK->Fill(RR1);
      }
      if (lhok1==0x3 && (acc&0x1)) {                     // K+ID: plot K-
	hR_phiPK[ipd2]->Fill(1/h2(5)); hR_phiAK->Fill(AR2*1000);
	hR_phiRK->Fill(RR2);
      }
    }

    // ********** RICH EFFICIENCY HISTOs and LIKELIHOOD PROFILES **********

    if (acc==0x3) {                          // ***** BOTH Hs W/IN RICH... *****
#  ifdef U3_USE_CHCHI2
      hC_phi2->Fill(m_KK,chiok_pm);               // ...COMBINED Eff ESTIMATE #1
#  else
      hL_phi2->Fill(m_KK,lhok_pm);
#  endif
    }

    if (KThr<mom1 && mom1<60) pDomain |= 0x1;
    if (KThr<mom2 && mom2<60) pDomain |= 0x2;
    if ((pDomain&0x1) && mom1<PRICHCut) pRange |= 0x1;
    if ((pDomain&0x2) && mom2<PRICHCut) pRange |= 0x2;
    if (acc&pDomain) {

      // ******************************************************************
      // *************** AT LEAST 1 HADRON w/in RICH, RANGE ***************
      // ******************************************************************

      // ***** Phase Space Binning *****
      int iPi, iP1, iP2;
      for (iPi=iP1 = 0; iPi<npB+1; iPi++) if (mom1-KThr<pBins[iPi]) {
	iP1 = iPi-1; break;
      }
      for (iPi=iP2 = 0; iPi<npB+1; iPi++) if (mom2-KThr<pBins[iPi]) {
	iP2 = iPi-1; break;
      }
      //int iP1 = (int)(mom1-KThr)/10, iP2 =  (int)(mom2-KThr)/10;
      if (iP1<0) iP1 = 0; /* to be sure */ if (iP1>=npB) iP1 = npB-1;
      if (iP2<0) iP2 = 0; /* to be sure */ if (iP2>=npB) iP2 = npB-1;
#  ifdef PS_USE_RADIUS
      int iA1 = (int)(RR1/RRBin), iA2 = (int)(RR2/RRBin);
#  else
      int iA1 = (int)(AR1/ARBin), iA2 = (int)(AR2/ARBin);
#  endif
      if (iA1>=nphiBins) iA1 = nphiBins-1;
      if (iA2>=nphiBins) iA2 = nphiBins-1;
      int iPA1 = 1+iP1+npB*iA1, iPA2 = 1+iP2+npB*iA2;
      time_t time = (time_t)e.UnixSeconds();
      int hour = localtime(&time)->tm_hour;
      int ihour = (hour+21)%24/6;                   // ***** AS A F(TIME of DAY)

      if ((acc&pRange)==0x3) {   // ***** BOTH Hs w/in RICH, w/in RANGE... *****
#  ifdef U3_USE_CHCHI2
	hC_phiM2->Fill(m_KK,chiok_pm);            // ...COMBINED Eff ESTIMATE #2
#  else
	hL_phiM2->Fill(m_KK,lhok_pm);
#  endif
      }

      // *************** chi2 Eff: OVERALL, f(PHASE SPACE,PD,T) ***************

#  ifdef U3_USE_CHCHI2
      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	int ok1 = 0, pv1 = 0;
	if (chiok1&0x1) ok1 = 1; if (chiok1&0x2) pv1 = 1;
	hC_phioP[0][0] ->Fill(m_KK,(double)ok1);
	hC_phioPV[0][0]->Fill(m_KK,(double)pv1);
	hC_phioPT[0][0]->Fill(m_KK,(double)(ok1&pv1));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hC_phioP[0][iPA1] ->Fill(m_KK,(double)ok1);
	hC_phioPV[0][iPA1]->Fill(m_KK,(double)pv1);
	hC_phioPT[0][iPA1]->Fill(m_KK,(double)(ok1&pv1));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hC_phiM[0]->Fill(m_KK,(double)(ok1&pv1));
	  hC_phiPD[0][ipd1]->Fill(m_KK,(double)(ok1&pv1));   // ***** AS A F(PD)
	  hC_phiH[0][ihour]->Fill(m_KK,(double)(ok1&pv1));// AS A F(TIME in DAY)
	}
      }

      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	hR_phioP[1][0]->Fill(RR2); hR_phioP[1][iPA2]->Fill(RR2);
	hR_phioA[1][0]->Fill(mom2-KThr); hR_phioA[1][iPA2]->Fill(mom2-KThr);
	int ok2 = 0, pv2 = 0;
	if (chiok2&0x1) ok2 = 1; if (chiok2&0x2) pv2 = 1;
	hC_phioP[1][0] ->Fill(m_KK,(double)ok2);
	hC_phioPV[1][0]->Fill(m_KK,(double)pv2);
	hC_phioPT[1][0]->Fill(m_KK,(double)(ok2&pv2));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hC_phioP[1][iPA2] ->Fill(m_KK,(double)ok2);
	hC_phioPV[1][iPA2]->Fill(m_KK,(double)pv2);
	hC_phioPT[1][iPA2]->Fill(m_KK,(double)(ok2&pv2));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hC_phiM[1]->Fill(m_KK,(double)(ok2&pv2));
	  hC_phiPD[1][ipd2]->Fill(m_KK,(double)(ok2&pv2));   // ***** AS A F(PD)
	  hC_phiH[1][ihour]->Fill(m_KK,(double)(ok2&pv2));// AS A F(TIME in DAY)
	}
      }
#  endif
#  ifdef U3_PLOT_dTheta
      // *************** dTheta Eff: OVERALL, f(PHASE SPACE, PD) ***************

      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	hT_phioP[0][0] ->Fill(m_KK,(double)(dthok1%2));
	hT_phioPV[0][0]->Fill(m_KK,(double)(dthok1/2));
	hT_phioPT[0][0]->Fill(m_KK,(double)(dthok1/3));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hT_phioP[0][iPA1] ->Fill(m_KK,(double)(dthok1%2));
	hT_phioPV[0][iPA1]->Fill(m_KK,(double)(dthok1/2));
	hT_phioPT[0][iPA1]->Fill(m_KK,(double)(dthok1/3));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hT_phiM[0]->Fill(m_KK,(double)(dthok1/3));
	  hT_phiPD[0][ipd1]->Fill(m_KK,(double)(dthok1/3));  // ***** AS A F(PD)
	}
      }
      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	hT_phioP[1][0] ->Fill(m_KK,(double)(dthok2%2));
	hT_phioPV[1][0]->Fill(m_KK,(double)(dthok2/2));
	hT_phioPT[1][0]->Fill(m_KK,(double)(dthok2/3));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hT_phioP[1][iPA2] ->Fill(m_KK,(double)(dthok2%2));
	hT_phioPV[1][iPA2]->Fill(m_KK,(double)(dthok2/2));
	hT_phioPT[1][iPA2]->Fill(m_KK,(double)(dthok2/3));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hT_phiM[1]->Fill(m_KK,(double)(dthok2/3));
	  hT_phiPD[1][ipd2]->Fill(m_KK,(double)(dthok2/3));  // ***** AS A F(PD)
	}
      }
#  endif
      // ********** LIKELIHOODS **********

      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	hR_phioP[0][0]->Fill(RR1);       hR_phioP[0][iPA1]->Fill(RR1);
	hR_phioA[0][0]->Fill(mom1-KThr); hR_phioA[0][iPA1]->Fill(mom1-KThr);
	//#define DEBUG_SELECTION
#  ifdef DEBUG_SELECTION
	static FILE *f = fopen("sel.log","w");
	if (fabs(m_KK-M_phi)<0.100) fprintf(f,"%d %4d %5d   %d %d %.3f\n",
		Run,e.SpillNum(),e.EvInSpill(),iET1,iET2,m_KK-M_phi);
#  endif
	int ok1 = 0, pv1 = 0;
	if (lhok1&0x1) ok1 = 1; if (lhok1&0x2) pv1 = 1;
	hL_phioP[0][0] ->Fill(m_KK,(double)ok1);
	hL_phioPV[0][0]->Fill(m_KK,(double)pv1);
	hL_phioPT[0][0]->Fill(m_KK,(double)(ok1&pv1));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hL_phioP[0][iPA1] ->Fill(m_KK,(double)ok1);
	hL_phioPV[0][iPA1]->Fill(m_KK,(double)pv1);
	hL_phioPT[0][iPA1]->Fill(m_KK,(double)(ok1&pv1));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hL_phiM[0]->Fill(m_KK,(double)(ok1&pv1));
	  hL_phiPD[0][ipd1]->Fill(m_KK,(double)(ok1&pv1));   // ***** AS A F(PD)
	  hL_phiH[0][ihour]->Fill(m_KK,(double)(ok1&pv1));// AS A F(TIME in DAY)
	}
	if (PIDs[iET1]&0x10) {                        // ...and bckLH available
	  double KLH = richLHs[1][iET1]; double piLH = richLHs[0][iET1];
	  double KpiLH = KLH/piLH;
	  if (pRange&0x1) {
	    hR_LHK[0]->Fill(m_KK,KLH);               // LH profile
	    if (PIDs[iET1]&0x20)                       // ...and pi LH available
	      hR_LHKpi[0]->Fill(m_KK,KpiLH);
	  }
	  hR_LHKvsP[0][iP1]->Fill(m_KK,KLH);           // ...vs P
	}
      }
      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	int ok2 = 0, pv2 = 0;
	if (lhok2&0x1) ok2 = 1; if (lhok2&0x2) pv2 = 1;
	hL_phioP[1][0] ->Fill(m_KK,(double)ok2);
	hL_phioPV[1][0]->Fill(m_KK,(double)pv2);
	hL_phioPT[1][0]->Fill(m_KK,(double)(ok2&pv2));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hL_phioP[1][iPA2] ->Fill(m_KK,(double)ok2);
	hL_phioPV[1][iPA2]->Fill(m_KK,(double)pv2);
	hL_phioPT[1][iPA2]->Fill(m_KK,(double)(ok2&pv2));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hL_phiM[1]->Fill(m_KK,(double)(ok2&pv2));
	  hL_phiPD[1][ipd2]->Fill(m_KK,(double)(ok2&pv2));   // ***** AS A F(PD)
	  hL_phiH[1][ihour]->Fill(m_KK,(double)(ok2&pv2));// AS A F(TIME in DAY)
	}
	if (PIDs[iET2]&0x10) {                        // ...and bckLH available
	  double KLH = richLHs[1][iET2]; double piLH = richLHs[0][iET2];
	  double KpiLH = KLH/piLH;
	  if (pRange&0x2) {
	    hR_LHK[1]->Fill(m_KK,KLH);               // LH profile
	    if (PIDs[iET2]&0x20)                       // ...and pi LH available
	      hR_LHKpi[1]->Fill(m_KK,KpiLH);
	  }
	  hR_LHKvsP[1][iP2]->Fill(m_KK,KLH);           // ...vs P
	}
      }

      if (isphi) {                // *************** KK = phi ***************
	//             ********** RICH RESOLUTION: RESIDUAL...
	if ((acc&pRange&0x1) &&                      // >0 W/IN ACC, p RANGE
	    (thid&0x1) &&                            // ...dtheta of >0
	    thMom1) {                                // >KThr+margin
	  hR_dvsr[2]->Fill(Run,thRich1-thMom1);   // ... vs. Run
	  hR_dKPD[0][ipd1]->Fill(thRich1-thMom1); // ... f(pseudoPhotoDetector)
	  int ia = fabs(YR1)<10?0:1;     // Resol. vs. ThetaC as f(angle,+/-)
	  hR_dKvst[0][ia]->Fill(thMom1,thRich1-thMom1);
#  if defined UserEvent_NTuple && ( UserEvent_NTuple & 2 )
	  double ax1 = h1.Mom3().X()/h1.Mom3().Z();
	  double ay1 = h1.Mom3().Y()/h1.Mom3().Z();
	  ntDth->Fill(thRich1-thMom1,mom1-KThr,
		      XR1,YR1,ax1,ay1,phD1,2);
#  endif
	}
	if ((acc&pRange&0x2) &&                      // <0 W/IN ACC, p RANGE
	    (thid&0x2) &&                            // ...dtheta of <0
	    thMom2) {                                // >KThr+margin
	  hR_dvsr[2]->Fill(Run,thRich2-thMom2);   // ... vs. Run
	  hR_dKPD[1][ipd2]->Fill(thRich2-thMom2); // ... f(pseudoPhotonDetector)
	  int ia = fabs(YR2)<10?0:1;     // Resol. vs. theta as f(angle,+/-)
	  hR_dKvst[1][ia]->Fill(thMom2,thRich2-thMom2);
#  if defined UserEvent_NTuple && ( UserEvent_NTuple & 2 )
	  double ax2 = h2.Mom3().X()/h2.Mom3().Z();
	  double ay2 = h2.Mom3().Y()/h2.Mom3().Z();
	  ntDth->Fill(thRich2-thMom2,mom2-KThr,
		      XR2,YR2,ax2,ay2,phD2,-2);
#  endif
	}
      }

      int rRange = 0;  // Radius cut
      if (RR1>rRICHCut) rRange |= 0x1; if (RR2>rRICHCut) rRange |= 0x2;
      if ((acc&pRange&rRange)==0x3) {// ***** BOTH Hs w/in RICH, pR RANGEs *****
#  ifdef U3_USE_CHCHI2
	hC_phiMA2->Fill(m_KK,chiok_pm);           // ...COMBINED Eff ESTIMATE #3
#  else
	hL_phiMA2->Fill(m_KK,lhok_pm);
#  endif
      }
    }  // End one hadron w/in momentum range

    if (piThr*1.5<mom1 && mom1<KThr) subRange |= 0x1;
    if (piThr*1.5<mom2 && mom2<KThr) subRange |= 0x2;
    if (acc&subRange) {

      // *********************************************************************
      // *************** AT LEAST 1 HADRON w/in RICH, subRANGE ***************
      // *********************************************************************

      // ********** SUBTHRESHOLD Chi2 and LIKELIHOODS **********

      if (acc&subRange&0x1) {    // ***** >0 w/in RICH, sub RANGE *****
#  ifdef U3_PLOT_CHi2
	hC_phiS[0]->Fill(m_KK,(double)(chiSubok1));
#  endif
	hL_phiS[0]->Fill(m_KK,(double)(lhSubok1));
	if (PIDs[iET1]&0x10) {                        // ...and bckLH available
	  double piLH = richLHs[0][iET1];
	  hR_LHSpi[0]->Fill(m_KK,piLH);                 // => piLH Profile
	}
      }
      if (acc&subRange&0x2) {    // ***** <0 W/IN RICH, SUB RANGE *****
#  ifdef U3_PLOT_CHi2
	hC_phiS[1]->Fill(m_KK,(double)(chiSubok2));
#  endif
	hL_phiS[1]->Fill(m_KK,(double)(lhSubok2));
	if (PIDs[iET2]&0x10) {                        // ...and bckLH available
	  double piLH = richLHs[0][iET2];
	  hR_LHSpi[1]->Fill(m_KK,piLH);                 // => piLH Profile
	}
      }
    }  // End one hadron sub momentum range

  }  // End one hadron w/ RICH acceptance
}
// *******************************************************************
// ************************* bookK0_RICHPerf *************************
// *******************************************************************
void bookK0_RICHPerf(double dM_K0, double K0MassCut)
{
  // ******************** K0 * RICH PERFS ********************
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
    gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();
  TDirectory *dRICH; if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs")))
    gDirectory->mkdir("RICHPerfs","RICH Performances from K0 decay");
  if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs"))) {
    printf("\n** U3:\a No creating subdir \"RICHPerfs\" in TDirectory \"%s/K0\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dRICH->cd();

#  ifdef U3_USE_CHCHI2
  char cChi2Cut[] = "#pi#chi^{2}<5.0,<0.95K,p  ";
  sprintf(cChi2Cut,"#pi#chi^{2}<%.1f,<K,p",CHChi2Mx); string sChi2piCut(cChi2Cut);
  // Cuts w/ pi veto
  sprintf(cChi2Cut,"#pi#chi^{2}<%.1f,<%.2fK,p",CHChi2Mx,CHChi2Cut);
  string sChi2Cut(cChi2Cut);
  sprintf(cChi2Cut,"#pi#chi^{2}<%.2fK,p",CHChi2Cut);  string sChi2VCut(cChi2Cut);
#  endif

  char cLHCut[] = "#piLH>1.005bck,>0.999K,p  ";
  int precLHBC = LHBckCut-1 ? int(-log10(fabs(LHBckCut-1))+.999) : 0;
  int precLHC  = LHCut-1    ? int(-log10(fabs(LHCut   -1))+.999) : 0;
  sprintf(cLHCut,"#piLH>%#.*fbck,>K,p",precLHBC,LHBckCut);
  string sLHpiCut(cLHCut);
  sprintf(cLHCut,"#piLH>%#.*fbck,>%#.*fK,p",precLHBC,LHBckCut,precLHC,LHCut);
  string sLHCut(cLHCut);
  sprintf(cLHCut,"#piLH>%#.*fK,p",precLHC,LHCut);
  string sLHVCut(cLHCut);
#  ifdef piS_e_REJECTION
  char cLHeVeto[] = ",eVeto(1.5)  ";
  sprintf(cLHeVeto,",eVeto(%.1f)",piLHeVeto); string sLHeVeto(cLHeVeto);
#  endif

#  ifdef U3_PLOT_dTheta
  char cdTpiCut[] = "#delta_{#pi}#theta<2.5#sigma  ";
  sprintf(cdTpiCut,"#delta_{#pi}#theta<%.1f#sigma",nSigmas);
  char cKVCut[] = "#delta_{K}#theta>2.5#sigma  ";
  sprintf(cKVCut,"#delta_{K}#theta>%.1f#sigma",nVeto);
  string sdTpiCut(cdTpiCut), sKVCut(cKVCut);
  string sdTCut(cdTpiCut+string("#times")+cKVCut);
#  endif
#  ifdef RICH_ANGLE_CUT
  char cRCut[] = ",A>25"; sprintf(cRCut,",R>%.0f",1000*aRICHCut);
#  else    // ***** Radius cut *****
  char cRCut[] = ",R>12"; sprintf(cRCut,",R>%.0f",rRICHCut);
#  endif

  string sRCut(cRCut);
  char cPRange[] = ",#piThr<P<100GeV "; sprintf(cPRange,",#piThr<P<%.0fGeV",PRICHCut);
  string pipi("#pi+#pi- - "), sPRange(cPRange);
  string spiID("#pi#pm:");

  //   ********** COMBINED +/- HISTOS FOR CORRELATED EFFICIENCIES **********
#  ifdef U3_PLOT_CHi2
  //           *************** Chi2 ***************
  // Raw efficiency
  hC_K02   = new TH2D("hC_K0",  string(pipi+spiID+sChi2Cut).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  // Momenta over threshold
  hC_K0M2  = new TH2D("hC_K0M", string(pipi+spiID+sChi2Cut+sPRange).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  // Momenta over threshold AND radius(angle) < limit
  hC_K0MA2 = new TH2D("hC_K0MA",string(pipi+spiID+sChi2Cut+sPRange+sRCut).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
#  else
  //            *************** LH ***************
#    ifdef piS_e_REJECTION
  hL_K02   = new TH2D("hL_K0",  string(pipi+spiID+sLHCut+sLHeVeto).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
#    else
  hL_K02   = new TH2D("hL_K0",  string(pipi+spiID+sLHCut).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
#    endif
  // Momenta over threshold
  hL_K0M2  = new TH2D("hL_K0M", string(pipi+spiID+sLHCut+sPRange).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
  // Momenta over threshold AND radius(angle) < limit
  hL_K0MA2 = new TH2D("hL_K0MA",string(pipi+spiID+sLHCut+sPRange+sRCut).c_str(),
		      150,M_K0-dM_K0,M_K0+dM_K0,4,-.5,3.5);
#  endif

  // *************** Chi2, dTheta, LikeliHood ***************
  // (Successively pi+ and pi-)
  char cpipm[] = "#pi+:", hName[] = "hC_K0Mp";
  for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON pi+ pi- **********
    sprintf(cpipm,"#pi%c:",pm?'-':'+'); string spm(cpipm);
    // ***** piID Eff AND Purity FOR MOMENTA OVER THRESHOLD *****
#  ifdef U3_PLOT_CHi2
    // chi2
    sprintf(hName,"hC_K0M%c",pm?'m':'p');
    hC_K0M[pm] = new TH2D(hName,string(pipi+spm+sChi2Cut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
    // dTheta
    sprintf(hName,"hT_K0M%c",pm?'m':'p');
    hT_K0M[pm] = new TH2D(hName,string(pipi+spm+sdTCut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
    // Likelihood
    sprintf(hName,"hL_K0M%c",pm?'m':'p');
    hL_K0M[pm] = new TH2D(hName,string(pipi+spm+sLHCut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
  }
  for (int ipa = 0; ipa<nK0B; ipa++) {

    // *************** piID EFF FOR P>THRES AS A F(PHASE SPACE) ***************

    char hN[] = "hC_K0oPTp10";
    int PMn, PMx; char cABin[] = ",12<A<24";
    if (ipa==0) { PMn = 0;           PMx = 60; }
    else { 
      //PMn = (ipa-1)%3*10; PMx = PMn+10;
      int ip = (ipa-1)%npB; PMn = pBins[ip]; PMx = pBins[ip+1];
    }
#  ifdef PS_USE_RADIUS
    char V = 'R'; double RBin = RRBin;
#  else
    char V = 'A'; double RBin = ARBin*1000;
#  endif
    if      (ipa==0)              sprintf(cABin,"%s","\0");
    else if (0<ipa && ipa<=npB)   sprintf(cABin,",%c<%.0f",V,RBin); 
    else if (npB*(nK0Bins-1)<ipa) sprintf(cABin,",%.0f<%c",(ipa-1)/npB*RBin,V);
    else                          sprintf(cABin,",%.0f<%c<%.0f",
				       (ipa-1)/npB*RBin,V,((ipa-1)/npB+1)*RBin);
    char cPBin[] = ",15<P-#piThr<60GeV";
    if      (ipa==0)             sprintf(cPBin,",#piThr<P<60GeV");
    else if ((ipa-1)%npB==npB-1) sprintf(cPBin,",#piThr+%d<P<60GeV",PMn);
    else                         sprintf(cPBin,",%d<P-#piThr<%dGeV",PMn,PMx);
    string sABin(cABin), sPBin(cPBin), sPABin(sPBin+sABin);

    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON pi+ pi- **********
      sprintf(cpipm,"#pi%c:",pm?'-':'+'); string spm(cpipm);

#  ifdef U3_PLOT_CHi2
      // ********** USING CHI2 CUTS **********
      sprintf(hN,"hC_K0oP%c%d",pm?'m':'p',ipa);
      hC_K0oP[pm][ipa] = new TH2D(hN,string(pipi+spm+sChi2piCut+sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: pi ID + K Veto
      sprintf(hN,"hC_K0oPT%c%d",pm?'m':'p',ipa);
      hC_K0oPT[pm][ipa]= new TH2D(hN,string(pipi+spm+sChi2Cut  +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Enhance K rejection (K Veto)
      sprintf(hN,"hC_K0oPV%c%d",pm?'m':'p',ipa);
      hC_K0oPV[pm][ipa]= new TH2D(hN,string(pipi+spm+sChi2VCut +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      // ********** USING dTheta CUTS **********
      sprintf(hN,"hT_K0oP%c%d",pm?'m':'p',ipa);
      hT_K0oP[pm][ipa] = new TH2D(hN,string(pipi+spm+sdTpiCut+sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: pi ID + K Veto
      sprintf(hN,"hT_K0oPT%c%d",pm?'m':'p',ipa);
      hT_K0oPT[pm][ipa]= new TH2D(hN,string(pipi+spm+sdTCut  +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // K Veto Only
      sprintf(hN,"hT_K0oPV%c%d",pm?'m':'p',ipa);
      hT_K0oPV[pm][ipa]= new TH2D(hN,string(pipi+spm+sKVCut  +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
      // ***** USING LH CUTS *****
      sprintf(hN,"hL_K0oP%c%d",pm?'m':'p',ipa);
      hL_K0oP[pm][ipa] = new TH2D(hN,string(pipi+spm+sLHpiCut+sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: pi ID + K Veto
      sprintf(hN,"hL_K0oPT%c%d",pm?'m':'p',ipa);
      hL_K0oPT[pm][ipa]= new TH2D(hN,string(pipi+spm+sLHCut  +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // K enhanced rejection (K veto)
      sprintf(hN,"hL_K0oPV%c%d",pm?'m':'p',ipa);
      hL_K0oPV[pm][ipa]= new TH2D(hN,string(pipi+spm+sLHVCut +sPABin).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    }
  }
  for (int ipd = 0; ipd<4; ipd++) {

    // *************** piID EFF FOR P > THRES AS A F(pseudoPD) ***************

    char hN[] = "hT_K0PDp0", cPD[] = ",PD0";
    sprintf(cPD,",PD%d",ipd); string sPD(cPD);
    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cpipm,"#pi%c:",pm?'-':'+'); string spm(cpipm);

#  ifdef U3_PLOT_CHi2
      sprintf(hN,"hC_K0PD%c%d",pm?'m':'p',ipd);
      hC_K0PD[pm][ipd] = new TH2D(hN,string(pipi+spm+sChi2Cut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      sprintf(hN,"hT_K0PD%c%d",pm?'m':'p',ipd);
      hT_K0PD[pm][ipd] = new TH2D(hN,string(pipi+spm+sdTCut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
      sprintf(hN,"hL_K0PD%c%d",pm?'m':'p',ipd);
      hL_K0PD[pm][ipd]= new TH2D(hN,string(pipi+spm+sLHCut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    }
  }

  for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON pi+ pi- **********
    char hName[] = "hR_dpivstp01";
    char hTitle[] =   "RICH #pi+#delta#theta vs. #theta - K0#pm20Mev,PD0,#pi-ID,D_{#nu}<.3,A<100  ";
    //char hTitle[] = "#pi+ (X,Y) @ Z=750 - K0#pm20Mev  ";  // XY
    char cpm = pm?'-':'+';

    //    *************** CERENKOV ANGLE RESOLUTIONS... ***************
    //    ***** (N.B.: w/ ID of -/+ counterpart) *****

    for (int ipd = 0; ipd<16; ipd++) {         // ***** ...as a f(pseudoPD, +/-)
      sprintf(hName,"hR_dpiPD%c%d",pm?'m':'p',ipd);
      sprintf(hTitle,
	      "RICH #pi%c#delta#theta - K0#pm%.0fMev,PD%d,#pi%cID",
	      cpm,K0MassCut*1000,ipd,pm?'+':'-');
      hR_dpiPD[pm][ipd] = new TH1D(hName,hTitle,100,-10,10);
    }
    //                   ***** vs. Theta as a f(photonDensity, +/-, SM1/2) *****
    for (int ia = 0; ia<2; ia++) {
      for (int iph = 0; iph<2; iph++) {
	sprintf(hName,"hR_dpivst%c%d%d",pm?'m':'p',iph,ia);
	sprintf(hTitle,"RICH #pi%c#delta#theta vs. #theta - K0#pm%.0fMev,#pi%cID,D_{#nu}%c.3,%s",
	 cpm,K0MassCut*1000,pm?'+':'-',iph?'>':'<',ia?"Y>10":"Y<10");
	hR_dpivst[pm][iph][ia] = new TH2D(hName,hTitle,16,0,60,120,-6,6);
      }
    }
    sprintf(hName,"hR_nphvsp%c",pm?'m':'p');    // ***** # of photon vs. p *****
    sprintf(hTitle,
	    "RICH #pi%c##nu vs. P-Thr - K0#pm%.0fMev,#pi%cID",
	    cpm,K0MassCut*1000,pm?'+':'-');
    hR_nphvsp[pm] = new TH2D(hName,hTitle,16,0,40,50,0,1);
    sprintf(hName,"hR_chivsp%c",pm?'m':'p');           // ***** chi2 vs. p *****
    sprintf(hTitle,
	    "RICH #pi%c#chi^{2} vs. P-Thr - K0#pm%.0fMev,#pi%cID",
	    cpm,K0MassCut*1000,pm?'m':'p');
    hR_chivsp[pm] = new TH2D(hName,hTitle,16,0,40,50,0,5);

    // *************** LIKELIHOOD PROFILES ***************
    // w/ RICH and P range
    sprintf(hName,"hR_LHpi%c",pm?'m':'p');
    sprintf(hTitle,"RICH #pi%cLH - #piThr<p<%.0fGeV",pm?'-':'+',PRICHCut);
    hR_LHpi[pm]  = new TProfile(hName,hTitle,50,M_K0-dM_K0,M_K0+dM_K0,"S");
  }
}
//*******************************************************************
// ************************* bookK0_RICHPur *************************
//*******************************************************************
void bookK0_RICHPur(double dM_K0)
{
  // ******************** K0 * RICH PURITY ********************
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
    gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();
  TDirectory *dRICH; if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs")))
    gDirectory->mkdir("RICHPerfs","RICH Performances from K0 decay");
  if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs"))) {
    printf("\n** U3:\a No creating subdir \"RICHPerfs\" in TDirectory \"%s/K0\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dRICH->cd();

#  ifdef U3_USE_CHCHI2
  char cChi2Cut[] = "K#chi^{2}<5.0,<0.95#pi,p  ";
  sprintf(cChi2Cut,"K#chi^{2}<%.1f,<#pi,p",CHChi2Mx); string sChi2KCut(cChi2Cut);
  // Cuts w/ pi veto
  sprintf(cChi2Cut,"K#chi^{2}<%.1f,<%.2f#pi,p",CHChi2Mx,CHChi2Cut);
  string sChi2Cut(cChi2Cut);
  sprintf(cChi2Cut,"K#chi^{2}<%.2f#pi,p",CHChi2Cut);  string sChi2VCut(cChi2Cut);
  sprintf(cChi2Cut,"#pi#chi^{2}>%.1f",subKThrChi2piVeto);
  string sChi2SCut(cChi2Cut);
#  endif
  char cLHCut[] = "KLH>1.005bck,>0.999#pi,p  ";
  int precLHBC = LHBckCut-1 ? int(-log10(fabs(LHBckCut-1))+.999) : 0;
  int precLHC  = LHCut-1    ? int(-log10(fabs(LHCut   -1))+.999) : 0;
  sprintf(cLHCut,"KLH>%#.*fbck,>#pi,p",precLHBC,LHBckCut);
  string sLHKCut(cLHCut);
  sprintf(cLHCut,"KLH>%#.*fbck,>%#.*f#pi,p",precLHBC,LHBckCut,precLHC,LHCut);
  string sLHCut(cLHCut);
  sprintf(cLHCut,"KLH>%#.*f#pi,p",precLHC,LHCut);
  string sLHVCut(cLHCut);
  int precLHSC = subKThrLHpiVeto-1 ? int(-log10(fabs(subKThrLHpiVeto-1))+.999)
    : 0;
  sprintf(cLHCut,"#piLH<%#.*fbck",precLHSC,subKThrLHpiVeto);
  string sLHSCut(cLHCut);

#  ifdef U3_PLOT_dTheta
  char cdTKCut[] = "#delta_{K}#theta<2.5#sigma  ";
  sprintf(cdTKCut,"#delta_{K}#theta<%.1f#sigma",nSigmas);
  char cpiVCut[] = "#delta_{#pi}#theta>2.5#sigma  ";
  sprintf(cpiVCut,"#delta_{#pi}#theta>%.1f#sigma",nVeto);
  string sdTKCut(cdTKCut), spiVCut(cpiVCut);
  string sdTCut(cdTKCut+string("#times")+cpiVCut);
#  endif
#  ifdef RICH_ANGLE_CUT
  char cRCut[] = ",A>25"; sprintf(cRCut,",R>%.0f",1000*aRICHCut);
#  else    // ***** Radius cut *****
  char cRCut[] = ",R>12"; sprintf(cRCut,",R>%.0f",rRICHCut);
#  endif
  string sRCut(cRCut);

  char cPRange[] = ",#piThr<P<100GeV"; sprintf(cPRange,",#piThr<P<%.0fGeV",PRICHCut);
  string pipi("#pi+#pi- - "), sPRange(cPRange);
  char cPSub[] = ",1.50#piThr<P<KThr"; sprintf(cPSub,",%4.2f#piThr<P<KThr",thrMargin);
  string sPSub(cPSub);

  // *************** Chi2, dTheta, LikeliHood ***************
  // (Successively K+ and K-)
  char cKpm[] = "K+:", hName[] = "hC_K0Mp";
  for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
    sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);
    // ***** K misID Eff FOR MOMENTA OVER THRESHOLD *****
#  ifdef U3_PLOT_CHi2
    // chi2
    sprintf(hName,"hD_K0M%c",pm?'m':'p');
    hD_K0M[pm] = new TH2D(hName,string(pipi+spm+sChi2Cut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    // SubThreshold chi2
    sprintf(hName,"hD_K0S%c",pm?'m':'p');
    hD_K0S[pm] = new TH2D(hName,string(pipi+spm+sChi2SCut+sPSub).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
    // dTheta
    sprintf(hName,"hU_K0M%c",pm?'m':'p');
    hU_K0M[pm] = new TH2D(hName,string(pipi+spm+sdTCut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
    // Likelihood
    sprintf(hName,"hM_K0M%c",pm?'m':'p');
    hM_K0M[pm] = new TH2D(hName,string(pipi+spm+sLHCut+sPRange).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    // SubThreshold Likelihood
    sprintf(hName,"hM_K0S%c",pm?'m':'p');
    hM_K0S[pm] = new TH2D(hName,string(pipi+spm+sLHSCut+sPSub).c_str(),
			  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
  }
  for (int ipa = 0; ipa<nphB; ipa++) {

    // *************** KmisID EFF FOR P>THR AS A F(PHASE SPACE) ***************

    char hN[] = "hC_K0oPTp10";
    int PMn, PMx; char cABin[] = ",12<A<24";
    if (ipa==0) { PMn = 0;           PMx = 60; }
    else { 
      //PMn = (ipa-1)%3*10; PMx = PMn+10;
      int ip = (ipa-1)%npB; PMn = pBins[ip]; PMx = pBins[ip+1];
    }
#  ifdef PS_USE_RADIUS
    char V = 'R'; double RBin = RRBin;
#  else
    char V = 'A'; double RBin = ARBin*1000;
#  endif
    if      (ipa==0)               sprintf(cABin,"%s","\0");
    else if (0<ipa && ipa<=npB)    sprintf(cABin,",%c<%.0f",V,RBin); 
    else if (npB*(nphiBins-1)<ipa) sprintf(cABin,",%.0f<%c",(ipa-1)/npB*RBin,V);
    else                           sprintf(cABin,",%.0f<%c<%.0f",
				       (ipa-1)/npB*RBin,V,((ipa-1)/npB+1)*RBin);
    char cPBin[] = ",15<P-KThr<60GeV";
    if      (ipa==0)             sprintf(cPBin,",KThr<P<60GeV");
    else if ((ipa-1)%npB==npB-1) sprintf(cPBin,",KThr+%d<P<60GeV",PMn);
    else                         sprintf(cPBin,",%d<P-KThr<%dGeV",PMn,PMx);
    string sABin(cABin), sPBin(cPBin), sPABin(sPBin+sABin);

    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_PLOT_CHi2
      // ********** USING CHI2 CUTS **********
      sprintf(hN,"hD_K0oP%c%d",pm?'m':'p',ipa);
      hD_K0oP[pm][ipa]  = new TH2D(hN,string(pipi+spm+sChi2KCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hD_K0oPT%c%d",pm?'m':'p',ipa);
      hD_K0oPT[pm][ipa] = new TH2D(hN,string(pipi+spm+sChi2Cut +sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Enhance pi rejection (pi Veto)
      sprintf(hN,"hD_K0oPV%c%d",pm?'m':'p',ipa);
      hD_K0oPV[pm][ipa] = new TH2D(hN,string(pipi+spm+sChi2VCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      // ********** USING dTheta CUTS **********
      sprintf(hN,"hU_K0oP%c%d",pm?'m':'p',ipa);
      hU_K0oP[pm][ipa]  = new TH2D(hN,string(pipi+spm+sdTKCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hU_K0oPT%c%d",pm?'m':'p',ipa);
      hU_K0oPT[pm][ipa] = new TH2D(hN,string(pipi+spm+sdTCut +sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // pi Veto Only
      sprintf(hN,"hU_K0oPV%c%d",pm?'m':'p',ipa);
      hU_K0oPV[pm][ipa] = new TH2D(hN,string(pipi+spm+spiVCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
      // Radius Distribution
      sprintf(hN,"hR_K0oP%c%d",pm?'m':'p',ipa);
      hR_K0oP[pm][ipa]= new TH1D(hN,string(pipi+spm+sPABin).c_str(),
				 80,0,80);
      // Momentum Distribution
      sprintf(hN,"hR_K0oA%c%d",pm?'m':'p',ipa);
      hR_K0oA[pm][ipa]= new TH1D(hN,string(pipi+spm+sPABin).c_str(),
				 80,0,40);
      // ********** USING LH CUTS **********
      sprintf(hN,"hM_K0oP%c%d",pm?'m':'p',ipa);
      hM_K0oP[pm][ipa]  = new TH2D(hN,string(pipi+spm+sLHKCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Total: K ID + pi Veto
      sprintf(hN,"hM_K0oPT%c%d",pm?'m':'p',ipa);
      hM_K0oPT[pm][ipa] = new TH2D(hN,string(pipi+spm+sLHCut +sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
      // Enhanced pi rejection (pi Veto)
      sprintf(hN,"hM_K0oPV%c%d",pm?'m':'p',ipa);
      hM_K0oPV[pm][ipa] = new TH2D(hN,string(pipi+spm+sLHVCut+sPABin).c_str(),
				   150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    }
  }
  for (int ipd = 0; ipd<4; ipd++) {

    // *************** KmisID EFF FOR P > THRES AS A F(pseudoPD) ***************

    char hN[] = "hT_K0PDp0", cPD[] = ",PD0";
    sprintf(cPD,",PD%d",ipd); string sPD(cPD);
    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON K+ K- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_PLOT_CHi2
      sprintf(hN,"hD_K0PD%c%d",pm?'m':'p',ipd);
      hD_K0PD[pm][ipd] = new TH2D(hN,string(pipi+spm+sChi2Cut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
#  ifdef U3_PLOT_dTheta
      sprintf(hN,"hU_K0PD%c%d",pm?'m':'p',ipd);
      hU_K0PD[pm][ipd] = new TH2D(hN,string(pipi+spm+sdTCut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
      sprintf(hN,"hM_K0PD%c%d",pm?'m':'p',ipd);
      hM_K0PD[pm][ipd] = new TH2D(hN,string(pipi+spm+sLHCut+sPD).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    }
  }
  for (int ih = 0; ih<4; ih++) {

    // ********** KmisID EFF FOR P > THRES AS A F(time of day) **********

    char hN[] = "hT_K0Hp0";
    char h1[] = "03:09", h2[] = "09:15", h3[] = "15:21", h4[] = "21:03";
    char *hours[4] = { h1, h2, h3, h4}; string sH(hours[ih]);
    for (int pm = 0; pm<2; pm++) {   // ********** LOOP ON pi+ pi- **********
      sprintf(cKpm,"K%c:",pm?'-':'+'); string spm(cKpm);

#  ifdef U3_PLOT_CHi2
      sprintf(hN,"hD_K0H%c%d",pm?'m':'p',ih);
      hD_K0H[pm][ih] = new TH2D(hN,string(pipi+spm+sChi2Cut+sH).c_str(),
				150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
#  endif
      sprintf(hN,"hM_K0H%c%d",pm?'m':'p',ih);
      hM_K0H[pm][ih] = new TH2D(hN,string(pipi+spm+sLHCut+sH).c_str(),
				  150,M_K0-dM_K0,M_K0+dM_K0,2,-.5,1.5);
    }
  }
}
// ********************************************************************
// *************************  fillK0_RICHPerf *************************
// ********************************************************************
void fillK0_RICHPerf(PaEvent &e, int iET1, int iET2,
		     const TLorentzVector &lvpi1, const TLorentzVector &lvpi2,
		     double m_pipi, bool isK0, bool cleanV0, int Run)
{
  // ********** DETERMINE TRACKS' ACCEPTANCE by RICH **********

  const PaTrack &trk1 = e.vTrack()[iET1], &trk2 = e.vTrack()[iET2];
  // Start hlx
  const PaTPar &h1 = trk1.vTPar()[0], &h2 = trk2.vTPar()[0];
  double XR1, YR1, RR1, XR2, YR2, RR2;
  double *XRs[2] = {&XR1,&XR2}, *YRs[2] = {&YR1,&YR2}, *RRs[2] = {&RR1,&RR2};
#  ifndef PS_USE_RADIUS
  double AR1, AR2, *ARs[2] = {&AR1,&AR2};
#  endif
  int acc = 0, pRange = 0, pDomain = 0; // Acceptance and p Range (+/-) Flags...
  int pm, iET; const PaTrack *trk; for (pm = 0, iET = iET1, trk = &trk1; pm<2;
					pm++) {
    if (aRs[iET]==0) {
 printf("** U3/K0 RICHPerf:\a Evt %d,%d Track %d RICH not yet checked\n\n",
	e.RunNum(),(int)e.UniqueEvNum(),iET); assert(false);
    }
    if (aRs[iET]>0) acc |= 0x1<<pm;
    const vector<Float_t> &aux = trk->vAux();
    double XR = aux[0], YR = aux[1], RR = sqrt(XR*XR+YR*YR);
    *XRs[pm] = XR; *YRs[pm] = YR; *RRs[pm] = RR;
#  ifndef PS_USE_RADIUS
    float tgx = aux[2], tgy = aux[3];
    *ARs[pm] = acos(1/sqrt(1.+ tgx*tgx + tgy*tgy));
#  endif
    iET = iET2; trk = &trk2;
  }
  // More precision: which part of acceptance is hit...
  int ipd1  = 0, ipd2 = 0;   // .. PD = pseudo Photon Detector
  if (XR1>0) ipd1 |= 0x1; if (YR1>0) ipd1 |= 0x2;
  if (XR2>0) ipd2 |= 0x1; if (YR2>0) ipd2 |= 0x2;
  int jpd1 = ipd1, jpd2 = ipd2;
  if (fabs(XR1)>45) jpd1 += 4; if (fabs(YR1)>45) jpd1 += 8;
  if (fabs(XR2)>45) jpd2 += 4; if (fabs(YR2)>45) jpd2 += 8;
  double mom1 = h1.Mom(), mom2 = h2.Mom();

#  ifdef DISCARD_MUONS             // ********** DISCARD MUONs **********
  if (trk1.XX0()>15) acc &= 0x2; if (trk2.XX0()>15) acc &= 0x1;
#  endif

  if (acc) {
    // *********************************************************************
    // ******************** AT LEAST 1 HADRON W/IN RICH ********************
    // *********************************************************************

    // ***** SET A NUMBER OF RICH-RELATED VARIABLES *****
    int thid = 0;                // "th" = theta availability of theta
    if (PIDs[iET1]&0x8) thid |= 0x1; if (PIDs[iET2]&0x8) thid |= 0x2;
    double thRich1 = 0, thRich2 = 0, thMom1 = 0, thMom2 = 0;
#  ifdef U3_PLOT_dTheta  // deltaTheta pi(K)ID: Retain only 2 LS bits
    unsigned int dthKok1 = rich->KID (trk1,thRich1,thMom1)&0x3;
    unsigned int dthKok2 = rich->KID (trk2,thRich2,thMom2)&0x3;
    // N.B.: The 2 instructions infra set "thMom"
    unsigned int dthok1  = rich->piID(trk1,thRich1,thMom1)&0x3;
    unsigned int dthok2  = rich->piID(trk2,thRich2,thMom2)&0x3;
#  else
    if (PIDs[iET1]&0x8) thRich1 = trk1.RichInf(RICH_ESTIMATOR);
    if (PIDs[iET2]&0x8) thRich2 = trk2.RichInf(RICH_ESTIMATOR);
    if (piThr+.001<mom1)
      thMom1 = acos(sqrt(M2_pi/mom1/mom1+1)/rich->CurrentIndAPV)*1000;
    if (piThr+.001<mom2)
      thMom2 = acos(sqrt(M2_pi/mom2/mom2+1)/rich->CurrentIndAPV)*1000;
    /*
    double th1 = thMom1, th2 = thMom2;
    rich->piID(trk1,thRich1,thMom1)&0x3; rich->piID(trk2,thRich2,thMom2)&0x3;
    if (th1 && (PIDs[iET1]&0x8) && fabs(thMom1-th1)>.00001) 
      printf("Inconsistency: %.2f => %.4f != %.4f\n",mom1,thMom1,th1);
    if (th2 && (PIDs[iET2]&0x8) && fabs(thMom2-th2)>.00001) 
      printf("Inconsistency: %.2f => %.4f != %.4f\n",mom2,thMom2,th2);
    */
#  endif
#  ifdef U3_USE_CHCHI2
    // ***** pi(K)ID: chi2 method *****
    int chiok1 = 0, chiok2 = 0;
    int chiKok1 = 0, chiKok2 = 0, chiSubKok1 = 0, chiSubKok2 = 0;
    if (PIDs[iET1]&0x10) {
      double pichi2 = richChi2s[0][iET1], Kchi2  = richChi2s[1][iET1];
      double pchi2  = richChi2s[2][iET1];
      if (pichi2<          Kchi2 && pichi2<          pchi2 &&
	  pichi2<CHChi2Mx)                                  chiok1 |= 0x1;
      if (pichi2<CHChi2Cut*Kchi2 && pichi2<CHChi2Cut*pchi2) chiok1 |= 0x2;
      if (Kchi2<          pichi2 && Kchi2<           pchi2 &&
	  Kchi2<CHChi2Mx)                                   chiKok1 |= 0x1;
      if (Kchi2<CHChi2Cut*pichi2 && Kchi2 <CHChi2Cut*pchi2) chiKok1 |= 0x2;
      if (pichi2>subKThrChi2piVeto) {
#    ifdef REJECT_SubKThr_e
	bool eid = richChi2s[3][iET1]<subKThrChi2piVeto;
#    else
	bool eid = false;
#    endif
	if (!eid) chiSubKok1 = 1;
      }
    }
    if (PIDs[iET2]&0x10) {
      double pichi2 = richChi2s[0][iET2], Kchi2  = richChi2s[1][iET2];
      double pchi2  = richChi2s[2][iET2];
      if (pichi2<          Kchi2 && pichi2<          pchi2 &&
	  pichi2<CHChi2Mx)                                  chiok2 |= 0x1;
      if (pichi2<CHChi2Cut*Kchi2 && pichi2<CHChi2Cut*pchi2) chiok2 |= 0x2;
      if (Kchi2<          pichi2 && Kchi2<           pchi2 &&
	  Kchi2<CHChi2Mx)                                   chiKok2 |= 0x1;
      if (Kchi2<CHChi2Cut*pichi2 && Kchi2 <CHChi2Cut*pchi2) chiKok2 |= 0x2;
      if (pichi2>subKThrChi2piVeto) {
#    ifdef REJECT_SubKThr_e
	bool eid = richChi2s[3][iET2]<subKThrChi2piVeto;
#    else
	int eid = 0;
#    endif
	if (!eid) chiSubKok2 = 1;
      }
    }
#  endif
    // ***** pi(K)ID: LH method *****
    int lhok1 = 0, lhok2 = 0;
    int lhKok1 = 0, lhKok2 = 0, lhSubKok1 = 0, lhSubKok2 = 0;
#  ifdef piS_e_REJECTION
    int lhpiSok1 = 0, lhpiSok2 = 0;
#  endif
    if (PIDs[iET1]&0x10) {
      double piLH = richLHs[0][iET1], KLH  = richLHs[1][iET1];
      double pLH  = richLHs[2][iET1], eLH  = richLHs[3][iET1];
      if (piLH>      KLH && piLH>      pLH && piLH>LHBckCut) lhok1 |= 0x1;
      if (KLH>      piLH && KLH>       pLH &&  KLH>LHBckCut) lhKok1 |= 0x1;
      if (piLH>LHCut*KLH && piLH>LHCut*pLH)                  lhok1 |= 0x2;
      if (KLH>LHCut*piLH && KLH >LHCut*pLH)                  lhKok1 |= 0x2;
      if (piLH<subKThrLHpiVeto) {
#  ifdef REJECT_SubKThr_e
	bool eid = eLH>subKThrLHeVeto;
#  else
	bool eid = false;
#  endif
	if (!eid) lhSubKok1 = 1;
      }
#  ifdef piS_e_REJECTION
      lhpiSok1 = lhok1;
      if (eLH>piLHeVeto*piLH) {
	/* */                                  lhpiSok1 = 0; lhok1 = 0;
      }
#  endif
    }
#  if K_BELOW_THR > 1
    else if (acc&0x1) {
      // SubThreshold KID: no tighter requirement on momentum than >piThr here,
      // but later on the "lhSubKok1" will be enabled only in a resticted range.
      if      (piThr<mom1)                                   lhSubKok1 = 1;
      //else... in fact, no real need to ID pion below threshold any further
      //       than what's achieved by eVeto.
    }
#  endif
    if (PIDs[iET2]&0x10) {
      double piLH = richLHs[0][iET2], KLH  = richLHs[1][iET2];
      double pLH  = richLHs[2][iET2], eLH  = richLHs[3][iET2];
      if (piLH>      KLH && piLH>      pLH && piLH>LHBckCut) lhok2 |= 0x1;
      if (KLH>      piLH && KLH>       pLH &&  KLH>LHBckCut) lhKok2 |= 0x1;
      if (piLH>LHCut*KLH && piLH>LHCut*pLH)                  lhok2 |= 0x2;
      if (KLH>LHCut*piLH && KLH >LHCut*pLH)                  lhKok2 |= 0x2;
      if (piLH<subKThrLHpiVeto) {
#  ifdef REJECT_SubKThr_e
	bool eid = eLH>subKThrLHeVeto;
#  else
	bool eid = false;
#  endif
	if (!eid) lhSubKok2 = 1;
      }
#  ifdef piS_e_REJECTION
      lhpiSok2 = lhok2;
      if (eLH>piLHeVeto*piLH) {
	/* */                                  lhpiSok2 = 0; lhok2 = 0;
      }
#  endif
    }
#  if K_BELOW_THR > 1
    else if (acc&0x2) {
      // SubThreshold KID: no tighter requirement on momentum than >piThr here
      if      (piThr<mom2)                                   lhSubKok2 = 1;
      //else... cf. supra.
    }
#  endif
    double phD1 = 0, phD2 = 0;                       // Photon density
    if (PIDs[iET1]&0x8) phD1 = trk1.RichInf(9)/thRich1;
    if (PIDs[iET2]&0x8) phD2 = trk2.RichInf(9)/thRich2;

    //            ***** COMBINED ID FLAGS *****
#  ifdef U3_USE_CHCHI2
    int chiok_pm = 0;
    if (chiok1==0x3) chiok_pm |= 0x1; if (chiok2==0x3) chiok_pm |= 0x2;
#  else
    int lhok_pm = 0;
    if (lhok1==0x3)  lhok_pm  |= 0x1; if (lhok2==0x3)  lhok_pm  |= 0x2;
#  endif

    // ********** RICH EFFICIENCY/PURITY HISTOs, LIKELIHOOD PROFILES **********

#  ifdef U3_REQUIRE_CLEANV0 // Was uspposed to allow to monitor S/B, as ``K0 selection becomes more and more demanding''
    bool v0ok = cleanV0;
#  else
    bool v0ok = true;
#  endif
    if (acc==0x3 &&                          // ***** BOTH Hs W/IN RICH... *****
	v0ok) {
#  ifdef U3_USE_CHCHI2
      hC_K02->Fill(m_pipi,chiok_pm);              // ...COMBINED Eff ESTIMATE #1
#  else
      hL_K02->Fill(m_pipi,lhok_pm);
#  endif
    }

    if (piThr<mom1 && mom1<60) pDomain |= 0x1;
    if (piThr<mom2 && mom2<60) pDomain |= 0x2;
    if (pDomain && mom1<PRICHCut) pRange |= 0x1;
    if (pDomain && mom2<PRICHCut) pRange |= 0x2;
    if (acc&pDomain) {

      // ******************************************************************
      // *************** AT LEAST 1 HADRON w/in RICH, RANGE ***************
      // ******************************************************************

      // ***** Phase Space Binning *****
      int iPi, iP1, iP2;
      for (iPi=iP1 = 0; iPi<npB+1; iPi++) if (mom1-piThr<pBins[iPi]) {
	iP1 = iPi-1; break;
      }
      for (iPi=iP2 = 0; iPi<npB+1; iPi++) if (mom2-piThr<pBins[iPi]) {
	iP2 = iPi-1; break;
      }
      //int iP1 = (int)(mom1-piThr)/10, iP2 =  (int)(mom2-piThr)/10;
      if (iP1<0) iP1 = 0; /* to be sure */ if (iP1>=npB) iP1 = npB-1;
      if (iP2<0) iP2 = 0; /* to be sure */ if (iP2>=npB) iP2 = npB-1;
#  ifdef PS_USE_RADIUS
      int iA1 = (int)(RR1/RRBin), iA2 = (int)(RR2/RRBin);
#  else
      int iA1 = (int)(AR1/ARBin), iA2 = (int)(AR2/ARBin);
#  endif
      if (iA1>=nK0Bins)  iA1 = nK0Bins-1;  if (iA2>=nK0Bins)  iA2 = nK0Bins-1;
      int iPA1 = 1+iP1+npB*iA1, iPA2 = 1+iP2+npB*iA2;
      // Binning for K misID
      int jP1, jP2;
      for (iPi=jP1 = 0; iPi<npB+1; iPi++) if (mom1-KThr<pBins[iPi]) {
	jP1 = iPi-1; break;
      }
      for (iPi=jP2 = 0; iPi<npB+1; iPi++) if (mom2-KThr<pBins[iPi]) {
	jP2 = iPi-1; break;
      }
      //int jP1 = (int)(mom1-KThr)/10, jP2 =  (int)(mom2-KThr)/10;
      if (jP1<0) jP1 = 0; /* to be sure */ if (jP1>2) jP1 = 2;
      if (jP2<0) jP2 = 0; /* to be sure */ if (jP2>2) jP2 = 2;
      int jA1 = iA1, jA2 = iA2; // w/ lesser extent 'cause few stats out there
      if (jA1>=nphiBins) jA1 = nphiBins-1; if (jA2>=nphiBins) jA2 = nphiBins-1;
      int jPA1 = 1+jP1+npB*jA1, jPA2 = 1+jP2+npB*jA2;
      time_t time = (time_t)e.UnixSeconds();
      int hour = localtime(&time)->tm_hour;
      int ihour = (hour+21)%24/6;                   // ***** AS A F(TIME of DAY)
      if (acc==0x3 &&            // ***** BOTH Hs w/in RICH, w/in RANGE... *****
	  v0ok) {
#  ifdef U3_USE_CHCHI2
	hC_K0M2->Fill(m_pipi,chiok_pm);           // ...COMBINED Eff ESTIMATE #2
#  else
	hL_K0M2->Fill(m_pipi,lhok_pm);
#  endif
      }
      // *************** chi2 Eff: OVERALL, f(PHASE SPACE,PD)   ***************
      // *************** chi2 Pur: OVERALL, f(PHASE SPACE,PD,T) ***************

#  ifdef U3_PLOT_CHi2
      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	int ok1 = 0, pv1 = 0;
	if (chiok1&0x1) ok1 = 1; if (chiok1&0x2) pv1 = 1;
	hC_K0oP[0][0] ->Fill(m_pipi,(double)ok1);
	hC_K0oPV[0][0]->Fill(m_pipi,(double)pv1);
	hC_K0oPT[0][0]->Fill(m_pipi,(double)(ok1&pv1));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hC_K0oP[0][iPA1] ->Fill(m_pipi,(double)ok1);
	hC_K0oPV[0][iPA1]->Fill(m_pipi,(double)pv1);
	hC_K0oPT[0][iPA1]->Fill(m_pipi,(double)(ok1&pv1));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hC_K0M[0]->Fill(m_pipi,(double)(ok1&pv1));
	  hC_K0PD[0][ipd1]->Fill(m_pipi,(double)(ok1&pv1));  // ***** AS A F(PD)
	}
	if (KThr<mom1 &&              // ***** ...AND w/in K RANGE... *****
	    lhok2==0x3 &&             // ***** ...AND pi-ID to get a pi+ beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  hR_K0oP[0][0]->Fill(RR1); hR_K0oP[0][jPA1]->Fill(RR1); 
	  hR_K0oA[0][0]->Fill(mom1-KThr); hR_K0oA[0][jPA1]->Fill(mom1-KThr); 
	  int ok1 = 0, pv1 = 0;       // ***** ...K misID Eff *****
	  if (chiKok1&0x1) ok1 = 1; if (chiKok1&0x2) pv1 = 1;
	  hD_K0oP[0][0] ->Fill(m_pipi,(double)ok1);
	  hD_K0oPV[0][0]->Fill(m_pipi,(double)pv1);
	  hD_K0oPT[0][0]->Fill(m_pipi,(double)(ok1&pv1));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hD_K0oP[0][jPA1] ->Fill(m_pipi,(double)ok1);
	  hD_K0oPV[0][jPA1]->Fill(m_pipi,(double)pv1);
	  hD_K0oPT[0][jPA1]->Fill(m_pipi,(double)(ok1&pv1));
	  if (pRange&0x1) {               // **********...w/in pRANGE **********
	    hD_K0M[0]->Fill(m_pipi,(double)(ok1&pv1));
	    hD_K0PD[0][ipd1]->Fill(m_pipi,(double)(ok1&pv1));// ***** AS A F(PD)
	    hD_K0H[0][ihour]->Fill(m_pipi,(double)(ok1&pv1));// AS A F(TIME o D)
	  }
	}
	else if (piThr*1.5<mom1<KThr &&   // ***** ...AND w/ sub KRANGE... *****
		 lhok2==0x3 && cleanV0)
	  hD_K0S[0]->Fill(m_pipi,(double)(chiSubKok1));
      }
      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	int ok2 = 0, pv2 = 0;
	if (chiok2&0x1) ok2 = 1; if (chiok2&0x2) pv2 = 1;
	hC_K0oP[1][0] ->Fill(m_pipi,(double)ok2);
	hC_K0oPV[1][0]->Fill(m_pipi,(double)pv2);
	hC_K0oPT[1][0]->Fill(m_pipi,(double)(ok2&pv2));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hC_K0oP[1][iPA2] ->Fill(m_pipi,(double)ok2);
	hC_K0oPV[1][iPA2]->Fill(m_pipi,(double)pv2);
	hC_K0oPT[1][iPA2]->Fill(m_pipi,(double)(ok2&pv2));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hC_K0M[1]->Fill(m_pipi,(double)(ok2&pv2));
	  hC_K0PD[1][ipd2]->Fill(m_pipi,(double)(ok2&pv2));  // ***** AS A F(PD)
	}
	if (KThr<mom2 &&              // ***** ...AND w/in K RANGE... *****
	    lhok1==0x3 &&             // ***** ...AND pi+ID to get a pi- beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  hR_K0oP[1][0]->Fill(RR2); hR_K0oP[1][jPA2]->Fill(RR2); 
	  hR_K0oA[1][0]->Fill(mom2-KThr); hR_K0oA[1][jPA2]->Fill(mom2-KThr); 
	  int ok2 = 0, pv2 = 0;       // ***** ...K misID Eff *****
	  if (chiKok2&0x1) ok2 = 1; if (chiKok2&0x2) pv2 = 1;
	  hD_K0oP[1][0] ->Fill(m_pipi,(double)ok2);
	  hD_K0oPV[1][0]->Fill(m_pipi,(double)pv2);
	  hD_K0oPT[1][0]->Fill(m_pipi,(double)(ok2&pv2));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hD_K0oP[1][jPA2] ->Fill(m_pipi,(double)ok2);
	  hD_K0oPV[1][jPA2]->Fill(m_pipi,(double)pv2);
	  hD_K0oPT[1][jPA2]->Fill(m_pipi,(double)(ok2&pv2));
	  if (pRange&0x2) {               // **********...w/in pRANGE **********
	    hD_K0M[1]->Fill(m_pipi,(double)(ok2&pv2));
	    hD_K0PD[1][ipd2]->Fill(m_pipi,(double)(ok2&pv2));// ***** AS A F(PD)
	    hD_K0H[1][ihour]->Fill(m_pipi,(double)(ok2&pv2));// AS A F(TIME o D)
	  }
	}
	else if (piThr*1.5<mom2<KThr &&   // ***** ...AND w/ sub KRANGE... *****
		 lhok1==0x3 && cleanV0)
	  hD_K0S[1]->Fill(m_pipi,(double)(chiSubKok2));
      }
#  endif
#  ifdef U3_PLOT_dTheta
      // ********** dTheta Perf: OVERALL, f(PHASE SPACE), f(PD) **********

      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	hT_K0oP[0][0] ->Fill(m_pipi,(double)(dthok1%2));
	hT_K0oPV[0][0]->Fill(m_pipi,(double)(dthok1/2));
	hT_K0oPT[0][0]->Fill(m_pipi,(double)(dthok1/3));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hT_K0oP[0][iPA1] ->Fill(m_pipi,(double)(dthok1%2));
	hT_K0oPV[0][iPA1]->Fill(m_pipi,(double)(dthok1/2));
	hT_K0oPT[0][iPA1]->Fill(m_pipi,(double)(dthok1/3));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hT_K0M[0]->Fill(m_pipi,(double)(dthok1/3));
	  hT_K0PD[0][ipd1]->Fill(m_pipi,(double)(dthok1/3)); // ***** AS A F(PD)
	}
	if (KThr<mom1 &&              // ***** ...AND w/in K RANGE... *****
	    lhok2==0x3 &&             // ***** ...AND pi-ID to get a pi+ beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  //                             ***** ...K misID Eff *****
	  hU_K0oP[0][0] ->Fill(m_pipi,(double)(dthKok1%2));
	  hU_K0oPV[0][0]->Fill(m_pipi,(double)(dthKok1/2));
	  hU_K0oPT[0][0]->Fill(m_pipi,(double)(dthKok1/3));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hU_K0oP[0][jPA1] ->Fill(m_pipi,(double)(dthKok1%2));
	  hU_K0oPV[0][jPA1]->Fill(m_pipi,(double)(dthKok1/2));
	  hU_K0oPT[0][jPA1]->Fill(m_pipi,(double)(dthKok1/3));
	  if (pRange&0x1)                 // **********...w/in pRANGE **********
	    hU_K0PD[0][ipd1]->Fill(m_pipi,(double)(dthKok1/3));// ***** AS A F(PD)
	}
      }
      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	hT_K0oP[1][0] ->Fill(m_pipi,(double)(dthok2%2));
	hT_K0oPV[1][0]->Fill(m_pipi,(double)(dthok2/2));
	hT_K0oPT[1][0]->Fill(m_pipi,(double)(dthok2/3));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hT_K0oP[1][iPA2] ->Fill(m_pipi,(double)(dthok2%2));
	hT_K0oPV[1][iPA2]->Fill(m_pipi,(double)(dthok2/2));
	hT_K0oPT[1][iPA2]->Fill(m_pipi,(double)(dthok2/3));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hT_K0M[1]->Fill(m_pipi,(double)(dthok2/3));
	  hT_K0PD[1][ipd2]->Fill(m_pipi,(double)(dthok2/3)); // ***** AS A F(PD)
	}
	if (KThr<mom2 &&              // ***** ...AND w/in K RANGE... *****
	    lhok1==0x3 &&             // ***** ...AND pi+ID to get a pi- beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  //                             ***** ...K misID Eff *****
	  hU_K0oP[1][0] ->Fill(m_pipi,(double)(dthKok2%2));
	  hU_K0oPV[1][0]->Fill(m_pipi,(double)(dthKok2/2));
	  hU_K0oPT[1][0]->Fill(m_pipi,(double)(dthKok2/3));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hU_K0oP[1][jPA2] ->Fill(m_pipi,(double)(dthKok2%2));
	  hU_K0oPV[1][jPA2]->Fill(m_pipi,(double)(dthKok2/2));
	  hU_K0oPT[1][jPA2]->Fill(m_pipi,(double)(dthKok2/3));
	  if (pRange&0x2)                 // **********...w/in pRANGE **********
	    hU_K0PD[1][ipd2]->Fill(m_pipi,(double)(dthKok2/3));// ***** AS A F(PD)
	}
      }
#  endif
      // ********** LIKELIHOODS **********


#  ifdef DEBUG_K0
      if (cleanV0) {
	if (dbK0p1p2==0x3) { dbK0Event |= 0x20; dbK0v1 = m_pipi; }
	if (.484672<=m_pipi && m_pipi<.509672) { // K0 +/- 2 sigma=6MeV
	  if (dbK0p1p2==0x3) dbK0Event = 0xff;
	  static FILE *dbK0Refs_o = NULL; // To compare 2 processings w/ Ediff
	  if (!dbK0Refs_o) dbK0Refs_o = fopen("Refs.out","w");
	  fprintf(dbK0_o,
 "%5d %3d %5d %10Ld 0xff  %.6f  %3d %6.3f %6.3f %6.3f  %3d %6.3f %6.3f %6.3f\n",
		    Run,e.SpillNum(),e.EvInSpill(),e.UniqueEvNum(),
		    m_pipi,
		  iET1,mom1,h1.Theta(),trk1.Chi2tot()/trk1.Ndf(),
		  iET2,mom2,h2.Theta(),trk2.Chi2tot()/trk2.Ndf());
	  fprintf(dbK0Refs_o,"%5d %2d %5d\n",Run,e.SpillNum(),e.EvInSpill());
	}
      }
#  endif
      int subThr = 0;
      if (piThr*thrMargin<mom1 && mom1<KThr) subThr |= 0x1;
      if (piThr*thrMargin<mom2 && mom2<KThr) subThr |= 0x2;
      if (acc&pDomain&0x1) {// ********** >0 w/in RICH w/in pDOMAIN...**********
	int ok1 = 0, pv1 = 0;
	if (lhok1&0x1) ok1 = 1; if (lhok1&0x2) pv1 = 1;
	hL_K0oP[0][0] ->Fill(m_pipi,(double)ok1);
	hL_K0oPV[0][0]->Fill(m_pipi,(double)pv1);
	hL_K0oPT[0][0]->Fill(m_pipi,(double)(ok1&pv1));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hL_K0oP[0][iPA1] ->Fill(m_pipi,(double)ok1);
	hL_K0oPV[0][iPA1]->Fill(m_pipi,(double)pv1);
	hL_K0oPT[0][iPA1]->Fill(m_pipi,(double)(ok1&pv1));
	if (pRange&0x1) {                 // **********...w/in pRANGE **********
	  hL_K0M[0]->Fill(m_pipi,(double)(ok1&pv1));        // pi+ Likelihood ID
	  hL_K0PD[0][ipd1]->Fill(m_pipi,(double)(ok1&pv1));  // ***** AS A F(PD)
	}
	if (KThr<mom1 &&              // ***** ...AND w/in K RANGE... *****
	    lhok2==0x3 &&             // ***** ...AND pi-ID to get a pi+ beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  //                             ***** ...K misID Eff *****
	  int ok1 = 0, pv1 = 0;       // ***** ...K misID Eff *****
	  if (lhKok1&0x1) ok1 = 1; if (lhKok1&0x2) pv1 = 1;
	  hM_K0oP[0][0] ->Fill(m_pipi,(double)ok1);
	  hM_K0oPV[0][0]->Fill(m_pipi,(double)pv1);
	  hM_K0oPT[0][0]->Fill(m_pipi,(double)(ok1&pv1));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hM_K0oP[0][jPA1] ->Fill(m_pipi,(double)ok1);
	  hM_K0oPV[0][jPA1]->Fill(m_pipi,(double)pv1);
	  hM_K0oPT[0][jPA1]->Fill(m_pipi,(double)(ok1&pv1));
	  if (pRange&0x1) {               // **********...w/in pRANGE **********
	    hM_K0M[0]->Fill(m_pipi,(double)(ok1&pv1));
	    hM_K0PD[0][ipd1]->Fill(m_pipi,(double)(ok1&pv1));// ***** AS A F(PD)
	    hM_K0H[0][ihour]->Fill(m_pipi,(double)(ok1&pv1));// AS A F(TIME o D)
	    if (PIDs[iET1]&0x10)   // Likelihood available => Profile
	      hR_LHpi[0]->Fill(m_pipi,richLHs[0][iET1]);
	  }
	}
	else if ((subThr&0x1) &&          // ***** ...AND w/ sub KRANGE... *****
		 lhok2==0x3 && cleanV0)
	  hM_K0S[0]->Fill(m_pipi,(double)(lhSubKok1));
      }
      if (acc&pDomain&0x2) {// ********** <0 w/in RICH w/in pDOMAIN...**********
	int ok2 = 0, pv2 = 0;
	if (lhok2&0x1) ok2 = 1; if (lhok2&0x2) pv2 = 1;
	hL_K0oP[1][0] ->Fill(m_pipi,(double)ok2);
	hL_K0oPV[1][0]->Fill(m_pipi,(double)pv2);
	hL_K0oPT[1][0]->Fill(m_pipi,(double)(ok2&pv2));
	//                                ***** w/,w/o VETO, AS A F(PHASE SPACE)
	hL_K0oP[1][iPA2] ->Fill(m_pipi,(double)ok2);
	hL_K0oPV[1][iPA2]->Fill(m_pipi,(double)pv2);
	hL_K0oPT[1][iPA2]->Fill(m_pipi,(double)(ok2&pv2));
	if (pRange&0x2) {                 // **********...w/in pRANGE **********
	  hL_K0M[1]->Fill(m_pipi,(double)(ok2&pv2));
	  hL_K0PD[1][ipd2]->Fill(m_pipi,(double)(ok2&pv2));  // ***** AS A F(PD)
	}
	if (KThr<mom2 &&              // ***** ...AND w/in K RANGE... *****
	    lhok1==0x3 &&             // ***** ...AND pi+ID to get a pi- beam
	    cleanV0) {   // ***** ... AND cleanV0: 'cause p>KThr gets it dirtier
	  //                             ***** ...K misID Eff *****
	  int ok2 = 0, pv2 = 0;       // ***** ...K misID Eff *****
	  if (lhKok2&0x1) ok2 = 1; if (lhKok2&0x2) pv2 = 1;
	  hM_K0oP[1][0] ->Fill(m_pipi,(double)ok2);
	  hM_K0oPV[1][0]->Fill(m_pipi,(double)pv2);
	  hM_K0oPT[1][0]->Fill(m_pipi,(double)(ok2&pv2));
	  //                              ***** w/,w/o VETO, AS A F(PHASE SPACE)
	  hM_K0oP[1][jPA2] ->Fill(m_pipi,(double)ok2);
	  hM_K0oPV[1][jPA2]->Fill(m_pipi,(double)pv2);
	  hM_K0oPT[1][jPA2]->Fill(m_pipi,(double)(ok2&pv2));
	  if (pRange&0x2) {               // **********...w/in pRANGE **********
	    hM_K0M[1]->Fill(m_pipi,(double)(ok2&pv2));
	    hM_K0PD[1][ipd2]->Fill(m_pipi,(double)(ok2&pv2));// ***** AS A F(PD)
	    hM_K0H[1][ihour]->Fill(m_pipi,(double)(ok2&pv2));// AS A F(TIME o D)
	    if (PIDs[iET2]&0x10)   // Likelihood available => Profile
	      hR_LHpi[1]->Fill(m_pipi,richLHs[0][iET2]);
	  }
	}
	else if ((subThr&0x2) &&          // ***** ...AND w/ sub KRANGE... *****
		 lhok1==0x3 && cleanV0)
	  hM_K0S[1]->Fill(m_pipi,(double)(lhSubKok2));
      }

      if (isK0) {               // *************** pipi = K0 ***************
	// ********** ...RICH CHARACTERISTICS: Resolution, etc... **********
	// id not asked for: 'cause include cut on beta<1
	if ((acc&pRange&0x1) && (thid&0x1) &&          // IF theta AVAILABLE....
	    thMom1 &&                                  // ...>piThr+margin
#  ifdef piS_e_REJECTION
	    lhpiSok1 &&                                // ...piS e-Rejection
#  endif
	    lhok2==0x3) {                              // ...IF pi-ID...
	  // ***** RESIDUALS...
	  hR_dpiPD[0][jpd1]->Fill(thRich1-thMom1);// ...f(pseudoPhotonDetector)
	  hR_dvsr[0]->Fill(Run,thRich1-thMom1);   // ...vs. Run
	  int iph = phD1<.3?0:1;// Resolution vs. p-Thr as a f(phDensity,+/-)...
	  int ia = fabs(YR1)<10?0:1;                 // ...and angle
	  hR_dpivst[0][iph][ia]->Fill(thMom1,thRich1-thMom1);
	  // ***** ...ELSE
	  hR_nphvsp[0]->Fill(mom1-piThr,phD1);       // # of photons
	  hR_chivsp[0]->Fill(mom1-piThr,trk1.RichInf(12));  // chi2 vs. p
#  if defined UserEvent_NTuple && ( UserEvent_NTuple & 1 )
	  if (cleanV0) {
	    double ax1 = h1.Mom3().X()/h1.Mom3().Z();
	    double ay1 = h1.Mom3().Y()/h1.Mom3().Z();
	    ntDth->Fill(thRich1-thMom1,mom1-piThr,
			XR1,YR1,ax1,ay1,phD1,1);
	  }
#  endif
	}
	if ((acc&pRange&0x2) && (thid&0x2) &&          // IF theta AVAILABLE....
	    thMom2 &&                                  // ...>piThr+margin
#  ifdef piS_e_REJECTION
	    lhpiSok2 &&                                // ...piS e-Rejection
#  endif
	    lhok1==0x3) {                              // ...IF pi+ID...
	  // ***** RESIDUALS...
	  hR_dpiPD[1][jpd2]->Fill(thRich2-thMom2);// ...f(pseudoPhotonDetector)
	  hR_dvsr[1]->Fill(Run,thRich2-thMom2);   // ...vs. Run
	  int iph = phD2<.3?0:1;// Resolution vs. p-Thr as a f(phDensity,+/-)...
	  int ia = fabs(YR2)<10?0:1;                 // ...and angle
	  hR_dpivst[1][iph][ia]->Fill(thMom2,thRich2-thMom2);
	  // ***** ...ELSE
	  hR_nphvsp[1]->Fill(mom2-piThr,phD2);       // # of photons
	  hR_chivsp[1]->Fill(mom2-piThr,trk2.RichInf(12));  // chi2 vs. p
#  if defined UserEvent_NTuple && ( UserEvent_NTuple & 1 )
	  if (cleanV0) {
	    double ax2 = h2.Mom3().X()/h2.Mom3().Z();
	    double ay2 = h2.Mom3().Y()/h2.Mom3().Z();
	    ntDth->Fill(thRich2-thMom2,mom2-piThr,
			XR2,YR2,ax2,ay2,phD2,-1);
	  }
#  endif
	}
      }
#  define U3_dTheta_K0_SELECTION
#  ifdef U3_dTheta_K0_SELECTION
      static TH2D *hR_K0 = 0;
      if (!hR_K0)
	hR_K0 = new TH2D("hR_K0","#pi+#pi- - #delta#theta selection",
			 150,M_K0-0.075,M_K0+0.075,2,-.5,1.5);
      if ((acc&pRange&0x1) && (thid&0x1) &&
#    ifdef piS_e_REJECTION
	  lhpiSok1 &&
#    endif
	  lhok2==0x3) hR_K0->Fill(m_pipi,0);
      if ((acc&pRange&0x2) && (thid&0x2) &&
#    ifdef piS_e_REJECTION
	  lhpiSok2 &&
#    endif
	  lhok1==0x3) hR_K0->Fill(m_pipi,1);
#  endif

      if (RR1<rRICHCut || (pRange&0x1)==0) acc &= 0x2;  // Radius cut
      if (RR2<rRICHCut || (pRange&0x2)==0) acc &= 0x1;
      if (acc==0x3 &&         // ***** BOTH Hs w/in RICH, RANGE, RADIUS... *****
	  v0ok) {
#  ifdef U3_USE_CHCHI2
	hC_K0MA2->Fill(m_pipi,chiok_pm);          // ...COMBINED Eff ESTIMATE #3
#  else
	hL_K0MA2->Fill(m_pipi,lhok_pm);
#  endif
      }
    }  // End one hadron w/in momentum range
  }  // End one hadron w/in RICH 
}
// **********************************************************************
// ********************      bookLambda_RICHPerf     ********************
// **********************************************************************
void bookLambda_RICHPerf(double dM_Lambda, double LMassCut)
{
  // ******************** Lambda * RICH PERFS ********************
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dLambda; if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda")))
    gDirectory->mkdir("Lambda","Lambda Histos");
  if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda"))) {
    printf("\n** U3:\a No creating subdir \"Lambda\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dLambda->cd();
  TDirectory *dRICH; if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs")))
    gDirectory->mkdir("RICHPerfs","RICH Performances from Lambda decay");
  if (!(dRICH = (TDirectory*)gDirectory->Get("RICHPerfs"))) {
    printf("\n** U3:\a No creating subdir \"RICHPerfs\" in TDirectory \"%s/Lambda\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dRICH->cd();

#  ifdef U3_USE_CHCHI2
  char cChi2Cut[] = "p#chi^{2}<5.0,<0.95#pi,K  ";
  sprintf(cChi2Cut,"p#chi^{2}<%.1f,<#pi,K",CHChi2Mx);
  string sChi2piCut(cChi2Cut);
  // Cuts w/ pi veto
  sprintf(cChi2Cut,"p#chi^{2}<%.1f,<%.2f#pi,K",CHChi2Mx,CHChi2Cut);
  string sChi2Cut(cChi2Cut);
#  endif
  char cLHCut[] = "p>1.001bck,>0.999#piKe,P>1.05pThr  ";
  int precLHBC = LHBckCut-1 ? int(-log10(fabs(LHBckCut-1))+.999) : 0;
  int precLHC  = LHCut-1    ? int(-log10(fabs(LHCut   -1))+.999) : 0;
  sprintf(cLHCut,"p>%#.*fbck,>%#.*f#piKe,P>%.2fpThr",precLHBC,LHBckCut,precLHC,LHCut,pIDPmin);
  string sLHCut(cLHCut);
  char cLHVeto[] = ",#piKeVeto(1.50,1.40,>1.10#piThr)  ";
  sprintf(cLHVeto,",#piKeVeto(%.2f,%.2f)",subpThrLHpiVeto,subpThrLHVeto);
  string sLHVeto(cLHVeto);
#  ifdef piS_e_REJECTION
  char cLHeReject[] = ",#pi:eRejection(1.5)  ";
  sprintf(cLHeReject,",#pi:eRejection(%.1f)",piLHeVeto); string sLHeReject(cLHeReject);
#  endif

#  ifdef U3_PLOT_dTheta
  char cdTpCut[] = "#delta_{p}#theta<2.5#sigma  ";
  sprintf(cdTpCut,"#delta_{p}#theta<%.1f#sigma",nSigmas);
  char cKVCut[] = "#delta_{K}#theta>2.5#sigma  ";
  sprintf(cKVCut,"#delta_{K}#theta>%.1f#sigma",nVeto);
  string sdTpCut(cdTpCut), sKVCut(cKVCut);
  string sdTCut(cdTpCut+string("#times")+cKVCut);
#  endif
#  ifdef RICH_ANGLE_CUT
  char cRCut[] = ",A>25"; sprintf(cRCut,",R>%.0f",1000*aRICHCut);
#  else    // ***** Radius cut *****
  char cRCut[] = ",R>12"; sprintf(cRCut,",R>%.0f",rRICHCut);
#  endif
  string sRCut(cRCut);

  string ppi("p#pi- - "), pip("#bar{p}#pi+ - ");
  char cPRange[] = ",1.10#piThr<P<100GeV";
  sprintf(cPRange,",%.2f#piThr<P<%.0fGeV",subpThrPmin,pIDPmax);
  string sPRange(cPRange), spID("pID:");

  char hName[] = "hR_dpPD00";
  int nBinsX = int((M_Lam+dM_Lambda-M_p-M_pi+.005)/.0005+.5);
  double xMn = M_Lam+dM_Lambda-nBinsX*.0005;
#  if U3_USE_CHCHI2
  //           *************** Chi2 ***************
  // piD...
  // ...Raw efficiency
  // ...Momenta over threshold
  hC_L   = new TH2D("hC_L",  string(ppi+spID+sChi2Cut).c_str(),
		      nBinsX,xMn,M_Lam+dM_Lambda,2,-.5,1.5);
  hC_LM  = new TH2D("hC_LM", string(ppi+spID+sChi2Cut+sPRange).c_str(),
		      nBinsX,xMn,M_Lam+dM_Lambda,2,-.5,1.5);
  // ...Momenta over threshold AND radius(angle) < limit
  hC_LMA = new TH2D("hC_LMA",string(ppi+spID+sChi2Cut+sPRange+sRCut).c_str(),
		      nBinsX,xMn,M_Lam+dM_Lambda,2,-.5,1.5);
#  else
  //        *************** Likelihoods  ***************  
  const char *LaL[] = { "L", "aL" }; for (int iLaL = 0; iLaL<2; iLaL++) {
    sprintf(hName,"hL_%s",LaL[iLaL]);         // For I) Lambda, I)) anti-Lambda
    string &sp = iLaL ? pip : ppi;
#    ifdef piS_e_REJECTION
    hL_L[iLaL]   = new TH2D(hName,
  string(sp+spID+sLHCut+sLHVeto+sLHeReject).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,6,-1.5,4.5);
    sprintf(hName,"hL_%sM",LaL[iLaL]);        // ...Momenta over threshold
    hL_LM[iLaL]  = new TH2D(hName,
  string(sp+spID+sLHCut+sLHVeto+sLHeReject+sPRange).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,6,-1.5,4.5);
    sprintf(hName,"hL_%sMA",LaL[iLaL]);       // ...AND radius(angle) < limit
    hL_LMA[iLaL] = new TH2D(hName,
  string(sp+spID+sLHCut+sLHVeto+sLHeReject+sPRange+sRCut).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,6,-1.5,4.5);
#    else
    hL_L[iLaL]   = new TH2D(hName,
      string(sp+spID+sLHCut+sLHVeto).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,4,-1.5,2.5);
    sprintf(hName,"hL_%sM",LaL[iLaL]);         // ...Momenta over threshold
    hL_LM[iLaL]  = new TH2D(hName,
      string(sp+spID+sLHCut+sLHVeto+sPRange).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,4,-1.5,2.5);
    sprintf(hName,"hL_%sMA",LaL[iLaL]);        // ...AND radius(angle) < limit
    hL_LMA[iLaL] = new TH2D(hName,
      string(sp+spID+sLHCut+sLHVeto+sPRange+sRCut).c_str(),
			     nBinsX,xMn,M_Lam+dM_Lambda,4,-1.5,2.5);
#    endif
  }
#  endif

#  ifdef U3_PLOT_dTheta
  // *************** dTheta: pID Eff FOR MOMENTA OVER THRESHOLD *****
  hT_LM  = new TH2D("hT_LM", string(ppi+spID+sdTCut+sPRange).c_str(),
		     nBinsX,xMn,M_Lam+dM_Lambda,2,-.5,1.5);
#  endif

  char hTitle[]   = "RICH p#delta#theta -  #Lambda#pm5.2Mev,PD0  ";
  //char hTitle[] = "RICH #theta vs. Pp - #Lambda#pm20Mev ";
  for (int ipd = 0; ipd<4; ipd++) {
    // *************** CERENKOV ANGLE RESOLUTION ***************
    sprintf(hName,"hR_dpPD%d",ipd);
    sprintf(hTitle,"RICH p#delta#theta - #Lambda#pm%.1fMev,PD%d",
	    LMassCut*1000,ipd);
    hR_dpPD[ipd] = new TH1D(hName,hTitle,100,-10,10);
  }

  // *************** LIKELIHOOD PROFILES ***************
  // w/in RICH and P range
  hR_LHp   = new TProfile("hR_LHp","RICH pLH - pThr<p<45",
			  nBinsX,xMn,M_Lam+dM_Lambda);
  hR_LHppi = new TProfile("hR_LHppi","RICH p/#piLH - pThr<p<45",
			  nBinsX,xMn,M_Lam+dM_Lambda);
}
// ********************************************************************
// ********************     fillLambda_RICHPerf    ********************
// ********************************************************************
void fillLambda_RICHPerf(PaEvent &e, int iET1, int iET2, int iETp,
			 const TLorentzVector &lvp1,
			 double m_ppi, bool &isLambda, int Run)
{
  // ********** DETERMINE TRACKS' ACCEPTANCE by RICH **********

  const PaTrack &trkp = e.vTrack()[iETp];
#  ifdef piS_e_REJECTION
  int iETpi, iLaL; if (iETp==iET1) { iETpi = iET2; iLaL = 0; }
  else                             { iETpi = iET1; iLaL = 1; }
#  else
  int        iLaL; if (iETp==iET1)                 iLaL = 0;
  else                                             iLaL = 1;
#  endif
  int acc = 0; // Acceptance Flag...
  if (aRs[iETp]==0) {
 printf("** U3/Lambda RICHPerf:\a Evt %d,%d Track %d RICH not yet checked\n\n",
	e.RunNum(),(int)e.UniqueEvNum(),iETp); assert(false);
  }
  if (aRs[iETp]>0) acc |= 0x1;
  const vector<Float_t> &aux = trkp.vAux();
  if (aRs[iETp]<-1) { // Ensure then "aRs" are meaningful
    printf("** U3:\a Evt %d,%d Track %d(<-Lambda): Inconsistency: aRs = %d\n",
	   e.RunNum(),(int)e.UniqueEvNum(),iETp,aRs[iETp]);
    abort();
  }
  double XR = aux[0], YR = aux[1], RR = sqrt(XR*XR+YR*YR);
  // Start hlx
  const PaTPar &hp = trkp.vTPar()[0];
  // More precision: which part of acceptance is hit...
  int ipdp  = 0;  // .. PD = pseudo Photon Detector
  if (XR>0) ipdp |= 0x1; if (YR>0) ipdp |= 0x2;
  double momp = hp.Mom();
  int subDomain = 0;

  if (acc) {    // ********** p(anti-p) W/IN RICH ACCEPTANCE **********


    if (piThr<momp &&  // SubThreshold domain >piThr. Alternative could be >KThr
	momp<pIDPmin*pThr) subDomain = 1;
    // ***** SET A NUMBER OF RICH-RELATED VARIABLES *****
    int thid = 0;              // "th" = availability of theta
    if (PIDs[iETp]&0x8) thid |= 0x1;
    double thRich = 0, thMom = 0;
#  ifdef U3_PLOT_dTheta    // deltaTheta pID: Retain only 2 LS bits
    unsigned int dthok = rich->pID(trkp,thRich,thMom)&0x3;
#  else
    if (pThr+.001<momp)
      thMom = acos(sqrt(M2_p/momp/momp+1)/rich->CurrentIndAPV)*1000;
    if (PIDs[iETp]&0x8) thRich = trkp.RichInf(RICH_ESTIMATOR);
#  endif
#  ifdef U3_USE_CHCHI2
    // ***** pID: chi2 method *****
    int chiok = 0;
    if (PIDs[iETp]&0x10) {
      double pichi2 = richChi2s[0][iETp], Kchi2  = richChi2s[1][iETp];
      double pchi2  = richChi2s[2][iETp];
      if (pchi2<          Kchi2 && pchi2<          pichi2 &&
	  pchi2<CHChi2Mx)                                      chiok |= 0x1;
      if (pchi2<CHChi2Cut*Kchi2 && pchi2<CHChi2Cut*pichi2)     chiok |= 0x2;
    }
#  endif
    // ***** p(pi)ID: LH method *****
    int lhok = 0, lhSubok = 0;
#  ifdef piS_e_REJECTION
    int lhpiok = 1;
#  endif
    if (PIDs[iETp]&0x10) {                     // Likelihood for p (or anti-p)
      double piLH = richLHs[0][iETp], KLH = richLHs[1][iETp];
      double pLH  = richLHs[2][iETp], eLH = richLHs[3][iETp];
      if (!subDomain && pLH>LHCut*KLH && pLH>LHCut*piLH && pLH>LHCut*eLH &&
	  pLH>LHBckCut)                                    lhok = 1;
      if (subDomain && piLH<subpThrLHpiVeto && (momp<KThr || KLH<subpThrLHVeto))  {
#    ifdef REJECT_SubKThr_e
	bool eid = eLH>subpThrLHVeto;
#    else
	bool eid = false;
#    endif
	if (!eid)                                          lhSubok = 1;
      }
    }
#  if K_BELOW_THR > 1
    else if (acc&0x1) {
      // SubThreshold pID: no tighter requirement on momentum than >piThr here,
      // but later on the "lhSubKok" will be enabled only in a resticted range.
      if (subDomain)                                       lhSubok = 1;
    }
#  endif
#  ifdef piS_e_REJECTION
    if (PIDs[iETpi]&0x10) {                    // Likelihood for pi
      double piLH = richLHs[0][iETpi], eLH = richLHs[3][iETpi];
      if (eLH>piLHeVeto*piLH)                              lhpiok = 0;
    }
#  endif

    // ********** RICH EFFICIENCY HISTOs and LIKELIHOOD PROFILES **********

#  ifdef U3_USE_CHCHI2                                         // Eff ESTIMATE #1
    if (!iLaL) hC_L->Fill(m_ppi,(double)(chiok==0x3?1:0));
#  else
    int idok; if (momp<pThr*pIDPmin) idok = lhSubok ? 1 : -1;
    else                             idok = lhok    ? 2 :  0;
#    ifdef piS_e_REJECTION
    if (idok>0 && lhpiok) idok += 0x2;
#    endif
    hL_L[iLaL]->Fill(m_ppi,(double)idok);
#  endif

    if (piThr*subpThrPmin<momp && momp<pIDPmax) {
      // ***** p w/in RICH, pRANGE (high upper bound (=pIDPmax) for pID) *****

#  ifdef U3_USE_CHCHI2                                         // Eff ESTIMATE #2
      if (!iLaL) hC_LpM->Fill(m_ppi,(double)(chiok==0x3?1:0));
#  else
      hL_LM[iLaL]->Fill(m_ppi,idok);
#  endif
#  ifdef U3_PLOT_dTheta
      if (!iLaL) hT_LpM->Fill(m_ppi,(double)(dthok1/3));
#  endif

      if (pThr<momp && PIDs[iETp]&0x10) {             // Likelihood for p/anti-p
	double pLH = richLHs[2][iETp]; double piLH = richLHs[0][iETp];
	double ppiLH = pLH/piLH;
	hR_LHp->Fill(m_ppi,pLH);
	if (PIDs[iETp]&0x20)              // LH divided by LH for pi
	  hR_LHppi->Fill(m_ppi,ppiLH);
      }

      if (isLambda) {                               // ***** ppi == Lambda *****
	if (thid&0x1) {                                // if theta available
	  hR_dpPD[ipdp]->Fill(thRich-thMom);  // Resol. as a f(pseudoPD)
	}
      }

      if (RR>rRICHCut) {                     // ***** RADIUS CUT...
#  ifdef U3_USE_CHCHI2                                        // Eff ESTIMATE #3
	hC_LMA->Fill(m_ppi,(double)(chiok==0x3?1:0));
#  else
	hL_LMA[iLaL]->Fill(m_ppi,idok);
#  endif
      }
#  ifndef U3_USE_CHCHI2
      if (idok<=0) isLambda = 0;    // ***** UPDATE "isLambda" W/ PID *****
#    ifdef piS_e_REJECTION
      if (!lhpiok) isLambda = 0;
#    endif
#  endif
    }  // End p w/in p Range
  }  // End p (>0 particle in fact) w/in RICH acceptance
}

#  ifdef DEBUG_K0
void dbK0GetEvent(const PaEvent &e)
{
  // Search input "K0s.in" file for next event where a K0 is expected (from
  // what was determined by the processing that produced the input file)

  // Sets "dbK0Event" flag and "dbK0" track characteristics if found.

  if (dbK0Event==-3) return;          // ==-3 => No input available
  static FILE *dbK0_i = 0;
  if (dbK0Event==-2) {                 // ==0 Means: Init
    dbK0_o  = fopen("K0s.out","w");          // ...=> Open output
    if (!(dbK0_i = fopen("K0s.in","r"))) {   // ...=> Try and open input
      dbK0Event = -3;                 // =-3: No input ever
      printf("\n * DEBUG_K0:\a No input list of K0s \"%s\"\n","K0s.in"); return;
    }
  }

  dbK0Event = -1;    // Default is no match

  long long curEvent = e.UniqueEvNum();
  static int curSpill = 998; int prvSpill = curSpill; curSpill = e.SpillNum();
  static long long fileEvent = -1; static int fileSpill = 999;
  if (fileEvent!=-1) {
    if      (curEvent==fileEvent) {
      dbK0Event = 0; return;    // Input file event reached: set "dbK0Event"
    }
    else if (curSpill>=prvSpill) {
      if (curEvent<fileEvent)
	return;  // Input file event not yet reached: return
    }
    //...else We have (probably) siwtched to a new chunk...
  }

  // Search the file
  char line[121]; static int iline = 0; char *ioError; int nItems, nChars;
  int run, evInSp, iET; unsigned int score; float m_pipi, chi2;
  while (fileEvent<curEvent && fileSpill>=prvSpill) {

    while ((ioError = fgets(line,120,dbK0_i))>0 &&
	   (nItems = sscanf(line,"%d %d %d %Ld 0x%x %n",
			    &run,&fileSpill,&evInSp,&fileEvent,
			    &score,&nChars))==5 && score!=0xff) iline++;
    if (nItems<5) {
      printf("\n** U3:\a DEBUG_K0: Error syntax \"%s\" line %d\n%s\n\n",
	     "K0s.in",iline,line); assert(false);
    }
    if (ioError==0) {         // EoF 
      printf("\n** U3:\a DEBUG_K0: EoF \"%s\"\n\n","K0s.in");
      dbK0Event = -3; return;                  // => =-3: No more any input
    }
    if ((nItems = sscanf(line+nChars,"%f %d %f %f %f %d %f %f %f",&m_pipi,
			 &iET,dbK0p+0,dbK0a+0,&chi2,
			 &iET,dbK0p+1,dbK0a+1,&chi2))!=9) {
      printf("\n** U3:\a DEBUG_K0: Error syntax \"%s\" line %d\n%s\n\n",
	   "K0s.in",iline,line); assert(false);
    }
  }
  if (fileEvent==curEvent)
    dbK0Event = 0;                 // Input file event reached: set "dbK0Event"
}
#  endif
#endif  // End "U3_FILL_RICH_PERFS
// **********************************************************************
// ****************************** bookK0X *******************************
// **********************************************************************
void bookK0X(double dM_phi, double K0MassCut, double KSMassCut)
{
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dK0; if (!(dK0 = (TDirectory*)gDirectory->Get("K0")))
		     gDirectory->mkdir("K0","K0 Histos");
  if (!(dK0 = (TDirectory*)gDirectory->Get("K0"))) {
    printf("\n** U3:\a No creating subdir \"K0\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dK0->cd();
  TDirectory *dK0X; if (!(dK0X = (TDirectory*)gDirectory->Get("K0X")))
		      gDirectory->mkdir("K0X","D+/-->K0pi, D0->K*pi,K0phi");
  if (!(dK0X = (TDirectory*)gDirectory->Get("K0X"))) {
    printf("\n** U3:\a No creating subdir \"K0X\" in TDirectory \"%s/K0\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dK0X->cd();

  //           ********** D+/- -> K0pi **********
  int nBinsX = int((M_D+.4-M_KS+.2)/.005+.5); double xMn = M_D+.4-nBinsX*.005;
  hm_K0pi  = new TH1D("hm_K0pi",  "K0#pi",            nBinsX,xMn,M_D+.4);
  hm_K0piC = new TH1D("hm_K0piC", "K0#pi - KinCuts",nBinsX,xMn,M_D+.4);

  //            ********** D0 -> K* pi **********
  char KSpiName[] = "hm_KSppiC";
  char KSpiTitle[] = "K*+#pi- - -2.8<D0#pi+-D*<3.2MeV,#pi+#pi-=K0#pm15.2MeV,K0#pi+=K*#pm100MeV,#piID,zD>0.25,#SigmapT>1.6GeV  ";
  for (int pm = 0; pm<2; pm++) {
    char sK = pm ? '-' : '+', spi = pm ? '+' : '-';
    sprintf(KSpiName,"hm_KS%cpi",pm?'m':'p');
    sprintf(KSpiTitle,"K*%c#pi%c - #pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID",sK,spi,K0MassCut*1000,sK,KSMassCut*1000);
    hm_KSpi[pm] = new TH1D(KSpiName,KSpiTitle,160,M_D0-.4,M_D0+.4);
    sprintf(KSpiName,"hm_KS%cpiC",pm?'m':'p');
    sprintf(KSpiTitle,"K*%c#pi%c - #pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,KinCuts",sK,spi,K0MassCut*1000,sK,KSMassCut*1000);
    hm_KSpiC[pm] = new TH1D(KSpiName,KSpiTitle,160,M_D0-.4,M_D0+.4);
    sprintf(KSpiName,"hm_KS%cpit",pm?'m':'p');
    sprintf(KSpiTitle,"K*%c#pi%c - #pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,zD>%.2f,#SigmapT>1.6GeV",sK,spi,K0MassCut*1000,sK,KSMassCut*1000,zCut0);
    hm_KSpit[pm] = new TH1D(KSpiName,KSpiTitle,160,M_D0-.4,M_D0+.4);
    sprintf(KSpiName,"hm_KS%cpiT",pm?'m':'p');
    sprintf(KSpiTitle,"K*%c#pi%c - #pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,zD>%.2f,#SigmapT>2GeV",sK,spi,K0MassCut*1000,sK,KSMassCut*1000,zCut0);
    hm_KSpiT[pm] = new TH1D(KSpiName,KSpiTitle,160,M_D0-.4,M_D0+.4);
    for (int wr = 0; wr<2; wr++) {
      char spiwr = wr ? sK : spi; int pmwr = wr ? pm+2 : pm;
      sprintf(KSpiName,"h%c_KS%cpi",wr?'W':'S',pm?'m':'p');
      sprintf(KSpiTitle,"K*%c#pi%c - %.1f<D0#pi%c-D*<%.1fMeV,#pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID",sK,spi,D0piLow*1000,spiwr,D0piUp*1000,K0MassCut*1000,sK,KSMassCut*1000);
      hS_KSpi[pmwr] = new TH1D(KSpiName,KSpiTitle,240,M_D0-.6,M_D0+.6);
      sprintf(KSpiName,"h%c_KS%cpiC",wr?'W':'S',pm?'m':'p');
      sprintf(KSpiTitle,"K*%c#pi%c - %.1f<D0#pi%c-D*<%.1fMeV,#pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,KinCuts",sK,spi,D0piLow*1000,spiwr,D0piUp*1000,K0MassCut*1000,sK,KSMassCut*1000);
      hS_KSpiC[pmwr] = new TH1D(KSpiName,KSpiTitle,240,M_D0-.6,M_D0+.6);
      sprintf(KSpiName,"h%c_KS%cpit",wr?'W':'S',pm?'m':'p');
      sprintf(KSpiTitle,"K*%c#pi%c - %.1f<D0#pi%c-D*<%.1fMeV,#pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,zD>%.2f,#SigmapT>1.6GeV",sK,spi,D0piLow*1000,spiwr,D0piUp*1000,K0MassCut*1000,sK,KSMassCut*1000,zCutS);
      hS_KSpit[pmwr] = new TH1D(KSpiName,KSpiTitle,240,M_D0-.6,M_D0+.6);
      sprintf(KSpiName,"h%c_KS%cpiT",wr?'W':'S',pm?'m':'p');
      sprintf(KSpiTitle,"K*%c#pi%c - %.1f<D0#pi%c-D*<%.1fMeV,#pi+#pi-=K0#pm%.1fMeV,K0#pi%c=K*#pm%.0fMeV,#piID,zD>%.2f,#SigmapT>2GeV",sK,spi,D0piLow*1000,spiwr,D0piUp*1000,K0MassCut*1000,sK,KSMassCut*1000,zCutS);
      hS_KSpiT[pmwr] = new TH1D(KSpiName,KSpiTitle,240,M_D0-.6,M_D0+.6);
    }
  }
  char KSpipiName[] = "hm_KSpipiC";
  char KSpipiTitle[] = "D0(K*+#pi-)#pi- - K*#pi=D0#pm30MeV,#piID,KinCuts  ";
  for (int pm = 0; pm<2; pm++) for (int wr = 0; wr<2; wr++) {
    char sK = pm ? '-' : '+', spi = pm ? '+' : '-', spiwr = wr ? sK : spi;
    double xMn = M_D0+M_pi-.01, xMx = M_D0+M_pi+.05;
    int pmwr = pm +2*wr;
    sprintf(KSpipiName,"h%c_KS%cpipi",wr?'w':'m',pm?'m':'p');
    sprintf(KSpipiTitle,"D0(K*%c#pi%c)#pi%c - K*#pi=D0#pm30MeV,#piID",sK,spi,spiwr);
    hm_KSpipi[pmwr] = new TH1D(KSpipiName,KSpipiTitle,120,xMn,xMx);
    sprintf(KSpipiName,"h%c_KS%cpipiC",wr?'w':'m',pm?'m':'p');
    sprintf(KSpipiTitle,"D0(K*%c#pi%c)#pi%c - K*#pi=D0#pm30MeV,#piID,KinCuts",sK,spi,spiwr);
    hm_KSpipiC[pmwr] = new TH1D(KSpipiName,KSpipiTitle,120,xMn,xMx);
  }

  //             ********** D0->K0phi **********
  hm_K0phi  = new TH1D("hm_K0phi",  "K0#phi",160,M_D0-.4,M_D0+.4);
  hm_K0phiC = new TH1D("hm_K0phiC", "K0#phi - KinCuts",
		       160,M_D0-.4,M_D0+.4);
  hm_K0phiA = new TH1D("hm_K0phiA", "K0#phi - KThr<pK",
		       160,M_D0-.4,M_D0+.4);
  hm_K0phiAC= new TH1D("hm_K0phiAC","K0#phi - KThr<pK,KinCuts",
		       160,M_D0-.4,M_D0+.4);
  //           ***** K0phi conditioned by D* *****
  hS_K0phi  = new TH1D("hS_K0phi", "K0#phi - D0pi=D*#pm3MeV",
		       240,M_D0-.6,M_D0+.6);
  hS_K0phiC = new TH1D("hS_K0phiC","K0#phi - D0pi=D*#pm3MeV,KinCuts",
		       240,M_D0-.6,M_D0+.6);
  // D*->K0phipi
  hs_K0phipi  = new TH1D("hs_K0phipi" ,"D0#pi - K0#phi=D0#pm30MeV",
		      120,M_D0+M_pi-.01,M_D0+M_pi+.05);
  hs_K0phipiC = new TH1D("hs_K0phipiC","D0#pi - K0#phi=D0#pm30MeV,KinCuts",
		      120,M_D0+M_pi-.01,M_D0+M_pi+.05);

  for (int iKp = 0; iKp<2; iKp++) {
    char hN[] = "hm_K0pCD2", hID[] = "e,#pi,KVeto,R>12  ";
    char hT[] = "K0p - subThr pID,R>12,KinCuts,MissM2<10  ";
    for (int ipID = 0; ipID<2; ipID++) {
      char Kp = iKp ? 'p' : 'K';
      // ********** K0 + K or p: 2 different pID **********
      //#define U3_K0X_RICH_RADIUS_CUT
#ifdef U3_K0X_RICH_RADIUS_CUT
#  ifdef RICH_ANGLE_CUT
      char cRCut[] = ",A>25"; sprintf(cRCut,",R>%.0f",1000*aRICHCut);
#  else    // ***** Radius cut *****
      char cRCut[] = ",R>12"; sprintf(cRCut,",R>%.0f",rRICHCut);
#  endif
#else
      char cRCut[] = "";
#endif
      if (ipID) sprintf(hID,"%cID%s",Kp,cRCut);
      else      sprintf(hID,"subThr %cID%s",Kp,cRCut);
      char cID = ipID ? '2' : '\0';
      sprintf(hN,"hm_K0%c%c",Kp,cID);   sprintf(hT,"K0%c - %s",Kp,hID);
      if (iKp) hm_K0p[ipID]   = new TH1D(hN,hT,325,1.4,2.9);
      else     hm_K0K[ipID]   = new TH1D(hN,hT,160,M_Ds-.4,M_Ds+.4);
      sprintf(hN,"hm_K0%cC%c",Kp,cID);  sprintf(hT,"K0%c - %s,KinCuts",Kp,hID);
      if (iKp) hm_K0pC[ipID]  = new TH1D(hN,hT,325,1.4,2.9);
      else     hm_K0KC[ipID]  = new TH1D(hN,hT,160,M_Ds-.4,M_Ds+.4);
      sprintf(hN,"hm_K0%cCM%c",Kp,cID); sprintf(hT,"K0%c - %s,KinCuts,MissM2<10",Kp,hID);
      if (iKp) hm_K0pCM[ipID] = new TH1D(hN,hT,325,1.4,2.9);
    }
  }
  hc_missM2 = new TH1D("hc_missM2", "MissMass2 - K0p=L_{c}#pm40MeV",
		       150,-50,100);
  //         ***** KINEMATICAL CUT *****
  hk_tK0phi  = new TH1D("hk_tK0phi","cos#theta* - K0#phi",200,-1,1);
  hk_tKSpi[0]= new TH1D("hk_tKSppi","cos#theta* - K*+#pi-",
			200,-1,1);
  hk_tKSpi[1]= new TH1D("hk_tKSmpi","cos#theta* - K*-#pi+",
			200,-1,1);

  //           ********** INCLUSIVE phi **********
  // - hm_iphi??: TH2D, in conjunction w/ a K0 candidate
  // - hm_Iphi??: TH1D, purely inclusive
  gDirectory->cd("/");
  TDirectory *dIphi; if (!(dIphi = (TDirectory*)gDirectory->Get("Iphi")))
    gDirectory->mkdir("Iphi","Inclusive phi Histos");
  if (!(dIphi = (TDirectory*)gDirectory->Get("Iphi"))) {
    printf("\n** U3:\a No creating subdir \"Iphi\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dIphi->cd();

  char iphiTitle[] = "K+K- - KThr<P:KID(1.005,1.01)|p<KThr:e(0.96),#pi(0.96)Veto,pT>100MeV ";
#ifdef REJECT_SubKThr_e
  sprintf(iphiTitle,"K+K- - P<KThr:e(%.2f),#pi(%.2f)Veto,pT>%d",
	  subKThrLHeVeto,subKThrLHpiVeto,Iphi_pTCUT);
#else
  sprintf(iphiTitle,"K+K- - P<KThr:#piVeto(%.2f),pT>%d",
	  subKThrLHpiVeto,Iphi_pTCUT);
#endif
  hm_iphibb = new TH2D("hm_iphibb",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi,3,.5,3.5);
  sprintf(iphiTitle,"K+K- - P<KThr:e#piVeto(%.2f),pT>%d",
	  subKThrLHpiVeto,Iphi_pTCUT);
  hm_Iphibb = new TH1D("hm_Iphibb",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi);  
#ifdef REJECT_SubKThr_e
  sprintf(iphiTitle,
	  "K+K- - KThr<P:KID(%.3f,%.2f)|P<KThr:e(%.2f),#pi(%.2f)Veto,pT>%d",
	  LHCut,LHBckCut,subKThrLHeVeto,subKThrLHpiVeto,Iphi_pTCUT);
#else
  sprintf(iphiTitle,
	  "K+K- -  KThr<P:KID(%.3f,%.2f)|P<KThr:#piVeto(%.2f),pT>%d",
	  LHCut,LHBckCut,subKThrLHpiVeto,Iphi_pTCUT);
#endif
  hm_iphiab = new TH2D("hm_iphiab",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi,3,.5,3.5);
  sprintf(iphiTitle,
	  "K+K- - KThr<P:KID(%.3f,%.2f)|P<KThr:e#piVeto(%.2f),pT>%d",
	  LHCut,LHBckCut,subKThrLHpiVeto,Iphi_pTCUT);
  hm_Iphiab = new TH1D("hm_Iphiab",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi);
  sprintf(iphiTitle,"K+K- -  KThr<P:KID(%.3f,%.2f),pT>%d",
	  LHCut,LHBckCut,Iphi_pTCUT);
  hm_iphiaa = new TH2D("hm_iphiaa",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi,3,.5,3.5);
  hm_Iphiaa = new TH1D("hm_Iphiaa",iphiTitle,
		       200,2*M_K-.02,M_phi+dM_phi);
}
// **********************************************************************
// ****************************** fillK0X *******************************
// **********************************************************************
void fillK0X(PaEvent &e, int ipV, int imuS,
	     const TLorentzVector &lvq,
	     int iET1, int iET2, const TLorentzVector &lvK0,
	     double KSMassCut, double dM_phi, double phiMassCut)
{
  double nu = lvq.E();
  const PaVertex &pV = e.vVertex(ipV); int nTrksPV = pV.NOutParticles();
  int ip3; double spT; for (ip3 = 0, spT = 0; ip3<nTrksPV; ip3++) {
    //    *************** DETERMINE SUM OF pT ***************
    int iEP3 = pV.iOutParticle(ip3);
    if (iEP3==imuS) continue;                               // ***** EXCLUDE mu'
    const PaParticle &pa3 = e.vParticle(iEP3);
    int iET3 = pa3.iTrack(); PaTrack &trk3 = e.vTrack()[iET3];
    int zones3 = tZones[iET3];
    //         ***** 3RD PARTICLE: REJECTIONS *****
    if (zones3==0x1) continue;                     // ***** EXCLUDE FRINGE FIELD
    if (trk3.Chi2tot()/trk3.Ndf()>chi2Cut) continue;    // ***** CUT ON CHI2/NDF
    const PaTPar &hv3 = pa3.ParInVtx(ipV); TVector3 v33 = hv3.Mom3();
    spT += fabs(v33.Perp(lvq.Vect()));
  }
  for (ip3 = 0; ip3<nTrksPV; ip3++) { 

    // *************************************************************
    // *************** SEARCH FOR K0 + pi,K,p,phi... ***************
    // *************************************************************

    int iEP3 = pV.iOutParticle(ip3);
    const PaParticle &pa3 = e.vParticle(iEP3);

    // ********** LOOP on PARTICLES in PRIMARY VERTEX **********

    if (iEP3==imuS) continue;                               // ***** EXCLUDE mu'
    int iET3 = pa3.iTrack();
    if (iET3==iET1 || iET3==iET2) continue;          // ***** EXCLUDE pi from K0

    PaTrack &trk3 = e.vTrack()[iET3]; int zones3 = tZones[iET3];

    //      ********** 3RD PARTICLE: REJECTIONS **********
    if (zones3==0x1) continue;                     // ***** EXCLUDE FRINGE FIELD
#ifdef DISCARD_MUONS
    if (trk3.XX0()>15) continue;                          // ***** EXCLUDE MUONS
#endif
    if (trk3.Chi2tot()/trk3.Ndf()>chi2Cut) continue;    // ***** CUT ON CHI2/NDF
    // Prior to, possibly re-scale momentum, determine acceptance by RICH and
    // exit angle for SM1 and incidence/exit angle for SM2.
    // (Cf. also comment supra, in V0 block, concenring Marcin's rescaling.)
    if (aRs[iET3]==0) /* Not done? */transport(e,trk3,aRs[iET3],zones3);
    if (aRs[iET3]<-1) continue; // Skip track starting downstream of RICH
    const PaTPar &hv3 = pa3.ParInVtx(ipV); TVector3 v33 = hv3.Mom3();
    // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
    const PaTPar &hi3 = trk3.vTPar()[0]; double mom3;
    mom3 = fabs(1/hi3(5));

    int Kok3 = 0, pok3 = 0;        // ***** RICH PID **********
    int abovepi = 0, aboveK = 0, abovep = 0;
    if (piThr*thrMargin<mom3) {
      abovepi |= 0x1;
      if (KThr*thrMargin<mom3) {
	aboveK |= 0x1;
	if (pThr*thrMargin<mom3) abovep |= 0x1;
      }
    }
    int piok3 = abovepi ? 1 : 2;
    if (PIDs[iET3]&0x10) {
      //                           ***** KID: above Thr = 1, below = 2
      double piLH = richLHs[0][iET3], KLH = richLHs[1][iET3];
      double pLH  = richLHs[2][iET3], eLH = richLHs[3][iET3];
#ifdef piS_e_REJECTION
      if (eLH>piLHeVeto*piLH)                               piok3 = 0;
#endif
      if (piLH<KLH || piLH<pLH || (abovepi && piLH<LHBckCut)) {
	/* */                                               piok3 = 0;
	if (aboveK) {
	  if (KLH>LHCut*piLH && KLH>LHCut*pLH && KLH>LHBckCut)
	    /* */                                            Kok3 = 1;
	}
	else if (abovepi) {
#ifdef REJECT_SubKThr_e
	  if (piLH<subKThrLHpiVeto && eLH<subKThrLHeVeto)    Kok3 = 2;
	  if (piLH<subpThrLHpiVeto && eLH<subpThrLHVeto)     pok3 = 2;
#else
	  if (piLH<subKThrLHpiVeto)                          Kok3 = 2;
	  if (piLH<subpThrLHpiVeto)                          pok3 = 2;
#endif
	}
	//else No telling K from pi

	// pID: above pThr = 1 = No piID nor KID subThreshold ID
	// Warning: infra not updated after definition of eID changed
	if (abovep) {
	  if (pLH>LHCut*piLH && pLH>LHCut*KLH && pLH>LHBckCut)
	    /* */                                            pok3 = 1;
	}
	else if (aboveK) {
#ifdef REJECT_SubKThr_e
	  if (piLH<subKThrLHpiVeto && KLH<subKThrLHpiVeto &&
	      eLH<subKThrLHeVeto)                            pok3 = 2;
#else
	  if (piLH<subKThrLHpiVeto && KLH<subKThrLHpiVeto)   pok3 = 2;
#endif
	}
	//else no telling p from K
      }
    }
#if defined K_BELOW_THR && K_BELOW_THR > 1
    bool checkAcceptance = piok3 || Kok3 || pok3 || !(PIDs[iET3]&0x10);
#else
    bool checkAcceptance = piok3 || Kok3 || pok3;
#endif
    if (checkAcceptance && aRs[iET3]==0) {
  printf("\n\n** U3:\a Run %d Evt %d Track %d RICH not yet checked!\n\n",
	 e.RunNum(),(int)e.UniqueEvNum(),iET3); assert(false);
    }
#if defined K_BELOW_THR && K_BELOW_THR > 1
    if (!(PIDs[iET3]&0x10) && aRs[iET3]>0 && abovepi) {
      /* */                                                 piok3 = 0;
      if      (!aboveK)                                      Kok3 = 2;
      else if (!abovep)                                      pok3 = 2;
    }
#endif
    if (piok3 || Kok3 || pok3) {    // ***** ID CONDITIONED by...*****
#ifdef U3_K0X_RICH_RADIUS_CUT
      if (aRs[iET3]<2) {                // ...RADIUS/ANGLE @ RICH 
	Kok3 = 0; pok3 = 0; piok3 = 0;
      }
#else
      if (aRs[iET3]<0) {                // ...RICH ACCEPTANCE *****
	Kok3 = 0; pok3 = 0; piok3 = 0;
      }
#endif
    }
    if ((piok3&&Kok3) || (pok3&&piok3)) {
  printf("\n** U3:\a Evt %d,%d Track %d(<-K0+X) RICH not yet checked\n\n",
	 e.RunNum(),(int)e.UniqueEvNum(),iET3); assert(false);
    }

    TLorentzVector *lvK0pi3 = 0, *lvpi3 = 0;
    if (piok3) {
      //  ***** piID: search for D+/- -> K0pi+/- or K*(->K0pi)pi *****
      double Epi3 = sqrt(v33.Mag2()+M2_pi);
      lvpi3 = new TLorentzVector(v33,Epi3);
      lvK0pi3 = new TLorentzVector(*lvpi3); *lvK0pi3 += lvK0;
      double m_K0pi = lvK0pi3->M(); hm_K0pi->Fill(m_K0pi);
      TLorentzVector lvK02 = lvK0;	       // Boost K0 to K0pi frame
      lvK02.Boost(-lvK0pi3->BoostVector());
      double cthstar = cos(lvK02.Angle(lvK0pi3->Vect()));
      if (nu!=0 && lvK0pi3->Energy()/nu>.25 && fabs(cthstar)<.75) {
	hm_K0piC->Fill(m_K0pi);
      }
      if (fabs(m_K0pi-M_KS)>KSMassCut) {
	delete lvK0pi3; lvK0pi3 = 0;
      }
    }
    else if (pok3 && pa3.Q()>0) {
      //        ********** p+ID: search for Lc and 5q **********
      double Ep3 = sqrt(v33.Mag2()+M2_p);
      TLorentzVector lvp3(v33,Ep3);
      TLorentzVector lvK0p = lvK0 + lvp3; double m_K0p = lvK0p.M();
      hm_K0p[pok3-1]->Fill(m_K0p);
      TLorentzVector lvK02 = lvK0;
      lvK02.Boost(-lvK0p.BoostVector());	// Boost K0 to K0p frame
      double cthstar = cos(lvK02.Angle(lvK0p.Vect()));
      if (nu!=0 && lvK0p.Energy()/nu>.25 && fabs(cthstar)<.75) {
	hm_K0pC[pok3-1]->Fill(m_K0p);
      }
      TLorentzVector lvp(0,0,0,M_p), lvMiss = lvp+lvq-lvK0p;
      double missM2 = lvMiss.M2();
      if (fabs(m_K0p-M_Lc)<.04) hc_missM2->Fill(missM2);
      if (missM2<10) hm_K0pCM[pok3-1]->Fill(m_K0p);
    }
    //      ********** Search for phi -> KK **********
    // Here we have in view to look for phi w/ possibly only one KID
    double EK3 = sqrt(v33.Mag2()+M2_K);
    TLorentzVector lvK3(v33,EK3);
    if (Kok3) {
      //     ********** KID: search for Ds+/- -> K0K+/- **********
      TLorentzVector lvK0K(lvK3); lvK0K += lvK0;
      double m_K0K = lvK0K.M(); if (fabs(m_K0K-M_Ds)<.4) {
	hm_K0K[Kok3-1]->Fill(m_K0K);
	TLorentzVector lvK02 = lvK0;
	lvK02.Boost(-lvK0K.BoostVector());  // Boost K0 to K0K frame
	double cthstar = cos(lvK02.Angle(lvK0K.Vect()));
	if (nu!=0 && lvK0K.Energy()/nu>.25 && fabs(cthstar)<.75) {
	  hm_K0KC[Kok3-1]->Fill(m_K0K);
	}
      }
    }

    if (pa3.Q()<0 ||        // ***** KEEP SEARCHING ONLY if 3 > 0 OR...
	(!Kok3 && !piok3)) {// ...CANDIDATES D0 -> K0 phi OR -> K0* pi
      if (lvK0pi3) delete lvK0pi3; if (lvpi3) delete lvpi3;
      continue;
    }

    vector<UsParticle> KSpis, K0phis;
    for (int ip4 = 0; ip4<nTrksPV; ip4++) {
      int iEP4 = pV.iOutParticle(ip4);
      const PaParticle &pa4 = e.vParticle(iEP4);
      if (pa4.Q()>0) continue;             // Particle 4 is the <0

      // ***** LOOP on <0 PARTICLES in PRIMARY VERTEX *****

      if (iEP4==imuS) continue;                   // ***** EXCLUDE mu'
      int iET4 = pa4.iTrack();
      if (iET4==iET2) continue;           // ***** EXCLUDE pi- from K0

      PaTrack &trk4 = e.vTrack()[iET4];
      int zones4 = tZones[iET4];

      // ***** 4th Particle: REJECTIONS *****
      if (zones4==0x1) continue;                   // ***** EXCLUDE FRINGE FIELD
#ifdef DISCARD_MUONS
      if (trk4.XX0()>15) continue;                        // ***** EXCLUDE MUONS
#endif
      if (trk4.Chi2tot()/trk4.Ndf()>chi2Cut) continue;  // ***** CUT ON CHI2/NDF

      // Prior to, possibly re-scale momentum, determine acceptance by RICH and
      // exit angle for SM1 and incidence/exit angle for SM2.
      // (Cf. also comment supra, in V0 block, concenring Marcin's rescaling.)
      if (aRs[iET4]==0) /* !done? */transport(e,trk4,aRs[iET4],zones4);
      if (aRs[iET4]<-1) continue; // Skip track starting downstream of RICH
      const PaTPar &hv4 = pa4.ParInVtx(ipV); TVector3 v43 = hv4.Mom3();
      // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
      const PaTPar &hi4 = trk4.vTPar()[0]; double mom4;
      mom4 = fabs(1/hi4(5));

      if (piThr*thrMargin<mom4) {           // ***** RICH KID *****
	abovepi |= 0x2;
	if (KThr*thrMargin<mom4) aboveK |= 0x2;
	else                     aboveK &= 0x1;
      }
      else abovepi &= 0x1;
      int piok4 = (abovepi&0x2) ? 1 : 2, Kok4 = 0;
      if (PIDs[iET4]&0x10) {
	double piLH = richLHs[0][iET4], KLH = richLHs[1][iET4];
	double pLH  = richLHs[2][iET4], eLH = richLHs[3][iET4];
#ifdef piS_e_REJECTION
	if (eLH>piLHeVeto*piLH)                             piok4 = 0;
#endif
	if (piLH<KLH || piLH<pLH || ((abovepi&0x2) && piLH<LHBckCut)) {
	  /* */                                             piok4 = 0;
	  if (aboveK&0x2) {          // KID: above Thr = 1, below = 2
	    if (KLH>LHCut*piLH && KLH>LHCut*pLH && KLH>LHBckCut)
	      /* */                                          Kok4 = 1;
	  }
	  else if (abovepi&0x2) {
#ifdef REJECT_SubKThr_e
	    if (piLH<subKThrLHpiVeto && eLH<subKThrLHeVeto)  Kok4 = 2;
#else
	    if (piLH<subKThrLHpiVeto)                        Kok4 = 2;
#endif
	  }
	}
      }
#if defined K_BELOW_THR && K_BELOW_THR > 1
      checkAcceptance =  piok4 || Kok4 || !(PIDs[iET4]&0x10);
#else
      checkAcceptance =  piok4 || Kok4;
#endif
      if (checkAcceptance && aRs[iET4]==0) {
  printf("\n** U3:\a Evt %d,%d Track %d(<-K0+XX') RICH not yet checked\n\n",
	 e.RunNum(),(int)e.UniqueEvNum(),iET4); assert(false);
      }
#if defined K_BELOW_THR && K_BELOW_THR > 1
      if (!(PIDs[iET4]&0x10) && aRs[iET4]>0 && (abovepi&0x2)) {
	/* */                                               piok4 = 0;
	if      (!(aboveK&0x2))                              Kok4 = 2;
      }
#endif
      if (piok4 || Kok4) {        // ***** ID CONDITIONED by... *****
#ifdef U3_K0X_RICH_RADIUS_CUT
	if (aRs[iET4]<2) {        // ...RADIUS/ANGLE @ RICH
	  piok4 = 0 ; Kok4 = 0;
	}
#else
	if (aRs[iET4]<0) {        // ...RICH ACCEPTANCE *****
	  piok4 = 0; Kok4 = 0;
	}
#endif
      }
      if (piok4&&Kok4) {
   printf("\n** U3:\a Run %d Evt %d 4-Track %d: PID inconsistency\n\n",
	 e.RunNum(),(int)e.UniqueEvNum(),iET4); assert(false);
      }

      if (piok3 && piok4) {
	// ********** piID: search for D0 -> K*(->K0pi)pi **********
	double Epi4 = sqrt(v43.Mag2()+M2_pi);
	TLorentzVector lvpi4(v43,Epi4);
	TLorentzVector lvK0pi4(lvpi4); lvK0pi4 += lvK0;
	double m_K0pi = lvK0pi4.M();
	int pm = -1; if (fabs(m_K0pi-M_KS)<KSMassCut) {
	  pm = 1; lvK0pi4 += *lvpi3;
	  KSpis.push_back(UsParticle(lvK0pi4,ipV,iET3,iET4));
	  KSpis.back().type = 2; // UsParticle.type==1: K*-
	}
	else if (lvK0pi3) {
	  pm = 0; TLorentzVector lvKSpi(*lvK0pi3); lvKSpi += lvpi4;
	  KSpis.push_back(UsParticle(lvKSpi,ipV,iET3,iET4));
	  KSpis.back().type = 1; // UsParticle.type==1: K*+
	}
	if (pm>=0) {
	  TLorentzVector &lvKSpi = KSpis.back().lv;
	  double m_KSpi = lvKSpi.M(); hm_KSpi[pm]->Fill(m_KSpi);
	  if (spT>1.6) {
	    hm_KSpit[pm]->Fill(m_KSpi);
	    if (spT>2) hm_KSpiT[pm]->Fill(m_KSpi);
	  }
	}
      }

      int okphi = 0; if (Kok3) okphi |= 0x1; if (Kok4) okphi |= 0x2;
      if (okphi) {
	// ********** KID: search for D0-> K0phi(->KK) **********
	double EK4 = sqrt(v43.Mag2()+M2_K);
	TLorentzVector lvK4(v43,EK4);
	TLorentzVector lvKK(lvK3); lvKK += lvK4;
	double pT34 = lvK4.Perp(lvK3.Vect());
	if (pT34<Iphi_pTCUT*.001) continue;
	if      (aboveK==0)   hm_iphibb->Fill(lvKK.M(),okphi);
	else if (aboveK==0x3) hm_iphiaa->Fill(lvKK.M(),okphi);
	else                  hm_iphiab->Fill(lvKK.M(),okphi);

	if (fabs(lvKK.M()-M_phi)>phiMassCut || // phi mass cut...
	    okphi!=0x3) continue;              // K+ AND K- ID

	//                        ***** phi SELECTION *****
	TLorentzVector lvK0phi = lvK0 + lvKK;
	double m_K0phi = lvK0phi.M();
	if (fabs(m_K0phi-M_D0)>.6) continue;
	K0phis.push_back(UsParticle(lvK0phi,ipV,iET3,iET4));
	K0phis.back().type = aboveK ? 0xc : 0x4; 
	hm_K0phi->Fill(m_K0phi); if (aboveK) hm_K0phiA->Fill(m_K0phi);
      }
    }

    int nKSpis = (int)KSpis.size();
    int nKXs = nKSpis+(int)K0phis.size();
    for (int iKX = 0; iKX<nKXs; iKX++) {
      UsParticle &KX = iKX<nKSpis ? KSpis[iKX] : K0phis[iKX-nKSpis];
      TLorentzVector &lvKX = KX.lv; double m_KX = lvKX.M();
      if (fabs(m_KX-M_D0)>.6) continue;
      int type = KX.type, iET4 = KX.i2;
      TLorentzVector lvK02 = lvK0;
      lvK02.Boost(-lvKX.BoostVector());  // Boost K0 to K*pi frame
      double cthstar = cos(lvK02.Angle(lvKX.Vect()));
      if (type==1 || type==2) hk_tKSpi[type-1]->Fill(cthstar);
      else if (type&0x4)             hk_tK0phi->Fill(cthstar);
      else {
	printf("** U3:\a Run %d Evt %d: KX type(=%d) unknown\n\n",
	 e.RunNum(),(int)e.UniqueEvNum(),type); assert(false);
      }
      int kinOK = 0; if (nu!=0) {
	double zKX = lvKX.Energy()/nu;
	if (zKX>zCutS && cLowCutS<cthstar && cthstar<cUpCutS)
	  kinOK |= 0x1;
	if (zKX>zCut0 && cLowCut0<cthstar && cthstar<cUpCut0)
	  kinOK |= 0x2;
	if (kinOK&0x2) {
	  if (type==1 || type==2) hm_KSpiC[type-1]->Fill(m_KX);
	  else if (type&0x4) {
	    hm_K0phiC->Fill(m_KX);
	    if (type&0x8) hm_K0phiAC->Fill(m_KX);
	  }
	}
      }
      double ED0 = sqrt(lvKX.Vect().Mag2()+M_D0*M_D0);
      TLorentzVector lvD0(lvKX.Vect(),ED0);

      for (int ip5 = 0; ip5<nTrksPV; ip5++) {
	//   ********** SEARCH FOR D*->D0(->K0phi | K*pi)pi **********
	int iEP5 = pV.iOutParticle(ip5);
	const PaParticle &pa5 = e.vParticle(iEP5);
	if (iEP5==imuS) continue;                 // ***** EXCLUDE mu'
	int iET5 = pa5.iTrack();
	if (iET5==iET1 || iET5==iET2 ||       // Exclude pi+/- from K0
	    iET5==iET3 || iET5==iET4)         // ...and from phi/K*pi
	  continue;

	PaTrack &trk5 = e.vTrack()[iET5];//double zL5 = trk5.ZLast();
	int zones5 = tZones[iET5];

	// *************** REJECTIONS ***************
	if (zones5==0x1) continue;             // Exclude fringe field
#ifdef DISCARD_MUONS
	if (trk5.XX0()>15) continue;                  // Exclude muons
#endif
	if (trk5.Chi2tot()/trk5.Ndf()>chi2Cut)
	  continue;                           // Cut on chi2/ndf
#ifdef piS_e_REJECTION
	if (PIDs[iET5]&0x10) {
	  double eLH = richLHs[3][iET5], piLH = richLHs[2][iET5];
	  if (eLH>piLHeVeto*piLH) continue;
	}
#endif

	// Prior to, possibly re-scale momentum, determine acceptance by RICH
	// and exit angle for SM1 and incidence/exit angle for SM2.
	// (Cf. also comment supra, in V0 block, concenring Marcin's rescaling.)
	if (aRs[iET5]==0)/* !done? */transport(e,trk5,aRs[iET5],zones5);
	if (aRs[iET5]<-1) continue; // Skip track starting downstream of RICH
	const PaTPar &hv5 = pa5.ParInVtx(ipV); TVector3 v53 = hv5.Mom3();
	double Epi5 = sqrt(v53.Mag2()+M2_pi);
	TLorentzVector lv53(v53,Epi5);

	TLorentzVector lvD0pi = lvD0 + lv53;
	double m_D0pi = lvD0pi.M();
	if (fabs(m_KX-M_D0)<.03) {
	  if (type==1 || type==2) {
	    int pm = type-1, pmwr = pm + (pa5.Q()!=2*pm-1 ? 2 : 0); 
	    hm_KSpipi[pmwr]->Fill(m_D0pi);
	    if (kinOK&0x1) hm_KSpipiC[pmwr]->Fill(m_D0pi);
	  }
	  else if (type&0x4) {
	    hs_K0phipi->Fill(m_D0pi);
	    if (kinOK&0x1) hs_K0phipiC->Fill(m_D0pi);
	  }
	}
	if (D0piLow<m_D0pi-M_DS && m_D0pi-M_DS<D0piUp) {
	  if (type==1 || type==2) {
	    int pm = type-1, pmwr = pm + (pa5.Q()!=2*pm-1 ? 2 : 0);
	    hS_KSpi[pmwr]->Fill(m_KX);
	    if (kinOK&0x1) hS_KSpiC[pmwr]->Fill(m_KX);
	    if (spT>1.6) {
	      hS_KSpit[pmwr]->Fill(m_KX);
	      if (spT>2) hS_KSpiT[pmwr]->Fill(m_KX);
	    }
	  }
	  else if (type&0x4) {
	    hS_K0phi->Fill(m_KX);
	    if (kinOK) hS_K0phiC->Fill(m_KX);
	  }
	}
      }  // End loop on pa5
    }  // End loop on K*pi/K0phi's
    if (lvK0pi3) delete lvK0pi3; if (lvpi3) delete lvpi3;
  }  // End loop on pa3
}
// **********************************************************************
// *************************    fillLambdaX    **************************
// **********************************************************************
void fillLambdaX(PaEvent &e, int ipV, int imuS,
		 const TLorentzVector &lvq,
		 int iET1, int iET2, int iLaL, const TLorentzVector *lvL)
{
  const PaVertex &pV = e.vVertex(ipV); int nTrksPV = pV.NOutParticles();
  for (int ip3 = 0; ip3<nTrksPV; ip3++) { 
    // *************************************************************
    //    ********** SEARCH FOR Sigma* -> Lambda pi **********
    // *************************************************************
    int iEP3 = pV.iOutParticle(ip3);
    const PaParticle &pa3 = e.vParticle(iEP3);

    // ********** LOOP on PARTICLES in PRIMARY VERTEX **********

    if (iEP3==imuS) continue;                               // ***** EXCLUDE mu'
    int iET3 = pa3.iTrack();
    if (iET3==iET1 || iET3==iET2) continue;    // ***** EXCLUDE p,pi from Lambda

    PaTrack &trk3 = e.vTrack()[iET3]; int zones3 = tZones[iET3];

    //      ********** 3RD PARTICLE: REJECTIONS **********
    if (zones3==0x1) continue;                     // ***** EXCLUDE FRINGE FIELD
#  ifdef DISCARD_MUONS
    if (trk3.XX0()>15) continue;                          // ***** EXCLUDE MUONS
#  endif
    if (trk3.Chi2tot()/trk3.Ndf()>chi2Cut) continue;    // ***** CUT ON CHI2/NDF

    // Prior to, possibly re-scale momentum, determine acceptance by RICH and
    // exit angle for SM1 and incidence/exit angle for SM2.
    // (Cf. also comment supra, in V0 block, concenring Marcin's rescaling.)
    if (aRs[iET3]==0) /* Not done? */transport(e,trk3,aRs[iET3],zones3);
    if (aRs[iET3]<-1) continue; // Skip track starting downstream of RICH
    const PaTPar &hv3 = pa3.ParInVtx(ipV); TVector3 v33 = hv3.Mom3();
    // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
    const PaTPar &hi3 = trk3.vTPar()[0]; double mom3;
    mom3 = fabs(1/hi3(5));
#  ifdef piS_e_REJECTION
    if (PIDs[iET3]&0x10) {
      double piLH = richLHs[0][iET3], eLH = richLHs[3][iET3];
      if (eLH>piLHeVeto*piLH) continue;
    }
#  endif

    double Epi3 = sqrt(v33.Mag2()+M2_pi);
    TLorentzVector lvpi3(v33,Epi3);
    TLorentzVector lvLpi(*lvL); lvLpi += lvpi3;
    int pm = pa3.Q()>0 ? 0 : 1;
    double m_Lpi = lvLpi.M(); hm_Lpi[iLaL+2*pm]->Fill(m_Lpi);

    if (!(PIDs[iET3]&0x10)) continue;                     // ***** Sigma w/ piID
    double piLH = richLHs[0][iET3], KLH = richLHs[1][iET3];
    double pLH  = richLHs[2][iET3];
    if (rich->piThr<mom3 && mom3<PRICHCut &&
	piLH>LHCut*KLH && piLH>LHCut*pLH && piLH>LHBckCut) {
      hm_LpiID[iLaL+2*pm]->Fill(m_Lpi);
    }
  }  // End search fo Sigma* -> L + pi
}
// **********************************************************************
// *************************   fillLCascades   **************************
// **********************************************************************
void fillLCascades(PaEvent &e, int ipV, TMatrixD &CovP, int imuS,
		   int isV, TMatrixD &CovS, int iET1, int iET2, int iLaL,
		   const TLorentzVector *lvL, unsigned short ppiID,
		   bool &uDSTSelection)
{
  const PaVertex &pV = e.vVertex(ipV);
  double Xp = pV.Pos(0), Yp = pV.Pos(1), Zp = pV.Pos(2);
  const PaVertex &sV = e.vVertex(isV);
  double Xs = sV.Pos(0), Ys = sV.Pos(1), Zs = sV.Pos(2);
  const TVector3 &vL = lvL->Vect();
  // ********** Lambda PARAMETRISATION : "hL"
  PaTPar hL(Zs,Xs,Ys,vL[0]/vL[2],vL[1]/vL[2],1/vL.Mag());
  int iEP1 = sV.iOutParticle(0), iEP2 = sV.iOutParticle(1);
  // Get the Pa3Par's of the two particles in sVertex
  const PaParticle &pa1 = e.vParticle(iEP1), &pa2 = e.vParticle(iEP2);
  int jV1, jV2, jvtx; for (jvtx = 0, jV1 = -1; jvtx<pa1.NVertex(); jvtx++) {
    if (pa1.iVertex(jvtx)==isV) { jV1 = jvtx; break; }
  }
  /* */               for (jvtx = 0, jV2 = -1; jvtx<pa2.NVertex(); jvtx++) {
    if (pa2.iVertex(jvtx)==isV) { jV2 = jvtx; break; }
  }
  if (jV1<0 || jV2<0) {
    unsigned int paPat = 0; if (jV1<0) paPat |= 0x1; if (jV2<0) paPat |= 0x2; 
    printf("** U3: Inconsistency: Evt %d#%d Vtx %d Particles 0x%x w/ %d in vertex list\n",
	   Run,EvNum,isV,paPat,isV);
    abort();
  }
  const Pa3Par &p31 = pa1.vVtxPars()[jV1], &p32 = pa2.vVtxPars()[jV2];
  double P1 = 1/fabs(p31(3)), P2 = 1/fabs(p32(3));
  double w1 = sqrt(1.+ p31(1)*p31(1)+p31(2)*p31(2)), c1 = P1/w1;
  double w2 = sqrt(1.+ p32(1)*p32(1)+p32(2)*p32(2)), c2 = P2/w2;
  //#define U3_DEBUG_CDA
#ifdef U3_DEBUG_CDA
  static int jdebug = 0;
  if (jdebug) {
    Pa3Par p3L;
    p3L(1) = (c1*p31(1)+c2*p32(1))/(c1+c2);
    p3L(2) = (c1*p31(2)+c2*p32(2))/(c1+c2);
    p3L(3) = hL(5);
    TVector3 vLp = p3L.Mom3();
    printf("Evt #%d:",EvNum);
    for (int i = 0; i<3; i++) printf(" %.3f %.3f ",vL[i],vLp[i]);
    printf("\n");
  }
#endif
  // ***** COVARIANCE
  // Simplification: the two coefficients "c(1|2)" are fixed.
  hL(1,1) = sV.Cov(0); hL(2,1) = sV.Cov(1); hL(2,2) = sV.Cov(2);
  hL(3,3) = (c1*p31(1,1)+c2*p32(1,1))/(c1+c2);
  hL(4,3) = (c1*p31(2,1)+c2*p32(2,1))/(c1+c2);
  hL(4,4) = (c1*p31(2,2)+c2*p32(2,2))/(c1+c2);

  // ***** V0 SELECTION: RECALCULATE
  float dZp = CovP(2,2), dZs = CovS(2,2), Z0 = (Zp+Zs)/2;
  TVector3 vv3(Xs-Xp,Ys-Yp,Zs-Zp);	// Vertex Line
  double vvv0 = vv3*lvL->Vect(), dist = vv3.Mag();
  double ctheta = vvv0/dist/lvL->Vect().Mag();
  // ***** RICH
  int PID = 0; // 0x1: piID, 0x2: pID,
  if (ppiID&0x30) PID += 1;
  if (ppiID&0x6)  PID += 2;

  int nTrks = e.vTrack().size(); for (int iET3 = 0; iET3<nTrks; iET3++) {

    // ********** LOOP on PARTICLES in PRIMARY VERTEX **********

    if (iET3==iET1 || iET3==iET2) continue;    // ***** EXCLUDE p,pi FROM Lambda
		
    PaTrack &trk3 = e.vTrack()[iET3]; int zones3 = tZones[iET3];
    PaTPar hi3 = trk3.vTPar()[0];    
    if (!hi3.HasMom()) continue;                       // ***** REQUIRE MOMENTUM
    int iEP3 = trk3.iParticle(); if (iEP3<0) continue;
    if (iEP3==imuS) continue;                               // ***** EXCLUDE mu'
    const PaParticle &pa3 = e.vParticle(iEP3);
    if (pa3.IsBeam()) continue;                            // ***** EXCLUDE BEAM
    if ((!iLaL && pa3.Q()>0) || (iLaL && pa3.Q()<0))         // ***** Lpi-/aLpi+
      continue;
    //                      ***** 3RD PARTICLE: REJECTIONS *****
    if (zones3==0x1) continue;                     // ***** EXCLUDE FRINGE FIELD
#  ifdef DISCARD_MUONS
    if (trk3.XX0()>15) continue;                          // ***** EXCLUDE MUONS
#  endif
    if (trk3.Chi2tot()/trk3.Ndf()>chi2Cut) continue;    // ***** CUT ON CHI2/NDF

    //                ********** REQUIRE CASCADE **********
    // - Cut on distance @ CDA
    // - Set "s_d3h" flag.
    //#define U03_DEBUG_CDA
#  ifdef U3_DEBUG_CDA
    static int debug = 0;
    if (debug) {
      PaTPar H(hL), Hi3 = hi3;
      printf("DEBUG_CDA: %.2f,%.2f,%.2f / %.2f,%.2f,%.2f / %.2f[%.2f,%.2f]",
	     Hi3(1),Hi3(2),Hi3(0),H(1),H(2),H(0),Z0,Zp+3*dZp,Zs-3*dZs);
      if (H.FindCDA(Hi3,Z0,Zp+3*dZp,Zs-3*dZs,true)) {
	printf(" -> %.2f,%.2f,%.2f / %.2f,%.2f,%.2f\n",
	       Hi3(1),Hi3(2),Hi3(0),H(1),H(2),H(0));
      }
      else printf("\n");
    }
#  endif
    PaTPar h(hL) /* Keep "hL" unchanged */;  double dZCDA;
    // Double-precision covariance is used in the extrapolations performed in
    // "FindCDA". Somehow we need explict call to "f2d" to init it correctly.
    h.f2d();
    if (!h.FindCDA(hi3,Z0,Zp+3*dZp,Zs-3*dZs,true,false,dZCDA)) continue;
    // Cascade origin: we have two different estimates => Average? although
    // could be that the one from Lambda, w/ its two tracks and vertex fit is
    // more accurate...
    //static double a = .5;
    //double b = 1-a;
    //double Xc = (a*hi3(1)+b*h(1))/2, Yc = (a*hi3(2)+b*h(2))/2, Zc = h(0);
    // => Turned out that a = 0 yield highest "cthpcCa" for evt 105951061 of
    //   2016P09#275479 (volume .002 of slot5.1).
    double Xc = h(1),   Yc = h(2), Zc = h(0);
    double X3 = hi3(1), Y3 = hi3(2), dX3h = X3-Xc, dY3h = Y3-Yc;
    double d3h = sqrt(dX3h*dX3h+dY3h*dY3h);
    double dd3h = dX3h*dX3h*(hi3(1,1)+h(1,1))+dY3h*dY3h*(hi3(2,2)+h(2,2))+
      2*dX3h*dY3h*(hi3(2,1)+h(2,1));
    int s_d3h; if (dd3h<0) {
      printf("** U3: Inconsistency: Evt %d#%d Vtx,Particle %d,%d: dd3h error %f\n",
	     Run,EvNum,isV,iET3,dd3h);
      s_d3h = -1;
    }
    else {
      dd3h = sqrt(dd3h)/d3h;
      if      (d3h/dd3h>16) s_d3h = -1;
      else if (d3h/dd3h>8)  s_d3h = 0;
      else                  s_d3h = 1;
    }
    if (s_d3h<0) // To speed up processing:          *****  Loose cut on D @ CDA
      continue;
#  ifdef U3_DEBUG_CDA
    if (debug)
      printf(" -> %.2f,%.2f,%.2f / %.2f,%.2f,%.2f => %.2f+/-%.2f\n",
	     hi3(1),hi3(2),hi3(0),h(1),h(2),h(0),d3h,dd3h);
#  endif
    // ***** REQUIRE EXTRA MARGIN...
    // ...w.r.t. FindCDA, which takes into account the uncertainty on ZCDA (a
    // precaution that can be useful in the context of the PV search, but can
    // lead here to CDA upstream of pV, i.e. something which will later on
    // never pass the cascade selection cuts, wichsoever.
    if (Zc<Zp+3*dZp) continue;

    //       ********** RICH and PID **********
    // Determine acceptance by RICH and exit angle for SM1 and incidence/exit
    // angle for SM2.
    // (Cf. also comment supra, in V0 block, concerning Marcin's rescaling.)
    if (aRs[iET3]==0) /* Not done? */transport(e,trk3,aRs[iET3],zones3);
    if (aRs[iET3]<-1) continue; // Skip track starting downstream of RICH
    TVector3 v33 = hi3.Mom3(); // 3rd particle 3-vector @ point of CDA
    // ***** PARTICLE ID
    int Ppi3ID = PID; if (Ppi3ID>=1) {
      if (PIDs[iET3]&0x10) {         // ***** SEARCH for Xi *****
	double mom3 = hi3.Mom();     // ***** SEARCH for Omega *****
	double piLH = richLHs[0][iET3], eLH = richLHs[3][iET3];
	if (thrMargin*piThr<mom3) { if (eLH<piLHeVeto*piLH) Ppi3ID += 3; }
	else if (eLH<subpiThrLHeVeto) Ppi3ID += 3;
      }
    }
    int PK3ID = PID; if (PK3ID>=1) {
      double mom3 = hi3.Mom();      // ***** SEARCH for Omega *****
      if (rich->piThr*subKThrPmin<mom3 && mom3<PRICHCut) {
	int KID = 0; if (PIDs[iET3]&0x10) {
	  double piLH = richLHs[0][iET3], KLH = richLHs[1][iET3];
	  double pLH =  richLHs[2][iET3], eLH = richLHs[3][iET3];
	  if (rich->KThr*thrMargin<mom3) {
	    if (KLH>LHCut*piLH && KLH>LHCut*pLH && KLH>LHBckCut)
	      KID = 0x1;
	  }
	  else {
	    if (piLH<subKThrLHpiVeto) KID = 0x2;
#  ifdef REJECT_SubKThr_e
	    if (eLH>subKThrLHeVeto) KID = 0;
#  endif
	  }
	}
	else if (mom3<rich->KThr) {
	  if (aRs[iET3]==0) {
	    printf("\n** U3:\a Evt %d,%d Track %d(<-Omega) RICH not yet checked\n\n",
		   e.RunNum(),(int)e.UniqueEvNum(),iET2); assert(false);
	  }
	  if (aRs[iET3]) KID = 0x2;
	}
	if (KID) PK3ID += 3;
      }
    }

    //              ********** INVARIANT MASS **********
    double Epi3 = sqrt(v33.Mag2()+M2_pi);
    TLorentzVector lvpi3(v33,Epi3);
    TLorentzVector lvLpi(*lvL); lvLpi += lvpi3;
    double m_Lpi = lvLpi.M();                    // ***** Lambda+pi
    double EK3 = sqrt(v33.Mag2()+M2_K);
    TLorentzVector lvK3(v33,EK3);
    TLorentzVector lvLK(*lvL); lvLK += lvK3;
    double m_LK = lvLK.M();                      // ***** Lambda+K

    //             ********** CASCADE SELECTION CUTS **********
    // ***** NO CUT
    hm_LpiCa[iLaL]->Fill(m_Lpi,Ppi3ID);
    hm_LKCa[iLaL]->Fill(m_LK,PK3ID);
    // ***** PRIMARY<->CASCADE
    TVector3 vpc(Xc-Xp,Yc-Yp,Zc-Zp);  // pV-Cascade Vertex Line
    double pcCa = vpc*lvLpi.Vect(), dipc = vpc.Mag();
    double cthpcCa = pcCa/dipc/lvLpi.Vect().Mag();
    dipc = cthpcCa>0 ? dipc : -dipc;
    TVectorD Vpc(0,2,Xc-Xp,Yc-Yp,Zc-Zp,"END");
    TVectorD tmp = Vpc; tmp *= CovP;
    double ddipc = tmp*Vpc+Vpc[2]*Vpc[2]*dZCDA*dZCDA;
    ddipc = sqrt(ddipc)/fabs(dipc);
    // ***** CASCADE<->SECONDARY (Lambda decay)
    TVector3 vcs(Xs-Xc,Ys-Yc,Zs-Zc);  // Cascade-Secondary Vertex Line
    double csV0 = vcs*vL, dics = vcs.Mag();
    dics = csV0>0 ? dics : -dics;
    TVectorD Vcs(0,2,Xs-Xc,Ys-Yc,Zs-Zc,"END");
    tmp = Vcs; tmp *= CovS;
    double ddics = tmp*Vcs+Vcs[2]*Vcs[2]*dZCDA*dZCDA;
    ddics = sqrt(ddics)/fabs(dics);
    unsigned int s_dipc = dipc>3 ? 0x1 : 0;
    if (dipc>8*ddipc) s_dipc |= 0x2;
    int s_pcCa; if (cthpcCa>.9998) { s_pcCa = cthpcCa>.99995 ? 2 : 1; }
    else                             s_pcCa = 0;
    unsigned int s_dics = dics>3 ? 0x1 : 0;
    if (dics>12*ddics) s_dics |= 0x2;
    // ***** INV.MASS HISTOS
    unsigned int cascade = 0;
    if (s_d3h && (s_dipc&0x1) && s_pcCa) {
      cascade |= 0x1;
      hm_LpiCa3[iLaL]->Fill(m_Lpi,Ppi3ID);
      if (s_pcCa>1) hm_LKCa3 [iLaL]->Fill(m_LK, PK3ID);
    }
    if (s_d3h && (s_dipc&0x2) && s_pcCa) {
      cascade |= 0x2;
      hm_LpiCaD[iLaL]->Fill(m_Lpi,Ppi3ID);
      if (s_pcCa>1) hm_LKCaD [iLaL]->Fill(m_LK, PK3ID);
    }
    if ((cascade&0x1) && (s_dics&0x1)) {
      hm_LpiCa6[iLaL]->Fill(m_Lpi,Ppi3ID);
      if (s_pcCa>1) hm_LKCa6 [iLaL]->Fill(m_LK, PK3ID);
    }
    if ((cascade&0x2) && (s_dics&0x2)) {
      hm_LpiCaC[iLaL]->Fill(m_Lpi,Ppi3ID);
      if (s_pcCa>1) hm_LKCaC [iLaL]->Fill(m_LK, PK3ID);
    }
    // ********** CASCADE SELECTION HISTOS
    //                                             ***** Xi SIGNAL and SIDEBANDS
    if      (fabs(m_Lpi-M_Xi)<.008) {
      // The cut on M_Lpi is somwhat arbitrary, not based on a Gaussian fit
      if ((s_dipc&0x1) && s_pcCa && (s_dics&0x1)) hs_d3hXi->Fill(dd3h);
      if ((s_dipc&0x2) && s_pcCa && (s_dics&0x2)) hs_D3hXi->Fill(d3h/dd3h);
      if (s_pcCa && (s_dics&0x1) && s_d3h)        hs_dpcXi->Fill(dipc);
      if (s_pcCa && (s_dics&0x2) && s_d3h)        hs_DpcXi->Fill(dipc/ddipc);
      if ((s_dipc&0x1) && (s_dics&0x1) && s_d3h)  hs_cpcXi->Fill(cthpcCa);
      if ((s_dipc&0x2) && (s_dics&0x2) && s_d3h)  hs_CpcXi->Fill(cthpcCa);
      if ((s_dipc&0x1) && s_pcCa && s_d3h)        hs_dcsXi->Fill(dics);
      if ((s_dipc&0x2) && s_pcCa && s_d3h)        hs_DcsXi->Fill(dics/ddics);
      if ((s_dipc&0x1) && s_pcCa && (s_dics&0x1) && s_d3h)
	hs_cpsLa->Fill(ctheta);
      if ((s_dipc&0x2) && s_pcCa && (s_dics&0x2) && s_d3h)
	hs_CpsLa->Fill(ctheta);
    }
    else if (fabs(m_Lpi-(M_Xi-2.5*.008))<.004 ||
	     fabs(m_Lpi-(M_Xi+2.5*.008))<.004) {
      if ((s_dipc&0x1) && s_pcCa && (s_dics&0x1)) hn_d3hXi->Fill(dd3h);
      if ((s_dipc&0x2) && s_pcCa && (s_dics&0x2)) hn_D3hXi->Fill(d3h/dd3h);
      if (s_pcCa && (s_dics&0x1) && s_d3h)        hn_dpcXi->Fill(dipc);
      if (s_pcCa && (s_dics&0x2) && s_d3h)        hn_DpcXi->Fill(dipc/ddipc);
      if ((s_dipc&0x1) && (s_dics&0x1) && s_d3h)  hn_cpcXi->Fill(cthpcCa);
      if ((s_dipc&0x2) && (s_dics&0x2) && s_d3h)  hn_CpcXi->Fill(cthpcCa);
      if ((s_dipc&0x1) && s_pcCa && s_d3h)        hn_dcsXi->Fill(dics);
      if ((s_dipc&0x2) && s_pcCa && s_d3h)        hn_DcsXi->Fill(dics/ddics);
      if ((s_dipc&0x1) && s_pcCa && (s_dics&0x1) && s_d3h)
	hn_cpsLa->Fill(ctheta);
      if ((s_dipc&0x2) && s_pcCa && (s_dics&0x2) && s_d3h)
	hn_CpsLa->Fill(ctheta);
    }
    //                                             ********** uDST SELECTION: Xi
    if (m_Lpi<M_Xi+.41) // I.e. (little more than) "hm_Lpi" range
      uDSTSelection = true;

    //                                          ********** uDST SELECTION: Omega
    // (Not yet included, for the're incompatible w/ angle cut.) (!?)
    if (m_LK<M_Omega+.41) // I.e. (little more than) "hm_LK" range
      uDSTSelection = true;
  } // End loop over 3rd particle for Ksi and Omega search
  if (lvL) delete lvL;
}
// **********************************************************************
// ****************************** transport *******************************
// **********************************************************************
void transport(PaEvent &e,
	       PaTrack &trk,    // 
	       int &aR,         // RICH acceptance flag
	       int zones)
{
  // Transport argument track to
  // - RICH => Determine acceptance by RICH (if starting downstream of it).
  //  -> Fill argument "aR":
  //     =0: On input, meaning not yet evaluated
  //     =1: W/in RICH acceptance
  //     =2: Tighter cuts
  //     =-1: W/out acceptance
  //     =-2: Downstream of RICH
  //  -> Fill PaTrack's "vAux" [0,1] w/ (X,Y)
  //  -> Fill PaTrack's "vAux" [2,3] w/ (tgx,tgy)
  // - Downstream of SM2
  // => Determine exit (resp. incidence) angles for SM1 (resp. SM2), equated to
  // above @ RICH) and horizontal exit angle for SM2 (if involved).
  //  -> Track's "vAux": FIVE components if a 0x4-track, THREE otherwise.
  vector<float> &aux = trk.vAux(); aux.clear();
  static PaTPar Hout; // Generic helix used in extrapolations
  vector<PaTPar> &hs = trk.vTPar();
  if (hs[0](0)<ZRICH) { // Since starting downstream of RICH => not processed by RICH software...
    int ih, ihBest; double best; for (ih = 0, ihBest = 0, best = 10000;
				    ih<(int)hs.size(); ih++) {
      PaTPar &hR = hs[ih]; double Z0 = hR(0); if (Z0>ZSM2) continue;
      double dist = fabs(Z0-ZRICH); if (dist<best) {
	ihBest = ih; best = dist;
      }
    }
    const PaTPar &hR = hs[ihBest]; hR.Extrapolate(ZRICH,Hout,false);
    float XR = Hout(1), YR = Hout(2), RR = sqrt(XR*XR+YR*YR);
    if (rRICH<RR && fabs(XR)<XRICH && fabs(YR)<YRICH)
      // Fill all of "aRs", even if not useful, for it is computed only once,
      // cf. "aRs[iET]==0" conditioning the call to "transport".
      aR = RR>rRICHCut? 2 : 1;
    else
      aR = -1;
    aux.push_back(XR); aux.push_back(YR);
    aux.push_back(Hout(3)); aux.push_back(Hout(4));
  }
  else {
    aR = -2;
    aux.push_back(0); aux.push_back(0); aux.push_back(0); aux.push_back(0);
  }
  
  if (!(zones&0x4)) return;

  int ih, ihBest; double best; for (ih = 0, ihBest = -1, best = 10000;
				    ih<(int)hs.size(); ih++) {
    PaTPar &hS = hs[ih]; double Z0 = hS(0); if (Z0<ZSM2) continue;
    double dist = fabs(Z0-ZSM2+300); if (dist<best) {
      ihBest = ih; best = dist;
    }
  }
  if (ihBest<0) {
    // For &0x4-zone tracks, may happen if only the first helix is available...
    hs[0].Extrapolate(ZSM2+300/* i.e. 1m beyond SM2 exit */,Hout,false);
    aux.push_back(Hout(3));
  }
  else
    aux.push_back(hs[ihBest](3));
}
#ifdef U3_FILL_Lambda
// **********************************************************************
// *************************   bookLambdaMisc  **************************
// **********************************************************************
void bookLCascades(double LMassCut, int PIDScheme)
{
  const char *LaL[] = {"L","aL"};
  char LName[] =
    "hm_aLpiCac ";
  const char *LbL[] = {"#Lambda","#bar{#Lambda}"};
  char LTitle[] =
    "#bar{#Lambda}#pi+ - #bar{#Lambda}#pm5.2MeV,p:pID,e#piKVeto,#pi:eRejection,#pi:eRejection,Cascade  ";
  //"P#pi - p#pi=#bar{#Lambda}#pm5.2MeV,pID,eVeto)  ";
  size_t sT = strlen(LTitle)+1;
  char pID[]  = ",p:pID,e#piKVeto";
#  ifdef piS_e_REJECTION
  char piID[] = ",#pi:eRejection";
#  else
  char piID[] = "\0";
#endif

  for (int iLaL = 0; iLaL<2; iLaL++) {
    //                                                ***** hm_Lpi,KCa: Cascades
    int nBinsX = int((M_SigS+.4-M_Lam-M_pi+.05)/.002+.5);
    double xMn = M_SigS+.4-nBinsX*.002;
    sprintf(LName,"hm_%spiCa",LaL[iLaL]);
    snprintf(LTitle,sT,"%s#pi%c - %s#pm%.1fMeV%s%s%s%s,Cascade",
	    LbL[iLaL],iLaL?'+':'-',LbL[iLaL],LMassCut*1000,
	    PIDScheme?",#Lambda:":"\0",pID,piID,piID);
    hm_LpiCa[iLaL] =  new TH2D(LName,LTitle,nBinsX,xMn,M_SigS+.4,6,.5,6.5);
    sprintf(LName,"hm_%spiCa3",LaL[iLaL]);
    hm_LpiCa3[iLaL] = new TH2D(LName,LTitle,nBinsX,xMn,M_SigS+.4,6,.5,6.5);
    sprintf(LName,"hm_%spiCa6",LaL[iLaL]);
    hm_LpiCa6[iLaL] = new TH2D(LName,LTitle,nBinsX,xMn,M_SigS+.4,6,.5,6.5);
    sprintf(LName,"hm_%spiCaD",LaL[iLaL]);
    hm_LpiCaD[iLaL] = new TH2D(LName,LTitle,nBinsX,xMn,M_SigS+.4,6,.5,6.5);
    sprintf(LName,"hm_%spiCaC",LaL[iLaL]);
    hm_LpiCaC[iLaL] = new TH2D(LName,LTitle,nBinsX,xMn,M_SigS+.4,6,.5,6.5);
    nBinsX = int((M_Omega+.4-M_Lam-M_K+.05)/.002+.5);
    xMn = M_Omega+.4-nBinsX*.002;
    sprintf(LName,"hm_%sKCa",LaL[iLaL]);
    snprintf(LTitle,sT,"%sK%c - %s#pm%.1fMeV%s%s%s,K:ID,Cascade",
	    LbL[iLaL],iLaL?'+':'-',LbL[iLaL],LMassCut*1000,
	    PIDScheme?",#Lambda:":"\0",pID,piID);
    hm_LKCa[iLaL]  =  new TH2D(LName,LTitle,nBinsX,xMn,M_Omega+.4,6,.5,6.5);
    sprintf(LName,"hm_%sKCa3",LaL[iLaL]);
    hm_LKCa3[iLaL]  = new TH2D(LName,LTitle,nBinsX,xMn,M_Omega+.4,6,.5,6.5);
    sprintf(LName,"hm_%sKCa6",LaL[iLaL]);
    hm_LKCa6[iLaL]  = new TH2D(LName,LTitle,nBinsX,xMn,M_Omega+.4,6,.5,6.5);
    sprintf(LName,"hm_%sKCaD",LaL[iLaL]);
    hm_LKCaD[iLaL]  = new TH2D(LName,LTitle,nBinsX,xMn,M_Omega+.4,6,.5,6.5);
    sprintf(LName,"hm_%sKCaC",LaL[iLaL]);
    hm_LKCaC[iLaL]  = new TH2D(LName,LTitle,nBinsX,xMn,M_Omega+.4,6,.5,6.5);
  }
  // ********** CASCADE SELECTION HISTOS
  hs_d3hXi = new TH1D("hs_d3hXi","#Lambda#pi=#Xi+/-8MeV;d3h (cm)",   64, 0,1);
  hs_D3hXi = new TH1D("hs_D3hXi","#Lambda#pi=#Xi+/-8MeV;d3h/dd3h",   64, 0,16);
  hs_dpcXi = new TH1D("hs_dpcXi","#Lambda#pi=#Xi+/-8MeV;dipc (cm)",  128,0,160);
  hs_DpcXi = new TH1D("hs_DpcXi","#Lambda#pi=#Xi+/-8MeV;dipc/ddipc", 128,0,16);
  hs_cpcXi = new TH1D("hs_cpcXi","#Lambda#pi=#Xi+/-8MeV;cthpcCa",    100,.9995,1);
  hs_CpcXi = new TH1D("hs_CpcXi","#Lambda#pi=#Xi+/-8MeV;cthpcCa",    100,.9995,1);
  hs_dcsXi = new TH1D("hs_dcsXi","#Lambda#pi=#Xi+/-8MeV;dics (cm)",  128,0,160);
  hs_DcsXi = new TH1D("hs_DcsXi","#Lambda#pi=#Xi+/-8MeV;dics/ddics", 128,0,16);
  hs_cpsLa = new TH1D("hs_cpsLa","#Lambda#pi=#Xi+/-8MeV;cthpsLa",    100,.9995,1);
  hs_CpsLa = new TH1D("hs_CpsLa","#Lambda#pi=#Xi+/-8MeV;cthpsLa",    100,.9995,1);
  hn_d3hXi = new TH1D("hn_d3hXi","#Lambda#pi=#XiSideBand;d3h (cm)",  64, 0,1);
  hn_D3hXi = new TH1D("hn_D3hXi","#Lambda#pi=#XiSideBand;d3h/dd3h",  64, 0,16);
  hn_dpcXi = new TH1D("hn_dpcXi","#Lambda#pi=#XiSideBand;dipc (cm)", 128,0,160);
  hn_DpcXi = new TH1D("hn_DpcXi","#Lambda#pi=#XiSideBand;dipc/ddipc",128,0,16);
  hn_cpcXi = new TH1D("hn_cpcXi","#Lambda#pi=#XiSideBand;cthpcCa",   100,.9995,1);
  hn_CpcXi = new TH1D("hn_CpcXi","#Lambda#pi=#XiSideBand;cthpcCa",   100,.9995,1);
  hn_dcsXi = new TH1D("hn_dcsXi","#Lambda#pi=#XiSideBand;dics (cm)", 128,0,160);
  hn_DcsXi = new TH1D("hn_DcsXi","#Lambda#pi=#XiSideBand;dics/ddics",128,0,16);
  hn_cpsLa = new TH1D("hn_cpsLa","#Lambda#pi=#XiSideBand;cthpsLa",   100,.9995,1);
  hn_CpsLa = new TH1D("hn_CpsLa","#Lambda#pi=#XiSideBand;cthpsLa",   100,.9995,1);
}
void bookLambdaMisc(const double *V0DsdD, const double *V0cthCut, const double *V0pTCut,
		    double dM_Lambda, double LMassCut, double *xFbins,
		    int PIDScheme)
{
  double cthCut = V0cthCut[1]; int pTCut = int(V0pTCut[1]*1000+.5);
  gDirectory->cd("/"); TDirectory *baseDir = gDirectory;
  TDirectory *dLambda; if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda")))
		     gDirectory->mkdir("Lambda","Lambda Histos");
  if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda"))) {
    printf("\n** U3:\a No creating subdir \"Lambda\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dLambda->cd();

  const char *LaL[] = {"L","aL"};
  char LName[] =
    "hm_aLcthS0 ";
  //"hm_aLVsID ";
  //"hm_aLpiCac ";
  const char *LbL[] = {"#Lambda","#bar{#Lambda}"};
  char LTitle[] =
    "p+#pi- - D>2.0#deltaD,c#theta>0.99995,pT>100MeV,p:pID,e#piKVeto,#pi:eRejection  vs. c(#theta*) - .05<xF<.17  ";
  //"#Lambda#pm5.2MeV,D>2.0#deltaD,c#theta>0.99990,pT>100MeV,#piID;Pp (GeV);#ThetaC+d#ThetaC-#ThetaCp (mrd)  ";
  //"#bar{#Lambda}#pi+ - #bar{#Lambda}#pm5.2MeV,p:pID,e#piKVeto,#pi:eRejection  ";
  //"P#pi - p#pi=#bar{#Lambda}#pm5.2MeV,pID,eVeto)  ";
  size_t sT = strlen(LTitle)+1;
  char xFRange[] = " - .05<xF<.17";
  char pID[]  = ",p:pID,e#piKVeto";
#  ifdef piS_e_REJECTION
  char piID[] = ",#pi:eRejection";
#  else
  char piID[] = "\0";
#endif
  double DSdDCut = V0DsdD[1];

  // ***** Q2 and xBj BINING
  double sq16_10; sq16_10 = sqrt(sqrt(sqrt(sqrt(10.))));
  double Q2bin; int bin; double Q2bins[104];
  for (bin = 103, Q2bin = 100*sq16_10*sq16_10; bin>=0; bin--) {
    Q2bins[bin] = Q2bin; Q2bin /= sq16_10;
  }
  double xBbin; double xBbins[104];
  for (bin = 103, xBbin = 1*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10*sq16_10;
       bin>=0; bin--) {
    xBbins[bin] = xBbin; xBbin /= sq16_10;
  }

  // ***** BOOK
  // - hk_Pp|piLp  (hk_(p|pi)[ID])[a]L: p,pi MOMENTA for Lambda's W/ and W/O PID
  // - hm_LcthS    (hm_[a]LcthSi: cos(theta*) vs. Mass as a f(xF)
  // - hQ2,xB,y,xF (h(Q2|xB|yB|xF)_[a]L: Kinematics of Lambda's
  int nBinsX = int((M_Lam+dM_Lambda-M_p-M_pi+.005)/.0005+.5);
  double xMn = M_Lam+dM_Lambda-nBinsX*.0005;
  const char *decays[] = {"p","#pi"};
  for (int iLaL = 0; iLaL<2; iLaL++) {
    const char *ppi = decays[iLaL], *pip = decays[1-iLaL];
    for (int PID = 0; PID<2; PID++) {
      //                ***** hp_Lp,pi: p,pi MOMENTA for Lambda's W/ and W/O PID
      sprintf(LName,"hk_Pp%s%s",LaL[iLaL],PID?"ID":"\0");
      snprintf(LTitle,sT,"Pp - p#pi=%s#pm%.1fMeV%s%s",
	      LbL[iLaL],LMassCut*1000,PID?pID:"\0",PID?piID:"\0");
      hk_PpL[iLaL+PID*2]  = new TH1D(LName,LTitle,100,0,100);
      sprintf(LName,"hk_Ppi%s%s",LaL[iLaL],PID?"ID":"\0");
      snprintf(LTitle,sT,"P#pi - p#pi=%s#pm%.1fMeV%s%s",
	      LbL[iLaL],LMassCut*1000,PID?pID:"\0",PID?piID:"\0");
      hk_PpiL[iLaL+PID*2] = new TH1D(LName,LTitle,100,0,100);
      if (PID) {
	sprintf(LName,"hm_%sLID",LaL[iLaL]);
	snprintf(LTitle,sT,"D>%.1f/#deltaD,c#theta>%.5f,pT>%dMeV%s%s;M%s%s (GeV)",
		V0DsdD[1],cthCut,pTCut,pID,piID,ppi,pip);
	hm_LID[iLaL] = new TH1D(LName,LTitle,nBinsX,xMn,M_Lam+dM_Lambda);
      }
    }
    //               ***** hm_LcthS: COSTHETA* vs. ppi(pibar{p}) MASS as a f(xF)
    for (int ixF = 0; ixF<5; ixF++) {
      if (ixF==0) sprintf(xFRange,"%c",'\0');
      else        sprintf(xFRange," - .%02d<xF<.%02d",
			  int(xFbins[ixF-1]*100+.5),int(xFbins[ixF]*100+.5));
      if (iLaL==0) {
	if (ixF==0) sprintf(LName,"hm_LcthS%c",'\0');
	else        sprintf(LName,"hm_LcthS%d",ixF-1);
	snprintf(LTitle,sT,"%s - D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV%s%s;M%s%s (GeV);(#theta*)",
		xFRange,DSdDCut,cthCut,pTCut,pID,piID,ppi,pip);
	hm_LcthS[ixF]  = new TH2D
	  (LName,LTitle,nBinsX,xMn,M_Lam+dM_Lambda,16,-1,1);
      }
      else {
	if (ixF==0) sprintf(LName,"hm_aLcthS%c",'\0');
	else        sprintf(LName,"hm_aLcthS%d",ixF-1);
	snprintf(LTitle,sT,"%s - D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV%s%s;M%s%s (GeV);c(#theta*)",
		xFRange,DSdDCut,cthCut,pTCut,pID,piID,ppi,pip);
	hm_aLcthS[ixF] = new TH2D
	  (LName,LTitle,nBinsX,xMn,M_Lam+dM_Lambda,16,-1,1);
      }
    }

    // "Q2 - #bar{#Lambda}#pm5.2MeV,..."   // ***** KINEMATICS for Lambda's
    sprintf(LName,"hQ2_%s",LaL[iLaL]);
    sprintf(LTitle,"Q2 - %s#pm%.1fMeV%s%s",
	    LbL[iLaL],LMassCut*1000,pID,piID);
    hQ2_L[iLaL] = new TH1D(LName,LTitle,   103,Q2bins);
    sprintf(LName,"hxB_%s",LaL[iLaL]);
    sprintf(LTitle,"xB - %s#pm%.1fMeV%s%s",
	    LbL[iLaL],LMassCut*1000,pID,piID);
    hxB_L[iLaL] = new TH1D(LName,LTitle,   103,xBbins);
    sprintf(LName,"hyB_%s",LaL[iLaL]);
    sprintf(LTitle,"y - %s#pm%.1fMeV%s%s",
	    LbL[iLaL],LMassCut*1000,pID,piID);
    hyB_L[iLaL] = new TH1D(LName,LTitle,   110,-.05,1.05);
    sprintf(LName,"hxF_%s",LaL[iLaL]);
    sprintf(LTitle,"x_{F} - %s#pm%.1fMeV%s%s",
	    LbL[iLaL],LMassCut*1000,pID,piID);
    hxF_L[iLaL] = new TH1D(LName,LTitle,   200, -1, 1);
    sprintf(LName,"hv_pZ_%s",LaL[iLaL]);
    sprintf(LTitle,"pVertex Z - %s#pm%.1fMeV%s%s",
	    LbL[iLaL],LMassCut*1000,pID,piID);
    getFIMMAbscissae();
    double ZMn, ZMx; int nZbins; getTargetBinning(nZbins,ZMn,ZMx);
    hv_pZ_L[iLaL] = new TH1D(LName,LTitle, nZbins,ZMn,ZMx);
  }

  // ***** BOOK
  // - hm_LVsID    (hm_[a]LVsID: Mass vs. p/pi-ID
  char cID[] =
    "pID:0x1=1.05#piThr<P<100GeV 0x2=P<1.05pTHr,Else<0.92Back 0x4=p>1.20Back,1.01Else,"
    "#piID:0x8=#piThr<P<60GeV 0x10=e<1.50#pi 0x20=#pi>1.05Back,1.01Hadron;Mp#pi";
  size_t size = strlen(cID)+1;
  for (int iLaL = 0; iLaL<2; iLaL++) {
    const char *ppi = decays[iLaL], *pip = decays[1-iLaL];
    sprintf(LName,"hm_%sVsID",LaL[iLaL]);
    if (snprintf(cID,size,
    "pID:0x1=%.2f#piThr<P<%.0fGeV 0x2=P<%.2fpTHr,Else<%.2fBack 0x4=p>%.2fBack,%.2fElse,"
    "#piID:0x8=#piThr<P<%.0fGeV 0x10=e<%.2f#pi 0x20=#pi>%.2fBack,%.2fHadron;M%s%s",
       subpThrPmin,pIDPmax,pIDPmin,subpThrLHVeto,LHBckCut,LHCut,
		 PRICHCut,piLHeVeto,LHBckCut,LHCut,ppi,pip)>=(int)size) {
      printf("\n** bookLambdaMisc:\a Error writing string \"%s...\"\n\n",cID);
      assert(false);
    }
    hm_LVsID[iLaL] = new TH2D(LName,cID,nBinsX,xMn,M_Lam+dM_Lambda,64,-0.5,63.5);
  }

  // (X,Y), THETA vs. P @ RICH
  // - hk_X|YL     (hk_(X|Y)(p|pi)L): X/Y @ RICH (-/+ along the verticalaxis)
  // - hk_PThL     (hk_PTh(p|pi)L): theta vs. P @ RICH
  for (int ippi = 0; ippi<2; ippi++) {
    const char *ppi = decays[ippi], *pip = decays[1-ippi];
    sprintf(LName,"hk_X%sL",ippi?"pi":"p");
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,%s#mpID;XRICH (cm);%s#pm",LMassCut*1000,DSdDCut,cthCut,pTCut,pip,ppi);
    hk_XL[ippi] = new TH2D(LName,LTitle,160,-160,160,2,-2,2);
    LName[3] = 'Y';
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,%s#mpID;YRICH (cm);%s#pm",LMassCut*1000,DSdDCut,cthCut,pTCut,pip,ppi);
    hk_YL[ippi] = new TH2D(LName,LTitle,160,-160,160,2,-2,2);
    sprintf(LName,"hk_PTh%sL",ippi?"pi":"p");
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,%s#mpID;P%s#pm;#thetaRICH (rd)",LMassCut*1000,DSdDCut,cthCut,pTCut,pip,ppi);
    hk_PThL[ippi] = new TH2D(LName,LTitle,100,0,100,80,0,0.8);
    if (ippi) continue;
    //                                               ***** THETA CHERENKOV vs. P
    // N.B.: Stricter cthCut
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,#piID;Pp (GeV);#ThetaC (mrd)",         LMassCut*1000,V0DsdD[2],V0cthCut[2],pTCut);
    hR_thCPL = new TH2D("hR_thCPL",LTitle,320,-80,80,240,0,60);
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,#piID;Pp (GeV);#ThetaC+d#ThetaC (mrd)",LMassCut*1000,V0DsdD[2],V0cthCut[2],pTCut);
    hR_ThCPL = new TH2D("hR_ThCPL",LTitle,320,-80,80,240,0,60);
    snprintf(LTitle,sT,"#Lambda#pm%.1fMeV,D>%.1f#deltaD,c#theta>%.5f,pT>%dMeV,#piID;Pp (GeV);#ThetaC+d#ThetaC (mrd)",LMassCut*1000,V0DsdD[2],V0cthCut[2],pTCut);
    hR_THCPL = new TH2D("hR_THCPL",LTitle,320,-80,80,128,-4,4);
  }

  // ***** BOOK
  // - hm_Lpi     (hm_[a]Lpi(m|p)[ID]): Sigma*
  // - hm_Lpi,KCa (hm_[a]LpiCa): Cascade (Xi, Omega)
  TDirectory *dLX; if (!(dLX = (TDirectory*)gDirectory->Get("LambdaX")))
    gDirectory->mkdir("LambdaX","Sigma*, Xsi, Omega");
  if (!(dLX = (TDirectory*)gDirectory->Get("LambdaX"))) {
    printf("\n** U3:\a No creating subdir \"K0X\" in TDirectory \"%s/Lambda\"\n\n",
	   baseDir->GetName()); assert(false);
  }
  dLX->cd();

  for (int iLaL = 0; iLaL<2; iLaL++) {
    //                                                      ***** hm_Lpi: Sigma*
    int nBinsX = int((M_SigS+.4-M_Lam-M_pi+.05)/.002+.5);
    double xMn = M_SigS+.4-nBinsX*.002;
    for (int pm = 0; pm<2; pm++) {
      sprintf(LName,"hm_%spi%c",LaL[iLaL],pm?'m':'p');
      sprintf(LTitle,"%s#pi%c - %s#pm%.1fMeV%s%s%s",
	      LbL[iLaL],pm?'-':'+',LbL[iLaL],LMassCut*1000,
	      PIDScheme?",#Lambda:":"\0",pID,piID);
      hm_Lpi[iLaL+2*pm] = new TH1D(LName,LTitle,nBinsX,xMn,M_SigS+.4);
      sprintf(LName,"hm_%spi%cID",LaL[iLaL],pm?'m':'p');
      sprintf(LTitle,"%s#pi%c - %s#pm%.1fMeV%s%s%s,#pi:ID",
	      LbL[iLaL],pm?'-':'+',LbL[iLaL],LMassCut*1000,
	      PIDScheme?",#Lambda:":"\0",pID,piID);
      hm_LpiID[iLaL+2*pm] = new TH1D(LName,LTitle,nBinsX,xMn,M_SigS+.4);
    }
  }
  bookLCascades(LMassCut,PIDScheme);
}
// **********************************************************************
// *************************  printLambdaCuts  **************************
// **********************************************************************
void printLambdaCuts(double DsdDCut, double cthCut, double pTCut, double vChi2Cut)
{
  printf("\n * U3_FILL_Lambda:\n V0: D/D,cos(theta),pT,vChi2 Event: Q2,yMn,yMx,Tgt");
#  ifdef Lambda_p_ID
  printf(" pID");
#  endif
#  ifdef piS_e_REJECTION
  printf(" piSeVeto: eLH/pi,bckLH");
#  endif
  printf("\n");
  printf(" %.1f,%.5f,%.3f,%.1f   %.3f",DsdDCut,cthCut,pTCut,vChi2Cut,Q2Cut);
#  ifdef U3_y_CUT
  printf(",%4.2f,%4.2f", yLowCut,yUpCut);
#  endif
# if U3_TARGET_CUT > 1
  printf(",Target#%d",U3_TARGET_CUT);
#  endif
  printf("\n");
#  ifdef Lambda_p_ID
  printf(" pID: pMn/pThr,pMx,LH,bckLH Veto: pMn/piThr,piLH,eKLH ");
  printf("   %4.2f,%.0f,%.2f,%.2f",pIDPmin,pIDPmax,LHCut,LHBckCut);
  printf(" %4.2f,%4.2f,%4.2f",
	 subpThrPmin,subpThrLHpiVeto,subpThrLHVeto);
#  endif
#  ifdef piS_e_REJECTION
  printf(" %4.2f",piLHeVeto);
#  endif
  printf("\n");
}
#endif
#if defined U3_K0_HIGHLEVEL || defined U3_L_HIGHLEVEL
// **********************************************************************
// *************************      getGM0709WB  **************************
// *************************   fillK0GM0709WB  **************************
// *************************    fillLGM0709WB  **************************
// **********************************************************************
void getGM0709WB(int *gm0709Ws, unsigned int *gm0709Bs)
{
  // ***** SPECIAL HIGH LEVEL HISTO for (L?)SAS *****
  // ***** => DETERMINE (WORD,BITS) ASSOCIATED TO GM07:09
  //  - SAS SAT = GEM07:09: require all 3 stations
  //  - SAS LAT = GEM07:09: require there be none.

  // "gm0709" = words storing the bitwise description of GM07:09's. (Since these
  // are less than 32 units apart in the detector list, 2 words must suffice.)
  static int gm0709Ws[2]; static unsigned int gm0709Bs[2];

  static bool first = true; if (first) {
    first = false;
    const PaSetup &setup = PaSetup::Ref();
    for (int i = 0; i<2; i++) { gm0709Ws[i] = -1; gm0709Bs[i] = 0; }
    char coords[] = "XYVU", charName[] = "GM07X1__";
    int iw, idet; for (iw=idet = 0; iw<HIT_MAP_SIZE; iw++) {
      for (int ib = 0; ib<32; ib++, idet++) {
	if (idet>=setup.NDetectors()) continue;
	const string &name = setup.Detector(idet).Name();
	if (name.find("GM")) continue;
	for (int coord = 0; coord<4; coord++)
	  for (int station = 7; station<=9; station++) {
	    sprintf(charName,"GM0%d%c1__",station,coords[coord]);
	    string planeName(charName);
	    if (name==planeName) {
	      int igmW, match; for (igmW=match = 0; igmW<2; igmW++) {
		if (gm0709Ws[igmW]<0 || gm0709Ws[igmW]==iw) {
		  gm0709Ws[igmW] = iw; gm0709Bs[igmW] |= 1<<ib;
		  match = 1; break;
		}
	      }
	      if (!match) {
  printf("\n** U3: No match for word %d(=%s) in gm0709 words = (%d,%d)\n\n",
	 iw,charName,gm0709Ws[0],gm0709Ws[1]); assert(false);
	      }
	    }
	  }
      }
    }
    if (gm0709Ws[0]<0 || gm0709Ws[1]!=-1 && gm0709Ws[1]!=gm0709Ws[0]+1) {
      printf("\n** U3: Inconsistent gm0709 words = (%d,%d)\n\n",
	     gm0709Ws[0],gm0709Ws[1]); assert(false);
    }
  }
}
#  ifdef U3_K0_HIGHLEVEL
void fillK0GM0709(double m_pipi, const PaTrack &trk1, const PaTrack &trk2,
		  const int *gm0709Ws, const unsigned int *gm0709Bs)
{
  const PaTrack *t; int pm, gm; for (pm=gm = 0, t = &trk1; pm<2; pm++) {
    const unsigned int *found = t->FoundHitsBitMap();
    int iw, nGMs; for (iw = 0, nGMs = 0; iw<2; iw++) {
      int word = gm0709Ws[iw]; if (word==-1) continue;
      unsigned int gms = found[word]&gm0709Bs[iw];
      for (unsigned int ib = 1; ib<=0x10000000; ib <<= 1)
	if (ib&gms) nGMs++;
    }
    if (nGMs>7) gm |= 1<<pm; t = &trk2;
  }
  if      (gm!=0) hm_K0cssS->Fill(m_pipi);
  else if (gm==0) hm_K0cssL->Fill(m_pipi);
}
#  endif
#  ifdef U3_L_HIGHLEVEL
void fillLGM0709(double m_ppi, const PaTrack &trk1, const PaTrack &trk2,
		 const int *gm0709Ws, const unsigned int *gm0709Bs)
{
  const PaTrack *t; int pm, gm; for (pm=gm = 0, t = &trk1; pm<2; pm++) {
    const unsigned int *found = t->FoundHitsBitMap();
    int iw, nGMs; for (iw = 0, nGMs = 0; iw<2; iw++) {
      int word = gm0709Ws[iw]; if (word==-1) continue;
      unsigned int gms = found[word]&gm0709Bs[iw];
      for (unsigned int ib = 1; ib<=0x10000000; ib <<= 1)
	if (ib<gms) nGMs++;
    }
    if (nGMs>9) gm |= 1<<pm; t = &trk2;
  }
  if      (gm!=0) FillLambda(m_ppi,0x10);
  else if (gm==0) FillLambda(m_ppi,0x20);
}
#  endif
#endif
#ifdef U3_K0_HIGHLEVEL  // K0 vs. Incidence and related
// **********************************************************************
// *************************  fillK0HighLevel  **************************
// **********************************************************************
void fillK0HighLevel(const PaTrack &trk1, const PaTrack &trk2,
		     double m_pipi, const TLorentzVector &lvpipi,
		     double Theta, double Thetx, double Thety, double ThetX,
		     double alpha, double beta)
{
  double zL1 = trk1.ZLast(), zL2 = trk2.ZLast();

  if      (zL1<ZSM2 && zL2<ZSM2)                                  // ***** LAS^2
    fillK0IncidenceLAS(m_pipi,lvpipi,
		       Theta,Thetx,Thety,ThetX,alpha,beta);
  else if (zL1<ZSM2 || zL2<ZSM2<0) {                            // ***** LASxSAS
    if (zL1>ZSM2) {
      double alpha1 = trk1.vAux()[2], bin = floor(alpha1/.020)+2;
      if (bin<0) bin = 0; if (bin>3) bin = 3;
      hm_K0cslx->Fill(m_pipi,bin);
    }
    else {
      double alpha2 = trk2.vAux()[2], bin = floor(alpha2/.020)+1;
      if (bin<0) bin = 0; if (bin>3) bin = 3;
      hm_K0clsx->Fill(m_pipi,bin);
    }
  }
  else                                                            // ***** SAS^2
    fillK0IncidenceSAS(m_pipi,lvpipi,hv1,hv2,
		       Theta,Thetx,Thety,ThetX,alpha,beta);

  //                             No LAS: tracking only downstream of SM1
  const PaTPar &hi1 = trk1.vTPar()[0], &hi2 = trk2.vTPar()[0];
  const PaTPar *h; int pm, las;
  for (pm=las = 0, h = &hi1; pm<2; pm++) {
    if ((*h)(0)<350) las |= 1<<pm; h = &hi2;
  }
  if (!las) {
    const int *gm0709Ws; const unsigned int *gm0709Bs;
    getGM0709WB(gm0709Ws,gm0709Bs);
    hm_K0cnn->Fill(m_pipi);
    if (hi1(5)<.020)  hm_K0cNNp->Fill(m_pipi);
    if (hi2(5)>-.020) hm_K0cNNm->Fill(m_pipi);
    // SAT or LAT at the exit of SM2?
    fillK0GM0709(m_pipi,trk1,trk2,gm0709Ws,gm0709Bs);
  }
}
// **********************************************************************
// *************************   getK0Incidence  **************************
// *************************  fillK0Incidence  **************************
// **********************************************************************
void getK0Incidence(const PaTrack &trk1, const PaTrack &trk2,
		    int zones1, int zones2,
		    const TVector3 &v13, const TVector3 &v23,
		    const PaTPar &hi1, const PaTPar &hi2,
		    double &Theta, double &Thetx, double &Thety, double &ThetX,
		    double &alpha, double &beta)
{
  //             ***** INCIDENCE ANGLE (and related) *****
  double cTh = v13*v23/v13.Mag()/v23.Mag();
  if (cTh>1.) cTh = 1.;
  Theta = acos(cTh);
  double alpha1 = (zones1&0x4) ? trk1.vAux()[2] : hi1(3);
  double alpha2 = (zones2&0x4) ? trk2.vAux()[2] : hi2(3);
  Thetx = alpha2-alpha1, alpha = (alpha2+alpha1)/2;
  alpha1 = (zones1&0x4) ? trk1.vAux()[3] : hi1(4);
  alpha2 = (zones2&0x4) ? trk2.vAux()[3] : hi2(4);
  Thety = alpha2-alpha1;
  alpha1 = (zones1&0x4) ? trk1.vAux()[4] : trk1.vAux()[2];
  alpha2 = (zones2&0x4) ? trk2.vAux()[4] : trk2.vAux()[2];
  ThetX = (alpha2-alpha1)/2, beta = (alpha2+alpha1)/2;
}
void fillK0IncidenceLAS(double m_pipi, const TLorentzVector &lvpipi,
			double Theta, double Thetx, double Thety, double ThetX,
			double alpha, double beta)
{
  if      (Theta<.05) hm_K0cllT->Fill(m_pipi,0);        // vs. Theta
  else if (Theta<.1)  hm_K0cllT->Fill(m_pipi,1);
  else                hm_K0cllT->Fill(m_pipi,2);
  double bin = floor(Thetx/.05)+2;
  if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0cllx->Fill(m_pipi,bin);                          // vs. Theta X
  bin = floor(Thety/.05)+2; if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0clly->Fill(m_pipi,bin);                          // vs. Theta Y
  bin = floor(alpha/.025)+2; if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0clla->Fill(m_pipi,bin);                          // vs. alpha
  bin = floor(ThetX/.05); if (bin<0) bin = 0; if (bin>5) bin = 5;
  hm_K0cllX->Fill(m_pipi,bin);                          // vs. Beta X
  bin = floor(beta/.08)+3; if (bin<0) bin = 0; if (bin>5) bin = 5;
  hm_K0cllb->Fill(m_pipi,bin);                          // vs. Beta X
  /*
    printf("\nEvt %d %02d,%02d 0x%x,0x%x %5.1f,%5.1f %.4f  %6.3f,%6.3f  %6.2f => %.0f    %.4f %.4f\n",
    EvNum,iET1,iET2,zones1,zones2,1/hv1(5),1/hv2(5),m_pipi-M_K0,alpha1,alpha2,Thetx,floor(Thetx/.05)+2,fSM1*(1-dSM1dx*(hi1(3)+.025)),fSM1*(1-dSM1dx*(hi2(3)+.025)));
  */
  double pK0 = lvpipi.P(); bin = floor(pK0/12); if (bin>5) bin = 5;
  hm_K0cllP->Fill(m_pipi,bin);
}
void fillK0IncidenceSAS(double m_pipi, const TLorentzVector &lvpipi,
			const PaTPar &hv1, const PaTPar &hv2,
			double Theta, double Thetx, double Thety, double ThetX,
			double alpha, double beta)
{
  if      (Theta<.01)  hm_K0cssT->Fill(m_pipi,0);       // vs. Theta
  else if (Theta<.025) hm_K0cssT->Fill(m_pipi,1);
  else                 hm_K0cssT->Fill(m_pipi,2);
  double bin = floor(Thetx/.025)+1;
  if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0cssx->Fill(m_pipi,bin);                          // vs. Theta X
  bin = floor(Thety/.015)+2; if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0cssy->Fill(m_pipi,bin);                          // vs. Theta Y
  bin = floor(alpha/.01)+2; if (bin<0) bin = 0; if (bin>3) bin = 3;
  hm_K0cssa->Fill(m_pipi,bin);                          // vs. alpha
  bin = floor(ThetX/.025); if (bin<0) bin = 0; if (bin>5) bin = 5;
  hm_K0cssX->Fill(m_pipi,bin);                          // vs. Beta X
  bin = floor(beta/.03)+3; if (bin<0) bin = 0; if (bin>5) bin = 5;
  hm_K0cssb->Fill(m_pipi,bin);                          // vs. Beta X
  /*
    if (fabs(m_pipi-M_K0)<.03)
    printf("\nEvt %d %02d,%02d 0x%x,0x%x %5.1f,%5.1f %.4f   %6.3f,%6.3f  %6.3f => %.0f    %.4f %.4f\n",
    EvNum,iET1,iET2,zones1,zones2,1/hv1(5),1/hv2(5),m_pipi-M_K0,alpha1,alpha2,Thetx,floor(Thetx/.025)+3,fSM2*(1-dSM2dx*trk1.vAux()[2]),fSM2*(1+dSM2dx*trk2.vAux()[2]));
  */
  double pK0 = lvpipi.P(); bin = floor(pK0/12)-2; if (bin>5) bin = 6;
  hm_K0cssP->Fill(m_pipi,bin);
  if (hv1(5)<.020)  hm_K0cSSp->Fill(m_pipi);
  if (hv2(5)>-.020) hm_K0cSSm->Fill(m_pipi);
  if (hv1(5)<.0167)  hm_K0cSSP->Fill(m_pipi);
  if (hv2(5)>-.0167) hm_K0cSSM->Fill(m_pipi);
}
#endif
#ifdef U3_OUTPUT_TREE
// **********************************************************************
// ************************* copyCSEventHeader **************************
// ************************* copyCSHadronData  **************************
// ************************* copyCSResonanceData ************************
// ************************* getHadronIndex    **************************
// **********************************************************************
void getRICHMultiplicity(PaEvent &e, const PaVertex &pV,
			 int &nTrksRIt, int &nTrksRIb);
void copyCSEventHeader(PaEvent &e,
		       unsigned short primaryPat, unsigned short spillOKPat,
		       const PaVertex &pV,
		       double E0, double Q2, double xB, double y,
		       double index, double prodIndex)
{
  fCSEvt->runNo      = e.RunNum();
  fCSEvt->spillNo    = e.SpillNum();
  fCSEvt->evtNo      = (int)e.UniqueEvNum();
  fCSEvt->trigMask   = e.TrigMask();
  fCSEvt->primaryPat = primaryPat;
  fCSEvt->spillOKPat = spillOKPat;

  fCSEvt->Xp = pV.X(); fCSEvt->Yp = pV.Y(); fCSEvt->Zp = pV.Z();

  fCSEvt->E0 = E0; fCSEvt->Q2 = Q2; fCSEvt->y = y;

  fCSEvt->index = index-1; fCSEvt->prodIndex = prodIndex-1; fCSEvt->piThr = piThr;
  fCSEvt->dEK = 100;  // A large number, that won't pass any exclusivity
  fCSEvt->nOuts = pV.NOutParticles();
  getRICHMultiplicity(e,pV,fCSEvt->nTrksRIt,fCSEvt->nTrksRIb);
}
void copyCSHadronData(PaEvent &e, const PaParticle *pa, int iV)
{
  int iET = pa->iTrack(); PaTrack &trk = e.vTrack()[iET];
  const vector<Float_t> &aux = trk.vAux();
  if (aRs[iET]<-1) { // Ensure then "aRs" are meaningful
    printf("** U3:\a Evt %d,%d Track %d(<-Hadron): Inconsistency: aRs = %d\n",
	   e.RunNum(),(int)e.UniqueEvNum(),iET,aRs[iET]);
    abort();
  }
  CSHadronData h;
  // 3-momentum @ pV
  const PaTPar &hv = pa->ParInVtx(iV); const TVector3 v3 = hv.Mom3();
  h.Px = v3.X(); h.Py = v3.Y(); h.Pz = v3.Z();
  // P is taken as close to RICH as possible, yet upstream of SM2
  const vector<PaTPar> &hs = trk.vTPar();
  int i, iBest; double best; for (i = 0, iBest = 0, best = 10000;
				  i<(int)hs.size(); i++) {
    const PaTPar &hR = hs[i]; double Z0 = hR(0); if (Z0>ZSM2) continue;
    double dist = fabs(Z0-ZRICH); if (dist<best) { iBest = i; best = dist; }
  }
  const PaTPar &hR = hs[iBest]; // At worst, "iBest" is =0, i.e. 1rst measured point
  h.qP = 1/hR(5); h.phiR = atan2(hR(4),hR(3));
  h.XX0 = trk.XX0();
  h.ZFirst = hs[0](0);
  h.ZLast = trk.ZLast(); // Get ZLast from PaTrack dedicated method, even if it's not the fastest way, given that not that many tracks are concerned. 
  h.chi2 = trk.Chi2tot()/trk.Ndf();
  for (int i = 0; i<6; ++i) {
    float LH = fPid->GetLike(i,trk);
    if (fabs(LH)<1e-6)// It can happen that the LHs come close to 0...
      // ... this behaviour is prone to make the analysis depend upo context.
      // => Let's set =0 all such cases (as is done in "../GetPIDS.cc")
      LH = 0;
    h.LH[i] = LH;
  }
  if (trk.NRichInf()>0) {
    for (int i = 0; i<3; ++i) h.dLHdI[i] = trk.RichInf(i+4);
    if (trk.RichInf(0)) {
      h.hasR = true; h.thC = trk.RichInf(7);
    }
  }
  h.XR = aux[0]; h.YR = aux[1]; h.tgXR = aux[2]; h.tgYR = aux[3];

  // ***** Get electromagnetic calorimeter signals
  int iC; for (iC = 0, h.ECAL = 0; iC<pa->NCalorim(); iC++) {
    const PaCaloClus &clus = e.vCaloClus(pa->iCalorim(iC));
    if (clus.CalorimName()[0]=='E') // Only ECALs
      h.ECAL += clus.E();
  }
  if (e.IsMC()) {
    int iMCT = trk.iMCtrack(); if (iMCT>=0) h.MCpid = e.vMCtrack(iMCT).Pid();
  }
  fHadrons.push_back(h);
}
void copyCSResonanceData(const PaVertex &v,
			 double m, double D, double dD, double cth,
			 double pT, double alpha)
{
  CSResonanceData r;
  r.Xs = v.Pos(0); r.Ys = v.Pos(1); r.Zs = v.Pos(2);
  r.chi2 = v.Chi2()/v.Ndf();
  r.m = m; r.D = D; r.dD = dD; r. cth = cth; r.pT = pT; r.alpha = alpha;
  fResonances.push_back(r);
}
int getHadronIndex(const PaParticle *pa) {
  int paInd = pa->MyIndex();
  map<int,int>::iterator it = fPa2fH.find(paInd);
  int fHInd = it!=fPa2fH.end() ? it->second : -1;
  if (fHInd<0) {
    fHInd = fHadrons.size();
    if (fPa2fH.size()!=fHadrons.size()) {
      printf("** U3::getHadronIndex: Evt %d#%d fHadrons(%d) != fPa2fH(%d)\n",
	     Run,EvNum,(int)fHadrons.size(),(int)fPa2fH.size());
      abort();
    }
    fPa2fH[paInd] = fHInd;
  }
  //#define U3_DEBUG_fPa2fH
#ifdef U3_DEBUG_fPa2fH
  if (fHInd!=(int)fHadrons.size()) {
    int paInd = pa->MyIndex(), iET = pa->iTrack(), nfHs = fHadrons.size();
    const PaTrack &trk = PaEvent::Ptr()->vTrack()[iET];
    const PaTPar &hi = trk.vTPar()[0]; double qP = 1/hi(5);
    printf("U3_DEBUG_fPa2fH: Evt %d#%d pa#%d(%.2fGeV) fHInd %d/%d\n",
	   Run,EvNum,paInd,qP,fHInd,nfHs);
    for (int i = 0; i<nfHs; i++) {
      const CSHadronData &h = fHadrons[i]; printf("fH %d(%.2fGeV) ",i,h.qP);
    }
    printf("\n");
  }
#endif
  return fHInd;
}
#endif
// **********************************************************************
// *************************  checkExclusive   **************************
// **********************************************************************
void checkExclusive(PaEvent &e, int ipV, int imuS,
		    const PaParticle *&pa1, const PaParticle *&pa2, const PaParticle *&pa3,
		    unsigned short &exclPhiPat)
{
  // - Check # of 2ndaries: must be =3, allowance for an extra, slow, 4th.
  // - If =3, set "exclPhiPat = 0x1".
  // - Check there's a scattered mu, else abort.
  // - Returns the two candidate hadrons. Ordered (h+,h-1), if opposite signs.
  const PaVertex &pV = e.vVertex(ipV);
  int nTrksPV = pV.NOutParticles(), iVP, iVPSlow = -1;
#ifdef U3_QUASIEXCL_PHI
  if (nTrksPV==4) {
    // Case: 3 h-tracks: single out slowest one (! mu'); to be later disregarded.
    double pMn = 0; for (iVP = 0; iVP<nTrksPV; iVP++) {
      int iEP = pV.iOutParticle(iVP); if (iEP==imuS) continue;
      const PaParticle &pa = e.vParticle(iEP);
      const PaTPar &hv = pa.ParInVtx(ipV); double p = hv.Mom();
      if (iVPSlow<0 || p<pMn) { pMn = p; iVPSlow = iVP; }
    }
    pa3 = &e.vParticle(pV.iOutParticle(iVPSlow));
  }
  else if (nTrksPV!=3) {
    printf("\n** U3:\a Run %d Evt %d: Exclusivity Inconsistency:\"U3_QUASIEXCL_PHI\" defined, but #tracks in pV = %d\n",e.RunNum(),(int)e.UniqueEvNum(),nTrksPV);
      assert(false);
  }
  else exclPhiPat |= 0x1;
#else
  if (nTrksPV!=3) {
    printf("\n** U3:\a Run %d Evt %d: Exclusivity Inconsistency: #tracks in pV =%d, while \"U3_QUASIEXCL_PHI\" undefined\n",e.RunNum(),(int)e.UniqueEvNum(),nTrksPV);
    assert(false);
  }
  pa3 = 0;
  exclPhiPat |= 0x1;
#endif
  int iVPM, iVP1, iVP2; for (iVP = 0, iVPM=iVP1=iVP2 = -1; iVP<nTrksPV; iVP++) {
    if (iVP==iVPSlow) continue;  // Skip slow, extra, particle
    int iEP = pV.iOutParticle(iVP);
    if      (iEP==imuS) iVPM = iVP;
    else if (iVP1<0) iVP1 = iVP;
    else iVP2 = iVP;
  }
  if (iVPM<0 || iVP1<0 || iVP2<0) {
    printf("\n** U3:\a Run %d Evt %d: Exclusivity Inconsistency:\n",
	   e.RunNum(),(int)e.UniqueEvNum());
    printf(  "  mu' = %d",imuS);
    if (iVPSlow>=0) printf(", slow = %d",pV.iOutParticle(iVPSlow));
    printf("\n  %d OutParticles =",nTrksPV);
    for (iVP = 0; iVP<nTrksPV; iVP++)
      printf(" %d",pV.iOutParticle(iVP));
    printf("\n\n");
    assert(false);
  }
  pa1 = &e.vParticle(pV.iOutParticle(iVP1));
  pa2 = &e.vParticle(pV.iOutParticle(iVP2));
  int pairSigns = pa1->Q()+pa2->Q();
  if (pairSigns<-2 || 2<pairSigns || abs(pairSigns)==1) {
    printf("\n** U3:\a Run %d Evt %d: Exclusivity Inconsistency: h charges = %d,%d\n\n",
	   e.RunNum(),(int)e.UniqueEvNum(),pa1->Q(),pa2->Q()); assert(false);
  }
  if (pairSigns==0) {
    exclPhiPat |= 0x2;
    if (pa1->Q()<pa2->Q()) swap(pa1,pa2);       // First particle = positive
  }
}
// **********************************************************************
// *************************  detachedHadrons  **************************
// **********************************************************************
bool detachedHadrons(PaEvent &e, const PaVertex &pV)
{
  // Returns true if any, good, 2ndary track out of vertex, yet compatible w/ it.
  // Requirements for the tracks:
  // - Cut on chi2 and NDF: to filter out unreliable tracks
  // - Timed (since tracks w/o time measurement, and thus w/ hits only in drift
  //  detectors) are typically either halo tracks (and hence off-time) or low P
  //  tracks (and hence, typically, off-RICH).
  // - Time compatible w/ event Time.
  // - Not fully downstream of SM1
  // - Not paraxial
  // - Compatible w/ pVertex

  const PaParticle &beam = e.vParticle(pV.InParticle());
  const PaTrack &beamTrk = e.vTrack()[beam.iTrack()];
  double t0 = beamTrk.MeanTime(), dt0 = beamTrk.SigmaTime();

  int nTrks = e.vTrack().size(), nPats = nTrks/32+1, nTrksPV = pV.NOutParticles();
  unsigned int *tPats = new unsigned int[nPats];
  memset((void*)tPats,0,nPats*sizeof(unsigned int));
  for (int iEP = 0; iEP<nTrksPV; iEP++) {
      const PaParticle &p = e.vParticle(pV.iOutParticle(iEP));
      int iET = p.iTrack(); tPats[iET/32] |= 1<<iET%32;
  }
  double Zp = pV.Z(), Xp = pV.X(), Yp = pV.Y(); int iET, match; PaTPar H;
  for (iET=match = 0; iET<nTrks; iET++) {
    if ((1<<iET%32)&tPats[iET/32]) continue;// Skip tracks belonging to pV
    const PaTrack &trk = e.vTrack()[iET];

    if (trk.Ndf()<TRACK_NDF_MIN)            // Skip few hit tracks
      continue;
    if (trk.Chi2tot()/trk.Ndf()>chi2Cut)    // Skip bad chi2 tracks
      continue;

    double diffT = abs(trk.MeanTime()-t0), dt = trk.SigmaTime();
    if (!dt) continue;                      // Skip tracks w/o time measurement
    if (diffT>3*sqrt(dt0*dt0+dt*dt))        // Skip tracks offset in time
      continue;

    const PaTPar &hi = trk.vTPar()[0];
    if (hi(0)<ZTarget) continue;            // Skip beam tracks
    if (trk.vTPar().size()>1) {
      // First and last helices may have been swapped => try both.
      const PaTPar &hf = trk.vTPar()[1]; if (hf(0)<ZTarget) continue;
    }

    if (hi(0)>ZSM1)                         // Skip if fully downstream of SM1
      // - For simplicity's sake.
      // - Such tracks, when w/o momentum, would have to be rejected: their
      //  compatibility w/ the target cannot be checked.
      // - Else, those extending downstream of SM2 have a high P, which is then
      //  missing in the energy conservation budget of the pV, are likely not to
      //  pass the elasticity cut which the events we're here interested w/
      //  must have passed (since this routine is called from "fillExclusive").
      continue;
    double xp = hi(3), yp = hi(4);
    if (sqrt(xp*xp+yp*yp)<1e-3)             // Skip, paraxial, halo tracks
      // 10-3 is thought to be the divergence of the beam. I checked it on
      // a few MC events, at least...
      continue;

    //                                         Extrapolate to Zp
    trk.Extrap(Zp,H); double dX = H(1)-Xp, dY = H(2)-Yp;
    if (dX*dX+dY*dY>0.04) continue;         // Skip tracks far away from pVertex

    match = 1; break;
  }
  delete tPats;
  return match;
}
// **********************************************************************
// *************************  nInTimeNeutrals  **************************
// **********************************************************************
int nInTimeNeutrals(PaEvent &e, const PaVertex &pV)
{
  // Returns # of neutrals, i.e. ECAL non associated clusters, in time w/ pV. 

  double eCuts[2] = {0.6,1.2}; // Cf. http://wwwcompass.cern.ch/compass/software/analysis/transparencies/2013/am_130115/uhl_130116.pdf

  const PaParticle &beam = e.vParticle(pV.InParticle());
  const PaTrack &beamTrk = e.vTrack()[beam.iTrack()];
  double t0 = beamTrk.MeanTime(), dt0 = beamTrk.SigmaTime();

  vector<PaParticle> &particles = e.vParticle();

  int i, n; for (i=n = 0; i<e.NCaloClus(); i++) { // Loop over calo clusters
    const PaCaloClus &cl = e.vCaloClus(i);
    const string &name = cl.CalorimName(); if (name[0]!='E') continue;
    int ec; if (name=="EC01P1") ec = 0;                // ***** ENERGY LOWER CUT
    else    if (name=="EC02P1") ec = 1;
    else continue;
    if (cl.E()<eCuts[ec]) continue;

    //                                                            ***** TIME CUT
    if (!cl.HasTime()) continue; // No time? Skip, for what else can one do?
    double diffT = fabs(cl.Time()-t0), dt = cl.SigmaT();
    if (diffT>3*sqrt(dt*dt+dt*dt0)) continue;

    //                        ***** REJECT CLUSTERS ASSOCIATED TO CHARGED TRACKS
    int j, associated, nTot = cl.NParticles();
    for (j=associated = 0; j<nTot; j++) {
      if (particles[cl.iParticle(j)].iTrack()>=0) { associated = 1; break; }
    }
    if (associated) continue;

    n++;
  }
  return n;
}
#ifdef U3_OUTPUT_TREE
// **********************************************************************
// *************************  nInTimeNeutrals  **************************
// **********************************************************************
int nChargedKsInPV(PaEvent &e, const PaVertex &pV)
{
  int ip, nKs = 0; for (ip=nKs = 0; ip<pV.NOutParticles(); ip++) {
    int iET = e.vParticle(pV.iOutParticle(ip)).iTrack();
    PaTrack &trk = e.vTrack()[iET];
    if (!(PIDs[iET]&0x10)) continue;
    // Prior to, possibly re-scale momentum, determine exit angle for SM1
    // and incidence/exit angle for SM2. May need to extrapolate, if so
    // take advantage to set also acceptance by RICH.
    if (!aRs[iET]) /* If not yet done */ transport(e,trk,aRs[iET],tZones[iET]);
    if (aRs[iET]<0) continue;
    const PaTPar &hi = trk.vTPar()[0]; TVector3 v3 = hi.Mom3();
    double mom = v3.Mag(); if (mom<thrMargin*piThr || PRICHCut<mom) continue;
    double piLH = richLHs[0][iET], KLH = richLHs[1][iET];
    double pLH  = richLHs[2][iET], eLH = richLHs[3][iET];
    double hadronLH = pLH>piLH ? pLH : piLH;
    double elseLH = eLH>hadronLH ? eLH : hadronLH;
    if ((mom<KThr && elseLH<subKThrLHpiVeto) ||           // SubKThr else veto
	(KThr<mom && KLH>LHBckCut && KLH>hadronLH*LHCut)) // Actual ID
      nKs++;
  }
  return nKs;
}
#endif
// **********************************************************************
// ****************************** fillIphi *******************************
// **********************************************************************
void fillIphi(PaEvent &e, int ipV, int imuS, double dM_phi,
	      bool &uDSTSelection)
{
  const PaVertex &pV = e.vVertex(ipV); int nTrksPV = pV.NOutParticles();
  for (int ip1 = 0; ip1<nTrksPV; ip1++) {

    // *************************************************************
    // ***************   SEARCH FOR INCLUSIVE phi    ***************
    // - Let's disregard wrong sign pairs (at variance w/ the exlc. phi case):
    //  What we have in mind w/ these incl. phi's is to X-check the RICH perf,
    //  determined on excl. phi's, in the extended kinematical domain they allow
    //  to access. For this exercise, we do not foresee any acrobatic exercise
    //  that would benefit from a background evaluation based on wrong sign (
    //  contrary to what may happen w/ the excl. phi's where we want to also
    //  look at the misidentifications...)
    // *************************************************************

    int iEP1 = pV.iOutParticle(ip1); const PaParticle &pa1 = e.vParticle(iEP1);
    if (pa1.Q()<0) continue;                         // ***** 1ST PARTICLE IS >0

    //      ********** LOOP on >0 PARTICLE in PRIMARY VERTEX **********

    if (iEP1==imuS) continue;                               // ***** EXCLUDE mu'
    int iET1 = pa1.iTrack();
    PaTrack &trk1 = e.vTrack()[iET1]; int zones1 = tZones[iET1];

    //                                       ********** 1ST PARTICLE: REJECTIONS
    if (zones1==0x1) continue;                     // ***** EXCLUDE FRINGE FIELD
#ifdef DISCARD_MUONS
    if (trk1.XX0()>15) continue;            // ***** (UPON OPTION) EXCLUDE MUONS
#endif
    if (trk1.Chi2tot()/trk1.Ndf()>chi2Cut) continue;    // ***** CUT ON CHI2/NDF
    // Determine acceptance by RICH and exit angle for SM1 and incidence/exit
    // angle for SM2.
    if (aRs[iET1]==0) /* Not done? */transport(e,trk1,aRs[iET1],zones1);
    if (aRs[iET1]<-1) continue; // Skip track starting downstream of RICH
    const PaTPar &hv1 = pa1.ParInVtx(ipV); TVector3 v13 = hv1.Mom3();
    // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
    const PaTPar &hi1 = trk1.vTPar()[0]; double mom1;
    mom1 = fabs(1/hi1(5));
    double EK1 = sqrt(v13.Mag2()+M2_K); TLorentzVector lvK1(v13,EK1);
    for (int ip2 = 0; ip2<nTrksPV; ip2++) {
      int iEP2 = pV.iOutParticle(ip2); const PaParticle &pa2 = e.vParticle(iEP2);
      if (pa2.Q()>0) continue;                       // ***** 2ND PARTICLE IS <0

      //     ********** LOOP on <0 PARTICLES in PRIMARY VERTEX **********

      if (iEP2==imuS) continue;                              // ***** EXCLUDE mu'
      int iET2 = pa2.iTrack();
      PaTrack &trk2 = e.vTrack()[iET2]; int zones2 = tZones[iET2];

      //                                     ********** 2ND PARTICLE: REJECTIONS
      if (zones2==0x1) continue;                   // ***** EXCLUDE FRINGE FIELD
#ifdef DISCARD_MUONS
      if (trk2.XX0()>15) continue;          // ***** (UPON OPTION) EXCLUDE MUONS
#endif
      if (trk2.Chi2tot()/trk2.Ndf()>chi2Cut) continue;  // ***** CUT ON CHI2/NDF

      // Determine acceptance by RICH and exit angle for SM1 and incidence/exit
      // angle for SM2.
      if (aRs[iET2]==0) /* !done? */transport(e,trk2,aRs[iET2],zones2);
      if (aRs[iET2]<-1) continue; // Skip track starting downstream of RICH
      const PaTPar &hv2 = pa2.ParInVtx(ipV); TVector3 v23 = hv2.Mom3();
      // For P @ RICH, best would be a smoothing point in zone 0x2. Next best...
      const PaTPar &hi2 = trk2.vTPar()[0]; double mom2;
      mom2 = fabs(1/hi2(5));
      double EK2 = sqrt(v23.Mag2()+M2_K); TLorentzVector lvKK(v23,EK2);
      lvKK += lvK1; double m_KK = lvKK.M();
      if (fabs(m_KK-M_phi)>dM_phi) continue;
      double pT = lvK1.Perp(lvKK.Vect());
#ifdef Iphi_pTCUT
      if (pT<Iphi_pTCUT*.001) continue;
#endif
      double PL1 = sqrt(v13.Mag2()-pT*pT), PL2 = sqrt(v23.Mag2()-pT*pT);
      double alpha = (PL1-PL2)/(PL1+PL2);

      uDSTSelection = true;              // ********** uDST SELECTION: Incl. phi

      unsigned short id = 0; // 0x9: piTHr<P<PMax, 0x12: subKthr, 0x24: KID
      unsigned short idKpm = 0;// ID above or below K threshold
      unsigned short RICHOK = 0;
      int iET, pm; double mom; for (pm = 0, iET = iET1, mom = mom1; pm<2;
				    pm++) { 
	if (aRs[iET]>0 && thrMargin*piThr<mom && mom<PRICHCut) {
	  RICHOK |= 0x1<<pm; unsigned short idpm = 0;
	  if (PIDs[iET]&0x10) {
	    double piLH = richLHs[0][iET], KLH = richLHs[1][iET];
	    double pLH  = richLHs[2][iET], eLH = richLHs[3][iET];
	    double hadronLH = pLH>piLH ? pLH : piLH;
	    double elseLH = eLH>hadronLH ? eLH : hadronLH;
	    if (mom<thrMargin*KThr && // "thrMargin" used here and not infra...
		// ...so as to maximize efficiency and leave the problem of
		// purity to TTree projection time.
		elseLH<subKThrLHpiVeto) idpm |= 0x2;
	    if (KThr<mom && KLH>LHBckCut && KLH>hadronLH*LHCut) idpm |= 0x4;
	  }
	  else if (mom<thrMargin*KThr)
	    // Absence of RICH block or LH =0: cf. comment in PID for Lambda
	    idpm |= 0x2;
	  if ((idpm&0x4) ||     // Actual KID
	      (idpm&0x2))       // SubKThreshold piVeto
	    idKpm |= 1<<pm;
	  id |= idpm<<(3*pm);
	}
	iET = iET2; mom = mom2;
      }

      if (id) { // Let's require the KID of one of the hadrons
	if      ((id&0x24)==0x24) hm_Iphibb->Fill(m_KK); // Direct ID
	else if (id&0x24)         hm_Iphiab->Fill(m_KK); // Mixed case
	else                      hm_Iphiaa->Fill(m_KK); // SubKThr ID

      }
#ifdef U3_OUTPUT_TREE
      //                                                  ***** FILL OUTPUT TREE
      copyCSResonanceData(pV,m_KK,0,0,0,pT,alpha);
      CSResonanceData &r = fResonances.back();
      r.phiPat = 0x6|RICHOK<<4|idKpm<<6|0x400;
      const PaParticle *pa; for (pm = 0, pa = &pa1; pm<2; pm++) {
	int fHInd = getHadronIndex(pa);
	if (pm==0) r.h1 = fHInd; //adrons.size();
	else       r.h2 = fHInd; //adrons.size();
	if (fHInd==(int)fHadrons.size()) copyCSHadronData(e,pa,ipV);
	pa = &pa2;
      }
#endif
    } // End loop on pa2
  }  // End loop on pa1
}
// **********************************************************************
// ********************     getRICHMultiplicity     *********************
// **********************************************************************
void getRICHMultiplicity(PaEvent &e, const PaVertex &pV,
			 int &nTrksRIt, int &nTrksRIb)
{
  // Returns top/bottom (derived from incidence on RICH) track multiplicities
  // Requirements for the tracks:
  // - Reco'd dowstream of SM1 (checked as i) in zone 0x2, ii) upstream of RICH).
  // - Cut on chi2 and NDF: to filter out unreliable tracks.
  // - Timed (since tracks w/o time measurement, and thus w/ hits only in drift
  //  detectors) are typically either halo tracks (and hence off-time) or low P
  //  tracks (and hence, typically, off-RICH).
  // - Time compatible w/ event Time.

  const PaParticle &beam = e.vParticle(pV.InParticle());
  const PaTrack &beamTrk = e.vTrack()[beam.iTrack()];
  double t0 = beamTrk.MeanTime(), dt0 = beamTrk.SigmaTime();

  int nTrks = e.vTrack().size();
  int iET; for (iET = 0, nTrksRIt=nTrksRIb = 0; iET<nTrks; iET++) {
    int zones = tZones[iET]; if (!(zones&0x2)) continue;
    PaTrack &trk = e.vTrack()[iET]; int nDF = trk.Ndf();

    if (trk.Chi2tot()/nDF>chi2Cut || nDF<TRACK_NDF_MIN) continue;

    double diffT = abs(trk.MeanTime()-t0), dt = trk.SigmaTime();
    if (!dt) continue;                      // Skip tracks w/o time measurement
    if (diffT>3*sqrt(dt0*dt0+dt*dt))        // Skip tracks offset in time
      continue;

    const PaTPar &hi = trk.vTPar()[0]; float Z0 = hi(0);
    if (Z0>ZRICH) continue;
    float Yp; if (!hi(5)) {
      Yp = hi(4); float XR = hi(1)+hi(3)*(ZRICH-Z0), YR = hi(2)+Yp*(ZRICH-Z0);
      if (XR*XR+YR*YR<r2RICH || fabs(XR)>XRICH || fabs(YR)>YRICH) continue;
    }
    else {
      if (aRs[iET]==0) /* If not yet done */ transport(e,trk,aRs[iET],zones);
      if (aRs[iET]<0) continue;
      Yp = trk.vAux()[3];
    }
    if (Yp>0) nTrksRIt++;
    else      nTrksRIb++;
  }
}
void getFIMMAbscissae()
{
#ifdef U3_FIMP_VTCS
  for (int co = 0; co<3; co++) ZFI03s[co] = -1000;
  for (int co = 0; co<2; co++) ZMP01UVs[co] = -1000; 
#endif
  const PaSetup &setup = PaSetup::Ref();
  int idet, ndets = setup.NDetectors();
  for (idet = 0, ZFI02X=ZMM02X = -1000; idet<ndets; idet++) {
    const PaDetect &det = setup.Detector(idet);
    if      (det.Name().find("FI02X")==0) ZFI02X = det.Z();
    else if (det.Name().find("MM02X")==0) ZMM02X = det.Z();
    else if (det.Name().find("MP02X")==0) ZMM02X = det.Z();
#ifdef U3_FIMP_VTCS
    else if (det.Name().find("FI03X")==0) ZFI03s[0]= det.Z();
    else if (det.Name().find("FI03Y")==0) ZFI03s[1]= det.Z();
    else if (det.Name().find("FI03U")==0) ZFI03s[2]= det.Z();
    else if (det.Name().find("MP01U")==0) ZMP01UVs[0]= det.Z();
    else if (det.Name().find("MP01V")==0) ZMP01UVs[1]= det.Z();
#endif
  }
#ifdef U3_FIMP_VTCS
  for (int co = 0; co<3; co++) ZFI03s[co] += corrFI03;
  for (int co = 0; co<2; co++) ZMP01UVs[co] += corrMP01UV; 
#endif
}
void getTargetBinning(int &nZbins, double &ZMn, double &ZMx) {
  // Return [zMn,ZMx] rounded in 10cm.
  // Corresponding to zone large enough to encompass all kinds of target cuts.
  ZMn = ZFI02X<-999?-200:ZFI02X; ZMx = ZMM02X<-999?300:ZMM02X;
  double dZ = ZMx-ZMn; ZMn -= dZ*.05; ZMx += dZ*.05; dZ = ZMx-ZMn;
  // Have a binning in 
  ZMn = int(ZMn*10+.5); ZMx = int(ZMx*10+.5); nZbins = (ZMx-ZMn)*4;
  ZMn /= 10; ZMx /=10;
}
////////////////////          RCS HISTORY          ////////////////////
// $Log: UserEvent103.cc,v $
// Revision 1.13  2021/04/05 13:26:00  ybedfer
// Cosmetics...
//
// Revision 1.12  2021/03/31 17:33:58  ybedfer
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
// Revision 1.8  2018/04/26 09:11:21  ybedfer
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
//   - New subroutines: "checkExclusive,detachedHadrons,nInTimeNeutrals".
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
// Revision 1.5  2009/07/03 13:48:50  ybedfer
//  - Bug fix: negative pID by absence of photons (RICH block) is when
//   P < pThreshold!
//  - pT cut now "V0pTCut" (an array, for both Lambda and K0) conditioned by
//   "U3_V0_CUT".
//  - MC specifics:
//    - "MCLambda" class, for MC acceptance*efficiency.
//    - "U3_L_MC_ANALYSIS" = "U3_Lambda_ANALYSIS" except for PID.
//    - "U3_IDEAL_RICH".
//    - "histoLambdaLevel" taken into account in Lambda histogramming.
//    - Cut on trigger pattern (rejection of zero patterns) moved downstream (
//     of "MCInfo::SortLambdaOut").
//  - "U3_TARGET_y_CUTS" enabled by default.
//  - RICH acceptance flags "aRs" made a global variable.
//  - "fillExclusive" in an independent routine and simplified: no independent
//   histo of "dE" upon mu'ID type.
//  - "U3_DEBUG_W38_3" cancelled.
//
// Revision 1.4  2009/06/07 00:14:27  ybedfer
//  - "U3_V0_WO_PID"
//  - "U3_SELECTION_SCHEME": new option 2008, w/ RICH cuts = 2006's,
//   except for "piS_e_REJECTION" (why!?).
//  - Hadron
//    - "U_HADRON_SETUP": W/ impact on:
//      - "InitPIDs" (scattered mu, else..?).
//      - Beam selection.
//    - "U3_SELECTION_SCHEME==2008": target cut.
//    - Switch between LAS an SAS based upon "tZones", no longer "zLast" (
//     in fact not systematically, but then threshold is "zLast=1950").
//  - Primary vertex:
//    - No more condition on the sign of the incident particle.
//    - Meaning of "primary" modified.
//      Among others, no longer distinguishing between PaVertex and InitPID ID.
//     => When excluding muon later on: simplify.
//     => Change the meaning of histos: "hQ2p", "h_dE",...
//  - "U_SMU_FROM_CORAL": for debugging the impact of resting on coral for
//   the sc.mu ID. Not yet checked is the impact w/ coral ID is supplemented
//   by some recipe (supposed to make the unreliability of ID based on the
//   sole extrapolation).
//  - Momentum re-scaling:
//    - "SCALE_MOMENTA" cpp macro w/ options.
//    - "ScaleMomenta" class.
//    - Prior to re-scaling, determine incidence/exit angles using newly
//     introduced "transport" routine. (Note: "transport" modifies PaTrack.)
//  - V0 cuts:
//    - "U3_Lambda_V0_CUT" -> more appropriate name "U3_V0_CUT".
//    - "V0DsdD": [2] -> [3].
//    - "V0_STRICTER_SELECTION" macro cancelled.
//    - A priori impact on uDST selection.
//    - Chi2 cut < 10, instead of 5.
//  - piS_e_ALL_REJECTION: Veto all K0 histos if any of +/- is e. Note: yet
//   exclude veto from "K0LambdaphiXi_SELECTION".
//  - K0:
//    - Evaluation histos: Consider also incidence angle.
//    - Upon "U3_V0_HIGH_LEVEL": histos for evaluating mass shifts vs.
//     incidence
//  - Lambda:
//    - xF: New calculation. Supposed to be correct one, this time.
//    - Lambda selection w/ pT>25, instead of 50.
//    - "U3_Lambda_ANALYSIS":
//      - In 2008: #undef "U3_TARGET_y_CUTS" temporarily (why? guess it's
//       for maximising statistics).
//    - "Book/FillLambda":
//      - Arg. PaTrack's instead of double "zLast".
//      - Conditioned by "U3_V0_HIGH_LEVEL" instead of "U3_FILL_Lambda", which
//       is V0- and no longer only Lambda-specific.
//    - Upon "U3_V0_HIGH_LEVEL": histos for evaluating mass shifts vs.
//     incidence
//    - Lambda+K: New binning, centered on Omega.
//    - Cascades:
//      - Relaxed cut on collinearity.
//      - Stricter cut on CDA.
//      - "K0Lambdaphi" -> "K0LambdaphiXi", w/ special Xi uDST selection.
//  - Block of cpp macros: add a 5th sub-block for check/ensuring consistency
//   among macros.
//  - "tZones" made a static pointer.
//  - "zSM1,zSM2" Z absissae of magnets (to be retrieved from PaSetup)
//  - Bug fixes:
//    - Upon target/y/Q2 cuts, do not exit right away, but flag event for later
//     cleanly exiting.
//    - Subthreshold pID for K0+X.
//    - RICH perfs no longer condioned by (lately introduced) "isLambda".
//    - New RICH acceptance also in "U3_FILL_RICH_PERFS" conditioned code.
//  - GM07:09.
//  - Trigger type: a unique type for I and SI Middle.
//
// Revision 1.3  2009/06/06 19:36:40  ybedfer
//  - Cleaner set of cpp macros:
//    - Reminder: all macros here after mentioned are "object-like".
//    - Clearer rationale for the overall macro scheme, w/ a clear separation
//     between high-level and low-level macros, w/ high-level maros being sets
//     of low-level ones corresponding to:
//      - either a particular configuration of processing (e.g. final analysis
//       or uDST selection),
//      - or a particular (2002,..,2008) context of data taking in the case of
//       "U3_SELECTION_SCHEME".
//    - The macro definition block is organised in 4 sub-blocks:
//      I) Block enabling (w/ a #define statement) one, and only one, of the
//        high-level macros of the ''processing'' kind, but also listing the
//        rest of them, commented out.
//     II) Block of "U3_SELECTION_SCHEME": w/ an assignment (to 2002,..,2008).
//    III) Block w/ the list of all low-level definitions, some of which being
//        uncommented (the objective is to have them all uncommented).
//     IV) Block specifying the contents of all the high-level definitions, of
//        the ''processing'' kind, i.e. sequences of "#undef#" or "#define"
//        statements, where "#define"'s overwrite default assignments.
//    - Take advantage of this re-shuffling to opt for more adequate macro
//      names: U3_microDST_SELECTION -> U3_uDST_PRESELECTION
//  - New high level macros:
//    - U3_Lambda_uDST_PRESEL = uDST selection having Lambda analysis in
//     mind, although it does cater for K0 based physics channels. The older
//     "U3_uDST_PRESELECTION" is kept, as a backup.
//    - U3_Lambda_ANALYSIS.
//    - U3_ERASE_K0s: to be used in conjunction w/ the D0(->Kpi)-dedicated
//     "UserEvent5": it removes pions from K0 decay from the pVertex.
//  - "Lamda_p_ID": Low-level cpp macro conditioning a new algorithm of Lambda
//   selection, resorting to RICH-based p(roton) ID.
//  - A set of macros allowing to by-pass part of the processing steps:
//   "U3_FILL_RICH_PERFS", "U3_FILL_Lambda", "U3_GENERALPURPOSE_HISTOS", taking
//   care (in part) not to declare dedicated variables.
//  - RICH cuts, conditioned by several cpp macros:
//    - U3_SELECTION_SCHEME
//    - U3_Lambda_LOOSE (typically "LHCut = 0.95, LHBckCut = 0.95",
//     "subpThr[pi]LHVeto = 1.80".
//    - Acceptance from 2008/8/6 Giulia Pesaro <gpesaro@cern.ch>.
//  - Event cuts:
//    - "U3_Q2_CUT".
//    - Target cut depends upon "U3_SELECTION_SCHEME".
//  - V0 cuts: playing a bit dependence upon high-level macro, but no fully
//   satisfying solution.
//  - "Book|FillLambda" w/ extra arguments for PID and some extra histos.
//  - "xF" bins.
//  - RICH, MCInfo: Instantiate it only once.
//  - More histos: Lambda+pi, Lambda+K, Q2 and else for Lambda.
//  - Definition of histos: collinearity range.
//  - "K0Lambdaphi_SELECTION": moved from (1,2) to (1,1).
//  - K0 selection:
//    - "DEBUG_K0".
//    - "U3_DEBUG_e_REJECTION"
//    - K0+X PID: introduce "subpThrLHpiVeto" // to "subKThrLHpiVeto".
//  - Lambda selection:
//    - Macro "Lambda_p_ID".
//    - Selection evaluation: moved behing cut on pT and vChi2.
//    - ...
//  - RICH perfs: discard "LHDIFF_SELECTION".
//
// Revision 1.2  2008/04/17 20:58:05  ybedfer
//  Initialisation of "tZones, aRs, PIDs" (w/ "memset"): moved to after all allocations
// have been done.
//
// Revision 1.2  2004/05/20 22:56:45  ybedfer
//  Reshuffling:
//   - Subroutines for booking and filling RICH perfs.
//   - Header file for the definition of histo's
//////////////////////////////////////////////////////////////////////

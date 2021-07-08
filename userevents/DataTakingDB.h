// $Id: DataTakingDB.h,v 1.12 2020/10/14 14:54:47 ybedfer Exp $

// Class DataTakingDB
// - defining some cuts, correction factors and other parameters,
// - based on a built-in data base.
// - to be used to set those quantities (cuts, etc...) given run#

#ifndef DataTakingDB_h
#define DataTakingDB_h

#include "PaAlgo.h"
#include "PaEvent.h"
#include "DISTarget.h"

class DataTakingDB
{
public:
  DataTakingDB();

  void Init(PaEvent &e,
	    // =1: Get best pV from PHAST
	    // =2: Get it from CORAL, if latter is guessed to be up-to-date
	    // =3: Get it fomr CORAL in any case
	    int bestPVOption,
	    // =1: Loose chi2 cut (typically for microDST preselection)
	    // =2: Chi2 cut set depending upon period
	    int chi2CutOption,
	    // =0: No RICH settings requested
	    // =1: Default (as a f(year))
	    // =2: Special loose setting for Lambda
	    int RICHOption,
	    // =0: No cut
	    // =1: Clone "PaAlgo::GetTargetLocation"
	    // =2: Home made
	    // =3: DISTarget if relevant, else as for =1
	    // =4: Between bT and spectro
	    int targetOption); // !=0: Include FIMP (i.e. FI03/MP01UV vertices))

  //! Returns pointer to DataTakingDB
  static DataTakingDB *Ptr() { return address; }

  int Year;
  const string *Period;  // E.g. "P1C" or "W42"
  int SubPeriod;// Some periods are split, cf. "periodDB" in "./DataTakingDB.cc"
  const string *PhysicsType; 
  int StrtRun, StopRun;
  double MeanIndex; // Mean refractive index
  int ThCIndex;     // Index of estimator for theta Cerenkov
  double PBeam, DPBMS;  // DPBMS =0 means no BMS, i.e. not a muon beam
  double DqP2BMSCut;   // Cut on cov(1/P,1/P) to reject beams w/o BMS
  double yLowCut, yUpCut;
  int  BeamRecMethod;// 0: old, 1: Pawel's, -1 = old, but backpropag unreliable
  bool BadBackProp(const PaParticle &muI) {
    if (BeamRecMethod==0) return muI.Chi2CutFlag();
    if (BeamRecMethod==-1) return false;
    double lh = muI.BackPropLH();
    return lh<0.005 && 1<lh;
  }

  // Settings
  // - These are determined on view on the sole data taking.
  // - Reconstruction (coral version, production slot#...) is not taken into
  //  account, while it certainly can have an impact on such setting as
  //   - "Chi2Cut" (which depends on e.g. the quality of alignment available to
  //    the production),
  //   - "GetBestPVfromCORAL" (i.e. use "iBestCoralPrimaryVertex" instead of
  //    "iBestPrimaryVertex", w/ the meaning of the former varying w/ coral
  //    version.
  // => It's implicitly assumed the processed DST data come from the latest
  //   production, or rather from the latest known to the author of this class.
  double Chi2Cut;
  bool GetBestPVfromCORAL;

  // ***** RICH
  // Basics
  double PRICHCut;	// Max. P for pi vs. K separation (case of p is different)
  double LHCut;		// Cut on p1LH/p2LH, for particles p1,p2 = pi, K
  double LHBckCut;	// Cut on p1LH/bckLH for particles p1 = pi, K. Also used for for eVeto: eLH1>LHeVeto*LHBckCut

  // pID (e.g. for Lambda physics and K0+p search)
  double pIDPmax;	// In GeV
  // Thresholds (e.g. pThr) can vary (just as refractive index varies):
  // therefore lower bounds are expressed as a factor of threshold (upper could
  // be so just as well, only that then it does not matter from a statistics
  // point of view.
  double pIDPmin;	// Times pThr
  double SubpThrPmin;	// Time piThr. To be used to veto pions in Lambda physics
  // piVeto in pID.
  // - Can be demending in term of cut piLH/bckLH<cut: in Lambda physics e.g.,
  //  we want to remove the bulk of the K0's mimicking Lambda's, the rest being
  //  removed by a mass cut around the (very narrow: 1-2 MeV) Lambda peak.
  // - In a hadron multiplicity analysis, on the contrary, where we would want
  //  to separate p from K or pi over the P range P>KThr, we need be stricter
  //  => another cut is used.
  double SubpThrLHpiVeto;	// For piLH/bckLH
  double SubpThrLHVeto;		// For e,KLH/bckLH

  // ***** TARGET
  int TargetOption; // =1: Clone "PaAlgo::GetTargetLocation", =2: Home made, =3: DISTarget id relevant, else as for =1, =4: Between bT and spectro
  int TargetFIMP;   // !=0: Include FIMP
  int TargetType;
  // - Following "PaAlgo::GetTargetLocation" notation, but for (x,y,z)->(X,Y,Z)
  // - Some definition are made more universal, applying also when there are
  //  fewer than 3 cells
  // ZU_1  z position of the most upstream part of the upstream cell
  //    -> z position of the most upstream part of the target
  // ZD_2  z position of the most downstream part of the downstream cell
  //    -> z position of the most downstream part of the target
  // - For the rest, let's recall that
  // XU    shift in x of the most upstream part of the target
  // YU    shift in y of the most upstream part of the target
  // etc...
  // XD    shift in x of the most downstream part of the target
  // YD    shift in y of the most downstream part of the target
  // RTCut  the radial cut
  // YTCut  Y cut for top part of the cell
  // If there are fewer than 3 cells. 
  double XU, YU, ZU_1, ZU_2, ZC_1, ZC_2, XD, YD, ZD_1, ZD_2;
  double ZTMn, ZTMx, R2TMx; // Target zone cuts
  double RTCut, R2TCut, YTCut;
  // "WinTarget" depends upon "TargetOption/TargetFIMP"
  unsigned int WinTarget(double X, double Y, double Z);
  // "WinTargetZone" depends upon "TargetFIMP"
  bool WinTargetZone(double X, double Y, double Z);
  bool CellsCrossed(const PaTPar &par);
  DISTarget *fDISTarget;

  // ***** SCATTERED MU
  bool CheckYokeSM2;


  // ***** BAD SPILLS
  set<unsigned> BadSpills, BadSpillsCalo, BadSpillsRICH, BadSpillsRPD;
  void GetBadSpills();
  void GetBadSpills(set<unsigned> &badSpills,
		    set<unsigned> &badSpillsCalo, set<unsigned> &badSpillsRICH);
  unsigned int BadSpillPat(int run, int spill);

private:
  static DataTakingDB *address;        //! (-) address of this object
};
#endif

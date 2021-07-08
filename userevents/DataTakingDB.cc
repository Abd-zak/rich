// $Id: DataTakingDB.cc,v 1.12 2020/10/14 14:54:47 ybedfer Exp $

// DataTakingDB contructor
// - defining some cuts, correction factors and other parameters,
// - based on a built-in data base.
// - to be used to set those quantities (cuts, etc...) given run#


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "TSystem.h"
#include "TFile.h"

#include "Phast.h"
#include "PaAlgo.h"

#include "DataTakingDB.h"

// Interface
bool periodDB(int run,
	      int &year, const string *&period, int &subPeriod,
	      int &strtRun, int &stopRun,
	      double &meanIndex,
	      const string *&physicsType);
void settingDB(int year, const string &period, const string &physicsType,
	       int &thCIndex,
	       double &pBeam, double &dpBMS, int &beamRecMethod,
	       double &dqP2BMSCut,
	       double &yLowCut, double &yUpCut, double &chi2Cut);
bool XCheckYearSetup(int year);

DataTakingDB *DataTakingDB::address = 0; // Init static pointer

DataTakingDB::DataTakingDB()
{
  if (address) {
    printf("\n** DataTakingDB:\a DataTakingDB already instantiated\n\n");
    abort();
  }
  address = this;

  Year = 0; Period = 0; SubPeriod = 0;
}

void DataTakingDB::Init(PaEvent &e,
			int bestPVOption, // =1: CORAL if newer coral, =2: CORAL in any case, else: PHAST
			int chi2CutOption, int RICHOption,
			// "targetOption": =0: NO CUT, =1: SATNDARD, =2: HOME MADE, =3: DISTarget (if relevant, else as for =1), =4: LOOSE (between bT and spectro)
			int targetOption)
{
  const PaSetup &setup = PaSetup::Ref();

  // Backup previous values
  int prvYear = Year;
  const string *prvPeriod = Period ? new string(*Period) : 0;
  int prvSubPeriod = SubPeriod;

  // Year and period from mDST
  Year = e.Year();                     if (Year==1900) Year = 0;
  Period = new string(e.PeriodName()); if (*Period=="XXXX") Period = 0;

  // "PeriodDB" OK else abort => "Year/Period" receive a setting in any case!
  int year = 0; const string *period;
  if (!periodDB(e.RunNum(),year,period,SubPeriod,StrtRun,StopRun,
		MeanIndex,
		PhysicsType)) abort();

  // ***** X-CHECK "PeriodDB" against PHAST *****
  if (Year && year!=Year) {
    printf("\n** DataTakingDB: DataTaking year retrieved by PHAST from DST file differ from that derived from DST's run#(=%d) based on built-in DB:\n   PHAST(%d) != DB(%d)\n\n",e.RunNum(),Year,year); abort();
  }
  if (Period && *Period==*period) {
    printf("\n** DataTakingDB: DataTaking period retrieved by PHAST from DST file differ from that derived from DST's run#(=%d) based on built-in DB:\n   PHAST(%s) != DB(%s)\n\n",e.RunNum(),Period->c_str(),period->c_str()); abort();
  }
  Year = year; Period = period;

  // ***** ADDITIONAL X-CHECK
  if (!XCheckYearSetup(Year)) abort();

  // ***** NO CHANGE? RETURN
  if (Year==prvYear && *prvPeriod==*Period && SubPeriod==prvSubPeriod) {
    printf("\n * DataTakingDB: Run %d %d/%s => No change\n",
	   e.RunNum(),Year,Period->c_str());
    return;
  }

  // ***** DERIVE SETTINGS from Year and PhysicsType
  settingDB(Year,*Period,*PhysicsType,
	    ThCIndex,PBeam,DPBMS,BeamRecMethod,DqP2BMSCut,yLowCut,yUpCut,Chi2Cut);


  // ***** SUMMARY of PARAMETERS and OPTION *****
  printf("\n\n * DataTakingDB:\n");

  printf("  Year = %d",Year);
  printf(" => %s Physics, Chi2 Cut < %.1f\n",PhysicsType->c_str(),Chi2Cut);
  if (chi2CutOption==1) {
    double looseCut = 10; 
    if (looseCut-Chi2Cut>.1)
      printf("   Chi2 cut from DB (<%.1f) downgraded to < %.1f  (because option \"Loose chi2 cut\" requested)\n",Chi2Cut,looseCut);
    else
      printf("   Chi2 cut from DB (<%.1f) does not contradict requested \"Loose chi2 cut\" option\n",Chi2Cut);
    Chi2Cut = looseCut;
  }

  // ***** BEST PRIMARY VERTEX STRATEGY
  // - CORAL nice algorithm is only available in recent versions.
  // - Let's guess if current DST data could have been produced w/ it.
  // - We assume in any case that processed data are from the latest production.
  // - And as of 2012/10: the production status seems to be: [2002,2006] and
  //  2007 longitudinal are too old.
  bool newerCoral = Year<=2006 ? false : true;
  if (Year==2007) { // Is it longitudinal? Look at target field
    PaMagInfo *m = setup.PtrMagField()->getMagInfo();
    int nms = setup.PtrMagField()->getNumOfMags();
    if (nms!=3) {
      printf("\n** DataTakingDB:\a # of magnets in PaSetup = %d, while 3 are expected in 2007\n\n",nms); abort();
    }
    if (m[0].flg1==4) /* transverse */ newerCoral = true;
  }
  if (newerCoral)
    printf("\n   DST data were (most probably) produced w/ a newer CORAL => best pV from CORAL is an interesting alternative...\n");
  else
    printf("\n   DST data were (most probably) produced w/ an older CORAL => best pV from CORAL is NOT an interesting alternative...\n");
  // "bestPVOption": =1: CORAL if newer coral, =2: CORAL in any case, else: PHAST
  GetBestPVfromCORAL = bestPVOption==2 || bestPVOption==1 && newerCoral;
  if (GetBestPVfromCORAL) {
    if (newerCoral)
      printf("  ... => Requested \"Best pV from CORAL\" option granted.\n");
    else if (bestPVOption==1)
      printf("  ... => Requested \"Best pV from CORAL\" option turned down.\n");
    else
      printf("  ... Nevertheless \"Best pV from CORAL\" was requested imperatively!\n");
  }
  else if (newerCoral)
    printf(" ... Option \"Best pV from CORAL\" was NOT requested, though.\n");
  else
    printf(" ... Option \"Best pV from CORAL\" was NOT requested, anyway.\n");

  // ***** BEAM

  // ***** SCATTERED MUON
  if (Year<=2006) CheckYokeSM2 = true;
  else            CheckYokeSM2 = false;
  if (!GetBestPVfromCORAL)
    printf("\n   Option \"Best pV from PHAST\": CheckYokeSM2 %srequested.\n",
	   CheckYokeSM2?"":"NOT ");

  // ***** RICH
  if (RICHOption) {
    // Settings independent of "RICHOption"
    PRICHCut = Year>=2006 ? 60 : 45;
    // Settings dependent on "RICHOption"
    if      (RICHOption==1) { // Default
      if (Year>=2006) {
	LHCut = 1.01; LHBckCut = 1.20;
      }
      else {
	LHCut = 1.01; LHBckCut = 1.20;
      }
      // pID (for Lambda physics and K0+p search)
      pIDPmax     = 80   /* GeV */;
      pIDPmin     = 1.05 /* times pThr*/;
      SubpThrPmin = 1.1  /* times piThr */;
      SubpThrLHpiVeto = 1.40; SubpThrLHVeto = 1.50;
    }
    else if (RICHOption==2) { // U3_Lambda_LOOSE
      if (Year>=2006) {
	LHCut = 0.95; LHBckCut = 0.95;
      }
      else {
	printf("\n Special RICHOption = 2 requested, while year of data taking != 2006\n\n");
	abort();
      }
      pIDPmax     = 100   /* GeV */;
      pIDPmin     = 1.05  /* times pThr*/;
      SubpThrPmin = 1.05  /* times piThr */;
      SubpThrLHpiVeto = 1.80; SubpThrLHVeto = 1.80;
    }
    //else...
//static double LHCut = 1.01, LHBckCut = 1.02, subKThrLHpiVeto = 0.92;  // 04W39
//static double LHCut = 1, LHBckCut = 1.2, subKThrLHpiVeto = 0.92;  // 04W39
//static double LHCut = 1.01, LHBckCut = 1.1, subKThrLHpiVeto = 0.98;
//static double LHCut = 1.01, LHBckCut = 1.02, subKThrLHpiVeto = 1.01;
//static double LHCut = 1.01, LHBckCut = 1.05, subKThrLHpiVeto = 1.00;
//static double LHCut = 1.0001, LHBckCut = 1.01, subKThrLHpiVeto = 1.01;  // last

//static double LHCut = 1.0001, LHBckCut = 1, subKThrLHpiVeto = 1;  // Std
//static double LHCut = 1.0001, LHBckCut = 1.01, subKThrLHpiVeto = 1.2;  // Stf
//static double LHCut = 1.0001, LHBckCut = 1.05, subKThrLHpiVeto = 1.01; // Stf
//static double LHCut = 1.0001, LHBckCut = 1.20, subKThrLHpiVeto = 1.01; // Stf
  }

  //                         ********** TARGET **********

  // "targetOption":
  // =0: No cut
  // =1: Clone "PaAlgo::GetTargetLocation",
  // =2: Home made,
  // =3: DISTarget if relevant, else as for =1,
  // =4: Between bT and spectro, or (MP01U+MP01X)/2 in 2016/17
  TargetOption = targetOption;
  if (TargetOption) {
    if (TargetOption==3) {
      //         ***** DISTarget: RELEVANT? IF NOT FALL BACK ON =1 *****
      printf("\n** DataTakingDB: \"DISTarget\" option requested");
      if (!e.IsMC()) {
	TargetOption = 1;
	printf(", while !MC\n => Falling back to \"PaAlgo\" option.\n");
      }
      else {
	TFile *in_file = Phast::Ref().in_file;
	// => let's look at the file name, search for a ".r<version#>." field.
	int svnRevision = 0;
	const string spath (in_file->GetName()); int nidx = spath.rfind("/");
	string sname = spath.substr(nidx+1,string::npos);
	const char *fname = sname.c_str();
	printf(", while analyzing run #%d \"%s\"\n",e.RunNum(),fname);
	int ridx = sname.find(".r"); if (ridx) {
	  const char *rfield = fname+(ridx+2)*sizeof(char); char *end, **endptr = &end;
	  int i = strtol(rfield,endptr,10);
	  if (*endptr!=rfield && **endptr=='.') {
	    printf(" => SVN version# = %d\n",i); svnRevision = i;
	    fDISTarget = new DISTarget(Year,svnRevision);
	  }
	}
	if (!svnRevision) {
	  TargetOption = 1;
	  printf(" => No retrieving SVN version# => Falling back to \"PaAlgo\" option.\n");
	}
      }
    }

    //            ***** SWITCH DEPENDING ON "TargetOption" *****
    if (TargetOption==3) { // I.e. TargetOption #3 requested and granted
      // Data members are not used in "DataTakingDB::WinTarget/CellsCrossed"
      // Nevertheless  
      fDISTarget->getCuts(RTCut,YTCut);
      fDISTarget->getPos(XU,YU,ZU_1,ZU_2,
			 ZC_1,ZC_2,    
			 XD,YD,ZD_1,ZD_2);
    }
    else if (TargetOption==4) { // Loose target: delimited by last bT and
      // first spectro detector, enlarged by +/-1 cm
      RTCut=YTCut = 0;
    }
    else {
      if      (2008<=Year && Year<=2009 ||
	       Year==2012 && *PhysicsType!="DVCS")
	TargetType = 2008;
      else if (Year==2012 || Year==2016 || Year==2017)
	TargetType = 2012;
      else if (Year>=2007)
	TargetType = 2007;
      else if (Year==2006)
	TargetType = 2006;
      else if (Year<=2004)
	TargetType = 2002;
      else
	TargetType = 0;

      if (TargetOption==1) {
	if (!TargetType) {
	  printf("\n** DataTakingDB:\a Target cut requested for Year = %d: no method available\n\n",Year);
	  abort();
	}
	if (TargetType==2008) {
	  printf("\n** DataTakingDB:\a Target cut alla \"PaAlgo\" requested for Year = %d: no method available\n\n",Year);
	  abort();
	}

	if (TargetType==2012) {
	  double dummy; vector<double> vDummy;
	  PaAlgo::GetDVCSTargetLocationCenter
	    (e.RunNum(),dummy,dummy,dummy,dummy,0,RTCut,dummy,YTCut,vDummy);
	  // - "ZU_1" and "ZD_2" were originally taken from:
	  // "/afs/cern.ch/compass/dvcs/Production/Analysis/Students/avidon/phast/user/FluxUtils.h"
	  //  =>
 	  //  if      (Year==2012) { ZU_1 = -311.19; ZD_2 = -71.19; }
	  //  else if (Year==2016) { ZU_1 = -318.5;  ZD_2 = -78.5;  }
	  // - Since this file does not exist any longer, and I don't know where
	  // to find the up-to-date values
	  //  =>
	  // - 2012: Recycle "avidon" 
	  // - 2016: Assign Nicolas Z range, see "nicolas/analySIDIS_split.cc"
	  //  NOTA BENE: Would still have to check whether values correspond
	  // to the flux calculation.
	  // - 2017: Don't know yet => Large values, as a reminder that
	  //  something needs to be tuned.
	  if (Year==2012)      { ZU_1 = -311.19; ZD_2 = -71.19; }
	  else if (Year==2016) { ZU_1 = -325;    ZD_2 = -71; }
	  else { // As of 2020/07, it's 2017
	    /* */              { ZU_1 = -500;    ZD_2 = -50; }
	  }
	}
	else {
	  bool ok;
	  if (Year>=2006)
	    ok = 
	      PaAlgo::GetTargetLocation(e.RunNum(),XU,YU,ZU_1,ZU_2,ZC_1,ZC_2,
					XD,YD,ZD_1,ZD_2,RTCut,YTCut);
	  else
	    ok =
	      PaAlgo::GetTargetLocation(e.RunNum(),XU,YU,ZU_1,ZU_2,
					XD,YD,ZD_1,ZD_2,RTCut,YTCut);
	  if (!ok) {
	    printf("\n** DataTakingDB:\a Target cut alla \"PaAlgo\" requested for run %d: unsuccessful!\n\n",e.RunNum());
	    abort();
	  }
	}
      }
      else if (TargetOption==2) { // =2: Home made
	if      (TargetType==2012) {
	  // http://wwwcompass.cern.ch/compass/results/2014/october_dvcs_camera/note_DVCS_2012.pdf
	  ZU_1 = -311.2; ZD_2 = -71.2; XU=XD = 0; YU = -.45; YD = -.8;
	  RTCut = 1.9;
	}
	else if (TargetType==2008) {
	// From "$COMPASS_FILES/2009/detectors.*.dat"
	/*
	  targ    1  TMAT     15      5    0.000   3.500  40.000   -48.5200     0.0000     0.0000
	*/
	  ZU_1 = -68.52; ZD_2 = -48.52; XU=XD=YU=YD = 0;
	  RTCut = 1.75;
	}
	else {
	  printf("\n** DataTakingDB:\a TargetOption =2 requested for TargetType =%d: no method available\n\n",TargetType); abort();
	}
	YTCut = RTCut;
	R2TCut = RTCut*RTCut;
      }
    }
    R2TCut = RTCut*RTCut;
  }

  //                 ***** TARGET ZONE *****
  double ZTarget = setup.TargetCenterZ();
  R2TMx = 10000;  // Large initialisation value.
  double ZMP01U, ZMP01X;
  int imap; for (imap = 0, ZTMn = ZTarget+1, ZTMx = ZTarget-1,
		   ZMP01U=ZMP01X = 0;
		 imap<(int)setup.NDetectors(); imap++) {
    const PaDetect &d = setup.Detector(imap); double Z = d.Z();
    if      (Z<ZTarget)    ZTMn = Z;
    else if (ZTMx<ZTarget) ZTMx = Z;
    if      (d.Name().find("MP01U")==0) ZMP01U = d.Z();
    else if (d.Name().find("MP01X")==0) ZMP01X = d.Z();
  }
  if (ZTMn>ZTarget) {
    printf("\n** DataTakingDB:\a Error assigning target zone lower bound!\n\n"); abort();
  }
  if (ZTMx<ZTarget) {
    printf("\n** DataTakingDB:\a Error assigning target zone upper bound!\n\n"); abort();
  }
  // Add/subtract a margin, not to mix in vertices from detectors' material.
  // Make it large enough to avoid cases where tracking confuses hits from
  // incident and outgoing track and subsequent vertexing biases momentum.
  // So that gain obtained by cuting on target-zone instead of target proper is
  // not at the expense of event quality.
  ZTMn += 2; ZTMx -= 2;
  bool targetFIMP = // !=0: Include FIMP (i.e. FI03/MP01UV vertices)
    2016<=Year && Year<=2017;
  if (targetFIMP) { // Upper bound of target zone half-way between MPU and MPX
    if (!ZMP01U || !ZMP01X) {
      printf("\n** DataTakingDB:\a Error assigning FIMP upper bound!\n\n"); abort();
    }
    ZTMx = (ZMP01U+ZMP01X)/2;
    // Radius cut: Larger than the radius of target proper. Should be used to bin the Y vs. X vertex histo
    double RTMx = 2.5; R2TMx = RTMx*RTMx; 
  }

  printf("\n");
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   WinTarget   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int DataTakingDB::WinTarget(double X, double Y, double Z)
{
  unsigned int winTarget;  // 0x1: w/in Z cut, 0x2: w/in R cut
  // "TargetOption"
  // =0: NO CUT,
  // =1: STANDARD: Clone "PaAlgo::GetTargetLocation",
  // =2: HOME MADE,
  // =3: DISTarget if relevant, else as for =1,
  // =4: LOOSE: Between bT and spectro
  if     (TargetOption==3) {                        // ***** DISTarget
    winTarget = fDISTarget->Zcut(Z) ? 0x1 : 0; 
    if (fDISTarget->Rcut(X,Y,Z)) winTarget |= 0x2;
    return winTarget;
  }
  else if (TargetOption==2) {                       // ***** HOME MADE
    if (TargetType==2012) {                   // 2012
      winTarget = ZU_1<Z && Z<ZD_2 ? 0x1 : 0;
      double Y_Z = YU+(Z-ZU_1)*(YD-YU)/(ZD_2-ZU_1);
      if (X*X+(Y-Y_Z)*(Y-Y_Z)<R2TCut) winTarget |= 0x2;
    }
    else {
      if (TargetType==2008) {
	winTarget = ZU_1<Z && Z<ZD_2 ? 0x1 : 0;
	if (X*X+Y*Y<R2TCut) winTarget |= 0x2;
      }
    }
    return winTarget;
  }
  else if (TargetOption==4)                         // ***** LOOSE
    return WinTargetZone(X,Y,Z);
  else if (TargetOption==1) {                       // ***** STANDARD...
    if (TargetType==2012) {                           // ...2012,17-17: PaAlgo
      winTarget = ZU_1<Z && Z<ZD_2 ? 0x1 : 0;
      int run = PaEvent::Ptr()->RunNum();
      if (PaAlgo::InTarget(X,Y,Z,'O',run,RTCut,YTCut,ZU_1,ZD_2,RTCut))
	winTarget |= 0x2;
      return winTarget;
    }

    // Extracted from "PaAlgo::InTarget" mutatis mutandis.
  
    int udc; for (udc=winTarget = 0; udc<(TargetType==2002?2:3); udc++) {
      if      (udc==0) { //if(        Cell == 'U' ) {
	if( Z > ZU_1 && Z < ZU_2 ) { winTarget = 0x1; break; }
      }
      else if (udc==2) { //} else if( Cell == 'D' ) {
	if( Z > ZD_1 && Z < ZD_2 ) { winTarget = 0x1; break; }
      }
      else {             //} else if( Cell == 'C' ) {
	if( Z > ZC_1 && Z < ZC_2 ) { winTarget = 0x1; break; }
      }
      /*
	} else {
	cout<<"inTarget PROBLEM: no info for cell "<<Cell<<endl;
	return false;
      */
    }
    double XC = (XD-XU) * (ZU_1-Z) / (ZU_1-ZD_2) + XU;
    double YC = (YD-YU) * (ZU_1-Z) / (ZU_1-ZD_2) + YU;
    double R2 = (X-XC)*(X-XC) + (Y-YC)*(Y-YC);
    if(        R2 > R2TCut    ) return winTarget;
    if(      Y-YC > YTCut     ) return winTarget;
    winTarget |= 0x2;
    return winTarget;
  }
  else // I.e. "TargetOption==0"
    return true;
}
bool DataTakingDB::WinTargetZone(double X, double Y, double Z)
{
  return ZTMn<Z && Z<ZTMx && X*X+Y*Y<R2TMx;
  
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  CellsCrossed ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool DataTakingDB::CellsCrossed(const PaTPar &par)
{
  if (TargetOption==1 && TargetType==2012) {
    int run = PaEvent::Ptr()->RunNum();
    return PaAlgo::CrossCells(par,run,RTCut,YTCut,ZU_1,ZD_2,RTCut);
  }
  if (TargetOption==3) return fDISTarget->cellsCrossed(par);

  // Derived from "PaAlgo::CrossCells"
  PaTPar parE;


  par.Extrapolate(ZU_1,parE,0);
  double XL = parE(1), YL = parE(2), R2L = (XL-XU)*(XL-XU)+(YL-YU)*(YL-YU);
  if (R2L>R2TCut) return false;
  if (YL-YU>YTCut) return false;

  par.Extrapolate(ZD_2,parE,0);
  double XR = parE(1), YR = parE(2), R2R = (XR-XD)*(XR-XD)+(YR-YD)*(YR-YD);
  if (R2R>R2TCut) return false;
  if (YR-YD>YTCut) return false;
  return true;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    periodDB   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool periodDB(int run,
	      int &year, const string *&period, int &subPeriod,
	      int &strtRun, int &stopRun,
	      double &meanIndex,
	      const string *&physicsType)
{
  // ***** DEFAULT VALUES
  year = 0;
  int iPeriod = -1; period = 0; subPeriod = 0;
  strtRun = run; stopRun = run+100;
  meanIndex = 0;
  physicsType = 0;

  // ******************** DETERMINE PERIOD ********************
  const int nY2s = 6;
  const int Y2Runs[nY2s][2] = { {20019,20351},   // P1C
                                {20437,21134},   // P2A
                                {22018,22341},   // P2D
                                {22373,22698},   // P2E
                                {22754,22972},   // P2F
                                {23016,23449} }; // P2G

  const int nY3s = 10;
  const int Y3Runs[nY3s][2] = { {27599,28178},   // P1A 
                                {28232,28433},   // P1B
                                {28568,28879},   // P1C
                                {28881,29366},   // P1D
                                {29963,30379},   // P1E 
                                {30448,30717},   // P1F 
           {0,0},{0,0},         {31581,31855},   // P1I 
                                {31930,32236} }; // P1J 

  const int nY4s = 27;
  const int Y4Runs[nY4s][2] = { {    0,    0},
				{34930,35228},   // W21
				{35230,35605},   // W22
				{35606,35841},   // W23			      
           {0,0},{0,0},         {36672,36949},   // W26
				{37007,37145},   // W27
				{37183,37521},   // W28
				{37522,37767},   // W29
				{37846,38061},   // W30
				{38065,38520},   // W31
				{38521,38920},   // W32
				{38921,39169},   // W33
				{39173,39545},   // W34
				{39546,39818},   // W35
				{39843,39988},   // W36
				{39999,40450},   // W37
				{40453,40951},   // W38
				{40952,41165},   // W39
				{41220,41397},   // W40
           {0,0},{0,0},{0,0},   {42282,42680},   // W44
				{42704,42945},   // W45
				{42970,43350} }; // W46

  const int nY6s = 12;
  const int Y6Runs[nY6s][2] = { {50936,51098},   // W35
				{51264,51608},   // W36
				{51750,51986},   // W37
	   {0,0},{0,0},		{52565,52931},   // W40
				{52959,53231},   // W41
				{53347,53468},   // W42
				{53534,53760},   // W43
				{53780,53925},   // W44
				{54126,54356},   // W45
				{54506,54663} }; // W46

  const int nY7s = 18;
  const int Y7Runs[nY7s][2] = { {57898,58191},   // W25 T #0
				{58263,58580},   // W26 T
				{58780,59036},   // W27 T
				{58263,59399},   // W28 T
				{59963,60087},   // W30 T
				{60147,60332},   // W31 T #5

				{60340,60834},   // W32 L #6
				{60839,61122},   // W33 L
				{61168,61252},   // W34 L
				{61459,61678},   // W35 L
				{61717,62040},   // W36 L
				{62095,62367},   // W37 L
				{62398,62564},   // W38 L #12

				{62701,62899},   // W39 T #13
				{62997,63122},   // W40 T
				{63194,63356},   // W41 T
				{63456,63672},   // W42 T
				{63746,63807} }; // W43 T #17

  const int nY8s = 1;
  const int Y8Runs[nY8s][2] = { {67474,71359} }; // W26/39

  const int nY9s = 18;
  //http://wwwcompass.cern.ch/twiki/bin/view/HadronAnalysis/Periods2009
  const int Y9Runs[nY9s][2] = { {74777,75501},   //  0 W24    pi-
				{75511,76162},   //  1 W25    pi-
				{76192,76399},   //  2 W27d   pi-
				{76435,76554},   //  3 W27c   p
				{76556,77213},   //  4 W29    p
				{77239,77734},   //  5 W31    p
				{77759,78530},   //  6 W33    p
				{78531,79147},   //  7 W35    pi-
				{79485,79764},   //  8 W38    mu+ DVCS
				{79773,79793},   //  9 W39    mu- DVCS
				{79794,79964},   // 10 W40    mu+ DVCS
				{80034,80521},   // 11 W42    pi-
				{80608,81093},   // 12 W43    pi-
				{81176,81452},   // 13 W44    pi-
				{81473,81963},   // 14 W45f   pi- Prim
				{81969,82001},   // 15 W45m   mu- Prim
				{82006,82153},   // 16 W45l   pi- Prim
				{82223,82350} }; // 17 W47    DY

  const int nYas = 1;
  const int YaRuns[nYas][2] = { {85197,89276} }; // W24/44

  const int nYbs = 12; // Cf. http://na58dst1.home.cern.ch/na58dst1/dstprod.html
  const int YbRuns[nYbs][2] = { {91737,92168}, 	 //  0 W25
				{92169,92676}, 	 //  1 W27 	
				{92683,92869}, 	 //  2 W30
				{92871,93161}, 	 //  3 W31
				{93162,93325}, 	 //  4 W32
				{93326,93496}, 	 //  5 W33
				{93507,93735}, 	 //  6 W34
				{94089,94577}, 	 //  7 W36
				{94579,94844}, 	 //  8 W38
				{94864,95287}, 	 //  9 W39
				{95326,95766}, 	 // 10 W41
				{95789,96118} }; // 11 W43

  const int nYcs = 31;
  // Cf. https://twiki.cern.ch/twiki/bin/view/Compass/HadronAnalysis/Chronology2012
  const int YcRuns[nYcs][2] = { {102550,102861},   // pi1
				{102878,102942},   // pi2
				{102951,102988},   // mu1
				{103161,103280},   // pi3
				{103514,103530},   // mu2
				{103587,103606},   // mu3
				{103637,103678},   // mu4
				{103682,103804},   // pi4
				{103844,103981},   // pi5
				{103984,104037},   // mu5
				{104105,104177},   // mu6
				{104179,104318},   // pi6
				{104365,104515},   // pi7
				{104518,104528},   // mu7
				{104760,104806},   // mu8
				{104811,104919},   // pi8
				{104957,105091},   // pi9
				{105095,105126},   // mu9
				{105201,105250},   // mua
				{105254,105334},   // pia
				{105351,105415},   // pib
				{105437,105475},   // mub
				{105481,105486},   // mup
				{105493,105596},   // piD in W37
				{105609,105688},   // piC in W38
				{105691,105692},   // piE in W38
				{107881,108205},   // W44 : DVCS
				{108213,108431},   // W45
				{108440,108567},   // W46
				{108581,108902},   // W47
				{108910,109081} }; // W48

  const int nYes = 3;
  // From "generalprod/testcoral/dy2014t(1|2|3)"
  const int YeRuns[nYes][2] = { {254440,254619},   // T05
				{254682,254899},   // T06
				{255037,255214} }; // T07
  const int nYfs = 9;
  // From "generalprod/testcoral/dy15W??t3"
  const int YfRuns[nYfs][2] = { {259363,260016},   // W07
				{260074,260565},   // W08
				{260627,261496},   // W09
				{261515,262221},   // W10
				{262370,263090},   // W11
				{263143,263603},   // W12
				{263655,264134},   // W13
				{264170,264562},   // W14
				{264619,264857} }; // W15

  const int nY10s = 13;
  // From "/castor/.../2016/raw"
  // (and https://twiki.cern.ch/twiki/bin/view/Compass/GPD/Runs2016
  const int Y10Runs[nY10s][2]={ {269874,270099},   // W04
				{269874,270399},   // W05
				{270400,271149},   // W06
				{271150,271699},   // W07
				{271700,272329},   // W08
				{272340,273116},   // W09
				{273122,273654},   // W10
				{273655,274329},   // W11
				{274330,274909},   // W12
				{274910,275430},   // W13
				{275475,275908},   // W14
				{275910,276368},   // W15
				{276381,276524} }; // W16

  const int nY11s = 9;
  // From "http://wwwcompass.cern.ch/compass/run/run2017/periods.html"
  const int Y11Runs[nY11s][2]={ {278474,278706},   // W04
				{278723,279158},   // W05
				{279209,279429},   // W06
				{279522,279796},   // W07
				{279860,280203},   // W08
				{280310,280619},   // W09
				{280630,281065},   // W10
				{281099,281365},   // W11
				{281380,281722} }; // W12

  if (Y2Runs[0][0]<=run && run<=Y2Runs[nY2s-1][1]) {
    year = 2002;                                   // ********** 2002 **********
    for (int i = 0; i<nY2s; i++)
      if      (Y2Runs[i][0]<=run && run<=Y2Runs[i][1]) {
	iPeriod = i;
	strtRun = Y2Runs[iPeriod][0]; stopRun = Y2Runs[iPeriod][1];
	switch (iPeriod) {
	case 0:  meanIndex = 1.001124; break;
	case 1:  meanIndex = 1.001124; break;
	case 2:  meanIndex = 1.001144; break;
	case 3:  meanIndex = 1.001263; break;
	case 4:  meanIndex = 1.001415; break;
	case 5:  meanIndex = 1.001417; break;
	default: meanIndex = 1.001460;
	}
	static const char *Y2Periods[nY2s] = {"P1C","P2A","P2D","P2E","P2F","P2G"};
	period = new string(Y2Periods[iPeriod]);

      }
    if (*period=="P12B" || *period=="P2C" || *period=="P2H")
      physicsType = new string("Transverse");
    else
      physicsType = new string("Longitudinal");
  }
  else if (Y3Runs[0][0]<=run && run<=Y3Runs[nY3s-1][1]) {
    year = 2003;                                  //  ********** 2003 **********
    for (int i = 0; i<nY3s; i++)
      if  (Y4Runs[i][0]<=run && run<=Y4Runs[i][1]) {
	iPeriod = i;
	meanIndex = 1.001460;
	strtRun = Y3Runs[iPeriod][0]; stopRun = Y3Runs[iPeriod][1];
	static const char *Y3Periods[nY3s] = {"P1A","P1B","P1C","P1D","P1E","P1F","  ","  ","P1I","P1J"};
	period = new string(Y3Periods[iPeriod]);
	break;
      }
    if (*period=="P1G" || *period=="P1H")
      physicsType = new string("Transverse");
    else
      physicsType = new string("Longitudinal");
  }
  else if (Y4Runs[0][0]<=run && run<=Y4Runs[nY4s-1][1]) {
    year = 2004;                                  //  ********** 2004 **********
    for (int i = 0; i<nY4s; i++)
      if  (Y4Runs[i][0]<=run && run<=Y4Runs[i][1]) {
	iPeriod = i;
	meanIndex = 1.001460;
	strtRun = Y4Runs[iPeriod][0]; stopRun = Y4Runs[iPeriod][1];
	static const char *Y4Periods[nY4s] = {
	  "   ","W21","W22","W23","   ","   ","W26","W27","W28","W29",
	  "W30","W31","W32","W33","W34","W35","W36","W37","W38","W39",
	  "W40","   ","   ","   ","W44","W45","W46"};
	period = new string(Y4Periods[iPeriod]);
	break;
      }
    int periodN = atoi(period->c_str()+1);
    if (33<=periodN && periodN<=36)
      physicsType = new string("Transverse");
    else
      physicsType = new string("Longitudinal");
  }
  else if (Y6Runs[0][0]<=run && run<=Y6Runs[nY6s-1][1]) {
    year = 2006;                                  //  ********** 2006 **********

    // - Infra is valid for mass production of mid-2007
    // - W40 has to be split into 2 sub-periods, W44 into 3 => arrays[nY6s*3]
    int W40Period = 5, W40Runs[2][2] = { {52565,52635}, {52671,52931}};
    int W44Period = 9, W44Runs[3][2] = { {53780,53805}, {53815,53846}, {53860,53925}};

    for (int i = 0; i<nY6s; i++)
      if  (Y6Runs[i][0]<=run && run<=Y6Runs[i][1]) {
	iPeriod = i;
	if      (iPeriod==W40Period) {
	  if      (W40Runs[0][0]<=run && run<=W40Runs[0][1]) subPeriod = 0;
	  else if (W40Runs[1][0]<=run && run<=W40Runs[1][1]) subPeriod = 1;
	  else {
 printf("** periodDB: #%d (06W40) in neither [#%d,#%d] nor [#%d,#%d]\n",
	run,W40Runs[0][0],W40Runs[0][1],W40Runs[1][0],W40Runs[1][1]); exit(1);
	  }
	}
	else if (iPeriod==W44Period) {
	  if      (W44Runs[0][0]<=run && run<=W44Runs[0][1]) subPeriod = 0;
	  else if (W44Runs[1][0]<=run && run<=W44Runs[1][1]) subPeriod = 1;
	  else if (W44Runs[2][0]<=run && run<=W44Runs[2][1]) subPeriod = 2;
	  else {
 printf("** periodDB: #%d (06W44) in neither [#%d,#%d] nor [#%d,#%d]\n",
	run,W44Runs[0][0],W44Runs[0][1],W44Runs[1][0],W44Runs[1][1]); exit(1);
	  }
	}
	else subPeriod = 0;
	meanIndex = 1.001350; // For PMTs
	strtRun = Y6Runs[iPeriod][0]; stopRun = Y6Runs[iPeriod][1];
	static const char *Y6Periods[nY6s] = {
 "W35","W36","W37","W38","W39","W40","W41","W42","W43","W44","W45","W46"};
	period = new string(Y6Periods[iPeriod]);
	break;
      }
    physicsType = new string("Longitudinal");
  }
  else if (Y7Runs[0][0]<=run && run<=Y7Runs[nY7s-1][1]) {
    year = 2007;                                  //  ********** 2007 **********
    for (int i = 0; i<nY7s; i++)
      if  (Y7Runs[i][0]<=run && run<=Y7Runs[i][1]) {
	iPeriod = i;
	strtRun = Y7Runs[iPeriod][0]; stopRun = Y7Runs[iPeriod][1];
	static const char *Y7Periods[nY7s] = {      
	  "W25","W26","W27","W28","W30","W31",
	  "W32","W33","W34","W35","W36","W37","W38",
	  "W39","W40","W41","W42","W43" };
	period = new string(Y7Periods[iPeriod]);
	break;
      }
    if (6<=iPeriod && iPeriod<=12)
      physicsType = new string("Longitudinal");
    else
      physicsType = new string("Transverse");
  }
  else if (Y8Runs[0][0]<=run && run<=Y8Runs[nY8s-1][1]) {
    year = 2008;                                  //  ********** 2008 **********
    iPeriod = 0;
    strtRun = Y8Runs[iPeriod][0]; stopRun = Y8Runs[iPeriod][1];
    physicsType = new string("Hadron");
  }
  else if (Y9Runs[0][0]<=run && run<=Y9Runs[nY9s-1][1]) {
    year = 2009;                                  //  ********** 2009 **********
    for (int i = 0; i<nY9s; i++)
      if  (Y9Runs[i][0]<=run && run<=Y9Runs[i][1]) {
	iPeriod = i;
	strtRun = Y9Runs[iPeriod][0]; stopRun = Y9Runs[iPeriod][1];
	static const char *Y9Periods[nY9s] = {      
	  "W24", "W25","W27d","W27c","W29", "W31", "W33", "W35","W38","W39",
	  "W40", "W42","W43", "W44", "W45f","W45m","W45l","W47"}; 
	period = new string(Y9Periods[iPeriod]);
	break;
      }
    if (*period=="W45m") {
      physicsType = new string("Primakoff_mu");
    }
    else if (iPeriod<8 || 10<iPeriod) {
      physicsType = new string("Hadron");
    }
  }
  else if (YaRuns[0][0]<=run && run<=YaRuns[nYas-1][1]) {
    year = 2010;                                  //  ********** 2010 **********
    iPeriod = 0;
    strtRun = YaRuns[iPeriod][0]; stopRun = YaRuns[iPeriod][1];
    physicsType = new string("Transverse");
  }
  else if (YbRuns[0][0]<=run && run<=YbRuns[nYbs-1][1]) {
    year = 2011;                                  //  ********** 2011 **********
    for (int i = 0; i<nYbs; i++)
      if  (YbRuns[i][0]<=run && run<=YbRuns[i][1]) {
	iPeriod = i;
	strtRun = YbRuns[iPeriod][0]; stopRun = YbRuns[iPeriod][1];
	static const char *YbPeriods[nYbs] = {
	  "W25","W27","W30","W31","W32","W33","W34","W36","W38","W39",
	  "W41","W43"};
	period = new string(YbPeriods[iPeriod]);
	// Mean indices evaluated by averaging w/ mock fluxes, cf. ~/prods/2011
	static const double YbMeanIndex[nYbs] = {
	  1.001297,1.001293,1.001296,1.001296,1.001294,1.001292,0,1.001303,1.001313,1.001320,
	  1.001325,1.001317};
	meanIndex = YbMeanIndex[iPeriod];
	break;
      }
    physicsType = new string("Longitudinal");
  }
  else if (YcRuns[0][0]<=run && run<=YcRuns[nYcs-1][1]) {
    year = 2012;                                  //  ********** 2012 **********
    int i, isDVCS; for (i=isDVCS = 0; i<nYcs; i++)
      if  (YcRuns[i][0]<=run && run<=YcRuns[i][1]) {
	iPeriod = i;
	strtRun = YcRuns[iPeriod][0]; stopRun = YcRuns[iPeriod][1];
	static const char *YcPeriods[nYcs] = {
	  "pi1","pi2","mu1","pi3","mu2","mu3","mu4","pi4","pi5","mu5",
	  "mu6","pi6","pi7","mu7","mu8","pi8","pi9","mu9","mua","pia",
	  "pib","mub","mup","piD","piC","piE","W44","W45","W46","W47",
	  "W48"};
	period = new string(YcPeriods[iPeriod]);
	break;
      }
    if (period->c_str()[0]=='W')
      physicsType = new string("DVCS");
    else {
      if      (period->find("pi")==0)
	physicsType = new string("Hadron");
      else if (period->find("mu")==0)
	physicsType = new string("Primakoff_mu");
    }
  }
  else if (YeRuns[0][0]<=run && run<=YeRuns[nYes-1][1]) {
    year = 2014;                                  //  ********** 2014 **********
    for (int i = 0; i<nYes; i++)
      if  (YeRuns[i][0]<=run && run<=YeRuns[i][1]) {
	iPeriod = i;
	strtRun = YeRuns[iPeriod][0]; stopRun = YeRuns[iPeriod][1];
	static const char *YePeriods[nYes] = {      
	  "T05","T06","T07"};
	period = new string(YePeriods[iPeriod]);
	break;
      }
    physicsType = new string("DY");
  }
  else if (YfRuns[0][0]<=run && run<=YfRuns[nYfs-1][1]) {
    year = 2015;                                  //  ********** 2015 **********
    for (int i = 0; i<nYfs; i++)
      if  (YfRuns[i][0]<=run && run<=YfRuns[i][1]) {
	iPeriod = i;
	strtRun = YfRuns[iPeriod][0]; stopRun = YfRuns[iPeriod][1];
	static const char *YfPeriods[nYfs] = {
	  "W07","W08","W09","W10","W11","W12","W13","W14","W15"};
	period = new string(YfPeriods[iPeriod]);
	break;
      }
    physicsType = new string("DY");
  }
  else if (Y10Runs[0][0]<=run && run<=Y10Runs[nY10s-1][1]) {
    year = 2016;                                  //  ********** 2016 **********
    for (int i = 0; i<nY10s; i++)
      if  (Y10Runs[i][0]<=run && run<=Y10Runs[i][1]) {
	iPeriod = i;
	strtRun = Y10Runs[iPeriod][0]; stopRun = Y10Runs[iPeriod][1];
	switch (iPeriod) { // From sumProd and UserEvent109
	  // ~/phast.utils/sumProd.csh Phi_8910/h-27.....log $rMn $rMx
	  // ./phast -u 109 -U 1 -T 16 -T P09 ...
	case 8:   meanIndex = 1.001434; break; // rMn,rMx 274509 274901
	case 9:   meanIndex = 1.001451; break; // rMn,rMx 274946 275393 
	case 10:  meanIndex = 1.001452; break; // P09: rMn,rMx 275478 275907
	case 11:  meanIndex = 1.001458; break; // rMn,rMx 276014 276318
	default:  meanIndex = 1.001460;
	}
	static const char *Y10Periods[nY10s] = {
	  "W04","W05","W06","W07","W08","W09","P05",
	  "P06","P07","P08","P09","P10","P11"};
	period = new string(Y10Periods[iPeriod]);
	break;
      }
    physicsType = new string("DVCS");
  }
  else if (Y11Runs[0][0]<=run && run<=Y11Runs[nY11s-1][1]) {
    year = 2017;                                  //  ********** 2017 **********
    for (int i = 0; i<nY11s; i++)
      if  (Y11Runs[i][0]<=run && run<=Y11Runs[i][1]) {
	iPeriod = i;
	strtRun = Y11Runs[iPeriod][0]; stopRun = Y11Runs[iPeriod][1];
	static const char *Y11Periods[nY11s] = {
	  "W04","W05","W06","W07","W08","W09","W10",
	  "W11","W12"};
	period = new string(Y11Periods[iPeriod]);
	break;
      }
    physicsType = new string("DVCS");
  }

  if (!year) {
    // Many (not all) of the settings managed by "DataTakingDB" rely on "Year".
    // => Let's abort if it cannot be defined
    printf("\n** periodDB:\a Year of data taking cannot be determined from built-in DB!\n\n");
    return false;
  }
  if (!period || !physicsType) {
    // If "year" is set, these two variables should also be set.
    printf("\n** periodDB:\a Inconsistency: no \"%s\"\n\n",
	   period?"physicsType":"period");
    return false;
  }
  return true;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   settingDB   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void settingDB(int year, const string &period, const string &physicsType,
	       int &thCIndex,
	       double &pBeam, double &dpBMS, int &beamRecMethod,
	       double &dqP2BMSCut,
	       double &yLowCut, double &yUpCut, double &chi2Cut)
{
  // Defines a number of settings, depending upon input "year" and "physicsType"
  thCIndex // "7" is the best option, but...
    // ..it's not been enabled in latest mass productions: to gain in speed
    = year<2008 ? 7 : 8;
  // Override: Following Marcin's line of thought: ring fit offers the
  // opportunity to be independent of the Likelihood calculation.
  thCIndex = 10;
  pBeam = 160; dpBMS = 20;
  beamRecMethod = 1; // Default beam reco method is new, Pawel's, one
  // y cut:
  // - These are minimal cuts, guaranteeing a reliable reconstruction.
  //   Even looser would be something intended to explore the peak at y=1. Then,
  //  better cancel cut alltogether.
  // - Special, hadron data, cases: values below are updated.
  yLowCut = .05; yUpCut = 0.95;
  // Beams w/o BMS:
  // - Can be singled out by looking at the uncertainty on their momentum: it's
  //  to be larger than usual. How large is this? This is assigned by TraFDic
  //  ("TEv::TracksFit2"), according to TraF option dCut[5]. This, in turn, is
  //  expected to correspond to the beam momentum spread, which is typically
  //  15GeV/sqrt(12)
  //    => d(1/P) = 15/sqrt(12.)/160/160 ~= 1.69e-4, cov(1/P,1/P) ~= 2.86e-8
  //  Whereas BMS precision on P is i) variable, ii) typically .5%
  //    => d(1/P) ~= .005/160            ~= 3.13e-5, cov(1/P,1/P) ~= 9.77e-10
  //   This is consistent w/ what SIDIS analysis uses, as reported by Ana Sofia:
  //   On 11/03/2012 05:07 PM, Ana Sofia Da Silva Nunes wrote:
  //  > In the latest SIDIS release note, a cut of \sigma^2_{q/|P|}<20.e-9 GeV^-2
  // - But in some cases, the actual cut on cov(1/P,1/P) may have to be tigher,
  //  e.g. in 07W45. But how far thighter can we go, w/o cutting away also some
  //  genuine BMS, those w/ worst precision?
  //   On 11/03/2012 05:07 PM, Ana Sofia Da Silva Nunes wrote:
  //  > In any case, what removes the spike in the beam momentum of 07W45 (at
  //  > least for this run) is actually the cut \sigma^2_{q/|P|}>0.5e-9 GeV^-2
  //   And this seems to be definitely too tight
  // => Have to investigate.
  // => In the mean time, let's set the cut @ 20.e-9, except for 2007W45
  dqP2BMSCut = 20.e-9;
  if      (2002<=year && year<=2006) {
    beamRecMethod = 0; // 2002,6: last mass prod. uses old beam
  }
  else if (year==2007) {
    // Longitudinal 2007: Last mass prod. uses old beam, and is most probably
    // unreliable, given that it's similar to 2007 transverse penultimate mass
    // prod., where this bad reliability was evidenced (by the observation of
    // a high rate of badly negative missing energy events). Also, the
    // backpropagation info (PaParticle::Chi2CutFlag, in the present case) is
    // unreliable: let's then not use it, cf. "DatatakingDB::BadBackProp".
    if (physicsType!="Transverse") beamRecMethod = -1;
    if (period=="W45") {
      dqP2BMSCut = 0.5e-9;
      printf(" * settingDB:\a 2007W45: exceptionally set q(1/P)² = %.2e\n",
	     dqP2BMSCut);
    }
  }
  else if (year==2011) {
    pBeam = 200; dpBMS = 25;              // 200 GeV mu beam in 2011
  }
  else if (physicsType=="Hadron") {
    pBeam = 190; dpBMS = 0; yLowCut=yUpCut = 0;       // Not a muon beam
  }
  else if ((year==2009 || year==2012) && physicsType=="Primakoff_mu") {
    pBeam = 190; yLowCut = 0.3; yUpCut = 0.95; // 190 GeV mu beam in Primakoff
  }

  // CHI2 CUT and PHYSICS TYPE
  // - Chi2 Cut may be allowed to be tighter for some particular years, as
  //  precribed by the info collected in "periodDB".
  // - This tight setting can be bypassed by "chi2CutOption=1".
  chi2Cut = 10;
  if      (2002<=year && year<=2004) chi2Cut = 8; // 2002:4 newest productions
  else if (2008<=year && year<=2009 ||
	   year==2012) {
    chi2Cut = 8; // No target field => easier alignment
  }

}
#include "TSystemDirectory.h"
#include "TList.h"
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GetBadSpills ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  BadSpillPat  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
unsigned int DataTakingDB::BadSpillPat(int run, int spill)
{
  unsigned int rs = (run<<12)+spill;
  // There is room for THREE bad spill lists. Yet, most of the time only one
  // is available. Then, by default, the corresponding OK bits are set.
  set<unsigned> *bs[] = {&BadSpills, &BadSpillsCalo, &BadSpillsRICH};
  int ibs; unsigned int spillPat; for (ibs = 0, spillPat = 0x7; ibs<3; ibs++) { 
    set<unsigned> *b = bs[ibs];
    if (b->find(rs)!=b->end()) spillPat ^= 1<<ibs;
  }
  return spillPat;
}
void DataTakingDB::GetBadSpills()
{
  GetBadSpills(BadSpills,BadSpillsCalo,BadSpillsRICH);
}
void DataTakingDB::GetBadSpills(set<unsigned> &badSpills,
				set<unsigned> &badSpillsCalo,
				set<unsigned> &badSpillsRICH)
{
  badSpills.clear(); badSpillsCalo.clear(); badSpillsRICH.clear();
  // BAD SPILL LIST EXPECTED TO BE STORE in "~/phast.utils/badSpills/<year>"
  string sDir(gSystem->ExpandPathName("~/phast.utils/badSpills/"));
  char cYear[] = "2006"; sprintf(cYear,"%d",Year); sDir += string(cYear);
  sDir += string("/");
  printf(" * DataTakingDB: Searching \"%s\" for Bad Spill Lists relevant for period \"%s\"...",sDir.c_str(),Period?Period->c_str():"?");
  TSystemDirectory d(".",sDir.c_str()); // Working directory
  TList *liste = d.GetListOfFiles(); if (!liste) {
    printf("\n \"%s\" does not exist or is empty\n\n",sDir.c_str());
    abort();
  }
  // I) badSpills: LOOKING FOR FILES OF TYPE, e.g., "BadSpill_06W45..."
  //II) badSpillsCalo/RICH: LOOKING FOR, e.g., "BadSpill_06W45...Calo/RICH..."
  string sBadSpill("BadSpill_");
  if (Period) { sBadSpill += string(cYear+2); sBadSpill += *Period; }
  else          sBadSpill += string(cYear);
  printf("\n                targeting \"%s*\" files...",sBadSpill.c_str());
  const char *cBadSpill = sBadSpill.c_str(); int nChars = strlen(cBadSpill);
  char newline = '\n';
  vector<string> infiles[3]; TObject *obj = liste->First(); do {
    const char *cFName = obj->GetName();
    if (!strncmp(cFName,cBadSpill,nChars)) {
      int ib = 0; if (strstr(cFName,"Calo")) ib = 1;
      else if        (strstr(cFName,"RICH")) ib = 2;
      if (!infiles[ib].empty()) {
	printf("%c  ... Found TWO different files \"%s\" and \"%s\". Which is valid?\n\n",
	       newline,infiles[ib][0].c_str(),cFName); abort();
      }
      string sFPath = sDir+string(cFName); infiles[ib].push_back(sFPath);
      printf("%c  ... Found \"%s\"\n",newline,infiles[ib][0].c_str()); newline = '\0';
      bool not_badspl = false;
      switch (ib) {
      case 0: PaUtils::GetBadSpills(infiles[0],badSpills,not_badspl); break;
      case 1: PaUtils::GetBadSpills(infiles[1],badSpillsCalo,not_badspl); break;
      case 2: PaUtils::GetBadSpills(infiles[2],badSpillsRICH,not_badspl); break;
      }
    }
  } while ((obj = liste->After(obj)));
  if (infiles[0].empty()) {
    printf("%c  ... Found no plain \"%s*\" file.\n\n",newline,cBadSpill);
  }
  else printf("\n");
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~     XCheckYearSetup     ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool XCheckYearSetup(int year)
{
  // X-check what we can guess as to the year of dataking from the target
  // position in PaSetup against argument "year".
  const PaSetup &setup = PaSetup::Ref();
  bool error = false; string tgtCenterType;
  if (fabs(setup.TargetCenterZ())<5) {             // 2006/7...
    // ...Design value is 0. Real value may be shifted by 2 to 3 cm, cf.
    // "PaAlgo::GetTargetLocation"
    tgtCenterType = "2006";
    error = year!=2006 && year!=2007 && year!= 2010 && year != 2011;
  }
  else if (fabs(setup.TargetCenterZ()+35)<5) {     // 2002:2004
    tgtCenterType = "2002"; error = year<2002 || 2004<year;
  }
  else if (fabs(setup.TargetCenterZ()+48.52)<5) {  // LH2 2008,9
    tgtCenterType = "LH2"; error = year!=2008 && year!=2009;
  }
  else if (fabs(setup.TargetCenterZ()+192.43)<5) { // LH2 2012
    tgtCenterType = "DVCS"; error = year!=2012;
  }
  else if (fabs(setup.TargetCenterZ()+198.50)<5) { // LH2 2016
    tgtCenterType = "DVCS"; error = year!=2016;
  }
  else if (fabs(setup.TargetCenterZ()+230.00)<5) { // DY
    tgtCenterType = "DY"; error = year!=2014 && year!=2015;
  }
  else {
 printf("\n** DataTakingDB: PaSetup::TargetCenter =%.3fcm: none of \"2002\", \"2006\", \"LH2\", \"DVCS\"\n",
	setup.TargetCenterZ()); return false;
  }
  if (error) {
 printf("\n** DataTakingDB: PaSetup::TargetCenter =%.3fcm, i.e. \"%s\"-like, inconsistent w/ year = %d\n",
	setup.TargetCenterZ(),tgtCenterType.c_str(),year); return false;
  }
  return true;
}
// **********************************************************************
// *************************    GetTrigType    **************************
// **********************************************************************
int GetTrigType(unsigned int evTrig)
{
  int trigType = -1; // Type=ilmoLc; M=any middle, C=pure C, precedence=ilmoL
  DataTakingDB *DTDB = DataTakingDB::Ptr();
  if      (DTDB->Year<2010) { // <2010: no LAST yet
    if      (evTrig&0x1)                trigType = 0; // useless for transverse
    else if (evTrig&0x4)                trigType = 1;
    else if ((evTrig&0x2)   /* SI */ ||
	     (evTrig&0x100) /* I  */)   trigType = 2;
    else if (evTrig&0x8)                trigType = 3;
    else if (evTrig&0x10) /* Pure C */  trigType = 5;
  }
  else if (DTDB->Year<2012) { // 2010:2011: w/ LAST
    if      (evTrig&0x1)                trigType = 0; // useless for transverse
    else if (evTrig&0x4)                trigType = 1;
    else if ((evTrig&0x2)   /* SI */ ||
	     (evTrig&0x100) /* I  */)   trigType = 2;
    else if (evTrig&0x8)                trigType = 3;
    else if (evTrig&0x200)              trigType = 4;
    else if (evTrig&0x10) /* Pure C */  trigType = 5;
  }
  else if (DTDB->Year==2012 || // 2012: no IT, nor I MT
	   DTDB->Year==2016 || DTDB->Year==2017) {
    if      (evTrig&0x4)                trigType = 1;
    else if (evTrig&0x2)                trigType = 2;
    else if (evTrig&0x8)                trigType = 3;
    else if (evTrig&0x200)              trigType = 4;  // LAST
    else if (evTrig&0x10) /* Pure C */  trigType = 5;
    else if (evTrig&0xc0)               trigType = 6;
  }
  return trigType;
}

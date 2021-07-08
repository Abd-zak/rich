// $Id: GetPIDs.cc,v 1.4 2021/04/05 12:50:38 ybedfer Exp $

// Fill a number of arrays (of size = to # of tracks in Evt) for PID

// $Log: GetPIDs.cc,v $
// Revision 1.4  2021/04/05 12:50:38  ybedfer
// Cosmetics...
//
// Revision 1.3  2020/07/29 17:12:31  ybedfer
// Use last helix in vTPar to determine Zlast.
//
// Revision 1.2  2013/08/08 21:15:52  ybedfer
// - Added input argument "thCIndex" = index of estimator for theta Cerenkov.
// - Returned values.
//   - Array "PIDs": 0, 0x8 or 0x10. And no longer 0x20.
//   - LHs: if small => set it =0.
// - Zones.
//   - Get ZLast from last helix.
//   - Check that the array of helices is such that last helix can be taken as
//    array[1].
//
// Revision 1.1  2012/11/05 20:26:13  ybedfer
// Initial revision
//

// The working of the routine depends upon the following cpp macros
#define PaTrack_RICHDATASIZE 21

#include "Phast.h"
#include "RICH.h"
#include "PaSetup.h"

void GetPIDs(PaEvent &e, RICH *rich,
	     double runIndex,   // Best estimate of index
	     double prodIndex,  // Index value used in production
	     double deltaIndex, // Diff w.r.t. mean index in period
	     int thCIndex,      // Index of estimator for theta Cerenkov
	     // Pointers to arrays expected to be allocated size(#tracks)
	     int *PIDs, // Expected = 0 on input; 0x8: RICH blcok exists, 0x10: bckLH != 0
	     double **richLHs, double **richChi2s,
	     int *tZones,
	     double *richDths)  // dtheta/dn/d1/beta * DeltaN
{

#ifdef DEBUG_PIDs
  const double rRICH = 5, ZRICH = 750, r2RICH = rRICH*rRICH, RRICH = 300;
  char cIni = '\n';
#endif

  const PaSetup &setup = PaSetup::Ref();
  static double ZSM1, ZSM2;

  static bool first = true; if (first) {
    //                ***** GET MAGNET ABSCISSAE *****
    PaMagInfo *m = setup.PtrMagField()->getMagInfo();
    int nMags = setup.PtrMagField()->getNumOfMags(), iSM2 = nMags-1;
    int iSM1 = iSM2-1; ZSM1 = m[iSM1].zcm/10; ZSM2 = m[iSM2].zcm/10;
    first = false;
  }

  for (int iET = 0; iET<(int)e.vTrack().size(); iET++) {

    // ********** LOOP on ALL TRACKS => RICH pID

    const PaTrack &trk = e.vTrack()[iET];
    const PaTPar &hi = trk.vTPar()[0];   // Initial helix
    double mom = hi.Mom();

    // ********** RICH **********
    //0) BACKGROUND Likelihood
    //1) PION Likelihood
    //2) KAON Likelihood
    //3) PROTON Likelihood
    //4) PION Likelihood Derivative
    //5) KAON Likelihood Derivative
    //6) PROTON Likelihood Derivative
    //7) THETAL Max. Likelihood Cherenkov angle
    //8) THETAR Reconstructed Cherenkov angle
    //9) N Photons in the Ring
    //10) THETAF Cherenkov angle from Ring Fit
    //11) Ring Chisquare
    //12) PION Chisquare
    //13) KAON Chisquare
    //14) PROTON Chisquare
    //15] electron Like
    //16] muon Like
    //17] electron Like derivative
    //18] muon Like derivative
    //19] electron chisquare
    //20] muon chisquare

    if (trk.NRichInf()>=PaTrack_RICHDATASIZE) {
      int iInf, ipiKp;
      PIDs[iET] |= 0x8; // 0x8 <-> Something about RICH exists

      if (richDths) {
	double theta = trk.RichInf(thCIndex), sth = sin(theta*0.001);
	richDths[iET] = 1000*deltaIndex/runIndex/runIndex/sth;
      }
      double bckLH = trk.RichInf(0);
      if (fabs(bckLH)<1e-6) // It can happen that the LHs come close to 0...
	// ... this behaviour is prone to make the analysis depend upon context.
	// => Let's set =0 all such cases
	bckLH = 0;
      if (bckLH>0) {

	// Correct Likehoods and Chi2's from error on index
	// Store result in arrays "richLHs" and "richChi2s"

	PIDs[iET] |= 0x10;                     // Bckgrd likelihood available
	
	for (iInf = 1, ipiKp = 0; iInf<=3; iInf++, ipiKp++) {
	  // ***** RICH pi, K AND p LH *****
	  // Likelihood pi (resp. K,p) corrected and divided
	  // by likelihood background
	  double LH = trk.RichInf(iInf);
	  if (fabs(LH)<1e-6) LH =0; // Cf. comment supra in the bckLH case.
	  if (LH) {
	    LH -= (runIndex-prodIndex)*trk.RichInf(iInf+3);
	    richLHs[ipiKp][iET] = LH/bckLH;
	  }
	  else richLHs[ipiKp][iET] = 0;
	  if (richChi2s) richChi2s[ipiKp][iET] = trk.RichInf(iInf+11);
	}
	if (richChi2s)
	  rich->CorrChiSq(richChi2s[0][iET],richChi2s[1][iET],richChi2s[2][iET],
			  trk.RichInf(11),trk.RichInf(10),mom);
#if PaTrack_RICHDATASIZE == 21
	int iemu; for (iInf = 15, iemu = 3; iInf<=16; iInf++, iemu++) {
	  // ***** RICH e AND mu *****
	  double LH = trk.RichInf(iInf);
	  if (fabs(LH)<1e-6) LH =0; // Cf. comment supra in the bckLH case.
	  if (LH) {
	    LH -= (runIndex-prodIndex)*trk.RichInf(iInf+2);
	    richLHs[iemu][iET] = LH/bckLH;
	  }
	  else richLHs[iemu][iET] = 0;
	  if (richChi2s) richChi2s[iemu][iET] = trk.RichInf(iInf+4);
	}
#endif
      }  // End loop on infos
    }  // End RICH info
#ifdef DEBUG_PIDs
    else {
      // Extrapolate to RICH in view of checking whether w/in RICH acceptance
      // (Note: this is also done in "UserEvent3/transport".)
      PaTPar Hout; hi.Extrapolate(ZRICH,Hout,false);
      double XR = Hout.Pos(0), YR = Hout.Pos(1);
      if (XR*XR+YR*YR>r2RICH && fabs(XR)<RRICH*1.5 && fabs(YR)<RRICH*1.5) {
	int win = XR*XR+YR*YR>r2RICH*1.5 && fabs(XR)<RRICH && fabs(YR)<RRICH;
	if (trk.ZLast()>setup.TargetCenterZ()) 
	  printf("%cDEBUG_PIDs: Run %d Evt %7d Track %2d %6.1fGeV %4.0f<->%4.0fcm %6.1f,%6.1fcm => %s: No (%d) RICH\n",cIni,e.RunNum(),(int)e.UniqueEvNum(),iET,mom,hi(0),trk.ZLast(),XR,YR,win?"IN":"??",trk.NRichInf());
	cIni = '\0';
      }
    }
#endif

    int nHelices = trk.NTPar(); double zL;
    if (nHelices>1) {
      //  If #helices >1, last helix must be available => use it to access
      // ZLast, for it's faster and hence interesting since we are doing for all tracks.
      const vector<PaTPar> &hs = trk.vTPar();
      /*
	// Check which is the very last helix: is it [1] or last in vTPar?
	// => Conclusion: last helix in vTPar, which is implemented infra.
      static int nHsCheck = 2; if (nHelices>nHsCheck) {
	//  Not sure that PTPar [0] and [1] are first and last (or the reverse
	// for those strict beam tracks that ends updtream of target) => check
	// it upon 1st occurence.
	double zMx = hs[0](0)>hs[1](0) ? hs[0](0) : hs[1](0);
	for (int i = 2; i<nHelices; i++) {
	  if (hs[i](0)>zMx+.1) {
 printf("\n ** GetPIDs:\a Evt %d,%d Track %d: %d PaTPar's h%d(=%.2f) > h0(=%.2f),h1(=%.2f)\n\n",
	    e.RunNum(),(int)e.UniqueEvNum(),iET,
	    nHelices,i,hs[i](0),hs[0](0),hs[1](0));
            assert(false);
	  }
	  nHsCheck = nHelices;
	}
      }
      */
      zL = hs[nHelices-1](0);
    }
    else zL = trk.ZLast();
    // Zones: 0x1 -SM1- 0x2 -SM2- 0x4 -muWall- 0x8
    int zones = 0;
    if      (hi(0)<ZSM1) { zones |= 0x1; if (zL>ZSM1) zones |= 0x2; }
    else if (hi(0)<ZSM2)   zones |= 0x2;
    if      (zL>ZSM2) {    zones |= 0x4; if (zL>3800) zones |= 0x8; }
    tZones[iET] = zones;
    if (zones==0x1) continue;                            // Exclude fringe field

  } // ********** END LOOP ON PARTICLES IN EVENT **********

}

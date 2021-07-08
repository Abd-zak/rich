// $Id: RICH.cc,v 1.2 2021/04/05 13:01:48 ybedfer Exp ybedfer $

// $Log: RICH.cc,v $
// Revision 1.2  2021/04/05 13:01:48  ybedfer
// No change.
//
// Revision 1.1  2012/11/05 20:28:13  ybedfer
// Initial revision
//
// Revision 1.1  2004/02/11 14:36:10  ybedfer
// Initial revision
//

#include <stdio.h>
#include <vector>
#include <unistd.h>
#ifndef __CINT__
#  include <wordexp.h>
#endif
#include "PaMetaDB.h"
#include "PaSetup.h"
#include "RICH.h"
#include "Masses.h"

void updateThresholds();

RICH *RICH::address = 0; // Init static pointer

RICH::RICH(double nSigmas, double nVeto)
{

  if (address) {
    printf("** RICH::RICH:\a RICH already instantiated\n");
    assert(false);
  }
  address = this;

  // ********** pID **********

  // ***** SET pID CUTS *****
  //! dtheta(p) = [0]+(x<[1])*[2]*(x-[1])*(x-[1]); x = p-threshold.
  const double ResolutionvsP[4] = {.73, 11,.0090, .0053};
  //const double ResolutionvsP[4] = {.73, 11,.0090, .0025};
  for (int i = 0; i<4; i++) if (i!=1) {
    dthvsP[i] = ResolutionvsP[i];
    dthCut[i] = ResolutionvsP[i]*nSigmas; dthVeto[i] = ResolutionvsP[i]*nVeto; 
  }
  dthVeto[1]=dthCut[1]=dthvsP[1] = ResolutionvsP[1];
   

  //M_pi  = 0.13957018; M2_pi = M_pi*M_pi;
  //M_p   = 0.9382723;  M2_p  = M_p *M_p ;
  //M_K   = 0.493677;   M2_K  = M_K *M_K ;

}

void RICH::GetNRICH(double &index, double &indAPV)
{
  index = CurrentIndex; indAPV = CurrentIndAPV;
}

void RICH::UpdateIndices(int run, bool isMC,
			 double &runIndex,  double &prodIndex,
			 double &runIndAPV, double &prodIndAPV) {
  // Retrieve indices (production one from setup and up-to-date one from DB).
  // Update RICH class member indices (and thresholds as well).
  // Return the indices

  //  The following has been cancelled. But it was not checked that it's
  // always useless:
  //#define RICH_USE_BUILTIN_INDEX
  //#if defined RICH_USE_BUILTIN_INDEX
  // Working around a former bug in PHAST failing to convey prod index to mDST
  // Most of the time th mDSTs used to be produced w/ "./pkopt/rich1.2002.opt"
  // => RICHONE   C4F10RefrIndex              1.00129       // mixed N2-C4F10 (20-80%) refraction index  (def.=1.00153)-data           !!!
  //if (e.IsMC()) {
  //  runIndex=prodIndex = 1.00129;
  //  rich->SetNRICH(runIndex);   // In order to update thresholds
  //}
  //else {
  //  printf("** UserEvent1:\a \"RICH_USE_BUILTIN_INDEX\" defined while !MC\n");
  //  exit(1);
  //}

  // Last run of 04W46. Used to switch from the 1-index case (2002:2004) to
  // the 2-index case (2006...) . =43350 = last run in DB. =41397 = last run
  // in PHAST's RICH index DB.
  const int last2004Run = 43400;
  const PaSetup &setup = PaSetup::Ref();
  if (setup.vNminusOne().size()==2) {
    prodIndAPV = setup.vNminusOne()[0]; prodIndex = setup.vNminusOne()[1];
  }
  else if (setup.NminusOne()<0) {
    printf("** RICH:\a No RICH production index for run #%d!\n",run);
    if (isMC) {
      prodIndAPV = 1.00153 /* cf. rich1.2004.opt */; prodIndex = -1;
      printf(" * RICH: Set it = %f!\n",prodIndAPV);
      prodIndAPV = (prodIndAPV-1)*1e6;
    }
    else exit(1);
  }
  else {
    if (run>last2004Run) {
#define RICH_IGNORE_MISSING_PMT_prodIndex
#ifdef RICH_IGNORE_MISSING_PMT_prodIndex
      prodIndAPV = setup.NminusOne(); prodIndex = -1;
#else
      printf("** RICH:\a No 2 RICH production indices for run #%d!\n",run);
      exit(1);
#endif
    }
    else {
      // APV indices is equated to PMT one => So that PID methods, which are
      // based on "CurrentIndAPV" work also for 2002:4. Same thing for the
      // chi2 correction.
      prodIndex=prodIndAPV = setup.NminusOne();
    }
  }
  if (prodIndex>0)  prodIndex = 1+(double)prodIndex/1e6;
  if (prodIndAPV>0) prodIndAPV = 1+(double)prodIndAPV/1e6;

  if (isMC) {
    runIndex = prodIndex; runIndAPV = prodIndAPV;
    printf(" * RICH:\a MC => run Index = Prod Index\n");
  }
  else {
    const PaMetaDB &db = PaMetaDB::Ref();
    if (db.vNminusOne(run).size()!=2) {
      runIndex = (double)db.NminusOne(run);
      if (runIndex>0) runIndex = 1+runIndex/1000000;
      runIndAPV = runIndex;   // APV indices is equated to PMT one, cf. supra.

    }
    else {
      runIndex  = 1+(double)db.vNminusOne(run)[1]/1000000;
      runIndAPV = 1+(double)db.vNminusOne(run)[0]/1000000;
    }
    if (runIndex<=0) {
      printf("** RICH:\a No RICH DB index for run #%d!\n",run); exit(1);
    }
    if (prodIndex<=0 && run>last2004Run) {
      // This can only happen if "RICH_IGNORE_MISSING_PMT_prodIndex" supra.
      // => Make up 2nd, PMT, production index from DB index
      printf("** RICH:\a No 2 RICH production indices for run #%d!\n",run);
      printf("  => Make up UV one from DB\n");
      prodIndex = runIndex;
    }
  }

  CurrentProdIndex = prodIndex; CurrentProdIndAPV = prodIndAPV;
  CurrentIndex = runIndex; CurrentIndAPV = runIndAPV; UpdateThresholds();

  if (runIndAPV>0)
    printf(" * RICH: RICH indices V,UV = %.5f,%.5f (prod. %.5f,%.5f)\n",
	   runIndex,runIndAPV,prodIndex,prodIndAPV);
  else
    printf(" * RICH: RICH index =%.5f (prod. %.5f)\n",runIndex,prodIndex);
}

void RICH::UpdateThresholds() {

  double *Thrs[5] =    {&eThr,   &muThr,   &piThr,   &KThr,   &pThr};
  double *ThrAPVs[5] = {&eThrAPV,&muThrAPV,&piThrAPV,&KThrAPV,&pThrAPV};
  for (int i = 0; i<5; i++) {
    static double mass;  // "static" for compiler not to complain 
    switch (i) {
    case 0: mass = M_e; break;
    case 1: mass = M_mu; break;
    case 2: mass = M_pi; break;
    case 3: mass = M_K; break;
    case 4: mass = M_p; break;
    }
    double beta = 1/CurrentIndex, gamma = 1/sqrt(1-beta*beta);
    *Thrs[i] = mass*beta*gamma;
    if (CurrentIndAPV>0) {
      beta = 1/CurrentIndAPV, gamma = 1/sqrt(1-beta*beta);
      *ThrAPVs[i] = mass*beta*gamma;
    }
  }
}

void RICH::SetNRICH(double index, double indAPV) {
  CurrentIndex = index; CurrentIndAPV = indAPV; UpdateThresholds();
}
void RICH::SetNRICHProd(double index, double indAPV) {
  CurrentProdIndex = index; CurrentProdIndAPV = indAPV;
}

void RICH::GetNProd(double &index, double &indAPV)
{
  index = CurrentProdIndex; indAPV = CurrentProdIndAPV;
}

#define PaTrack_RICHDATASIZE 21

// eID: Return LS bit if eID, NLS bit if eVeto, NNLS bit if crude eID
unsigned int RICH::eID(const PaTrack &trk)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-eThr;
    if (dp>1.e-3) {  // Positivity + a margin to ensure acos's arg < 1 infra
      double idcut = dthCut[0], pvcut =  dthVeto[0];
      double thRich = trk.RichInf(RICH_ESTIMATOR);
      double phDensity = trk.RichInf(9)/thRich;
      double p2 = p*p;
      if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
	if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
	else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
      }
      dp = p-piThr;
      if (dp<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dp*dp;
	else              pvcut += dthVeto[3]*dp*dp;
      }
      double thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (pvcut<thRich-thMom)       pID |= 0x2;          //       ?Veto: 0x2
      thMom = acos(sqrt(M2_e/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                      // Crude eID:    0x4
	if (fabs(thRich-thMom)<idcut) pID |= 0x1;        //       eID:    0x1
      }
    }
  }
  return pID;
}
// muid: Return LS bit if muid, NLS bit if eVeto, NNLS bit if crude muid
unsigned int RICH::muID(const PaTrack &trk)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-muThr;
    if (dp>1.e-3) {  // Positivity + a margin to ensure acos's arg < 1 infra
      double idcut = dthCut[0], pvcut =  dthVeto[0];
      double thRich = trk.RichInf(RICH_ESTIMATOR);
      double phDensity = trk.RichInf(9)/thRich;
      double p2 = p*p;
      if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
	if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
	else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
      }
      dp = p-piThr;
      if (dp<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dp*dp;
	else              pvcut += dthVeto[3]*dp*dp;
      }
      double thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (pvcut<thRich-thMom)       pID |= 0x2;          //       ?Veto: 0x2
      thMom = acos(sqrt(M2_mu/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                      // Crude muid:    0x4
	if (fabs(thRich-thMom)<idcut) pID |= 0x1;        //       muid:    0x1
      }
    }
  }
  return pID;
}
// KID: Return LS bit if KID and NLS bit if piVeto (N.B. other bits also set)
unsigned int RICH::KID(const PaTrack &trk)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-KThr;
    if (dp>1.e-3) {  // Positivity + a margin to ensure acos's arg < 1 infra
      double idcut = dthCut[0], pvcut =  dthVeto[0];
      double thRich = trk.RichInf(RICH_ESTIMATOR);
      double phDensity = trk.RichInf(9)/thRich;
      double p2 = p*p;
      if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
	if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
	else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
      }
      dp = p-piThr;
      if (dp<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dp*dp;
	else              pvcut += dthVeto[3]*dp*dp;
      }
      double thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (pvcut/1.5<thMom-thRich) {
	pID |= 0x8;                                         // Crude piVeto: 0x8
	if (pvcut<thMom-thRich) pID |= 0x2;                 //       piVeto: 0x2
      }
      thMom = acos(sqrt(M2_K/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude KID:    0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       KID:    0x1
      }
    }
  }
  return pID;
}
// pID: Return LS bit if pID and NLS bit if KVeto (N.B. other bits also set)
unsigned int RICH::pID(const PaTrack &trk)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-pThr;
    if (dp>1.e-3) {  // Positivity + a margin to ensure acos's arg < 1 infra
      double idcut = dthCut[0], pvcut =  dthVeto[0];
      double thRich = trk.RichInf(RICH_ESTIMATOR);
      double phDensity = trk.RichInf(9)/thRich;
      double p2 = p*p;
      if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
	if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
	else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
      }
      dp = p-KThr;
      if (dp<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dp*dp;
	else              pvcut += dthVeto[3]*dp*dp;
      }
      double thMom = acos(sqrt(M2_K/p2+1)/CurrentIndAPV)*1000;  // Use K for Veto
      if (pvcut/1.5<thMom-thRich) {
	pID |= 0x8;                                         // Crude KVeto:  0x8
	if (pvcut<thMom-thRich) pID |= 0x2;                 //       KVeto:  0x2
      }
      thMom = acos(sqrt(M2_p/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude pID:    0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       pID:    0x1
      }
    }
  }
  return pID;
}
// piID: Return LS bit if piID and NLS bit if KVeto (N.B. other bits also set)
unsigned int RICH::piID(const PaTrack &trk)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-piThr;
    if (dp>1.e-3) {  // Positivity + a margin to ensure acos's arg < 1 infra
      double idcut = dthCut[0], pvcut =  dthVeto[0];
      double thRich = trk.RichInf(RICH_ESTIMATOR);
      double phDensity = trk.RichInf(9)/thRich;
      double p2 = p*p;
      if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
	if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
	else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
      }
      dp = p-KThr;
      if (dp<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dp*dp;
	else              pvcut += dthVeto[3]*dp*dp;
      }
      double arg = sqrt(M2_K/p2+1)/CurrentIndAPV;               // Use K for Veto
      if (arg>=1) pID |= 0xa;
      else {
	double thMom = acos(arg)*1000;
	if (pvcut/1.5<thRich-thMom) {
	  pID |= 0x8;                                       // Crude KVeto:  0x8
	  if (pvcut<thRich-thMom) pID |= 0x2;               //       KVeto:  0x2
	}
      }
      double thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude piID:   0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       piID:   0x1
      }
    }
  }
  return pID;
}
unsigned int RICH::eID(const PaTrack &trk, double &thRich, double &thMom)
{
  unsigned int pID = 0; thRich = -1; thMom = -1;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-eThr;
    double idcut = dthCut[0], pvcut =  dthVeto[0];
    thRich = trk.RichInf(RICH_ESTIMATOR);
    double phDensity = trk.RichInf(9)/thRich;
    double p2 = p*p;
    if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
      if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
      else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
    }
    double dppi = p-piThr;
    if (dppi>1.e-3) {
      if (dppi<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dppi*dppi;
	else              pvcut += dthVeto[3]*dppi*dppi;
      }
      thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (pvcut<thRich-thMom)       pID |= 0x2;          //       ?Veto: 0x2
    }
    if (dp>1.e-3) {      // Positivity + a margin to ensure acos's arg < 1 infra
      thMom = acos(sqrt(M2_e/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                      // Crude eID:    0x4
	if (fabs(thRich-thMom)<idcut) pID |= 0x1;        //       eID:    0x1
      }
    }
  }
  return pID;
}
unsigned int RICH::KID(const PaTrack &trk, double &thRich, double &thMom)
{
  unsigned int pID = 0; thRich = -1; thMom = -1;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-KThr;
    double idcut = dthCut[0], pvcut =  dthVeto[0];
    thRich = trk.RichInf(RICH_ESTIMATOR);
    double phDensity = trk.RichInf(9)/thRich;
    double p2 = p*p;
    if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
      if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
      else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
    }
    double dppi = p-piThr;
    if (dppi>1.e-3) {
      if (dppi<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dppi*dppi;
	else              pvcut += dthVeto[3]*dppi*dppi;
      }
      double thMompi = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (pvcut/1.5<thMompi-thRich) {
	pID |= 0x8;                                         // Crude piVeto: 0x8
	if (pvcut<thMompi-thRich) pID |= 0x2;               //       piVeto: 0x2
      }
    }
    if (dp>1.e-3) {      // Positivity + a margin to ensure acos's arg < 1 infra
      thMom = acos(sqrt(M2_K/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude KID:    0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       KID:    0x1
      }
    }
  }
  return pID;
}
unsigned int RICH::pID(const PaTrack &trk, double &thRich, double &thMom)
{
  unsigned int pID = 0;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-pThr;
    double idcut = dthCut[0], pvcut =  dthVeto[0];
    thRich = trk.RichInf(RICH_ESTIMATOR);
    double phDensity = trk.RichInf(9)/thRich;
    double p2 = p*p;
    if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
      if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
      else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
    }
    double dpK = p-KThr;
    if (dpK>1.e-3) {
      if (dpK<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dpK*dpK;
	else              pvcut += dthVeto[3]*dpK*dpK;
      }
      double thMomK = acos(sqrt(M2_K/p2+1)/CurrentIndAPV)*1000; // Use K for Veto
      if (pvcut/1.5<thMomK-thRich) {
	pID |= 0x8;                                         // Crude KVeto:  0x8
	if (pvcut<thMomK-thRich) pID |= 0x2;                //       KVeto:  0x2
      }
    }
    if (dp>1.e-3) {      // Positivity + a margin to ensure acos's arg < 1 infra
      thMom = acos(sqrt(M2_p/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude pID:    0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       pID:    0x1
      }
    }
  }
  return pID;
}
unsigned int RICH::piID(const PaTrack &trk, double &thRich, double &thMom)
{
  unsigned int pID = 0; thRich = -1; thMom = -1;
  if (trk.NRichInf()==PaTrack_RICHDATASIZE) {
    double p = (trk.vTPar()[0]).Mom(), dp = p-piThr;
    double idcut = dthCut[0], pvcut =  dthVeto[0];
    thRich = trk.RichInf(RICH_ESTIMATOR);
    double phDensity = trk.RichInf(9)/thRich;
    double p2 = p*p;
    if (dp<dthCut[1]) {            // Quadratic dependence upon p below [1]...
      if (phDensity<.0) idcut += dthCut[2]*dp*dp;   // ...depending in turn...
      else              idcut += dthCut[3]*dp*dp;   // ...upon photon density
    }
    double dpK = p-KThr;
    if (dpK>1.e-3) {
      if (dpK<dthVeto[1]) {
	if (phDensity<.0) pvcut += dthVeto[2]*dpK*dpK;
	else              pvcut += dthVeto[3]*dpK*dpK;
      }
      thMom = acos(sqrt(M2_K/p2+1))*1000;                    // Use K for Veto
      if (pvcut/1.5<thRich-thMom) {
	pID |= 0x8;                                       // Crude KVeto:  0x8
	if (pvcut<thRich-thMom) pID |= 0x2;               //       KVeto:  0x2
      }
    }
    else pID |= 0xa;
    if (dp>1.e-3) {      // Positivity + a margin to ensure acos's arg < 1 infra
      thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
      if (fabs(thRich-thMom)<idcut*1.5) {
	pID |= 0x4;                                         // Crude piID:   0x4
	if (fabs(thRich-thMom)<idcut)  pID |= 0x1;          //       piID:   0x1
      }
    }
  }
  return pID;
}

// piKID: Return bit pattern: LS: ID exists
//                           NLS: p>piThr
//                          NNLS: p>KThr...
//       ...AND diffs (measure-expectation)/resol
unsigned int RICH::piKID(const PaTrack &trk, double &piDth, double &KDth)
{
  unsigned int pID = trk.NRichInf()==PaTrack_RICHDATASIZE ? 1 : 0;
  double p = (trk.vTPar()[0]).Mom(), dp = p-piThr;
  if (dp>1.e-3)   // Positivity + a margin to ensure acos's arg < 1 infra
    pID |= 0x2;
  if (pID==0x3) {
    double picut = dthvsP[0], Kcut =  dthvsP[0];
    double thRich = trk.RichInf(RICH_ESTIMATOR);
    double phDensity = trk.RichInf(9)/thRich;
    double p2 = p*p;
    if (dp<dthvsP[1]) {            // Quadratic dependence upon p below [1]...
      if (phDensity<.0) picut += dthvsP[2]*dp*dp;   // ...depending in turn...
      else              picut += dthvsP[3]*dp*dp;   // ...upon photon density
    }
    dp = p-KThr;
    if (dp<dthvsP[1]) {
      if (phDensity<.0) Kcut += dthvsP[2]*dp*dp;
      else              Kcut += dthvsP[3]*dp*dp;
    }
    double arg = sqrt(M2_K/p2+1)/CurrentIndAPV;               // Use KID
    if (arg<1) {	
      double thMom = acos(arg)*1000;
      pID |= 0x4; KDth = (thRich-thMom)/Kcut;
    }
    else KDth = 1.e6;
    double thMom = acos(sqrt(M2_pi/p2+1)/CurrentIndAPV)*1000;
    piDth = (thRich-thMom)/picut;
  }
  else {
    piDth = 1.e6; KDth = 2.e6;
    if (p>KThr) pID |= 0x4;
  }
  return pID;
}

// *************** CORRECT chi2's FROM ERROR ON INDEX ***************

bool RICH::CorrChiSq(double& chisqPi, double& chisqKa, double& chisqPr,
		     double  chisqRi, double theRi, double mom)
{

//   From Paolo   03/06/10

  static float massPi = 0.13957;
  static float massKa = 0.49369;
  static float massPr = 0.93827;

  bool condChi = false;
  if( chisqPi > 0.  &&  chisqKa > 0.  &&  chisqRi > 0.  &&
      chisqPi < 1000.  &&  chisqKa < 1000.  &&  chisqRi < 1000. )
    condChi = true;
  if( !condChi ) return  condChi;

  mom = fabs( mom );

  double ePi = sqrt(mom*mom + massPi*massPi);
  double tThePi = ePi / (CurrentIndAPV * mom);
  if( tThePi > 1. ) tThePi = 1.;
  tThePi = acos( tThePi ) * 1000.;
  double pThePi = ePi / (CurrentProdIndAPV * mom);
  if( pThePi > 1. ) pThePi = 1.;
  pThePi = acos( pThePi ) * 1000.;
  double eKa = sqrt(mom*mom + massKa*massKa);
  double tTheKa = eKa / (CurrentIndAPV * mom);
  if( tTheKa > 1. ) tTheKa = 1.;
  tTheKa = acos( tTheKa ) * 1000.;
  double pTheKa = eKa / (CurrentProdIndAPV * mom);
  if( pTheKa > 1. ) pTheKa = 1.;
  pTheKa = acos( pTheKa ) * 1000.;
  double ePr = sqrt(mom*mom + massPr*massPr);
  double tThePr = ePr / (CurrentIndAPV * mom);
  if( tThePr > 1. ) tThePr = 1.;
  tThePr = acos( tThePr ) * 1000.;
  double pThePr = ePr / (CurrentProdIndAPV * mom);
  if( pThePr > 1. ) pThePr = 1.;
  pThePr = acos( pThePr ) * 1000.;

  double dChiKaPi = chisqKa - chisqPi;
  double dChiRiKa = chisqRi - chisqKa;
  double dTheKaRi = pTheKa - theRi;
  double dThePiKa = pThePi - pTheKa;
  double dTheRiKa2 = theRi*theRi - pTheKa*pTheKa;
  double dTheKaPi2 = pTheKa*pTheKa - pThePi*pThePi;

  double sumoS2 = (dChiKaPi * dTheKaRi - dChiRiKa * dThePiKa) /
                  (dTheKaPi2 * dTheKaRi - dTheRiKa2 * dThePiKa);
  double sumThoS2 = (dChiKaPi - dTheKaPi2* sumoS2) / (2.* dThePiKa);
  double sumTh2oS2 = chisqPi + 2.* pThePi* sumThoS2 - pThePi*pThePi* sumoS2;

  double chisqPiC = sumTh2oS2 - 2.* tThePi* sumThoS2 + tThePi*tThePi* sumoS2;
  double chisqKaC = sumTh2oS2 - 2.* tTheKa* sumThoS2 + tTheKa*tTheKa* sumoS2;
  double chisqPrC = sumTh2oS2 - 2.* tThePr* sumThoS2 + tThePr*tThePr* sumoS2;

  chisqPi = chisqPiC;
  chisqKa = chisqKaC;
  chisqPr = chisqPrC;

  return  condChi;

}

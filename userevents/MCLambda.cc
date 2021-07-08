// $Id: MCLambda.cc,v 1.3 2020/07/30 10:55:31 ybedfer Exp $

// $Log: MCLambda.cc,v $
// Revision 1.3  2020/07/30 10:55:31  ybedfer
// Compilation warning fix.
//
// Revision 1.2  2012/11/05 19:59:30  ybedfer
// No change w.r.t. v1.2.
//
// Revision 1.1  2009/07/03 14:05:11  ybedfer
// Initial revision
//

#include "MCLambda.h"

MCLambda* MCLambda::address = 0; // Init static pointer

MCLambda::MCLambda(double yLow, double yUp, double pTCut)
{
  if (address) {
  printf("** MCLambda::MCLambda:\a MCLambda already instantiated\n"); exit(1);
  }
  address = this;

  if ((MC = MCInfo::Ptr())==0) {
  printf("** MCLambda::MCLambda:\a No MCInfo yet instantiated!\n"); exit(1);
  }

  yLowCut = yLow; yUpCut = yUp;                         // ***** INITIALIZATIONS
  V0pTCut = pTCut; V0cthCut = 0;

  // *************** HISTOS BOOKING **********

  char title[] =
    "Reco'ible/Reco'd #bar{#Lambda}* vs. c#theta* - pT>25MeV,VV.V0>0.99990,RICH      ";
  //"Reco'ible/Reco'd #bar{#Lambda}* vs. c#theta* (Q2>1,.35<y<.85)      ";
  for (int iLaL = 0; iLaL<2; iLaL++) {  // Loop on Lambda/anti-Lambda
    // Successively all Ls
    //              reco'ible
    //              reco'ible RICH'able
    char hName[] = "hc_ArALki"; char LaL[] = "#bar{#Lambda}", LAL[] = "AL";
    if (iLaL) { sprintf(LaL,"bar{#Lambda}"); sprintf(LAL,"AL"); }
    else      { sprintf(LaL,"#Lambda");      sprintf(LAL,"L"); }
    const char *vars[] = {"c#theta*","Pp"}, cps[] = "cp"; 
    for (int icp = 0; icp<2; icp++) {
      const char *var = vars[icp], cp = cps[icp];
      int nBins; double xMn, xMx;
      if (icp==0) { nBins = 16; xMn = -1; xMx = 1; }
      else        { nBins = 40; xMn =  0; xMx = 80; }
      // "(Q2,y)-Accepted/Reconstructed" histos: hc_Ar..
      sprintf(hName,"h%c_Ar%s",cp,LAL);
      sprintf(title,"(Q2,y)OK/Reco'd %s vs. %s (Q2>1,%.2f<y<%.2f)",
	      LaL,var,yLow,yUp);
      hc_ArL[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Ar%sk",cp,LAL);
      sprintf(title,"(Q2,y)OK/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_ArLk[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Ar%si",cp,LAL);
      sprintf(title,"(Q2,y)OK/Reco'd %s vs. %s - RICH",LaL,var);
      hc_ArLi[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Ar%ski",cp,LAL);
      sprintf(title,"(Q2,y)OK/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f,RICH",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_ArLki[icp][iLaL] = new TH1D(hName,title,nBins,xMn,xMx);
      // "Reconstructed outside (Q2,y) acceptance" histos: hc_r..
      sprintf(hName,"h%c_r%s",cp,LAL);
      sprintf(title,"(Q2,y)!OK/Reco'd %s vs. %s (Q2>1,%.2f<y<%.2f)",
	      LaL,var,yLow,yUp);
      hc_rL[icp][iLaL]    = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_r%sk",cp,LAL);
      sprintf(title,"(Q2,y)!OK/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_rLk[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_r%si",cp,LAL);
      sprintf(title,"(Q2,y)!OK/Reco'd %s vs. %s - RICH",LaL,var);
      hc_rLi[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_r%ski",cp,LAL);
      sprintf(title,"(Q2,y)!OK/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f,RICH",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_rLki[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      // "Reconstructible and reconstructed" histos: hc_Rr
      sprintf(hName,"h%c_Rr%s",cp,LAL);
      sprintf(title,"Reco'ible/Reco'd %s vs. %s (Q2>1,%.2f<y<%.2f)",
	      LaL,var,yLow,yUp);
      hc_RrL[icp][iLaL]   = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Rr%sk",cp,LAL);
      sprintf(title,"Reco'ible/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_RrLk[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Rr%si",cp,LAL);
      sprintf(title,"Reco'ible/Reco'd %s vs. %s - RICH",LaL,var);
      hc_RrLi[icp][iLaL]  = new TH1D(hName,title,nBins,xMn,xMx);
      sprintf(hName,"h%c_Rr%ski",cp,LAL);
      sprintf(title,"Reco'ible/Reco'd %s vs. %s - pt>%.0fMeV,V0.VV>%.5f,RICH",
	      LaL,var,V0pTCut*1000,V0cthCut);
      hc_RrLki[icp][iLaL] = new TH1D(hName,title,nBins,xMn,xMx);
    }

    // ***** Sigma* -> Lpi: ACCEPTANCE, RECONSTRUCTIBILITY *****
    // Successively all Sigma*
    //              Reco'ible w/in kin cuts and RICH'able
    if (iLaL) { sprintf(LAL,"AS"); }
    else      { sprintf(LAL,"S"); }
    // "(Q2,y)-Accepted/Reconstructed" histo: hc_Ar..
    sprintf(hName,"hc_Ar%s",LAL);
    sprintf(title,"(Q2,y)OK/Reco'd %s* vs. c#theta* (Q2>1,%.2f<y<%.2f)",
	    LaL,yLow,yUp);
    hc_ArS[iLaL]   = new TH1D(hName,title,16,-1,1);
    // "(Q2,y)-Accepted/Reconstructed" histo: hc_r..
    sprintf(hName,"hc_r%s",LAL);
    sprintf(title,"!(Q2,y)OK/Reco'd %s* vs. c#theta* (Q2>1,%.2f<y<%.2f)",
	    LaL,yLow,yUp);
    hc_rS[iLaL]    = new TH1D(hName,title,16,-1,1);
    // "Reconstructible and reconstructed" histos: hc_Rr
    sprintf(hName,"hc_Rr%ski",LAL);
    sprintf(title,"Reco'ible %s* vs. c#theta* - pt>%.0fMeV,V0.VV>%.5f,RICH",
	    LaL,V0pTCut*1000,V0cthCut);
    hc_RrSki[iLaL] = new TH1D(hName,title,16,-1,1);
  }

}

void MCLambda::Update()  // Import the MCInfo flags of the current PaEvent's MC tracks and check (and store) MC Q2 and y are w/in range
{
  MCTypes = MC->GetTypes();
  MCDIS = 0;
  if (MC->Q2>1)                        MCDIS |= 0x1;
  if (yLowCut<MC->yB && MC->yB<yUpCut) MCDIS |= 0x2;
}

int MCLambda::Fill(PaEvent &e, // Only used to document the error messages
		   PaTrack &tp, PaTrack &tpi, // p and pi reco'd tracks
		   int iLaL,      // Lambda or anti-Lambda
		   bool winMCut)  // Reco passed mass cut
  // Retuned value = bit pattern: LS 2 bits = 2^iLaL, if argument (p,pi) corresponds to a genuine MC Lambda
{
  int iMCTp  = tp.iMCtrack(); if (iMCTp<0) return 0;
  int typep  = MCTypes[iMCTp];
  if (!(typep&0x01000000))  return 0; // p not from (a)Lambda
  int iMCTpi = tpi.iMCtrack(); if (iMCTpi<0) return 0;
  int typepi = MCTypes[iMCTpi];
  if (!(typepi&0x02000000)) return 0; // pi not from (a)Lambda

  // X-check that p and pi correspond to Lambda w/ C = "iLaL"
  if (iMCTp!=MC->iMCL_p[iLaL] || iMCTpi!=MC->iMCL_pi[iLaL]) return 0;
  // Consistency check: Lambda w/ C = "iLaL" does exist
  int iMCTL = MC->iMCT_L[iLaL]; if (iMCTL<0) {
 printf("** MCLambda::MCLambda:\a Run %d Evt %d: Inconsistency: Lambda(C=%d) -> p+pi(=%d+%d) does not exist!\n",
	e.RunNum(),(int)e.UniqueEvNum(),iLaL,iMCTp,iMCTpi); exit(1);
  }

  // Input (p,pi) does correspond to MC w/ C = "iLaL"

  if (winMCut) {
    int typeL = MCTypes[iMCTL];
    double cthstar = MC->cthstarL[iLaL], pp = MC->PL_p[iLaL];
    if ((MCDIS&0x3)==0x3) {
      hc_ArL[0][iLaL]->Fill(cthstar); hc_ArL[1][iLaL]->Fill(pp);
      if (typeL&0x00400000) {
	hc_ArLk[0][iLaL]->Fill(cthstar); hc_ArLk[1][iLaL]->Fill(pp);
      }
      if (typep&0x10000000) {
	hc_ArLi[0][iLaL]->Fill(cthstar); hc_ArLi[1][iLaL]->Fill(pp);
	if (typeL&0x00400000) {
	  hc_ArLki[0][iLaL]->Fill(cthstar); hc_ArLki[1][iLaL]->Fill(pp);
	}
      }
    }
    else {
      hc_rL[0][iLaL]->Fill(cthstar);  hc_rL[1][iLaL]->Fill(pp);
      if (typeL&0x00400000) {
	hc_rLk[0][iLaL]->Fill(cthstar); hc_rLk[1][iLaL]->Fill(pp);
      }
      if (typep&0x10000000) {
	hc_rLi[0][iLaL]->Fill(cthstar); hc_rLi[1][iLaL]->Fill(pp);
	if (typeL&0x00400000) {
	  hc_rLki[0][iLaL]->Fill(cthstar); hc_rLki[1][iLaL]->Fill(pp);
	}
      }
    }

    if (typeL&0x0080) {
      // Reconstruction efficiency = Rr / R (from MCInfo)
      // Note: In order to make it here the reco'd (p,pi) will have had to pass
      // the RICH ID, while it would be interesting to know the reco efficiency
      // independent of the RICH efficiency. This one can achieve by replaying
      // the MC data w/ 100% efficient RICH.
      hc_RrL[0][iLaL]->Fill(cthstar); hc_RrL[1][iLaL]->Fill(pp);
      if (typeL&0x00400000) {
	hc_RrLk[0][iLaL]->Fill(cthstar); hc_RrLk[1][iLaL]->Fill(pp);
      }
      if (typep&0x10000000) {
	hc_RrLi[0][iLaL]->Fill(cthstar); hc_RrLi[1][iLaL]->Fill(pp);
	if (typeL&0x00400000) {
	  hc_RrLki[0][iLaL]->Fill(cthstar); hc_RrLki[1][iLaL]->Fill(pp);
	}
      }
    }
  }
  return (iLaL? 0x2 : 0x1);
}

// $Id: HistoLambda.cc,v 1.3 2013/08/08 21:19:43 ybedfer Exp $

// Booking for Lambda selection

// $Log: HistoLambda.cc,v $
// Revision 1.3  2013/08/08 21:19:43  ybedfer
// Bypass the filling "hm_LcllX". Because it aborts the execution. Haven't understood why...
//
// Revision 1.2  2012/11/05 20:23:57  ybedfer
// - Histos now created in sub-TDirectory "Lambda" of base directory.
// - More systematics in the nomenclature of histo names and titles.
// - "FillLambdaMC": hm_[a]Lambda now 2D, w/ Y-binning distinguishing among
//  fake or genuine, pi mu-ID'd or not. .
//
// Revision 1.1  2009/07/03 14:31:15  ybedfer
// Initial revision
//

#include "TH1.h"
#include "TH2.h"
#include "Masses.h"
#include "Phast.h"

// *************** Lambda MASS conditioned by VERTEX LINE ***************
static TH1D *hm_Lc0, *hm_aLc0;
static TH1D *hm_Lc11, *hm_Lc12, *hm_Lc21, *hm_Lc22, *hm_aLc22;
static TH1D *hm_LC22;                      // Requiring mu/mu'
static TH1D *hm_Lcll, *hm_Lcsl, *hm_Lcss;  // In LAS, L+SAS, SAS
static TH1D *hm_Lcnn;                      // No LAS: tracking downstream of SM1
static TH1D *hm_LcssS, *hm_LcssL;          // In SAS, SAT and LAT
static TH2D *hm_LcllX, *hm_LcssX;          // Vs. Theta X
static TH2D *hm_LcllB, *hm_LcssB;          // Vs. Beta X
static TH1D *hm_Li, *hm_Lo;
// ***** High level histos (to be enabled explicitly in the calling routine) 
static TH1D *hm_LV22;                      // Either decay belongs to best pV
static TH1D *hm_LR32;                      // p outside RICH, w/ tight V0 cut
static TH1D *hm_LK32;                      // hm_LR32, w/ pipi outside K0
static TH2D *hm_Lambda, *hm_aLambda;       // Genuine/Fake MC (a-)Lambda
static int hLevel_Lambda;

void BookLambda(const double *V0DsdD, const double *V0cthCut, double dM_Lambda,
		double pTCut,
		int PIDScheme,  // 0x1=p: pID|e,pi,KVeto, 0x2=pi: e-rejection
		int histoLevel)
{
  gDirectory->cd("/");
  TDirectory *dLambda; if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda")))
    gDirectory->mkdir("Lambda","Lambda Histos");
  if (!(dLambda = (TDirectory*)gDirectory->Get("Lambda"))) {
    printf("\n** US3:\a No creating subdir \"Lambda\" in base TDirectory \"%s\"\n\n",
	   gDirectory->GetName()); assert(false);
  }
  dLambda->cd();

  int nBinsX = int((M_Lam+dM_Lambda-M_p-M_pi+.005)/.0005+.5);
  double xMn = M_Lam+dM_Lambda-nBinsX*.0005;
  char title[] =
    "p^{+}#pi^{-} - D/#deltaD>2.0#sigma,cos#theta>.99990,p_{T}<100MeV,p:pID,e#piKVeto,#pi:eVeto - (L)SAS^{2} vs. #Theta_{X}   ";
  //"p^{+}#pi^{-} - D/#deltaD>2.0#sigma,cos#theta>.99990,p_{T}<100MeV,p:pID,e#piKVeto,#pi:eVeto - !w/in RICH|K0   ";
  //"p^{+}#pi^{-} - D/#deltaD>2.0#sigma,cos#theta>.99990,p_{T}<100MeV,p:pID,e#piKVeto,#pi:eVeto - MC #bar{#Lambda}  ";
  char LFormat[]  =
    "p#pi^{-} - D/#deltaD>%.1f#sigma,cos#theta>%.5f,p_{T}>%.0fMeV%s%s%s";
  char aLFormat[] =
    "p^{-}#pi^{+} - D/#deltaD>%.1f#sigma,cos#theta>%.5f,p_{T}>%.0fMeV%s%s%s";

  char pScheme[] = ",p:pID,e#piKVeto", piScheme[] = ",#pi:eVeto";
  const char *pID  = (PIDScheme&0x1) ? pScheme  : "\0";
  const char *piID = (PIDScheme&0x2) ? piScheme : "\0";

  pTCut *= 1000;
  sprintf(title,LFormat,V0DsdD[0],V0cthCut[0],pTCut,pID,piID," ");
  hm_Lc11  = new TH1D("hm_Lc11",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[0],V0cthCut[1],pTCut,pID,piID," ");
  hm_Lc12  = new TH1D("hm_Lc12",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[0],pTCut,pID,piID," ");
  hm_Lc21  = new TH1D("hm_Lc21",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," ");
  hm_Lc22  = new TH1D("hm_Lc22",title,nBinsX,xMn,M_Lam+dM_Lambda);

  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - #mu#mu'");
  hm_LC22  = new TH1D("hm_LC22",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - OUT");
  hm_Lo    = new TH1D("hm_Lo",  title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - IN");
  hm_Li    = new TH1D("hm_Li",  title,nBinsX,xMn,M_Lam+dM_Lambda);

  //                                  Lambda as a function of spectrometer
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - (L)SAS^{2}");
  hm_Lcss  = new TH1D("hm_Lcss",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - LAS#times(L)SAS");
  hm_Lcsl  = new TH1D("hm_Lcsl",title,nBinsX,xMn,M_Lam+dM_Lambda);
  sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - LAS^{2}");
  hm_Lcll  = new TH1D("hm_Lcll",title,nBinsX,xMn,M_Lam+dM_Lambda);

  hLevel_Lambda = histoLevel;
  if (histoLevel&0x10) {
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - SAS^{2}");
    hm_Lcnn  = new TH1D("hm_Lcnn",title,nBinsX,xMn,M_Lam+dM_Lambda);
    //                                  Vs. Theta X
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - LAS^{2} vs. #Theta_{X}");
    hm_LcllX  = new TH2D("hm_LcllX",title,nBinsX,xMn,M_Lam+dM_Lambda,4,-.5,3.5);
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - (L)SAS^{2} vs. #Theta_{X}");
    hm_LcssX  = new TH2D("hm_LcssX",title,nBinsX,xMn,M_Lam+dM_Lambda,4,-.5,3.5);
    //                                  Vs. Beta X
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - LAS^{2} vs. #Beta_{X}");
    hm_LcllB  = new TH2D("hm_LcllB",title,nBinsX,xMn,M_Lam+dM_Lambda,4,-.5,3.5);
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - (L)SAS^{2} vs. #Beta_{X}");
    hm_LcssB  = new TH2D("hm_LcssB",title,nBinsX,xMn,M_Lam+dM_Lambda,4,-.5,3.5);
  }

  if (histoLevel&0x10) {
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - SAS/SAT");
    hm_LcssS  = new TH1D("hm_LcssS",title,nBinsX,xMn,M_Lam+dM_Lambda);
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - SAS/LAT");
    hm_LcssL  = new TH1D("hm_LcssL",title,nBinsX,xMn,M_Lam+dM_Lambda);
  }
  else {
    hm_LcssS=hm_LcssL = 0;
  }

  //                                 \bar{Lambda}
  sprintf(title,aLFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," ");
  hm_aLc22 = new TH1D("hm_aLc22",title,nBinsX,xMn,M_Lam+dM_Lambda);

  if (histoLevel&0x100) {         // Monte-Carlo genuine Lambda's
    // MC would have to be 3D, w/ kin and cos(theta*) binning
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID," - MC #Lambda");
    hm_Lambda  = new TH2D("hm_Lambda", title,nBinsX,xMn,M_Lam+dM_Lambda,
			  5,-.5,4.5);
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - MC #bar{#Lambda}");
    hm_aLambda = new TH2D("hm_aLambda",title,nBinsX,xMn,M_Lam+dM_Lambda,
			  4,-.5,3.5);
  }

  // Special histos used for debugging
  if (histoLevel&0x2) {
    sprintf(title,LFormat,V0DsdD[1],V0cthCut[1],pTCut,pID,piID,
	    " - p|#pi in pV");
    hm_LV22  = new TH1D("hm_LV22",title,nBinsX,xMn,M_Lam+dM_Lambda);
  }
  if (histoLevel&0x4) {
    sprintf(title,LFormat,V0DsdD[2],V0cthCut[1],pTCut,"\0","\0",
	    " - !w/in RICH");
    hm_LR32  = new TH1D("hm_LR32",title,nBinsX,xMn,M_Lam+dM_Lambda);
  }  
  if (histoLevel&0x8) {
    sprintf(title,LFormat,V0DsdD[2],V0cthCut[1],pTCut,"\0","\0",
	    " - !w/in RICH|K0");
    hm_LK32  = new TH1D("hm_LRK2",title,nBinsX,xMn,M_Lam+dM_Lambda);
  }  
}

void FillLambda(double m_ppi, double m_pip,
		int s_dist, int s_angle,  // Vertex line cut indices
		bool hasmumup,            // p-V has mu/mu'
		const PaTrack &trk1,      // >0 track
		const PaTrack &trk2,      // <0 track
		bool downstreamOfTarget,  // s-Vertex downstream of target
		int pID)      // 2-bit pattern: LSB = Lambda, MSB = \bar{Lambda}
{
  // The routine is expecting that a cut on pT has already been performed.
  // Further cuts on vertex line, viz.:...

  if (s_dist>=1 &&          // ...Vertex ditance
      s_angle>=1) {         // ...Angle between vertex line and V0 momentum
    if (pID&0x1) {
      hm_Lc11->Fill(m_ppi);
      if (s_angle>=2) hm_Lc12->Fill(m_ppi);
      if (s_dist>=2) hm_Lc21->Fill(m_ppi);
    }
    if (s_dist>=2 && s_angle>=2) { // ********** STRICTEST CUTS:... **********
      if (pID&0x1) {
	hm_Lc22->Fill(m_ppi);
	if (hasmumup) hm_LC22->Fill(m_ppi);        // ...Requiring mu/mu'

	double zL1 = trk1.ZLast(), zL2 = trk2.ZLast();
	if (zL1<1950) {                            // ...LAS and SAS Lambdas
	  if (zL2<1950) {
	    hm_Lcll->Fill(m_ppi);
	    if (hLevel_Lambda&0x10) {
	      const PaTPar &hi1 = trk1.vTPar()[0], &hi2 = trk2.vTPar()[0];
	      double ThetX = hi2(3)-hi1(3), bin = floor(ThetX/.04)+2;
	      if (bin<0) bin = 0; if (bin>3) bin = 3;
	      hm_LcllX->Fill(m_ppi,bin);                // vs. Theta X
	      double BetaX = trk2.vAux()[0]-trk1.vAux()[0];
	      bin = floor(BetaX/.04)+3; if (bin<0) bin = 0; if (bin>5) bin = 5;
	      hm_LcllB->Fill(m_ppi,bin);                // vs. Beta X
	    }
	  }
	  else          hm_Lcsl->Fill(m_ppi);
	}
	else {
	  if (zL2<1950) hm_Lcsl->Fill(m_ppi);
	  else {
	    hm_Lcss->Fill(m_ppi);                  // SAS times SAS
	    if (hLevel_Lambda&0x10) {
	      double ThetX = trk2.vAux()[0]-trk1.vAux()[0];
	      double bin = floor(ThetX/.025)+2;
	      if (bin<0) bin = 0; if (bin>3) bin = 3;
	      hm_LcssX->Fill(m_ppi,bin);                // vs. Theta X
	      double BetaX = trk2.vAux()[2]-trk1.vAux()[2];
	      bin = floor(BetaX/.025)+3; if (bin<0) bin = 0; if (bin>5) bin = 5;
	      hm_LcssB->Fill(m_ppi,bin);                // vs. Beta X

	      //                        No LAS: tracking only downstream of SM1
	      const PaTrack *t; int pm, las;
	      for (pm=las = 0, t = &trk1; pm<2; pm++) {
		if (t->vTPar()[0](0)<350) las |= 1<<pm; t = &trk2;
	      }
	      if (!las) hm_Lcnn->Fill(m_ppi);
	    }
	  }
	}

	if (downstreamOfTarget) hm_Lo->Fill(m_ppi);// ...In/outside target
	else                    hm_Li->Fill(m_ppi);
      }
      if (pID&0x2) hm_aLc22->Fill(m_pip);          // ...Anti-Lambda's
    }
  }
}
// ***** FILL HIGH LEVEL HISTOS
void FillLambda(double m_ppi, int histoLevel)
{
  if (histoLevel&0x2) {
    if (!hm_LV22) {
      printf("** FillLambda:\a High level \"hm_LV22\" accessed, yet not booked!\n");
      exit(1);
    }
    hm_LV22->Fill(m_ppi);
  }
  if (histoLevel&0x4) {
    if (!hm_LR32) {
      printf("** FillLambda:\a High level \"hm_LR32\" accessed, yet not booked!\n");
      exit(1);
    }
    hm_LR32->Fill(m_ppi);
  }
  if (histoLevel&0x8) {
    if (!hm_LK32) {
      printf("** FillLambda:\a High level \"hm_LK32\" accessed, yet not booked!\n");
      exit(1);
    }
    hm_LK32->Fill(m_ppi);
  }
  if (histoLevel&0x10) {        // SAS w/ SAT at the exit of SM2?
    if (!hm_LcssS) {
      printf("** FillLambda:\a High level \"hm_LcssS\" accessed, yet not booked!\n");
      exit(1);
    }
    hm_LcssS->Fill(m_ppi);
  }
  if (histoLevel&0x20) {        // SAS w/ LAT at the exit of SM2?
    if (!hm_LcssL) {
      printf("** FillLambda:\a High level \"hm_LcssL\" accessed, yet not booked!\n");
      exit(1);
    }
    hm_LcssL->Fill(m_ppi);
  }
}
void FillLambdaMC(double m,
		  int pID,      // 2^iLaL
		  int genuine,  // LS 2 bits: 2^iLaL, NLS 2 bits: HH decay
		  bool kinMCOK,
		  bool piIsMu)
{
  int fake;
  if      (genuine&0x3) {
    if (kinMCOK)        fake = 0;
    else                fake = 1;
  }
  else if (genuine)     fake = 2;
  else if (piIsMu)      fake = 3;
  else                  fake = 4;
  if      (pID&0x1) hm_Lambda ->Fill(m,fake);    // MC Lambda
  else if (pID&0x2) hm_aLambda->Fill(m,fake);    // MC a-Lambda
}

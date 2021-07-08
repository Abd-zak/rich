// $Id$

// Plot CSEventData filtered by "fit_table plots"
// - Plot histos of the "a" (i.e. "all") series, such as "h_K0_pip_a_14_3"...
// - ...for given <particleName>, viz.: K0 or Lambda or phi,
// - ...for pip or pim, depending on <decayCharge> (kp or km in the "phi" case).
// - Looping on momentum and angle bins.

#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

bool parsePCSArgs(const char *particleName, int decayCharge, int armenteros,
		  int binning, int nAs,
		  int &iPar, int &iA0, int &iAF, char &cpm, int &Julien_3);
bool getPCSColours(int nAs, int *&cols, int &ncols);
void setPCSStyle();
bool getPCSHist2(const char *pN, const char *decay, int cpm, int iP, int iA,
		 TH2D *&h2);
bool getPCSHisto(const char *pN, const char *decay, int cpm, int iP, int iA,
		 TH1D *&h);
void getPCSCanvas(const char *pN, int cpm, int iP, int armenteros, TCanvas *&c);
void setPCSPaveTexts(TH1 *h, int rank, int optStat, int col, int mode,
		     int rankMx = 3);
bool PlotCSEvent(const char *particleName, int decayCharge,
		 int armenteros = -1 /* Armenteros' angle bin */,
		 int binning = 0)
{
  int iPar;                                   // Particle index
  int nPs = 15, iA0 = 0, nAs = 4, iAF = nAs;  // Momentum and Angle binning
  char cpm;                                   // Charge
  int Julien_3;    // Case of Julien_3: only 2 angle bins
  // ********** PARSE ARG.S
  if (!parsePCSArgs(particleName,decayCharge,armenteros,binning,nAs,
		    iPar,iA0,iAF,cpm,Julien_3)) return false;
  // ***** ASSOCIATED DECAY
  const char *decay[] = { "pi", "pi",     "k",   "k",    "pi"};
  // ***** LOCAL VARIABLES
  int iP, iA; const char *pN = particleName;
  static TH1D *h;
  int *cols, ncols; if (!getPCSColours(nAs,cols,ncols)) return false;
  setPCSStyle();  // ***** DRAWING STYLE
  // ***** LOOP ON P...
  // ***** IN TWO SUCCESSIVE CANVASES
  static TCanvas *c; static int iPad;
  static double yMx; TH1D *h0; // In each pad, remember maximum maximorum

  for (iP = 0; iP<nPs; iP++ ) {  // ********** LOOP ON P BINS
    iPad++; if (c) c->cd(iPad);
    for (iA = iA0, h0 = 0; iA<iAF; iA++) {  // ********** LOOP ON A BINS
      int col = Julien_3 ? cols[iA+1] : cols[iA];
      if (armenteros>=0) {       // ***** PLOT TH2D Armenteros
	TH2D *h2 = 0;
	if (!getPCSHist2(pN,decay[iPar],cpm,iP,iA,h2)) return false;
	if ((iP==0 || iP==8) && iA==iA0) {
	  getPCSCanvas(pN,cpm,iP,armenteros,c); iPad = 1; c->cd(iPad);
	}
	h2->GetXaxis()->SetNdivisions(505);
	h2->Draw("colz");
	gPad->Update(); setPCSPaveTexts(h2,iA,1000000,col,3);
      }
      else {                                 // ***** PLOT TH1D INVARIANT MASS
	if (!getPCSHisto(pN,decay[iPar],cpm,iP,iA,h)) return false;
	if ((iP==0 || iP==8) && iA==iA0) {
	  getPCSCanvas(pN,cpm,iP,armenteros,c); iPad = 1; c->cd(iPad);
	}
	h->SetLineColor(col);
	TAxis *ax;
	ax = h->GetXaxis(); ax->SetNdivisions(505); ax->SetLabelSize(.045);
	ax = h->GetYaxis(); ax->SetNdivisions(505); ax->SetLabelSize(.045);
	if (iA!=iA0) h->Draw("sames");
	else         h->Draw();
	// ***** GET MAXIMUM MAXIMORUM
	// ***** TITLE OF 0th HISTO -> OVERALL TITLE, W/ ALL theta RANGES
	double y = h->GetMaximum(); if (!h0) {
	  h0 = h; yMx = y;
	  const char *titre0 = h->GetTitle(), *cTheta;
	  if ((cTheta = strstr(titre0,"0.00<#theta<0.01"))) {
	    const char *AllA = "#theta: 0,#color[618]{.01},#color[601]{.04},#color[418]{.12},#color[801]{.3}";
	    size_t head = cTheta-titre0;
	    char *titreA = new char[head+strlen(AllA)+1];
	    sprintf(titreA,"%s",titre0);
	    sprintf(titreA+head,"%s",AllA);
	    h->SetTitle(titreA);
	  }
	}
	else {
	  int nBins = h->GetNbinsX(), bin; for (bin = 1; bin<=nBins; bin++) {
	    double y = h->GetBinContent(bin); if (y>yMx) yMx = y;
	  }
	}
	gPad->Update(); setPCSPaveTexts(h,iA,1000000,col,3);
      }
    }
    if (armenteros<0) h0->SetMaximum(yMx*1.05);
  }
  return true;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  parsePCSArgs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  getPCSHisto  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  getPCSHist2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ getPCSCanvas  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ getPCSColours ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bool parsePCSArgs(const char *particleName, int decayCharge, int armenteros,
		  int binning, int nAs,
		  int &iPar, int &iA0, int &iAF, char &cpm, int &Julien_3)
{
  if (armenteros<-1 || nAs<=armenteros) {
    printf("** PlotCSEvent: Bad arg. <armenteros> = %d. Must be in [-1,%d] => exit...\n",
	   armenteros,nAs-1);
    return false;
  }
  if (armenteros>=0) {  iA0 = armenteros; iAF = armenteros+1; }
  if (decayCharge!=1 && decayCharge!=-1) {
    printf("** PlotCSEvent: Bad arg. <decayCharge> = %d. Must be +1 or -1 => exit...\n",
	   decayCharge);
    return false;
  }
  int pm = (1-decayCharge)/2; cpm = pm?'m':'p';
  const char *pN = particleName; if (!pN) {
    printf("** PlotCSEvent: Unspecified arg. <particleName> => exit...\n");
    return false;
  }
  // ***** CHECK "particleName" is either of:
  const char *names[] = { "K0", "Lambda", "phi", "ephi", "k0" };
  int nPars = sizeof(names)/sizeof(char*);
  for (int iter = 0; iter<2; iter++) {
    int match; for (iPar = 0, match = -1; iPar<nPars; iPar++) {
      if (!strcmp(pN,names[iPar])) { match = iPar; break; }
      else if (iter) {
	if (!iPar) printf("** PlotCSEvent: Bad arg. <particleName> =\"%s\". Must be either of:",pN);
	printf(" \"%s\"",names[iPar]);
      }
    }
    if (match>=0) { iPar = match; break; }
    else if (iter) { printf("\n"); return false; }
  }
  if (armenteros<0 && iPar==4 || binning) { // Case of Julien_3: only 2 angles
    Julien_3 = 1; iA0 = 0; iAF = 2;
  }
  else
    Julien_3 = 0;
  return true;
}
bool getPCSColours(int nAs, int *&cols, int &ncols)
{
  // Create enough colours to cover the "nAs" angle bins
  int violet = kMagenta+2;
  int bleu   = kBlue+1;
  int vert   = kGreen+2;
  int orange = kOrange+1;
  static int colours[] = { violet, bleu, vert, orange };
  ncols = sizeof(colours)/sizeof(int); if (ncols!=nAs) {
    printf("** PlotCSEvent: Inconsistency # of cols(%d) != nAs(%d)\n",ncols,nAs);
    return false;
  }
  cols = colours;
  return true;
}
bool getPCSHist2(const char *pN, const char *decay, int cpm, int iP, int iA,
		 TH2D *&h)
{
  char hN[] = "am_Lambda_pip_pi_10_10";
  sprintf(hN,"am_%s_%s%c_a_%d_%d",pN,decay,cpm,iP,iA);
  if (!(h = (TH2D*)gDirectory->Get(hN))) {
    printf("** PlotCSEvent: TH2D \"%s\" does not exist => exit...\n",hN);
    return false;
  }
  return true;
}
bool getPCSHisto(const char *pN, const char *decay, int cpm, int iP, int iA,
		 TH1D *&h)
{
  char hN[] = "h_Lambda_pip_pi_10_10";
  sprintf(hN,"h_%s_%s%c_a_%d_%d",pN,decay,cpm,iP,iA);
  if (!(h = (TH1D*)gDirectory->Get(hN))) {
    printf("** PlotCSEvent: TH1D \"%s\" does not exist => exit...\n",hN);
    return false;
  }
  return true;
}
void getPCSCanvas(const char *pN, int cpm, int iP, int armenteros, TCanvas *&c)
{
  // Create a new TCanvas w/ as specific a name as possible, i.e.:
  // - based on particle name and decay charge,
  // - P bin (in fact first P bin in the sequence of 8 covered by the canvas,
  // - 'M' for mass histos, "armenteros" index otherwise.
  // - W/ appended 'p' if name already exists, so that one can draw w/in a
  //  same root session histos obtained in two different conditions.
  char cN[] = "cMLambdap00";
  if (armenteros<0) {
    sprintf(cN,"cM%s%c%d",            pN,cpm,iP<8?0:1);
    if (gROOT->FindObject(cN))
      sprintf(cN,"cM%s%c%dp",            pN,cpm,iP<8?0:1);
  }
  else {
    sprintf(cN,"c%d%s%c%d",armenteros,pN,cpm,iP<8?0:1);
    if (gROOT->FindObject(cN))
      sprintf(cN,"c%d%s%c%dp",armenteros,pN,cpm,iP<8?0:1);
  }
  int wtop = 64+(iP==8?256:0);
  c = new TCanvas(cN,cN,wtop,wtop,1024,512); c->Divide(4,2);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~     setPCSPaveTexts     ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~     setPCSStyle         ~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "TStyle.h"
void setPCSStyle()
{
  // Limit stats info to the minimum: not to clutter up the picture where we
  // already have many histos.
  gStyle->SetOptStat(1000000);
  // Setting a fix length format, and hence same for everybody, gives a more
  // uniform touch to the ensemble.
  // We expect counts in the ten thousands at most => "%7.0f" should be enough.
  gStyle->SetStatFormat("7.0f");
}
#include "TPaveStats.h"
#include "TClass.h"
#include "TH2D.h"
void setPCSPaveTexts(TH1 *h, int rank, int optStat, int col, int mode, int rankMx)
{
  // - "rank": 0 = rightmost, 1,2,.. shifted to the left
  // - "mode": 0 = TRatioPlot, 1 = standard horiz., 2 = horiz. w/ fit, 3 = vert.
  TPaveStats *st; if ((st = (TPaveStats*)h->GetFunction("stats"))) {
    if (optStat) st->SetOptStat(optStat);
    if (mode==3) {
      double dY; if (optStat==10) {
	if (h->IsA()->InheritsFrom(TH2::Class()))
	  dY = .08;
	else
	  dY = .16; // Assuming it's w/ fit
      }
      else
	dY = .12;
      if (rank<=rankMx) {
	st->SetY2NDC(.995-rank*dY); st->SetY1NDC(.995-(rank+1)*dY);
	st->SetX1NDC(.645); st->SetX2NDC(.995);
      }
      else {
	rank -= rankMx;
	st->SetY2NDC(.995-rank*dY); st->SetY1NDC(.995-(rank+1)*dY);
	st->SetX1NDC(.050); st->SetX2NDC(.350);
      }
    }
    else {
      st->SetX2NDC(.995-rank*.30); st->SetX1NDC(.695-rank*.30);
      st->SetY2NDC(.995);
      if      (mode==2)  st->SetY1NDC(.75);
      else if (mode==1) {
	if (optStat==10) st->SetY1NDC(.95);
	else             st->SetY1NDC(.875);
      }
      else               st->SetY1NDC(.85);
    }
    if (col) st->SetTextColor(col);
    st->Draw();
  }
}

// $Id: PlotFitTable.C,v 1.1 2021/03/15 17:28:43 ybedfer Exp ybedfer $

// Plot output of "rich_matrices/*/fit_table" (let's call it "rich.root").

/*
  // - "rich.root" contains TGraphErrors of performance:
  //   - efficiency Eh = hadron h ID'd as h
  //   - misidentification Mhh' = h ID'd ad h' != h.  
  // - This, for h(h') = pi, K, p, u(meaning ambiguous response)...
  // - ...in a multidimensional binning B = [charge] x [incidence] x [momentum].
  //  E.g. B can be = [+,-] x [0-3] x [0-14].
  // - "PlotRICHPerf" can:
  //     I) Plot Eh or Mh for [a subset of] B.
  //    II) Plot several such in a same TCanvas.
  //   III) Plot two such plus a comparison of the two.

  // Load and compile "PlotFitTable.C":
  .L PlotFitTable.C++

  // Online help:
  PlotFitTable();

  // Optional: set global title
  SetTitre("2016/P78910.sl71");

  // I) Plot Epi for B = [-] x [0-3]
  PlotRICHPerf("pi",0,0x1,0x2,0xf,0);

  //  II) Plot Epi, MpiK, Mpip for B = [+,-] x [1-2] in 3 pads
  PlotRICHPerf("pi",0,0xd,0x3,0x6,0);
  // II') Plot Epi, MpiK, Mpip for B = [+,-] x [1-2] in 2 pads 
  PlotRICHPerf("pi",0,0x3,0x3,0x6,0);

  // III) Compare EK for B = [-,+] x [1-2]: CSEvent vs. Julien in 2+1/2 pads...
  TFile *_file0 = TFile::Open("rich.root");
  TDirectory *dCS = gDirectory;
  PlotRICHPerf("K",0,0x1,0x3,0x6,2.5);
  TFile *_file1 = TFile::Open("../Julien/rich.root");
  TDirectory *dJu = gDirectory;
  PlotRICHPerf("K",0,0x31,0x3,0x3,-2); // "-2" is here somehow redundant... but nevertheless required
  //     ...Refinement: restrict P range (N.B.: click on pads to update them)
  hKCom->GetXaxis()->SetRange(16,78); SetFTRange();
  //     ...Compare also MKpi and MKp
  dCS->cd();
  PlotRICHPerf("K",0,0x2,0x3,0x6,2.5);
  dJu->cd();
  PlotRICHPerf("K",0,0x32,0x3,0x3,-2);
  hKpCom2->GetXaxis()->SetRange(16,78); SetFTRange();

  // III') Compare EK + vs. -
  dCS->cd();
  PlotRICHPerf("K",0,0x1, 0x1,0x6,2.5);
  PlotRICHPerf("K",0,0x11,0x2,0x6,-2);
  hKmCom2->GetXaxis()->SetRange(16,78); SetFTRange();
  PlotRICHPerf("K",0,0x2, 0x1,0x6,2.5);
  PlotRICHPerf("K",0,0x12,0x2,0x6,-2);
  hKmpCom2->GetXaxis()->SetRange(16,78); SetFTRange();
 */

#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TAxis.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"

static int black  = 1;
static int violet = kMagenta+2;
static int bleu   = kBlue+1;
static int cyan   = kCyan+1;
static int vert   = kGreen+2;
static int jaune  = kOrange+0;
static int orange = kOrange+1;
static int rouge  = kRed+0;

vector<TGraphErrors*> vGs;
static char *titre = 0;

bool getFTCanvas(const char *particle, unsigned int mode,
		 int &nPads, int hPad,
		 unsigned int ts, unsigned int signs, int binning,
		 TPad **&pads, char &tag);
void setFTAxes(TH2D *h, int nPads, int hPad, double paMn = 0, double paMx = 0);
void setFTLegend(TH2D *h, int effOrPur);
TGraphErrors *getFTGraph(const char *gN);
void compareFT(TGraphErrors *gi, TGraphErrors *gj);
void lowerK(char *name);
void upperK(char *name);
void setFTPaveText(TPaveText *pTxt,
		   unsigned int ts, unsigned int signs, int binning);
TPaveText *getFTPaveText(TPad *pad);

void PlotFitTable()
{
  printf("PlotRICHPerf: Plot output of \"rich_matrices/*/fit_table\"\n");
  printf("Usage: PlotRICHPerf(\n");
  printf("  const char *particle,    // pi, K,    p\n");
  printf("  const char *sample = 0,  // K0, ephi, Lambda, ...; D=pi:K0,K:ephi,p:Lambda\n");
  printf("  unsigned int mode = 0x3, // 0x1: Eff, 0x2: All Purities,\n");
  printf("                           // 0x4: Purity par->(par+1)%%3\n");
  printf("                           // 0x8:        par->(par+2)%%3\n");
  printf("                           // 0x10: Compare\n");
  printf("                           // 0x20: Julien's, 2-theta binning\n");
  printf("  unsigned int signs = 0x3,// 0x1 = +, 0x2 = -\n");
  printf("  unsigned int ts = 0x6,   // theta bins\n");
  printf("  int nPads = 0,           // =0: Derive # of pads from <mode>\n");
  printf("                           // >0 = # of pads, integer or half integer\n");
  printf("                           // <0: Use pad #(|nPads|-1); only for single pad <mod>,\n");
  printf("                           // >0 = #pads; has to accommodate <mode> at least\n");
  printf("\n");
}
int PlotRICHPerf(const char *particle,    // pi, K,    p
		 const char *sample = 0,  // K0, ephi, Lambda, ...
		 unsigned int mode = 0x3, // 0x1: Eff, 0x2: All Purities, 0x4: Purity pa->par+1%3, 0x10: Compare, 0x20: Julien's, 2-theta binning
		 unsigned int signs = 0x3,// 0x1 = +, 0x2 = -
		 unsigned int ts = 0x6,   // theta bins
		 double nPads = 0) // =0: Derive # of pads from <mode>
		                   // >0  = # of pad, integer or half integer
		                   // <0: Use pad #(|nPads|-1); only for single pad <mod>
{
  // Plot perf (efficiency,purity) from "rich.root"
  // - particles:         pi, K,    p
  // - default samples:   K0, ephi, Lamda
  const char *particles[] = {"pi",      "K",       "p"};    int npas = sizeof(particles)/sizeof(char*);
  const char *paNs[] =      {"#bf{#pi}","K",       "p"};
  const char *samples[] =   {"K0",      "iphi",    "Lambda"}; int nsas = sizeof(samples)/sizeof(char*);
  const char *saNs[] =      {"K0",      "i#varphi","#Lambda"};
  // ##### SIGN
  const char sCs[] = "pm", sNs[] = "+-"; char sC = '\0', sN = '\0';
  char **pmNs = new char*[npas];
  for (int ipa = 0; ipa<npas; ipa++) {
    size_t lP = strlen(paNs[ipa])+1;
    if (signs==0x3) {
      pmNs[ipa] = new char[lP];   strcpy(pmNs[ipa],paNs[ipa]);
    }
    else {
      if (signs==0x1) { sC = sCs[0]; sN = sNs[0]; }
      else            { sC = sCs[1]; sN = sNs[1]; }
      pmNs[ipa] = new char[lP+1]; snprintf(pmNs[ipa],lP+1,"%s%c",paNs[ipa],sN);
    }
  }

  // ########## BINNING
  double piMn = .855, piMx = 1.025, pimn = -.025, pimx = .08;
  double KMn =  .3,   KMx = 1.025,  Kmn =  -.025, Kmx =  .2;
  double paMn, paMx, pamn, pamx; if (!strcmp(particle,"pi")) {
    paMn = piMn; paMx = piMx; pamn = pimn; pamx = pimx;
  }
  else {
    paMn = KMn;  paMx = KMx;  pamn = Kmn;  pamx = Kmx;
  }

  // ########## CHECK ARGS
  // ##### <particle> and <sample> have to belong to known sets
  int ipa, jpa; for (ipa = 0, jpa = -1; ipa<npas; ipa++) {
    if (!strcmp(particle,particles[ipa])) { jpa = ipa; break; }
  }
  if (jpa<0) {
    printf("** PlotRICHPerf: Bad <particle> arg. (=\"%s\"): not in:",particle);
    for (ipa = 0; ipa<npas; ipa++) printf(" \"%s\"",particles[ipa]);
    printf("\n");
    return 1;
  }
  int jsa; if (sample) {
    int isa; for (isa = 0, jsa = -1; isa<nsas; isa++) {
      if (!strcmp(sample,samples[isa])) { jsa = isa; break; }
    }
    if (jsa<0) {
      printf("** PlotRICHPerf: Bad <sample> arg. (=\"%s\"): not in:",sample);
      for (isa = 0; isa<nsas; isa++) printf(" \"%s\"",samples[isa]);
      printf("\n");
      return 1;
    }
  }
  else
    jsa = jpa; // ##### DEFAULT SAMPLE

  // ##### MODE
  // 0xf0 bits: Derive "compare" and "binning", then clean away
  int compare = (mode&0x10) ? 1 : 0;
  int binning = (mode&0x20) ? 1 : 0;
  mode &= 0xf;
  if (compare && nPads>=0) {
    printf("** fit_table: \"compare\" option requires <0 \"nPads\" (as of v1.1)\n");
    return 1;
  }
 
  // ########## #PADS
  // "nPads", if >0, can be integer or half integer.
  int nPads2 = (int)round(2*nPads), iPads = nPads2/2, hPad = nPads2%2;
  if (nPads<0 && hPad) {
    printf("** PlotRICHPerf: Arg. <nPads> (=\"%.1f\") can be non-integer only if >0\n",nPads);
    return 1;
  }
  if (iPads>3 || iPads==3 && hPad) {
    printf("** PlotRICHPerf: Arg. <nPads> = %.1f: should be <3 \n",nPads);
    return 1;
  }
  int mPads, im; for (im=mPads = 0; im<4; im++) if (0x1<<im&mode) mPads++; 
  // <nPads><0: Accept only single pad mode
  if (nPads<0) {
    if (mPads>1) {
      printf("** PlotRICHPerf: Arg. <mode>(=0x%x) requires %d pads => cannot accept arg. <nPads>=0\n",mode,mPads);
      return 1;
    }
  }
  if (iPads==0) iPads = mPads;
  // <nPads> ENOUGH TO ACCOMODATE <mode>
  if (nPads>0 && nPads<mPads) {
    printf("** PlotRICHPerf: Arg. <mode>(=0x%x) requires %d pads > arg. <nPads>(=%.1f)\n",mode,mPads,nPads);
    return 1;
  }

  // ########## CANVAS
  TPad **pads; char tag; int kPads = iPads;
  if (getFTCanvas(particle,mode,kPads,hPad,ts,signs,binning,pads,tag))
    return 1;
  int ipad; if (nPads<0) { // If "nPads" <0, get pad from getFTCanvas
    ipad = -nPads-1; // nPads=-1 means 1st pad, which is that w/ highest index (since w/ build pad bottom->up, see "stackPads").
  }
  else
    ipad = 0;

  // ########## EFFICIENCY
  if (mode&0x1) {
    char hN[] = "hpimEff00"; size_t lN = strlen(hN)+1;
    char hT[] = "#bar{#Lambda}: #pi^{-} #rightarrow #pi^{-};#it{P} (GeV)  "; size_t lT = strlen(hT)+1;
    // ##### HISTO
    if (signs==0x3) snprintf(hN,lN,"h%sEff%d%c",  particle,   ipad,tag);
    else            snprintf(hN,lN,"h%s%cEff%d%c",particle,sC,ipad,tag);
    snprintf(hT,lT,"%s: %s #rightarrow %s;#it{P} (GeV)  ",
	     saNs[jsa],pmNs[jpa],pmNs[jpa]);
    TH2D *hEff = new TH2D(hN,hT,100,3,50,100,-.025,1.025); hEff->SetStats(0);
    pads[ipad]->cd(); gPad->SetGridy();
    hEff->Draw();
    setFTAxes(hEff,iPads+hPad,0,paMn,paMx); setFTLegend(hEff,0);
    // ##### GRAPH
    TGraphErrors *gi; char giN[] = "pi_m_pi_0"; size_t lGi = strlen(giN)+1;
    int cols[] = {black,violet,bleu,cyan,vert,jaune,orange,rouge};
    int firstT, t; for (t = 0, firstT=1; t<4; t++) {
      if (!(0x1<<t&ts)) continue;
      int jt = binning ? t+1 : t; 
      for (int s = 0; s<2; s++) {
	if (!(0x1<<s&signs)) continue;
	pads[ipad]->cd();
	snprintf(giN,lGi,"%s_%c_%s_%d",particle,sCs[s],particle,t);
	if (binning) lowerK(giN);
	if (!(gi = (TGraphErrors*)gDirectory->Get(giN))) {
	  printf("** PlotRICHPerf: No \"%s\" TGraphErrors in \"%s\" TDirectory\n",
		 giN,gDirectory->GetName());
	}
	else {
	  // Build a new TGraph w/ abscissae shifted so that data points do not
	  // hide each other.
	  TGraphErrors *g = new TGraphErrors();
	  g->SetMarkerStyle(20); g->SetMarkerSize(.75);
	  char gN[] = "gpipT99_99"; size_t lG = strlen(gN)+1;
	  snprintf(gN,lG,"g%s%cT%d_%d%c",particle,sCs[s],jt,ipad,tag);
	  upperK(gN);
	  g->SetName(gN); vGs.push_back(g);
	  for (int pt = 0; pt<gi->GetN(); pt++) {
	    double x, y, dy; gi->GetPoint(pt,x,y); dy = gi->GetErrorY(pt);
	    // If both signs, <0 on the left, >0 on the right
	    // And t=0 and t=4 also oon either side but farther away.
	    double dx; if (signs==0x3) {
	      if (ts==0x6 || ts==0x3) {
		double dxs[] = {.040,.120,.2,.28}; dx = (1-2*s)*dxs[jt];
	      }
	      else {
		double dxs[] = {.020,.060,.1,.14}; dx = (1-2*s)*dxs[jt];
	      }
	    }
	    else {
	      if (ts==0x6 || ts==0x3) {
		double dxs[] = {-.150,-.050,.050,.150};  dx = dxs[jt];
	      }
	      else {
		double dxs[] = {-.075,-.025,.025,.075};  dx = dxs[jt];
	      }
	    }
	    g->SetPoint(pt,x+dx,y); g->SetPointError(pt,0,dy);
	  }
	  int col = cols[2*jt+s];
	  g->SetMarkerColor(col); g->SetLineColor(col);
	  g->Draw("e1p");
	  if (compare) {        // ***** COMPARISON
	    TGraphErrors *gj; for (int is = 0; is<2; is++) {
	      int sp = is ? 1-s : s; // Sign? Try first "s", then "1-s"
	      snprintf(gN,lG,"g%s%cT%d_%d%c",particle,sCs[sp],jt,ipad-1,tag);
	      upperK(gN);
	      gj = getFTGraph(gN); if (gj) {
		if (is && firstT) {
		  TPaveText *pTxt = getFTPaveText(pads[0]);
		  if (pTxt) {
		    pTxt->Clear(); setFTPaveText(pTxt,ts,0x3,binning);
		    TCanvas *c = pads[0]->GetCanvas(); c->cd(); pTxt->Draw();
		  }
		  else
		    printf("** PlotRICHPerf: No updating global TPaveText title\n");
		}
		break;
	      }
	    }
	    if (!gj) {
	      for (int is = 0; is<2; is++) {
		int sp; if (is) { printf(" nor ");               sp = 1-s; }
		else { printf("** fit_table: No TGraphErrors "); sp = s; }
		snprintf(gN,lG,"g%s%cT%d_%d%c",particle,sCs[sp],jt,ipad-1,tag);
		printf("\"%s\"",gN);
	      }
	      printf(" in \"vGs\" => No comparison\n");
	    }
	    else  {
	      pads[ipad+1]->cd();
	      if (firstT) {
		char HN[] = "hpipCom0"; size_t LN = strlen(HN)+1;
		if (signs==0x3)
		  snprintf(HN,LN,"h%sCom%c",  particle,   tag);
		else
		  snprintf(HN,LN,"h%s%cCom%c",particle,sC,tag);
		TH2D *hCom = new TH2D(HN,hT,100,3,50,100,-.15,.15); hCom->SetStats(0);
		pads[ipad+1]->cd(); gPad->SetGridy();
		hCom->Draw();
		setFTAxes(hCom,iPads+hPad,1);
		// Y axis range is set wide in ctor. It's here narrowed...
		// ...but the user can always widen it again.
		hCom->GetYaxis()->SetRange(36,66);
		firstT = 0;
	      }
	      compareFT(g,gj);
	    }
	  }
	}
      }
    }
    ipad++;
  }
  // ########## PURITY
  if (mode&0xe) {
    char hN[] = "hpimPur00"; size_t lN = strlen(hN)+1;
    char hT[] = "#bar{#Lambda}: #pi^{-} #rightarrow #pi^{-}/#pi^{-};#it{P} (GeV)  "; size_t lT = strlen(hT)+1;
    int kpa = (jpa+1)%3, lpa = (jpa+2)%3; bool allAtOnce = (mode&0xe)==0x2;

    for (int ipa = 0; ipa<2; ipa++) {
      int misID = 0x4<<ipa; int mpa = ipa?lpa:kpa;
      if (allAtOnce) { if (ipa>0) break; }
      else if (!(mode&misID)) continue;
      // ##### HISTOS
      if (allAtOnce) {
	if (signs==0x3)
	  snprintf(hN,lN,"h%sPur%d%c",  particle,   ipad,tag);
	else
	  snprintf(hN,lN,"h%s%cPur%d%c",particle,sC,ipad,tag);
	snprintf(hT,lT,"%s: %s #rightarrow %s/%s;#it{P} (GeV)  ",
		 saNs[jsa],pmNs[jpa],pmNs[kpa],pmNs[lpa]);
      }
      else {
	if (signs==0x3)
	  snprintf(hN,lN,"%s%s%d%c",  particle,   particles[mpa],ipad,tag);
	else
	  snprintf(hN,lN,"%s%c%s%d%c",particle,sC,particles[mpa],ipad,tag);
	snprintf(hT,lT,"%s: %s #rightarrow %s;#it{P} (GeV)  ",
		 saNs[jsa],pmNs[jpa],pmNs[mpa]);
      }
      double dy = .001;// Avoid top label to be eaten by TPad or TPaveText above
      TH2D *hPur = new TH2D(hN,hT,100,3,50,100,-.025,1-dy); hPur->SetStats(0);
      pads[ipad]->cd(); gPad->SetGridy();
      hPur->Draw();
      setFTAxes(hPur,iPads+hPad,0,pamn,pamx); setFTLegend(hPur,1);
      // ##### GRAPH
      TGraphErrors *gi; char gN[] = "pi_m_pi_0"; size_t lG = strlen(gN)+1;
      int cols[] = {black,violet,bleu,cyan,vert,jaune,orange,rouge};
      int firstT, t; for (t = 0, firstT = 1; t<4; t++) {
	if (!(0x1<<t&ts)) continue;
	int jt = binning ? t+1 : t; 
	for (int s = 0; s<2; s++) {
	  if (!(0x1<<s&signs)) continue;
	  for (int kl = 0; kl<2; kl++) {
	    pads[ipad]->cd();
	    if (!allAtOnce) {
	      if (kl!=ipa) continue;
	    }
	    const char *partjcle = kl?particles[lpa]:particles[kpa];
	    snprintf(gN,lG,"%s_%c_%s_%d",particle,sCs[s],partjcle,t);
	    if (binning) lowerK(gN);
	    if (!(gi = (TGraphErrors*)gDirectory->Get(gN))) {
	      printf("** PlotRICHPerf: No \"%s\" TGraphErrors in \"%s\" TDirectory\n",
		     gN,gDirectory->GetName());
	    }
	    else {
	      // Build a new TGraph w/ abscissae shifted so that data points do
	      // not hide each other.
	      TGraphErrors *g = new TGraphErrors();
	      g->SetMarkerStyle(20+4*kl); g->SetMarkerSize(.75);
	      char gN[] = "gpipKT99_99"; size_t lG = strlen(gN)+1;
	      snprintf(gN,lG,"g%s%c%sT%d_%d%c",particle,sCs[s],partjcle,jt,ipad,tag);
	      upperK(gN);
	      g->SetName(gN); vGs.push_back(g);
	      //printf("___ pushing <%s> => %d\n",gN,(int)vGs.size());
	      for (int pt = 0; pt<gi->GetN(); pt++) {
		double x, y, dy; gi->GetPoint(pt,x,y); dy = gi->GetErrorY(pt);
		// If both signs, <0 on the left, >0 on the right
		// And t=0 and t=4 also oon either side but farther away.
		double dx; if (signs==0x3) {
		  if (ts==0x6 || ts==0x3) {
		    double dxs[] = {.040,.120,.2,.28}; dx = (1-2*s)*dxs[jt];
		  }
		  else {
		    double dxs[] = {.020,.060,.1,.14}; dx = (1-2*s)*dxs[jt];
		  }
		}
		else {
		  if (ts==0x6 || ts==0x3) {
		    double dxs[] = {-.150,-.050,.050,.150};  dx = dxs[jt];
		  }
		  else {
		    double dxs[] = {-.075,-.025,.025,.075};  dx = dxs[jt];
		  }
		}
		g->SetPoint(pt,x+dx,y); g->SetPointError(pt,0,dy);
	      }
	      int col = binning==1 ? cols[2*t+2+s] : cols[2*t+s];
	      g->SetMarkerColor(col); g->SetLineColor(col);
	      g->Draw("e1p");
	      if (compare) {        // ***** COMPARISON
		TGraphErrors *gj; for (int is = 0; is<2; is++) {
		  int sp = is ? 1-s : s; // Sign? Try first "s", then "1-s"
		  snprintf(gN,lG,"g%s%c%sT%d_%d%c",particle,sCs[sp],partjcle,jt,ipad-1,tag);
		  upperK(gN);
		  gj = getFTGraph(gN); if (gj) {
		    if (is && firstT) {
		      TPaveText *pTxt = getFTPaveText(pads[0]);
		      if (pTxt) {
			pTxt->Clear(); setFTPaveText(pTxt,ts,0x3,binning);
			TCanvas *c = pads[0]->GetCanvas(); c->cd(); pTxt->Draw();
		      }
		      else
			printf("** PlotRICHPerf: No updating global TPaveText title\n");
		    }
		    break;
		  }
		}
		if (!gj) {
		  for (int is = 0; is<2; is++) {
		    int sp; if (is) { printf(" nor ");               sp = 1-s; }
		    else { printf("** fit_table: No TGraphErrors "); sp = s; }
		    snprintf(gN,lG,"g%s%cT%d_%d%c",particle,sCs[sp],jt,ipad-1,tag);
		    printf("\"%s\"",gN);
		  }
		  printf(" in \"vGs\" => No comparison\n");
		}
		else  {
		  pads[ipad+1]->cd();
		  if (firstT) {
		    char HN[] = "hpipKCom0"; size_t LN = strlen(HN)+1;
		    if (signs==0x3)
		      snprintf(HN,LN,"h%s%sCom%c",  particle,   partjcle,tag);
		    else
		      snprintf(HN,LN,"h%s%c%sCom%c",particle,sC,partjcle,tag);
		    upperK(HN);
		    TH2D *hCom = new TH2D(HN,hT,100,3,50,100,-.15,.15); hCom->SetStats(0);
		    pads[ipad+1]->cd(); gPad->SetGridy();
		    hCom->Draw();
		    setFTAxes(hCom,iPads+hPad,1);
		    // Y axis range is set wide in ctor. It's here narrowed...
		    // ...but the user can always widen it again.
		    hCom->GetYaxis()->SetRange(36,66);
		    firstT = 0;
		  }
		  compareFT(g,gj);
		}
	      }
	    }
	  }
	}
      }
      if ((mode&0xc)==0xc) ipad++;
    }
  }
  if (mode&0x10) // mode&0xà1, i.e. comparison requested. => Clear "vGs", now
    // that it fulfilled its purpose of storing the first leg of the comparison.
    vGs.clear();
  return 0;
}
string getFitTabTitle(unsigned int ts, unsigned int signs, int binning)
{
  if (binning==1) // Julien's binning: 0x1 means [0.01,0.04], i.e. CSEvent 0x2
    ts = ts<<1;
  string title;
  const char *cThs[4][2] = {{"[0,.01]",               "#color[618]{[0,.01]}"},
			    {"#color[601]{[.01,.04]}","#color[433]{[.01,.04]}"},
			    {"#color[418]{[.04,.12]}","#color[800]{[.04,.12]}"},
			    {"#color[801]{[.12,.3]}", "#color[632]{[.12,.3]}"}};
  for (int s = 0; s<2; s++) {
    if (!(0x1<<s&signs)) continue;
    if (!s) title += "  #bf{#it{h}}^{#plus}: #scale[.8]{";
    else    title += "  #bf{#it{h}}^{#minus}: #scale[.8]{";    
    for (int t = 0; t<4; t++)
      if (0x1<<t&ts) title += string(cThs[t][s]);
    title += "}";
  }
  return title;
}
void stackPads(const char *genericPadName, int nPads, int hPad, TPad **pads)
{

  gStyle->SetOptTitle(0);

  Double_t W  = 0.95, H;
  if (hPad==0) {
    if      (nPads==3) H = 0.29; // Pad Width, Height
    else if (nPads==2) H = 0.43;
    else               H = 0.86;
  }
  else {
    if      (nPads==2) H = 0.35;
    else               H = 0.57;
  }

  Int_t    Nx = 1; double Ny = nPads+.5*hPad;  // Number of pads along X, Y
  Double_t Xm = (1-(Nx*W))/2; // X Margin
  Double_t Ym = (1-(Ny*H))/2; // Y Margin   Double_t dw = (W*0.1)/4;
  Double_t dw = (W*0.1)/4;
  Double_t dh; if (nPads==3) dh = H*0.5;
  else         if (nPads==2) {
    if (hPad)                dh = H*0.52;
    else                     dh = H*0.3;
  }
  else                       dh = H*0.2;
   
  char tag[] = "_1_i";
  double Xlow = Xm, Xup = Xlow+W+dw;
  double Ylow = Ym, Yup;
  if (hPad) Yup = Ylow+H/2+dh/3;
  else      Yup = Ylow+H+dh/3;
  int nPadsTot = nPads+hPad; for (int i = 0; i<nPadsTot; i++) {
    sprintf(tag,"_1_%d",i);
    string pS(genericPadName); pS += string(tag); const char *pN = pS.c_str();
    int j = nPadsTot-1-i;
    //printf("%d,%d => %d %.4f,%.4f %.4f,%.4f,%.4f\n",
    //i,nPadsTot,j,Xlow,Xup,Ylow,Yup,H);
    pads[j] = new TPad(pN,pN,Xlow,Ylow,Xup,Yup, 0, 0, 0);
    pads[j]->SetTopMargin(0);   
    if (i==0) pads[j]->SetBottomMargin(dh);
    else      pads[j]->SetBottomMargin(0);
    /*
      if (i<=1) pads[i]->SetTopMargin(0);
      if (i>=1) pads[i]->SetBottomMargin(0);
      else      pads[i]->SetBottomMargin(0.18);
    */
    pads[j]->Draw();
    Ylow = Yup;
    Yup += H;
    /*
      if (i==0) Yup -= 2*dh;
      else      Yup += dh;
    */
  }
}
void setFTAxes(TH2D *h, int nPads, int hPad, double paMn, double paMx) {
  TAxis *ax = h->GetXaxis(); ax->SetNdivisions(505);
  double siz; if (hPad==0) {
    if      (nPads) siz = 0.053;
    else if (nPads) siz = 0.05;
    else            siz = 0.04;
  }
  else              siz = 0.08;
  ax->SetLabelSize(siz); ax->SetTitleSize(siz*1.2); ax->SetTitleOffset(.9);
  TAxis *ay = h->GetYaxis(); ay->SetNdivisions(505);
  ay->SetLabelSize(siz);
  if (paMn) {
    int nBins = ay->GetNbins(); double xMn = ay->GetXmin(), xMx = ay->GetXmax();
    double dy = xMx-xMn; int first = nBins*(paMn-xMn)/dy, last = nBins*(paMx-xMn)/dy;
    ay->SetRange(first,last);
  }
}
void setFTLegend(TH2D *h, int effOrPur)
{  
  double x1, x2, y1, y2; // Has to be ~centre, since you may want to cut away
  // the left, low P, part of the histogram
  if (effOrPur==0) {
    x1 = .4; x2 = .6; y1 = .2; y2 = .3;  // Centre-bottom
  }
  else {
    x1 = .4; x2 = .6; y1 = .8; y2 = .9;  // Centre-top
  }
  TPaveText *pTxt = new TPaveText(x1,y1,x2,y2,"NDC");
  pTxt->InsertText(h->GetTitle());
  pTxt->Draw();
}
#include "TClass.h"
#include "TList.h"
bool getFTCanvas(const char *particle, unsigned int mode,
		 int &nPads, // <0: use current, pre-existing TCanvas; upon return total number of pads
		 int hPad,   // half-pad
		 unsigned int ts, unsigned int signs, int binning,
		 TPad **&pads, char &tag)
{
  if (nPads<0) {
    if (!gPad) {
      printf("** PlotRICHPerf: Arg. <nPads> is <0 while no TCanvas pre-exists. => Exiting...\n");
      return 1;
    }
    TCanvas *c = gPad->GetCanvas();
    const char *cPerf; if (!(cPerf = strstr(c->GetName(),"Perf"))) {
      printf("** PlotRICHPerf: Arg. <nPads> is <0 while no current TCanvas(=\"%s\") is not a PlotFitTable one\n",c->GetName());
      return 1;
    }
    if (strlen(cPerf)>4) tag = *(cPerf+4);
    else                 tag = '\0';
    TList *l = c->GetListOfPrimitives(); TObject *o = l->First(); TPad *pad = 0;
    do {
      if (o->IsA()->InheritsFrom(TPad::Class())) { pad = (TPad*)o; break; }
    } while ((o = l->After(o)));
    if (!pad) {
      printf("** PlotRICHPerf: Arg. <nPads> is <0 while current TCanvas(=\"%s\") contains no TPad\n",c->GetName());
      return 1;
    }
    l = pad->GetListOfPrimitives(); o = l->First();
    int nPadsTot = 0;  do {
      if (o->IsA()->InheritsFrom(TPad::Class())) nPadsTot++;
    } while ((o = l->After(o)));
    if (-nPads>nPadsTot) {
      printf("** PlotRICHPerf: Arg. <nPads> is %d while current TCanvas(=\"%s\") contains only %d TPad's\n",nPads,c->GetName(),nPadsTot);
      return 1;
    }
    pads = new TPad*[nPadsTot];
    int ipad = 0; o = l->First(); do {
      if (o->IsA()->InheritsFrom(TPad::Class()))
	pads[nPadsTot-ipad-1] = (TPad*)o; ipad++;
    } while ((o = l->After(o)));
    nPads = nPadsTot;
    return 0;
  }
  // If TCanvas already existing, create new one w/ same name + '2'
  char cN[] = "cpimPerf0"; size_t lC = strlen(cN)+1; tag = '\0';
  int mPads, im; for (im=mPads = 0; im<4; im++) if (0x1<<im&mode) mPads++; 
  if (signs==0x3 ||
      // If nPads>mPads and sign is given, we can expect next execution to be
      // the complementary sign => Don't specify the charge in TCanvas name
      signs!=0x3 && (nPads>mPads || nPads==0))
    snprintf(cN,lC,"c%sPerf",particle);
  else
    snprintf(cN,lC,"c%s%cPerf",particle,signs==0x1?'p':'m');
  if (gROOT->FindObject(cN)) {
    tag = '2'; sprintf(cN+strlen(cN),"%c",tag);
  }
  double nPadsEff = nPads+.5*hPad;
  TCanvas *c = new TCanvas(cN,cN,896,256+nPadsEff*176);
  // TPaveText w/ global title
  c->Divide(1,1); c->cd(1); TPad *pad = (TPad*)gPad;
  double x1NDC, y1NDC, x2NDC, y2NDC; pad->GetPadPar(x1NDC,y1NDC,x2NDC,y2NDC);
  double dy;
  if      (nPads==3) dy = 0.02;
  else if (nPads==2) dy = 0.02;
  else               dy = 0.02;
  y2NDC -= .02 + (4-nPadsEff)*dy;
  pad->SetPad(x1NDC,y1NDC,x2NDC,y2NDC);
  TPaveText *pTxt = new TPaveText(.1,y2NDC+.005,.9,.995,"NDC");
  cN[0] = 'p'; pTxt->SetName(cN);
  setFTPaveText(pTxt,ts,signs,binning);
  c->cd(); pTxt->Draw();
  pad->cd(); pads = new TPad*[nPads]; stackPads(cN,nPads,hPad,pads);
  return 0;
}
TGraphErrors *getFTGraph(const char *gN)
{
  // Get TGraphErrors w/ name <gN> from global "vGs".
  TGraphErrors *g; int ig; for (ig = 0, g = 0; ig<(int)vGs.size(); ig++) {
    TGraphErrors *gi = vGs[ig]; const char *giN = gi->GetName();
    //printf("==== <%s>\n",giN);
    if (!strcmp(gN,giN)) { g = gi; break; }
  }
  return g;
}
void compareFT(TGraphErrors *gi, TGraphErrors *gj)
{
  printf(" * compareFT: \"%s\" - \"%s\"\n",gi->GetName(),gj->GetName());
  TGraphErrors *g;
  int ni = gi->GetN(), nj = gj->GetN(); int i, k;
  for (i=k = 0, g = 0; i<ni; i++) {
    double xi, yi, dyi; gi->GetPoint(i,xi,yi); dyi = gi->GetErrorY(i);
    for (int j = 0; j<nj; j++) {
      double xj, yj, dyj; gj->GetPoint(j,xj,yj); dyj = gj->GetErrorY(j);
      if (xj>xi+.45) // I.e. larger than any of the artificial shifts (see "dx"
	// supra) and smaller than any of the the p bin widths
	break;
      else if (xj<xi-.45)
	continue;
      double y = yi-yj, dy = sqrt(dyi*dyi+dyj*dyj);
      if (!g) {
	g = new TGraphErrors(*gi);
	const char *nami = gi->GetName();
	const char *namj = gj->GetName()+/* exclude leading 'g' */1;
	string name = string(nami)+string(namj); g->SetName(name.c_str());
      }
      g->SetPoint(k,xi,y); g->SetPointError(k,0,dy); k++;
    }
  }
  //g->Print();
  g->Draw("e1p");
  return;
}
void lowerK(char *name)
{
  char *c = name; do {
    if (*c=='K') *c = 'k';
  } while (*(++c)!='\0');
}
void upperK(char *name)
{
  char *c = name; do {
    if (*c=='k') *c = 'K';
  } while (*(++c)!='\0');
}
int SetFTRange()
{
  // Set the range of all pads = to that of bottom pad.
  // Bottom pad is the first in the list (see how pads are built in "stackPads").

  // ########## GET TCanvas
  TCanvas *c = gPad->GetCanvas(); if (!c) {
    printf("** SetRange: No parent TCanvas to current TPad \"%s\"\n",
	   gPad->GetName());
    return 1;
  }
  // ########## GET BOTTOM TPad
  // ##### BACKGROUND PAD (i.e. pad obtained w/ "Divide(1,1)")
  TList *l = c->GetListOfPrimitives(); TObject *o = l->First(); TPad *pad0 = 0;
  do {
    if (o->IsA()->InheritsFrom(TPad::Class())) { pad0 = (TPad*)o; break; }
  } while ((o = l->After(o)));
  if (!pad0) {
    printf("** SetRange: Parent TCanvas \"%s\" contains no TPad\n",
	   c->GetName());
    return 1;
  }
  // ##### BOTTOM PAD
  l = pad0->GetListOfPrimitives(); o = l->First(); TPad *pad = 0;
  do {
    if (o->IsA()->InheritsFrom(TPad::Class())) { pad = (TPad*)o; break; }
  } while ((o = l->After(o)));
  if (!pad) {
    printf("** SetRange: Background TPad \"%s\" contains no TPad\n",
	   pad0->GetName());
    return 1;
  }
  printf(" * SetRange: All pads of TCanvas \"%s\" aligned on TPad \"%s\"\n",
	 c->GetName(),pad->GetName());
  // ########## HISTO in BOTTOM PAD
  TList *lp = pad->GetListOfPrimitives(); TObject *op = lp->First(); TH2D *h2 = 0;
  do {
    if (op->IsA()->InheritsFrom(TH2D::Class())) { h2 = (TH2D*)op; break; }
  } while ((op = lp->After(op)));
  if (!h2) {
    printf("** SetRange: Bottom TPad \"%s\" contains no TH2D\n",pad->GetName());
    return 1;
  }
  TAxis *ax = h2->GetXaxis(); int first = ax->GetFirst(), last = ax->GetLast();
  // ########## LOOP ON REST OF PADS
  int nPads = 0; while ((o = l->After(o))) {
    if (o->IsA()->InheritsFrom(TPad::Class())) {
      pad = (TPad*)o;
      lp = pad->GetListOfPrimitives(); op = lp->First(); TH2D *h2 = 0;
      do {
	if (op->IsA()->InheritsFrom(TH2D::Class())) {
	  h2 = (TH2D*)op; ax = h2->GetXaxis(); ax->SetRange(first,last);
	  pad->Update(); nPads++; break;
	}
      } while ((op = lp->After(op)));
      if (!h2) {
	printf("** SetRange: TPad \"%s\" contains no TH2D\n",pad->GetName());
      }
    }
  }
  if (nPads) {
    printf(" * SetRange: %d pad's found and aligned\n",nPads); return 0;
  }
  else {
    printf(" * SetRange: No pad found\n"); return 1;
  }
}
void SetTitre(const char *text)
{
  size_t len = strlen(text)+1; titre = new char[len];
  strcpy(titre,text);
}
void setFTPaveText(TPaveText *pTxt,
		   unsigned int ts, unsigned int signs, int binning)
{
  string title = getFitTabTitle(ts,signs,binning);
  if (titre) {
    string longTitle = string(titre)+" - "+title;
    pTxt->InsertText(longTitle.c_str());
  }
  else
    pTxt->InsertText(title.c_str());
}
TPaveText *getFTPaveText(TPad *pad)
{
  TCanvas *c = gPad->GetCanvas(); if (c) {
    TList *l = c->GetListOfPrimitives();
    TObject *o = l->First(); do {
      if (o->IsA()->InheritsFrom(TPaveText::Class())) {
	return (TPaveText*)o;
      }
    } while ((o = l->After(o)));
  }
  return 0;
}

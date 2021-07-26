#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <cmath>

#include <TROOT.h>

#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixD.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>
#include <TPDF.h>

#include <RooGlobalFunc.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooAddition.h>
#include <RooCategory.h>
#include <RooChebychev.h>
#include <RooConstVar.h>
#include <RooDataHist.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooMinuit.h>
//#include <RooMy.h>
#include <RooPlot.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooBreitWigner.h>
#include <RooRelBreitWigner.hh>
#include <RooSimultaneous.h>
#include <RooVoigtian.h>


/*
 *
 *	const TMatrixDSym& cor = r->correlationMatrix() ;
 *	const TMatrixDSym& cov = r->covarianceMatrix() ;
 *	// Print correlation, covariance matrix
 *	cor.Print() ;
 *	cov.Print() ;
 *
 *
 *	// Make list of model parameters
 *	RooArgSet* params = model.getParameters(x) ;
 *	// Write LaTex table to file
 *	params->printLatex(Sibling(*initParams),OutputFile("rf407_latextables.tex")) ;
 *
 *
 *
 */



// ******************************************************************************************

using namespace std;

using std::cout;
using std::cerr;
using std::endl;
using std::map;
using std::stringstream;
using std::vector;
using std::pair;
using std::string;
using std::make_pair;
using std::ios;
using std::setw;
using std::setprecision;


using namespace RooFit ;

// ******************************************************************************************
// const int Np = 6;
// const int Nt = 7;
// const double p_bins[Np+1]={10. ,15. ,20. ,25. ,30. ,40. ,50.};
// const double t_bins[Nt+1]={0.00,0.01,0.02,0.03,0.04,0.06,0.09,0.12};

// const int Np = 13;
// const int Nt = 4;
// const double p_bins[Np+1]={10.,11.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.,50.};
// const double t_bins[Nt+1]={0.00,0.01,0.04,0.12,0.3};

const int Np = 15;
const int Nt = 4;
const double p_bins[Np+1]={3.,5.,7.,10.,12.,13.,15.,17.,19.,22.,25.,27.,30.,35.,40.,50.};
const double t_bins[Nt+1]={0.00,0.01,0.04,0.12,0.3};

const string chan[8] = {"K0_pip","K0_pim","phi_kp","phi_km","Lambda_pip","Lambda_pim","ephi_kp","ephi_km"};

// const int Np = 1;
// const int Nt = 1;
// const double p_bins[Np+1]={10.,50.};
// const double t_bins[Nt+1]={0.01,0.12};

double N_id[8][6][Np][Nt]; // [channels][a,pi,K,p,u and background]
void initCounts();
//double R_id[6][4][Np][Nt];

double last_para[30];
double last_para2[Np][Nt][30];

string analysis;
string data_file;
string data_template;
int data_ff_nb;
int data_lf_nb;
int data_nb;
string fit_type;
string hist_file_K0 =   "hist_K0.root";
string hist_file_Lam =  "hist_Lambda.root";
string hist_file_iphi = "hist.iphi.root";
string hist_file_ephi = "hist.ephi.root";
string out_file = "rich.root";
int id_lst[5]; double lh_cut[5][6]; // LikeliHood cuts
TH1D* h[8][5][Np][Nt];
TH2D* h2[8][5][Np][Nt];
// Kinematics histos
TH2D *am_all, *am_K0, *am_L;
TH2D *am_K0p, *am_K0m;
TH1D *DdD_K0, *DdD_L, *cth_K0, *cth_L, *pT_K0, *pT_L;
TH1D *Z_all,  *Z_K0,  *Z_L;
TH2D *XY_all, *XY_K0, *XY_L;
TH1D *Tr_K0,  *Tr_L,  *Rb_K0,  *Rb_L,  *Rt_K0, *Rt_L;
TH2D *am_Incl, *am_Excl, *am_Iphi, *am_Ephi;
TH1D *pT_Incl, *pT_Excl, *dE_Incl, *dE_Excl;
TH1D *Z_Iphi,  *Z_Ephi;
TH2D *XY_Iphi, *XY_Ephi;
TH1D *pT_Iphi, *pT_Ephi, *dE_Iphi, *dE_Ephi;
TH1D *Tr_Iphi, *Tr_Ephi, *Rb_Iphi, *Rb_Ephi, *Rt_Iphi, *Rt_Ephi;
// Kinematics cuts
double DdD_cuts[2], cth_cuts[2]; // 0: K0, 1: Lambda.
double pT_cuts[4];               // 0: K0, 1: Lambda, 2: Incl. phi, 3: Excl. phi 
double dE_cuts[2];               // 0: !Incl., 1: Excl.
RooFitResult *r[6][Np][Nt];
double Ns[6][Np][Nt], Bs[6][Np][Nt];
int lw = 1;
double thr_diff = 0.;
int retry = 20;
bool rpipe = false;

bool use_improve = false;
bool use_hesse = true;
bool use_minos = false;
bool use_sidebins = true;
TFile* input;
TFile* input_K0;
TFile* input_iphi;
TFile* input_Lam;
int  verbose;

// CSEvenData TTree
#include "CSEventData.h"
CSEventData *ev;
vector<CSHadronData> *hadrons;
vector<CSResonanceData> *resonances;

// ******************************************************************************************

int main(int, char**);
bool get_inputFile(int pi);
bool read_options(string optFile);
void get_input_data_t1();
void get_input_data_t2();
void get_input_data2();
void bookKineHistos();
void writeKineHistos(const char *particleName);

RooDataHist* gen_K0(int, int , int, int);
void write_hist();
void create_hist();
void get_plots();
void fit_table_K0(int);
void fit_table_phi(int);
void fit_table_Lambda(int);
void print_table();
void set_plot_style();

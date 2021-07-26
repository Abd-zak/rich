#include "fit_table.h"

// ******************************************************************************************
//
// Using the selected data for K0, Lambda and phi in a simultaneous fit to
// produce the RICH Table with pi, p and K.
//
// ******************************************************************************************

// OUTLINE OF fit_table CODE (as of 2021/05)
// - plots
//   - "get_inputFile"
//   - "get_input_data_t1": Read in input TTree's for K0/Lambda
//   - "get_input_data_t2": Read in input TTree's for phi
//   - "(create|write)_hist" to handle output histos.
//   - HISTOS "h" of INV.MASS are indexed by [channel][id][Pbin][Tbin], where
//    "id" id a,pi,K,p,u 'a' means 'all', i.e. noID, and 'u' means 'unknown'.
//     Same for "h2's", which contain Armenteros plots.
//   - CHANNELS "chan[]" specify:
//     i) The DECAYING, neutral, PARTICLE, p0 = K0, Lambda or phi.
//    ii) The SPECTATOR decay PARTICLE, pS, i.e. the counterpart of the decay
//       PARTICLE UNDER EXAM, pE.
//    "K0_pip" means p0 = K0, pS = pi+ and therefore pE = pi-. And mutatis
//    mutandis for the rest.
//     Channels are defined in the header file.
//   - "getPID" returns particle's ID, base on LikeliHood cuts defined in
//    options file.
//   - The spectator is subjected to an ID-BASED SELECTION. The selection
//    criterion is the same for + and -. But changes from channel to channel
//    otherwise.
//      - For p0 = K0, pi-ID is required for pS. May not be optimum, given that:
//        - Pions are the most numerous of all particles.
//        => Rejecting then the rest of the particles does not improve much
//          Signal/Background.
//        - At high momentum P, p/K cannot be resolved and RICH response tends
//         to be 'unknown'.
//        => Pion selection cuts significantly into overall statistics.
//      - For p0 = phi, pi-veto is required for pS (or so I(Y.B.) understand).
//       Again, this cuts into the statistics. But this time, it does improve
//       S/B. Remains to estimate whether this improves the statistical
//       significance (i.e. S/sqrt(S+B)) and hence improve the precision on RICH
//       performance. A trade-off could be envisioned, whereby pi-veto is only
//       required for lower P of pS. This could improve precision at high P of
//       pE since P of K+ and K- from phi decay tend to be close to one another
//       given the small mass diff. Mphi-2*MK.
//      - For p0 = Lambda, pi-ID is required (or so I(Y.B.) understand). Again
//       that may not be helpful (for the same reason as in the p0 = K0 case).
// - fit
//   - "get_plots": Read in histos to be fitted
//   - "fit_table_(K0|Lambda|phi)": Fit histos. The excl. phi channel is not
//    yet processed.
//   - The fitting is a simultaneous fit of a,pi,K,pi,u w/ the constraint that
//    a = pi+K+p+u.
//   - Output is, inter alia, in the shape of TGraphErrors indexed by Pbin:
//     - TGraphErrors named "particle_p/m_ID_Tbin" is Efficiency (ID==particle)
//      or misidentification probability (ID!=particle), where particle =
//      (pi||K||p), "p/m" = "+/-" is the charge of pE and ID=a,pi,K,p,u (where
//      'u' means 'unknown'). (Nota bene that here "p/m" concerns pE, contrary
//      to what we have in "chan[]", where it concerns pS.)
//     - Sub-TDirectory "Stats" contains TGraphErrors of S/B and statistical
//      significance (S/sqrt(S+B)).
// - "read_options" reads options from options file (D="options_fit.dat") for
//  both "plots" and "fit".

// TODO:
// - K signal selection: trade-off pi-Veto/K-ID w/ strict pi-Veto, i.e.
//  loose pi-LH cut (e.g. pi-LH not the largest, w/o any margin factor).
// - p+/- signal selection: try to loosen the requirement on h-/+ from pi-ID
//  (see "id_(p|m)==0" in "get_input_data_t1"), which restricts h-/+ to w/in
//  the pi/K-separation range (i.e. something like P<50 GeV), to (pi|K)-ID.
// - Several fit attempts: to investigate systematics of fit unstability
// - Apply event cuts:
// + BMS,
// + Extra track in exclusivity selection

const int start = 0;
const int stop = 6;

const double M_pi   = 0.13957018, M2_pi = M_pi*M_pi;
const double M_p    = 0.9382723,  M2_p  = M_p* M_p;
const double M_K    = 0.493677;
const double M_K0   = 0.497672;
const double M_Lam  = 1.115684;
const double M_phi  = 1.019456;

#define CSEVENTDATA 3

/**********************************************************************/
void usage() {
  printf(" * fit_table: RICH Table with pi,K,p via fit of Lambda,K0,phi\n");
  printf("Usage: fit_table [-f <optFile>] [-v] <mode>\n");
  printf("  <mode> = plots: Create the histograms for fitting.\n");
  printf("  <mode> = fit  : Do the fit and produce the table.\n");
  printf("  <mode> = test : Read options file and exit\n");
  printf("  -f: <optFile> specified on command line.\n");
  printf("  -h: Print this message and exit.\n");
  printf("  -v: Verbose.\n");
  printf("Default options file = \"./options_fit.dat\"\n");
  exit(1);
}
int main(int argc, char *argv[]){

  //   ********** PARSE COMMAND LINE
  verbose = 0; string optFile("options_fit.dat"); bool badCommandLine = false;
  int iarg = 1; while (iarg<argc && argv[iarg][0]=='-') {
    if      (string(argv[iarg])=="-f") {
      if (++iarg<argc) optFile = string(argv[iarg]);
      else badCommandLine = true;
    }
    else if (string(argv[iarg])=="-h") usage();
    else if (string(argv[iarg])=="-v") verbose = 1;
    iarg++;
  }
  if (!badCommandLine) badCommandLine = argc!=1+iarg || argv[iarg][0]=='-';
  string mode; if (!badCommandLine) {
    mode = string(argv[iarg]);
    badCommandLine = mode!="fit" && mode!="plots" && mode!="test";
  }
  if (badCommandLine) {
    cerr << "** fit_table: Ill formed command line: \"" << argv[0];
    for (int i = 1; i<argc; i++) cerr << " " << argv[i];
    cerr << "\"\n\n";
    usage();
  }

  if (!read_options(optFile)) return 1;    // ***** READ OPTIONS

  if (verbose)
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  else
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gErrorIgnoreLevel = kWarning;

  if      (mode=="plots") {             // ***** plots
    create_hist();
    ev = new CSEventData;
    hadrons = new vector<CSHadronData>;
    resonances = new vector<CSResonanceData>;
    for(int i=0; i<data_nb; i++) {
      if(get_inputFile(i+data_ff_nb)) {
	if(analysis=="K0L") get_input_data_t1(); // t1: K0 Lambda
	if(analysis=="phi") get_input_data_t2(); // t2: phi, both incl. and excl.
	input->Close();
      }
    }
    write_hist();
    return 0;
  }
  else if (mode=="fit") {               // ***** fit
    initCounts();
    get_plots();
    if (fit_type=="all") {
      fit_table_K0(0);     fit_table_K0(1);     // pi-, pi+
      fit_table_phi(0);    fit_table_phi(1);    // k-,  k+  from incl. phi 
      fit_table_Lambda(0); fit_table_Lambda(1); // p-,  p+
    }
    else if (fit_type=="pi") {
      fit_table_K0(0);     fit_table_K0(1);     // pi-, pi+
    }
    else if (fit_type=="k") {
      fit_table_phi(0);    fit_table_phi(1);    // k-,  k+  from incl. phi
    }
    else if (fit_type=="p") {
      fit_table_Lambda(0); fit_table_Lambda(1); // p-,  p+
    }
    print_table();
    return 0;
  }
  else if (mode=="test") return 0;      // ***** test
  return 0;
}

/**********************************************************************/
bool get_inputFile(int pi)
{
  string fileString = Form("%s/%s-%d.root",data_file.c_str(),data_template.c_str(),pi);
  const char *fileName = fileString.c_str();
  input = TFile::Open(fileName);
  if(input) {
    if (verbose) printf("=> \"%s\"\n",fileName);
    return 1;
  }
  else
    return 0;
}

/**********************************************************************/
bool read_options(string optFile) {
  //                  ********** INITIALISE LIKELIHOOD CUTS
  // ***** lh_cut[i][j]
  //  See (in particular to understand the above/below-threshold switch for p):
  // http://wwwcompass.cern.ch/compass/results/2015/april_kaon_mult/Multiplicities_note.pdf
  //  lh_cut[i][j] = Cut on LHk/LHj to select k = pi,K,bkg,p
  // i = 0,1,2,3 = pi,K,pbthr(below-threshold),pathr(above-threshold) and
  // i = 4 has a special meaning of upper cut on LH_(pi|K)(+|-)/LH_bg for
  //  below-threshold pID.
  // j = 0,1,2,3,4,5 = pi,K,p,e,mu,bkg
  // [10]: LH_K_pi,                          [12]: LH_K_p  [15]: LH_K_bg
  // [20]: LH_bthr_p_pi    [21]: LH_bthr_p_K
  // [30]: LH_athr_p_pi    [31]: LH_athr_p_K               [35]: LH_athr_p_bg
  // [40]: LH_bthr_p_pi_bg [41]: LH_bthr_p_K_bg
  for(int i = 0; i<6; i++) for(int j = 0; j<5; j++) lh_cut[j][i] = -1;
  const int ID_list[5] = {0,1,5,2,3};
  for (int i = 0; i<5; i++) id_lst[i] = ID_list[i];

  ifstream stream; string var1, var2; stringstream nn;

  map< string ,pair<int,int> > opt;
  // 	opt["LH_pi_pi:"] = make_pair(0,0);
  opt["LH_pi_K:"]  = make_pair(0,1);
  opt["LH_pi_p:"]  = make_pair(0,2);
  opt["LH_pi_e:"]  = make_pair(0,3);
  opt["LH_pi_mu:"] = make_pair(0,4);
  opt["LH_pi_bg:"] = make_pair(0,5);

  opt["LH_K_pi:"] = make_pair(1,0);
  // 	opt["LH_K_K:"]  = make_pair(1,1);
  opt["LH_K_p:"]  = make_pair(1,2);
  opt["LH_K_e:"]  = make_pair(1,3);
  opt["LH_K_mu:"] = make_pair(1,4);
  opt["LH_K_bg:"] = make_pair(1,5);

  opt["LH_bthr_p_pi:"] = make_pair(2,0);
  opt["LH_bthr_p_K:"]  = make_pair(2,1);
  opt["LH_bthr_p_p:"]  = make_pair(2,2);
  opt["LH_bthr_p_e:"]  = make_pair(2,3);
  opt["LH_bthr_p_mu:"] = make_pair(2,4);
  // 	opt["LH_bthr_p_bg:"] = make_pair(2,5);

  opt["LH_athr_p_pi:"] = make_pair(3,0);
  opt["LH_athr_p_K:"]  = make_pair(3,1);
  // 	opt["LH_athr_p_p:"]  = make_pair(3,2);
  opt["LH_athr_p_e:"]  = make_pair(3,3);
  opt["LH_athr_p_mu:"] = make_pair(3,4);
  opt["LH_athr_p_bg:"] = make_pair(3,5);

  opt["LH_bthr_p_pi_bg:"] = make_pair(4,0);
  opt["LH_bthr_p_K_bg:"]  = make_pair(4,1);
  opt["LH_bthr_m_pi_bg:"] = make_pair(4,2);
  opt["LH_bthr_m_K_bg:"]  = make_pair(4,3);

  stream.open(optFile.c_str(), ifstream::in);
  if (stream.fail()) {
    cerr << "Error opening <optFile> \""<<optFile<<"\" => Exiting...\n";
    return false;
  }

  map<string, pair<int,int> >::iterator it;

  string line;

  if(!stream.good()) { cerr << "No option file" << endl; std::exit(0);}
  while(true){
    if (stream.eof()) break;
    getline(stream,line);
    if( line[0] != '#' && !line.empty() ) {

      nn.str("");
      nn.clear();
      nn.str(line);
      nn >> var1 >> var2;
      cout << var1 << " " << var2 << endl;
      if (var1 == "analysis:")		analysis = var2;
      if (var1 == "data_file:")		data_file = var2;
      if (var1 == "data_template:")		data_template = var2;
      if (var1 == "data_firstfile_nb:")		data_ff_nb = stoi(var2.c_str());
      if (var1 == "data_lastfile_nb:")		data_lf_nb = stoi(var2.c_str());
      if (var1 == "fit_type:")		fit_type = var2;
      if (var1 == "hist_file_K0:")	hist_file_K0 = var2;
      if (var1 == "hist_file_iphi:")	hist_file_iphi = var2;
      if (var1 == "hist_file_ephi:")	hist_file_ephi = var2;
      if (var1 == "hist_file_Lam:")	hist_file_Lam = var2;
      if (var1 == "out_file:")		out_file = var2;
      if (var1 == "line_width:")		stringstream ( var2 ) >> lw;
      if (var1 == "remove_richpipe:"){if(var2=="true") rpipe = true; else rpipe = false;}
      if (var1 == "max_retry:")		stringstream ( var2 ) >> retry;
      if (var1 == "thr_diff:")		stringstream ( var2 ) >> thr_diff;
      if (var1 == "minuit_improve:")	{if(var2=="true") use_improve = true; else use_improve = false;}
      if (var1 == "minuit_hesse:")	{if(var2=="true") use_hesse = true; else use_hesse = false;}
      if (var1 == "minuit_minos:")	{if(var2=="true") use_minos = true; else use_minos = false;}
      if (var1 == "sidebins:")		{if(var2=="true") use_sidebins = true; else use_sidebins = false;}
      else{
	it = opt.find( var1 );
	if(it != opt.end()){
	  stringstream ( var2 ) >> lh_cut[(*it).second.first][(*it).second.second];
	}
      }
    }
  }
  data_nb = data_lf_nb-data_ff_nb+1;
  // cout << rpipe << endl;
  stream.close();

#define DEBUG_PID
#ifdef DEBUG_PID
  // - 0 - pion
  // - 1 - kaon
  // - 2 - proton
  // - 3 - electron
  // - 4 - muon
  // - 5 - background
    printf("=====================\nlh_cut[i][j]:\n");
    for (int i = 0; i<5; i++) {
      for (int j = 0; j<6; j++) {
	printf("[%d][%d] %5.2f ",i,j,lh_cut[i][j]);
      }
      printf("\n");
    }
#endif

  return true;
}

/**********************************************************************/
void create_hist(){
  stringstream nn;

  //const string chan[8] = {"K0_pip","K0_pim","phi_kp","phi_km","Lambda_pip","Lambda_pim","ephi_kp","ephi_km"};
  const char  *tags[8] = {"K0+",   "K0-",   "#phi+", "#phi-", "#Lambda+",  "#Lambda-",  "#phi+",  "#phi-"};
  const string id[5]   = {"a","pi","K","p","u"};
  // 	const Int_t Nbins[6]   = {120,  120,    33,    33,    70,   70};
  // 	const Double_t min[6]  = {0.44, 0.44, 0.98,  0.98,  1.09, 1.09};
  // 	const Double_t max[6]  = {0.56, 0.56, 1.046, 1.046, 1.16, 1.16};
  // 	const Int_t Nbins[6]   = {120,  120,     25,    25,   70,   70};
  // 	const Double_t min[6]  = {0.44, 0.44,   1.0,   1.0, 1.09, 1.09};
  // 	const Double_t max[6]  = {0.56, 0.56, 1.042, 1.042, 1.16, 1.16};
  const Int_t Nbins[8]   = { 120,  120,    30,    30,   70,   70,    30,    30};
  const Double_t min[8]  = {0.44, 0.44, 0.995, 0.995,  1.1,  1.1, 0.995, 0.995};
  const Double_t max[8]  = {0.56, 0.56, 1.042, 1.042, 1.13, 1.13, 1.042, 1.042};
  // 	const Double_t min[6]  = {0.44,0.44,0.98,0.98,1.10,1.10};
  // 	const Double_t max[6]  = {0.56,0.56,1.12,1.12,1.13,1.13};

  // ********** KINEMATICAL CUTS
  // ***** K0L: 0: K0, 1: Lambda
  /*
  DdD_cuts[0] = 2;      DdD_cuts[1] = 2;  //99999
  cth_cuts[0] = .99999; cth_cuts[1] = .99999;
  pT_cuts[0] =   .020;  pT_cuts[1] =  .020;
  */
  DdD_cuts[0] = 4;      DdD_cuts[1] = 4;  // tighter
  cth_cuts[0] = .99999; cth_cuts[1] = .99999;
  pT_cuts[0] =   .020;  pT_cuts[1] =  .020;
  // ***** phi
  dE_cuts[0] = 2.5; dE_cuts[1] = 2.5;
  pT_cuts[2] = .020; pT_cuts[3] = .020;
  // Avoid incl./excl. double counting
  if (pT_cuts[3]>pT_cuts[2]+.001) {
    printf("** fit_table: Overlapping pT cuts for Excl./Incl. phi\n");
    abort();
  }
  bookKineHistos();

  char hT[] = "#Lambda+ pi- 17<P<50 0.00<#theta<.004"; size_t sT = strlen(hT)+1;
  for(int i = 0; i<8; i++){      // K0 iphi Lambda ephi * h+/-ID-of-counterpart
    for(int j = 0; j<5; j++){      // All pi K p unID'd
      for(int p = 0; p<Np;p++){
	for(int t = 0; t<Nt; t++){
	  nn.str("");
	  nn.clear();
	  nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
	  snprintf(hT,sT,"%s %s%c %.0f<P<%.0f %.2f<#theta<%.2f",
		   tags[i],id[j].c_str(),(i%2)?'+':'-',
		   p_bins[p],p_bins[p+1],t_bins[t],t_bins[t+1]);
	  h[i][j][p][t] = new TH1D(nn.str().c_str(),hT,Nbins[i],min[i],max[i]);
	  h[i][j][p][t]->Sumw2();

	  nn.str("");
	  nn.clear();
	  nn << "am_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
	  h2[i][j][p][t] = new TH2D(nn.str().c_str(),"",100,-1.,1.,80,0.,0.4);

	  nn.str("");
	  nn.clear();
	  nn << "all" << (i%2?'+':'-') << ": " << p_bins[p] << " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
	  h2[i][j][p][t]->SetTitle(nn.str().c_str());
	  h2[i][j][p][t]->GetXaxis()->SetTitle("#alpha");
	  h2[i][j][p][t]->GetYaxis()->SetTitle("p_{t}");
	}
      }
    }
  }
}

/**********************************************************************/
void write_hist(){
  stringstream nn;
  
  //const string chan[8] = {"K0_pip","K0_pim","phi_kp","phi_km","Lambda_pip","Lambda_pim","ephi_kp","ephi_km"};
  const string id[5]   = {"a","pi","K","p","u"};


  for(int i =0; i<6; i++){

    cout << std::left
	 << setw(3) << i
	 << setw(15) << h[i][0][0][0]->GetEntries()
	 << setw(15) << h[i][1][0][0]->GetEntries()
	 << setw(15) << h[i][2][0][0]->GetEntries()
	 << setw(15) << h[i][3][0][0]->GetEntries()
	 << setw(15) << h[i][4][0][0]->GetEntries()
	 << setw(15) << h[i][0][0][0]->GetEntries() - h[i][1][0][0]->GetEntries() - h[i][2][0][0]->GetEntries() - h[i][3][0][0]->GetEntries() - h[i][4][0][0]->GetEntries()  << endl;

  }

  if(analysis=="K0L") {
    TFile* output = new TFile(hist_file_K0.c_str(),"RECREATE");
    for(int i = 0; i<2; i++) { // K0 +/-
      for(int j = 0; j<5; j++) { // a pi k p u
	for(int p = 0; p<Np;p++){
	  for(int t = 0; t<Nt; t++){
	    ((TH1D*)h[i][j][p][t])->Write();
	    ((TH2D*)h2[i][j][p][t])->Write();
	  }
	}
      }
    }
    writeKineHistos("K0");
    output->Close();
    delete output;

    output = new TFile(hist_file_Lam.c_str(),"RECREATE");
    for(int i = 4; i<6; i++) { // Lamda +/-
      for(int j = 0; j<5; j++) { // a pi k p u
	for(int p = 0; p<Np;p++){
	  for(int t = 0; t<Nt; t++){
	    ((TH1D*)h[i][j][p][t])->Write();
	    ((TH2D*)h2[i][j][p][t])->Write();
	  }
	}
      }
    }
    writeKineHistos("Lambda");
    output->Close();
  }
  if(analysis=="phi") {
    TFile* output = new TFile(hist_file_iphi.c_str(),"RECREATE");
    for(int i = 2; i<4; i++) { // Incl. phi +/-
      for(int j = 0; j<5; j++) { // a pi k p u
	for(int p = 0; p<Np;p++){
	  for(int t = 0; t<Nt; t++){
	    ((TH1D*)h[i][j][p][t])->Write();
	    ((TH2D*)h2[i][j][p][t])->Write();
	  }
	}
      }
    }
    writeKineHistos("Iphi");
    output->Close();
    delete output;

    output = new TFile(hist_file_ephi.c_str(),"RECREATE");
    for(int i = 6; i<8; i++) { // Excl. phi +/-
      for(int j = 0; j<5; j++) { // a pi k p u
	for(int p = 0; p<Np;p++){
	  for(int t = 0; t<Nt; t++){
	    ((TH1D*)h[i][j][p][t])->Write();
	    ((TH2D*)h2[i][j][p][t])->Write();
	  }
	}
      }
    }
    writeKineHistos("Ephi");
    output->Close();
  }

}
/**********************************************************************/
int getPID(double PR,     // Momentum @ RICH
	   double *p_lh,  // Array of LikeliHoods
	   double pi_thr, double p_thr,
	   int charge)
{
  // Returns id = 0:pi,1:K,2:p,3:e,5:bkg
  // (Notes: In its state of 2021/05, the method is the result of some dirty
  // adjustement, via "id_lst". As a consequence, it never returns id = 5 (i.e.
  // 'background'). Instead 3 (i.e. 'e') is returned where 5 would be expected,
  // e.g. when LH_bkg is the largest. This has no real impact on the working of
  // "fit_table".)
  int id = 3;  // Init. w/ 3 (see supra).
  // ***** p(m|p)_lh[j]
  // - 0 - pion
  // - 1 - kaon
  // - 2 - proton
  // - 3 - electron
  // - 4 - muon
  // - 5 - background  
  // ***** lh_cut[k][j]:   LH_k/LH_j > lh_cut[k][j]
  // k = pi,K,bkg,p = id_lst[i], i = 0,1,2,3 = pi,K,pbthr,pathr
  //                       [01]: LH_pi_K,    [02]: LH_pi_p [05]: LH_pi_bg
  // [10]: LH_K_pi,                          [12]: LH_K_p  [15]: LH_K_bg
  // [20]: LH_bthr_p_pi    [21]: LH_bthr_p_K
  // [30]: LH_athr_p_pi    [31]: LH_athr_p_K               [35]: LH_athr_p_bg
  // [40]: LH_bthr_p_pi_bg [41]: LH_bthr_p_K_bg

  if (PR>pi_thr){
    for (int i = 0; i<4; i++) { // i = pi,K,???,p-above-threshold (??? dirty trick 
      if(lh_cut[i][0]==-1 || p_lh[id_lst[i]]>lh_cut[i][0] * p_lh[0]){
	if(lh_cut[i][1]==-1 || p_lh[id_lst[i]]>lh_cut[i][1] * p_lh[1]){
	  if(lh_cut[i][2]==-1 || p_lh[id_lst[i]]>lh_cut[i][2] * p_lh[2]){
	    if(lh_cut[i][3]==-1 || p_lh[id_lst[i]]>lh_cut[i][3] * p_lh[3]){
	      if(lh_cut[i][4]==-1 || p_lh[id_lst[i]]>lh_cut[i][4] * p_lh[4]){
		if(lh_cut[i][5]==-1 || p_lh[id_lst[i]]>lh_cut[i][5] * p_lh[5]){
		  if(i != 2 || (PR<=p_thr+thr_diff && id !=0 && id !=1)){
		    if(i != 3 || PR>=p_thr-thr_diff){
		      id = id_lst[i];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //if (PR<=p_thr+thr_diff && p_lh[0]==-1 && p_lh[1]==-1 && p_lh[2]==-1 && p_lh[3]==-1 && p_lh[4]==-1 && p_lh[5]==-1) id = 5;
  if (p_lh[0]==-1 && p_lh[1]==-1 && p_lh[2]==-1 && p_lh[3]==-1 && p_lh[4]==-1 && p_lh[5]==-1) id = 3;
  // Below threshold
  if (charge==1) {
    if (PR<=p_thr+thr_diff && ((p_lh[0]==0 && p_lh[1]==0 && p_lh[2]==0
				&& p_lh[3]==0 && p_lh[4]==0 && p_lh[5]==0)
			       || (p_lh[0]<lh_cut[4][0] * p_lh[5]
				   && p_lh[1]<lh_cut[4][1] * p_lh[5])) ) id = 2;
  }
  else {
    if (PR<=p_thr+thr_diff && ((p_lh[0]==0 && p_lh[1]==0 && p_lh[2]==0
				&& p_lh[3]==0 && p_lh[4]==0 && p_lh[5]==0)
			       || (p_lh[0]<lh_cut[4][2] * p_lh[5]
				   && p_lh[1]<lh_cut[4][3] * p_lh[5])) ) id = 2;
  }
  return id;
}
/**********************************************************************/
void get_input_data_t1(){
  TTree *tree = (TTree*)input->Get("CSEvtTree");
  if (!tree) {
    printf("** get_input_data_t1: No \"CSEvtTree\" TTree in TFile \"%s\"\n",
	   gFile->GetName());
    exit(1);
  }
  tree->SetBranchAddress("CSEvt",&ev);
  tree->SetBranchAddress("Hs",&hadrons);
  tree->SetBranchAddress("Rs",&resonances);
  Long64_t nentries1 = tree->GetEntries();

  stringstream cut;

  double pp_lh[7],pm_lh[7];
  for(int i = 0; i<7; i++) {
    pm_lh[i] = -1;
    pp_lh[i] = -1;
  }

  cout << nentries1  << endl;

  // ********** PREPARE SELECTIONS
  // ***** RUN DEPENDENT VARIABLES
  int prvRun = 0;
  static double pi_thr, k_thr, p_thr;
  // ***** EVENT DEPENDENT
  int prvEvt = 0, prvIhp = -1, prvIhm = -1;
  // ***** K0/Lambda
  unsigned short K0Required =     0x3f;
  unsigned short LambdaRequired = 0x3f;
  
  for (Long64_t jentry=0; jentry<nentries1;jentry++) {

    tree->GetEntry(jentry);

    if (ev->runNo!=prvRun) {     // ***** RUN DEPENDENT VARIABLES
      prvRun = ev->runNo;
      pi_thr = ev->piThr; k_thr = pi_thr*M_K/M_pi; p_thr = pi_thr*M_p/M_pi;
    }
    double Xp = ev->Xp, Yp = ev->Yp, Zp = ev->Zp;
    int evt = ev->evtNo;
    
    vector<CSHadronData> &hdrns= *hadrons;
    vector<CSResonanceData> &vRes = *resonances; int nRes = vRes.size();
    for (int iRes = 0; iRes<nRes; iRes++) {
      int p_bin_m = -1;
      int p_bin_p = -1;
      int t_bin_m = -1;
      int t_bin_p = -1;
      CSResonanceData &res = vRes[iRes];
      unsigned short K0Pat = res.K0Pat, LambdaPat = res.LambdaPat;
      if ((K0Pat&K0Required)!=K0Required &&
	  (LambdaPat&LambdaRequired)!=LambdaRequired) continue;
      short ihp = res.h1, ihm = res.h2;
      CSHadronData &hp = hdrns[ihp], &hm = hdrns[ihm];

      // ***** P AND theta BINNING
      // P and theta are taken @ RICH
      double PRp = hp.qP, PRm = hm.qP;
      if (PRp<0 || PRm>0) { // PRp>0 PRm<0 by construction: double check it.
	printf("** fit_table: Evt %d#%d, CsRes %d: PRp,PRm = %.2f,%.2f\n",
	       ev->runNo,ev->evtNo,iRes,PRp,PRm);
	abort();
      }
      PRp = fabs(PRp); PRm = fabs(PRm);
      float tgXR, tgYR;
      tgXR = hp.tgXR; tgYR = hp.tgYR;
      double thRp = acos(1/sqrt(1.+ tgXR*tgXR + tgYR*tgYR));
      tgXR = hm.tgXR; tgYR = hm.tgYR;
      double thRm = acos(1/sqrt(1.+ tgXR*tgXR + tgYR*tgYR));
      // ***** (P,thR) => BINNING
      for (int i = 0 ; i<Np; i++) {
	if (PRm>=p_bins[i] && PRm<p_bins[i+1]) {
	  //if(PRp> pi_thr){ // granted by (K0|Lambda)Pat &= 0x10
	  p_bin_m = i;
	}
	if (PRp>=p_bins[i] && PRp<p_bins[i+1]) {
	  //if(PRm> pi_thr){ // granted by (K0|Lambda)Pat &= 0x20
	  p_bin_p = i;
	}
      }
      for (int i = 0 ; i<Nt; i++) {
	if (thRm>=t_bins[i] && thRm<t_bins[i+1]) t_bin_m = i;
	if (thRp>=t_bins[i] && thRp<t_bins[i+1]) t_bin_p = i;
      }
      if (p_bin_m==-1 && p_bin_p==-1 && t_bin_m==-1 && t_bin_p==-1)
	continue;   // ***** BOTH m AND p OUT OF SCOPE

      // ***** REJECT RICH PIPE: could be revisited: why both + and -?
      bool pipe = false;
      float pp_x = hp.XR, pp_y = hp.YR, pm_x = hm.XR, pm_y = hm.YR;
      if(rpipe){
	if(pp_x*pp_x + pp_y*pp_y >=25. && pm_x*pm_x + pm_y*pm_y >=25.) pipe = true;
      }
      if(!( pipe || !rpipe)) continue;

      double alpha = res.alpha, pT = res.pT;
      if (evt!=prvEvt || (prvIhp!=ihp && prvIhm!=ihm)) {// Avoid double counting
	prvEvt = evt; prvIhp = ihp; prvIhm = ihm;
	am_all->Fill(alpha,pT);  // ***** OVERALL ARMENTEROS/KINEMATICS
	Z_all->Fill(Zp); XY_all->Fill(Xp,Yp);
      }

      // ***** V0 CUTS
      double D = res.D, dD = res.dD, cth = res.cth;
      if (K0Pat) {
	if (D/dD<DdD_cuts[0]) continue;
	if (cth<cth_cuts[0]) continue;
	if (pT<pT_cuts[0]) continue;
      }
      else {
	if (D/dD<DdD_cuts[1]) continue;
	if (cth<cth_cuts[1]) continue;
	if (pT<pT_cuts[1]) continue;
      }
      if (K0Pat) {
	if (dD) DdD_K0->Fill(D/dD);
	else // Consistency check: "dD" is !=0 by construction
	  printf("** get_input_data_t1: Evt %d#%d, CsRes %d: D = %.2f, dD =0\n",
		    ev->runNo,ev->evtNo,iRes,D);
	cth_K0->Fill(cth); pT_K0->Fill(pT);
      }
      else {
	if (dD) DdD_L->Fill(D/dD);
	else // Consistency check: "dD" is !=0 by construction
	  printf("** get_input_data_t1: Evt %d#%d, CsRes %d: D = %.2f, dD =0\n",
		    ev->runNo,ev->evtNo,iRes,D);
	cth_L->Fill(cth);  pT_L->Fill(pT);
      }

      // ********** PID
      for (int i = 0; i<6; i++) { pp_lh[i] = hp.LH[i]; pm_lh[i] = hm.LH[i]; }
      int id_p = getPID(PRp,pp_lh,pi_thr,p_thr, 1);
      int id_m = getPID(PRm,pm_lh,pi_thr,p_thr,-1);

      if( (pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) ||
	  (pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1)) continue;

      // ***** ARMENTEROS FOR K0 W/ pi-ID
      if (K0Pat) {
	if (id_m==0) am_K0p->Fill(alpha,pT);
	if (id_p==0) am_K0m->Fill(alpha,pT);
      }

      // ***** REJECT K0/Lambda REFLECTION
      int Lambda;
      if      (LambdaPat&0x800) Lambda = -1;
      else if (LambdaPat)       Lambda =  1;
      else                      Lambda =  0;
      // ***** MASSES: CONVERT M_RESONANCE -> M_REFLEXION
      double Pp = sqrt(hp.Px*hp.Px+hp.Py*hp.Py+hp.Pz*hp.Pz);
      double Pm = sqrt(hm.Px*hm.Px+hm.Py*hm.Py+hm.Pz*hm.Pz);
      double P2m = Pm*Pm, P2p = Pp*Pp, m2 = res.m*res.m;
      static double mpipi, mppi, mpip;
      double Epipi =  sqrt(P2p+M2_pi )+sqrt(P2m+M2_pi );
      if      (Lambda==1) {
	double Eppi = sqrt(P2p+M2_p )+sqrt(P2m+M2_pi);
	mpipi = sqrt(m2+Epipi*Epipi-Eppi*Eppi);
      }
      else if (Lambda==-1) {
	double Epip = sqrt(P2p+M2_pi)+sqrt(P2m+M2_p );
	mpipi = sqrt(m2+Epipi*Epipi-Epip*Epip);
      }
      else {
	double Eppi = sqrt(P2p+M2_p )+sqrt(P2m+M2_pi);
	double Epip = sqrt(P2p+M2_pi)+sqrt(P2m+M2_p );
	mpip =  sqrt(m2-Epipi*Epipi+Epip*Epip);
	mppi =  sqrt(m2-Epipi*Epipi+Eppi*Eppi);
      }
      bool fillHisto;
      if (Lambda) {
	fillHisto = fabs(mpipi-M_K0)>0.02;
      }
      else {
	if (alpha>0) fillHisto = fabs(mppi-M_Lam)>0.01;
	else         fillHisto = fabs(mpip-M_Lam)>0.01;
      }
      if (!fillHisto) continue;

      if (K0Pat) {  // ***** ARMENTEROS/KINEMATICS AFTER CUTS
	am_K0->Fill(alpha,pT); Z_K0->Fill(Zp); XY_K0->Fill(Xp,Yp);
#if CSEVENTDATA == 3
	Tr_K0->Fill(ev->nOuts);
#endif
	Rb_K0->Fill(ev->nTrksRIb); Rt_K0->Fill(ev->nTrksRIt);
      }
      else {
	am_L->Fill(alpha,pT);  Z_L->Fill(Zp);  XY_L->Fill(Xp,Yp);
#if CSEVENTDATA == 3
	Tr_L->Fill(ev->nOuts);
#endif
	Rb_L->Fill(ev->nTrksRIb);  Rt_L->Fill(ev->nTrksRIt);
      }

      if (p_bin_m!=-1 && t_bin_m!=-1  &&  // ***** FILLING NEGATIVE pE- *****
	  id_p==0 /* ID-BASED SPECTATOR pS+ SELECTION */) {
	if      (Lambda==0) {
	  h[0][0][p_bin_m][t_bin_m]->Fill(res.m);
	  h2[0][0][p_bin_m][t_bin_m]->Fill(alpha,pT);
	}
	else if (Lambda==-1) {
	  h[4][0][p_bin_m][t_bin_m]->Fill(res.m);
	  h2[4][0][p_bin_m][t_bin_m]->Fill(alpha,pT);
	}
	for(int i = 0; i<5; i++){
	  if(id_m == id_lst[i]){
	    if      (Lambda==0) {
	      if(id_m!=5){
		h[0][id_m+1][p_bin_m][t_bin_m]->Fill(res.m);
		h2[0][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pT);
	      } else{
		h[0][3][p_bin_m][t_bin_m]->Fill(res.m);
		h2[0][3][p_bin_m][t_bin_m]->Fill(alpha,pT);
	      }
	    }
	    else if (Lambda==-1) {
	      if(id_m!=5){
		h[4][id_m+1][p_bin_m][t_bin_m]->Fill(res.m);
		h2[4][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pT);
	      } else{
		h[4][3][p_bin_m][t_bin_m]->Fill(res.m);
		h2[4][3][p_bin_m][t_bin_m]->Fill(alpha,pT);
	      }
	    }
	  }
	}
      } // End filling negative pE-
      if (p_bin_p!=-1 && t_bin_p!=-1 &&  // ***** FILLING POSITIVE pE+ *****
	  id_m==0 /* ID-BASED SPECTATOR pS- SELECTION = pi-ID */) {
	if      (Lambda==0) {
	  h[1][0][p_bin_p][t_bin_p]->Fill(res.m);
	  h2[1][0][p_bin_p][t_bin_p]->Fill(alpha,pT);
	}
	else if (Lambda==1) {
	  h[5][0][p_bin_p][t_bin_p]->Fill(res.m);
	  h2[5][0][p_bin_p][t_bin_p]->Fill(alpha,pT);
	}
	for(int i = 0; i<5; i++){
	  if(id_p == id_lst[i]){
	    if      (Lambda==0) {
	      if(id_p!=5){
		h[1][id_p+1][p_bin_p][t_bin_p]->Fill(res.m);
		h2[1][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pT);
	      } else{
		h[1][3][p_bin_p][t_bin_p]->Fill(res.m);
		h2[1][3][p_bin_p][t_bin_p]->Fill(alpha,pT);
	      }
	    }
	    else if (Lambda==1) {
	      if(id_p!=5){
		h[5][id_p+1][p_bin_p][t_bin_p]->Fill(res.m);
		h2[5][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pT);
	      } else{
		h[5][3][p_bin_p][t_bin_p]->Fill(res.m);
		h2[5][3][p_bin_p][t_bin_p]->Fill(alpha,pT);
	      }
	    }

	  }
	}
      } // End filling positive pE+
    } // End loop on resonances p0
  } // End loop on entries

  delete tree;
}
/**********************************************************************/
void get_input_data_t2(){
  TTree *tree2 = (TTree*)input->Get("CSEvtTree");
  if (!tree2) {
    printf("** get_input_data_t2: No \"CSEvtTree\" TTree in TFile \"%s\"\n",
	   gFile->GetName());
    exit(1);
  }
  tree2->SetBranchAddress("CSEvt",&ev);
  tree2->SetBranchAddress("Hs",&hadrons);
  tree2->SetBranchAddress("Rs",&resonances);
  Long64_t nentries2 = tree2->GetEntries();

  cout << nentries2  << endl;
  // cout << "tree 2" << endl;

  stringstream cut;

  double pp_lh[7],pm_lh[7];
  for(int i = 0; i<7; i++) { pm_lh[i] = -1; pp_lh[i] = -1; }

  //if(lv_phi->Mag() >  1.042 || lv_phi->Mag() < 0.995) continue;

  // ********** PREPARE SELECTIONS
  // ***** RUN DEPENDENT VARIABLES
  int prvRun = 0;
  static double pi_thr, k_thr, p_thr;
  // ***** phi
  // Incl.: Also: 0x1 (3 outs) to be rejected? 0x30: couldn't it be too strong?
  unsigned short IphiRequired = 0x436;
  // Excl.: Also 0x100 = No detached track? 0x200 = No ECAL?
  unsigned short EphiRequired = 0x37; // 0x8 (dEK w/in cuts applied later)

  for (Long64_t jentry=0; jentry<nentries2;jentry++) {
    tree2->GetEntry(jentry);

    if (ev->runNo!=prvRun) {     // ***** RUN DEPENDENT VARIABLES
      prvRun = ev->runNo;
      pi_thr = ev->piThr; k_thr = pi_thr*M_K/M_pi; p_thr = pi_thr*M_p/M_pi;
    }
    double Xp = ev->Xp, Yp = ev->Yp, Zp = ev->Zp;

    // ***** EXCLUSIVE phi: THREE CRITERIA: phiPat&0x1 !&0x400 dE<cut
    // MISSING ENERGY
    double dEK = ev->dEK;
    bool inclEMiss = dEK>dE_cuts[0];
    bool exclEMiss = fabs(dEK)<dE_cuts[1];

    /*
      // K* 4-vector
    TLorentzVector lv_pip, lv_pim, lv_ks1, lv_ks2;
    lv_pip.SetVectM(lv_kp->Vect(),m_pi);
    lv_pim.SetVectM(lv_km->Vect(),m_pi);
    lv_ks1 = lv_pim + *lv_kp;
    lv_ks2 = lv_pip + *lv_km;
    // cout << lv_pip.Mag() << " " << lv_pim.Mag() << " " << lv_ks.Mag() << endl;
    */
#define LOOSE_K_SELECTION
#ifdef LOOSE_K_SELECTION
    // Selection of K+/-: Don't require K-/+ID, but veto pion (and p as much a possible)
    // i) h-/+ not pi nor p (to which one should have added not e)
    // ii) Only far enough above K threshold (i.e. P > KThr*1.2): KID
    double piKMx = 50; // Max. P for pi/K resolution
    double pieMx = 40; // Max. P for e/pi resolution
#endif
    vector<CSHadronData> &hdrns= *hadrons;
    vector<CSResonanceData> &vRes = *resonances; int nRes = vRes.size();
    for (int iRes = 0; iRes<nRes; iRes++) {
      int p_bin_m = -1;
      int p_bin_p = -1;
      int t_bin_m = -1;
      int t_bin_p = -1;
      CSResonanceData &res = vRes[iRes];
      unsigned short phiPat = res.phiPat;
      // ***** INCL/EXCL SELEC>TION ...BUT FOR EMISS
      // (Note: possibly redundant...)
      bool isIncl = (phiPat&IphiRequired)==IphiRequired;
      bool isExcl = (phiPat&EphiRequired)==EphiRequired;
      if ((phiPat&0x1) && (phiPat&0x400)) {
	printf("** fit_table: Evt %d/%d, CsRes %d: phiPat 0x400 AND 0x1\n",
	       ev->runNo,ev->evtNo,iRes);
	abort();
      }
      int ie; if (isIncl) ie = 0;
      else    if (isExcl) ie = 1;
      else continue;
      short ihp = res.h1, ihm = res.h2;
      CSHadronData &hp = hdrns[ihp], &hm = hdrns[ihm];

      // ***** P AND theta BINNING
      // P and theta are taken @ RICH
      double PRp = hp.qP, PRm = hm.qP;
      if (PRp<0 || PRm>0) {// PRp>0 PRm<0: like-sign pairs may have been allowed
	printf(" * fit_table: Evt %d/%d, CsRes %d: PRp,PRm = %.2f,%.2f\n",
	       ev->runNo,ev->evtNo,iRes,PRp,PRm);
	continue;
      }
      PRp = fabs(PRp); PRm = fabs(PRm);
      float tgXR, tgYR;
      tgXR = hp.tgXR; tgYR = hp.tgYR;
      double thRp = acos(1/sqrt(1.+ tgXR*tgXR + tgYR*tgYR));
      tgXR = hm.tgXR; tgYR = hm.tgYR;
      double thRm = acos(1/sqrt(1.+ tgXR*tgXR + tgYR*tgYR));
      for(int i = 0 ; i<Np; i++){
	if(PRm >= p_bins[i] && PRm < p_bins[i+1]){
#ifdef LOOSE_K_SELECTION
	  // K- loose selection: require P+ > piThr instead of P+ > Kthr
	  if (PRp>pi_thr) p_bin_m = i; 
#else
	  if (PRp>k_thr)  p_bin_m = i;
#endif
	}
	if(PRp >= p_bins[i] && PRp < p_bins[i+1]){
#ifdef LOOSE_K_SELECTION
	  // K- loose selection: require P+ > piThr instead of P+ > Kthr
	  if (PRm>pi_thr) p_bin_p = i; 
#else
	  if (PRm>k_thr)  p_bin_p = i;
#endif
	}
      }
      for(int i = 0 ; i<Nt; i++){
	if(thRm >= t_bins[i] && thRm < t_bins[i+1]) t_bin_m = i;
	if(thRp >= t_bins[i] && thRp < t_bins[i+1]) t_bin_p = i;
      }
      if (p_bin_m==-1 && p_bin_p==-1 && t_bin_m==-1 && t_bin_p==-1)
	continue;   // ***** BOTH p AND m OUT OF SCOPE

      // ***** REJECT RICH PIPE: could be revisited: why both + and -?
      bool pipe = false;
      float pp_x2 = hp.XR, pp_y2 = hp.YR, pm_x2 = hm.XR, pm_y2 = hm.YR;
      if(rpipe){
	if(pp_x2*pp_x2 + pp_y2*pp_y2 >=25. && pm_x2*pm_x2 + pm_y2*pm_y2 >=25.) pipe = true;
      }
      if(!( pipe || !rpipe)) continue;      

      // ********** PID
      for (int i = 0; i<6; i++) { pp_lh[i] = hp.LH[i]; pm_lh[i] = hm.LH[i]; }
      int id_p = getPID(PRp,pp_lh,pi_thr,p_thr, 1);
      int id_m = getPID(PRm,pm_lh,pi_thr,p_thr,-1);

      if( (pp_lh[0] == -1 && pp_lh[1] == -1 && pp_lh[2] == -1 && pp_lh[3] == -1 && pp_lh[4] == -1 && pp_lh[5] == -1) || (pm_lh[0] == -1 && pm_lh[1] == -1 && pm_lh[2] == -1 && pm_lh[3] == -1 && pm_lh[4] == -1 && pm_lh[5] == -1)) continue;

      bool counterpartID = id_p==1 || id_m==1;
      
      double alpha = res.alpha, pT = res.pT;
      //              ***** OVERALL KINEMATICS (W/ COUNTERPART ID)
      bool EMissOK = ie==0 && inclEMiss || ie!=0 && exclEMiss;
      bool pTOK; if   (ie==0) pTOK = pT>pT_cuts[2];
      else                    pTOK = pT>pT_cuts[3];
      if (counterpartID) {
	if (pTOK) {
	  if (ie==0) dE_Incl->Fill(dEK);
	  else       dE_Excl->Fill(dEK);
	}
	if (EMissOK) {
	  if (ie==0) { pT_Incl->Fill(pT); am_Incl->Fill(alpha,pT); }
	  else       { pT_Excl->Fill(pT); am_Excl->Fill(alpha,pT); }
	}
      }

      if (!EMissOK) continue;                           // ***** EMiss CUT
      if (!pTOK) continue;                              // ***** pT CUT
      /*
      // Exclude K*
      if(fabs(lv_ks1.Mag()-0.89166)>0.05){
	if(fabs(lv_ks2.Mag()-0.89166)>0.05){
      */
      
      // ***** ARMENTEROS/KINEMATICS AFTER CUTS
      if (counterpartID) {
	if (ie==0) {
	  am_Iphi->Fill(alpha,pT); Z_Iphi->Fill(Zp); XY_Iphi->Fill(Xp,Yp);
	  pT_Iphi->Fill(pT); dE_Iphi->Fill(dEK); // To double-check kine. cuts
#if CSEVENTDATA == 3
	  Tr_Iphi->Fill(ev->nOuts);
#endif
	  Rb_Iphi->Fill(ev->nTrksRIb); Rt_Iphi->Fill(ev->nTrksRIt);
	}
	else {
	  am_Ephi->Fill(alpha,pT); Z_Ephi->Fill(Zp); XY_Ephi->Fill(Xp,Yp);
	  pT_Ephi->Fill(pT); dE_Ephi->Fill(dEK); // To double-check kine. cuts
#if CSEVENTDATA == 3
	  Tr_Ephi->Fill(ev->nOuts);
#endif
	  Rb_Ephi->Fill(ev->nTrksRIb); Rt_Ephi->Fill(ev->nTrksRIt);
	}
      }

      // ***** FILL HISTOS
      int IE =  ie ? 4 : 0;
      bool selectKm = p_bin_m!=-1 && t_bin_m!=-1;
      if (selectKm) {                  // ***** FILLING POSITIVE pE- = K- *****
	// ID-BASED SPECTATOR pS+ REJECTION
#ifdef LOOSE_K_SELECTION
	if (id_p==0 ||                   // pS = pi
	    PRp>p_thr*1.1 && id_p==2 ||  // pS = p above p-threshold
	    PRp>piKMx ||                 // pS above piK resolution
	    PRp<pieMx &&                 // pS = e w/in pie resolution
	    pp_lh[3]>2*pp_lh[5] && pp_lh[3]>1.5*pp_lh[1]) 
	  selectKm = false;
#else
	if (id_p!=1) selectKm = false;
#endif
      }
      if (selectKm) {
      // 				if(TMath::Abs(lv_ks1.Mag()-0.89166)>0.05){
      // 				if(lv_ks1.Mag()-0.89166<-0.05){
      //					if(TMath::Abs(lv_ks2.Mag()-0.89166)>0.03){
	h[2+IE][0][p_bin_m][t_bin_m]->Fill(res.m);
	h2[2+IE][0][p_bin_m][t_bin_m]->Fill(alpha,pT);
	for(int i = 0; i<5; i++){
	  if(id_m == id_lst[i]){
	    if(id_m!=5){
	      h[2+IE][id_m+1][p_bin_m][t_bin_m]->Fill(res.m);
	      h2[2+IE][id_m+1][p_bin_m][t_bin_m]->Fill(alpha,pT);
	    } else{
	      h[2+IE][3][p_bin_m][t_bin_m]->Fill(res.m);
	      h2[2+IE][3][p_bin_m][t_bin_m]->Fill(alpha,pT);
	    }
	  }
	}
	// 					}
	// 				}
      }
      bool selectKp = p_bin_p!=-1 && t_bin_p!=-1;
      if (selectKp) {                  // ***** FILLING POSITIVE pE+ = K+ *****
	// ID-BASED SPECTATOR pS- REJECTION
#ifdef LOOSE_K_SELECTION
	if (id_m==0 ||                   // pS = pi
	    PRm>p_thr*1.1 && id_m==2 ||  // pS = p above p-threshold
	    PRm>piKMx ||                 // pS above piK resolution
	    PRm<pieMx &&                 // pS = e w/in pie resolution
	    pm_lh[3]>2*pm_lh[5] && pm_lh[3]>1.5*pm_lh[1]) 
	  selectKp = false;
#else
	if (id_m!=1) selectKp = false;
#endif
      }
      if (selectKp) {
	//				if(TMath::Abs(lv_ks1.Mag()-0.89166)>0.03){
	// 					if(TMath::Abs(lv_ks2.Mag()-0.89166)>0.05){
	// 					if(lv_ks2.Mag()-0.89166<-0.05){
	h[3+IE][0][p_bin_p][t_bin_p]->Fill(res.m);
	if (p_bin_p==7 && t_bin_p==1 && .995<res.m && res.m<1.042 && ie==0) {
	  static FILE *fp = 0; static int nevts = 0;
	  if (!fp) fp = fopen("phip7_1.txt","w");
	  if (!fp) {
	    printf("No opening \"phip7_1.txt\"\n"); abort();
	  }
	  double integral = h[3][0][7][1]->Integral(1,30);
	  printf("%d#%9d %4d,%4.0f %.3f %.3f,%.3f %.3f,%.3f\n",
		 /**/  ev->runNo,ev->evtNo,++nevts,integral,res.m,PRp,thRp,PRm,thRm);
	  fprintf(fp,"%d#%9d %4d,%4.0f %.3f %.3f,%.3f %.3f,%.3f\n",
		  /**/ ev->runNo,ev->evtNo,  nevts,integral,res.m,PRp,thRp,PRm,thRm);
	}
	h2[3+IE][0][p_bin_p][t_bin_p]->Fill(alpha,pT);
	for(int i = 0; i<5; i++){
	  if(id_p == id_lst[i]){
	    if(id_p!=5){
	      h[3+IE][id_p+1][p_bin_p][t_bin_p]->Fill(res.m);
	      h2[3+IE][id_p+1][p_bin_p][t_bin_p]->Fill(alpha,pT);
	    } else{
	      h[3+IE][3][p_bin_p][t_bin_p]->Fill(res.m);
	      h2[3+IE][3][p_bin_p][t_bin_p]->Fill(alpha,pT);
	    }
	  }
	}
      }
    } // End loop on resonances p0
  } // End loop on entries

  delete tree2;

}
/**********************************************************************/
void initCounts()
{
  for(int i = 0; i<8; i++)    // Channels
    for(int j = 0; j<6; j++)    // a,pi,K,p,u  and background
      for(int t = 0; t< Nt; t++)
	for(int p = 0; p<Np; p++) N_id[i][j][p][t] = -1;
}
/**********************************************************************/
void get_plots(){
  //const string chan[8] = {"K0_pip","K0_pim","phi_kp","phi_km","Lambda_pip","Lambda_pim","ephi_kp","ephi_km"};
  const string id[5]   = {"a","pi","K","p","u"};
  stringstream nn;
  // cout << 1 << endl;
  input_K0 = new TFile(hist_file_K0.c_str(),"READ");
  if (!input_K0 || input_K0->IsZombie()) {
    printf("** fit_table: Error opening TFile \"%s\"\n",hist_file_K0.c_str());
    abort();
  }
  for(int i = 0; i<2; i++){      //6
    for(int j = 0; j<5; j++){  //4
      for(int p = 0; p<Np;p++){
	for(int t = 0; t<Nt; t++){
	  nn.str("");
	  nn.clear();
	  nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
	  h[i][j][p][t] = (TH1D*)input_K0->Get(nn.str().c_str());
	  if (!input_K0->Get(nn.str().c_str())) {
	    printf("** fit_table: TH1D \"%s\" does not exist\n",
		   nn.str().c_str());
	    abort();
	  }
	}
      }
    }
  }

  // cout << 2 << endl;
  input_iphi = new TFile(hist_file_iphi.c_str(),"READ");
  if (!input_iphi || input_iphi->IsZombie()) {
    printf("** fit_table: Error opening TFile \"%s\"\n",hist_file_iphi.c_str());
    abort();
  }
  for(int i = 2; i<4; i++){      //6
    for(int j = 0; j<5; j++){  //4
      for(int p = 0; p<Np;p++){
	for(int t = 0; t<Nt; t++){
	  nn.str("");
	  nn.clear();
	  nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
	  h[i][j][p][t] = (TH1D*)input_iphi->Get(nn.str().c_str());
	}
      }
    }
  }

  // cout << 3 << endl;
  input_Lam = new TFile(hist_file_Lam.c_str(),"READ");
  for(int i = 4; i<6; i++){      //6
    for(int j = 0; j<5; j++){  //4
      for(int p = 0; p<Np;p++){
	for(int t = 0; t<Nt; t++){
	  nn.str("");
	  nn.clear();
	  nn << "h_" << chan[i] << "_" << id[j] << "_" << p <<"_" <<t;
	  h[i][j][p][t] = (TH1D*)input_Lam->Get(nn.str().c_str());
	}
      }
    }
  }




}
/**********************************************************************/
// - "fit_table_(K0|phi|Lambda)":
//   - Fit model for simultaneous fit of "all,pi,K,p,u" histos.
//   - Enforcing "N_a_s" and "N_a_b" for "all" to be the sum of other four.
//   - Fit is iterated until:
//     - "RooFitResult::covQual" =3,
//     - "RooMinuit::migrad/improve" be OK,
//     - "RooFitResult" uncertainty on efficiency does not depart too much
//       from the binomial formula "sqrt(eff(1-eff)/N)".
// - "getEffDEff": compute efficiency/purities and uncertainties and
//  comparison w/ binomial.
// - "root2PDF": convert "RooPlot's" into PDFs.
// - "getSB": compute Signal and Background, store in "N|Bs[channel][p][t]".
/**********************************************************************/
int getEffDEff(const RooFitResult *r, int h2h, double *effs, double *dEffs) 
{
  // - Computes efficiency (h->h) and purities (h->!h) and their uncertainties,
  //  given input "RooFitResult".
  // - Compares efficiency uncertainty w/ binomial formula "sqrt(eff(1-ff)/N)",
  //  and returns a status word.
  // - Upon return, efficiency and purities are stored in arg. "eff/dEffs".
  int status = 0;
  int ip, jp, kp, idcs[4]; // pi, K, p, u
  const char *pNames[] = {"N_pi_s","N_k_s","N_p_s","N_u_s"};
  const RooArgList &pars = r->floatParsFinal();
  for (ip = 0; ip<4; ip++) idcs[ip] = pars.index(pNames[ip]);
  double nA; for (ip = 0, nA = 0; ip<4; ip++) {
    nA += dynamic_cast<const RooRealVar&>(pars[idcs[ip]]).getVal();
  }
  for (ip = 0; ip<4; ip++) {
    if (nA) {
      double nS = dynamic_cast<const RooRealVar&>(pars[idcs[ip]]).getVal();
      double dAdpi[4], dEff, eff = nS/nA; effs[ip] = eff;
      for (jp = 0; jp<4; jp++) {
	dAdpi[jp] = -nS/nA/nA; if (jp==ip) dAdpi[ip] += 1/nA;
      }
      for (jp = 0, dEff = 0; jp<4; jp++) {
	int jdx = idcs[jp];
	for (kp = jp+1; kp<4; kp++) {
	  int kdx = idcs[kp];
	  dEff += 2*dAdpi[jp]*dAdpi[kp]*r->covarianceMatrix()(jdx,kdx);
	}
	dEff += dAdpi[jp]*dAdpi[jp]*r->covarianceMatrix()(jdx,jdx);
      }
      dEffs[ip] = sqrt(dEff);
      double bDEff = sqrt(eff*(1-eff)/nA);
      if (ip==h2h && bDEff>.005 && dEff>1.5*bDEff) status = 1;
    }
    else {
      effs[ip]=dEffs[ip] = 0;
    }
  }
  return status;
}
void root2PDF(const char *outFName, int cc)
{
  // RooPlot's to PDF from ROOT file complete w/ all plots for current <cc>.
  char cN[] = "cpP99T99", mp[] = "mp";
  int lastT, lastP, t; for (t = 0, lastT=lastP = -1; t<Nt; t++) {
    for (int p = 0; p<Np; p++) {
      sprintf(cN,"c%cP%dT%d",mp[cc],p,t);
      TCanvas *c = (TCanvas*)gDirectory->Get(cN); if (c) {
	lastT = t; lastP = p;
      }
    }
  }
  int first; for (t = 0, first = 1; t<Nt; t++) {
    for (int p = 0; p<Np; p++) {
      sprintf(cN,"c%cP%dT%d",mp[cc],p,t);
      TCanvas *c = (TCanvas*)gDirectory->Get(cN); if (c) {
	if (first) {
	  string pdfOpen =  string(outFName)+string(".pdf(");
	  c->Print(pdfOpen.c_str());
	  first = 0;
	}
	else if (t==lastT && p==lastP) {
	  string pdfClose = string(outFName)+string(".pdf)");
	  c->Print(pdfClose.c_str());
	}
	else {
	  string pdfName =  string(outFName)+string(".pdf");
	  c->Print(pdfName.c_str());
	}
      }
    }
  }
}
/**********************************************************************/
void getSB(RooAddPdf &model_all, RooRealVar &x, const char *bgn, int i, int p, int t)
{
  // Integrals
  // - of S and B from "sig" and <bgn> PDFs,
  // - over "SBRange" range assigned to arg. <x> variable.
  // See: https://root-forum.cern.ch/t/integrating-a-rooabspdf-compoment-from-rooaddpdf/5361
  // S and B stored in array "Ns" and "Bs" (of [channel][p][t]).
  const RooArgList &list_all = model_all.pdfList();
  RooAbsPdf *sigAll = (RooAbsPdf*)list_all.find("sig");
  RooAbsPdf *bgnAll = (RooAbsPdf*)list_all.find(bgn);
  if (sigAll) {
    RooAbsReal *SAll = sigAll->createIntegral(x,x,"SBRange");
    Ns[i][p][t] = SAll->getVal();
  }
  else
    Ns[i][p][t] = 0;
  if (bgnAll) {
    RooAbsReal *BAll = bgnAll->createIntegral(x,x,"SBRange");
    Bs[i][p][t] = BAll->getVal();
  }
  else
    Bs[i][p][t] = 0;
}
/**********************************************************************/
void fit_table_K0(int cc){
  if( cc != 0 && cc != 1) return;
  cout << endl;
  if(cc == 0) cout << "Fits of K0 sample for pi- efficiency:" << endl;
  if(cc == 1) cout << "Fits of K0 sample for pi+ efficiency:" << endl;

  // Output ROOT and PDF files
  // (In two steps: i) to ROOT file, ii) from ROOT file to PDF, via "root2PDF").
  TFile *fOut = 0; TDirectory *rootOutDir = 0;
  char outFName[] = "test_K0_0"; sprintf(outFName,"test_K0_%d",cc);
  TDirectory *rootApp = gDirectory;

  for(int t = 0; t<Nt; t++){
    for(int p = 0; p<Np; p++){      // ********** LOOP ON [p][t] BINS
      if(t==3 && p>6) continue;
      Int_t ent = h[0+cc][0][p][t]->GetEntries(); //[p][t]
      cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << setw(7) << "Ent:" << ent << endl;
      if (ent<25)                   // ***** SKIP IF #ENTRIES TOO LOW
	// (Not much thinking into it: numerical value retained merely skips
	// most obvioulsy problematic cases...)
	continue;
      // ***** MODEL FOR K0
      // Signal
      RooRealVar x("x","M",0.44,0.56,"GeV");
      x.setBins(120);

      RooRealVar mean("mean","mean",0.5,0.47,0.52,"GeV") ;
      RooRealVar sigma1("sigma1","sigma1",0.01,0.,0.1,"GeV");
      RooRealVar sigma2("sigma2","sigma2",0.005,0.,0.1,"GeV");

      RooGaussian gauss1("gauss1","gauss1",x,mean,sigma1) ;
      RooGaussian gauss2("gauss2","gauss2",x,mean,sigma2) ;

      // Background
      RooRealVar k0a("k0a","k0a",	-0.2,	-1.,	1.) ;
      RooRealVar k1a("k1a","k1a",	-0.3,	-1.,	1.) ;
      RooRealVar k2a("k2a","k2a",	0.,		-1.,	1.) ;
      RooChebychev bgna("bgna","bgna",x,RooArgSet(k0a,k1a,k2a));//,k3a,k4a)) ;

      RooRealVar k0p("k0p","k0p",	-0.6,	-1.,	1.) ;
      RooRealVar k1p("k1p","k1p",	-0.4,	-1.,	1.) ;
      RooRealVar k2p("k2p","k2p",	0.08,	-1.,	1.) ;
      RooChebychev bgnp("bgnp","bgnp",x,RooArgSet(k0p,k1p,k2p));//,k3p,k4p)) ;

      //Construct composite pdf
      //RooRealVar N_a_b("N_a_b","N_a_b",0.05*ent,0.,1.05*ent) ;

      ent = h[0+cc][1][p][t]->GetEntries(); //[p][t]
      RooRealVar N_pi_s("N_pi_s","N_pi_s",0.9*ent,0.,1.05*ent) ;
      RooRealVar N_pi_b("N_pi_b","N_pi_b",0.1*ent,0.,1.05*ent) ;

      ent = h[0+cc][2][p][t]->GetEntries(); //[p][t]
      RooRealVar N_k_s("N_k_s","N_k_s",0.3*ent,0.,1.05*ent) ;
      RooRealVar N_k_b("N_k_b","N_k_b",0.7*ent,0.,1.05*ent) ;

      ent = h[0+cc][3][p][t]->GetEntries(); //[p][t]
      RooRealVar N_p_s("N_p_s","N_p_s",0.2*ent,0.,1.05*ent) ;
      RooRealVar N_p_b("N_p_b","N_p_b",0.8*ent,0.,1.05*ent) ;

      ent = h[0+cc][4][p][t]->GetEntries(); //[p][t]
      RooRealVar N_u_s("N_u_s","N_u_s",0.8*ent,0.,1.05*ent) ;
      RooRealVar N_u_b("N_u_b","N_u_b",0.05*ent,0.,1.05*ent) ;

      /*
	// Don't remember what this is?...
	if(cc == 1 && t==2 && p >Np-4){
	cout << cc << " " << t << " " << p << endl;
	N_p_s.setVal(0.);
	}
      */

      RooRealVar frac_s("frac_s","frac_s",0.65,0.5,0.7) ; //
      RooAddPdf sig("sig","sig",RooArgList(gauss1,gauss2),frac_s) ;

      RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));
      RooFormulaVar N_a_b("N_a_b","N_pi_b + N_k_b + N_p_b + N_u_b",RooArgSet(N_pi_b,N_k_b,N_p_b,N_u_b));

      RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgna),RooArgList(N_a_s,N_a_b)) ;
      RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgna),RooArgList(N_pi_s,N_pi_b)) ; // checked with bgnpi
      RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgna),RooArgList(N_k_s,N_k_b)) ; // checked with bgnk
      RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgnp),RooArgList(N_p_s,N_p_b)) ;
      RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgna),RooArgList(N_u_s,N_u_b)) ;

      RooDataHist data_all("data_all","data_all",x,Import(*h[0+cc][0][p][t])); //[p][t]
      RooDataHist data_pi("data_pi","data_pi",x,Import(*h[0+cc][1][p][t])); //[p][t]
      RooDataHist data_k("data_k","data_k",x,Import(*h[0+cc][2][p][t])); //[p][t]
      RooDataHist data_p("data_p","data_p",x,Import(*h[0+cc][3][p][t])); //[p][t]
      RooDataHist data_unk("data_unk","data_unk",x,Import(*h[0+cc][4][p][t])); //[p][t]

      RooCategory sample("sample","sample") ;
      sample.defineType("all") ;
      sample.defineType("pi") ;
      sample.defineType("k") ;
      sample.defineType("p") ;
      sample.defineType("unk") ;

      RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));

      RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
      simPdf.addPdf(model_all,"all") ;
      simPdf.addPdf(model_pi,"pi") ;
      simPdf.addPdf(model_k,"k") ;
      simPdf.addPdf(model_p,"p") ;
      simPdf.addPdf(model_unk,"unk") ;

      if(!use_sidebins){
	//RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
	RooAbsReal* nll = simPdf.createNLL(combData,Extended(true)) ;
	//RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
	//RooMinuit minu(nll_r) ;
	RooMinuit minu(*nll) ;
	minu.setPrintLevel(-1);
	minu.setNoWarn();
	int iter = 0, status, improved; double nLLMn = 0; do {
	  if(iter>0){
	    N_pi_s.setVal(1.01*N_pi_s.getVal());
	    N_k_s.setVal( 0.98*N_k_s.getVal());
	    N_p_s.setVal( 0.98*N_p_s.getVal());
	    N_u_s.setVal( 0.98*N_u_s.getVal());
	  }

	  minu.hesse();
	  minu.simplex();
	  int ret = minu.migrad();
	  if (ret) {
	    ret = minu.migrad();
	    if (!ret) printf("=== Sucessful réessai\n");
	  }
	  if (!ret) {
	    minu.improve(); improved = 1;
	  }
	  else              improved = 0;
	  minu.hesse();

	  r[0+cc][p][t] = minu.save();
	  double minNll = r[0+cc][p][t]->minNll();
	  if (!nLLMn || minNll<nLLMn) nLLMn = minNll;
	  status = 1-improved;
	  if (status==0) {
	    double effs[4], dEffs[4]; status = getEffDEff(r[0+cc][p][t],0,effs,dEffs);
	    if (status) {
	      double val = effs[0], err = dEffs[0], bDEff = sqrt(val*(1-val)/N_a_s.getVal());
	      printf("Eff+/-dEff = %.2f+/-%.2f >> %.2f\n",val,err,bDEff);
	    }
	  }
	  iter++;
	} while ((r[0+cc][p][t]->covQual()!=3 || status) && iter<=retry);
	printf("=> %2d iters (%.2f) -> %d",iter,nLLMn,r[0+cc][p][t]->covQual());
	double minNll = r[0+cc][p][t]->minNll();
	if (minNll>nLLMn*0.9) printf(" %.02f!\n",minNll);
	else                  printf("\n");
	// ***** GET S and B in K0 range
	// (Range is fixed and set to minmise fluctuation in peak position.)
	x.setRange("SBRange",M_K0-.03,M_K0+.03); getSB(model_all,x,"bgna",0+cc,p,t);
      }else {

	N_id[0+cc][0][p][t] = h[0+cc][0][p][t]->Integral(33,85) - h[0+cc][0][p][t]->Integral(3,29) - h[0+cc][0][p][t]->Integral(89,115);
	N_id[0+cc][1][p][t] = ( h[0+cc][1][p][t]->Integral(33,85) - h[0+cc][1][p][t]->Integral(3,29) - h[0+cc][1][p][t]->Integral(89,115) )/N_id[0+cc][0][p][t];
	N_id[0+cc][2][p][t] = ( h[0+cc][2][p][t]->Integral(33,85) - h[0+cc][2][p][t]->Integral(3,29) - h[0+cc][2][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];
	N_id[0+cc][3][p][t] = ( h[0+cc][3][p][t]->Integral(33,85) - h[0+cc][3][p][t]->Integral(3,29) - h[0+cc][3][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];
	N_id[0+cc][4][p][t] = ( h[0+cc][4][p][t]->Integral(33,85) - h[0+cc][4][p][t]->Integral(3,29) - h[0+cc][4][p][t]->Integral(89,115) ) /N_id[0+cc][0][p][t];

      }

      // ***** SAVE SIGNALS AND OVERFALL BACKGROUND
      N_id[0+cc][0][p][t] = N_a_s.getVal();
      N_id[0+cc][1][p][t] = N_pi_s.getVal();
      N_id[0+cc][2][p][t] = N_k_s.getVal() ;
      N_id[0+cc][3][p][t] = N_p_s.getVal() ;
      N_id[0+cc][4][p][t] = N_u_s.getVal() ;
      N_id[0+cc][5][p][t] = N_a_b.getVal() ;

      TGaxis::SetMaxDigits(3);
      stringstream nn;
      nn.str("");
      nn << "all" << (cc?'+':'-') << ": " << p_bins[p] << " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
      RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;

      RooPlot* frame2 = x.frame(Title("pi")) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;

      RooPlot* frame3 = x.frame(Title("k")) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;

      RooPlot* frame4 = x.frame(Title("p")) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgnp"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;

      RooPlot* frame5 = x.frame(Title("noID")) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;

      TCanvas* c = new TCanvas("c","",1280,384);
      RooPlot *frames[] = {frame1,frame2,frame3,frame4,frame5};
      for (int iF = 0; iF<5; iF++) {
	RooPlot *&frame = frames[iF];
	frame->GetXaxis()->SetTitleOffset(1.01);
	frame->GetXaxis()->SetTitleSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetYaxis()->SetTitleSize(0.045);
	frame->GetYaxis()->SetLabelSize(0.045);
	frame->GetYaxis()->SetNdivisions(505);

      }
      TPad* pa[5];
      pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
      pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
      pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
      pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
      pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
      for(int dr = 0; dr<5;dr++){
	pa[dr]->SetLeftMargin(0.12);
	pa[dr]->SetRightMargin(0.08);
	pa[dr]->Draw();
      }
      pa[0]->cd() ; frame1->Draw() ;
      pa[1]->cd() ; frame2->Draw() ;
      pa[2]->cd() ; frame3->Draw() ;
      pa[3]->cd() ; frame4->Draw() ;
      pa[4]->cd() ; frame5->Draw() ;

      char cN[] = "cpP99T99", mp[] = "mp"; sprintf(cN,"c%cP%dT%d",mp[cc],p,t);
      char cT[] = "cpP99T99_K0"; sprintf(cT,"c%cP%dT%d_K0",mp[cc],p,t);
      c->SetName(cN); c->SetTitle();

      // Output to ROOT and PDF files
      if (!fOut) {
	string rootName = string(outFName)+string(".root");
	fOut = TFile::Open(rootName.c_str(),"RECREATE");
	c->Write(); rootOutDir = gDirectory;
      }
      else {
	rootOutDir->cd(); c->Write();
      }
      delete c; rootApp->cd();
    }
  }
  if (fOut) { // Write to PDF and close ROOT file
    fOut->cd(); root2PDF(outFName,cc); fOut->Close(); rootApp->cd();
  }
}

/**********************************************************************/
void fit_table_phi(int cc){
  if( cc != 0 && cc != 1) return;
  cout << endl;
  if(cc == 0) cout << "Fits of PHI sample for K- efficiency:" << endl;
  if(cc == 1) cout << "Fits of PHI sample for K+ efficiency:" << endl;

  // Output ROOT and PDF files
  // (In two steps: i) to ROOT file, ii) from ROOT file to PDF, via "root2PDF").
  TFile *fOut = 0; TDirectory *rootOutDir = 0;
  char outFName[] = "test_iphi_0"; sprintf(outFName,"test_iphi_%d",cc);
  TDirectory *rootApp = gDirectory;

  for(int t = 0; t<Nt; t++){
    for(int p = 0; p<Np; p++){      // ********** LOOP ON [p][t] BINS
      if (p<2) continue;            // ***** Skip if p < 7 GeV, see "p_bins" in header file
      Int_t ent = h[2+cc][0][p][t]->GetEntries(); //[p][t]
      cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << setw(7) << "Ent:" << ent << endl;
      if (ent<25)                   // ***** SKIP IF #ENTRIES TOO LOW
	// (Not much thinking into it: numerical value retained merely skips
	// most obvioulsy problematic cases...)
	continue;
      //*****  MODEL FOR phi
      // Signal
      RooRealVar x("x","M",0.995,1.042,"GeV");

      RooRealVar mass("mass","mass #phi",1.01945,"GeV") ;
      RooRealVar width("width","width #phi",0.00426,"GeV");
      RooRealVar mean("mean","mean",0,-.04,.04,"GeV") ;
      RooRealVar sigma("sigma","sigma",0.0031,0.0025,0.005,"GeV");
      RooRealVar d1("d1","d1",1.);
      RooRealVar d2("d2","d2",0.4937);

      RooRelBreitWigner sb("sb","relBW",x,mass,width,d1,d2,d2,d1);
      RooGaussian sg("sg","gauss",x,mean,sigma);


      RooFFTConvPdf sig("sig","sig",x,sb,sg);
      //sig.setShift(0,0);

      // Background
      RooRealVar thr("thr","threshold",0.987354,0.985,0.99,"GeV");
      RooRealVar thr2("thr2","threshold",0.987354,0.985,0.99,"GeV");
      RooRealVar thr3("thr3","threshold",0.987354,0.985,0.99,"GeV");
      RooRealVar thr4("thr4","threshold",0.987354,0.985,0.99,"GeV");
      RooRealVar thr5("thr5","threshold",0.987354,0.985,0.99,"GeV");

      RooRealVar n("n","n",  0.4,  0.,   1.5) ;
      RooRealVar a("a","a",  4.,   0.,    100.) ;
      RooGenericPdf bgn1("bgn1","x<thr ? 0 :TMath::Power(x - thr,n)*TMath::Exp(-a*(x - thr))",RooArgList(x,thr,n,a));

      RooRealVar n2("n2","n2",  0.4,  0.,   1.5) ;
      RooRealVar a2("a2","a2",  4.,   0.,    100.) ;
      RooGenericPdf bgn2("bgn2","x<thr2 ? 0 :TMath::Power(x - thr2,n2)*TMath::Exp(-a2*(x - thr2))",RooArgList(x,thr2,n2,a2));

      RooRealVar n3("n3","n3",  0.4,  0.,   1.5) ;
      RooRealVar a3("a3","a3",  4.,   0.,    100.) ;
      RooGenericPdf bgn3("bgn3","x<thr3 ? 0 :TMath::Power(x - thr3,n3)*TMath::Exp(-a3*(x - thr3))",RooArgList(x,thr3,n3,a3));
      /*
      RooRealVar n4("n4","n4",  0.4,  0.,   1.5) ;
      RooRealVar a4("a4","a4",  4.,   0.,    100.) ;
      RooGenericPdf bgn4("bgn4","x<thr4 ? 0 :TMath::Power(x - thr4,n4)*TMath::Exp(-a4*(x - thr4))",RooArgList(x,thr4,n4,a4));

      RooRealVar n5("n5","n5",  0.4,  0.,   1.5) ;
      RooRealVar a5("a5","a5",  4.,   0.,    100.) ;
      RooGenericPdf bgn5("bgn5","x<thr5 ? 0 :TMath::Power(x - thr5,n5)*TMath::Exp(-a5*(x - thr5))",RooArgList(x,thr5,n5,a5));
      */

      //Construct composite pdf
      //RooRealVar N_a_b("N_a_b","N_a_b",0.2*ent,0.,1.05*ent) ;

      ent = h[2+cc][1][p][t]->GetEntries(); //[p][t]
      RooRealVar N_pi_s("N_pi_s","N_pi_s",0.12*ent,0.,1.05*ent) ;
      RooRealVar N_pi_b("N_pi_b","N_pi_b",0.5*ent,0.,1.05*ent) ;
      RooRealVar N_pi_b2("N_pi_b2","N_pi_b2",0.38*ent,0.,1.05*ent) ;

      ent = h[2+cc][2][p][t]->GetEntries(); //[p][t]
      RooRealVar N_k_s("N_k_s","N_k_s",0.87*ent,0.,1.05*ent) ;
      RooRealVar N_k_b("N_k_b","N_k_b",0.13*ent,0.,1.05*ent) ;

      ent = h[2+cc][3][p][t]->GetEntries(); //[p][t]
      RooRealVar N_p_s("N_p_s","N_p_s",0.87*ent,0.,1.05*ent) ;
      RooRealVar N_p_b("N_p_b","N_p_b",0.13*ent,0.,1.05*ent) ;

      ent = h[2+cc][4][p][t]->GetEntries(); //[p][t]
      RooRealVar N_u_s("N_u_s","N_u_s",0.77*ent,0.,1.05*ent) ;
      RooRealVar N_u_b("N_u_b","N_u_b",0.23*ent,0.,1.05*ent) ;
      RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));
      RooFormulaVar N_a_b("N_a_b","N_pi_b + N_k_b + N_p_b + N_u_b",RooArgSet(N_pi_b,N_k_b,N_p_b,N_u_b));


      //bedingung nur für t = 2

      if(p>0){
	double offset1 = 1;
	double offset2 = 1;
	thr.setVal(last_para[0]);
	thr2.setVal(last_para[1]);
	thr3.setVal(last_para[2]);
	thr4.setVal(last_para[3]);
	thr5.setVal(last_para[4]);
	n.setVal(last_para[5]);
	n2.setVal(last_para[6]);
	n3.setVal(last_para[7]);
	/*
	n4.setVal(last_para[8]);
	n5.setVal(last_para[9]);
	*/
	a.setVal(last_para[10]);
	a2.setVal(last_para[11]);
	a3.setVal(last_para[12]);
	/*
	a4.setVal(last_para[13]);
	a5.setVal(last_para[14]);
	*/
	N_pi_s.setVal(offset2*last_para[15]);
	N_pi_b.setVal(offset2*last_para[16]);
	N_k_s.setVal(offset1*last_para[17]);
	N_k_b.setVal(offset1*last_para[18]);
	N_p_s.setVal(offset2*last_para[19]);
	N_p_b.setVal(offset2*last_para[20]);
	N_u_s.setVal(offset2*last_para[21]);
	N_u_b.setVal(offset2*last_para[22]);
	//N_a_b.setVal(last_para[23]);
	mean.setVal(last_para[24]);
	sigma.setVal(last_para[25]);
	mean.setVal(last_para[26]);
      }


      RooArgList* lst_k_pdf = new RooArgList;
      lst_k_pdf->add(sig);
      lst_k_pdf->add(bgn3);

      RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgn1),RooArgList(N_a_s,N_a_b)) ;
      RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgn2),RooArgList(N_pi_s,N_pi_b)) ;
      RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgn3),RooArgList(N_k_s,N_k_b)) ;
      RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgn1),RooArgList(N_p_s,N_p_b)) ;
      RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgn1),RooArgList(N_u_s,N_u_b)) ;

      x.setBins(30);

      RooDataHist data_all("data_all","data_all",x,Import(*h[2+cc][0][p][t])); //[p][t]
      RooDataHist data_pi("data_pi","data_pi",x,Import(*h[2+cc][1][p][t])); //[p][t]
      RooDataHist data_k("data_k","data_k",x,Import(*h[2+cc][2][p][t])); //[p][t]
      RooDataHist data_p("data_p","data_p",x,Import(*h[2+cc][3][p][t])); //[p][t]
      RooDataHist data_unk("data_unk","data_unk",x,Import(*h[2+cc][4][p][t])); //[p][t]

      RooCategory sample("sample","sample") ;
      sample.defineType("all") ;
      sample.defineType("pi") ;
      sample.defineType("k") ;
      sample.defineType("p") ;

      RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));
      RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
      simPdf.addPdf(model_all,"all") ;
      simPdf.addPdf(model_pi,"pi") ;
      simPdf.addPdf(model_k,"k") ;
      simPdf.addPdf(model_p,"p") ;
      simPdf.addPdf(model_unk,"unk") ;

      if(!use_sidebins){

	//RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
	RooAbsReal* nll = simPdf.createNLL(combData,Extended(true));
	//RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
	//RooMinuit minu(nll_r) ;
	RooMinuit minu(*nll) ;
	minu.setPrintLevel(-1);
	minu.setNoWarn();
	int iter = 0, status, improved; double nLLMn = 0; do {
	  if(iter>0){
	    if(t == 0){
	      N_pi_s.setVal(0.95*N_pi_s.getVal());
	      N_k_s.setVal( 1.05*N_k_s.getVal() );
	      N_p_s.setVal( 0.90*N_p_s.getVal() );
	      N_u_s.setVal( 0.95*N_u_s.getVal() );
	    }else if(t == 1){
	      N_pi_s.setVal(0.94*N_pi_s.getVal());
	      N_k_s.setVal( 1.06*N_k_s.getVal() );
	      N_p_s.setVal( 0.91*N_p_s.getVal() );
	      N_u_s.setVal( 0.94*N_u_s.getVal() );
	    }else if(t == 2){
	      N_pi_s.setVal(0.91*N_pi_s.getVal());
	      N_k_s.setVal( 1.09*N_k_s.getVal() );
	      N_p_s.setVal( 0.91*N_p_s.getVal() );
	      N_u_s.setVal( 0.92*N_u_s.getVal() );
	    }else if(t == 3){
	      N_pi_s.setVal(0.90*N_pi_s.getVal());
	      N_k_s.setVal( 1.08*N_k_s.getVal() );
	      N_p_s.setVal( 0.90*N_p_s.getVal() );
	      N_u_s.setVal( 0.90*N_u_s.getVal() );
	    }
	  }

	  minu.hesse();
	  minu.simplex();
	  int ret = minu.migrad();
	  if (ret) {
	    ret = minu.migrad();
	    if (!ret) printf("=== Sucessful réessai\n");
	  }
	  if (!ret) {
	    minu.improve(); improved = 1;
	  }
	  else              improved = 0;
	  status = minu.hesse();

	  r[2+cc][p][t] = minu.save();
	  double minNll = r[2+cc][p][t]->minNll();
	  if (!nLLMn || minNll<nLLMn) nLLMn = minNll;
	  status = 1-improved;
	  if (status==0) {
	    double effs[4], dEffs[4]; status = getEffDEff(r[2+cc][p][t],1,effs,dEffs);
	    if (status) {
	      double val = effs[1], err = dEffs[1], bDEff = sqrt(val*(1-val)/N_a_s.getVal());
	      printf("Eff+/-dEff = %.2f+/-%.2f >> %.2f\n",val,err,bDEff);
	    }
	  }
	  iter++;
	} while ((r[2+cc][p][t]->covQual()!=3 || status) && iter<=retry);
	// 					-1 "Unknown, matrix was externally provided"
	// 					 0 "Not calculated at all"
	// 					 1 "Approximation only, not accurate"
	// 					 2 "Full matrix, but forced positive-definite"
	// 					 3 "Full, accurate covariance matrix"
	printf("=> %2d iters (%.2f) -> %d",iter,nLLMn,r[2+cc][p][t]->covQual());
	double minNll = r[2+cc][p][t]->minNll();
	if (minNll>nLLMn*0.9) printf(" %.02f!\n",minNll);
	else                  printf("\n");
	// ***** GET S and B in K0 range
	// (Range is fixed and set to minmise fluctuation in peak position.)
	x.setRange("SBRange",M_phi-.01,M_phi+.01); getSB(model_all,x,"bgn1",2+cc,p,t);
      }else {

	N_id[2+cc][0][p][t] = h[2+cc][0][p][t]->Integral(26,55) - h[2+cc][0][p][t]->Integral(10,24) - h[2+cc][0][p][t]->Integral(57,72);
	N_id[2+cc][1][p][t] = ( h[2+cc][1][p][t]->Integral(26,55) - h[2+cc][1][p][t]->Integral(10,24) - h[2+cc][1][p][t]->Integral(57,72) )/N_id[2+cc][0][p][t];
	N_id[2+cc][2][p][t] = ( h[2+cc][2][p][t]->Integral(26,55) - h[2+cc][2][p][t]->Integral(10,24) - h[2+cc][2][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];
	N_id[2+cc][3][p][t] = ( h[2+cc][3][p][t]->Integral(26,55) - h[2+cc][3][p][t]->Integral(10,24) - h[2+cc][3][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];
	N_id[2+cc][4][p][t] = ( h[2+cc][4][p][t]->Integral(26,55) - h[2+cc][4][p][t]->Integral(10,24) - h[2+cc][4][p][t]->Integral(57,72) ) /N_id[2+cc][0][p][t];

      }

      // ***** SAVE SIGNALS AND OVERFALL BACKGROUND
      N_id[2+cc][0][p][t] = N_a_s.getVal();
      N_id[2+cc][1][p][t] = N_pi_s.getVal();
      N_id[2+cc][2][p][t] = N_k_s.getVal();
      N_id[2+cc][3][p][t] = N_p_s.getVal();
      N_id[2+cc][4][p][t] = N_u_s.getVal();
      N_id[2+cc][5][p][t] = N_a_b.getVal() ;

      // ***** REMEMBER FIT PARAMETERS
      last_para[0] = thr.getVal();
      last_para[1] = thr2.getVal();
      last_para[2] = thr3.getVal();
      last_para[3] = thr4.getVal();
      last_para[4] = thr5.getVal();
      last_para[5] = n.getVal();
      last_para[6] = n2.getVal();
      last_para[7] = n3.getVal();
      /*
      last_para[8] = n4.getVal();
      last_para[9] = n5.getVal();
      */
      last_para[10] = a.getVal();
      last_para[11] = a2.getVal();
      last_para[12] = a3.getVal();
      /*
      last_para[13] = a4.getVal();
      last_para[14] = a5.getVal();
      */
      last_para[15] = N_pi_s.getVal();
      last_para[16] = N_pi_b.getVal();
      last_para[17] = N_k_s.getVal();
      last_para[18] = N_k_b.getVal();
      last_para[19] = N_p_s.getVal();
      last_para[20] = N_p_b.getVal();
      last_para[21] = N_u_s.getVal();
      last_para[22] = N_u_b.getVal();
      //last_para[23] = N_a_b.getVal();
      last_para[24] = mean.getVal();
      last_para[25] = sigma.getVal();
      last_para[26] = mean.getVal();
      // 			N_id[2+cc][0][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_a_s"))->getVal();
      // 			N_id[2+cc][1][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_pi"))->getVal();
      // 			N_id[2+cc][2][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_k"))->getVal();
      // 			N_id[2+cc][3][p][t] = ((RooRealVar*)r->floatParsFinal().find("frac_p"))->getVal();

      // 			N_id[2+cc][0][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_a_s"))->getVal();
      // 			N_id[2+cc][1][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_pi_s"))->getVal();
      // 			N_id[2+cc][2][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_k_s"))->getVal();
      // 			N_id[2+cc][3][p][t] = ((RooRealVar*)r->floatParsFinal().find("N_p_s"))->getVal();

      //			R_id[2+cc][0][p][t] = TMath::Sqrt(N_a_s1.getError()*N_a_s1.getError() + N_a_s2.getError()*N_a_s2.getError());
      //			R_id[2+cc][1][p][t] = TMath::Sqrt(N_pi_s1.getError()*N_pi_s1.getError() + N_pi_s2.getError()*N_pi_s2.getError());
      //			R_id[2+cc][2][p][t] = TMath::Sqrt(N_k_s1.getError()*N_k_s1.getError() + N_k_s2.getError()*N_k_s2.getError());
      //			R_id[2+cc][3][p][t] = TMath::Sqrt(N_p_s1.getError()*N_p_s1.getError() + N_p_s2.getError()*N_p_s2.getError());
      stringstream nn;
      TGaxis::SetMaxDigits(3);
      nn.str("");
      nn << "all" << (cc?'+':'-') << ": " << p_bins[p] << " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
      RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;

      RooPlot* frame2 = x.frame(Title("pi")) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgn2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;

      RooPlot* frame3 = x.frame(Title("k")) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgn3"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;

      RooPlot* frame4 = x.frame(Title("p")) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;

      RooPlot* frame5 = x.frame(Title("noID")) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgn1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;

      TCanvas* c = new TCanvas("c","",1280,384);
      RooPlot *frames[] = {frame1,frame2,frame3,frame4,frame5};
      for (int iF = 0; iF<5; iF++) {
	RooPlot *&frame = frames[iF];
	frame->GetXaxis()->SetTitleOffset(1.01);
	frame->GetXaxis()->SetTitleSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetYaxis()->SetTitleSize(0.045);
	frame->GetYaxis()->SetLabelSize(0.045);
	frame->GetYaxis()->SetNdivisions(505);
	frame->SetStats(0);
      }

      c->Clear();
      TPad* pa[5];
      pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
      pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
      pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
      pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
      pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
      for(int dr = 0; dr<5;dr++){
	pa[dr]->SetLeftMargin(0.12);
	pa[dr]->SetRightMargin(0.08);
	pa[dr]->Draw();
      }
      pa[0]->cd() ; frame1->Draw() ;
      // TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
      // bb->SetY1(0.3);
      pa[1]->cd() ; frame2->Draw() ;
      pa[2]->cd() ; frame3->Draw() ;
      pa[3]->cd() ; frame4->Draw() ;
      pa[4]->cd() ; frame5->Draw() ;

      char cN[] = "cpP99T99", mp[] = "mp"; sprintf(cN,"c%cP%dT%d",mp[cc],p,t);
      char cT[] = "cpP99T99_iphi"; sprintf(cT,"c%cP%dT%d_iphi",mp[cc],p,t);
      c->SetName(cN); c->SetTitle();

      // Output to ROOT and PDF files
      if (!fOut) {
	string rootName = string(outFName)+string(".root");
	fOut = TFile::Open(rootName.c_str(),"RECREATE");
	c->Write(); rootOutDir = gDirectory;
      }
      else {
	rootOutDir->cd(); c->Write();
      }
      delete c; rootApp->cd();
    }
  }
  if (fOut) { // Write to PDF and close ROOT file
    fOut->cd(); root2PDF(outFName,cc); fOut->Close(); rootApp->cd();
  }
}

/**********************************************************************/
void fit_table_Lambda(int cc){
  retry = 1;

  if( cc != 0 && cc != 1) return;
  cout << endl;
  if(cc == 0) cout << "Fits of LAMBDA sample for pbar efficiency:" << endl;
  if(cc == 1) cout << "Fits of LAMBDA sample for p efficiency:" << endl;

  // Output ROOT and PDF files
  // (In two steps: i) to ROOT file, ii) from ROOT file to PDF, via "root2PDF").
  TFile *fOut = 0; TDirectory *rootOutDir = 0;
  char outFName[] = "test_Lambda_0"; sprintf(outFName,"test_Lambda_%d",cc);
  TDirectory *rootApp = gDirectory;

  for(int t = 0; t<Nt;t++){
    for(int p = 0; p<Np;p++){
      for(int k = 0; k<30;k++){
	last_para2[p][t][k] = 0.;
      }
    }
  }

  for(int t = 0; t<Nt; t++){
    for(int p = 0; p<Np; p++){      // ********** LOOP ON [p][t] BINS
      if(t==3 && p>6) continue;
      Int_t ent = h[4+cc][0][p][t]->GetEntries(); //[p][t]
      cout << setw(7) << "theta:" << setw(3)<< t   << setw(7) << "mom:" << setw(3) << p << setw(7) << "Ent:" << ent << endl;
      if (ent<25)                   // ***** SKIP IF #ENTRIES TOO LOW
	// (Not much thinking into it: numerical value retained merely skips
	// most obvioulsy problematic cases...)
	continue;
      // ***** MODEL FOR Lambda
      // Signal
      RooRealVar x("x","M",1.1,1.13,"GeV");
      x.setBins(70);
      RooRealVar mean("mean","mean",1.115,1.11,1.12,"GeV") ;
      RooRealVar sigma1("sigma1","sigma1",0.004,0.0020,0.005,"GeV");
      RooRealVar frac_s("frac_s","frac_s",0.15,0.1,0.7) ;

      RooFormulaVar sigma2("sigma2","sigma1*frac_s",RooArgSet(sigma1,frac_s));

      RooGaussian gauss1("gauss1","gauss1",x,mean,sigma1) ;
      RooGaussian gauss2("gauss2","gauss2",x,mean,sigma2) ;


      RooAddPdf sig("sig","sig",RooArgList(gauss1,gauss2),frac_s) ;
      // 			RooVoigtian sig("sig","sig",x,mean,sigma2,sigma1);
      // Background
      RooRealVar n("n","n",  4.1,  0.,   7.5) ;
      RooRealVar a("a","a",  60.,   0.,    300.) ;
      // 			RooRealVar thr("thr","threshold",1.0778423,1.,1.09,"GeV");
      RooRealVar thr("thr","threshold",1.0778423);
      RooGenericPdf bgna("bgna","x<thr ? 0 :TMath::Power(x - thr,n)*TMath::Exp(-a*(x - thr))",RooArgList(x,thr,n,a));
      /*
      RooRealVar n2("n2","n2",  2.7,  0.,   7.5) ;
      RooRealVar a2("a2","a2",  19.,   0.,    300.) ;
      // 			RooRealVar thr2("thr2","threshold2",1.0778423,1.,1.09,"GeV");

      RooGenericPdf bgn2("bgn2","x<thr ? 0 :TMath::Power(x - thr,n2)*TMath::Exp(-a2*(x - thr))",RooArgList(x,thr,n2,a2));

      RooRealVar n3("n3","n3",  1.7,  0.,   7.5) ;
      RooRealVar a3("a3","a3",  31.,   0.,    300.) ;
      // 			RooRealVar thr3("thr3","threshold3",1.0778423,1.,1.09,"GeV");
      RooGenericPdf bgn3("bgn3","x<thr ? 0 :TMath::Power(x - thr,n3)*TMath::Exp(-a3*(x - thr))",RooArgList(x,thr,n3,a3));

      RooRealVar n4("n4","n4",  1.7,  0.,   7.5) ;
      RooRealVar a4("a4","a4",  31.,   0.,    300.) ;
      // 			RooRealVar thr4("thr4","threshold4",1.0778423,1.,1.09,"GeV");
      RooGenericPdf bgn4("bgn4","x<thr ? 0 :TMath::Power(x - thr,n4)*TMath::Exp(-a4*(x - thr))",RooArgList(x,thr,n4,a4));

      RooRealVar n5("n5","n5",  1.7,  0.,   7.5) ;
      RooRealVar a5("a5","a5",  31.,   0.,    300.) ;
      // 			RooRealVar thr5("thr5","threshold5",1.0778423,1.,1.09,"GeV");
      RooGenericPdf bgn5("bgn5","x<thr ? 0 :TMath::Power(x - thr,n5)*TMath::Exp(-a5*(x - thr))",RooArgList(x,thr,n5,a5));
      */
      //Construct composite pdf
      //RooRealVar N_a_b("N_a_b","N_a_b",0.45*ent,0.,1.05*ent) ;

      ent = h[4+cc][1][p][t]->GetEntries(); //[p][t]
      RooRealVar N_pi_s("N_pi_s","N_pi_s",0.13*ent,0.,1.05*ent) ;
      RooRealVar N_pi_b("N_pi_b","N_pi_b",0.82*ent,0.,1.05*ent) ;

      ent = h[4+cc][2][p][t]->GetEntries(); //[p][t]
      RooRealVar N_k_s("N_k_s","N_k_s",0.45*ent,0.,1.05*ent) ;
      RooRealVar N_k_b("N_k_b","N_k_b",0.5*ent,0.,1.05*ent) ;

      ent = h[4+cc][3][p][t]->GetEntries(); //[p][t]
      RooRealVar N_p_s("N_p_s","N_p_s",0.72*ent,0.,1.05*ent) ;
      RooRealVar N_p_b("N_p_b","N_p_b",0.11*ent,0.,1.05*ent) ;

      ent = h[4+cc][4][p][t]->GetEntries(); //[p][t]
      RooRealVar N_u_s("N_u_s","N_u_s",0.77*ent,0.,1.05*ent) ;
      RooRealVar N_u_b("N_u_b","N_u_b",0.11*ent,0.,1.05*ent) ;

      RooFormulaVar N_a_s("N_a_s","N_pi_s + N_k_s + N_p_s + N_u_s",RooArgSet(N_pi_s,N_k_s,N_p_s,N_u_s));
      RooFormulaVar N_a_b("N_a_b","N_pi_b + N_k_b + N_p_b + N_u_b",RooArgSet(N_pi_b,N_k_b,N_p_b,N_u_b));

      if(t == 1 && p == Np-2) {
	N_k_s.setVal(0);
	N_pi_s.setVal(0);
      }

      RooAddPdf model_all("model_all","model_all",RooArgList(sig,bgna),RooArgList(N_a_s,N_a_b)) ;
      RooAddPdf model_pi("model_pi","model_pi",RooArgList(sig,bgna),RooArgList(N_pi_s,N_pi_b)) ;
      RooAddPdf model_k("model_k","model_k",RooArgList(sig,bgna),RooArgList(N_k_s,N_k_b)) ;
      RooAddPdf model_p("model_p","model_p",RooArgList(sig,bgna),RooArgList(N_p_s,N_p_b)) ;
      RooAddPdf model_unk("model_unk","model_unk",RooArgList(sig,bgna),RooArgList(N_u_s,N_u_b)) ;


      RooDataHist data_all("data_all","data_all",x,Import(*h[4+cc][0][p][t])); //[p][t]
      RooDataHist data_pi("data_pi","data_pi",x,Import(*h[4+cc][1][p][t])); //[p][t]
      RooDataHist data_k("data_k","data_k",x,Import(*h[4+cc][2][p][t])); //[p][t]
      RooDataHist data_p("data_p","data_p",x,Import(*h[4+cc][3][p][t])); //[p][t]
      RooDataHist data_unk("data_unk","data_unk",x,Import(*h[4+cc][4][p][t])); //[p][t]



      RooCategory sample("sample","sample") ;
      sample.defineType("all") ;
      sample.defineType("pi") ;
      sample.defineType("k") ;
      sample.defineType("p") ;
      sample.defineType("unk") ;

      RooDataHist combData("combData","combined data",x,Index(sample),Import("all",data_all),Import("pi",data_pi),Import("k",data_k),Import("p",data_p),Import("unk",data_unk));
      RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
      simPdf.addPdf(model_all,"all") ;
      simPdf.addPdf(model_pi,"pi") ;
      simPdf.addPdf(model_k,"k") ;
      simPdf.addPdf(model_p,"p") ;
      simPdf.addPdf(model_unk,"unk") ;


      if(p>0 && t==0){
	thr.setVal(last_para2[t][p-1][0]);
	// 				thr2.setVal(last_para2[t][p-1][1]);
	// 				thr3.setVal(last_para2[t][p-1][2]);
	// 				thr4.setVal(last_para2[t][p-1][3]);
	// 				thr5.setVal(last_para2[t][p-1][4]);
	n.setVal(last_para2[t][p-1][5]);
	/*
	n2.setVal(last_para2[t][p-1][6]);
	n3.setVal(last_para2[t][p-1][7]);
	n4.setVal(last_para2[t][p-1][8]);
	n5.setVal(last_para2[t][p-1][9]);
	*/
	a.setVal(last_para2[t][p-1][10]);
	/*
	a2.setVal(last_para2[t][p-1][11]);
	a3.setVal(last_para2[t][p-1][12]);
	a4.setVal(last_para2[t][p-1][13]);
	a5.setVal(last_para2[t][p-1][14]);
	*/
	N_pi_s.setVal(last_para2[t][p-1][15]);
	N_pi_b.setVal(last_para2[t][p-1][16]);
	N_k_s.setVal(last_para2[t][p-1][17]);
	N_k_b.setVal(last_para2[t][p-1][18]);
	N_p_s.setVal(last_para2[t][p-1][19]);
	N_p_b.setVal(last_para2[t][p-1][20]);
	N_u_s.setVal(last_para2[t][p-1][21]);
	N_u_b.setVal(last_para2[t][p-1][22]);
	//N_a_b.setVal(last_para2[t][p-1][23]);
	mean.setVal(last_para2[t][p-1][24]);
	sigma1.setVal(last_para2[t][p-1][25]);
	frac_s.setVal(last_para2[t][p-1][26]);
      }else if(t>0){
	thr.setVal(last_para2[t-1][p][0]);
	// 				thr2.setVal(last_para2[t][p-1][1]);
	// 				thr3.setVal(last_para2[t][p-1][2]);
	// 				thr4.setVal(last_para2[t][p-1][3]);
	// 				thr5.setVal(last_para2[t][p-1][4]);
	n.setVal(last_para2[t-1][p][5]);
	/*
	n2.setVal(last_para2[t-1][p][6]);
	n3.setVal(last_para2[t-1][p][7]);
	n4.setVal(last_para2[t-1][p][8]);
	n5.setVal(last_para2[t-1][p][9]);
	*/
	a.setVal(last_para2[t-1][p][10]);
	/*
	a2.setVal(last_para2[t-1][p][11]);
	a3.setVal(last_para2[t-1][p][12]);
	a4.setVal(last_para2[t-1][p][13]);
	a5.setVal(last_para2[t-1][p][14]);
	*/
	N_pi_s.setVal(last_para2[t-1][p][15]);
	N_pi_b.setVal(last_para2[t-1][p][16]);
	N_k_s.setVal(last_para2[t-1][p][17]);
	N_k_b.setVal(last_para2[t-1][p][18]);
	N_p_s.setVal(last_para2[t-1][p][19]);
	N_p_b.setVal(last_para2[t-1][p][20]);
	N_u_s.setVal(last_para2[t-1][p][21]);
	N_u_b.setVal(last_para2[t-1][p][22]);
	//N_a_b.setVal(last_para2[t-1][p][23]);
	mean.setVal(last_para2[t-1][p][24]);
	sigma1.setVal(last_para2[t-1][p][25]);
	frac_s.setVal(last_para2[t-1][p][26]);
      }

      if(!use_sidebins){
	//RooFormulaVar restriction("restriction","0",RooArgSet());
	// 				RooFormulaVar restriction("restriction","100000*((TMath::Abs(1 - (N_pi_s + N_k_s + N_p_s + N_u_s)/N_a_s) > 1e-4) )",RooArgSet(N_a_s,N_pi_s,N_k_s,N_p_s,N_u_s));
	RooAbsReal* nll = simPdf.createNLL(combData,Extended(true)) ;
	//RooAddition nll_r("nll_r","nll_r",RooArgSet(*nll,restriction)) ;
	//RooMinuit minu(nll_r) ;
	RooMinuit minu(*nll) ;

	minu.setPrintLevel(-1);
	minu.setNoWarn();
	int iter = 0, status, improved; double nLLMn = 0; do {
	  mean.setVal(1.115) ;
	  sigma1.setVal(0.004);
	  frac_s.setVal(0.15) ;

	  if(iter>0){
	    if(cc==0){
	      N_pi_s.setVal(0.41*N_pi_s.getVal());
	      N_k_s.setVal( 0.57*N_k_s.getVal() );
	      N_p_s.setVal( 1.72*N_p_s.getVal() );
	      N_u_s.setVal( 0.34*N_u_s.getVal() );

	      N_pi_b.setVal( 1.02*N_pi_b.getVal() );
	      N_k_b.setVal(  1.02*N_k_b.getVal() );
	      N_p_b.setVal(  0.98*N_p_b.getVal() );
	      N_u_b.setVal(  1.02*N_u_b.getVal() );
	    }else{
	      N_pi_s.setVal(0.33*N_pi_s.getVal());
	      N_k_s.setVal( 0.42*N_k_s.getVal() );
	      N_p_s.setVal( 1.84*N_p_s.getVal() );
	      N_u_s.setVal( 0.27*N_u_s.getVal() );

	      N_pi_b.setVal( 1.01*N_pi_b.getVal() );
	      N_k_b.setVal(  1.01*N_k_b.getVal() );
	      N_p_b.setVal(  0.99*N_p_b.getVal() );
	      N_u_b.setVal(  1.01*N_u_b.getVal() );
	    }

	    // if(i%2 == 1){
	    // N_pi_s.setVal((1.-i/((double)retry*1.07))*N_pi_s.getVal());
	    // N_k_s.setVal( (1.-i/((double)retry*1.07))*N_k_s.getVal() );
	    // N_p_s.setVal( (1.+i/((double)retry*1.07))*N_p_s.getVal() );
	    // N_u_s.setVal( (1.-i/((double)retry*1.07))*N_u_s.getVal() );
	    // }else if(i%2 == 0){
	    // N_pi_s.setVal((1.+i/((double)retry*1.05))*N_pi_s.getVal());
	    // N_k_s.setVal( (1.+i/((double)retry*1.05))*N_k_s.getVal() );
	    // N_p_s.setVal( (1.-i/((double)retry*1.05))*N_p_s.getVal() );
	    // N_u_s.setVal( (1.+i/((double)retry*1.05))*N_u_s.getVal() );
	    // }


	    /*
	      if(t == 0){
	      if(i%2 == 0){
	      N_pi_s.setVal((1.-i/((double)retry*1.3))*N_pi_s.getVal());
	      N_k_s.setVal( (1.-i/((double)retry*1.3))*N_k_s.getVal() );
	      N_p_s.setVal( (1.+i/((double)retry*1.3))*N_p_s.getVal() );
	      N_u_s.setVal( (1.-i/((double)retry*1.3))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(0.90*N_pi_s.getVal());
	      // 								N_k_s.setVal( 0.90*N_k_s.getVal() );
	      // 								N_p_s.setVal( 1.10*N_p_s.getVal() );
	      // 								N_u_s.setVal( 0.90*N_u_s.getVal() );
	      }else if(i%2 == 1){
	      // 								N_pi_s.setVal(1.10*N_pi_s.getVal());
	      // 								N_k_s.setVal( 1.10*N_k_s.getVal() );
	      // 								N_p_s.setVal( 0.90*N_p_s.getVal() );
	      // 								N_u_s.setVal( 1.10*N_u_s.getVal() );

	      N_pi_s.setVal((1.+i/((double)retry*1.35))*N_pi_s.getVal());
	      N_k_s.setVal( (1.+i/((double)retry*1.35))*N_k_s.getVal() );
	      N_p_s.setVal( (1.-i/((double)retry*1.35))*N_p_s.getVal() );
	      N_u_s.setVal( (1.+i/((double)retry*1.35))*N_u_s.getVal() );
	      }
	      }else if(t == 1){
	      if(i%2 == 0){
	      N_pi_s.setVal((1.-i/((double)retry*1.3))*N_pi_s.getVal());
	      N_k_s.setVal( (1.-i/((double)retry*1.3))*N_k_s.getVal() );
	      N_p_s.setVal( (1.+i/((double)retry*1.3))*N_p_s.getVal() );
	      N_u_s.setVal( (1.-i/((double)retry*1.3))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(0.90*N_pi_s.getVal());
	      // 								N_k_s.setVal( 0.90*N_k_s.getVal() );
	      // 								N_p_s.setVal( 1.10*N_p_s.getVal() );
	      // 								N_u_s.setVal( 0.90*N_u_s.getVal() );
	      }else if(i%2 == 1){
	      // 								N_pi_s.setVal(1.10*N_pi_s.getVal());
	      // 								N_k_s.setVal( 1.10*N_k_s.getVal() );
	      // 								N_p_s.setVal( 0.90*N_p_s.getVal() );
	      // 								N_u_s.setVal( 1.10*N_u_s.getVal() );

	      N_pi_s.setVal((1.+i/((double)retry*1.35))*N_pi_s.getVal());
	      N_k_s.setVal( (1.+i/((double)retry*1.35))*N_k_s.getVal() );
	      N_p_s.setVal( (1.-i/((double)retry*1.35))*N_p_s.getVal() );
	      N_u_s.setVal( (1.+i/((double)retry*1.35))*N_u_s.getVal() );
	      }
	      }else if(t == 2){
	      if(i%2 == 0){
	      N_pi_s.setVal((1.-i/((double)retry*1.4))*N_pi_s.getVal());
	      N_k_s.setVal( (1.-i/((double)retry*1.4))*N_k_s.getVal() );
	      N_p_s.setVal( (1.+i/((double)retry*1.4))*N_p_s.getVal() );
	      N_u_s.setVal( (1.-i/((double)retry*1.4))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(0.50*N_pi_s.getVal());
	      // 								N_k_s.setVal( 0.50*N_k_s.getVal() );
	      // 								N_p_s.setVal( 1.50*N_p_s.getVal() );
	      // 								N_u_s.setVal( 0.50*N_u_s.getVal() );
	      }else if(i%2 == 1){
	      N_pi_s.setVal((1.+i/((double)retry*1.2))*N_pi_s.getVal());
	      N_k_s.setVal( (1.+i/((double)retry*1.2))*N_k_s.getVal() );
	      N_p_s.setVal( (1.-i/((double)retry*1.2))*N_p_s.getVal() );
	      N_u_s.setVal( (1.+i/((double)retry*1.2))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(1.50*N_pi_s.getVal());
	      // 								N_k_s.setVal( 1.50*N_k_s.getVal() );
	      // 								N_p_s.setVal( 0.50*N_p_s.getVal() );
	      // 								N_u_s.setVal( 1.50*N_u_s.getVal() );
	      }
	      }else if(t == 3){
	      if(i%2 == 0){
	      N_pi_s.setVal((1.-i/((double)retry*1.4))*N_pi_s.getVal());
	      N_k_s.setVal( (1.-i/((double)retry*1.4))*N_k_s.getVal() );
	      N_p_s.setVal( (1.+i/((double)retry*1.4))*N_p_s.getVal() );
	      N_u_s.setVal( (1.-i/((double)retry*1.4))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(0.50*N_pi_s.getVal());
	      // 								N_k_s.setVal( 0.50*N_k_s.getVal() );
	      // 								N_p_s.setVal( 1.50*N_p_s.getVal() );
	      // 								N_u_s.setVal( 0.50*N_u_s.getVal() );
	      }else if(i%2 == 1){
	      N_pi_s.setVal((1.+i/((double)retry*1.3))*N_pi_s.getVal());
	      N_k_s.setVal( (1.+i/((double)retry*1.3))*N_k_s.getVal() );
	      N_p_s.setVal( (1.-i/((double)retry*1.3))*N_p_s.getVal() );
	      N_u_s.setVal( (1.+i/((double)retry*1.3))*N_u_s.getVal() );

	      // 								N_pi_s.setVal(1.50*N_pi_s.getVal());
	      // 								N_k_s.setVal( 1.50*N_k_s.getVal() );
	      // 								N_p_s.setVal( 0.50*N_p_s.getVal() );
	      // 								N_u_s.setVal( 1.50*N_u_s.getVal() );
	      }
	      }*/
	  }
	  minu.hesse();
	  minu.simplex();
	  int ret = minu.migrad();
	  if (ret) {
	    ret = minu.migrad();
	    if (!ret) printf("=== Sucessful réessai\n");
	  }
	  if (!ret) {
	    minu.improve(); improved = 1;
	  }
	  else              improved = 0;
	  status = minu.hesse();
	  r[4+cc][p][t] = minu.save();
	  double minNll = r[4+cc][p][t]->minNll();
	  if (!nLLMn || minNll<nLLMn) nLLMn = minNll;
	  status = 1-improved;
	  if (status==0) {
	    double effs[4], dEffs[4]; status = getEffDEff(r[4+cc][p][t],2,effs,dEffs);
	    if (status) {
	      double val = effs[2], err = dEffs[2], bDEff = sqrt(val*(1-val)/N_a_s.getVal());
	      printf("Eff+/-dEff = %.2f+/-%.2f >> %.2f\n",val,err,bDEff);
	    }
	  }
	  iter++;
	} while ((r[4+cc][p][t]->covQual()!=3 || status) && iter<=retry);
	printf("=> %2d iters (%.2f) -> %d",iter,nLLMn,r[4+cc][p][t]->covQual());
	double minNll = r[4+cc][p][t]->minNll();
	if (minNll>nLLMn*0.9) printf(" %.02f!\n",minNll);
	else                  printf("\n");
	// ***** GET S and B in K0 range
	// (Range is fixed and set to minmise fluctuation in peak position.)
	x.setRange("SBRange",M_Lam-.006,M_Lam+.006); getSB(model_all,x,"bgna",4+cc,p,t);
      }else {

	N_id[4+cc][0][p][t] = h[4+cc][0][p][t]->Integral(18,37) - h[4+cc][0][p][t]->Integral(6,15) - h[4+cc][0][p][t]->Integral(40,50);
	N_id[4+cc][1][p][t] = ( h[4+cc][1][p][t]->Integral(18,37) - h[4+cc][1][p][t]->Integral(6,15) - h[4+cc][1][p][t]->Integral(40,50) )/N_id[4+cc][0][p][t];
	N_id[4+cc][2][p][t] = ( h[4+cc][2][p][t]->Integral(18,37) - h[4+cc][2][p][t]->Integral(6,15) - h[4+cc][2][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];
	N_id[4+cc][3][p][t] = ( h[4+cc][3][p][t]->Integral(18,37) - h[4+cc][3][p][t]->Integral(6,15) - h[4+cc][3][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];
	N_id[4+cc][4][p][t] = ( h[4+cc][4][p][t]->Integral(18,37) - h[4+cc][4][p][t]->Integral(6,15) - h[4+cc][4][p][t]->Integral(40,50) ) /N_id[4+cc][0][p][t];

      }

      // ***** SAVE SIGNALS AND OVERFALL BACKGROUND
      N_id[4+cc][0][p][t] = N_a_s.getVal();
      N_id[4+cc][1][p][t] = N_pi_s.getVal();
      N_id[4+cc][2][p][t] = N_k_s.getVal() ;
      N_id[4+cc][3][p][t] = N_p_s.getVal() ;
      N_id[4+cc][4][p][t] = N_u_s.getVal() ;
      N_id[4+cc][5][p][t] = N_a_b.getVal() ;


      if(p>-1){
	// ***** REMEMBER FIT PARAMETERS
	last_para2[p][t][0] = thr.getVal();
	// 				last_para2[p][t][1] = thr2.getVal();
	// 				last_para2[p][t][2] = thr3.getVal();
	// 				last_para2[p][t][3] = thr4.getVal();
	// 				last_para2[p][t][4] = thr5.getVal();
	last_para2[p][t][5] = n.getVal();
	/*
	last_para2[p][t][6] = n2.getVal();
	last_para2[p][t][7] = n3.getVal();
	last_para2[p][t][8] = n4.getVal();
	last_para2[p][t][9] = n5.getVal();
	*/
	last_para2[p][t][10] = a.getVal();
	/*
	last_para2[p][t][11] = a2.getVal();
	last_para2[p][t][12] = a3.getVal();
	last_para2[p][t][13] = a4.getVal();
	last_para2[p][t][14] = a5.getVal();
	*/
	last_para2[p][t][15] = N_pi_s.getVal();
	last_para2[p][t][16] = N_pi_b.getVal();
	last_para2[p][t][17] = N_k_s.getVal();
	last_para2[p][t][18] = N_k_b.getVal();
	last_para2[p][t][19] = N_p_s.getVal();
	last_para2[p][t][20] = N_p_b.getVal();
	last_para2[p][t][21] = N_u_s.getVal();
	last_para2[p][t][22] = N_u_b.getVal();
	last_para2[p][t][23] = N_a_b.getVal();
	last_para2[p][t][24] = mean.getVal();
	last_para2[p][t][25] = sigma1.getVal();
	last_para2[p][t][26] = frac_s.getVal();
      }


      stringstream nn;
      TGaxis::SetMaxDigits(3);
      nn.str("");
      nn << "all" << (cc?'+':'-') << ": " << p_bins[p] << " < p < " << p_bins[p+1] <<" , "<< t_bins[t] << " < #theta < " << t_bins[t+1];
      RooPlot* frame1 = x.frame(Title(nn.str().c_str())) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;
      // 			simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      // 			simPdf.plotOn(frame1,Slice(sample,"all"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame1,Slice(sample,"all"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame1,Cut("sample==sample::all")) ;
      // simPdf.paramOn(frame1,Layout(0.7,0.99,0.97));
      // frame1->getAttText()->SetTextSize(0.02);

      RooPlot* frame2 = x.frame(Title("pi")) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;
      // 			simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      // 			simPdf.plotOn(frame2,Slice(sample,"pi"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame2,Slice(sample,"pi"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame2,Cut("sample==sample::pi")) ;

      RooPlot* frame3 = x.frame(Title("k")) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;
      // 			simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      // 			simPdf.plotOn(frame3,Slice(sample,"k"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame3,Slice(sample,"k"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame3,Cut("sample==sample::k")) ;

      RooPlot* frame4 = x.frame(Title("p")) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;
      // 			simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      // 			simPdf.plotOn(frame4,Slice(sample,"p"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame4,Slice(sample,"p"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame4,Cut("sample==sample::p")) ;


      RooPlot* frame5 = x.frame(Title("noID")) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;
      // 			simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw),VisualizeError(*r,1),FillColor(kCyan)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),ProjWData(sample,combData),LineWidth(lw)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss1"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kRed)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("gauss2"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      // 			simPdf.plotOn(frame5,Slice(sample,"unk"),Components("sig"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kGreen)) ;
      simPdf.plotOn(frame5,Slice(sample,"unk"),Components("bgna"),ProjWData(sample,combData),LineStyle(kDashed),LineWidth(lw),LineColor(kOrange)) ;
      combData.plotOn(frame5,Cut("sample==sample::unk")) ;

      TCanvas* c = new TCanvas("c","",1280,384);
      RooPlot *frames[] = {frame1,frame2,frame3,frame4,frame5};
      for (int iF = 0; iF<5; iF++) {
	RooPlot *&frame = frames[iF];
	frame->GetXaxis()->SetTitleOffset(1.01);
	frame->GetXaxis()->SetTitleSize(0.045);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetTitleOffset(1.3);
	frame->GetYaxis()->SetTitleSize(0.045);
	frame->GetYaxis()->SetLabelSize(0.045);
	frame->GetYaxis()->SetNdivisions(505);
	frame->SetStats(0);
      }

      c->Clear();
      TPad* pa[5];
      pa[0] = new TPad("p1","",0.0,	0.,	0.2 ,1.);
      pa[1] = new TPad("p2","",0.2,	0.,	0.4 ,1.);
      pa[2] = new TPad("p3","",0.4,	0.,	0.6 ,1.);
      pa[3] = new TPad("p4","",0.6,	0.,	0.8 ,1.);
      pa[4] = new TPad("p5","",0.8,	0.,	1.0 ,1.);
      for(int dr = 0; dr<5;dr++){
	pa[dr]->SetLeftMargin(0.12);
	pa[dr]->SetRightMargin(0.08);
	pa[dr]->Draw();
      }
      pa[0]->cd() ; frame1->Draw() ;
      // TPaveText* bb = (TPaveText*)c->GetPad(1)->FindObject("simPdf_paramBox");
      // bb->SetY1(0.3);
      pa[1]->cd() ; frame2->Draw() ;
      pa[2]->cd() ; frame3->Draw() ;
      pa[3]->cd() ; frame4->Draw() ;
      pa[4]->cd() ; frame5->Draw() ;

      char cN[] = "cpP99T99", mp[] = "mp"; sprintf(cN,"c%cP%dT%d",mp[cc],p,t);
      char cT[] = "cpP99T99_Lambda"; sprintf(cT,"c%cP%dT%d_Lambda",mp[cc],p,t);
      c->SetName(cN);

      // Output to ROOT and PDF files
      if (!fOut) {
	string rootName = string(outFName)+string(".root");
	fOut = TFile::Open(rootName.c_str(),"RECREATE");
	c->Write(); rootOutDir = gDirectory;
      }
      else {
	rootOutDir->cd(); c->Write();
      }
      delete c; rootApp->cd();
    }
  }
  if (fOut) { // Write to PDF and close ROOT file
    fOut->cd(); root2PDF(outFName,cc); fOut->Close(); rootApp->cd();
  }
}

/**********************************************************************/
void print_table(){
  input_K0->Close();
  input_Lam->Close();
  input_iphi->Close();

  TGraphErrors *gr[6][4][Nt]; // Efficiency/purity; [1-3]=[piKpu]
  TGraph* grC[6][Nt];

  stringstream nn;

  TFile* ff= new TFile(out_file.c_str(),"RECREATE");
  TDirectory *dGraphs = gDirectory;
  TDirectory *dStats = gDirectory->mkdir("Stats");
  TDirectory *dFitResults = gDirectory->mkdir("FitResults");

  string p_name[6] = {"pi_m","pi_p","K_m","K_p","p_m","p_p"};
  string id[4] = {"pi","K","p","u"};
  Int_t color[4] = {kOrange+7,kAzure+4,kTeal+4,kYellow+2};

  ofstream ofs_matrix("rich_mat.txt", std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_err("rich_err.txt", std::ofstream::out | std::ofstream::trunc);

  int color2[7] = {kRed, kOrange+7, kYellow+2, kSpring-6, kCyan-6, kAzure+3, kViolet-1};
  // 	double shift[3] = {-0.005,0,0.005};
  // double shift[4] = {0.2,0.,0.1,0.3};
  double shift[4] = {0.,0.1,0.2,0.3};

  //int cov_elem[4] = {6,2,4,8}; // pi,k,p,u
  int cov_elem[4] = {5,1,3,7}; // pi,k,p,u

  for(int p = 0; p< Np; p++)
    {
      for(int t = 1; t< 3; t++)
	{
	  ofs_matrix << p << "\t" << t;
	  ofs_err << p << "\t" << t;

	  for(int i = start; i<stop; i++)
	    {

	      if (!r[i][p][t])
		printf("p,t,i: %d,%d,%d %s\n",p,t,i,"Empty");

	      double tmp_jak[4];
	      double ggg = N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t];
	      for (int j = 1; j<4; j++) {
		  double aaa = N_id[i][j][p][t];
		  double val = (ggg ? aaa/ggg : 0);
		  ofs_matrix << "\t" << val;

		  for(int ll=0; ll<4;ll++){
		    tmp_jak[ll] = -aaa/(ggg*ggg);
		    if(ll == j-1) tmp_jak[ll] += 1./ggg;
		  }

		}
	      if (!r[i][p][t]) {
		ofs_err << "\t" << 0 << "\t" << 0 << "\t" << 0 
			<< "\t" << 0 << "\t" << 0 << "\t" << 0;
	      }
	      else {
		ofs_err << "\t" << tmp_jak[0]*tmp_jak[0]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[0],cov_elem[0])
			<< "\t" << tmp_jak[1]*tmp_jak[1]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[1],cov_elem[1])
			<< "\t" << tmp_jak[2]*tmp_jak[2]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[2],cov_elem[2]);
		ofs_err << "\t" << tmp_jak[0]*tmp_jak[1]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[0],cov_elem[1])
			<< "\t" << tmp_jak[0]*tmp_jak[2]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[0],cov_elem[2])
			<< "\t" << tmp_jak[1]*tmp_jak[2]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[1],cov_elem[2]);
	      }
	    }
	  ofs_matrix << endl;
	  ofs_err << endl;
	}
    }

  ofs_matrix.close();
  ofs_err.close();

  dFitResults->cd();
  for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    for(int t = 0; t< Nt; t++){
      grC[i][t] = new TGraph();
      grC[i][t]->SetMarkerStyle(33);
      grC[i][t]->SetMarkerSize(2);
      grC[i][t]->SetMarkerColor(color2[t]);
      grC[i][t]->SetLineColor(color2[t]);
      int ipt, p; for(p=ipt = 0; p< Np; p++){
	if (!r[i][p][t]) continue; // E.g. if (t==3 && p>6)
	// 				diff[p] = (N_id[i][1][p][t] + N_id[i][2][p][t] + N_id[i][3][p][t] +N_id[i][4][p][t]) + shift[t];
	double P = p_bins[p]/2. + p_bins[p+1]/2.;
	double diff = (r[i][p][t]->covQual()/3.) + shift[t];
	grC[i][t]->SetPoint(ipt++,P,diff);
      }
      grC[i][t]->Write();
    }
  }
  dGraphs->cd();
  for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    for (int t = 0; t< Nt; t++) {
      for (int j = 1; j<5; j++) {  // identified as (pi,k,p)
	cout << i << " " << j << " " << t << " ";
	gr[i][j-1][t] = new TGraphErrors();
	gr[i][j-1][t]->SetMarkerStyle(33);
	gr[i][j-1][t]->SetMarkerColor(color[j-1]);
	gr[i][j-1][t]->SetLineColor(color[j-1]);
	nn.str("");
	nn.clear();
	nn << p_name[i] << "_" << id[j-1] << "_" << t;
	gr[i][j-1][t]->SetName(nn.str().c_str());
	int ipt, p; for(p=ipt = 0; p< Np; p++){
	  if (!r[i][p][t]) continue;
	  double aaa = N_id[i][j][p][t];
	  double ggg = N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t];
	  double val = (ggg ? aaa/ggg : 0);
	  cout << val << " ";
	  double P = p_bins[p]/2. + p_bins[p+1]/2.;
	  gr[i][j-1][t]->SetPoint(ipt,P,val);

	  double tmp = 0.;
	  double tmp_jak[4];

	  for(int ll=0; ll<4;ll++){
	    tmp_jak[ll] = -aaa/(ggg*ggg);
	    if(ll == j-1) tmp_jak[ll] += 1./ggg;
	  }
	  for(int ll = 0; ll<4;ll++) tmp += tmp_jak[ll]*tmp_jak[ll]*r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[ll]);
	  for(int ll = 0; ll<3; ll++){
	    for(int hh = ll+1; hh<4; hh++){
	      tmp += tmp_jak[ll]*tmp_jak[hh]*2.*r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[hh]);
	    }
	  }
	  double err = TMath::Sqrt(tmp);
	  // err[p] = TMath::Sqrt((aaa+1)*(ggg-aaa+1)/( (ggg+2)*(ggg+2)*(ggg+3) ));
	  gr[i][j-1][t]->SetPointError(ipt++,0,err);	  
	  //#define DEBUG_EFFDEFF
#ifdef DEBUG_EFFDEFF
	  static int idebug = 0;
	  if (idebug) {
	    int h2h = i/2, status;
	    double effs[4], dEffs[4]; status = getEffDEff(r[i][p][t],h2h,effs,dEffs);
	    printf("#%d/%d [%d][%d] => %.2f+/-%.2f  %.2f+/-%.2f  %.2f %d\n",i,j-1,p,t,
		   effs[j-1],dEffs[j-1],val,err,sqrt(val*(1-val)/ggg),status);
	  }
#endif
	}
	cout << endl;
	gr[i][j-1][t]->Write();
      }
    }
  }
  dStats->cd();
  for (int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    for (int t = 0; t< Nt; t++) {
      TGraphErrors stS;
      stS.SetMarkerStyle(20);
      nn.str("");
      nn.clear();
      nn << "stS_" << p_name[i] << "_" << t;
      stS.SetName(nn.str().c_str());
      int ipt, p; for (p=ipt = 0; p< Np; p++) {
	if (!r[i][p][t]) continue;
	double fS = Ns[i][p][t], fB = Bs[i][p][t];
	double B = N_id[i][5][p][t]*fB;
	double S = (N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t])*fS;
	double statS, dStatS; if (S+B>0) {
	  // dfdSj = (S/2+B)/(S+B)^3/2
	  // df/dB = -S/2/(S+B)^3/2
	  double dfdSj = (S/2+B)/sqrt(S+B)/(S+B);
	  double dfdB = -S/2/sqrt(S+B)/(S+B);
	  double dS2 = 0, dSB = 0; for (int ll = 0; ll<4; ll++) {
	    dS2 += r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[ll]);
	    dSB += r[i][p][t]->covarianceMatrix()(cov_elem[ll],0);
	  }
	  double dB2 = r[i][p][t]->covarianceMatrix()(0,0);
	  statS = S/sqrt(S+B); dStatS = dfdSj*dfdSj*dS2+dfdSj*dfdB*dSB+dfdB*dfdB*dB2;
	}
	else {
	  statS=dStatS = 0;
	}
	double P = p_bins[p]/2. + p_bins[p+1]/2.;    
	stS.SetPoint(ipt,P,statS); stS.SetPointError(ipt++,0,dStatS);
      }
      stS.Write();
      TGraphErrors sB;
      sB.SetMarkerStyle(20);
      nn.str(""); nn.clear(); nn << "sB_" << p_name[i] << "_" << t;
      sB.SetName(nn.str().c_str());
      for (p=ipt = 0; p< Np; p++) {
	if (!r[i][p][t]) continue;
	double  fS = Ns[i][p][t], fB = Bs[i][p][t], B = N_id[i][5][p][t]*fB;
	double sigB, dSigB; if (B) {
	  double S = (N_id[i][1][p][t]+N_id[i][2][p][t]+N_id[i][3][p][t]+N_id[i][4][p][t])*fS;
	  // dfdSj = 1/B
	  // df/dB = -S/B ^2
	  double dfdSj = 1/B;
	  double dfdB = -S/2/B/B;
	  double dS2 = 0, dSB = 0; for (int ll = 0; ll<4; ll++) {
	    dS2 += r[i][p][t]->covarianceMatrix()(cov_elem[ll],cov_elem[ll]);
	    dSB += r[i][p][t]->covarianceMatrix()(cov_elem[ll],0);
	  }
	  double dB2 = r[i][p][t]->covarianceMatrix()(0,0);
	  sigB = S/B; dSigB = dfdSj*dfdSj*dS2+dfdSj*dfdB*dSB+dfdB*dfdB*dB2;
	}
	else {
	  sigB=dSigB = 0;
	}
	double P = p_bins[p]/2. + p_bins[p+1]/2.;
	sB.SetPoint(ipt,P,sigB); sB.SetPointError(ipt++,0,dSigB);
      }
      sB.Write();
    }
  }
  dFitResults->cd();
  for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    for(int t = 0; t< Nt; t++){
      for(int p = 0; p< Np; p++){
	if (!r[i][p][t]) continue;
	nn.str("");
	nn.clear();
	nn <<"fit_" << p_name[i] << "_" << t << "_" << p;
	r[i][p][t]->Write(nn.str().c_str());
      }
    }
  }
  ff->Close();
  gStyle->SetOptStat(0);

  for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    for(int j = 1; j<5; j++) {  // identified as (pi,k,p)
      for(int t = 0; t< Nt; t++){
	gr[i][j-1][t]->SetMarkerColor(color2[t]);
	gr[i][j-1][t]->SetLineColor(color2[t]);
      }
    }
  }

  TCanvas* c = new TCanvas("c","",800,600);
  TLegend *leg = new TLegend(0.13,0.25,0.33,0.55);
  leg->SetLineColor(0);
  leg->SetNColumns(1);
  leg->SetFillStyle(0);


  for(int t = 0; t< Nt; t++){
    nn.str("");
    nn << t_bins[t] <<" #leq #theta < " << t_bins[t+1];
    leg->AddEntry(gr[start][0][t],nn.str().c_str(),"pl");
  }

  string label[8] = {"#pi^{-}","#pi^{+}","K^{-}","K^{+}","#bar{p}","p","unknown","unknown"};

  double min[3][4];
  double max[3][4];


  max[0][0]=1.0;
  max[0][1]=1.0;
  max[0][2]=1.0;
  max[0][3]=1.0;

  max[1][0]=1.0;
  max[1][1]=1.0;
  max[1][2]=1.0;
  max[1][3]=1.0;

  max[2][0]=1.0;
  max[2][1]=1.0;
  max[2][2]=1.0;
  max[2][3]=1.0;

  min[0][0]=0.;
  min[0][1]=0.;
  min[0][2]=0.;
  min[0][3]=0.;

  min[1][0]=0.;
  min[1][1]=0.;
  min[1][2]=0.;
  min[1][3]=0.;

  min[2][0]=0.;
  min[2][1]=0.;
  min[2][2]=0.;
  min[2][3]=0.;

  TLine* line = new TLine(0,1,50,1);
  line->SetLineStyle(2);
  line->SetLineColor(kGray+1);

  TLine* li[4];
  // 	li[0] = new TLine(0,0.995,50,0.995);
  // 	li[1] = new TLine(0,1.000,50,1.000);
  // 	li[2] = new TLine(0,1.005,50,1.005);
  li[0] = new TLine(0,1.0,50,1.0);
  li[1] = new TLine(0,1.1,50,1.1);
  li[2] = new TLine(0,1.2,50,1.2);
  li[3] = new TLine(0,1.3,50,1.3);
  for(int l = 0 ; l<4; l++){
    li[l]->SetLineStyle(2);
    li[l]->SetLineColor(kGray+1);
  }

  for(int i = start; i<stop; i++) {  // particle (pi,K,p) 0,6
    for(int j = 1; j<5; j++) {  // identified as (pi,K,p)
      c->Clear();
      TH1F* hr = c->DrawFrame(0.,min[i/2][j-1],50.,max[i/2][j-1]);
      hr->GetXaxis()->SetTitle("p (GeV/c)");
      nn.str("");
      nn << label[i] <<" #rightarrow " << label[2*(j-1)+i%2];
      hr->SetTitle(nn.str().c_str());
      leg->Draw();
      for(int t = 0 ; t< Nt; t++) gr[i][j-1][t]->Draw("pl");
      nn.str("");
      if(i%2 == 0) nn << "table/" << id[i/2] << "m_" << id[j-1] <<".pdf";
      else         nn << "table/" << id[i/2] << "p_" << id[j-1] <<".pdf";
      c->Print(nn.str().c_str());
    }
  }



  c->Print("check_sum.pdf(");
  for(int i = start; i<stop; i++) {  // particle (pi,k,p) 0,6
    c->Clear();
    // 		TH1F* hr = c->DrawFrame(0.,0.98,50.,1.02);
    TH1F* hr = c->DrawFrame(0.,0.,50.,1.37);
    hr->GetXaxis()->SetTitle("p (GeV/c)");
    hr->GetXaxis()->SetLabelSize(0.045);
    hr->GetXaxis()->SetTitleSize(0.045);
    hr->GetXaxis()->SetTitleOffset(0.9);
    hr->GetYaxis()->SetLabelSize(0.045);
    nn.str("");
    nn << "Check " << label[i];
    hr->SetTitle(nn.str().c_str());
    for(int t = 0 ; t< Nt; t++){
      li[t]->Draw("same");
      grC[i][t]->Draw("pl");
    }
    c->Print("check_sum.pdf");
  }
  c->Print("check_sum.pdf)");


  /*
    string charge[2] = {"n","p"};
    for(int c = 0; c<2;c++){
    for(int t = 0; t<Nt;t++){
    for(int p=0; p<Np;p++){
    TMatrixD eff3(3,3);
    nn.str("");
    nn << "eff3_"<<charge[c] <<"_t"<<t<<"_p"<<TMath::Nint(p_bins[p]);
    //0 pi
    //1 K
    //2 p

    eff3[0][2] = N_id[0+c][3][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
    eff3[0][1] = N_id[0+c][2][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
    eff3[0][0] = N_id[0+c][1][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);

    eff3[1][2] = N_id[2+c][3][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
    eff3[1][1] = N_id[2+c][2][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
    eff3[1][0] = N_id[2+c][1][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);

    eff3[2][2] = N_id[4+c][3][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff3[2][1] = N_id[4+c][2][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff3[2][0] = N_id[4+c][1][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);

    eff3.Write(nn.str().c_str());
    TMatrixD eff4(3,4);
    nn.str("");
    nn << "eff4_"<<charge[c] <<"_t"<<t<<"_p"<<TMath::Nint(p_bins[p]);
    //0 pi
    //1 K
    //2 p
    //3 u

    eff4[0][3] = N_id[0+c][4][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
    eff4[0][2] = N_id[0+c][3][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
    eff4[0][1] = N_id[0+c][2][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);
    eff4[0][0] = N_id[0+c][1][p][t]/(N_id[0+c][1][p][t]+N_id[0+c][2][p][t]+N_id[0+c][3][p][t]+N_id[0+c][4][p][t]);

    eff4[1][3] = N_id[2+c][4][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
    eff4[1][2] = N_id[2+c][3][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
    eff4[1][1] = N_id[2+c][2][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);
    eff4[1][0] = N_id[2+c][1][p][t]/(N_id[2+c][1][p][t]+N_id[2+c][2][p][t]+N_id[2+c][3][p][t]+N_id[2+c][4][p][t]);

    eff4[2][3] = N_id[4+c][4][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff4[2][2] = N_id[4+c][3][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff4[2][1] = N_id[4+c][2][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff4[2][0] = N_id[4+c][1][p][t]/(N_id[4+c][1][p][t]+N_id[4+c][2][p][t]+N_id[4+c][3][p][t]+N_id[4+c][4][p][t]);
    eff4.Write(nn.str().c_str());
    }
    }
    }
  */


  ff->Close();
}

/**********************************************************************/
void set_plot_style(){
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

/**********************************************************************/
void bookKineHistos()
{
  double ZMn = -4968, ZMx = 2322; int nZbins = (ZMx-ZMn)*4;
  ZMn /= 10; ZMx /=10;
  char tag[] =
    "D/#deltaD>4,c#theta>0.99997,pT>20MeV";
  //"pT>20MeV,EMiss<2.5GeV";
  size_t len = strlen(tag); string title;
  if      (analysis=="K0L") {
    am_all = new TH2D("am_all","All;#alpha;pT (GeV)",500,-1,1,500,0.,0.3);
    Z_all =  new TH1D("Z_all", "All;ZpV (cm)",nZbins,ZMn,ZMx);
    XY_all = new TH2D("XY_all","All;XpV (cm);YpV(cm)",100,-2.5,2.5,100,-2.5,2.5);
    am_K0 =  new TH2D("am_K0", "K0,#LambdaVeto;#alpha;pT (GeV)",
		      500,-1,1,500,0.,0.3);
    am_K0p = new TH2D("am_K0p","K0#rightarrow#pi+,#pi-ID;#alpha;pT (GeV)",
		      500,-1,1,500,0.,0.3);
    am_K0m = new TH2D("am_K0m","K0#rightarrow#pi-,#pi+ID;#alpha;pT (GeV)",
		      500,-1,1,500,0.,0.3);
    am_L =   new TH2D("am_L",  "#Lambda,K0Veto;#alpha;pT (GeV)",
		      500,-1,1,500,0.,0.3);
    Z_K0 =  new TH1D("Z_K0", "K0;ZpV (cm)",     nZbins,ZMn,ZMx);
    Z_L =   new TH1D("Z_L",  "#Lambda;ZpV (cm)",nZbins,ZMn,ZMx);
    XY_K0 = new TH2D("XY_K0","K0;XpV (cm);YpV (cm)",
		     100,-2.5,2.5,100,-2.5,2.5);
    XY_L =  new TH2D("XY_L", "#Lambda;XpV (cm);YpV (cm)",
		     100,-2.5,2.5,100,-2.5,2.5);
    snprintf(tag,len,"D/#deltaD>%.0f,c#theta>%.5f,pT>%.0fMeV",
	     DdD_cuts[0],cth_cuts[0],pT_cuts[0]*1000);
    title = string("K0 - ")+     string(tag)+string(";D/#deltaD");
    DdD_K0 = new TH1D("DdD_K0",title.c_str(),100,0,100);
    title = string("K0 - ")+     string(tag)+string(";cos#theta");
    cth_K0 = new TH1D("cth_K0",title.c_str(),100,.9998,1);
    cth_K0->GetXaxis()->SetNdivisions(505);
    title = string("K0 - ")+     string(tag)+string(";pT (GeV)");
    pT_K0 =  new TH1D("pT_K0", title.c_str(),100,0.,.25);
    title = string("K0 - ")+string(tag)+string(";tracks in pV");
    Tr_K0 =  new TH1D("Tr_K0", title.c_str(),32,-.5,31.5);
    title = string("K0 - ")+string(tag)+string(";pTracks in bottom RICH");
    Rb_K0 =  new TH1D("Rb_K0", title.c_str(),32,-.5,31.5);
    title = string("K0 - ")+string(tag)+string(";pTracks in top RICH");
    Rt_K0 =  new TH1D("Rt_K0", title.c_str(),32,-.5,31.5);
    snprintf(tag,len,"D/#deltaD>%.0f,c#theta>%.5f,pT>%.0fMeV",
	     DdD_cuts[1],cth_cuts[1],pT_cuts[1]*1000);
    title = string("#Lambda - ")+string(tag)+string(";D/#deltaD");
    DdD_L =  new TH1D("DdD_L", title.c_str(),100,0,100);
    title = string("#Lambda - ")+string(tag)+string(";cos#theta");
    cth_L =  new TH1D("cth_L", title.c_str(),100,.9998,1);
    cth_L->GetXaxis()->SetNdivisions(505);
    title = string("#Lambda - ")+string(tag)+string(";pT (GeV)");
    pT_L =   new TH1D("pT_L",  title.c_str(),100,0.,.25);
    title = string("#Lambda - ")+string(tag)+string(";tracks in pV");
    Tr_L =   new TH1D("Tr_L",  title.c_str(),32,-.5,31.5);
    title = string("#Lambda - ")+string(tag)+string(";pTracks in bottom RICH");
    Rb_L =   new TH1D("Rb_L",  title.c_str(),32,-.5,31.5);
    title = string("#Lambda - ")+string(tag)+string(";pTracks in top RICH");
    Rt_L =   new TH1D("Rt_L",  title.c_str(),32,-.5,31.5);
  }
  else {
    snprintf(tag,len,"pT>%.0fMeV",pT_cuts[2]*1000);
    title = string("Incl. - ")+string(tag)+string(";EMiss (GeV)");
    dE_Incl = new TH1D("dE_Incl",title.c_str(),100,-5,10);
    snprintf(tag,len,"pT>%.0fMeV",pT_cuts[3]*1000);
    title = string("Excl. - ")+string(tag)+string(";EMiss (GeV)");
    dE_Excl = new TH1D("dE_Excl",title.c_str(),100,-5,10);
    snprintf(tag,len,"EMiss>%.1fGeV",dE_cuts[0]);
    title = string("Incl. - ")+string(tag)+string(";#alpha;pT (GeV)");
    am_Incl = new TH2D("am_Incl",title.c_str(),500,-1,1,500,0.,0.3);
    title = string("Incl. - ")+string(tag)+string(";pT (GeV)");
    pT_Incl = new TH1D("pT_Incl",title.c_str(),100,0.,.25);
    snprintf(tag,len,"EMiss>%.1fGeV",dE_cuts[1]);
    title = string("Excl. - ")+string(tag)+string(";#alpha;pT (GeV)");
    am_Excl = new TH2D("am_Excl",title.c_str(),500,-1,1,500,0.,0.3);
    title = string("Excl. - ")+string(tag)+string(";pT (GeV)");
    pT_Excl = new TH1D("pT_Excl",title.c_str(),100,0.,.25);

    snprintf(tag,len,"pT>%.0fMeV,EMiss>%.1fGeV",pT_cuts[2]*1000,dE_cuts[0]);
    title = string("Incl. #phi - ")+string(tag)+string(";#alpha;pT (GeV)");
    am_Iphi = new TH2D("am_Iphi",title.c_str(),500,-1,1,500,0.,0.3);
    title = string("Incl. #phi - ")+string(tag)+string(";ZpV (cm)");
    Z_Iphi =  new TH1D("Z_Iphi", title.c_str(),nZbins,ZMn,ZMx);
    title = string("Incl. #phi - ")+string(tag)+string(";XpV (cm);YpV (cm)");
    XY_Iphi = new TH2D("XY_Iphi",title.c_str(),100,-2.5,2.5,100,-2.5,2.5);
    title = string("Incl. #phi - ")+string(tag)+string(";pT (GeV)");
    pT_Iphi = new TH1D("pT_Iphi",title.c_str(),100,0.,.25);
    title = string("Incl. #phi - ")+string(tag)+string(";EMiss (GeV)");
    dE_Iphi = new TH1D("dE_Iphi",title.c_str(),100,-5,10);
    title = string("Incl. #phi - ")+string(tag)+string(";tracks in pV");
    Tr_Iphi = new TH1D("Tr_Iphi",title.c_str(),32,-.5,31.5);
    title = string("Incl. #phi - ")+string(tag)+string(";pTracks in bottom RICH");
    Rb_Iphi = new TH1D("Rb_Iphi",title.c_str(),32,-.5,31.5);
    title = string("Incl. #phi - ")+string(tag)+string(";pTracks in top RICH");
    Rt_Iphi = new TH1D("Rt_Iphi",title.c_str(),32,-.5,31.5);
    snprintf(tag,len,"pT>%.0fMeV,EMiss<%.1fGeV",pT_cuts[3]*1000,dE_cuts[1]);
    title = string("Excl. #phi - ")+string(tag)+string(";#alpha;pT (GeV)");
    am_Ephi = new TH2D("am_Ephi",title.c_str(),500,-1,1,500,0.,0.3);
    title = string("Excl. #phi - ")+string(tag)+string(";ZpV (cm)");
    Z_Ephi =  new TH1D("Z_Ephi", title.c_str(),nZbins,ZMn,ZMx);
    title = string("Excl. #phi - ")+string(tag)+string(";XpV (cm);YpV (cm)");
    XY_Ephi = new TH2D("XY_Ephi",title.c_str(),100,-2.5,2.5,100,-2.5,2.5);
    title = string("Excl. #phi - ")+string(tag)+string(";pT (GeV)");
    pT_Ephi = new TH1D("pT_Ephi",title.c_str(),100,0.,.25);
    title = string("Excl. #phi - ")+string(tag)+string(";EMiss (GeV)");
    dE_Ephi = new TH1D("dE_Ephi",title.c_str(),100,-5,10);
    title = string("Excl. #phi - ")+string(tag)+string(";tracks in pV");
    Tr_Ephi = new TH1D("Tr_Ephi",title.c_str(),32,-.5,31.5);
    title = string("Excl. #phi - ")+string(tag)+string(";pTracks in bottom RICH");
    Rb_Ephi = new TH1D("Rb_Ephi",title.c_str(),32,-.5,31.5);
    title = string("Excl. #phi - ")+string(tag)+string(";pTracks in top RICH");
    Rt_Ephi = new TH1D("Rt_Ephi",title.c_str(),32,-.5,31.5);
  }
}
void writeKineHistos(const char *particleName)
{
  if      (!strncmp(particleName,"K0",2)) {
    am_all->Write(); Z_all->Write();  XY_all->Write();
    am_K0->Write();  Z_K0->Write();   XY_K0->Write();
    am_K0p->Write(); am_K0m->Write();
    DdD_K0->Write(); cth_K0->Write(); pT_K0->Write();
    Tr_K0->Write();  Rb_K0->Write();  Rt_K0->Write();
  }
  else if (!strncmp(particleName,"Lambda",6)) {
    am_all->Write(); Z_all->Write();  XY_all->Write();
    am_L->Write();   Z_L->Write();    XY_L->Write();
    DdD_L->Write();  cth_L->Write();  pT_L->Write();
    Tr_L->Write();   Rb_L->Write();   Rt_L->Write();
  }
  else if (!strncmp(particleName,"Iphi",6)) {
    pT_Incl->Write(); dE_Incl->Write();
    am_Incl->Write();
    pT_Iphi->Write(); dE_Iphi->Write();
    am_Iphi->Write(); Z_Iphi->Write(); XY_Iphi->Write();
    Tr_Iphi->Write(); Rb_Iphi->Write(); Rt_Iphi->Write();
  }
  else if (!strncmp(particleName,"Ephi",6)) {
    pT_Excl->Write(); dE_Excl->Write();
    am_Excl->Write();
    pT_Ephi->Write(); dE_Ephi->Write();
    am_Ephi->Write(); Z_Ephi->Write(); XY_Ephi->Write();
    Tr_Ephi->Write(); Rb_Ephi->Write(); Rt_Ephi->Write();
  }

}

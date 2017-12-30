// Multichannel Double Recursive Frequentist-Bayesian Particle Identification (pion, kaon, proton, electron etc.)
//
// COMPILE with: root analyscep.cc+ -b -q
//
// TODO: write down a C++ class of this.
//
//
// mikael.mieskolainen@cern.ch, 2017


// C++
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>

// ROOT
#include "TFile.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TProfile.h"
#include "TF1.h"
#include "TObject.h"
#include "TSystem.h"
#include "TLegend.h"


// Set "nice" 2D-plot style
// Read here more about problems with the Rainbow
void set_plot_style() {

    // Set Smooth color gradients
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);

    // Black-Red palette
    gStyle->SetPalette(53); // 56 for inverted

    gStyle->SetTitleOffset(1.6,"x");  //X-axis title offset from axis
    gStyle->SetTitleOffset(1.6,"y");  //Y-axis title offset from axis
    gStyle->SetTitleSize(0.03,"x");   //X-axis title size
    gStyle->SetTitleSize(0.03,"y");   //Y-axis
    gStyle->SetTitleSize(0.03,"z");
    gStyle->SetLabelOffset(0.025);
}

void SetHStyle(TH1F*& h1) {
  h1->SetLineColor(kBlack);
  h1->SetMarkerColor(kBlack);
  h1->SetMarkerStyle(kFullCircle);
  h1->SetMarkerSize(0.2);
}

// Global Style Setup
void setROOTstyle() {

  gStyle->SetOptStat(1); // Statistics BOX OFF [0,1]
  gStyle->SetTitleSize(0.0475,"t"); // Title with "t" (or anything else than xyz)
  gStyle->SetStatY(1.0);
  gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.09);

  // See below
  set_plot_style();
}

// PID vectors
std::vector<std::vector<std::vector<double>>> nsigma;
std::vector<std::vector<double>> posterior;
std::vector<std::vector<double>> particleweight;

// Text Labels
std::vector<TString> pname_a;
std::vector<TString> pname_b;

std::vector<TString> slabels;
std::vector<TString> alabels;
std::vector<TString> blabels;
std::vector<TString> plabels;


// Final state setup
const double MPI = 0.13957018;
const double MK  = 0.493677;
const double MP  = 0.9382720813;
const double MASS[3] = {MPI, MK, MP};

const int NFINALSTATES = 2; // number of track
const int NSPECIES = 3;     // number of particle species
const int NCHANNEL = 9;     // number of 2-body channels, NSPECIES^NFINALSTATE

// Fixed indexing
const int Pi_ind = 0;
const int Ka_ind = 1;
const int Pr_ind = 2;

const int TPC_ind = 0;
const int TOF_ind = 1;

// Phase I PID binning setup
const double    PMIN   = 0.0;
const double    PMAX   = 3.5;
const int       PBINS  = 25;
const int       ITER   = 10;

// Phase II system binning setup
const double    SPMIN   = 0.0;
const double    SPMAX   = 2.0;
const int       SPBINS  = 15;
const int       SITER   = 10;

// Phase III mode
// 0 for Hard (Maximum Probability) or 1 for Soft (weighted)
int PROBmode = 0;


// ROOT Tree input
TFile* f;
TTree* tree2track;

Int_t   run = 0;
Float_t zVtx = 0;
Float_t px1 = 0;
Float_t py1 = 0;
Float_t pz1 = 0;
Float_t px2 = 0;
Float_t py2 = 0;
Float_t pz2 = 0;
Float_t sigmaPiTPC1 = 0;
Float_t sigmaKaTPC1 = 0;
Float_t sigmaPrTPC1 = 0;
Float_t sigmaPiTPC2 = 0;
Float_t sigmaKaTPC2 = 0;
Float_t sigmaPrTPC2 = 0;
Float_t sigmaPiTOF1 = 0;
Float_t sigmaKaTOF1 = 0;
Float_t sigmaPrTOF1 = 0;
Float_t sigmaPiTOF2 = 0;
Float_t sigmaKaTOF2 = 0;
Float_t sigmaPrTOF2 = 0;
Float_t pxMc1 = 0;
Float_t pyMc1 = 0;
Float_t pzMc1 = 0;
Float_t pxMc2 = 0;
Float_t pyMc2 = 0;
Float_t pzMc2 = 0;
Int_t   pidCode1 = 0;
Int_t   pidCode2 = 0;


// Histograms
std::vector<TH1F*> h1M(NCHANNEL, 0);
std::vector<TH1F*> h1Pt(NCHANNEL, 0);
std::vector<TH1F*> h1pt(NCHANNEL, 0);
std::vector<TH1F*> h1Y(NCHANNEL, 0);
std::vector<TH1F*> h1y(NCHANNEL, 0);
std::vector<TH1F*> h1eta(NCHANNEL, 0);
std::vector<TH1F*> h1dy(NCHANNEL, 0);

std::vector<TProfile*> hprMPt(NCHANNEL, 0);
std::vector<TProfile*> hprMpt(NCHANNEL, 0);

std::vector<TH2F*> h2MPt(NCHANNEL, 0);
std::vector<TH2F*> h2Mdphi(NCHANNEL, 0);
std::vector<TH2F*> h2Mpt(NCHANNEL, 0);
std::vector<TH2F*> h2Ptdphi(NCHANNEL, 0);
std::vector<TH2F*> h2ptdphi(NCHANNEL, 0);

std::vector<TH1F*> h1XF(NCHANNEL, 0);
std::vector<TH2F*> h2MXF(NCHANNEL, 0);
std::vector<TProfile*> hprMXF(NCHANNEL, 0);


std::vector<TH1F*> h1XT(NCHANNEL, 0);
std::vector<TH2F*> h2MXT(NCHANNEL, 0);
std::vector<TProfile*> hprMXT(NCHANNEL, 0);


std::vector<TH2F*> h2EtaPhi(NCHANNEL, 0);
std::vector<TH2F*> h2MCosTheta(NCHANNEL, 0);

std::vector<TH2F*> h2TPC_Pi(NCHANNEL, 0);
std::vector<TH2F*> h2TPC_Ka(NCHANNEL, 0);
std::vector<TH2F*> h2TPC_Pr(NCHANNEL, 0);

std::vector<TH2F*> h2TOF_Pi(NCHANNEL, 0);
std::vector<TH2F*> h2TOF_Ka(NCHANNEL, 0);
std::vector<TH2F*> h2TOF_Pr(NCHANNEL, 0);

TProfile* hPl[NCHANNEL][8];


// Track 4-momentum
std::vector<TLorentzVector> p_g; // generated
std::vector<TLorentzVector> p_r; // reconstructed

// System
TLorentzVector system_g;
TLorentzVector system_r;

// SANITY CUTS
const double TPCCUT = 6.0; // sigmas
const double TOFCUT = 1e6;

// Constants
const double z95 = 1.96; // 95% Gaussian CL
const double PI = 3.14159265359; // pi

// ASCII output
FILE* asciif;


// ----------------------------------------------------------------------
// Function prototypes

std::vector<std::vector<double>> EM(double PMIN, double PMAX, int PBINS, int ITER);
int    GetIdx(double value, double MINVAL, double MAXVAL, int NUMBINS);
void   GetProb(std::vector<double>& posterior, const std::vector<std::vector<double>>& nsigma, const std::vector<double>& prior);
double fG(double n);
void   InitTree();
void   CloseTree();
void   GetSignal();
bool   CheckPIDsignal();
void   InitHistogram();
void   FillHisto(std::vector<double> prob_);
void   MakePlots();
void   Analyzer(std::vector<std::vector<double>>& channel_prior, int mode);
void   HEframe(std::vector<TLorentzVector>& p);

bool DEBUG = false;



// MAIN FUNCTION
int analyscep() {

  // Open ASCII output
  asciif = fopen("pipiPIDLHC16.csv", "w");
  if (asciif == NULL) {
    printf("Error opening ascii output file!\n");
    exit(1);
  }
  // --------------------------------------------------------------------------
  // 0. Init vectors

  InitHistogram();

  std::vector<std::vector<double>> ttemp = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  posterior = ttemp;

  // [N final states] x [P particle species] x [D detectors]

  std::vector<std::vector<double>> tempvec = {{0.0, 0.0},{0.0, 0.0},{0.0, 0.0}};

  TLorentzVector temp4vec;

  for (int i = 0; i < NFINALSTATES; ++i) {
    nsigma.push_back(tempvec);
    p_g.push_back(temp4vec);
    p_r.push_back(temp4vec);
  }

  printf("Running... \n");

  // --------------------------------------------------------------------------
  // 1. Run Phase I analysis
  particleweight = EM(PMIN, PMAX, PBINS, ITER);

  // --------------------------------------------------------------------------
  // 2. Run Phase II analysis

  std::vector<std::vector<double>> channel_prior(SPBINS, std::vector<double> (NCHANNEL, 1.0/(double)(NCHANNEL)));

  // EM iterations
  int mode = 1;
  for (int i = 0; i < SITER; ++i) {
    Analyzer(channel_prior, mode);
  }
  
  // Save final histograms
  mode = 2;
  Analyzer(channel_prior, mode);

  // --------------------------------------------------------------------------
  // 3. Make plots
  setROOTstyle();
  MakePlots();


  return EXIT_SUCCESS;
}


// Legendre polynomials
double legendre_pl(int l, double x) {

  if (l == 0) {
    return 1.0;
  } else if (l == 1) {
    return x;

  } else if (l == 2) {
    return 0.5*(3.0*x*x - 1.0);

  } else if (l == 3) {
    return 0.5*(5.0*x*x*x - 3.0*x);

  } else if (l == 4) {
    return (35.0*x*x*x*x - 30.0*x*x + 3.0) / 8.0;

  } else if (l == 5) {
    
    return (63.0 * x*x*x*x*x - 70*x*x*x + 15 * x)/ 8.0;
  } else if (l == 6) {
    
    return (231.0 * x*x*x*x*x*x - 315.0*x*x*x*x + 105.0*x*x - 5.0)/ 16.0;
  } else if (l == 7) {

    return (429.0 * x*x*x*x*x*x*x - 693.0*x*x*x*x*x + 315.0*x*x*x - 35.0*x)/ 16.0;
  } else if (l == 8) {

    return (6435.0*x*x*x*x*x*x*x*x - 12012.0*x*x*x*x*x*x + 6930.0*x*x*x*x - 1260.0*x*x + 35.0)/ 128.0;
  }
}

// Initialize input ROOT tree
void InitTree() {

f = new TFile("tree2track_CUP13.root");
//  f = new TFile("tree2track_kPipmOrexp.root");
//  f = new TFile("tree2track_kKpkmOrexp.root");

  tree2track = (TTree*) f->Get("tree2track");

  tree2track->SetBranchAddress("run",&run);
  tree2track->SetBranchAddress("pxMc1",&pxMc1);
  tree2track->SetBranchAddress("pyMc1",&pyMc1);
  tree2track->SetBranchAddress("pzMc1",&pzMc1);
  tree2track->SetBranchAddress("pxMc2",&pxMc2);
  tree2track->SetBranchAddress("pyMc2",&pyMc2);
  tree2track->SetBranchAddress("pzMc2",&pzMc2);
  tree2track->SetBranchAddress("pidCode1",&pidCode1);
  tree2track->SetBranchAddress("pidCode2",&pidCode2);
  tree2track->SetBranchAddress("zVtx",&zVtx);
  tree2track->SetBranchAddress("px1",&px1);
  tree2track->SetBranchAddress("py1",&py1);
  tree2track->SetBranchAddress("pz1",&pz1);
  tree2track->SetBranchAddress("px2",&px2);
  tree2track->SetBranchAddress("py2",&py2);
  tree2track->SetBranchAddress("pz2",&pz2);
  tree2track->SetBranchAddress("sigmaPiTPC1",&sigmaPiTPC1);
  tree2track->SetBranchAddress("sigmaKaTPC1",&sigmaKaTPC1);
  tree2track->SetBranchAddress("sigmaPrTPC1",&sigmaPrTPC1);
  tree2track->SetBranchAddress("sigmaPiTPC2",&sigmaPiTPC2);
  tree2track->SetBranchAddress("sigmaKaTPC2",&sigmaKaTPC2);
  tree2track->SetBranchAddress("sigmaPrTPC2",&sigmaPrTPC2);
  tree2track->SetBranchAddress("sigmaPiTOF1",&sigmaPiTOF1);
  tree2track->SetBranchAddress("sigmaKaTOF1",&sigmaKaTOF1);
  tree2track->SetBranchAddress("sigmaPrTOF1",&sigmaPrTOF1);
  tree2track->SetBranchAddress("sigmaPiTOF2",&sigmaPiTOF2);
  tree2track->SetBranchAddress("sigmaKaTOF2",&sigmaKaTOF2);
  tree2track->SetBranchAddress("sigmaPrTOF2",&sigmaPrTOF2);

}

void CloseTree() {
  delete f;
}


// Set detector PID signals
void GetSignal() {

  nsigma.at(0).at(Pi_ind).at(TPC_ind) = sigmaPiTPC1;
  nsigma.at(0).at(Ka_ind).at(TPC_ind) = sigmaKaTPC1;
  nsigma.at(0).at(Pr_ind).at(TPC_ind) = sigmaPrTPC1;

  nsigma.at(0).at(Pi_ind).at(TOF_ind) = sigmaPiTOF1;
  nsigma.at(0).at(Ka_ind).at(TOF_ind) = sigmaKaTOF1;
  nsigma.at(0).at(Pr_ind).at(TOF_ind) = sigmaPrTOF1;

  nsigma.at(1).at(Pi_ind).at(TPC_ind) = sigmaPiTPC2;
  nsigma.at(1).at(Ka_ind).at(TPC_ind) = sigmaKaTPC2;
  nsigma.at(1).at(Pr_ind).at(TPC_ind) = sigmaPrTPC2;

  nsigma.at(1).at(Pi_ind).at(TOF_ind) = sigmaPiTOF2;
  nsigma.at(1).at(Ka_ind).at(TOF_ind) = sigmaKaTOF2;
  nsigma.at(1).at(Pr_ind).at(TOF_ind) = sigmaPrTOF2;
}


// PID signal sanity checks
bool CheckPIDsignal() {
  if (std::abs(sigmaPiTPC1) < TPCCUT ||  std::abs(sigmaKaTPC1) < TPCCUT || std::abs(sigmaPrTPC1) < TPCCUT) {
    return true;
  } else {
    return false;
  }
  if (std::abs(sigmaPiTPC2) < TPCCUT ||  std::abs(sigmaKaTPC2) < TPCCUT || std::abs(sigmaPrTPC2) < TPCCUT) {
    return true;
  } else {
    return false;
  }
  /*
  if (std::abs(sigmaPiTOF1) < TOFCUT ||  std::abs(sigmaKaTOF1) < TOFCUT || std::abs(sigmaPrTOF1) < TOFCUT) {
    return true;
  } else {
    return false;
  }
  if (std::abs(sigmaPiTOF2) < TOFCUT ||  std::abs(sigmaKaTOF2) < TOFCUT || std::abs(sigmaPrTOF2) < TOFCUT) {
    return true;
  } else {
    return false;
  }*/
}


// Initialize histograms
void InitHistogram() {

  printf("InitHistogram:: \n");

  int    BINS  = 0;

  // Mass limits
  double MIN_M = 0;
  double MAX_M = 3.5;

  // Rapidity limits
  double MIN_Y = -1.25;
  double MAX_Y = 1.25;

  BINS = 250; MIN_M = 0.0;
  

  pname_a.push_back("#pi^{+}"); pname_a.push_back("K^{+}"); pname_a.push_back("p");
  pname_b.push_back("#pi^{-}"); pname_b.push_back("K^{-}"); pname_b.push_back("#bar{p}");
  plabels.push_back("Pi"); plabels.push_back("Ka"); plabels.push_back("Pr");

  for (int i = 0; i < NSPECIES; ++i) {
    for (int j = 0; j < NSPECIES; ++j) {
      TString label = pname_a[i] + pname_b[j];

      slabels.push_back(label);
      alabels.push_back(pname_a.at(i));
      blabels.push_back(pname_b.at(j));
    }
  }

  for (int k = 0; k < NCHANNEL; ++k) {

    h1M[k]   = new TH1F("M(%d)",  Form(";M(%s) (GeV); Events / (%0.3f GeV)", slabels.at(k).Data(), (MAX_M - MIN_M)/(double)BINS), BINS, MIN_M, MAX_M);
    h1Pt[k]  = new TH1F("Pt(%d)", Form(";System p_{t}(%s) (GeV); Events / (%0.3f GeV)", slabels.at(k).Data(), (MAX_M - MIN_M)/(double)BINS), BINS, MIN_M, MAX_M);
    h1pt[k]  = new TH1F("pt(%d)", Form(";Particle p_{t}(%s) (GeV); Events / (%0.3f GeV)", alabels.at(k).Data(), (MAX_M - MIN_M)/(double)BINS), BINS, MIN_M, MAX_M);

    h1Y[k]   = new TH1F("Y(%d)", Form(";Rapidity Y(%s); Events / %0.3f", slabels.at(k).Data(), (MAX_Y - MIN_Y)/(double)BINS), BINS, MIN_Y, MAX_Y);
    h1y[k]   = new TH1F("y(%d)", Form(";Rapidity y(%s); Events / %0.3f", alabels.at(k).Data(), (MAX_Y - MIN_Y)/(double)BINS), BINS, MIN_Y, MAX_Y);
    h1eta[k] = new TH1F("eta(%d)", Form(";Pseudorapidity #eta(%s); Events / %0.3f", alabels.at(k).Data(), (MAX_Y - MIN_Y)/(double)BINS), BINS, MIN_Y, MAX_Y);
    h1dy[k]  = new TH1F("y(%d)", Form(";Pair #Deltay(%s); Events / %0.3f", slabels.at(k).Data(), 2*(MAX_Y - MIN_Y)/(double)BINS), BINS, 2*MIN_Y, 2*MAX_Y);
    
    hprMPt[k]  = new TProfile("prof MPt(%d)", Form(";M(%s) (GeV); System <p_{t}> (GeV)", slabels.at(k).Data()), BINS, MIN_M, MAX_M);
    hprMpt[k]  = new TProfile("prof Mpt(%d)", Form(";M(%s) (GeV); Particle <p_{t}(%s)> (GeV)", slabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_M, MAX_M);

    h2EtaPhi[k]   = new TH2F("eta,phi", Form(";Particle #eta(%s); Particle #phi(%s) (rad)", alabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_Y, MAX_Y, BINS, -PI, PI);
    h2MCosTheta[k] = new TH2F("M,costheta", Form(";M(%s) (GeV); cos #theta(%s)_{ r.f.} (rad)", slabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_M, MAX_M, BINS, -1, 1);

    //h1M_pipi->Sumw2();
    
    h2MPt[k]    = new TH2F("M,Pt()", Form(";M(%s) (GeV); System p_{T} (GeV)", slabels.at(k).Data()), BINS, MIN_M, MAX_M, BINS, 0, 2.5);
    h2Mdphi[k]  = new TH2F("M,dphi()", Form(";M(%s) (GeV); #Delta#phi (rad)", slabels.at(k).Data()), BINS, MIN_M, MAX_M, BINS, 0, 3.14159);
    h2Mpt[k]    = new TH2F("M,pt()", Form(";M(%s) (GeV); %s p_{T} (GeV)", slabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_M, MAX_M, BINS, 0, 2.5);
    h2Ptdphi[k] = new TH2F("Pt,dphi()", Form(";System p_{T}(%s) (GeV); #Delta#phi (rad)", slabels.at(k).Data()), BINS, 0.0, 2.5, BINS, 0, 3.14159);
    h2ptdphi[k] = new TH2F("pt,dphi()", Form(";Particle p_{T}(%s) (GeV); #Delta#phi (rad)", alabels.at(k).Data()), BINS, 0.0, 2.5, BINS, 0, 3.14159);


    // Feynman xF
    h1XF[k]    = new TH1F("x_F ", Form(";Feynman x_{F}(%s);Events", alabels.at(k).Data()), BINS, 0, 1);
    h2MXF[k]   = new TH2F("M x_F ", Form(";M(%s) (GeV); Feynman x_{F}(%s)", slabels.at(k).Data(), alabels.at(k).Data()), BINS, 0, 4, BINS, 0, 1);
    hprMXF[k]  = new TProfile("prof MXF(%d)", Form(";M(%s) (GeV); Feynman <x_{F}(%s)>", slabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_M, MAX_M);

    // xT
    h1XT[k]    = new TH1F("x_T ", Form(";x_{T}(%s);Events", slabels.at(k).Data()), BINS, 0, 1);
    h2MXT[k]   = new TH2F("M x_T ", Form(";M(%s) (GeV); x_{T}(%s)", slabels.at(k).Data(), alabels.at(k).Data()), BINS, 0, 4, BINS, 0, 1);
    hprMXT[k]  = new TProfile("prof MXT(%d)", Form(";M(%s) (GeV); <x_{T}(%s)>", slabels.at(k).Data(), alabels.at(k).Data()), BINS, MIN_M, MAX_M);

    // Save errors
    hprMPt[k]->Sumw2();
    hprMpt[k]->Sumw2();
    hprMXF[k]->Sumw2();
    hprMXT[k]->Sumw2();


    h2TPC_Pi[k] = new TH2F("(p,TPC_Pi)", ";p (GeV); #sigma #pi TPC", 200, 0, 3, 200, -50, 50);
    h2TPC_Ka[k] = new TH2F("(p,TPC_Ka)", ";p (GeV); #sigma K TPC",   200, 0, 3, 200, -50, 50);
    h2TPC_Pr[k] = new TH2F("(p,TPC_Pr)", ";p (GeV); #sigma p TPC",   200, 0, 3, 200, -50, 50);

    h2TOF_Pi[k] = new TH2F("(p,TPC_Pi)", ";p (GeV); #sigma #pi TOF", 200, 0, 3, 200, -50, 50);
    h2TOF_Ka[k] = new TH2F("(p,TPC_Ka)", ";p (GeV); #sigma K TOF",   200, 0, 3, 200, -50, 50);
    h2TOF_Pr[k] = new TH2F("(p,TPC_Pr)", ";p (GeV); #sigma p TOF",   200, 0, 3, 200, -50, 50);

    // Legendre polynomials, DO NOT CHANGE THE Y-RANGE [-1,1]
    for (int i = 0; i < 8; ++i) {
      hPl[k][i] = new TProfile(Form("hPl%d", i+1),"", 100, 0.0, MAX_M, -1, 1);

      hPl[k][i]->SetXTitle(Form("M(%s) (GeV)", slabels.at(k).Data())); 
      hPl[k][i]->SetYTitle(Form("#LTP_{l}(cos(#theta)#GT |_{ r.f.}"));
    }

  }
}

// Fill histograms
void FillHisto(std::vector<double> prob_) {

  // Find the maximum probability channel
  uint k_max = 999;
  double P_max = -1;
  for (uint i = 0; i < prob_.size(); ++i) {
    if (prob_.at(i) > P_max) {
      P_max = prob_.at(i);
      k_max = i;
    }
  }


  // Fill hard classification vector [0,0,1,0,..,,0]
  std::vector<double> weight_(NCHANNEL, 0.0);

  // Hard
  if (PROBmode == 0) {
    for (uint k = 0; k < weight_.size(); ++k) {
      if (k == k_max) {
        weight_.at(k) = 1.0;
      }
    }
  }
  // Weighted
  if (PROBmode == 1) {
    weight_ = prob_;
  }

  // Hard classification for each k-th channel
  int k = 0;
  for (int i = 0; i < NSPECIES; ++i) {
    for (int j = 0; j < NSPECIES; ++j) {

      // NOTE HERE THAT one must NOT fill events with zero weights in the hard
      // mode, otherwise the statistical errors are not properly calculated.
      if (PROBmode == 0 && weight_.at(k) > 1e-5 || PROBmode == 1) {

        // Reconstruct final state 4-momentum
        p_r.at(0).SetVectM(p_r.at(0).Vect(), MASS[i]);
        p_r.at(1).SetVectM(p_r.at(1).Vect(), MASS[j]);
        system_r = p_r.at(0) + p_r.at(1);

        if (k == 0) { // Write down ascii output for channel[0]
          fprintf(asciif, "%0.6f,%0.6f,%0.6f,%0.6f \n", system_r.M(), system_r.Px(), system_r.Py(), system_r.Pz());
        }
        
        // Observables
        double M   = system_r.M();
        double Pt  = system_r.Perp();
        double Y   = system_r.Rapidity();
        double pt  = p_r.at(0).Perp();
        double y   = p_r.at(0).Rapidity();
        double eta = p_r.at(0).Eta();
        double phi = p_r.at(0).Phi();

        // Fill histograms >>
        h1M[k]->Fill(M, weight_.at(k));
        h1Pt[k]->Fill(Pt, weight_.at(k));
        h1pt[k]->Fill(pt, weight_.at(k));

        h1Y[k]->Fill(Y, weight_.at(k));
        h1y[k]->Fill(y, weight_.at(k));
        h1eta[k]->Fill(eta, weight_.at(k));
        h1dy[k]->Fill(p_r.at(0).Rapidity() - p_r.at(1).Rapidity(), weight_.at(k));

        hprMPt[k]->Fill(M, Pt, weight_.at(k));
        hprMpt[k]->Fill(M, pt, weight_.at(k));

        h2EtaPhi[k]->Fill(eta, phi, weight_.at(k));


        // Calculate Feynman x = 2 |p_z*| / sqrt(shat)
        // Boost particles to the central system rest frame
        TVector3 betavec = -system_r.BoostVector(); // Note the minus sign
        TLorentzVector pboosted = p_r.at(0);
        pboosted.Boost(betavec);
        double xf = 2*std::abs(pboosted.Pz()) / system_r.M();
        double xt = 2*std::abs(pboosted.Perp()) / system_r.M();

        h1XF[k]->Fill(xf, weight_.at(k));
        h2MXF[k]->Fill(M, xf, weight_.at(k));
        hprMXF[k]->Fill(M, xf, weight_.at(k));

        h1XT[k]->Fill(xt, weight_.at(k));
        h2MXT[k]->Fill(M, xt, weight_.at(k));
        hprMXT[k]->Fill(M, xt, weight_.at(k));

        h2MCosTheta[k]->Fill(M, pboosted.CosTheta(), weight_.at(k));
        h2MPt[k]->Fill(M, Pt, weight_.at(k));
        h2Mdphi[k]->Fill(M, p_r.at(0).DeltaPhi(p_r.at(1)), weight_.at(k));

        h2Mpt[k]->Fill(M, p_r.at(0).Perp(), weight_.at(k));
        h2Mpt[k]->Fill(M, p_r.at(1).Perp(), weight_.at(k));
        
        h2Ptdphi[k]->Fill(Pt, p_r.at(0).DeltaPhi(p_r.at(1)), weight_.at(k));
        h2ptdphi[k]->Fill(pt, p_r.at(0).DeltaPhi(p_r.at(1)), weight_.at(k));

        // Calculate cos(theta) in the non-rotated rest frame
        double costheta = pboosted.CosTheta();

        std::vector<TLorentzVector> fs;
        fs.push_back(p_r.at(0)); fs.push_back(p_r.at(1));

        HEframe(fs);

        double HEcostheta = fs.at(0).CosTheta();

        // Legendre polynomials P_l cos(theta), l = 1,2,3,4,5,6,7,8
        for (int l = 0; l < 8; ++l) { // note l+1
          double value = legendre_pl((l+1), HEcostheta); // cos(theta)
          hPl[k][l]->Fill(system_r.M(), value);
        }
      }
      ++k;
    }
  }
}

// From lab to the Helicity frame
// Quantization z-axis as the direction of the resonance in the lab frame
//
// Input is a vector of final state 4-momentum
void HEframe(std::vector<TLorentzVector>& p) {

    // Sum to get the system 4-momentum
    TLorentzVector X(0,0,0,0);
    for (UInt_t i = 0; i < p.size(); ++i) {
        X += p.at(i);
    }

        // ********************************************************************
        if (DEBUG) {
        printf("\n\n ::HELICITY FRAME:: \n");
        printf("- Pions in LAB FRAME: \n");
        p.at(0).Print();
        p.at(1).Print();
        }
        // ********************************************************************

    // ACTIVE (cf. PASSIVE) rotation of final states by z-y-z (phi,theta,-phi)
    // (in the opposite direction -> minus signs) Euler angle sequence.
    double Z_angle = - X.Phi();
    double Y_angle = - X.Theta();

    for (UInt_t i = 0; i < p.size(); ++i) {
        p.at(i).RotateZ(Z_angle);
        p.at(i).RotateY(Y_angle);
    //    p.at(i).RotateZ(-Z_angle);    // Comment this one out for tests
    }

        // ********************************************************************
        if (DEBUG) {

        printf("- Pions in ROTATED LAB FRAME: \n");
        p.at(0).Print();
        p.at(1).Print();
        
        TVector3 ex(1,0,0);
        TVector3 ey(0,1,0);
        TVector3 ez(0,0,1);

        // x -> x'
        ex.RotateZ(Z_angle);
        ex.RotateY(Y_angle);
     //   ex.RotateZ(-Z_angle);  // Comment this one out for tests
        // y -> y'
        ey.RotateZ(Z_angle);
        ey.RotateY(Y_angle);
   //     ey.RotateZ(-Z_angle);  // Comment this one out for tests
        // z -> z'
        ez.RotateZ(Z_angle);
        ez.RotateY(Y_angle);
 //       ez.RotateZ(-Z_angle);  // Comment this one out for tests

        printf("- AXIS vectors after rotation: \n");

        ex.Print();
        ey.Print();
        ez.Print();
        }
        // ********************************************************************

    // Construct the central system 4-momentum in ROTATED FRAME
    TLorentzVector XNEW(0,0,0,0);
    for (UInt_t i = 0; i < p.size(); ++i) {
        XNEW += p.at(i);
    }

    // Boost particles to the central system rest frame
    // -> Helicity frame obtained
    TVector3 betavec = -XNEW.BoostVector(); // Note the minus sign
    for (UInt_t i = 0; i < p.size(); ++i) {
        p.at(i).Boost(betavec);
    }

        // ********************************************************************
        if (DEBUG) {

        printf("- Pions after boost in HELICITY FRAME: \n");
        p.at(0).Print();
        p.at(1).Print();

        printf("- Central system in LAB FRAME: \n");
        X.Print();

        printf("- Central system in ROTATED LAB FRAME: \n");
        XNEW.Print();

        printf("\n");
        }
        // ********************************************************************

    TLorentzVector sum;
    for (UInt_t i = 0; i < p.size(); ++i) {
        sum += p.at(i);
    }
    Double_t epsilon = 1e-6;
    if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon || std::abs(sum.Pz()) > epsilon) {
        printf("HEframe:: Not a rest frame!! \n");
    }
}

// Create plots
void MakePlots() {

  printf("MakePlots:: \n");

  // Title
  //h1M_pipi->SetTitle(Form("#Sigma N = %0.0f / %d / %d", channelweight[0], N_selected, N_tot));

  for (int k = 0; k < NCHANNEL; ++k) {

    // Create output directory in a case
    TString output_dir = "./figs/" + slabels.at(k);
    gSystem->Exec(Form("mkdir %s", output_dir.Data()));

    TCanvas c10("c10", "c10", 400, 300);
    TCanvas c10_log("c10_log", "c10_log", 400, 300);
    c10_log.SetLogy();

    c10.cd();
    h1M[k]->Draw();    
    c10.SaveAs(Form("./figs/%s/h1M_%d.pdf", slabels.at(k).Data(), k));
    c10_log.cd(); 
    h1M[k]->Draw(); h1M[k]->SetMinimum(0.1);  //   Y-axis minimum (for log only)
    c10_log.SaveAs(Form("./figs/%s/h1M_logy_%d.pdf", slabels.at(k).Data(), k));

    c10.cd();
    h1Pt[k]->Draw();
    c10.SaveAs(Form("./figs/%s/h1Pt_%d.pdf", slabels.at(k).Data(), k));
    c10_log.cd();
    h1Pt[k]->Draw();  h1Pt[k]->SetMinimum(0.1);  //   Y-axis minimum (for log only)
    c10_log.SaveAs(Form("./figs/%s/h1Pt_logy_%d.pdf", slabels.at(k).Data(), k));

    c10.cd();
    h1pt[k]->Draw();
    c10.SaveAs(Form("./figs/%s/h1pt_%d.pdf", slabels.at(k).Data(), k));
    c10_log.cd();
    h1pt[k]->Draw(); h1pt[k]->SetMinimum(0.1);  //   Y-axis minimum (for log only)
    c10_log.SaveAs(Form("./figs/%s/h1pt_logy_%d.pdf", slabels.at(k).Data(), k));

    c10.cd();
    h1Y[k]->Draw();    c10.SaveAs(Form("./figs/%s/h1Y_%d.pdf", slabels.at(k).Data(), k));
    h1y[k]->Draw();    c10.SaveAs(Form("./figs/%s/h1y_%d.pdf", slabels.at(k).Data(), k));
    h1eta[k]->Draw();  c10.SaveAs(Form("./figs/%s/h1eta_%d.pdf", slabels.at(k).Data(), k));
    h1dy[k]->Draw();   c10.SaveAs(Form("./figs/%s/h1dy_%d.pdf", slabels.at(k).Data(), k));

    hprMPt[k]->Draw(); c10.SaveAs(Form("./figs/%s/hprMPt_%d.pdf", slabels.at(k).Data(), k));
    hprMpt[k]->Draw(); c10.SaveAs(Form("./figs/%s/hprMpt_%d.pdf", slabels.at(k).Data(), k));


    TCanvas c100("c100", "c100", 400, 300);
    c100.cd();
    h2EtaPhi[k]->Draw("COLZ");    c100.SaveAs(Form("./figs/%s/h2EtaPhi_%d.pdf", slabels.at(k).Data(), k));
    c100.cd();
    h2MCosTheta[k]->Draw("COLZ"); c100.SaveAs(Form("./figs/%s/h2MCosTheta_%d.pdf", slabels.at(k).Data(), k));


    c100.cd();
    h2MPt[k]->Draw("COLZ");    c100.SaveAs(Form("./figs/%s/h2MPt_%d.pdf", slabels.at(k).Data(), k));
    c100.cd();
    h2Mdphi[k]->Draw("COLZ");  c100.SaveAs(Form("./figs/%s/h2Mdphi_%d.pdf", slabels.at(k).Data(), k));
    c100.cd();
    h2Mpt[k]->Draw("COLZ");    c100.SaveAs(Form("./figs/%s/h2Mpt_%d.pdf", slabels.at(k).Data(), k));
    c100.cd();
    h2Ptdphi[k]->Draw("COLZ"); c100.SaveAs(Form("./figs/%s/h2Ptdphi_%d.pdf", slabels.at(k).Data(), k));
    c100.cd();
    h2ptdphi[k]->Draw("COLZ"); c100.SaveAs(Form("./figs/%s/h2ptdphi_%d.pdf", slabels.at(k).Data(), k));

    // -------------------------------------------------------------------------------------
    TCanvas cX("cX", "cX", 400, 300);
    cX.cd();
    h1XF[k]->Draw();           cX.SaveAs(Form("./figs/%s/h1XF_%d.pdf", slabels.at(k).Data(), k));
    TCanvas cX2("cX2", "cX2", 400, 300);
    cX2.cd();
    h2MXF[k]->Draw("COLZ");    cX2.SaveAs(Form("./figs/%s/h2MXF_%d.pdf", slabels.at(k).Data(), k));
    TCanvas cX3("cX3", "cX3", 400, 300);
    cX3.cd();
    hprMXF[k]->Draw();         cX3.SaveAs(Form("./figs/%s/hprMXF_%d.pdf", slabels.at(k).Data(), k));


    cX.cd();
    h1XT[k]->Draw();           cX.SaveAs(Form("./figs/%s/h1XT_%d.pdf", slabels.at(k).Data(), k));
    cX2.cd();
    h2MXT[k]->Draw("COLZ");    cX2.SaveAs(Form("./figs/%s/h2MXT_%d.pdf", slabels.at(k).Data(), k));
    cX3.cd();
    hprMXT[k]->Draw();         cX3.SaveAs(Form("./figs/%s/hprMXT_%d.pdf", slabels.at(k).Data(), k));


    // -------------------------------------------------------------------------------------
    // Legendre polynomials in a Rest Frame

    const Int_t colors[4] = {48, 53, 98, 32};

    // 1...4
    {
      TCanvas* c115 = new TCanvas("c115","Legendre polynomials",600,400);
      TLegend* leg[4];
      c115->Divide(2,2, 0.001, 0.001);

      for (Int_t l = 0; l < 4; ++l) { 
        c115->cd(l+1); // note (l+1)

        leg[l] = new TLegend(0.15,0.75,0.4,0.85); // x1,y1,x2,y2
        hPl[k][l]->SetLineColor(colors[l]);
        hPl[k][l]->Draw();
        hPl[k][l]->SetMinimum(-0.4); // Y-axis minimum
        hPl[k][l]->SetMaximum( 0.4); // Y-axis maximum

        leg[l]->SetFillColor(0);  // White background
        leg[l]->SetBorderSize(0); // No box
        leg[l]->AddEntry(hPl[k][l], Form("l = %d", l+1), "l");
        leg[l]->Draw();
      }
      c115->SaveAs(Form("./figs/%s/hPl_1to4_%d.pdf", slabels.at(k).Data(), k));
    }

    // 5...8
    {
      TCanvas* c115 = new TCanvas("c115","Legendre polynomials",600,400);
      TLegend* leg[4];
      c115->Divide(2,2, 0.001, 0.001);

      const Int_t colors[4] = {48, 53, 98, 32};
      for (Int_t l = 4; l < 8; ++l) { 
        c115->cd(l-3); // note (l-3)

        leg[l] = new TLegend(0.15,0.75,0.4,0.85); // x1,y1,x2,y2
        hPl[k][l]->SetLineColor(colors[l-4]);
        hPl[k][l]->Draw();
        hPl[k][l]->SetMinimum(-0.4); // Y-axis minimum
        hPl[k][l]->SetMaximum( 0.4); // Y-axis maximum

        leg[l]->SetFillColor(0);  // White background
        leg[l]->SetBorderSize(0); // No box
        leg[l]->AddEntry(hPl[k][l], Form("l = %d", l+1), "l");
        leg[l]->Draw();
      }
      c115->SaveAs(Form("./figs/%s/hPl_5to8_%d.pdf", slabels.at(k).Data(), k));
    }

  }

}

// Phase II analyzer
void Analyzer(std::vector<std::vector<double>>& channel_prior, int mode) {

  printf("EM Phase II: \n");

  std::vector<std::vector<double>> channelweight(SPBINS, std::vector<double> (NCHANNEL, 1.0/(double)(NCHANNEL)));

  // Initialize data source
  InitTree();

  std::vector<double> posterior_sum(NCHANNEL, 0.0);

  // Loop over events
  int N_selected = 0;
  for (Int_t ev = 0; ev < tree2track->GetEntries(); ++ev) {
    
    tree2track->GetEntry(ev);

    // Get detector signal
    GetSignal();

    // Check PID signal sanity
    if (!CheckPIDsignal())
      continue;
    ++N_selected;

    // Get momentum of tracks
    p_r.at(0).SetXYZM(px1,py1,pz1,0);
    p_r.at(1).SetXYZM(px2,py2,pz2,0);
    system_r = p_r.at(0) + p_r.at(1);

    // Find out the momentum bin for each final state
    std::vector<int> pbin(NFINALSTATES, 0);
    for (uint f = 0; f < NFINALSTATES; ++f) {
      pbin.at(f) = GetIdx(p_r.at(f).Perp(), PMIN, PMAX, PBINS);
    }

    // Get posteriori probabilities for each final states
    for (int f = 0; f < NFINALSTATES; ++f) {
      GetProb(posterior.at(f), nsigma.at(f), particleweight[pbin.at(f)]);
    }

    // Find out the momentum bin of the system
    int sbin = GetIdx(system_r.Perp(), SPMIN, SPMAX, SPBINS);

    // Probabilities of different decay channels, by factorizing (independence) P_tot = P_1 x P_2
    std::vector<double> prob_(NCHANNEL, 0.0);

    int k = 0;
    for (int i = 0; i < NSPECIES; ++i) {
      for (int j = 0; j < NSPECIES; ++j) {
        prob_[k] = posterior.at(0).at(i) * posterior.at(1).at(j) * channel_prior.at(sbin).at(k); // pi+pi-, KK etc.
        ++k;
      }
    }

    // Normalize the channel probabilities sum to one
    double summ = 1e-12;
    for (uint i = 0; i < prob_.size(); ++i) {
      summ += prob_[i];
    }
    for (uint i = 0; i < prob_.size(); ++i) {
      prob_.at(i) /= summ;
    }

    if (mode == 2) {
      // Histograms
      FillHisto(prob_);
    }

    // Save decay channel weights
    for (uint i = 0; i < prob_.size(); ++i) {
      channelweight.at(sbin).at(i) += prob_.at(i);
    }
  } // Event loop

  // Update priors
  for (int i = 0; i < SPBINS; ++i) {
    for (uint j = 0; j < NCHANNEL; ++j) {
      channel_prior.at(i).at(j) = channelweight.at(i).at(j) / (double) N_selected;
    }
  }

  // Collect per channel probabilities
  for (int c = 0; c < NCHANNEL; ++c) {
    for (int i = 0; i < SPBINS; ++i) {
      posterior_sum.at(c) += channelweight.at(i).at(c);
    }
  }

  // Normalize per bin
  double NTOT = 0.0;
  for (int i = 0; i < SPBINS; ++i) {
    double sum = 1e-12;
    for (int j = 0; j < NCHANNEL; ++j) {
      sum += channelweight[i][j];
    }
    NTOT += sum;
    printf("B[%2d]=[%0.2f,%0.2f] GeV [N = %0.2E]: ", i, (SPMAX - 0.0) / SPBINS * i, (SPMAX - 0.0) / SPBINS * (i+1), sum);

    // Ratios
    for (int j = 0; j < NCHANNEL; ++j) {
        channelweight[i][j] /= sum;
        printf("%0.3f ", channelweight[i][j]);   
    }
    printf(" +- ");
    // Uncertainty
    for (int j = 0; j < NCHANNEL; ++j) {
        printf("%0.1E ", z95 * std::sqrt( 1 / sum * channelweight[i][j] * (1.0 - channelweight[i][j]) ));   
    }
    printf("\n");
  }
  printf("NTOT = %0.2E events \n", NTOT);
  printf("\n");

  // Calculate sum
  double Nsum = 0;
  for (int i = 0; i < NCHANNEL; ++i) {
    Nsum += posterior_sum.at(i);
  }

  // Print out particle ratios with Binomial proportion uncertainty
  printf("EM Phase II: Integrated channel fractions (CL95 statistical): \n");
  for (int i = 0; i < NCHANNEL; ++i) {
    double p = posterior_sum.at(i) / Nsum;
    printf("Channel[%d] = %0.4f +- %0.2E (%s) \n", i, p, z95 * std::sqrt(1/Nsum * p * (1.0 - p)), slabels.at(i).Data() );
  }
  printf("\n");

  // Close our data source
  CloseTree();
}


// Phase I analyzer
std::vector<std::vector<double>> EM(double PMIN, double PMAX, int PBINS, int ITER) {

  printf("EM Phase I: \n");

  // Init data
  InitTree();

  std::vector<double> posterior_sum;
  for (int i = 0; i < NSPECIES; ++i) {
    posterior_sum.push_back(0.0);
  }
  std::vector<std::vector<double> > particleweight(PBINS, std::vector<double>(NSPECIES, 1.0/(double)NSPECIES));

  if (ITER == 0) {
    printf("EM Phase I: no iteration, flat priors returned\n");
    return particleweight;
  }

  for (int it = 0; it < ITER; ++it) {

    printf("EM Phase I: Iteration = %d/%d (95CL statistical)\n", it+1, ITER);

    // Loop over events
    for (Int_t ev = 0; ev < tree2track->GetEntries(); ++ev) {

      tree2track->GetEntry(ev);

      // Get detector signal
      GetSignal();

      // Check PID signal sanity
      if (!CheckPIDsignal())
        continue;

      // Get momentum of tracks
      p_r.at(0).SetXYZM(px1,py1,pz1,0);
      p_r.at(1).SetXYZM(px2,py2,pz2,0);

      // Find out the momentum bin for each final state
      std::vector<int> pbin(NFINALSTATES, 0);
      for (int f = 0; f < NFINALSTATES; ++f) {
        pbin.at(f) = GetIdx(p_r.at(f).Perp(), PMIN, PMAX, PBINS);
      }

      // Get posteriori probabilities for each final states
      for (int f = 0; f < NFINALSTATES; ++f) {
        GetProb(posterior.at(f), nsigma.at(f), particleweight[pbin.at(f)]);
      }

      // Save |momentum| dependent weight for iterative mapping
      for (int i = 0; i < NSPECIES; ++i) {

        // Save weight from each final state
        for (int f = 0; f < NFINALSTATES; ++f) {
          particleweight[pbin.at(f)][i] += posterior.at(f).at(i);
        }
      }

      // LAST ITERATION: Save posteriori probabilities for particle ratios 
      if (it == ITER - 1) {
        for (int i = 0; i < NSPECIES; ++i) {

          // Loop over final states
          for (int f = 0; f < NFINALSTATES; ++f) {
            posterior_sum.at(i) += posterior.at(f).at(i);  
          }
        }
      }
    }

    // Normalize per bin
    double NTOT = 0.0;
    for (int i = 0; i < PBINS; ++i) {
      double sum = 1e-12;
      for (int j = 0; j < NSPECIES; ++j) {
        sum += particleweight[i][j];
      }
      NTOT += sum;
      printf("bin[%2d] = [%0.3f, %0.3f] GeV [N = %0.2E]: ", i, (PMAX - 0.0) / PBINS * i, (PMAX - 0.0) / PBINS * (i+1), sum);

      // Ratios
      for (int j = 0; j < NSPECIES; ++j) {
          particleweight[i][j] /= sum;
          printf("%0.3f ", particleweight[i][j]);   
      }
      printf(" +- ");
      // Uncertainty
      for (int j = 0; j < NSPECIES; ++j) {
          printf("%0.1E ", z95 * std::sqrt( 1 / sum * particleweight[i][j] * (1.0 - particleweight[i][j]) ));   
      }
      printf("\n");
    }
    printf("NTOT = %0.2E particles \n", NTOT);
    printf("\n");
  } // EM iteration

  // Calculate sum
  double Nsum = 0;
  for (int i = 0; i < NSPECIES; ++i) {
    Nsum += posterior_sum.at(i);
  }

  // Print out particle ratios with Binomial proportion uncertainty
  printf("EM Phase I: Integrated particle fractions (CL95 statistical): \n");
  for (int i = 0; i < NSPECIES; ++i) {
    double p = posterior_sum.at(i) / Nsum;
    printf("Particle[%d] = %0.4f +- %0.2E (%s) \n", i, p, z95 * std::sqrt(1/Nsum * p * (1.0 - p)), plabels.at(i).Data() );
  }
  printf("\n");

  // Close our data
  CloseTree();

  return particleweight;
}


// Get table/histogram index for linearly spaced bins
// Gives exact uniform filling within bin boundaries.
int GetIdx(double value, double MINVAL, double MAXVAL, int NUMBINS) {

  const double BINWIDTH = (MAXVAL - MINVAL) / NUMBINS;
  int idx = std::floor((value - MINVAL) / BINWIDTH);

  if (idx < 0) {           // Underflow
    idx = 0;
  }
  if (idx > NUMBINS - 1) { // Overflow
    idx = NUMBINS - 1;
  }
  return idx;
}


// Get posteriori probabilities
void GetProb(std::vector<double>& posterior, 
  const std::vector<std::vector<double>>& nsigma, const std::vector<double>& prior) {

  const double BOUND = 50; // No signal bound +-

  std::vector<double> fval(prior.size(), 1.0); // init with 1.0

  // Evaluate likelihood for different particle types
  for (uint i = 0; i < prior.size(); ++i) {

    // Loop over detectors (TPC,TOF)
    //bool check = false; 
    for (uint det = 0; det < nsigma.at(i).size(); ++det) {

      // Missing values (such as missing TOF are omitted from the product likelihood)
      if (-BOUND < nsigma.at(i).at(det) && nsigma.at(i).at(det) < BOUND) {

        // Independent detectors -> product likelihood (correlation information neglected)
        fval.at(i) *= fG(nsigma.at(i).at(det)); 
      }
    }
  }

  // Evaluate Posterior ~ likelihood x prior
  for (uint i = 0; i < posterior.size(); ++i) {
    posterior.at(i) = fval.at(i) * prior.at(i);
  }

  // Normalize by evidence (denominator of Bayes formula)
  double denominator = 0;
  for (uint i = 0; i < posterior.size(); ++i) {
    denominator += posterior.at(i);
  }
  // Now Posterior = likelihood x prior / sum_i likelihood x prior
  denominator += 1e-12; // protection
  for (uint i = 0; i < posterior.size(); ++i) {
    posterior.at(i) /= denominator;
  }
}

// Gaussian likelihood (density), based on asumption that the generation of nsigma
// values has been normalized along the different particle species (pion, kaon) etc.
// -> can be compared.
double fG(double n) {

  const double sigma = 1.0; 
  return 1.0/(std::sqrt(2.0*PI)*sigma) * std::exp(-0.5*n*n);

}


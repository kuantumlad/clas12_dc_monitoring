#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include "Riostream.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TMath.h"
#include <vector>
#include <TAxis.h>
#include <TLorentzRotation.h>
#include<vector>
#include <algorithm>
#include <functional>
#include <TChain.h>
#include <TCutG.h>
#include "TStyle.h"
#include "TROOT.h"

using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;


float returnPEFromThetaE(float Ebeam, TLorentzVector p4_electron);
float returnPEFromThetaP(double Ebeam, TLorentzVector p4_proton);
float returnPEFromPProton(float Ebeam, TLorentzVector p4_proton);
float returnThetaEFromPE(float Ebeam, TLorentzVector p4_electron);
float returnThetaEFromPProton(float Ebeam, TLorentzVector p4_proton);
float returnThetaEFromThetaProton(float Ebeam, TLorentzVector p4_proton);
float returnThetaPFromThetaE(float Ebeam, TLorentzVector p4_electron);
float returnThetaPFromPE(float Ebeam, TLorentzVector p4_electron);
float returnThetaPFromPP(float Ebeam, TLorentzVector p4_proton);
float returnPPFromThetaE(double Ebeam, TLorentzVector p4_electron);
float returnPPromThetaP(float Ebeam, TLorentzVector p4_proton);
float returnPPromEP(double Ebeam, TLorentzVector p4_electron);
float returnEBeamFromThetaEThetaP(TLorentzVector p4_electron, TLorentzVector p4_proton);
double kin_W(TLorentzVector ele, float Ebeam);
double kin_mismass2(TLorentzVector ele, TLorentzVector hadron, float Ebeam);
double kin_pmismass2(TLorentzVector hadron, float Ebeam);
double kin_Q2(TLorentzVector ele, float Ebeam);
double returnPhi(double phi);

TF1* fitHisto(TH1D *htemp);


int electronProtonFinal(int model){

	gStyle->SetErrorX(0);
	gStyle->SetPadGridY(kTRUE);
    //  gStyle->SetGridStyle(2);
	//gStyle->SetPalette(1  ,0);
	// /gStyle->SetPalette(47);
	gStyle->SetCanvasColor(0);
	gStyle->SetFrameFillColor(0);
	Int_t cancolor = 10;
	gStyle->SetNdivisions(505, "X");
	gStyle->SetNdivisions(505, "Y");
	gStyle->SetPadGridX(kFALSE);
	gStyle->SetPadGridY(kFALSE);

	gStyle->SetLabelFont(42,"x");
	gStyle->SetLabelSize(0.05,"x");
	gStyle->SetLabelFont(42,"y");
	gStyle->SetLabelSize(0.05,"y");
	gStyle->SetOptTitle(0);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.1);
	gStyle->SetOptStat(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->cd();
	gROOT->ForceStyle();

	//histogram parameters

	int nWBins = 250;
	int nQ2Bins = 250;
	float lowBorderW = 0.85;
	float highBorderW = 1.7;
	float lowBorderQ2 = 0.0;
	float highBorderQ2 = 1.0;
	float W, Q2;
	int sectorElectron;
	int sectorProton;

	int nOverlallBins = 250;
	float eBeamLow = 0;
	float eBeamHigh = 2.5;
	float thetaELow = 5;
	float thetaEHigh = 30;
	float thetaPLow = 50;
	float thetaPHigh = 90;

	float pELow = 1.6;
	float pEHigh = 2.2;
	float pPLow = 0.1;
	float pPHigh = 1.0;

	float deltaEenergyLow = -1.2;
	float deltaEenergyHigh = 1.2;

	float deltaThetaELow = -10;
	float deltaThetaEHigh = 10;

	float deltaThetaPLow = -10;
	float deltaThetaPHigh = 14;
	
	float deltaPELow = -0.2;
	float deltaPEHigh = 0.2;

	float deltaPPLow = -0.4;
	float deltaPPHigh = 0.4;



	float lowWValueCut[6];// = {0.997378, 0.997977, 0.999086, 0.992596, 0.976615, 0.986326};
	float highWValueCut[6];// = {1.0431, 1.06575, 1.04881, 1.04508, 1.03242, 1.03843};
	//delta phi cut
	float lowPhiValueCut[6];// = {179.614, 178.718, 176.232 , 176.765, 177.109, 174.926 };
	float highPhiValueCut[6];// = { 187.164, 187.166, 184.541, 184.831, 183.889, 182.214};
	
	
	
	if (model == 11){
	  float lowWValueCut[6] = {0.931152, 0.930877, 0.930978, 0.930868, 0.931364, 0.93098};
	  float highWValueCut[6] = {0.965404, 0.967999,0.967974, 0.966602, 0.965179, 0.965636};
	  
	  float lowPhiValueCut[6] = {-1799.614, -1778.718, -1776.232 , -1776.765, -1777.109, -1774.926 };
	  float highPhiValueCut[6] = { 1987.164, 1877.166, 1874.541, 1874.831, 1873.889, 1872.214};
	}
	//if( model == 66666 ){
	std::cout << " Using outbending w fit values. These values are taken from the histograms" << std::endl;
 
	int sector;
	double mean;
	double sigma;

	string line;
	ifstream readFromWCut("w_cut_limits_job"+std::to_string(model)+".txt");
	if( readFromWCut.is_open() ){
	  while(readFromWCut >> sector ) {//std::getline (readFromWCut, line) ){

	    readFromWCut >> mean >> sigma;
	    std::cout << " >> W CUT PARAMETERS: " << sector << " " << mean << " " << sigma << std::endl;

	    lowWValueCut[sector]= mean - 2*sigma;
	    highWValueCut[sector]= mean + 2*sigma;
	  }
	  //}


	  int sector2;
	  double mean2;
	  double sigma2;
    
	  string line2;
	  ifstream readFromPhiCut("phi_cut_limits_job"+std::to_string(model)+".txt");
	  if( readFromPhiCut.is_open() ){
	    while( readFromPhiCut >> sector2  ){

	      readFromPhiCut >> mean2 >> sigma2;
	      std::cout << " >> PHI CUT PARAMETERS: " << sector2 << " " << mean2 << " " << sigma2 << std::endl;
	      
	      lowPhiValueCut[sector2]= mean2 - 2*sigma2;
	      highPhiValueCut[sector2]= mean2 + 2*sigma2;
	    }
	  }

	}
	else{
	  float lowWValueCut[6] = {0.931152, 0.930877, 0.930978, 0.930868, 0.931364, 0.93098};
	  float highWValueCut[6] = {0.965404, 0.967999,0.967974, 0.966602, 0.965179, 0.965636};
	  
	  float lowPhiValueCut[6] = {-1799.614, -1778.718, -1776.232 , -1776.765, -1777.109, -1774.926 };
	  float highPhiValueCut[6] = { 1987.164, 1877.166, 1874.541, 1874.831, 1873.889, 1872.214};
	}

	
	TFile *resultsF = new TFile(Form("monjob%d.root", model), "recreate");
	resultsF->mkdir("overview");
	resultsF->mkdir("eBeam");
	resultsF->mkdir("electronInclusive");
	resultsF->mkdir("energyConservation");
	resultsF->mkdir("thetaElectron");
	resultsF->mkdir("thetaProton");
	resultsF->mkdir("pElectron");
	resultsF->mkdir("pProton");



	TFile f(Form("job%d.root", model));

	float eBeam = 2.2212;
	float mProton = 0.98327;
	TTree *anaTree=(TTree *) f.Get("out_tree");
	TTree *simTree=(TTree *) f.Get("mc_tree");


	float part_px[100], part_py[100], part_pz[100], part_E[100];
	float part_vx[100], part_vy[100], part_vz[100];

	float prot_px[100], prot_py[100], prot_pz[100], prot_E[100];
	float prot_vx[100], prot_vy[100], prot_vz[100];

	float ele_vx[100], ele_vy[100], ele_vz[100];
	int sectorE[100];

	vector<float>  *vpart_px      = 0;
	vector<float>  *vpart_py      = 0;
	vector<float>  *vpart_pz      = 0;

	vector<float>  *vpart_vx      = 0;
	vector<float>  *vpart_vy      = 0;
	vector<float>  *vpart_vz      = 0;
	vector<float>  *vpart_E       = 0;
	vector<int>    *v_sectorE     = 0;

	vector<float>  *vprot_px      = 0;
	vector<float>  *vprot_py      = 0;
	vector<float>  *vprot_pz      = 0;
	vector<float>  *vprot_vx      = 0;
	vector<float>  *vprot_vy      = 0;
	vector<float>  *vprot_vz      = 0;
	vector<float>  *vprot_E       = 0;



	TLorentzVector p4_ele[100], p4_ele_gen[100], initialE, initialP, p4_ele_final;
	TLorentzVector p4_proton[100];
	anaTree->SetBranchAddress("p4_ele_px", &vpart_px);
	anaTree->SetBranchAddress("p4_ele_py", &vpart_py);
	anaTree->SetBranchAddress("p4_ele_pz", &vpart_pz);
	anaTree->SetBranchAddress("p4_ele_E", &vpart_E);

	anaTree->SetBranchAddress("p4_ele_vx", &vpart_vx);
	anaTree->SetBranchAddress("p4_ele_vy", &vpart_vy);
	anaTree->SetBranchAddress("p4_ele_vz", &vpart_vz);

	anaTree->SetBranchAddress("p4_prot_px", &vprot_px);
	anaTree->SetBranchAddress("p4_prot_py", &vprot_py);
	anaTree->SetBranchAddress("p4_prot_pz", &vprot_pz);
	anaTree->SetBranchAddress("p4_prot_E",  &vprot_E);

	anaTree->SetBranchAddress("p4_prot_vx", &vprot_vx);
	anaTree->SetBranchAddress("p4_prot_vy", &vprot_vy);
	anaTree->SetBranchAddress("p4_prot_vz", &vprot_vz);

	anaTree->SetBranchAddress("sectorE", &v_sectorE);
	initialE.SetPxPyPzE(0, 0, eBeam, eBeam);
	initialP.SetPxPyPzE(0, 0, 0, mProton);

//defining histograms
	//overview
		//inclusive
	TH2D *wQ2Inclusive;
	TH2D *wQ2SectorInclusive[6];
	TH1D *wSectorInclusive[6];
	TH2D *thetaPElastic[6];
	//theta e inclusive
	TH2D *deltaThetaEFromEPAsThetaEInclusive[6];
	TH2D *deltaThetaEFromEPAsPEInclusive[6];
	TH1D *deltaThetaEFromEPInclusive[6];

	//p e inclusive
	TH2D *deltaPEFromEThetaAsThetaEInclusive[6];
	TH2D *deltaPEFromEThetaAsPEInclusive[6];
	TH1D *deltaPEFromEThetaInclusive[6];


		//exclusive
	TH2D *wQ2Exclusive;
	TH2D *wQ2SectorExclusive[6];
	TH1D *wSectorExclusive[6];
	TH1D *phiEphiP[6];
	TH1D *phiEPhiWCut[6];
	TH1D *wSectorExclusivePhiCut[6];
	TH2F *thetaEPWSThetaBins[6][6];
	TH2F *zVsPhiElectron;
	TH2F *zVsPhiProton;
	TH2F *deltaZVsPhiElectron;
	TH2F *deltaZVsPhiProton;
	TH1F *zElectronSector[6];
	TH1F *zProtonSector[6];
	TH1F *deltaZElectronProtonSector[6];
	TH2F *zVsThetaElectron[6];
	TH2F *zVsThetaProton[6];
	TH2F *wVsPhi;


	//beam energy histograms
	TH2D *eBeamFromAnglesAsThetaP[6];
	TH2D *eBeamFromAnglesAsThetaE[6];
	TH2D *eBeamFromAnglesAsPP[6];
	TH2D *eBeamFromAnglesAsPE[6];
	TH1D *eBeamFromAngles[6];

	//energy conservation

	TH2D *deltaPxPx[6];
	TH2D *deltaPyPy[6];
	TH2D *deltaPzPz[6];


	//theta E histograms
	TH2D *deltaThetaEFromPPAsThetaE[6];
	TH2D *deltaThetaEFromPPAsThetaP[6];		
	TH2D *deltaThetaEFromPPAsPE[6];
	TH2D *deltaThetaEFromPPAsPP[6];		
	TH1D *deltaThetaEFromPP[6];
	
	TH2D *deltaThetaEFromEPAsThetaE[6];
	TH2D *deltaThetaEFromEPAsThetaP[6];		
	TH2D *deltaThetaEFromEPAsPE[6];
	TH2D *deltaThetaEFromEPAsPP[6];		
	TH1D *deltaThetaEFromEP[6];


	TH2D *deltaThetaEFromEPAsThetaEEnergyConservation[6];
	TH2D *deltaThetaEFromEPAsThetaPEnergyConservation[6];		
	TH2D *deltaThetaEFromEPAsPEEnergyConservation[6];
	TH2D *deltaThetaEFromEPAsPPEnergyConservation[6];		
	TH1D *deltaThetaEFromEPEnergyConservation[6];


	TH2D *deltaThetaEFromPThetaAsThetaE[6];
	TH2D *deltaThetaEFromPThetaAsThetaP[6];		
	TH2D *deltaThetaEFromPThetaAsPE[6];
	TH2D *deltaThetaEFromPThetaAsPP[6];		
	TH1D *deltaThetaEFromPTheta[6];
	

	//theta p histograms

	TH2D *deltaThetaPFromPPAsThetaE[6];
	TH2D *deltaThetaPFromPPAsThetaP[6];		
	TH2D *deltaThetaPFromPPAsPE[6];
	TH2D *deltaThetaPFromPPAsPP[6];		
	TH1D *deltaThetaPFromPP[6];
	
	TH2D *deltaThetaPFromEPAsThetaE[6];
	TH2D *deltaThetaPFromEPAsThetaP[6];		
	TH2D *deltaThetaPFromEPAsPE[6];
	TH2D *deltaThetaPFromEPAsPP[6];		
	TH1D *deltaThetaPFromEP[6];
	
	TH2D *deltaThetaPFromEThetaAsThetaE[6];
	TH2D *deltaThetaPFromEThetaAsThetaP[6];		
	TH2D *deltaThetaPFromEThetaAsPE[6];
	TH2D *deltaThetaPFromEThetaAsPP[6];		
	TH1D *deltaThetaPFromETheta[6];


	//p E histograms
	TH2D *deltaPEFromPPAsThetaE[6];
	TH2D *deltaPEFromPPAsThetaP[6];		
	TH2D *deltaPEFromPPAsPE[6];
	TH2D *deltaPEFromPPAsPP[6];		
	TH1D *deltaPEFromPP[6];

	TH2D *deltaPEFromPThetaAsThetaE[6];
	TH2D *deltaPEFromPThetaAsThetaP[6];		
	TH2D *deltaPEFromPThetaAsPE[6];
	TH2D *deltaPEFromPThetaAsPP[6];		
	TH1D *deltaPEFromPTheta[6];
	
	TH2D *deltaPEFromEThetaAsThetaE[6];
	TH2D *deltaPEFromEThetaAsThetaP[6];		
	TH2D *deltaPEFromEThetaAsPE[6];
	TH2D *deltaPEFromEThetaAsPP[6];		
	TH1D *deltaPEFromETheta[6];

	TH2F *deltaPEFromPThetaAsPhiE;


	//p p histograms
	
	TH2D *deltaPPFromEPAsThetaE[6];
	TH2D *deltaPPFromEPAsThetaP[6];		
	TH2D *deltaPPFromEPAsPE[6];
	TH2D *deltaPPFromEPAsPP[6];		
	TH1D *deltaPPFromEP[6];

	TH2D *deltaPPFromPThetaAsThetaE[6];
	TH2D *deltaPPFromPThetaAsThetaP[6];		
	TH2D *deltaPPFromPThetaAsPE[6];
	TH2D *deltaPPFromPThetaAsPP[6];		
	TH1D *deltaPPFromPTheta[6];
	
	TH2D *deltaPPFromEThetaAsThetaE[6];
	TH2D *deltaPPFromEThetaAsThetaP[6];		
	TH2D *deltaPPFromEThetaAsPE[6];
	TH2D *deltaPPFromEThetaAsPP[6];		
	TH1D *deltaPPFromETheta[6];


	wQ2Inclusive = new TH2D("wQ2Inclusive", "wQ2Inclusive", nWBins, lowBorderW, highBorderW, nQ2Bins, lowBorderQ2, highBorderQ2);
	wQ2Exclusive = new TH2D("wQ2Exclusive", "wQ2Exclusive", nWBins, lowBorderW, highBorderW, nQ2Bins, lowBorderQ2, highBorderQ2);


	zVsPhiElectron = new TH2F("zVsPhiElectron", "zVsPhiElectron", 720, -30, 350, 500, -20, 15);
	zVsPhiElectron->GetYaxis()->SetTitle("Z_{e}, cm");
	zVsPhiElectron->GetXaxis()->SetTitle("#phi_{e}");
	zVsPhiElectron->SetTitle("Z_{e} vs #phi_{e}");
	zVsPhiProton = new TH2F("zVsPhiProton", "zVsPhiProton", 720, -30, 350, 400, -20, 15);
	zVsPhiProton->GetYaxis()->SetTitle("Z_{p}, cm");
	zVsPhiProton->GetXaxis()->SetTitle("#phi_{p}");
	zVsPhiProton->SetTitle("Z_{p} vs #phi_{p}");


	deltaZVsPhiElectron = new TH2F("deltaZVsPhiElectron", "deltaZVsPhiElectron", 720, -30, 350, 400, -20, 15);
	deltaZVsPhiElectron->GetYaxis()->SetTitle("Z_{e} - Z_{p}, cm");
	deltaZVsPhiElectron->GetXaxis()->SetTitle("#phi_{e}");
	deltaZVsPhiElectron->SetTitle("#DeltaZ vs #phi_{e}");
	deltaZVsPhiProton = new TH2F("deltaZVsPhiProton", "deltaZVsPhiProton", 720, -30, 350, 400, -20, 15);
	deltaZVsPhiProton->GetYaxis()->SetTitle("Z_{e} - Z_{p}, cm");
	deltaZVsPhiProton->GetXaxis()->SetTitle("#phi_{p}");
	deltaZVsPhiProton->SetTitle("#DeltaZ vs #phi_{p}");

	deltaPEFromPThetaAsPhiE = new TH2F("deltaPEFromPThetaAsPhiE", "deltaPEFromPThetaAsPhiE", 720, -30, 350, nOverlallBins,deltaPELow, deltaPEHigh);
	deltaPEFromPThetaAsPhiE->SetTitle("#DeltaP_{e} vs #phi_{e} From #theta_{e}");
	deltaPEFromPThetaAsPhiE->GetXaxis()->SetTitle("#phi_{e}");
	deltaPEFromPThetaAsPhiE->GetXaxis()->CenterTitle();
	deltaPEFromPThetaAsPhiE->GetYaxis()->SetTitle("#DeltaP_{e}");
	deltaPEFromPThetaAsPhiE->GetYaxis()->CenterTitle();
	
	
	wVsPhi = new TH2F("wVsPhi", "wVsPhi", 720, -30, 350, nWBins, lowBorderW, highBorderW - 0.55);//+ 0.1, highBorderW - 0.6);
	wVsPhi->GetXaxis()->SetTitle("#phi_{e}");
	wVsPhi->GetXaxis()->CenterTitle();	
	wVsPhi->GetYaxis()->SetTitle("W, GeV");
	wVsPhi->GetYaxis()->CenterTitle();

	wVsPhi->SetTitle("W vs #phi_{e}");
	for (int s = 0; s < 6; s++){
		//overview
		//inclusive
		wQ2SectorInclusive[s] = new TH2D(Form("wQ2InclusiveS%d", s + 1), Form("wQ2InclusiveS%d", s + 1), nWBins, lowBorderW, highBorderW, nQ2Bins, lowBorderQ2, highBorderQ2);
		wQ2SectorInclusive[s]->GetXaxis()->SetTitle("W, GeV");
		wQ2SectorInclusive[s]->GetYaxis()->SetTitle("Q^{2}, GeV^{2}");

		wSectorInclusive[s] = new TH1D(Form("wInclusiveS%d", s + 1), Form("wInclusiveS%d", s + 1), nWBins, lowBorderW, highBorderW);
		wSectorInclusive[s]->GetXaxis()->SetTitle("W, GeV");
		thetaPElastic[s]= new TH2D(Form("thetaPElasticS%d", s + 1), Form("thetaPElasticS%d", s + 1), 200, 0, 2.5, 200, 0, 35);
		//exclusive
		wQ2SectorExclusive[s] = new TH2D(Form("wQ2ExclusiveS%d", s + 1), Form("wQ2ExclusiveS%d", s + 1), nWBins, lowBorderW, highBorderW, nQ2Bins, lowBorderQ2, highBorderQ2);
		wQ2SectorExclusive[s]->GetXaxis()->SetTitle("W, GeV");
		wQ2SectorExclusive[s]->GetYaxis()->SetTitle("Q^{2}, GeV^{2}");
		wQ2SectorExclusive[s]->SetTitle(Form("Sector %d, Q^{2} vs W, e(FD) + positive(CD)", s + 1));
		wSectorExclusive[s] = new TH1D(Form("wExclusiveS%d", s + 1), Form("wExclusiveS%d", s + 1), nWBins, lowBorderW, highBorderW);
		wSectorExclusive[s]->GetXaxis()->SetTitle("W, GeV");
		wSectorExclusive[s]->SetTitle(Form("Sector %d, W, e(FD) + positive(CD)", s + 1));
		phiEphiP[s] = new TH1D(Form("phiEphiPS%d", s + 1), Form("phiEphiPS%d", s + 1), 360, 150, 210);
		phiEphiP[s]->GetXaxis()->SetTitle("#phi_e - #phi_p");
		phiEPhiWCut[s] = new TH1D(Form("phiEPhiWCutS%d", s + 1), Form("phiePhiWCutS%d", s + 1), 360, 150, 210);
		phiEPhiWCut[s]->GetXaxis()->SetTitle("#phi_e - #phi_p with W cut");
		wSectorExclusivePhiCut[s] = new TH1D(Form("wSectorExclusivePhiCutS%d", s + 1), Form("wSectorExclusivePhiCutS%d", s + 1), nWBins, lowBorderW, highBorderW);
		wSectorExclusivePhiCut[s]->GetXaxis()->SetTitle("W, GeV");
		wSectorExclusivePhiCut[s]->SetTitle(Form("Sector %d, Exclusive W with all cuts", s + 1));
		zElectronSector[s] = new TH1F(Form("zElectronSectorS%d", s + 1), Form("zElectronSectorS%d", s + 1),  400, -20, 15);
		zElectronSector[s]->SetTitle(Form("Sector %d, Z_{e}", s + 1));
		zElectronSector[s]->GetXaxis()->SetTitle("Z_{e}, cm");
		zProtonSector[s] = new TH1F(Form("zProtonSectorS%d", s + 1), Form("zProtonSectorS%d", s + 1),  400, -20, 15);
		zProtonSector[s]->SetTitle(Form("Sector %d, Z_{p}", s + 1));
		zProtonSector[s]->GetXaxis()->SetTitle("Z_{p}, cm");
		deltaZElectronProtonSector[s] = new TH1F(Form("deltaZElectronProtonSectorS%d", s + 1), Form("deltaZElectronProtonSectorS%d", s + 1),  400, -20, 15);
		deltaZElectronProtonSector[s]->SetTitle(Form("Sector %d, Z_{e} - Z_{p}", s + 1));
		deltaZElectronProtonSector[s]->GetXaxis()->SetTitle("Z_{e} - Z_{p}");

		zVsThetaElectron[s]= new TH2F(Form("zVsThetaElectronS%d", s + 1), Form("zVsThetaElectronS%d", s + 1), nOverlallBins, thetaELow, thetaEHigh,  400, -20, 15);
		zVsThetaElectron[s]->GetYaxis()->SetTitle("Z_{e}, cm");
		zVsThetaElectron[s]->GetXaxis()->SetTitle("#theta_{e}");
		zVsThetaElectron[s]->SetTitle(Form("Sector %d, Z_{e} vs #theta_{e}", s + 1));
		zVsThetaProton[s]= new TH2F(Form("zVsThetaProtonS%d", s + 1), Form("zVsThetaProtonS%d", s + 1),  nOverlallBins, thetaPLow, thetaPHigh, 400, -20, 15);
		zVsThetaProton[s]->GetYaxis()->SetTitle("Z_{p}, cm");
		zVsThetaProton[s]->GetXaxis()->SetTitle("#theta_{p}");
		zVsThetaProton[s]->SetTitle(Form("Sector %d, Z_{p} vs #theta_{p}", s + 1));

		for (int t = 0; t < 6; t++){
			thetaEPWSThetaBins[s][t] = new TH2F(Form("thetaEPWSThetaBins%dT%d", s +  1, t + 1), Form("thetaEPWSThetaBins%dT%d", s +  1, t + 1), 390, 0.6, 1.9, 300 , 50, 100);
			thetaEPWSThetaBins[s][t]->GetXaxis()->SetTitle("W, GeV");
			thetaEPWSThetaBins[s][t]->GetYaxis()->SetTitle("#theta_{e} + #theta_{p}");
			thetaEPWSThetaBins[s][t]->SetTitle(Form("%d < #theta_{e} < %d", 5 + 5*t, 10 + 5*t));
		}
		//beam energy
		eBeamFromAnglesAsThetaP[s] = new TH2D(Form("eBeamFromAnglesAsThetaPS%d", s + 1), Form("eBeamFromAnglesAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaEenergyLow, deltaEenergyHigh);
		eBeamFromAnglesAsThetaP[s]->GetYaxis()->SetTitle("#DeltaE_{beam}");
		eBeamFromAnglesAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		eBeamFromAnglesAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaE_{beam} From #theta_{e} and #theta_{p}", s + 1));
		eBeamFromAnglesAsThetaE[s] = new TH2D(Form("eBeamFromAnglesAsThetaES%d", s + 1), Form("eBeamFromAnglesAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaEenergyLow, deltaEenergyHigh);
		eBeamFromAnglesAsThetaE[s]->GetYaxis()->SetTitle("#DeltaE_{beam}");
		eBeamFromAnglesAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		eBeamFromAnglesAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaE_{beam} From #theta_{e} and #theta_{p}", s + 1));
		eBeamFromAnglesAsPP[s] = new TH2D(Form("eBeamFromAnglesAsPPS%d", s + 1), Form("eBeamFromAnglesAsPPS%d", s + 1), nOverlallBins, pPLow, pPHigh, nOverlallBins, deltaEenergyLow, deltaEenergyHigh);
		eBeamFromAnglesAsPP[s]->GetYaxis()->SetTitle("#DeltaE_{beam}");
		eBeamFromAnglesAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		eBeamFromAnglesAsPP[s]->SetTitle(Form("Sector %d, #DeltaE_{beam} From #theta_{e} and #theta_{p}", s + 1));
		eBeamFromAnglesAsPE[s] = new TH2D(Form("eBeamFromAnglesAsPES%d", s + 1), Form("eBeamFromAnglesAsPES%d", s + 1),  nOverlallBins, pELow, pEHigh, nOverlallBins, deltaEenergyLow, deltaEenergyHigh);
		eBeamFromAnglesAsPE[s]->GetYaxis()->SetTitle("#DeltaE_{beam}");
		eBeamFromAnglesAsPE[s]->GetXaxis()->SetTitle("E_p");
		eBeamFromAnglesAsPE[s]->SetTitle(Form("Sector %d, #DeltaE_{beam} From #theta_{e} and #theta_{p}", s + 1));
		eBeamFromAngles[s] = new TH1D(Form("eBeamFromAnglesS%d", s + 1), Form("eBeamFromAnglesS%d", s + 1), nOverlallBins, deltaEenergyLow, deltaEenergyHigh);
		eBeamFromAngles[s]->GetXaxis()->SetTitle("#DeltaE_{beam}");
		eBeamFromAngles[s]->SetTitle(Form("Sector %d, #DeltaE_{beam} From #theta_{e} and #theta_{p}", s + 1));

		//		energy conservation	
		deltaPxPx[s] = new TH2D(Form("deltaPxPxS%d", s + 1), Form("deltaPxPxS%d", s + 1), 100, -0.7, 0.7, 100, -0.2, 0.2);
		deltaPxPx[s]->GetXaxis()->SetTitle("P_x");
		deltaPxPx[s]->GetYaxis()->SetTitle("#DeltaP_x");
		deltaPyPy[s] = new TH2D(Form("deltaPyPyS%d", s + 1), Form("deltaPyPyS%d", s + 1), 100, -0.7, 0.7, 100, -0.2, 0.2);
		deltaPyPy[s]->GetXaxis()->SetTitle("P_y");
		deltaPyPy[s]->GetYaxis()->SetTitle("#DeltaP_y");
		deltaPzPz[s] = new TH2D(Form("deltaPzPzS%d", s + 1), Form("deltaPzPzS%d", s + 1), 100, 1.5, 2.5, 100, -0.3, 0.3);
		deltaPzPz[s]->GetXaxis()->SetTitle("P_z");
		deltaPzPz[s]->GetYaxis()->SetTitle("#DeltaP_z");

		//theta electron inclusive
		deltaThetaEFromEPAsThetaEInclusive[s] = new TH2D(Form("deltaThetaEFromEPAsThetaEInclusiveS%d", s + 1), Form("deltaThetaEFromEPAsThetaEInclusiveS%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsThetaEInclusive[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsThetaEInclusive[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaEFromEPAsPEInclusive[s] = new TH2D(Form("deltaThetaEFromEPAsPEInclusiveS%d", s + 1), Form("deltaThetaEFromEPAsPEInclusiveS%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsPEInclusive[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsPEInclusive[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaEFromEPInclusive[s] = new TH1D(Form("deltaThetaEFromEPInclusiveS%d", s + 1), Form("deltaThetaEFromEPInclusiveS%d", s + 1), nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPInclusive[s]->GetXaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPInclusive[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From #theta_{e}", s + 1));

	//p electron inclusive
		deltaPEFromEThetaAsThetaEInclusive[s] = new TH2D(Form("deltaPEFromEThetaAsThetaEInclusiveS%d", s + 1), Form("deltaPEFromEThetaAsThetaEInclusiveS%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsThetaEInclusive[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsThetaEInclusive[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPEFromEThetaAsPEInclusive[s] = new TH2D(Form("deltaPEFromEThetaAsPEInclusiveS%d", s + 1), Form("deltaPEFromEThetaAsPEInclusiveS%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins,deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsPEInclusive[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsPEInclusive[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPEFromEThetaInclusive[s] = new TH1D(Form("deltaPEFromEThetaInclusiveS%d", s + 1), Form("deltaPEFromEThetaInclusiveS%d", s + 1), nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaInclusive[s]->GetXaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaInclusive[s]->SetTitle(Form("Sector %d, #DeltaP_{e} From #theta_{e}", s + 1));



		//theta electron
		deltaThetaEFromPPAsThetaE[s] = new TH2D(Form("deltaThetaEFromPPAsThetaES%d", s + 1), Form("deltaThetaEFromPPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPPAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaEFromPPAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{e} From P_{p}", s + 1));
		deltaThetaEFromPPAsThetaP[s] = new TH2D(Form("deltaThetaEFromPPAsThetaPS%d", s + 1), Form("deltaThetaEFromPPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPPAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaEFromPPAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{p} From P_{p}", s + 1));
		deltaThetaEFromPPAsPE[s] = new TH2D(Form("deltaThetaEFromPPAsPES%d", s + 1), Form("deltaThetaEFromPPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPPAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaEFromPPAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{e} From P_{p}", s + 1));
		deltaThetaEFromPPAsPP[s] = new TH2D(Form("deltaThetaEFromPPAsPPS%d", s + 1), Form("deltaThetaEFromPPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPPAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaEFromPPAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{p} From P_{p}", s + 1));
		deltaThetaEFromPP[s] = new TH1D(Form("deltaThetaEFromPPS%d", s + 1), Form("deltaThetaEFromPPS%d", s + 1), nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPP[s]->GetXaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{p}", s + 1));

		deltaThetaEFromEPAsThetaE[s] = new TH2D(Form("deltaThetaEFromEPAsThetaES%d", s + 1), Form("deltaThetaEFromEPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaEFromEPAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{e} From P_{e}", s + 1));
		deltaThetaEFromEPAsThetaP[s] = new TH2D(Form("deltaThetaEFromEPAsThetaPS%d", s + 1), Form("deltaThetaEFromEPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaEFromEPAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{p} From P_{e}", s + 1));
		deltaThetaEFromEPAsPE[s] = new TH2D(Form("deltaThetaEFromEPAsPES%d", s + 1), Form("deltaThetaEFromEPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaEFromEPAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{e} From P_{e}", s + 1));
		deltaThetaEFromEPAsPP[s] = new TH2D(Form("deltaThetaEFromEPAsPPS%d", s + 1), Form("deltaThetaEFromEPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaEFromEPAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{p} From P_{e}", s + 1));
		deltaThetaEFromEP[s] = new TH1D(Form("deltaThetaEFromEPS%d", s + 1), Form("deltaThetaEFromEPS%d", s + 1), nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEP[s]->GetXaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));



		deltaThetaEFromEPAsThetaEEnergyConservation[s] = new TH2D(Form("deltaThetaEFromEPAsThetaEEnergyConservationS%d", s + 1), Form("deltaThetaEFromEPAsThetaEEnergyConservationS%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsThetaEEnergyConservation[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsThetaEEnergyConservation[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaEFromEPAsThetaEEnergyConservation[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));
		deltaThetaEFromEPAsThetaPEnergyConservation[s] = new TH2D(Form("deltaThetaEFromEPAsThetaPEnergyConservationS%d", s + 1), Form("deltaThetaEFromEPAsThetaPEnergyConservationS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsThetaPEnergyConservation[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPAsThetaPEnergyConservation[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaEFromEPAsThetaPEnergyConservation[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));
		deltaThetaEFromEPAsPEEnergyConservation[s] = new TH2D(Form("deltaThetaEFromEPAsPEEnergyConservationS%d", s + 1), Form("deltaThetaEFromEPAsPEEnergyConservationS%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsPEEnergyConservation[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaThetaEFromEPAsPEEnergyConservation[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaEFromEPAsPEEnergyConservation[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));
		deltaThetaEFromEPAsPPEnergyConservation[s] = new TH2D(Form("deltaThetaEFromEPAsPPEnergyConservationS%d", s + 1), Form("deltaThetaEFromEPAsPPEnergyConservationS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPAsPPEnergyConservation[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaThetaEFromEPAsPPEnergyConservation[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaEFromEPAsPPEnergyConservation[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));
		deltaThetaEFromEPEnergyConservation[s] = new TH1D(Form("deltaThetaEFromEPEnergyConservationS%d", s + 1), Form("deltaThetaEFromEPEnergyConservationS%d", s + 1), nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromEPEnergyConservation[s]->GetXaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromEPEnergyConservation[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From P_{e}", s + 1));



		deltaThetaEFromPThetaAsThetaE[s] = new TH2D(Form("deltaThetaEFromPThetaAsThetaES%d", s + 1), Form("deltaThetaEFromPThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPThetaAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaEFromPThetaAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{e} From #theta_{p}", s + 1));
		deltaThetaEFromPThetaAsThetaP[s] = new TH2D(Form("deltaThetaEFromPThetaAsThetaPS%d", s + 1), Form("deltaThetaEFromPThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPThetaAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaEFromPThetaAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs #theta_{p} From #theta_{p}", s + 1));
		deltaThetaEFromPThetaAsPE[s] = new TH2D(Form("deltaThetaEFromPThetaAsPES%d", s + 1), Form("deltaThetaEFromPThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPThetaAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaEFromPThetaAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{e} From #theta_{p}", s + 1));
		deltaThetaEFromPThetaAsPP[s] = new TH2D(Form("deltaThetaEFromPThetaAsPPS%d", s + 1), Form("deltaThetaEFromPThetaAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPThetaAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaEFromPThetaAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} vs P_{p} From #theta_{p}", s + 1));
		deltaThetaEFromPTheta[s] = new TH1D(Form("deltaThetaEFromPThetaS%d", s + 1), Form("deltaThetaEFromPThetaS%d", s + 1), nOverlallBins, deltaThetaELow, deltaThetaEHigh);
		deltaThetaEFromPTheta[s]->GetXaxis()->SetTitle("#Delta#theta_{e}");
		deltaThetaEFromPTheta[s]->SetTitle(Form("Sector %d, #Delta#theta_{e} From #theta_{p}", s + 1));

		//theta proton
		deltaThetaPFromPPAsThetaE[s] = new TH2D(Form("deltaThetaPFromPPAsThetaES%d", s + 1), Form("deltaThetaPFromPPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromPPAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromPPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaPFromPPAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{e} From P_{p}", s + 1));
		deltaThetaPFromPPAsThetaP[s] = new TH2D(Form("deltaThetaPFromPPAsThetaPS%d", s + 1), Form("deltaThetaPFromPPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromPPAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromPPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaPFromPPAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{p} From P_{p}", s + 1));
		deltaThetaPFromPPAsPE[s] = new TH2D(Form("deltaThetaPFromPPAsPES%d", s + 1), Form("deltaThetaPFromPPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromPPAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromPPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaPFromPPAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{e} From P_{p}", s + 1));
		deltaThetaPFromPPAsPP[s] = new TH2D(Form("deltaThetaPFromPPAsPPS%d", s + 1), Form("deltaThetaPFromPPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromPPAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromPPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaPFromPPAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{p} From P_{p}", s + 1));
		deltaThetaPFromPP[s] = new TH1D(Form("deltaThetaPFromPPS%d", s + 1), Form("deltaThetaPFromPPS%d", s + 1), nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromPP[s]->GetXaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} From P_{p}", s + 1));

		deltaThetaPFromEPAsThetaE[s] = new TH2D(Form("deltaThetaPFromEPAsThetaES%d", s + 1), Form("deltaThetaPFromEPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEPAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaPFromEPAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{e} From P_{e}", s + 1));
		deltaThetaPFromEPAsThetaP[s] = new TH2D(Form("deltaThetaPFromEPAsThetaPS%d", s + 1), Form("deltaThetaPFromEPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEPAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaPFromEPAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{p} From P_{e}", s + 1));
		deltaThetaPFromEPAsPE[s] = new TH2D(Form("deltaThetaPFromEPAsPES%d", s + 1), Form("deltaThetaPFromEPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEPAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaPFromEPAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{e} From P_{e}", s + 1));
		deltaThetaPFromEPAsPP[s] = new TH2D(Form("deltaThetaPFromEPAsPPS%d", s + 1), Form("deltaThetaPFromEPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEPAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaPFromEPAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{p} From P_{e}", s + 1));
		deltaThetaPFromEP[s] = new TH1D(Form("deltaThetaPFromEPS%d", s + 1), Form("deltaThetaPFromPPS%d", s + 1), nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEP[s]->GetXaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} From P_{e}", s + 1));

		deltaThetaPFromEThetaAsThetaE[s] = new TH2D(Form("deltaThetaPFromEThetaAsThetaES%d", s + 1), Form("deltaThetaPFromEThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEThetaAsThetaE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaThetaPFromEThetaAsThetaE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{e} From #theta_{e}", s + 1));
		deltaThetaPFromEThetaAsThetaP[s] = new TH2D(Form("deltaThetaPFromEThetaAsThetaPS%d", s + 1), Form("deltaThetaPFromEThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEThetaAsThetaP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaThetaPFromEThetaAsThetaP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs #theta_{p} From #theta_{e}", s + 1));
		deltaThetaPFromEThetaAsPE[s] = new TH2D(Form("deltaThetaPFromEThetaAsPES%d", s + 1), Form("deltaThetaPFromEThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEThetaAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaThetaPFromEThetaAsPE[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{e} From #theta_{e}", s + 1));
		deltaThetaPFromEThetaAsPP[s] = new TH2D(Form("deltaThetaPFromEThetaAsPPS%d", s + 1), Form("deltaThetaPFromEThetaAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromEThetaAsPP[s]->GetYaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromEThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaThetaPFromEThetaAsPP[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} vs P_{p} From #theta_{e}", s + 1));
		deltaThetaPFromETheta[s] = new TH1D(Form("deltaThetaPFromEThetaS%d", s + 1), Form("deltaThetaPFromEThetaS%d", s + 1), nOverlallBins, deltaThetaPLow, deltaThetaPHigh);
		deltaThetaPFromETheta[s]->GetXaxis()->SetTitle("#Delta#theta_{p}");
		deltaThetaPFromETheta[s]->SetTitle(Form("Sector %d, #Delta#theta_{p} From #theta_{e}", s + 1));

		//p electron
		deltaPEFromPPAsThetaE[s] = new TH2D(Form("deltaPEFromPPAsThetaES%d", s + 1), Form("deltaPEFromPPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPPAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPEFromPPAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{e} From P_{p}", s + 1));
		deltaPEFromPPAsThetaP[s] = new TH2D(Form("deltaPEFromPPAsThetaPS%d", s + 1), Form("deltaPEFromPPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPPAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPEFromPPAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{p} From P_{p}", s + 1));
		deltaPEFromPPAsPE[s] = new TH2D(Form("deltaPEFromPPAsPES%d", s + 1), Form("deltaPEFromPPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPPAsPE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPEFromPPAsPE[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPEFromPPAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{e} From P_{p}", s + 1));
		deltaPEFromPPAsPP[s] = new TH2D(Form("deltaPEFromPPAsPPS%d", s + 1), Form("deltaPEFromPPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPPAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPEFromPPAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPPAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{p} From P_{p}", s + 1));
		deltaPEFromPP[s] = new TH1D(Form("deltaPEFromPPS%d", s + 1), Form("deltaPEFromPPS%d", s + 1), nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPP[s]->GetXaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} From P_{p}", s + 1));

		deltaPEFromEThetaAsThetaE[s] = new TH2D(Form("deltaPEFromEThetaAsThetaES%d", s + 1), Form("deltaPEFromEThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPEFromEThetaAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{e} From #theta_{e}", s + 1));
		deltaPEFromEThetaAsThetaP[s] = new TH2D(Form("deltaPEFromEThetaAsThetaPS%d", s + 1), Form("deltaPEFromEThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPEFromEThetaAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{p} From #theta_{e}", s + 1));
		deltaPEFromEThetaAsPE[s] = new TH2D(Form("deltaPEFromEThetaAsPES%d", s + 1), Form("deltaPEFromEThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsPE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPEFromEThetaAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{e} From #theta_{e}", s + 1));
		deltaPEFromEThetaAsPP[s] = new TH2D(Form("deltaPEFromEThetaAsPPS%d", s + 1), Form("deltaPEFromEThetaAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromEThetaAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromEThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPEFromEThetaAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{p} From #theta_{e}", s + 1));
		deltaPEFromETheta[s] = new TH1D(Form("deltaPEFromEThetaS%d", s + 1), Form("deltaPEFromEThetaS%d", s + 1), nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromETheta[s]->GetXaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromETheta[s]->SetTitle(Form("Sector %d, #DeltaP_{e} From #theta_{e}", s + 1));

		deltaPEFromPThetaAsThetaE[s] = new TH2D(Form("deltaPEFromPThetaAsThetaES%d", s + 1), Form("deltaPEFromPThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPThetaAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPEFromPThetaAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{e} From #theta_{p}", s + 1));
		deltaPEFromPThetaAsThetaP[s] = new TH2D(Form("deltaPEFromPThetaAsThetaPS%d", s + 1), Form("deltaPEFromPThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPThetaAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPEFromPThetaAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs #theta_{p} From #theta_{p}", s + 1));
		deltaPEFromPThetaAsPE[s] = new TH2D(Form("deltaPEFromPThetaAsPES%d", s + 1), Form("deltaPEFromPThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPThetaAsPE[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPEFromPThetaAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{e} From #theta_{p}", s + 1));
		deltaPEFromPThetaAsPP[s] = new TH2D(Form("deltaPEFromPThetaAsPPS%d", s + 1), Form("deltaPEFromPThetaAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPThetaAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPEFromPThetaAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{e} vs P_{p} From #theta_{p}", s + 1));
		deltaPEFromPTheta[s] = new TH1D(Form("deltaPEFromPThetaS%d", s + 1), Form("deltaPEFromPThetaS%d", s + 1), nOverlallBins, deltaPELow, deltaPEHigh);
		deltaPEFromPTheta[s]->GetXaxis()->SetTitle("#DeltaP_{e}");
		deltaPEFromPTheta[s]->SetTitle(Form("Sector %d, #DeltaP_{e} From #theta_{p}", s + 1));

		//p proton
		deltaPPFromPThetaAsThetaE[s] = new TH2D(Form("deltaPPFromPThetaAsThetaES%d", s + 1), Form("deltaPPFromPThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromPThetaAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromPThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPPFromPThetaAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{e} From #theta_{p}", s + 1));
		deltaPPFromPThetaAsThetaP[s] = new TH2D(Form("deltaPPFromPThetaAsThetaPS%d", s + 1), Form("deltaPPFromPThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromPThetaAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromPThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPPFromPThetaAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{p} From #theta_{p}", s + 1));
		deltaPPFromPThetaAsPE[s] = new TH2D(Form("deltaPPFromPThetaAsPES%d", s + 1), Form("deltaPPFromPThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromPThetaAsPE[s]->GetYaxis()->SetTitle("#Delta#theta_{e}");
		deltaPPFromPThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPPFromPThetaAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{e} From #theta_{p}", s + 1));
		deltaPPFromPThetaAsPP[s] = new TH2D(Form("deltaPPFromPThetaAsPPS%d", s + 1), Form("deltaPThetaFromPPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromPThetaAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromPThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPPFromPThetaAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{p} From #theta_{p}", s + 1));
		deltaPPFromPTheta[s] = new TH1D(Form("deltaPPFromPThetaS%d", s + 1), Form("deltaPPFromPThetaS%d", s + 1), nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromPTheta[s]->GetXaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromPTheta[s]->SetTitle(Form("Sector %d, #DeltaP_{p} From #theta_{p}", s + 1));

		deltaPPFromEPAsThetaE[s] = new TH2D(Form("deltaPPFromEPAsThetaES%d", s + 1), Form("deltaPPFromEPAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEPAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEPAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPPFromEPAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{e} From P_{e}", s + 1));
		deltaPPFromEPAsThetaP[s] = new TH2D(Form("deltaPPFromEPAsThetaPS%d", s + 1), Form("deltaPPFromEPAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEPAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEPAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPPFromEPAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{p} From P_{e}", s + 1));
		deltaPPFromEPAsPE[s] = new TH2D(Form("deltaPPFromEPAsPES%d", s + 1), Form("deltaPPFromEPAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEPAsPE[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEPAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPPFromEPAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{e} From P_{e}", s + 1));
		deltaPPFromEPAsPP[s] = new TH2D(Form("deltaPPFromEPAsPPS%d", s + 1), Form("deltaPPFromEPAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEPAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEPAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPPFromEPAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{p} From P_{e}", s + 1));
		deltaPPFromEP[s] = new TH1D(Form("deltaPPFromEPS%d", s + 1), Form("deltaPPFromEPS%d", s + 1), nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEP[s]->GetXaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} From P_{e}", s + 1));

		deltaPPFromEThetaAsThetaE[s] = new TH2D(Form("deltaPPFromEThetaAsThetaES%d", s + 1), Form("deltaPPFromEThetaAsThetaES%d", s + 1), nOverlallBins, thetaELow, thetaEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEThetaAsThetaE[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEThetaAsThetaE[s]->GetXaxis()->SetTitle("#theta_{e}");
		deltaPPFromEThetaAsThetaE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{e} From #theta_{e}", s + 1));
		deltaPPFromEThetaAsThetaP[s] = new TH2D(Form("deltaPPFromEThetaAsThetaPS%d", s + 1), Form("deltaPPFromEThetaAsThetaPS%d", s + 1), nOverlallBins, thetaPLow, thetaPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEThetaAsThetaP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEThetaAsThetaP[s]->GetXaxis()->SetTitle("#theta_{p}");
		deltaPPFromEThetaAsThetaP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs #theta_{p}  From #theta_{e}", s + 1));
		deltaPPFromEThetaAsPE[s] = new TH2D(Form("deltaPPFromEThetaAsPES%d", s + 1), Form("deltaPPFromEThetaAsPES%d", s + 1), nOverlallBins, pELow, pEHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEThetaAsPE[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEThetaAsPE[s]->GetXaxis()->SetTitle("P_{e}");
		deltaPPFromEThetaAsPE[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{e} From #theta_{e}", s + 1));
		deltaPPFromEThetaAsPP[s] = new TH2D(Form("deltaPPFromEThetaAsPPS%d", s + 1), Form("deltaPPFromEThetaAsPPS%d", s + 1), nOverlallBins, pPLow,  pPHigh, nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromEThetaAsPP[s]->GetYaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromEThetaAsPP[s]->GetXaxis()->SetTitle("P_{p}");
		deltaPPFromEThetaAsPP[s]->SetTitle(Form("Sector %d, #DeltaP_{p} vs P_{p} From #theta_{e}", s + 1));
		deltaPPFromETheta[s] = new TH1D(Form("deltaPPFromEThetaS%d", s + 1), Form("deltaPPFromEThetaS%d", s + 1), nOverlallBins, deltaPPLow, deltaPPHigh);
		deltaPPFromETheta[s]->GetXaxis()->SetTitle("#DeltaP_{p}");
		deltaPPFromETheta[s]->SetTitle(Form("Sector %d, #DeltaP_{p} From #theta_{e}", s + 1));

	}

	int NPart = 0;
	int NPartP = 0;
	int nEvents = 1000000;
	float toRD = 57.2958;
	float thetaEFromThetaPValue, thetaEFromPPValue, thetaEFromThetaEValue, thetaEFromEPValue;
	float pEFromThetaPValue, pEFromPPValue, pEFromThetaEValue, pEFromPEValue;
	float thetaPFromThetaPValue, thetaPFromPPValue, thetaPFromThetaEValue, thetaPFromEPValue;
	float pPFromThetaPValue, pPFromPPValue, pPFromThetaEValue, pPFromEPValue;
	float beamEnergyCalc;
	float thetaEFromEPValueEnergyConserve;
	float thetaEMeasured, thetaPMeasured, pEMeasured, pPMeasured;
	float phiEMeasured, phiPMeasured;

	int thetaBin = 0;
	for(Int_t k=0; k<anaTree->GetEntriesFast();k++){
		anaTree->GetEntry(k);
		NPart = vpart_px->size();
		if (NPart > 0 && NPart < 2){
			for(Int_t i = 0; i < NPart; i++){
				part_px[i] = vpart_px->at(i);
				part_py[i] = vpart_py->at(i);
				part_pz[i] = vpart_pz->at(i);
				part_vx[i] = vpart_vx->at(i);
				part_vy[i] = vpart_vy->at(i);
				part_vz[i] = vpart_vz->at(i);
				part_E[i] = vpart_E->at(i);
				sectorE[i] = v_sectorE->at(i);
				sectorElectron = sectorE[i] - 1;
				if (sectorElectron > -1 && sectorElectron < 7){
					p4_ele[i].SetPxPyPzE(part_px[i], part_py[i], part_pz[i], part_E[i]);
					W = kin_W(p4_ele[i],  eBeam);
					Q2 = kin_Q2(p4_ele[i],  eBeam);
					wQ2Inclusive->Fill(W, Q2);
					wQ2SectorInclusive[sectorElectron]->Fill(W, Q2);
					wSectorInclusive[sectorElectron]->Fill(W);
					thetaEMeasured = p4_ele[i].Theta()*toRD;
					pEMeasured = p4_ele[i].P();
					if (W < highWValueCut[sectorElectron] && W >  lowWValueCut[sectorElectron]){
											//theta electron inclusive
						thetaEFromEPValue = returnThetaEFromPE(eBeam, p4_ele[i]);
						thetaPElastic[sectorElectron]->Fill(pEMeasured, thetaEMeasured);
						deltaThetaEFromEPAsThetaEInclusive[sectorElectron]->Fill(thetaEMeasured, thetaEFromEPValue - thetaEMeasured);
						deltaThetaEFromEPAsPEInclusive[sectorElectron]->Fill(pEMeasured ,thetaEFromEPValue - thetaEMeasured);
						deltaThetaEFromEPInclusive[sectorElectron]->Fill(thetaEFromEPValue - thetaEMeasured);
						//p electron inclusive
						pEFromThetaEValue = returnPEFromThetaE(eBeam, p4_ele[i]);
						deltaPEFromEThetaAsThetaEInclusive[sectorElectron]->Fill(thetaEMeasured, pEFromThetaEValue - pEMeasured);
						deltaPEFromEThetaAsPEInclusive[sectorElectron]->Fill(pEMeasured, pEFromThetaEValue -pEMeasured);
						deltaPEFromEThetaInclusive[sectorElectron]->Fill(pEFromThetaEValue - pEMeasured);
					}
					NPartP = vprot_px->size();
					if (NPartP > 0 && NPartP < 2){
						for(Int_t j = 0; j < NPartP; j++){
							prot_px[j] = vprot_px->at(j);
							prot_py[j] = vprot_py->at(j);
							prot_pz[j] = vprot_pz->at(j);
							prot_E[j]  = vprot_E->at(j);
							prot_vx[i] = vprot_vx->at(i);
							prot_vy[i] = vprot_vy->at(i);
							prot_vz[i] = vprot_vz->at(i);
							p4_proton[j].SetPxPyPzE(prot_px[j], prot_py[j], prot_pz[j], prot_E[j]);
							thetaPMeasured = p4_proton[j].Theta()*toRD;
							pPMeasured = p4_proton[j].P();
							wQ2Exclusive->Fill(W, Q2);
							wQ2SectorExclusive[sectorElectron]->Fill(W, Q2);
							wSectorExclusive[sectorElectron]->Fill(W);
							p4_ele_final = initialE + initialP -p4_proton[j];
							float deltaPhi = p4_ele[i].Phi()*toRD - p4_proton[j].Phi()*toRD;
							phiEMeasured = returnPhi(p4_ele[i].Phi()*toRD);
							phiPMeasured = returnPhi(p4_proton[i].Phi()*toRD);
							wVsPhi->Fill(phiEMeasured, W);
							zVsPhiElectron->Fill(phiEMeasured, part_vz[i]);
							zVsPhiProton->Fill(phiPMeasured, prot_vz[i]);
							deltaZVsPhiElectron->Fill(phiEMeasured, part_vz[i] - prot_vz[j]);
							deltaZVsPhiProton->Fill(phiPMeasured, part_vz[i] - prot_vz[j]);
							zVsThetaProton[sectorElectron]->Fill(thetaPMeasured, prot_vz[i]);
							zVsThetaElectron[sectorElectron]->Fill(thetaEMeasured, part_vz[i]);
							if (sectorElectron  > -1 && sectorElectron < 6){
								zElectronSector[sectorElectron]->Fill(part_vz[i]);
								zProtonSector[sectorElectron]->Fill(prot_vz[j]);
								deltaZElectronProtonSector[sectorElectron]->Fill(part_vz[i] - prot_vz[j]);
							}
							if (deltaPhi < 0) deltaPhi = -deltaPhi;
							phiEphiP[sectorElectron]->Fill(deltaPhi);
							thetaBin = (p4_ele[0].Theta()*57.2958 - 5)/5;
							if (sectorElectron  > -1 && sectorElectron < 6 && thetaBin > -1 && thetaBin < 6){
								thetaEPWSThetaBins[sectorElectron][thetaBin]->Fill(W, p4_proton[i].Theta()*57.2958 + p4_ele[j].Theta()*57.2958);
							}
							if (W < highWValueCut[sectorElectron] && W >  lowWValueCut[sectorElectron] && deltaPhi > lowPhiValueCut[sectorElectron] && deltaPhi < highPhiValueCut[sectorElectron]){
								phiEPhiWCut[sectorElectron]->Fill(deltaPhi);
								wSectorExclusivePhiCut[sectorElectron]->Fill(W);
								
								beamEnergyCalc = returnEBeamFromThetaEThetaP(p4_ele[i], p4_proton[j]);
								eBeamFromAnglesAsThetaP[sectorElectron]->Fill(thetaPMeasured, beamEnergyCalc - eBeam);
								eBeamFromAnglesAsThetaE[sectorElectron]->Fill(thetaEMeasured, beamEnergyCalc - eBeam);
								eBeamFromAnglesAsPP[sectorElectron]->Fill(pPMeasured, beamEnergyCalc - eBeam);
								eBeamFromAnglesAsPE[sectorElectron]->Fill(pEMeasured, beamEnergyCalc - eBeam);
								eBeamFromAngles[sectorElectron]->Fill(beamEnergyCalc - eBeam);

								//theta E From different variables
								thetaEFromEPValue = returnThetaEFromPE(eBeam, p4_ele[i]);									

								deltaPxPx[sectorElectron]->Fill(p4_ele[i].Px(), p4_ele_final.Px() - p4_ele[i].Px()); 
								deltaPyPy[sectorElectron]->Fill(p4_ele[i].Py(), p4_ele_final.Py() - p4_ele[i].Py()); 
								deltaPzPz[sectorElectron]->Fill(p4_ele[i].Pz(), p4_ele_final.Pz() - p4_ele[i].Pz()); 

								deltaThetaEFromEPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaEFromEPValue + thetaEMeasured); 
								deltaThetaEFromEPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaEFromEPValue + thetaEMeasured);
								deltaThetaEFromEPAsPE[sectorElectron]->Fill(pEMeasured, -thetaEFromEPValue + thetaEMeasured);
								deltaThetaEFromEPAsPP[sectorElectron]->Fill(pPMeasured, -thetaEFromEPValue + thetaEMeasured);
								deltaThetaEFromEP[sectorElectron]->Fill(-thetaEFromEPValue + thetaEMeasured);

								thetaEFromEPValueEnergyConserve = returnThetaEFromPE(eBeam, p4_ele_final);	
								deltaThetaEFromEPAsThetaEEnergyConservation[sectorElectron]->Fill(thetaEMeasured, -thetaEFromEPValueEnergyConserve + thetaEMeasured); 
								deltaThetaEFromEPAsThetaPEnergyConservation[sectorElectron]->Fill(thetaPMeasured, -thetaEFromEPValueEnergyConserve + thetaEMeasured);
								deltaThetaEFromEPAsPEEnergyConservation[sectorElectron]->Fill(pEMeasured, -thetaEFromEPValueEnergyConserve + thetaEMeasured);
								deltaThetaEFromEPAsPPEnergyConservation[sectorElectron]->Fill(pPMeasured, -thetaEFromEPValueEnergyConserve + thetaEMeasured);
								deltaThetaEFromEPEnergyConservation[sectorElectron]->Fill(-thetaEFromEPValueEnergyConserve + thetaEMeasured);

								thetaEFromPPValue = returnThetaEFromPProton(eBeam, p4_proton[j]);
								deltaThetaEFromPPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaEFromPPValue + thetaEMeasured); 
								deltaThetaEFromPPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaEFromPPValue + thetaEMeasured);
								deltaThetaEFromPPAsPE[sectorElectron]->Fill(pEMeasured, -thetaEFromPPValue + thetaEMeasured);
								deltaThetaEFromPPAsPP[sectorElectron]->Fill(pPMeasured, -thetaEFromPPValue + thetaEMeasured);
								deltaThetaEFromPP[sectorElectron]->Fill(-thetaEFromPPValue + thetaEMeasured);

								thetaEFromThetaPValue = returnThetaEFromThetaProton(eBeam, p4_proton[i]);
								deltaThetaEFromPThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaEFromThetaPValue + thetaEMeasured); 
								deltaThetaEFromPThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaEFromThetaPValue + thetaEMeasured);
								deltaThetaEFromPThetaAsPE[sectorElectron]->Fill(pEMeasured, -thetaEFromThetaPValue + thetaEMeasured);
								deltaThetaEFromPThetaAsPP[sectorElectron]->Fill(pPMeasured, -thetaEFromThetaPValue + thetaEMeasured);
								deltaThetaEFromPTheta[sectorElectron]->Fill(-thetaEFromThetaPValue + thetaEMeasured);


								//theta P From different variables
								thetaPFromEPValue = returnThetaPFromPE(eBeam, p4_ele[i]);
								deltaThetaPFromEPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaPFromEPValue + thetaPMeasured); 
								deltaThetaPFromEPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaPFromEPValue + thetaPMeasured);
								deltaThetaPFromEPAsPE[sectorElectron]->Fill(pEMeasured, -thetaPFromEPValue + thetaPMeasured);
								deltaThetaPFromEPAsPP[sectorElectron]->Fill(pPMeasured, -thetaPFromEPValue + thetaPMeasured);
								deltaThetaPFromEP[sectorElectron]->Fill(-thetaPFromEPValue + thetaPMeasured);

								thetaPFromPPValue = returnThetaPFromPP(eBeam, p4_proton[j]);
								deltaThetaPFromPPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaPFromPPValue + thetaPMeasured); 
								deltaThetaPFromPPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaPFromPPValue + thetaPMeasured);
								deltaThetaPFromPPAsPE[sectorElectron]->Fill(pEMeasured, -thetaPFromPPValue + thetaPMeasured);
								deltaThetaPFromPPAsPP[sectorElectron]->Fill(pPMeasured, -thetaPFromPPValue + thetaPMeasured);
								deltaThetaPFromPP[sectorElectron]->Fill(-thetaPFromPPValue + thetaPMeasured);

								thetaPFromThetaEValue = returnThetaPFromThetaE(eBeam, p4_ele[i]);
								deltaThetaPFromEThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured, -thetaPFromThetaEValue + thetaPMeasured); 
								deltaThetaPFromEThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured, -thetaPFromThetaEValue + thetaPMeasured);
								deltaThetaPFromEThetaAsPE[sectorElectron]->Fill(pEMeasured, -thetaPFromThetaEValue + thetaPMeasured);
								deltaThetaPFromEThetaAsPP[sectorElectron]->Fill(pPMeasured, -thetaPFromThetaEValue + thetaPMeasured);
								deltaThetaPFromETheta[sectorElectron]->Fill(-thetaPFromThetaEValue + thetaPMeasured);

								//p E From different variables
								pEFromThetaEValue = returnPEFromThetaE(eBeam, p4_ele[i]);
								deltaPEFromPThetaAsPhiE->Fill(phiEMeasured, -pEFromThetaEValue + pEMeasured);
								deltaPEFromEThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured,-pEFromThetaEValue + pEMeasured);
								deltaPEFromEThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured,-pEFromThetaEValue + pEMeasured);
								deltaPEFromEThetaAsPE[sectorElectron]->Fill(pEMeasured, -pEFromThetaEValue + pEMeasured);
								deltaPEFromEThetaAsPP[sectorElectron]->Fill(pPMeasured, -pEFromThetaEValue + pEMeasured);
								deltaPEFromETheta[sectorElectron]->Fill(-pEFromThetaEValue + pEMeasured);

								pEFromPPValue = returnPEFromPProton(eBeam, p4_proton[j]);
								deltaPEFromPPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -pEFromPPValue + pEMeasured); 
								deltaPEFromPPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -pEFromPPValue + pEMeasured);
								deltaPEFromPPAsPE[sectorElectron]->Fill(pEMeasured, -pEFromPPValue + pEMeasured);
								deltaPEFromPPAsPP[sectorElectron]->Fill(pPMeasured, -pEFromPPValue + pEMeasured);
								deltaPEFromPP[sectorElectron]->Fill(-pEFromPPValue + pEMeasured);

								pEFromThetaPValue = returnPEFromThetaP(eBeam, p4_proton[j]);
								deltaPEFromPThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured, -pEFromThetaPValue + pEMeasured); 
								deltaPEFromPThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured, -pEFromThetaPValue + pEMeasured);
								deltaPEFromPThetaAsPE[sectorElectron]->Fill(pEMeasured, -pEFromThetaPValue + pEMeasured);
								deltaPEFromPThetaAsPP[sectorElectron]->Fill(pPMeasured, -pEFromThetaPValue + pEMeasured);
								deltaPEFromPTheta[sectorElectron]->Fill(-pEFromThetaPValue + pEMeasured);

								//p P From different variables
								pPFromThetaEValue = returnPPFromThetaE(eBeam, p4_ele[i]);
								deltaPPFromEThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured,-pPFromThetaEValue + pPMeasured);
								deltaPPFromEThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured,-pPFromThetaEValue + pPMeasured);
								deltaPPFromEThetaAsPE[sectorElectron]->Fill(pEMeasured, -pPFromThetaEValue + pPMeasured);
								deltaPPFromEThetaAsPP[sectorElectron]->Fill(pPMeasured, -pPFromThetaEValue + pPMeasured);
								deltaPPFromETheta[sectorElectron]->Fill(-pPFromThetaEValue + pPMeasured);

								pPFromEPValue = returnPPromEP(eBeam, p4_ele[i]);
								deltaPPFromEPAsThetaE[sectorElectron]->Fill(thetaEMeasured, -pPFromEPValue + pPMeasured); 
								deltaPPFromEPAsThetaP[sectorElectron]->Fill(thetaPMeasured, -pPFromEPValue + pPMeasured);
								deltaPPFromEPAsPE[sectorElectron]->Fill(pEMeasured, -pPFromEPValue + pPMeasured);
								deltaPPFromEPAsPP[sectorElectron]->Fill(pPMeasured, -pPFromEPValue + pPMeasured);
								deltaPPFromEP[sectorElectron]->Fill(-pPFromEPValue + pPMeasured);

								pPFromThetaPValue = returnPPromThetaP(eBeam, p4_proton[j]);
								deltaPPFromPThetaAsThetaE[sectorElectron]->Fill(thetaEMeasured, -pPFromThetaPValue + pPMeasured); 
								deltaPPFromPThetaAsThetaP[sectorElectron]->Fill(thetaPMeasured, -pPFromThetaPValue + pPMeasured);
								deltaPPFromPThetaAsPE[sectorElectron]->Fill(pEMeasured, -pPFromThetaPValue + pPMeasured);
								deltaPPFromPThetaAsPP[sectorElectron]->Fill(pPMeasured, -pPFromThetaPValue + pPMeasured);
								deltaPPFromPTheta[sectorElectron]->Fill(-pPFromThetaPValue + pPMeasured);


							}

						}
					}
				}
			}
			if (k%10000 == 0) cout << k <<endl;
			if(k == nEvents) break;
		}
	}
	resultsF->cd("overview");
	wQ2Inclusive->Write();
	wQ2Exclusive ->Write();
	zVsPhiElectron->Write();
	zVsPhiProton->Write();
	deltaZVsPhiElectron->Write();
	deltaZVsPhiProton->Write();
	wVsPhi->Write();

	for (int s = 0; s < 6; s++){
		//inclusive
		wQ2SectorInclusive[s]->Write();
		wSectorInclusive[s]->Write();

		//exclusive
		wQ2SectorExclusive[s]->Write(); 
		wSectorExclusive[s]->Write();

		phiEphiP[s]->Write();
		phiEPhiWCut[s]->Write();
		wSectorExclusivePhiCut[s]->Write();
		zElectronSector[s]->Write();
		zProtonSector[s]->Write();
		deltaZElectronProtonSector[s]->Write();

		zVsThetaElectron[s]->Write();
		zVsThetaProton[s]->Write();
		for (int t = 0; t < 6; t++){
			thetaEPWSThetaBins[s][t]->Write();
		}
	}
	resultsF->cd("eBeam");
	for (int s = 0; s < 6; s++){
		eBeamFromAnglesAsThetaP[s]->Write();
		eBeamFromAnglesAsThetaE[s]->Write();
		eBeamFromAnglesAsPP[s]->Write();
		eBeamFromAnglesAsPE[s]->Write();
		eBeamFromAngles[s]->Write();
	}

	resultsF->cd("electronInclusive");
	for (int s = 0; s < 6; s++){
		deltaThetaEFromEPAsThetaEInclusive[s]->Write(); 
		deltaThetaEFromEPAsPEInclusive[s]->Write(); 
		deltaThetaEFromEPInclusive[s]->Write(); 
		deltaPEFromEThetaAsThetaEInclusive[s]->Write(); 
		deltaPEFromEThetaAsPEInclusive[s]->Write(); 
		deltaPEFromEThetaInclusive[s]->Write();
		thetaPElastic[s]->Write();
	}
	resultsF->cd("energyConservation");
	for (int s = 0; s < 6; s++){
		deltaPxPx[s]->Write(); 
		deltaPyPy[s]->Write(); 
		deltaPzPz[s]->Write(); 
	}
	resultsF->cd("thetaElectron");
	for (int s = 0; s < 6; s++){
		deltaThetaEFromPPAsThetaE[s]->Write(); 
		deltaThetaEFromPPAsThetaP[s]->Write();
		deltaThetaEFromPPAsPE[s]->Write();
		deltaThetaEFromPPAsPP[s]->Write();
		deltaThetaEFromPP[s]->Write();

		deltaThetaEFromEPAsThetaE[s]->Write();
		deltaThetaEFromEPAsThetaP[s]->Write();
		deltaThetaEFromEPAsPE[s]->Write();
		deltaThetaEFromEPAsPP[s]->Write();
		deltaThetaEFromEP[s]->Write();


		deltaThetaEFromEPAsThetaEEnergyConservation[s]->Write();
		deltaThetaEFromEPAsThetaPEnergyConservation[s]->Write();
		deltaThetaEFromEPAsPEEnergyConservation[s]->Write();
		deltaThetaEFromEPAsPPEnergyConservation[s]->Write();
		deltaThetaEFromEPEnergyConservation[s]->Write();

		deltaThetaEFromPThetaAsThetaE[s]->Write();
		deltaThetaEFromPThetaAsThetaP[s]->Write();
		deltaThetaEFromPThetaAsPE[s]->Write();
		deltaThetaEFromPThetaAsPP[s]->Write();
		deltaThetaEFromPTheta[s]->Write();
	}

	resultsF->cd("thetaProton");
	for (int s = 0; s < 6; s++){
		deltaThetaPFromPPAsThetaE[s]->Write(); 
		deltaThetaPFromPPAsThetaP[s]->Write();
		deltaThetaPFromPPAsPE[s]->Write();
		deltaThetaPFromPPAsPP[s]->Write();
		deltaThetaPFromPP[s]->Write();

		deltaThetaPFromEPAsThetaE[s]->Write();
		deltaThetaPFromEPAsThetaP[s]->Write();
		deltaThetaPFromEPAsPE[s]->Write();
		deltaThetaPFromEPAsPP[s]->Write();
		deltaThetaPFromEP[s]->Write();

		deltaThetaPFromEThetaAsThetaE[s]->Write();
		deltaThetaPFromEThetaAsThetaP[s]->Write();
		deltaThetaPFromEThetaAsPE[s]->Write();
		deltaThetaPFromEThetaAsPP[s]->Write();
		deltaThetaPFromETheta[s]->Write();
	}

	resultsF->cd("pElectron");
	deltaPEFromPThetaAsPhiE->Write();
	for (int s = 0; s < 6; s++){
		deltaPEFromPPAsThetaE[s]->Write(); 
		deltaPEFromPPAsThetaP[s]->Write();
		deltaPEFromPPAsPE[s]->Write();
		deltaPEFromPPAsPP[s]->Write();
		deltaPEFromPP[s]->Write();

		deltaPEFromPThetaAsThetaE[s]->Write();
		deltaPEFromPThetaAsThetaP[s]->Write();
		deltaPEFromPThetaAsPE[s]->Write();
		deltaPEFromPThetaAsPP[s]->Write();
		deltaPEFromPTheta[s]->Write();

		deltaPEFromEThetaAsThetaE[s]->Write();
		deltaPEFromEThetaAsThetaP[s]->Write();
		deltaPEFromEThetaAsPE[s]->Write();
		deltaPEFromEThetaAsPP[s]->Write();
		deltaPEFromETheta[s]->Write();
	}

	resultsF->cd("pProton");
	for (int s = 0; s < 6; s++){
		deltaPPFromPThetaAsThetaE[s]->Write(); 
		deltaPPFromPThetaAsThetaP[s]->Write();
		deltaPPFromPThetaAsPE[s]->Write();
		deltaPPFromPThetaAsPP[s]->Write();
		deltaPPFromPTheta[s]->Write();

		deltaPPFromEPAsThetaE[s]->Write();
		deltaPPFromEPAsThetaP[s]->Write();
		deltaPPFromEPAsPE[s]->Write();
		deltaPPFromEPAsPP[s]->Write();
		deltaPPFromEP[s]->Write();

		deltaPPFromEThetaAsThetaE[s]->Write();
		deltaPPFromEThetaAsThetaP[s]->Write();
		deltaPPFromEThetaAsPE[s]->Write();
		deltaPPFromEThetaAsPP[s]->Write();
		deltaPPFromETheta[s]->Write();
	}


	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
	// write out fit limits for electronProtonAna
	ofstream outputLimits1;
	std::string writeDeltapefrometheta = "deltapefrometheta_fit_limits_job"+std::to_string(model)+".txt";
	std::cout << " Creating File " << writeDeltapefrometheta << std::endl;
	outputLimits1.open(writeDeltapefrometheta);

	for( int s = 0; s < 6; s++ ){    
	  TF1 *fittemp=fitHisto(deltaPEFromETheta[s]);
	  std::cout << " Sector s: " << s << ", " << fittemp->GetParameter(1) << ", " << fittemp->GetParameter(2);    
	  outputLimits1 << s << " " << fittemp->GetParameter(1) << " " << fittemp->GetParameter(2) << std::endl;
	}

	ofstream outputLimits2;
	std::string writeDeltaPPFromETheta = "deltaPPFromETheta_fit_limits_job"+std::to_string(model)+".txt";
	std::cout << " Creating File " << writeDeltaPPFromETheta << std::endl;
	outputLimits2.open(writeDeltaPPFromETheta);

	for( int s = 0; s < 6; s++ ){    
	  TF1 *fittemp=fitHisto(deltaPPFromETheta[s]);
	  std::cout << " Sector s: " << s << ", " << fittemp->GetParameter(1) << ", " << fittemp->GetParameter(2);    
	  outputLimits2 << s << " " << fittemp->GetParameter(1) << " " << fittemp->GetParameter(2) << std::endl;
	}


	ofstream outputLimits3;
	std::string writeDeltaThetaPFromETheta = "deltaThetaPFromETheta_fit_limits_job"+std::to_string(model)+".txt";
	std::cout << " Creating File " << writeDeltaThetaPFromETheta << std::endl;
	outputLimits3.open(writeDeltaThetaPFromETheta);

	for( int s = 0; s < 6; s++ ){    
	  TF1 *fittemp=fitHisto(deltaThetaPFromETheta[s]);
	  std::cout << " Sector s: " << s << ", " << fittemp->GetParameter(1) << ", " << fittemp->GetParameter(2);    
	  outputLimits3 << s << " " << fittemp->GetParameter(1) << " " << fittemp->GetParameter(2) << std::endl;
	}


	resultsF->Close();


	return 0;
}
//electron momentum

float returnPEFromThetaE(float Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	return Ebeam/(1+(2*Ebeam*pow(sin(p4_electron.Theta())/2,2))/mProton);
}
float returnPEFromThetaP(double Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	double thetaE = 2*atan(mProton/((mProton+Ebeam)*tan(p4_proton.Theta())));
	return Ebeam/(1 + Ebeam*(1 - cos(thetaE))/mProton);
}
float returnPEFromPProton(float Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	return Ebeam + mProton - sqrt(p4_proton.P()*p4_proton.P() + mProton*mProton);
}

//electron theta
float returnThetaEFromPE(float Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	if (sqrt((mProton/(2*p4_electron.P()))-(mProton/(2*Ebeam)))< 1 && sqrt((mProton/(2*p4_electron.P()))-(mProton/(2*Ebeam))) > -1){
		return 57.2958*2 * asin(sqrt((mProton/(2*p4_electron.P()))-(mProton/(2*Ebeam))));
	}
	else return -1000;
}

float returnThetaEFromPProton(float Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	double pE = Ebeam + mProton - sqrt(p4_proton.P()*p4_proton.P() + mProton*mProton);
	if (sqrt((mProton/(2*pE))-(mProton/(2*Ebeam)))< 1 && sqrt((mProton/(2*pE))-(mProton/(2*Ebeam))) > -1){
		return 57.2958*2*asin(sqrt((mProton/(2*pE))-(mProton/(2*Ebeam))));
	}
	else return -1000;
}

float returnThetaEFromThetaProton(float Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	return 2*atan(mProton/((mProton+Ebeam)*tan(p4_proton.Theta())))*57.2958;
}

//proton theta
float returnThetaPFromThetaE(float Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	return atan(mProton/((mProton+Ebeam)*tan(0.5*p4_electron.Theta())))*57.2958;
}
float returnThetaPFromPE(float Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	double thetaE = 2 * asin(sqrt((mProton/(2*p4_electron.P()))-(mProton/(2*Ebeam))));
	return atan(mProton/((mProton+Ebeam)*tan(0.5*thetaE)))*57.2958;
}
float returnThetaPFromPP(float Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	double pE = Ebeam + mProton - sqrt(p4_proton.P()*p4_proton.P() + mProton*mProton);
	if (sqrt((1/(2*pE))*(mProton - pE*mProton/Ebeam))> -1 && sqrt((1/(2*pE))*(mProton - pE*mProton/Ebeam)) < 1){
		double thetaE = 2*asin(sqrt((mProton/(2*pE))-(mProton/(2*Ebeam))));
		return atan(mProton/((mProton+Ebeam)*tan(0.5*thetaE)))*57.2958;
	}
	else return -1000;
}


//proton momentum
float returnPPFromThetaE(double Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	double pE = Ebeam/(1 + 2*Ebeam/mProton*sin(p4_electron.Theta()/2)*sin(p4_electron.Theta()/2));
	return sqrt((Ebeam + mProton - pE)*(Ebeam + mProton - pE)- mProton*mProton);
}
float returnPPromThetaP(float Ebeam, TLorentzVector p4_proton){
	double mProton = 0.93827;
	double thetaE = 2*atan(mProton/((mProton+Ebeam)*tan(p4_proton.Theta())));
	double pE = Ebeam/(1 + 2*Ebeam/mProton*sin(thetaE/2)*sin(thetaE/2));
	return sqrt((Ebeam + mProton - pE)*(Ebeam + mProton - pE)- mProton*mProton);
}

float returnPPromEP(double Ebeam, TLorentzVector p4_electron){
	double mProton = 0.93827;
	double pE = p4_electron.P();
	return sqrt((Ebeam + mProton - pE)*(Ebeam + mProton - pE)- mProton*mProton);
}


//beam From angles
float returnEBeamFromThetaEThetaP(TLorentzVector p4_electron, TLorentzVector p4_proton){
	double mProton = 0.93827;
	return mProton/(tan(p4_electron.Theta()/2)*tan(p4_proton.Theta())) - mProton;
}





double kin_W(TLorentzVector ele, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	TLorentzVector fCM = fGamma + target;
	return fCM.M();
}


double kin_mismass2(TLorentzVector ele, TLorentzVector hadron, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector miss  = beam + target - ele - hadron;
	return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}


double kin_pmismass2(TLorentzVector hadron, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector miss   = beam + target - hadron;
	return (miss.E()*miss.E() - miss.Px()*miss.Px() - miss.Py()*miss.Py() - miss.Pz()*miss.Pz());
}


double kin_Q2(TLorentzVector ele, float Ebeam){
	TLorentzVector beam(0,0,Ebeam,Ebeam);
	TLorentzVector target(0,0,0,0.93827);
	TLorentzVector fGamma = beam - ele;
	return -fGamma.M2();
}

double returnPhi(double phi){
	float phiLocal = phi;
	if (phi < -20) phiLocal = phi + 360;
	return phiLocal;
}

TF1* fitHisto(TH1D* htemp){
  //////////////////////////////////////////////////
  // Start Andrew's fitting method:
  double xlow,xhigh,histmax;
  int binlow,binhigh,binmax;
  binmax = htemp->GetMaximumBin();
  histmax = htemp->GetMaximum();
  binlow=binmax;
  binhigh=binmax;

  // The 0.65 parameter can be changed, this basically means start at the peak and work your way left and right
  // until you've gotten to 65% of the max amplitude.
  while(htemp->GetBinContent(binhigh++) >= .65*histmax&&binhigh<=htemp->GetNbinsX()){};
  while(htemp->GetBinContent(binlow--) >= .65*histmax&&binlow>=1){};
    
  xlow = htemp->GetBinLowEdge(binlow);
  xhigh = htemp->GetBinLowEdge(binhigh+1);
    
  htemp->Fit("gaus","","",xlow,xhigh);

  TF1 *ftemp = (TF1*) htemp->GetListOfFunctions()->FindObject("gaus");

  return ftemp;

  // End
  /////////////////////////////////////////////

}

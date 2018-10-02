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
#include "TLatex.h"

using namespace std;

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

using namespace ROOT::Math;
void electronProtonAnaShort(int model){
       //1011 s-06T-0.6
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
	gStyle->SetOptTitle(1);
	gStyle->SetPadTopMargin(0.1);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.1);
	gStyle->SetOptStat(10);
	gStyle->SetFrameBorderMode(0);
	gStyle->cd();
	gROOT->ForceStyle();
	TLatex sCaption;
	sCaption.SetNDC();
	sCaption.SetTextFont(42);
	sCaption.SetTextSize(0.05);
	gStyle->SetOptFit(10);
	gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
	gStyle->SetStatW(0.15);
	gStyle->SetStatH(0.15); 
//fit functions
	TF1 *elasGausPeak[6]; 
	TF1 *beamGausPeak[6]; 
	TF1 *deltaPhiFit[6];
	TF1 *deltaPPFit[6];
	TF1 *deltaEPFit[6];
	TF1 *deltaPThetaFit[6];
	TF1 *deltaEThetaFit[6];

	float lowWValueFit[6];
	float highWValueFit[6];

	float lowDeltaPEValueFit[6];// = {-0.02, -0.02, -0.02, -0.02, -0.02, -0.02};
	float highDeltaPEValueFit[6];// = {0.002, 0.002, 0.002, 0.002, 0.002, 0.002};

	float lowDeltaPPFromThetaEValueFit[6];
	float highDeltaPPFromThetaEValueFit[6];

	float lowDeltaThetaPValueFit[6];
	float highDeltaThetaPValueFit[6];

 	int sector0, sector1, sector2, sector3;
	double mean0, mean1, mean2, mean3;
	double sig0, sig1, sig2, sig3;
	
	string line0;
	std::string parentDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/markovMonitor/scripts/quick_fit_parameters/";

	ifstream readFromWCut(parentDirectory+"w_cut_limits_job"+std::to_string(model)+".txt");
	if( readFromWCut.is_open() ){
	  while(readFromWCut >> sector0 ) {//std::getline (readFromWCut, line) ){
	    
	    readFromWCut >> mean0 >> sig0;
	    std::cout << " >> W CUT PARAMETERS: " << sector0 << " " << mean0 << " " << sig0 << std::endl;
	    
	    lowWValueFit[sector0]= mean0 - sig0;
	    highWValueFit[sector0]= mean0 + sig0;
	  }
	}
	
 	//string line;
	ifstream readFromdeltaPEFromETheta(parentDirectory+"deltapefrometheta_fit_limits_job"+std::to_string(model)+".txt");
	if( readFromdeltaPEFromETheta.is_open() ){
	  while(readFromdeltaPEFromETheta >> sector1 ) {
	    
	    readFromdeltaPEFromETheta >> mean1 >> sig1;
	    std::cout << " >> deltaPeFromETheta FIT PARAMETERS: " << sector1 << " " << mean1 << " " << sig1 << std::endl;
	    lowDeltaPEValueFit[sector1]= mean1 - sig1;
	    highDeltaPEValueFit[sector1]= mean1 + sig1;
	  }
	}

	//special sector fit range for outbending run 2587...
	//lowDeltaPEValueFit[3]= -0.025;
	//highDeltaPEValueFit[3]= 0.025;

	ifstream readFromdeltaPPFromETheta(parentDirectory+"deltaPPFromETheta_fit_limits_job"+std::to_string(model)+".txt");
	if( readFromdeltaPPFromETheta.is_open() ){
	  while(readFromdeltaPPFromETheta >> sector2 ) {
	    
	    readFromdeltaPPFromETheta >> mean2 >> sig2;
	    std::cout << " >> deltaPPFromETheta FIT PARAMETERS: " << sector2 << " " << mean2 << " " << sig2 << std::endl;
	    lowDeltaPPFromThetaEValueFit[sector2]= mean2 - sig2;
	    highDeltaPPFromThetaEValueFit[sector2]= mean2 + sig2;
	  }
	}

	//lowDeltaPPFromThetaEValueFit[3] = -0.05;//
	//highDeltaPPFromThetaEValueFit[3] = 0.1;//highDeltaPPFromThetaEValueFit[3]*1.1;

	ifstream readFromdeltaThetaPFromETheta(parentDirectory+"deltaThetaPFromETheta_fit_limits_job"+std::to_string(model)+".txt");
	if( readFromdeltaThetaPFromETheta.is_open() ){
	  while(readFromdeltaThetaPFromETheta >> sector3 ) {
	    
	    readFromdeltaThetaPFromETheta >> mean3 >> sig3;
	    std::cout << " >> deltaThetaPFromETheta FIT PARAMETERS: " << sector3 << " " << mean3 << " " << sig3 << std::endl;
	    lowDeltaThetaPValueFit[sector3]= mean3 - sig3;
	    highDeltaThetaPValueFit[sector3]= mean3 + sig3;
	  }
	}
	

//fit parameters
//fit parameters for job 6 (outbending) this needs to be automated later
	float meanw1 = 0.964894; float sigw1 =  0.0430203;
	float meanw2 = 0.944595; float sigw2 =  0.0301894;
	float meanw3 = 0.938689; float sigw3 =  0.0521689;
	float meanw4 = 0.941072; float sigw4 =  0.029401;
	float meanw5 = 1.00682; float sigw5 =  0.0302363;
	float meanw6 =  1.00228; float sigw6 =  0.0386132;

	float meanphi1 = 182.027; float sigphi1 = 2* 1.19124;
	float meanphi2 = 181.946; float sigphi2 = 2*1.44038;
	float meanphi3 = 179.568; float sigphi3 =  2*1.54789;
	float meanphi4 = 182.02; float sigphi4 =  2*1.15751;
	float meanphi5 = 181.214; float sigphi5 = 2*0.883448;
	float meanphi6 = 179.029; float sigphi6 =  2*0.939062;

	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	//uncomment below out for inbending runs ////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////


	//float lowWValueFit[6] = {1.005, 1.01, 1.01, 1.005, 0.99, 0.995};
	//float highWValueFit[6] = {1.04, 1.055, 1.045, 1.04, 1.02, 1.03};
	/*
	float lowEBeamValueFit[6] = {-0.4, -0.5, -0.4, -0.4, -0.3, -0.4};
	float highEBeamValueFit[6] = {0.05, -0, 0.0, 0.05, 0.1, 0.05};
		
	float lowPhiValueFit[6] = {177, 178, 175, 175, 176, 175};
	float highPhiValueFit[6] = {187, 187, 185, 185, 185, 184};
	
	float lowDeltaThetaPFromEPValueFit[6] = {5, 5, 5, 4.5, 3.5, 4};
	float highDeltaThetaPFromEPValueFit[6] = {10, 12, 11, 11, 9, 9};
	
	float lowDeltaThetaPValueFit[6] = {-1, -0.5, -0.5, 0, -1, -0.5};
	float highDeltaThetaPValueFit[6] = {3, 4, 3, 3, 2.5, 2.5};
	
	float lowDeltaThetaEValueFit[6] = {-0.5, -0, -0.5, -0.5, -1, -0.5};
	float highDeltaThetaEValueFit[6] = {2., 2.5, 2., 2., 1.8, 2};
	
	
	float lowDeltaThetaEFromEPValueFit[6] = {-6, -6, -6, -6, -5, -5.5};
	float highDeltaThetaEFromEPValueFit[6] = {-3.5, -3.5, -4, -4, -3, -3};

	
	float lowDeltaPEValueFromThetaPFit[6] = {-0.12, -0.15, -0.14, -0.14, -0.1, -0.12};
	float highDeltaPEValueFromThetaPFit[6] = {-0.06, -0.07, -0.07, -0.06, -0.04, -0.05};
	
	float lowDeltaPEValueFit[6] = {-0.1, -0.12, -0.11, -0.1, -0.08, -0.1};
	float highDeltaPEValueFit[6] = {-0.07, -0.07, -0.07, -0.07, -0.05, -0.05};		
	
	float lowDeltaPPValueFit[6] = {-0.1, -0.1, -0.07, -0.13, -0.06, -0.07};
	float highDeltaPPValueFit[6] = {0.05, 0.03, 0.05, 0.14, 0.08, 0.07};
	
	float lowDeltaPPFromPEValueFit[6] = {-0.3, -0.35, -0.28, -0.28, -0.23, -0.26};
	float highDeltaPPFromPEValueFit[6] = {-0.18, -0.22, -0.2, -0.2, -0.12, -0.14};
	
	
	float lowDeltaPPFromThetaEValueFit[6] = {-0.08, -0.12, -0.08, -0.07, -0.065, -0.075};
	float highDeltaPPFromThetaEValueFit[6] = {-0.02, -0.05, -0.01, -0.01, 0.03, 0.01};
	*/
	if (model == 11){
		float lowWValueFit[6] = {0.935, 0.935, 0.935, 0.935, 0.935, 0.935};
		float highWValueFit[6] = {0.96, 0.96, 0.96, 0.96, 0.96, 0.96};


		float lowEBeamValueFit[6] = {-0.15, -0.15, -0.15, -0.15, -0.15, -0.15};
		float highEBeamValueFit[6] = {0.22, 0.23, 0.23, 0.23, 0.22, 0.22};

		float lowPhiValueFit[6] = {177, 178, 175, 175, 176, 175};
		float highPhiValueFit[6] = {187, 187, 185, 185, 185, 184};

		float lowDeltaThetaPValueFit[6] = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
		float highDeltaThetaPValueFit[6] = {1.2, 1.2, 1.2, 1.2, 1.2, 1.2};


		float lowDeltaThetaEValueFit[6] = {-1., -1.2, -1.2, -1.2, -1.2, -1.2};
		float highDeltaThetaEValueFit[6] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8};


		float lowDeltaThetaEFromEPValueFit[6] = {-1.2, -1.2, -1.2, -1.2, -1.2, -1.2};
		float highDeltaThetaEFromEPValueFit[6] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
		
		float lowDeltaPEValueFit[6] = {-0.02, -0.02, -0.02, -0.02, -0.02, -0.02};
		float highDeltaPEValueFit[6] = {0.002, 0.002, 0.002, 0.002, 0.002, 0.002};


		float lowDeltaPEValueFromThetaPFit[6] = {-0.04, -0.04, -0.04, -0.04, -0.04, -0.04};
		float highDeltaPEValueFromThetaPFit[6] = {0.025, 0.02, 0.02, 0.02, 0.02, 0.02};

		float lowDeltaThetaPFromEPValueFit[6] = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
		float highDeltaThetaPFromEPValueFit[6] = {2., 2., 2., 2., 2., 2.};

		float lowDeltaPPValueFit[6] = {-0.06, -0.06, -0.06, -0.06, -0.06, -0.06};
		float highDeltaPPValueFit[6] = {0.07, 0.07, 0.07, 0.07, 0.07, 0.07};

		float lowDeltaPPFromPEValueFit[6] ={-0.07, -0.07, -0.07, -0.07, -0.07, -0.07};
		float highDeltaPPFromPEValueFit[6] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

		//		float lowDeltaPPFromThetaEValueFit[6] = {0.018, 0.018, 0.015, 0.017, 0.02, 0.02};
		//		float highDeltaPPFromThetaEValueFit[6] = {0.14, 0.14, 0.14, 0.14, 0.14, 0.14};
		float lowDeltaPPFromThetaEValueFit[6] = {-0.04, -0.04, -0.04, -0.04, -0.04, -0.04};                               
		float highDeltaPPFromThetaEValueFit[6] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04};        
	}
	//else if( model == 66666 ){

	  std::cout << " using outbending cut values " << std::endl;
	  

	  float lowEBeamValueFit[6] = {-0.15, -0.15, -0.15, -0.15, -0.15, -0.15};
	  float highEBeamValueFit[6] = {0.22, 0.23, 0.23, 0.23, 0.22, 0.22};

	  float lowPhiValueFit[6] = {meanphi1 - sigphi1, meanphi2 - sigphi2, meanphi3 - sigphi3, meanphi4 - sigphi4, meanphi5 - sigphi5, meanphi6 - sigphi6 };
	  float highPhiValueFit[6] = {meanphi1 + sigphi1, meanphi2 + sigphi2, meanphi3 + sigphi3, meanphi4 + sigphi4, meanphi5 + sigphi5, meanphi6 + sigphi6 };


	  //commenting out to instead read from txt file 
	  //float lowDeltaThetaPValueFit[6] = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
	  //float highDeltaThetaPValueFit[6] = {1.2, 1.2, 1.2, 1.2, 1.2, 1.2};


	  float lowDeltaThetaEValueFit[6] = {-1., -1.2, -1.2, -1.2, -1.2, -1.2};
	  float highDeltaThetaEValueFit[6] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8};


	  float lowDeltaThetaEFromEPValueFit[6] = {-1.2, -1.2, -1.2, -1.2, -1.2, -1.2};
	  float highDeltaThetaEFromEPValueFit[6] = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
		
	  //float lowDeltaPEValueFit[6] = {-0.02, -0.02, -0.02, -0.02, -0.02, -0.02};
	  //float highDeltaPEValueFit[6] = {0.002, 0.002, 0.002, 0.002, 0.002, 0.002};


	  float lowDeltaPEValueFromThetaPFit[6] = {-0.04, -0.04, -0.04, -0.04, -0.04, -0.04};
	  float highDeltaPEValueFromThetaPFit[6] = {0.025, 0.02, 0.02, 0.02, 0.02, 0.02};

	  float lowDeltaThetaPFromEPValueFit[6] = {-1.5, -1.5, -1.5, -1.5, -1.5, -1.5};
	  float highDeltaThetaPFromEPValueFit[6] = {2., 2., 2., 2., 2., 2.};

	  float lowDeltaPPValueFit[6] = {-0.06, -0.06, -0.06, -0.06, -0.06, -0.06};
	  float highDeltaPPValueFit[6] = {0.07, 0.07, 0.07, 0.07, 0.07, 0.07};

	  float lowDeltaPPFromPEValueFit[6] ={-0.07, -0.07, -0.07, -0.07, -0.07, -0.07};
	  float highDeltaPPFromPEValueFit[6] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};


	  //commenting these out for fits
	  //		float lowDeltaPPFromThetaEValueFit[6] = {0.018, 0.018, 0.015, 0.017, 0.02, 0.02};
	  //		float highDeltaPPFromThetaEValueFit[6] = {0.14, 0.14, 0.14, 0.14, 0.14, 0.14};
	  //float lowDeltaPPFromThetaEValueFit[6] = {-0.04, -0.04, -0.04, -0.04, -0.04, -0.04};                               
	  //float highDeltaPPFromThetaEValueFit[6] = {0.04, 0.04, 0.04, 0.04, 0.04, 0.04};        
	  //}

	  //float lowWValueFit[6] = {meanw1 - sigw1, meanw2 - sigw2, meanw3 - sigw3, meanw4 - sigw4, meanw5 - sigw5, meanw6 - sigw6 };
	  //float highWValueFit[6] = {meanw1 + sigw1, meanw2 + sigw2, meanw3 + sigw3, meanw4 + sigw4, meanw5 + sigw5, meanw6 + sigw6 };

	std::cout << " check low values is used " << lowWValueFit[0] << std::endl;

	for (int s = 0; s < 6; s++){
		elasGausPeak[s] = new TF1(Form("elasGausPeakS%d", s + 1), "gaus", lowWValueFit[s], highWValueFit[s]);
		elasGausPeak[s]->SetLineColor(2);

		beamGausPeak[s] = new TF1(Form("beamGausPeakS%d", s + 1), "gaus", -0.3, 0.3);
		beamGausPeak[s]->SetLineColor(2);
		deltaPhiFit[s] = new TF1(Form("deltaPhiFitS%d", s + 1), "gaus", 170, 190);
		deltaPhiFit[s]->SetLineColor(2);
		deltaPPFit[s] = new TF1(Form("deltaPPFitS%d", s + 1), "gaus", -0.5, 0.5);
		deltaPPFit[s]->SetLineColor(2);
		deltaPThetaFit[s] = new TF1(Form("deltaPThetaFitS%d", s + 1), "gaus", -2, 2);
		deltaPThetaFit[s]->SetLineColor(2);
		deltaEThetaFit[s] = new TF1(Form("deltaEThetaFitS%d", s + 1), "gaus", -2, 2);
		deltaEThetaFit[s]->SetLineColor(2);
		deltaEPFit[s] = new TF1(Form("deltaEPFitS%d", s + 1), "gaus", -0.1, 0.1);
		deltaEPFit[s]->SetLineColor(2);
	}

	std::string parentInDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/markovMonitor/scripts/mon_job_root/";
	TFile *resultsF = new TFile(Form("%smonjob%d.root", parentInDirectory.c_str() ,model));

	//overview
		//inclusive
	TH2D *wQ2Inclusive;
	TH2D *wQ2SectorInclusive[6];
	TH1D *wSectorInclusive[6];

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
	TH2F *zVsPhiElectron;
	TH2F *zVsPhiProton;
	TH2F *deltaZVsPhiElectron;
	TH2F *deltaZVsPhiProton;
	TH1D *zElectronSector[6];
	TH1D *zProtonSector[6];
	TH2F *zVsThetaProton[6];
	TH2F *zVsThetaElectron[6];
	TH2F *wVsPhi;


	TH1D *deltaZElectronProtonSector[6];
	//beam energy
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
	TH2D *deltaPEFromPThetaAsPhiE;

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



	//reading histograms
	wQ2Inclusive = (TH2D*) resultsF->Get("overview/wQ2Inclusive");
	wQ2Exclusive = (TH2D*) resultsF->Get("overview/wQ2Exclusive");


	zVsPhiElectron = (TH2F*) resultsF->Get("overview/zVsPhiElectron");
	zVsPhiProton = (TH2F*) resultsF->Get("overview/zVsPhiProton");


	deltaZVsPhiElectron = (TH2F*) resultsF->Get("overview/deltaZVsPhiElectron");
	deltaZVsPhiProton = (TH2F*) resultsF->Get("overview/deltaZVsPhiProton");
	
	wVsPhi= (TH2F*) resultsF->Get("overview/wVsPhi");
	deltaPEFromPThetaAsPhiE = (TH2D*) resultsF->Get("pElectron/deltaPEFromPThetaAsPhiE");
	for (int s = 0; s < 6; s++){
		//inclusive
		wSectorInclusive[s] = (TH1D*) resultsF->Get(Form("overview/wInclusiveS%d", s + 1));

		//exclusive
		wQ2SectorExclusive[s] = (TH2D*) resultsF->Get(Form("overview/wQ2ExclusiveS%d", s + 1));
		wSectorExclusive[s] = (TH1D*) resultsF->Get(Form("overview/wExclusiveS%d", s + 1));

	//	phiEphiP[s] = (TH1D*) resultsF->Get(Form("overview/phiEphiPS%d", s + 1));
	//	phiEPhiWCut[s] = (TH1D*) resultsF->Get(Form("overview/phiEPhiWCutS%d", s + 1));
		wSectorExclusivePhiCut[s] = (TH1D*) resultsF->Get(Form("overview/wSectorExclusivePhiCutS%d", s + 1));
		zElectronSector[s] = (TH1D*) resultsF->Get(Form("overview/zElectronSectorS%d", s + 1));
		zProtonSector[s] = (TH1D*) resultsF->Get(Form("overview/zProtonSectorS%d", s + 1));
		deltaZElectronProtonSector[s] = (TH1D*) resultsF->Get(Form("overview/deltaZElectronProtonSectorS%d", s + 1));

		zVsThetaProton[s] = (TH2F*) resultsF->Get(Form("overview/zVsThetaProtonS%d", s + 1));
		zVsThetaElectron[s] = (TH2F*) resultsF->Get(Form("overview/zVsThetaElectronS%d", s + 1));

		//beam energy

		eBeamFromAnglesAsThetaP[s] = (TH2D*) resultsF->Get(Form("eBeam/eBeamFromAnglesAsThetaPS%d", s + 1));
		eBeamFromAnglesAsThetaE[s] = (TH2D*) resultsF->Get(Form("eBeam/eBeamFromAnglesAsThetaES%d", s + 1));
		eBeamFromAnglesAsPP[s] = (TH2D*) resultsF->Get(Form("eBeam/eBeamFromAnglesAsPPS%d", s + 1));
		eBeamFromAnglesAsPE[s]  = (TH2D*) resultsF->Get(Form("eBeam/eBeamFromAnglesAsPES%d", s + 1));
		eBeamFromAngles[s]  = (TH1D*) resultsF->Get(Form("eBeam/eBeamFromAnglesS%d", s + 1));

		//energy conservation
		deltaPxPx[s] = (TH2D*) resultsF->Get(Form("energyConservation/deltaPxPxS%d", s + 1));
		deltaPyPy[s] = (TH2D*) resultsF->Get(Form("energyConservation/deltaPyPyS%d", s + 1));
		deltaPzPz[s] = (TH2D*) resultsF->Get(Form("energyConservation/deltaPzPzS%d", s + 1));



		//theta elecrton

		deltaThetaEFromPPAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPPAsThetaES%d", s + 1));
		deltaThetaEFromPPAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPPAsThetaPS%d", s + 1));
		deltaThetaEFromPPAsPE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPPAsPES%d", s + 1));
		deltaThetaEFromPPAsPP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPPAsPPS%d", s + 1));
		deltaThetaEFromPP[s] = (TH1D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPPS%d", s + 1));

		deltaThetaEFromEPAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsThetaES%d", s + 1));
		deltaThetaEFromEPAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsThetaPS%d", s + 1));
		deltaThetaEFromEPAsPE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsPES%d", s + 1));
		deltaThetaEFromEPAsPP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsPPS%d", s + 1));
		deltaThetaEFromEP[s] = (TH1D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPS%d", s + 1));


		deltaThetaEFromEPAsThetaEEnergyConservation[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsThetaEEnergyConservationS%d", s + 1));
		deltaThetaEFromEPAsThetaPEnergyConservation[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsThetaPEnergyConservationS%d", s + 1));
		deltaThetaEFromEPAsPEEnergyConservation[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsPEEnergyConservationS%d", s + 1));
		deltaThetaEFromEPAsPPEnergyConservation[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPAsPPEnergyConservationS%d", s + 1));
		deltaThetaEFromEPEnergyConservation[s] = (TH1D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromEPEnergyConservationS%d", s + 1));

		deltaThetaEFromPThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPThetaAsThetaES%d", s + 1));
		deltaThetaEFromPThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPThetaAsThetaPS%d", s + 1));
		deltaThetaEFromPThetaAsPE[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPThetaAsPES%d", s + 1));
		deltaThetaEFromPThetaAsPP[s] = (TH2D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPThetaAsPPS%d", s + 1));
		deltaThetaEFromPTheta[s] = (TH1D*) resultsF->Get(Form("thetaElectron/deltaThetaEFromPThetaS%d", s + 1));

		//theta proton

		deltaThetaPFromPPAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromPPAsThetaES%d", s + 1));
		deltaThetaPFromPPAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromPPAsThetaPS%d", s + 1));
		deltaThetaPFromPPAsPE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromPPAsPES%d", s + 1));
		deltaThetaPFromPPAsPP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromPPAsPPS%d", s + 1));
		deltaThetaPFromPP[s] = (TH1D*) resultsF->Get(Form("thetaProton/deltaThetaPFromPPS%d", s + 1));

		deltaThetaPFromEPAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEPAsThetaES%d", s + 1));
		deltaThetaPFromEPAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEPAsThetaPS%d", s + 1));
		deltaThetaPFromEPAsPE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEPAsPES%d", s + 1));
		deltaThetaPFromEPAsPP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEPAsPPS%d", s + 1));
		deltaThetaPFromEP[s] = (TH1D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEPS%d", s + 1));

		deltaThetaPFromEThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEThetaAsThetaES%d", s + 1));
		deltaThetaPFromEThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEThetaAsThetaPS%d", s + 1));
		deltaThetaPFromEThetaAsPE[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEThetaAsPES%d", s + 1));
		deltaThetaPFromEThetaAsPP[s] = (TH2D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEThetaAsPPS%d", s + 1));
		deltaThetaPFromETheta[s] = (TH1D*) resultsF->Get(Form("thetaProton/deltaThetaPFromEThetaS%d", s + 1));


     	//P electron

		deltaPEFromPThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPThetaAsThetaES%d", s + 1));
		deltaPEFromPThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPThetaAsThetaPS%d", s + 1));
		deltaPEFromPThetaAsPE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPThetaAsPES%d", s + 1));
		deltaPEFromPThetaAsPP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPThetaAsPPS%d", s + 1));
		deltaPEFromPTheta[s] = (TH1D*) resultsF->Get(Form("pElectron/deltaPEFromPThetaS%d", s + 1));

		deltaPEFromPPAsThetaE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPPAsThetaES%d", s + 1));
		deltaPEFromPPAsThetaP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPPAsThetaPS%d", s + 1));
		deltaPEFromPPAsPE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPPAsPES%d", s + 1));
		deltaPEFromPPAsPP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromPPAsPPS%d", s + 1));
		deltaPEFromPP[s] = (TH1D*) resultsF->Get(Form("pElectron/deltaPEFromPPS%d", s + 1));

		deltaPEFromEThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromEThetaAsThetaES%d", s + 1));
		deltaPEFromEThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromEThetaAsThetaPS%d", s + 1));
		deltaPEFromEThetaAsPE[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromEThetaAsPES%d", s + 1));
		deltaPEFromEThetaAsPP[s] = (TH2D*) resultsF->Get(Form("pElectron/deltaPEFromEThetaAsPPS%d", s + 1));
		deltaPEFromETheta[s] = (TH1D*) resultsF->Get(Form("pElectron/deltaPEFromEThetaS%d", s + 1));


		//P proton
		deltaPPFromPThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromPThetaAsThetaES%d", s + 1));
		deltaPPFromPThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromPThetaAsThetaPS%d", s + 1));
		deltaPPFromPThetaAsPE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromPThetaAsPES%d", s + 1));
		deltaPPFromPThetaAsPP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromPThetaAsPPS%d", s + 1));
		deltaPPFromPTheta[s] = (TH1D*) resultsF->Get(Form("pProton/deltaPPFromPThetaS%d", s + 1));

		deltaPPFromEPAsThetaE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEPAsThetaES%d", s + 1));
		deltaPPFromEPAsThetaP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEPAsThetaPS%d", s + 1));
		deltaPPFromEPAsPE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEPAsPES%d", s + 1));
		deltaPPFromEPAsPP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEPAsPPS%d", s + 1));
		deltaPPFromEP[s] = (TH1D*) resultsF->Get(Form("pProton/deltaPPFromEPS%d", s + 1));

		deltaPPFromEThetaAsThetaE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEThetaAsThetaES%d", s + 1));
		deltaPPFromEThetaAsThetaP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEThetaAsThetaPS%d", s + 1));
		deltaPPFromEThetaAsPE[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEThetaAsPES%d", s + 1));
		deltaPPFromEThetaAsPP[s] = (TH2D*) resultsF->Get(Form("pProton/deltaPPFromEThetaAsPPS%d", s + 1));
		deltaPPFromETheta[s] = (TH1D*) resultsF->Get(Form("pProton/deltaPPFromEThetaS%d", s + 1));


	}

	//short version


	std::string parentPDFDirectory = "/w/hallb-scifs17exp/clas12/bclary/CLAS12/validation/validationBrandon/markovMonitor/scripts/monitoring_job_pdf/";

	TCanvas *finalPresentationSixPlots = new TCanvas("finalPresentationSixPlots", "finalPresentationSixPlots", 10, 10, 600, 800);
	finalPresentationSixPlots->Divide(2, 3, 0.001, 0.001);


	TCanvas *finalPresentationOnePlots = new TCanvas("finalPresentationOnePlots", "finalPresentationOnePlots", 10, 10, 600, 800);
	finalPresentationOnePlots->Divide(1, 1, 0.001, 0.001);
//	w 1D 
	for (int s = 0; s < 6; s++){
		finalPresentationSixPlots->cd(s + 1);
		wSectorExclusive[s]->Draw();
		wSectorExclusive[s]->Fit(Form("elasGausPeakS%d", s + 1), "Q", "R", lowWValueFit[s], highWValueFit[s]);
		elasGausPeak[s]->Draw("same");
		if (s == 1) sCaption.DrawLatex(0.832, 0.95, Form("Page 1"));

	}
	finalPresentationSixPlots->SaveAs(Form("%smonitoringJob%d.pdf(",parentPDFDirectory.c_str(), model));

//w vs q2

		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s + 1);
			wQ2SectorExclusive[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.832, 0.95, Form("Page 2"));

		}
		finalPresentationSixPlots->SaveAs(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		finalPresentationOnePlots->cd(1);
		wVsPhi->Draw("COLZ");
		sCaption.DrawLatex(0.832, 0.95, Form("Page 3"));

		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf", parentPDFDirectory.c_str(), model));

//pe 1D from theta e
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaPEFromETheta[s]->Draw("");
			deltaPEFromETheta[s]->Fit(Form("deltaEPFitS%d", s + 1), "Q", "", lowDeltaPEValueFit[s], highDeltaPEValueFit[s]);
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 4"));

		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

//pe as theta E from thetaE
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaPEFromEThetaAsThetaE[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 5"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		finalPresentationOnePlots->cd(1);
		deltaPEFromPThetaAsPhiE->Draw("COLZ");
		sCaption.DrawLatex(0.82, 0.95, Form("Page 6"));
		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));


//p p 1D from theta e
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaPPFromETheta[s]->Draw("");
			deltaPPFromETheta[s]->Fit(Form("deltaPThetaFitS%d", s + 1), "Q", "", lowDeltaPPFromThetaEValueFit[s], highDeltaPPFromThetaEValueFit[s]);
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 7"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf", parentPDFDirectory.c_str(), model));

		//as theta p from thetae
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaPPFromEThetaAsThetaP[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 8"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

//pp as p p from thetae
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaPPFromEThetaAsPP[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 9"));

		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

	//as theta E from thetae

		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaThetaPFromETheta[s]->Draw("");
			deltaThetaPFromETheta[s]->Fit(Form("deltaPPFitS%d", s + 1), "Q", "R", lowDeltaThetaPValueFit[s], highDeltaThetaPValueFit[s]);
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 10"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf", parentPDFDirectory.c_str(), model));


	//as theta E from thetae
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaThetaPFromEThetaAsThetaE[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 11"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf", parentPDFDirectory.c_str(), model));

//as theta p from thetae
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaThetaPFromEThetaAsThetaP[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 12"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));


//vertex quantities
		for (int s = 0; s < 6; s++){
		  finalPresentationSixPlots->cd(s+1);
		  zVsThetaElectron[s]->Draw("COLZ");
		  if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 13"));
                }
                finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		finalPresentationOnePlots->cd(1);
                zVsPhiElectron->Draw("COLZ");
                sCaption.DrawLatex(0.832, 0.95, Form("Page 14"));


		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));
		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			zElectronSector[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 15"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		for (int s = 0; s < 6; s++){
                  finalPresentationSixPlots->cd(s+1);
                  zVsThetaProton[s]->Draw("COLZ");
                  if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 16"));
                }

                finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		
		finalPresentationOnePlots->cd(1);
		zVsPhiProton->Draw("COLZ");
		sCaption.DrawLatex(0.82, 0.95, Form("Page 17"));
		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));
		

		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			zProtonSector[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 18"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf", parentPDFDirectory.c_str(),model));


		finalPresentationOnePlots->cd(1);
		deltaZVsPhiElectron->Draw("COLZ");
		sCaption.DrawLatex(0.82, 0.95, Form("Page 19"));
		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));
		finalPresentationOnePlots->cd(1);
		deltaZVsPhiProton->Draw("COLZ");
		sCaption.DrawLatex(0.82, 0.95, Form("Page 20"));
		finalPresentationOnePlots->Print(Form("%smonitoringJob%d.pdf",parentPDFDirectory.c_str(), model));

		for (int s = 0; s < 6; s++){
			finalPresentationSixPlots->cd(s+1);
			deltaZElectronProtonSector[s]->Draw("COLZ");
			if (s == 1) sCaption.DrawLatex(0.82, 0.95, Form("Page 21"));
		}
		finalPresentationSixPlots->Print(Form("%smonitoringJob%d.pdf)",parentPDFDirectory.c_str(), model));


}		




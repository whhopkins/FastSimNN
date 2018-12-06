// ANL analysis tutorial. S.Chekanov (ANL)

#ifndef Ana_h
#define Ana_h

#include<vector>
#include<iostream>
#include"TROOT.h"
#include<TSystem.h>
#include"stdio.h"
#include"TH1D.h"
#include"TH2D.h"
#include<map>
#include"TFile.h"
#include"TTree.h"
#include"TProfile.h"
#include<TClonesArray.h>
#include"TRandom3.h"
#include"TLorentzVector.h"
#include <string>
#include"SystemOfUnits.h"
#include "TMath.h"
#include "TRandom.h"
#include<fstream>
#include<stdlib.h>
#include "LParticle.h"

// Local include(s):
#include "fann.h"
#include "parallel_fann.h"
#include <libconfig.h++>

using namespace libconfig; 
using namespace std;


// main analysis class. inherent analysis.h which should be rebuild each time
class Ana  {
    public:
          virtual ~Ana();
          Ana();
          int nevv; 
	  int  Init();
          int  Finish();
          int  AnalysisJets(vector<LParticle> JetsTrue, vector<LParticle> JetsReco);

   vector <string> ntup;       // ntuple list
   TH1D *h_debug ;             // global debug file
   int   MaxEvents;            // max number of events


   int  MuPileup; // number of pileup (mu)
   static const int num_threads = 16; // number of threads 
   static const int nBins=16;
   static const int nBatch=100000;  // number of events in batches for training
   static const int nEpoch=200;    // max number of epochs 
   static const int nBinsNN=400;   // number of bins for resolution plots
   static const int MinEntries=10; // min nr of entries in pT for NN training (per bunch);

   double DeltaR;  // parameter used to match true jets with reco
   double MSESTOP;  // when stop training..

   // NN structure 
   static const int num_input=4;
   static const int num_output=4*nBinsNN;


   bool firstTime[nBins-1]; // if false, continue training;
   // if ANN is found, we will read the old one. See src/Init.cxx
   double eBins[nBins];
   // for correction
   vector<vector<float>> finput_jets;
   vector<vector<float>> foutput_jets;
   // for efficiency
   vector<vector<float>> finput_jets_eff;
   vector<vector<int>>   foutput_jets_eff;

   
   double initialRME[nBins-1];
   double finalRME[nBins-1];
   int    eventsBins[nBins-1]; 
   struct fann *ann_jets[nBins-1];
   string ann_jets_name[nBins-1];
   struct fann *ann_jets_eff[nBins-1];
   string ann_jets_eff_name[nBins-1];
   // to deal with scale extensions
   float jet_escale;
   float jet_eshift;
   float jet_etascale;
   float jet_etashift;

   static const int num_layers = 3;
   //number in hidden layer 1
   static  const int num_neurons_hidden_1=40;
   TH1D *h_jetpt;
   TH1D *h_jetpt_truth;

   double minPT;
   double maxEta;
 
 
protected:

   TH1D *h_dR;

   TH1D *h_in1;
   TH1D *h_in2;
   TH1D *h_in3;
   TH1D *h_in4;

   TH1D *h_out1;
   TH1D *h_out2;
   TH1D *h_out3;
   TH1D *h_out4;
   TH1D *h_out5;

   string ffile; 
   TFile *RootFile;
   TTree * m_ntuple;

   // objects for output ntuple 
   int RunNumber; // run number
   int EventNumber; // event number
   TClonesArray *a_jets;

   double mcEventWeight;

   bool isMC;
   int systematics;
   int type;
   double weight;

};

#endif

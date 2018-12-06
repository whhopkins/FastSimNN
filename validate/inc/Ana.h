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
#include <stdlib.h>


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
          int nev; 
	  int  Init();
          int  Finish();
          int  AnalysisJets(vector<LParticle> JetsTrue, vector<LParticle> JetsReco);

   vector <string> ntup;       // ntuple list
   TH1D *h_debug ;             // global debug file
   int   MaxEvents;            // max number of events


 
  // NN structure 
   unsigned int num_layers = 3;
   unsigned int num_neurons_hidden_1=40;
   int nBinsNN;   // number of bins for resolution plots

    // NN structure 
   unsigned int num_input;
   unsigned int num_output;

   int  MuPileup; 
   int nBins;
   double DeltaR;  // parameter used to match true jets with reco
   double  *eBins;
   double *initialRME;
   double *finalRME;
   // jets
   struct fann **ann_jets; 
   struct fann **ann_jets_eff;
   string *ann_jets_name;
   string *ann_jets_eff_name;

   
   int *BinOverTrue; // bins for reco/true distributions 

   TH1D *h_jetpt;
   TH1D *h_jetpt_truth;
   TH1D *h_jetpt_nn;
   TH1D *h_jetcutflow;

   TH1D *h_ptcor;
   TH1D *h_etacor; 
   TH1D *h_phicor;

   double minPT;
   double maxEta;
 
protected:

   string ffile; 
   TFile *RootFile;

  TH1D *h_jetres100;
  TH1D *h_jetres1000;


   TH1D *h_in1;
   TH1D *h_in2;
   TH1D *h_in3;
   TH1D *h_in4;

   TH1D *h_out1;
   TH1D *h_out2;
   TH1D *h_out3;
   TH1D *h_out4;
   TH1D *h_out5_eff;

   //TTree*  m_jets; //!
   std::vector<Double32_t> m_jetpt; //!
   std::vector<Double32_t> m_jeteta; //!
   std::vector<Double32_t> m_jetphi; //!

   //TTree*  m_gjets; //!
   std::vector<Double32_t> m_gjetpt; //!
   std::vector<Double32_t> m_gjeteta; //!
   std::vector<Double32_t> m_gjetphi; //!

   //TTree*  m_nnjets; //!
   std::vector<Double32_t> m_nnjetpt; //!
   std::vector<Double32_t> m_nnjeteta; //!
   std::vector<Double32_t> m_nnjetphi; //!

   TTree* m_ntuple;

   // objects for output ntuple 
   int RunNumber; // run number
   int EventNumber; // event number
/*
   TClonesArray *a_jets;
   TClonesArray *t_jets;
   TClonesArray *f_jets;
*/

   double mcEventWeight;

 

   vector<LParticle> JetsReco; 
   vector<LParticle> JetsTrue;
   vector<LParticle> JetsFastNN;

   bool isMC;
   int systematics;
   int type;
   double weight;

};

#endif

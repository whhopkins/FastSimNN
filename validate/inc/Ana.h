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
   int nBinsNN;   // number of bins for resolution plots

    // NN structure 
   unsigned int num_input;
   unsigned int num_output;
   unsigned int num_input_eff;
   unsigned int num_output_eff;

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
   // to deal with scale extensions
   float jet_escale;
   float jet_eshift;
   float jet_etascale;
   float jet_etashift;
   float jet_mscale;
   float jet_mshift;

   
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
   TH1D *h_out6_btag;

   //random bins generated
   TH1D *h_rout1;
   TH1D *h_rout2;
   TH1D *h_rout3;
   TH1D *h_rout4;



   std::vector<Double32_t> m_jetpt; //!
   std::vector<Double32_t> m_jeteta; //!
   std::vector<Double32_t> m_jetphi; //!
   std::vector<Double32_t> m_jetm; //!
   std::vector<Int_t>      m_jetbtag; //!

   std::vector<Double32_t> m_gjetpt; //!
   std::vector<Double32_t> m_gjeteta; //!
   std::vector<Double32_t> m_gjetphi; //!
   std::vector<Double32_t> m_gjetm; //!
   std::vector<Int_t>      m_gjetbtag; //!

   std::vector<Double32_t> m_nnjetpt; //!
   std::vector<Double32_t> m_nnjeteta; //!
   std::vector<Double32_t> m_nnjetphi; //!
   std::vector<Double32_t> m_nnjetm; //!
   std::vector<Double32_t> m_nnjetbtag; //!

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

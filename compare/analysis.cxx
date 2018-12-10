// -----------------------------------------------
// An analysis program that shows how to read ProMC files from HepSim
// and perform an analysis of reconstructed objects from Delphes 3.X on the fly
// S.Chekanov (ANL) chekanov@anl.gov
/// -----------------------------------------------

#include <stdexcept>
#include <iostream>
#include<fstream>
#include <sstream>
#include <memory>
#include <deque>
#include <stdlib.h>
#include <signal.h>
#include <stdio.h>
#include <algorithm>
#include "TROOT.h"
#include "TApplication.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"
#include "LParticle.h"


using namespace std;
#include <vector>
#include <map>
#include <TClonesArray.h>
#include <TLorentzVector.h>


//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{

      int MaxEvents=200000;

// Choose input and output files from the command line
        string infile("-"), outfile("-");
        if (argc == 3) {
                infile = argv[1];
                outfile = argv[2];
        } else if (argc != 3) {
                cerr << "Usage: " << argv[0] << "[<input> <out>]" << endl;
                exit(1);
        }


   TFile *f = new TFile( infile.c_str());
   TTree *m_ntuple = (TTree*)f->Get("Ntuple");

   std::vector<Double32_t> *jetpt=0; //!
   std::vector<Double32_t> *jeteta=0; //!
   std::vector<Double32_t> *jetphi=0; //!
   std::vector<Double32_t> *jetm=0; //!
   std::vector<Int_t>      *jetbtag=0; //!

   std::vector<Double32_t> *gjetpt=0; //!
   std::vector<Double32_t> *gjeteta=0; //!
   std::vector<Double32_t> *gjetphi=0; //!
   std::vector<Double32_t> *gjetm=0; //!
   std::vector<Int_t>      *gjetbtag=0; //!

   std::vector<Double32_t> *nnjetpt=0; //!
   std::vector<Double32_t> *nnjeteta=0; //!
   std::vector<Double32_t> *nnjetphi=0; //!
   std::vector<Double32_t> *nnjetm=0; //!
   std::vector<Double32_t> *nnjetbtag=0; //!

   TBranch        *b_jetpt;
   TBranch        *b_jeteta;
   TBranch        *b_jetphi;
   TBranch        *b_jetm;
   TBranch        *b_jetbtag;

   TBranch        *b_gjetpt;
   TBranch        *b_gjeteta;
   TBranch        *b_gjetphi;
   TBranch        *b_gjetm;
   TBranch        *b_gjetbtag;

   TBranch        *b_nnjetpt;
   TBranch        *b_nnjeteta;
   TBranch        *b_nnjetphi;
   TBranch        *b_nnjetm;
   TBranch        *b_nnjetbtag;

   m_ntuple->SetBranchAddress("AntiKt4JetPt",   &jetpt, &b_jetpt);
   m_ntuple->SetBranchAddress("AntiKt4JetEta",  &jeteta, &b_jeteta);
   m_ntuple->SetBranchAddress("AntiKt4JetPhi",   &jetphi, &b_jetphi);
   m_ntuple->SetBranchAddress("AntiKt4JetM",     &jetm, &b_jetm);
   m_ntuple->SetBranchAddress("AntiKt4JetBtag",  &jetbtag, &b_jetbtag);

   m_ntuple->SetBranchAddress("AntiKt4TruthJetPt",   &gjetpt, &b_gjetpt);
   m_ntuple->SetBranchAddress("AntiKt4TruthJetEta",  &gjeteta, &b_gjeteta);
   m_ntuple->SetBranchAddress("AntiKt4TruthJetPhi",   &gjetphi, &b_gjetphi);
   m_ntuple->SetBranchAddress("AntiKt4TruthJetM",      &gjetm, &b_gjetm);
   m_ntuple->SetBranchAddress("AntiKt4TruthJetBtag",   &gjetbtag, &b_gjetbtag);

   m_ntuple->SetBranchAddress("AntiKt4NNJetPt",   &nnjetpt, &b_nnjetpt);
   m_ntuple->SetBranchAddress("AntiKt4NNJetEta",  &nnjeteta, &b_nnjeteta);
   m_ntuple->SetBranchAddress("AntiKt4NNJetPhi",   &nnjetphi, &b_nnjetphi);
   m_ntuple->SetBranchAddress("AntiKt4NNJetM",      &nnjetm, &b_nnjetm);
   m_ntuple->SetBranchAddress("AntiKt4NNJetBtag",   &nnjetbtag, &b_nnjetbtag);


 TFile * RootFile = new TFile(outfile.c_str(), "RECREATE", "ProMC record");
 if (!RootFile){
             std::cout << "Error: Cannot create ROOT file" << std::endl;
             exit(1);
        }



  Int_t nentries = (Int_t)m_ntuple->GetEntries();
   for (Int_t event=0; event<nentries; event++) {
     m_ntuple->GetEntry(event);
     if (event%1000==0 || event<100) cout << "Process=" <<  event << endl;

     cout << "jetpt=" << jetpt->size() << endl;


  }


        f->Close();


        RootFile->Write();
        RootFile->Print();
        RootFile->Close();

        cout << "Output ROOT file: " << outfile << endl;


   return 0;


}

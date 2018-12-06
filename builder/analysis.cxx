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

#include "modules/Delphes.h"
#include "classes/DelphesStream.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesModule.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootTask.h"
#include "ExRootAnalysis/ExRootConfReader.h"

#include "Ana.h"

#include "fann.h"
#include "parallel_fann.h"


using namespace std;
extern vector<string> getfiles(char * directory);

bool reorder(const TLorentzVector &a, const TLorentzVector &b)
{
	return a.Pt() > b.Pt();
}

double const twopi=6.28318530718;
double const pi=3.14159265359;

//---------------------------------------------------------------------------

static bool interrupted = false;
void SignalHandler(int sig)
{
	interrupted = true;
}

//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{
	char appName[] = "analysis";
	stringstream message;


        gSystem->Load("delphes/libDelphes");


        Ana ana; // initialaze class with calculations


// decide how many events to run
  if( argc > 1 ) {
      const Long64_t e = atoll( argv[ 1 ] );
         ana.MaxEvents=e;
      }





        ana.Init(); // initialize histogram and global variable
        const int Nfiles = ana.ntup.size();
        cout << " -> No of files to read:" << Nfiles << endl;


        bool stop=false;
        int nev=0;
	for(unsigned int m=0; m < Nfiles; m++){
                if (stop) break;

		string Rfile=ana.ntup[m];
		TChain chain("Delphes");
                cout << "reading file=" << Rfile << endl;

		chain.Add( Rfile.c_str() );
		ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);

		Long64_t numberOfEntries = treeReader->GetEntries();

		// Get pointers to branches used in this analysis
		TClonesArray *branchJet = treeReader->UseBranch("Jet");
		TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
		TClonesArray *branchElectron = treeReader->UseBranch("Electron");

		// Loop over all events
		for(Int_t entry = 0; entry < numberOfEntries; ++entry)
		{

                        nev++;
                        if (nev%100==0 || nev<100) cout << "Event=" << nev << endl;

                        if (ana.MaxEvents !=-1 && nev>ana.MaxEvents) { stop=true; break; }; 


			// Load selected branches with data from specified event
			treeReader->ReadEntry(entry);

                       vector<LParticle> JetsReco;
                       vector<LParticle> JetsTrue;

                     
                        for(int i = 0; i < branchJet->GetEntriesFast(); ++i) { 
				// Take first jet
				Jet *jet = (Jet*) branchJet->At(i);
                                if (abs(jet->Eta)>ana.maxEta) continue;
 
                                TLorentzVector l;
                                l.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
                                LParticle p;
                                p.SetP(l);
                                ana.h_jetpt->Fill(jet->PT);
                                JetsReco.push_back(p);
                                //cout << jet->PT << endl;
			}

                        for(int i = 0; i < branchGenJet->GetEntriesFast(); ++i) { 
                                // Take first jet
                                Jet *jet = (Jet*) branchGenJet->At(i);
                                if (abs(jet->Eta)>ana.maxEta) continue;
                                TLorentzVector l;
                                l.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
                                LParticle p;
                                p.SetP(l);
                                ana.h_jetpt_truth->Fill(jet->PT);
                                JetsTrue.push_back(p);
                        }



          // main event loop
          ana.AnalysisJets(JetsTrue,JetsReco);  // analysis at the end












		}



	}// end loop over files


        ana.Finish();

	return 0;


}

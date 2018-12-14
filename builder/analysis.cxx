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

const double PI   = TMath::Pi();
const double PI2  = 2*PI;

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
	//char appName[] = "analysis";
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
	for(int m=0; m < Nfiles; m++){
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
		TClonesArray *branchMuon = treeReader->UseBranch("Muon");
		TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
                TClonesArray *branchParticle = treeReader->UseBranch("Particle");

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
				if (abs(jet->Eta)>ana.maxEta)     continue; // use large Eta (needs for matching)
                                if (jet->Mass<1)                  continue; // cannot be zero 

				// overlap removal
				bool remove=false;
				for(int i1 = 0; i1 < branchMuon->GetEntriesFast(); ++i1) {
					Muon *muon = (Muon*) branchMuon->At(i1);
					double deta=jet->Eta  - muon->Eta;
					double dphi=jet->Phi  - muon->Phi;
					if (abs(dphi)>PI) dphi=PI2-abs(dphi);
					double dR=sqrt(deta*deta + dphi*dphi);
					if (dR<0.1) remove=true;
				}
				if (remove) continue;


				// overlap removal
				remove=false;
				for(int i1 = 0; i1 < branchElectron->GetEntriesFast(); ++i1) {
					Electron *ele = (Electron*) branchElectron->At(i1);
					double deta=jet->Eta  - ele->Eta;
					double dphi=jet->Phi  - ele->Phi;
					if (abs(dphi)>PI) dphi=PI2-abs(dphi);
					double dR=sqrt(deta*deta + dphi*dphi);
					if (dR<0.1) remove=true;
				}
				if (remove) continue;


				// overlap removal
				remove=false;
				for(int i1 = 0; i1 < branchPhoton->GetEntriesFast(); ++i1) {
					Photon *ph = (Photon*) branchPhoton->At(i1);
					double deta=jet->Eta  - ph->Eta;
					double dphi=jet->Phi  - ph->Phi;
					if (abs(dphi)>PI) dphi=PI2-abs(dphi);
					double dR=sqrt(deta*deta + dphi*dphi);
					if (dR<0.1) remove=true;
				}
				if (remove) continue;


				TLorentzVector l;
				l.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				LParticle p;
                                p.SetType(0);
                                Bool_t BtagOk = ( jet->BTag & (1 << i) );
                                if (BtagOk) p.SetType(1);
                                //if (p.GetType() >0) cout << "reco B-tagging! " << endl;

				p.SetP(l);
				ana.h_jetpt->Fill(jet->PT);
				JetsReco.push_back(p);
				//cout << jet->PT << endl;
			}

			// use proper kinematic cuts
			for(int i = 0; i < branchGenJet->GetEntriesFast(); ++i) {
				// Take first jet
				Jet *jet = (Jet*) branchGenJet->At(i);
				if (abs(jet->Eta)>ana.maxEta) continue;
				if (jet->PT      <ana.minPT)  continue;
                                if (jet->Mass<1)              continue; // cannot be zero 


                                // get fraction of b-quark inside this jet (in %) 
                                double btag_fracmom=0;
                                for(int i1 = 0; i1 < branchParticle->GetEntriesFast(); ++i1) {
                                        GenParticle *ph = (GenParticle*) branchParticle->At(i1);
                                        int    pdgCode = TMath::Abs(ph->PID);
                                        if (pdgCode != 5) continue;
                                        double deta=jet->Eta  - ph->Eta;
                                        double dphi=jet->Phi  - ph->Phi;
                                        if (abs(dphi)>PI) dphi=PI2-abs(dphi);
                                        double dR=sqrt(deta*deta + dphi*dphi);
                                        double rat=(ph->PT / jet->PT); 
                                        if (dR<0.4 && rat>0.25) btag_fracmom= 100 * rat;
                                }

                                //if (btag_fracmom>0) cout << btag_fracmom << endl;

// Loop over all jet's constituents
/*
      int btag=0;
      for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        TObject *object = jet->Constituents.At(j);
        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          GenParticle * particle = (GenParticle*) object;
          int status = particle->Status;
          int pdgCode = TMath::Abs(particle->PID);

          cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << endl;
          //momentum += particle->P4();
         }
         }

*/

				TLorentzVector l;
				l.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
				LParticle p;
                                p.SetType((int)(btag_fracmom));
                                //Bool_t BtagOk = ( jet->BTag & (1 << i) );
                                //if (BtagOk) p.SetType(1);
                                //if (p.GetType() >0) cout << "truth B-tagging! " << endl;
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

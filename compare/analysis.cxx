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

const double PI   = TMath::Pi();
const double PI2  = 2*PI;


using namespace std;
#include <vector>
#include <map>
#include <TClonesArray.h>
#include <TLorentzVector.h>


//---------------------------------------------------------------------------

int main(int argc, char *argv[])
{

	int MaxEvents=200000;
	double EtaMax=2.5;
	double PtMin=25;

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
	std::vector<Int_t>      *nnjetbtag=0; //!

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

	std::vector<Double32_t> *gmupt=0; //!
	std::vector<Double32_t> *gmueta=0; //!
	std::vector<Double32_t> *gmuphi=0; //!
	std::vector<Double32_t> *mupt=0; //!
	std::vector<Double32_t> *mueta=0; //!
	std::vector<Double32_t> *muphi=0; //!
	std::vector<Double32_t> *nnmupt=0; //!
	std::vector<Double32_t> *nnmueta=0; //!
	std::vector<Double32_t> *nnmuphi=0; //!
	TBranch        *b_mupt;
	TBranch        *b_mueta;
	TBranch        *b_muphi;
	TBranch        *b_gmupt;
	TBranch        *b_gmueta;
	TBranch        *b_gmuphi;
	TBranch        *b_nnmupt;
	TBranch        *b_nnmueta;
	TBranch        *b_nnmuphi;

	std::vector<Double32_t> *gphpt=0; //!
	std::vector<Double32_t> *gpheta=0; //!
	std::vector<Double32_t> *gphphi=0; //!
	std::vector<Double32_t> *phpt=0; //!
	std::vector<Double32_t> *pheta=0; //!
	std::vector<Double32_t> *phphi=0; //!
	std::vector<Double32_t> *nnphpt=0; //!
	std::vector<Double32_t> *nnpheta=0; //!
	std::vector<Double32_t> *nnphphi=0; //!
	TBranch        *b_phpt;
	TBranch        *b_pheta;
	TBranch        *b_phphi;
	TBranch        *b_gphpt;
	TBranch        *b_gpheta;
	TBranch        *b_gphphi;
	TBranch        *b_nnphpt;
	TBranch        *b_nnpheta;
	TBranch        *b_nnphphi;



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


	// muons
	m_ntuple->SetBranchAddress("muonPt",   &mupt, &b_mupt);
	m_ntuple->SetBranchAddress("muonEta",  &mueta, &b_mueta);
	m_ntuple->SetBranchAddress("muonPhi",   &muphi, &b_muphi);

	m_ntuple->SetBranchAddress("muonTruthPt",   &gmupt, &b_gmupt);
	m_ntuple->SetBranchAddress("muonTruthEta",  &gmueta, &b_gmueta);
	m_ntuple->SetBranchAddress("muonTruthPhi",   &gmuphi, &b_gmuphi);

	m_ntuple->SetBranchAddress("muonNNPt",   &nnmupt, &b_nnmupt);
	m_ntuple->SetBranchAddress("muonNNEta",  &nnmueta, &b_nnmueta);
	m_ntuple->SetBranchAddress("muonNNPhi",   &nnmuphi, &b_nnmuphi);

	// photons
	m_ntuple->SetBranchAddress("photonPt",   &phpt, &b_phpt);
	m_ntuple->SetBranchAddress("photonEta",  &pheta, &b_pheta);
	m_ntuple->SetBranchAddress("photonPhi",   &phphi, &b_phphi);

	m_ntuple->SetBranchAddress("photonTruthPt",   &gphpt, &b_gphpt);
	m_ntuple->SetBranchAddress("photonTruthEta",  &gpheta, &b_gpheta);
	m_ntuple->SetBranchAddress("photonTruthPhi",   &gphphi, &b_gphphi);

	m_ntuple->SetBranchAddress("photonNNPt",   &nnphpt, &b_nnphpt);
	m_ntuple->SetBranchAddress("photonNNEta",  &nnpheta, &b_nnpheta);
	m_ntuple->SetBranchAddress("photonNNPhi",   &nnphphi, &b_nnphphi);



	TFile * RootFile = new TFile(outfile.c_str(), "RECREATE", "ProMC record");
	if (!RootFile){
		std::cout << "Error: Cannot create ROOT file" << std::endl;
		exit(1);
	}

	const double DeltaR=0.2;
	TH1D *h_dR = new TH1D("jet_dR", "truth-jet distance", 500, 0, 7);

	// for b-tagging efficiency
	TH1D *h_jet1_btag = new TH1D("jet1_btag", "jet1_btag",20,0.0,500.);
	TH1D *h_jet2_btag =  new TH1D("jet2_btag", "jet2_btag",20,0.0,500.);
	TH1D *h_jet1_all = new TH1D("jet1_all", "jet1_all",20,0.0,500.);
	TH1D *h_jet2_all = new TH1D("jet2_all", "jet2_all",20,0.0,500.);
	h_jet1_btag->Sumw2();
	h_jet2_btag->Sumw2();
	h_jet1_all->Sumw2();
	h_jet2_all->Sumw2();

	double ptBins[] = {PtMin,40,60,80,100,140,180, 210, 240, 290, 340, 400, 500, 800, 1000, 1500,  2000, 2500};
	const int nBins=sizeof(ptBins)/sizeof(double);

	// for muon efficiency
	TH1D *h_mu_reco_pt = new TH1D("mu_reco_pt", "mu_pt for delphes",nBins-1, ptBins);
	TH1D *h_mu_nn_pt =  new TH1D("mu_nn_pt", "mu_pt for NN",nBins-1, ptBins);
	TH1D *h_mu_all_pt = new TH1D("mu_true_pt", "mu pT for truth",nBins-1, ptBins);
	h_mu_reco_pt->Sumw2();
	h_mu_nn_pt->Sumw2();
	h_mu_all_pt->Sumw2();

	// for photon efficiency
	TH1D *h_ph_reco_pt = new TH1D("ph_reco_pt", "ph_pt for delphes", nBins-1, ptBins);
	TH1D *h_ph_nn_pt =  new TH1D("ph_nn_pt", "ph_pt for NN",nBins-1, ptBins);
	TH1D *h_ph_all_pt = new TH1D("ph_true_pt", "ph pT for truth",nBins-1, ptBins);
	h_ph_reco_pt->Sumw2();
	h_ph_nn_pt->Sumw2();
	h_ph_all_pt->Sumw2();


         // jet mass
        TH1D *h_jet_reco_mass = new TH1D("jet_reco_mass", "mass for delphes", 50,0,800);
        TH1D *h_jet_nn_mass =  new TH1D("jet_nn_mass", "mass for NN",20,0,400);
        TH1D *h_jet_truth1_mass = new TH1D("jet_true1_mass", "mass for truth",50, 0, 800);
        TH1D *h_jet_truth2_mass = new TH1D("jet_true2_mass", "mass for truth",50, 0, 800);
 

	// jet pT resolution
	static int nmax_jet=22;
	TH1D *h_jet1_res[nmax_jet];
	TH1D *h_jet1_ptr[nmax_jet];
	TH1D *h_jet2_res[nmax_jet];
	TH1D *h_jet2_ptr[nmax_jet];

	// muon pT resolution
	TH1D *h_muon1_res[nmax_jet];
	TH1D *h_muon1_ptr[nmax_jet];
	TH1D *h_muon2_res[nmax_jet];
	TH1D *h_muon2_ptr[nmax_jet];

	for (int j=0; j<nmax_jet; j++){
		int nbins=50+40*j; // increasing number of bins
		// Delphes
		h_jet1_res[j] = new TH1D(Form("jet1_resolution_%02d",j), Form("jets1_res_%02d",j),nbins,0.0,4.0);
		h_jet1_ptr[j] = new TH1D(Form("jet1_resolution_pt_%02d",j), Form("jets1_res_pt_%02d",j),2500,0,5000);
		h_jet1_res[j]->Sumw2();
		h_jet1_ptr[j]->Sumw2();
		// NN
		h_jet2_res[j] = new TH1D(Form("jet2_resolution_%02d",j), Form("jets2_res_%02d",j),nbins,0.0,4.0);
		h_jet2_ptr[j] = new TH1D(Form("jet2_resolution_pt_%02d",j), Form("jets2_res_pt_%02d",j),2500,0,5000);
		h_jet2_res[j]->Sumw2();
		h_jet2_ptr[j]->Sumw2();

                // muons //
		// Delphes
		h_muon1_res[j] = new TH1D(Form("muon1_resolution_%02d",j), Form("muon1_res_%02d",j),nbins,0.6,1.4);
		h_muon1_ptr[j] = new TH1D(Form("muon1_resolution_pt_%02d",j), Form("muon1_res_pt_%02d",j),2500,0,5000);
		h_muon1_res[j]->Sumw2();
		h_muon1_ptr[j]->Sumw2();
		// NN
		h_muon2_res[j] = new TH1D(Form("muon2_resolution_%02d",j), Form("muon2_res_%02d",j),nbins,0.6,1.4);
		h_muon2_ptr[j] = new TH1D(Form("muon2_resolution_pt_%02d",j), Form("muon2_res_pt_%02d",j),2500,0,5000);
		h_muon2_res[j]->Sumw2();
		h_muon2_ptr[j]->Sumw2();


		double x1=5+pow(2,(0.35*(j+12)));
		double x2=10+pow(2,(0.35*(j+12+1)));
		cout << j << ") " << x1 << " - " << x2 << " GeV " << endl;


	}

	// jet Eta resolution
	TH1D *h_jeteta1_res[nmax_jet];
	TH1D *h_jeteta1_ptr[nmax_jet];
	TH1D *h_jeteta2_res[nmax_jet];
	TH1D *h_jeteta2_ptr[nmax_jet];

	int nbins=100;
	for (int j=0; j<nmax_jet; j++){
		// Delphes
		h_jeteta1_res[j] = new TH1D(Form("jeteta1_resolution_%02d",j), Form("jetseta1_res_%02d",j),nbins,0.5,1.5);
		h_jeteta1_ptr[j] = new TH1D(Form("jeteta1_resolution_eta_%02d",j), Form("jetseta1_res_eta_%02d",j),200,-3,3);
		h_jeteta1_res[j]->Sumw2();
		h_jeteta1_ptr[j]->Sumw2();
		// NN
		h_jeteta2_res[j] = new TH1D(Form("jeteta2_resolution_%02d",j), Form("jetseta2_res_%02d",j),nbins,0.5,1.5);
		h_jeteta2_ptr[j] = new TH1D(Form("jeteta2_resolution_eta_%02d",j), Form("jetseta2_res_eta_%02d",j),200,-3,3);
		h_jeteta2_res[j]->Sumw2();
		h_jeteta2_ptr[j]->Sumw2();
	}




	Int_t nentries = (Int_t)m_ntuple->GetEntries();
	for (Int_t event=0; event<nentries; event++) {
		m_ntuple->GetEntry(event);
		if (event%1000==0 || event<100) cout << "Process=" <<  event << endl;


		for(unsigned int j = 0; j<gjetpt->size(); j++){
			double phiT = gjetphi->at(j);
			double ptT =  gjetpt->at(j);
			double etaT = gjeteta->at(j);
			double massT =  gjetm->at(j);
			float  btagT=  gjetbtag->at(j); // get  b-quark in 100%

			if (abs(etaT)>EtaMax) continue;
			if (ptT<PtMin) continue;

			// reco jets
			double pt_matched1 =-1;
			double eta_matched1 =-1;
			double btag_matched1=0;
			for(unsigned int i = 0; i<jetpt->size(); i++){
				double phi = jetphi->at(i);
				double pt =  jetpt->at(i);
				double eta = jeteta->at(i);
				double mass =jetm->at(i);
				float  btag= jetbtag->at(i); // get  b-quark in 100%
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				h_dR->Fill(dR);
			        if (dR<DeltaR) {pt_matched1=pt; 
                                                eta_matched1=eta; 
                                                btag_matched1=btag;
                                                h_jet_reco_mass->Fill(mass); 
                                                h_jet_truth1_mass->Fill(massT); 
                                                }

			}

			// NN jets
			double pt_matched2 =-1;
			double eta_matched2 =-1;
			double btag_matched2=0;
			for(unsigned int i = 0; i<nnjetpt->size(); i++){
				double phi = nnjetphi->at(i);
				double pt =  nnjetpt->at(i);
				double eta = nnjeteta->at(i);
				double mass =nnjetm->at(i);
				float  btag= nnjetbtag->at(i); // get  b-quark in 100%
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				h_dR->Fill(dR);
			       if (dR<DeltaR) {pt_matched2=pt; 
                                               eta_matched2=eta;  
                                               btag_matched2=btag; 
                                               h_jet_nn_mass->Fill(mass); 
                                               h_jet_truth2_mass->Fill(massT);
                                             }
				//if (abs(dPhi)<0.15) {eta_matched2=eta;}

			}

			// step for Eta
			const double delta=2*EtaMax/(nmax_jet-1);

			// Delphes jets
			if (pt_matched1>0 && ptT>0) {

				h_jet1_all->Fill(ptT);
				if (btag_matched1 > 0) h_jet1_btag->Fill(ptT);

				for (int kk=0; kk<nmax_jet; kk++){
					double x1=5+pow(2,(0.35*(kk+12)));
					double x2=10+pow(2,(0.35*(kk+12+1)));
					//cout << kk << " " << x1 << " " << x2 << endl;
					if (ptT>x1 && ptT<x2) {
						h_jet1_res[kk]->Fill(pt_matched1/ptT);
						h_jet1_ptr[kk]->Fill(ptT); }
				}

				for (int kk=0; kk<nmax_jet-1; kk++){
					double x1=-1*EtaMax+kk*delta;
					double x2= x1+delta;
					//cout << kk << " " << x1 << " " << x2 << endl;
					if (etaT>x1 && etaT<x2 && etaT !=0) {
						h_jeteta1_res[kk]->Fill(abs(eta_matched1/etaT));
						h_jeteta1_ptr[kk]->Fill(etaT);
					}
				}

			}



			// NN jets
			if (pt_matched2>0 && ptT>0) {


				h_jet2_all->Fill(ptT);
				if (btag_matched2 > 0) h_jet2_btag->Fill(ptT);

				for (int kk=0; kk<nmax_jet; kk++){
					double x1=5+pow(2,(0.35*(kk+12)));
					double x2=10+pow(2,(0.35*(kk+12+1)));
					//cout << kk << " " << x1 << " " << x2 << endl;
					if (ptT>x1 && ptT<x2)
					{h_jet2_res[kk]->Fill(pt_matched2/ptT); h_jet2_ptr[kk]->Fill(ptT); }
				}

				for (int kk=0; kk<nmax_jet-1; kk++){
					double x1=-1*EtaMax+kk*delta;
					double x2= x1+delta;
					//cout << kk << " " << x1 << " " << x2 << endl;
					if (etaT>x1 && etaT<x2 && etaT !=0)  {
						h_jeteta2_res[kk]->Fill(abs(eta_matched2/etaT));
						h_jeteta2_ptr[kk]->Fill(etaT); }
				}



			}

		} // end loop over true jets




		// *********************** muon's efficiency **************** /
		for(unsigned int i1 = 0; i1<gmupt->size(); i1++){
			double phiT = gmuphi->at(i1);
			double ptT =   gmupt->at(i1);
			double etaT = gmueta->at(i1);
			if (abs(etaT)>EtaMax) continue;
			if (ptT<PtMin) continue;

			h_mu_all_pt->Fill(ptT);

			double pt_matched1=-1;
			for(unsigned int i2 = 0; i2<mupt->size(); i2++){
				double phi = muphi->at(i2);
				double pt =  mupt->at(i2);
				double eta = mueta->at(i2);
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				if (dR<0.15) pt_matched1=pt;
			}


			if (pt_matched1>0){
				h_mu_reco_pt->Fill(ptT);
				for (int kk=0; kk<nmax_jet; kk++){
					double x1=5+pow(2,(0.35*(kk+12)));
					double x2=10+pow(2,(0.35*(kk+12+1)));
					//cout << kk << " " << x1 << " " << x2 << endl;
					if (ptT>x1 && ptT<x2) {
						h_muon1_res[kk]->Fill(pt_matched1/ptT);
						h_muon1_ptr[kk]->Fill(ptT); }
				}


			}


			double pt_matched2=-1;
			for(unsigned int i2 = 0; i2<nnmupt->size(); i2++){
				double phi = nnmuphi->at(i2);
				double pt =  nnmupt->at(i2);
				double eta = nnmueta->at(i2);
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				if (dR<0.15) pt_matched2=pt;
			}


			if (pt_matched2>0){
				h_mu_nn_pt->Fill(ptT);
				for (int kk=0; kk<nmax_jet; kk++){
					double x1=5+pow(2,(0.35*(kk+12)));
					double x2=10+pow(2,(0.35*(kk+12+1)));
					if (ptT>x1 && ptT<x2) {
						h_muon2_res[kk]->Fill(pt_matched2/ptT);
						h_muon2_ptr[kk]->Fill(ptT); }
				}
			}



		}

		// photon efficiency
		for(unsigned int i1 = 0; i1<gphpt->size(); i1++){
			double phiT = gphphi->at(i1);
			double ptT =   gphpt->at(i1);
			double etaT = gpheta->at(i1);

			if (abs(etaT)>EtaMax) continue;
			if (ptT<PtMin) continue;

			h_ph_all_pt->Fill(ptT);

			for(unsigned int i2 = 0; i2<phpt->size(); i2++){
				double phi = phphi->at(i2);
				double pt =  phpt->at(i2);
				double eta = pheta->at(i2);
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				if (dR<0.15) h_ph_reco_pt->Fill(ptT);
			}

			for(unsigned int i2 = 0; i2<nnphpt->size(); i2++){
				double phi = nnphphi->at(i2);
				double pt =  nnphpt->at(i2);
				double eta = nnpheta->at(i2);
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				if (dR<0.15) h_ph_nn_pt->Fill(ptT);
			}
		}



	} // end loop over events


	f->Close();


	RootFile->Write();
	RootFile->Print();
	RootFile->Close();

	cout << "Output ROOT file: " << outfile << endl;


	return 0;


}

// S.Chekanov (ANL)

#include "Ana.h"

#include <TRandom2.h>

const double PI   = TMath::Pi();
const double PI2  = 2*PI;
const double PIover2   = 0.5*PI;

int findCeil(int  arr[], int r, int l, int h)
{
	int mid;
	while (l < h)
	{
		mid = l + ((h - l) >> 1); // Same as mid = (l+h)/2
		(r > arr[mid]) ? (l = mid + 1) : (h = mid);
	}
	return (arr[l] >= r) ? l : -1;
}

// The main function that returns a random number from arr[] according to
// distribution array defined by freq[]. n is size of arrays.
int myRand(int arr[], int freq[], int n)
{
	// Create and fill prefix array
	int prefix[n], i;
	prefix[0] = freq[0];
	for (i = 1; i < n; ++i)
		prefix[i] = prefix[i - 1] + freq[i];

	// prefix[n-1] is sum of all frequencies. Generate a random number
	// with value from 1 to this sum
	int r = (rand() % prefix[n - 1]) + 1;

	// Find index of ceiling of r in prefix arrat
	int indexc = findCeil(prefix, r, 0, n - 1);
	return arr[indexc];
}


Int_t Ana::AnalysisJets(vector<LParticle> JetsTrue, vector<LParticle> JetsReco) {



	JetsFastNN.clear();
	m_nnjetpt.clear();
	m_nnjeteta.clear();
	m_nnjetphi.clear();


	m_gjetpt.clear();
	m_gjeteta.clear();
	m_gjetphi.clear();

	m_jetpt.clear();
	m_jeteta.clear();
	m_jetphi.clear();


	const  double EtaMax=maxEta;
	const  double PhiMax=PI;
	const  double intmove=100000; // number for converting to integer frequencies
	const  double delta=2.0/nBinsNN;


	TRandom2 Reta;
	TRandom2 Rphi;


	// fill RECO jets
	for(int j = 0; j<JetsReco.size(); j++){
		LParticle tjet = (LParticle)JetsReco.at(j);
		TLorentzVector L2 = tjet.GetP();
		double phi = L2.Phi();
		double pt =  L2.Perp();
		double eta = L2.PseudoRapidity();
		double ee =  L2.E();
		if (pt>minPT) {
			m_jetpt.push_back(pt);
			m_jeteta.push_back( eta );
			m_jetphi.push_back( phi  );
		}

	}


	for(int j = 0; j<JetsTrue.size(); j++){
		LParticle tjet = (LParticle)JetsTrue.at(j);
		TLorentzVector L2 = tjet.GetP();

		double phiT = L2.Phi();
		double ptT =  L2.Perp();
		double etaT = L2.PseudoRapidity();
		double eeT =  L2.E();

		if (ptT>minPT) {
			m_gjetpt.push_back(ptT);
			m_gjeteta.push_back(etaT);
			m_gjetphi.push_back( phiT );
		}


		h_jetcutflow->Fill(1);

		// cout << "True jet=" << j << " " << ptT << " " << etaT << endl;

		// corrected kinematics
		float pt=ptT;
		float eta=etaT;
		float phi=phiT;
		float ee=eeT;

		//if (phiT<0) phiT=abs(phiT)+PI;
		//if (phi<0) phi=abs(phi)+PI;

		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			double width=dmax-dmin;
			//double dcenter=dmin+0.5*width;

			if (eeT>dmin && eeT<dmax) {
				// Range [0-1]. 0.5 - no difference; 0 - missed reco

				float eeIN=(eeT-dmin)/width;
				float etaIN=etaT/EtaMax;
				float phiIN=phiT/PhiMax;
				float ptIN=(((eeT-ptT)/width)-1);

				fann_type input[num_input];
				input[0] = ptIN;
				input[1] = etaIN;
				input[2] = phiIN;
				input[3] = eeIN;

				fann_type *output = fann_run(ann_jets[m], input);


				fann_type input_eff[3];
				input_eff[0] = eeIN;
				input_eff[1] = etaIN;
				input_eff[2] = phiIN;
				fann_type *output_eff = fann_run(ann_jets_eff[m], input_eff);
				float prob_efficiency=output_eff[0];
				//double R1=Reta.Rndm();
				//cout << "efficiency=" << prob_efficiency << endl;
				h_out5_eff->Fill(prob_efficiency);
				if (prob_efficiency<-0.5) continue; // skip event since low efficiency


				double escale=1.0;
				double eshift=0.0;
			        if (MuPileup>100) {escale=0.5; eshift=-0.5;}


				double scale=0;
				// output of NN can have negaive values. So we add scale to fix those.
				int freqPT[nBinsNN];
				for (int jjj=0; jjj<nBinsNN; jjj++) freqPT[jjj]=int((scale+output[jjj])*intmove);
				//debug:
				//for (unsigned int jjj=0; jjj<nBinsNN; jjj++) cout << freqPT[jjj] << " ";
				//cout << endl;

				srand(time_t(NULL)); // reset random numbers
				int BinSelected=myRand(BinOverTrue, freqPT, nBinsNN); // select random value (bin) assuming frequencies from freqPT
				//int BinSelected=25;
				// this bin defines PT(reco)/PT(true)
				// take into account 0.25 in the ratio when created ratio -> 4 factor
				double recoOvertrue=-1+BinSelected*delta;
				//double ptcor=( (1/escale)*recoOvertrue +1);
				double ptcor= ( (recoOvertrue - eshift)/escale ) +1;
				pt =  ptT *ptcor; // gain
				//cout << "true=" << ptT << " reco=" << pt <<  " Corr=" << ptcor << " selected bin=" << BinSelected << endl;
				h_ptcor->Fill( ptcor );


				// unpack the outputs for pT
				//cout << "\n New jet:" << endl;
				fann_type ptout[nBinsNN];
				for (int jjj=0; jjj<nBinsNN; jjj++) {
					double d1=-1.0+jjj*delta;
					double d2=-1.0+(jjj+1)*delta;
					//if (ptIN>d1 && ptIN<d2) ptout[jjj]= output[jjj];
					h_out1->Fill( d1+delta, output[jjj]);
					//cout << jjj << " " << output[jjj]*intmove << endl;
				}


				// Eta
				int freqETA[nBinsNN];
				for (int jjj=0; jjj<nBinsNN; jjj++) freqETA[jjj]=int((scale+output[jjj+nBinsNN])*intmove);
				BinSelected=myRand(BinOverTrue, freqETA, nBinsNN); // select random value (bin) assuming frequencies
				recoOvertrue=-1+BinSelected*delta;
				double etacor=(0.5*recoOvertrue +1);
				h_etacor->Fill( etacor );
				eta =  etaT * etacor; // gain. Factor 0.5 is due to factor 2 used to buld the ratio in training
				int iii=0;
				for (int jjj=nBinsNN; jjj<2*nBinsNN; jjj++) {
					double d1=-1.0+iii*delta;
					//double d2=-1.0+(iii+1)*delta;
					//if (ptIN>d1 && ptIN<d2) ptout[jjj]= output[jjj];
					//cout << iii << " " << output[jjj] << endl;
					h_out2->Fill( d1+delta, output[jjj]);
					iii++;
				}

				// Phi
				int freqPHI[nBinsNN];
				for (int jjj=0; jjj<nBinsNN; jjj++) freqPHI[jjj]=int( (scale+output[jjj+2*nBinsNN])*intmove);
				BinSelected=myRand(BinOverTrue, freqPHI, nBinsNN); // select random value (bin) assuming frequencies
				recoOvertrue=-1+BinSelected*delta;
				double phicor=(0.5*recoOvertrue +1);
				h_phicor->Fill( phicor );
				phi =  phiT * phicor; // gain. Factor 0.5 is due to factor 2 used to buld the ratio in training
				iii=0;
				for (int jjj=2*nBinsNN; jjj<3*nBinsNN; jjj++) {
					double d1=-1.0+iii*delta;
					//double d2=-1.0+(iii+1)*delta;
					//if (ptIN>d1 && ptIN<d2) ptout[jjj]= output[jjj];
					//cout << iii << " " << output[jjj] << endl;
					h_out3->Fill( d1+delta, output[jjj]);
					iii++;
				}


				
                                // energy
				int freqE[nBinsNN];
				for (int jjj=0; jjj<nBinsNN; jjj++) freqE[jjj]=int( (scale+output[jjj+3*nBinsNN])*intmove);
				BinSelected=myRand(BinOverTrue, freqE, nBinsNN); // select random value (bin) assuming frequencies
				recoOvertrue=-1+BinSelected*delta;

				// correction
				double ecorr=( (recoOvertrue - eshift)/escale ) +1;

				ee =  eeT * ecorr; // gain
				iii=0;
				for (int jjj=3*nBinsNN; jjj<4*nBinsNN; jjj++) {
					double d1=-1.0+iii*delta;
					//double d2=-1.0+(iii+1)*delta;
					//if (ptIN>d1 && ptIN<d2) ptout[jjj]= output[jjj];
					//cout << iii << " " << output[jjj] << endl;
					h_out4->Fill( d1+delta, output[jjj]);
					iii++;
				}



				h_in1->Fill(ptIN);
				h_in2->Fill(etaIN);
				h_in3->Fill(phiIN);
				h_in4->Fill(eeIN);

				



				if (phi>PI || phi<-PI) phi=phiT; // overcorrection!
			}
		}


		h_jetpt_nn->Fill(pt);

		TLorentzVector l;
		l.SetPtEtaPhiE(pt,eta,phi,ee);
		LParticle p;
		p.SetP(l);
		JetsFastNN.push_back(p);
		h_jetcutflow->Fill(2.);


		if (pt>minPT) {
			m_nnjetpt.push_back(pt);
			m_nnjeteta.push_back(eta);
			m_nnjetphi.push_back( phi );
		}


		if (ptT>50)  h_jetres100->Fill(pt/ptT);
		if (ptT>200)  h_jetres1000->Fill(pt/ptT);

	}


	m_ntuple->Fill();



	return 0;

}


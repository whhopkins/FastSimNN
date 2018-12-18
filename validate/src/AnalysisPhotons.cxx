// S.Chekanov (ANL)

#include "Ana.h"


Int_t Ana::AnalysisPhotons(vector<LParticle> True, vector<LParticle> Reco) {

	// NN 
	m_nnphpt.clear();
	m_nnpheta.clear();
	m_nnphphi.clear();

	// true 
	m_gphpt.clear();
	m_gpheta.clear();
	m_gphphi.clear();

	// reco 
	m_phpt.clear();
	m_pheta.clear();
	m_phphi.clear();

	const double EtaMax=EMmaxEta;
	const double PhiMax=PI;
	const double intmove=100000; // number for converting to integer frequencies
	const double delta=2.0/nBinsNN;
	const double slicesEta=(2*EtaMax)/slices_etaphi; // slices in eta
	const double slicesPhi=(2*PhiMax)/slices_etaphi; // slices in phi

	//TRandom2 Rphi;

	// fill RECO photons 
	for(unsigned int j = 0; j<Reco.size(); j++){
		LParticle rph = (LParticle)Reco.at(j);
		TLorentzVector L2 = rph.GetP();
		double phi = L2.Phi();
		double pt =  L2.Perp();
		double eta = L2.PseudoRapidity();

		if (pt>minPT) {
			m_phpt.push_back(pt);
			m_pheta.push_back( eta );
			m_phphi.push_back( phi  );
			//m_phcharge.push_back( charge  );
		}

	}


	for(unsigned int j = 0; j<True.size(); j++){
		LParticle tm = (LParticle)True.at(j);
		TLorentzVector L2 = tm.GetP();

		double phiT = L2.Phi();
		double ptT =  L2.Perp();
		double etaT = L2.PseudoRapidity();
                double isolationT= tm.GetType()/1000.; // isolation 

		if (ptT>minPT && abs(etaT)<EMmaxEta) {
			m_gphpt.push_back(ptT);
			m_gpheta.push_back(etaT);
			m_gphphi.push_back( phiT );
			//m_gphcharge.push_back( chargeT );
		}


		// reco 
		if (ptT>minPT && abs(etaT)<EMmaxEta) { // only for the pT cut
			int indexMatch=-1;
			for(unsigned int i = 0; i<Reco.size(); i++){
				LParticle rm = (LParticle)Reco.at(i);
				TLorentzVector L1 = rm.GetP();
				double phi = L1.Phi();
				double eta = L1.PseudoRapidity();
				double dEta=etaT-eta;
				double dPhi=phiT-phi;
				if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
				double dR=sqrt(dEta*dEta+dPhi*dPhi);
				if (dR<DeltaR) indexMatch=i;
			}

		} // end cut on pT and Eta

		// cout << "True photon=" << j << " " << ptT << " " << etaT << endl;

		// corrected kinematics
		float pt=ptT;
		float eta=etaT;
		float phi=phiT;
		float charge=0;

		//if (phiT<0) phiT=abs(phiT)+PI;
		//if (phi<0) phi=abs(phi)+PI;
                float prob_efficiency=1;

		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			double width=dmax-dmin;

			if (ptT>dmin && ptT<=dmax) {
				// Range [0-1]. 0.5 - no difference; 0 - missed reco
				float dminmax=dmin+0.5*width;
				float ptIN=((ptT-dminmax)/(0.5*width));
				float etaIN=etaT/EtaMax; // range -1 -1
				float phiIN=phiT/PhiMax; // range -1-1 from -pi - pi

				// sliced input for NN
				float etaINSlice[slices_etaphi-1];
				// bin the resolution plots
				for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
					float d1=-EtaMax+jjj*slicesEta;
					float d2=d1+slicesEta;
					//float dmm=d1+0.5*slicesEta;
					if (etaT>d1  && etaT<=d2) etaINSlice[jjj]=1.0f; //     etaINSlice[jjj]=(etaT-dmm)/(0.5*dmm);
					else etaINSlice[jjj]=0;
				}

				float phiINSlice[slices_etaphi-1];
				for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
					float d1=-PhiMax+jjj*slicesPhi;
					float d2=d1+slicesPhi;
					//float dmm=d1+0.5*slicesPhi;
					if (phiT>d1  && phiT<=d2)  phiINSlice[jjj]=1.0f; //    phiINSlice[jjj]=(phiT-dmm)/(0.5*dmm);
					else phiINSlice[jjj]=0;
				}

				fann_type uinput[num_input];

				uinput[0] = ptIN;
				uinput[1] = 0;
				uinput[2] = etaIN;
				uinput[3] = phiIN;

				// eta and phi are sliced for ANN
				// this is needed to reproduce spacial defects
				int shift=4;
				int kshift=0;
				for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
					uinput[shift+kshift] =  etaINSlice[jjj];
					kshift++;
				}
				shift=shift+slices_etaphi-1;
				kshift=0;
				for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
					uinput[shift+kshift] =  phiINSlice[jjj];
					kshift++;
				}


				fann_type * output1 = fann_run(ann1_photons[m], uinput);

				double scale=0;
				// output of NN can have negaive values. So we add scale to fix those.
				//int freqPT[nBinsNN];
				//for (int jjj=0; jjj<nBinsNN; jjj++) freqPT[jjj]=(int)((scale+output[jjj])*intmove);
				//debug:
				//for (unsigned int jjj=0; jjj<nBinsNN; jjj++) cout << freqPT[jjj] << " ";
				//cout << endl;


				// you can use Boost but it is slower
                                /*
				vector<double> freqPT;
                                // create probability distribution
                                double xsum=0; 
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) xsum=xsum+output[jjj];
				for (int jjj=0; jjj<nBinsNN-1; jjj++) freqPT.push_back( (double)(output[jjj]/xsum));
				boost::random::discrete_distribution<> dist1(freqPT);
				int BinSelected = dist1(gen);
				//cout << BinSelected << endl;
                                */

                                int freqPT[nBinsNN-1];
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) freqPT[jjj]=(int)((scale+output1[jjj])*intmove);
				int BinSelected=myRand(BinOverTrue, freqPT, nBinsNN-1); // select random value (bin) assuming frequencies from freqPT
				double recoOvertrue=-1+BinSelected*delta;
				double ptcor= ( (recoOvertrue - em_eshift)/em_escale )+1.0;
				pt =  ptT *ptcor; // gain
				//cout << "true=" << ptT << " reco=" << pt <<  " Corr=" << ptcor << " selected bin=" << BinSelected << endl;
				//h_ptcor->Fill( ptcor );
				//h_rout1->Fill( (float)BinSelected );

				// unpack the outputs for pT
				//cout << "\n New photon:" << endl;
				//for (int jjj=0; jjj<nBinsNN-1; jjj++) {
			        //		double d1=-1.0+jjj*delta;
					//h_out1->Fill( d1+0.5*delta, output1[jjj]);
				//}


				// Eta
				fann_type * output2 = fann_run(ann2_photons[m], uinput);
				int freqETA[nBinsNN-1];
				for (int jjj=0; jjj<nBinsNN-1; jjj++) freqETA[jjj]=(int)((scale+output2[jjj])*intmove);
				BinSelected=myRand(BinOverTrue, freqETA, nBinsNN-1); // select random value (bin) assuming frequencies
                               
                                /* slow boost 
                                vector<double> freqETA;
                                xsum=0;
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) xsum=xsum+output[jjj];
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) freqETA.push_back( (double)(output[jjj]/xsum));
                                boost::random::discrete_distribution<> dist2(freqETA);
                                BinSelected = dist2(gen);
                                */
				recoOvertrue=-1+BinSelected*delta;
				double etacor= ( (recoOvertrue - em_etashift)/em_etascale )+1.0;
				//h_etacor->Fill( etacor );
				eta =  etaT * etacor;
				//for (int jjj=0; jjj<nBinsNN-1; jjj++) {
			        //		double d1=-1.0+jjj*delta;
					//h_out2->Fill( d1+0.5*delta, output2[jjj]);
				//}

				// Phi
				fann_type * output3 = fann_run(ann3_photons[m], uinput);
				int freqPHI[nBinsNN-1];
				for (int jjj=0; jjj<nBinsNN-1; jjj++) freqPHI[jjj]=(int)((scale+output3[jjj])*intmove);
				BinSelected=myRand(BinOverTrue, freqPHI, nBinsNN-1); // select random value (bin) assuming frequencies
                                /*
                                vector<double> freqPHI;
                                xsum=0;
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) xsum=xsum+output[jjj];
                                for (int jjj=0; jjj<nBinsNN-1; jjj++) freqPHI.push_back( (double)(output[jjj]/xsum));
                                boost::random::discrete_distribution<> dist3(freqPHI);
                                BinSelected = dist3(gen);
                                */

				//h_rout3->Fill( (float)BinSelected );
				recoOvertrue=-1+BinSelected*delta;
				double phicor= ( (recoOvertrue - em_etashift)/em_etascale )+1.0;
				//h_phicor->Fill( phicor );
				phi =  phiT * phicor;
				//for (int jjj=0; jjj<nBinsNN; jjj++) {
				//	double d1=-1.0+jjj*delta;
					//h_out3->Fill( d1+0.5*delta, output3[jjj]);
				//}


				if (phi>PI || phi<-PI) phi=phiT; // overcorrection!

                                // ----------------- feature NN starts here --------------------
                                // structure of input the same but mass is replaced with b-tagging
                                uinput[0] = ptIN;
                                uinput[1] = etaIN;
                                uinput[2] = phiIN;
                                uinput[3] = isolationT;

				fann_type * output5 = fann_run(ann5_photons[m], uinput);
				prob_efficiency=output5[0];
				charge=output5[1];
                                //cout << "B-tag NN input=" << uinput[3] << " btagT=" << btagT << " NN out=" << btagFound << " " << endl;


			}
		}


                h_ph_out5_eff->Fill(prob_efficiency);

		if (pt>minPT && prob_efficiency>0) {
			m_nnphpt.push_back(pt);
			m_nnpheta.push_back(eta);
			m_nnphphi.push_back(phi);
			//m_nnelcharge.push_back(charge);
                        // cout << btagFound << " " << bt << endl;
		}



	}




	return 0;

}


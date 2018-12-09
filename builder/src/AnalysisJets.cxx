#include "Ana.h"


const double PI   = TMath::Pi();
const double PI2  = 2*PI;
const double PIover2   = 0.5*PI;


Int_t Ana::AnalysisJets(vector<LParticle> JetsTrue, vector<LParticle> JetsReco) {


	const double EtaMax=maxEta;
	const double PhiMax=PI;
	const double delta=2.0/nBinsNN;


	for(unsigned int j = 0; j<JetsTrue.size(); j++){
		LParticle tjet = (LParticle)JetsTrue.at(j);
		TLorentzVector L2 = tjet.GetP();
		double phiT = L2.Phi();
		double ptT =  L2.Perp();
		double etaT = L2.PseudoRapidity();
		double massT =  L2.M();
                float  btagT=(float)tjet.GetType(); // get  b-quark in 100%

		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			//double width=dmax-dmin;
			//double dcenter=dmin+0.5*width;
			if (ptT>dmin && ptT<dmax)  eventsBins[m]++;
		}


		// reco jets
		int indexMatch=-1;
		for(unsigned int i = 0; i<JetsReco.size(); i++){
			LParticle rjet = (LParticle)JetsReco.at(i);
			TLorentzVector L1 = rjet.GetP();
			double phi = L1.Phi();
			double eta = L1.PseudoRapidity();
			double dEta=etaT-eta;
			double dPhi=phiT-phi;
			if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
			double dR=sqrt(dEta*dEta+dPhi*dPhi);
			h_dR->Fill(dR);
			if (dR<DeltaR) indexMatch=i;
		}

                float btag=0;
		if (indexMatch>-1){
			LParticle rjet = (LParticle)JetsReco.at(indexMatch);
			TLorentzVector L1 = rjet.GetP();
			float phi = L1.Phi();
			float pt  = L1.Perp();
			float eta = L1.PseudoRapidity();
			float mass  = L1.M();
                        btag  = (float)rjet.GetType();

			vector<float> input;
			vector<float> output;
			input.push_back((float)ptT);
			input.push_back((float)etaT);
			input.push_back((float)phiT);
			input.push_back( (float)massT);

			output.push_back(pt);
			output.push_back(eta);
			output.push_back(phi);
			output.push_back(mass);

			// push vectors
			finput_jets.push_back(input);
			foutput_jets.push_back(output);
		}

		// statistics limitted for efficiency
		vector<float> input2;
		vector<float> output2;
		input2.push_back((float)ptT);
		input2.push_back((float)etaT);
		input2.push_back((float)phiT);
                input2.push_back((float)btagT); // fraction of b-quark momenta in % (100-0) 

                float ibb=-1.0f;
		float iout=-1.0f; // no match
		if (indexMatch>-1)  iout=1.0f;
                if (btag>0)         ibb=1.0f;  // b-tagged and matched

                output2.push_back(iout);
                output2.push_back(ibb);

		finput_jets_eff.push_back(input2);
		foutput_jets_eff.push_back(output2);


	}


	if (nevv%nBatch == 0  && nevv>0) {
		cout << "### Reset MSE after =" << nevv << " training " << nEpoch << " epoches " << endl;



		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			double width=dmax-dmin;
			if (eventsBins[m]<MinEntries) continue;
			cout << "-> Start training for energy bin = " << m << endl;

			// empty data
			fann_train_data *  dataset=  fann_create_train(eventsBins[m], num_input, num_output);

			// create a dataset for a given bin
			int nn=0;
			for (unsigned int i=0; i<finput_jets.size(); i++){
				vector<float> input = finput_jets[i];
				vector<float> output = foutput_jets[i];
				float ptT=input[0];
				float etaT=input[1];
				float phiT=input[2];
				float massT=input[3];

				float pt=output[0];
				float eta=output[1];
				float phi=output[2];
				float mass=output[3];


				if (ptT>dmin && ptT<dmax) {

                                        float dminmax=dmin+0.5*width;
					float ptIN=((ptT-dminmax)/(0.5*width));
					float etaIN=etaT/EtaMax; // range -1 -1
					float phiIN=phiT/PhiMax; // range -1-1 from -pi - pi
					float mIN=-1+(massT/(0.5*dmin));

	                   		// assume -1 - 1 range.
					// scale to avoid sharp behaviour near 0
					// for pileup we make the distribution narrower and shift to -0.5 to fit to [-1,1]
					float ptOUT= jet_escale*((pt/ptT)-1) + jet_eshift;
					//float etaOUT=norm*(eta-etaT);
					float etaOUT=0;
					// make broader rescale for a better efficiency for -1 and 1 range
					if (etaT !=0) etaOUT=jet_etascale*(abs(eta/etaT)-1) + jet_etashift;
					//float phiOUT=norm*(phi-phiT);
					float  phiOUT=0;
					if (phiT !=0) phiOUT=jet_etascale*(abs(phi/phiT)-1) +  jet_etashift; 
					//float eeOUT=abs(ee-eeT)/dmax;
					// shrink pT (important for pileup!)
					float mOUT=jet_mscale*((mass/massT) -1) + jet_mshift;
					//float shiftOUT=0.0f;
					//if (ee-eeT>0) shiftOUT=1.0f; // gain or positive shift

					fann_type uinput[num_input];
					fann_type uoutput[num_output];

					uinput[0] = ptIN;
					uinput[1] = etaIN;
					uinput[2] = phiIN;
					uinput[3] = mIN;
					// outputs
					uoutput[0] = ptOUT;
					uoutput[1] = etaOUT;
					uoutput[2] = phiOUT;
					uoutput[3] = mOUT;

					for (unsigned int jjj=0; jjj<num_output; jjj++) uoutput[jjj]=0.0;

					// bin the resolution plots
					for (unsigned int jjj=0; jjj<nBinsNN; jjj++) {
						double d1=-1.0+jjj*delta;
						double d2=-1.0+(jjj+1)*delta;
						if (ptOUT>d1  && ptOUT<=d2)    uoutput[jjj]=1.0;
						if (etaOUT>d1 && etaOUT<=d2)   uoutput[jjj+nBinsNN]=1.0;
						if (phiOUT>d1 && phiOUT<=d2)   uoutput[jjj+2*nBinsNN]=1.0;
						if (mOUT>d1  && mOUT<=d2)      uoutput[jjj+3*nBinsNN]=1.0;
					}

					for (unsigned int kk=0; kk<num_input; kk++)  dataset->input[nn][kk] =uinput[kk];
					for (unsigned int kk=0; kk<num_output; kk++) dataset->output[nn][kk] =uoutput[kk];

					h_in1->Fill(ptIN);
					h_in2->Fill(etaIN);
					h_in3->Fill(phiIN);
					h_in4->Fill(mIN);

					h_out1->Fill(ptOUT);
					h_out2->Fill(etaOUT);
					h_out3->Fill(phiOUT);
					h_out4->Fill(mOUT);

					nn++;
				} // end energy range


			} // end over data

			double mse=0;
			for (int e=0; e<nEpoch; e++) {
                          mse = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann_jets[m], dataset, num_threads) : fann_train_epoch(ann_jets[m], dataset);
    			  if (e%100==0 ||  (e<10) || (e%20==0)) cout << "bin=" << m << " epoch=" << e << " MSE=" << mse << endl;
			  if (mse<MSESTOP) break;
			}
			cout << "  Bin=" << m << " with " << eventsBins[m] << " events has MSE=" <<  mse << " after epoch Nr=" << nEpoch << endl;


			fann_destroy_train(dataset) ; // clear

			// empty data
			fann_train_data *  dataset_eff=  fann_create_train(eventsBins[m], num_input_eff, num_output_eff);
			nn=0;
			cout << "\n  -> Training for efficiency  = " << m << " sample size=" << finput_jets_eff.size() << endl;
			for (unsigned int i=0; i<finput_jets_eff.size(); i++){
				//cout << i << endl;
				vector<float> input2 = finput_jets_eff[i];
				vector<float> output2 =  foutput_jets_eff[i];
				float ptT=input2[0];
				float etaT=input2[1];
				float phiT=input2[2];
                                float btagT=(float)input2[3]; // 0 - 100 
                                // outputs
				float match=output2[0];
                                float btag=output2[1];
				h_out5->Fill(match);
                                h_out6->Fill(btag); // -1 or 1 

				if (ptT>dmin && ptT<dmax) {
				        float dminmax=dmin+0.5*width;
                                        float ptIN=((ptT-dminmax)/(0.5*width));
					float etaIN=etaT/EtaMax; // range -1 -1
					float phiIN=phiT/PhiMax; // range -1-1 from -pi - pi
					fann_type uinput[num_input_eff];
					fann_type uoutput[num_output_eff];
					uinput[0]=ptIN;
					uinput[1]=etaIN;
					uinput[2]=phiIN;
                                        uinput[3]=(float)(btagT/100.) - 1; // normalize  -1 - 0 
                                        if (uinput[3]>1) uinput[3]=1;

					uoutput[0]=(float)match;
                                        uoutput[1]=(float)btag;

                                        // debug
                                        //if (uoutput[0]>0 && btag>0) cout << "In " << uinput[3]  << "  out=" << uoutput[1] << endl;

					for (unsigned int kk=0; kk<4; kk++)  dataset_eff->input[nn][kk] =uinput[kk];
					for (unsigned int kk=0; kk<2; kk++)  dataset_eff->output[nn][kk] =uoutput[kk];

					nn++;
				} // end energy range
			}// end of dataset

			for (int e=0; e<nEpoch*4; e++) {
                            float mmse = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann_jets_eff[m], dataset_eff, num_threads) : fann_train_epoch(ann_jets_eff[m], dataset_eff);

			     if (e%100==0 || (e<10)) cout << "bin=" << m << " epoch=" << e << " MSE=" << mmse << endl;
		             if (mse<MSESTOP) break;
			}
			fann_destroy_train(dataset_eff) ; // clear


		} // end run over energy bins



		cout << "### MSE for this batch after " << nevv << endl;
		for (int m=0; m<nBins-1; m++) eventsBins[m] =0;


		double sumMSE=0;
		for (int m=0; m<nBins-1; m++){
			double mse=fann_get_MSE(ann_jets[m]);
			cout << "Event=" << nevv << " Bin=" << m << "  MSE=" << mse << endl;
			if (mse<1000) sumMSE=sumMSE+mse;
		};
		cout << "## Total MSE for all bins for this batch " << sumMSE << endl;
		// batch if
		finput_jets.clear();
		foutput_jets.clear();
		finput_jets_eff.clear();
		foutput_jets_eff.clear();

	}

	nevv++;

	return 0;

}


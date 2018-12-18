// S.Chekanov (ANL)

#include "Ana.h"

Int_t Ana::AnalysisMuons(vector<LParticle> True, vector<LParticle> Reco) {

	const double EtaMax=maxEta;
	const double PhiMax=PI;
	const double delta=2.0/(nBinsNN);                // slices for resolution 
        const double slicesEta=(2*EtaMax)/slices_etaphi; // slices in eta  
        const double slicesPhi=(2*PhiMax)/slices_etaphi; // slices in phi  

	for(unsigned int j = 0; j<True.size(); j++){
		LParticle tm = (LParticle)True.at(j);
		TLorentzVector L2 = tm.GetP();
		double phiT = L2.Phi();
		double ptT =  L2.Perp();
		double etaT = L2.PseudoRapidity();
		int    chargeT =  tm.GetCharge();
                double isolationT= tm.GetType()/1000.; // isolation 

		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			if (ptT>dmin && ptT<=dmax)  eventsMuonBins[m]++;
		}


		// reco 
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

                float charge=0;
		if (indexMatch>-1){
			LParticle rm = (LParticle)Reco.at(indexMatch);
			TLorentzVector L1 = rm.GetP();
			float phi = L1.Phi();
			float pt  = L1.Perp();
			float eta = L1.PseudoRapidity();
			charge  = rm.GetCharge();

			vector<float> input;
			vector<float> output;
			input.push_back((float)ptT);
			input.push_back((float)etaT);
			input.push_back((float)phiT);
			input.push_back((float)chargeT);

			output.push_back(pt);
			output.push_back(eta);
			output.push_back(phi);
			output.push_back(charge);

			// push vectors
			finput_muons.push_back(input);
			foutput_muons.push_back(output);
		}

		// statistics limitted for efficiency. No matching 
		vector<float> input2;
		vector<float> output2;
		input2.push_back((float)ptT);
		input2.push_back((float)etaT);
		input2.push_back((float)phiT);
                input2.push_back((float)chargeT);   //  some feature 
                input2.push_back((float)isolationT); //  some feature 

                float ibb=-1.0f;
		float iout=-1.0f; // no match
		if (indexMatch>-1)  {iout=1.0f;
                                     if (charge >0) ibb=1.0f;
                                    }

                output2.push_back(iout);
                output2.push_back(ibb);

		finput_muons_eff.push_back(input2);
		foutput_muons_eff.push_back(output2);


	}


	if (nevv%nBatch == 0  && nevv>0) {
                cout << "\033[1;31m Training of muons \033[0m\n";
		cout << "Muons: ### Reset MSE after =" << nevv << " training " << nEpoch << " epoches " << endl;
		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			double width=dmax-dmin;
			if (eventsMuonBins[m]<1) continue;
			cout << "-> Muons: Training for energy bin = " << m << endl;

			// empty data
			fann_train_data *  dataset1=  fann_create_train(eventsMuonBins[m], num_input, num_output);
                        fann_train_data *  dataset2=  fann_create_train(eventsMuonBins[m], num_input, num_output);
                        fann_train_data *  dataset3=  fann_create_train(eventsMuonBins[m], num_input, num_output);
                        fann_train_data *  dataset4=  fann_create_train(eventsMuonBins[m], num_input, num_output);

			// create a dataset for a given bin
			int nn=0;
			for (unsigned int i=0; i<finput_muons.size(); i++){
				vector<float> input = finput_muons[i];
				vector<float> output = foutput_muons[i];
				float ptT=input[0];
				float etaT=input[1];
				float phiT=input[2];
				float chargeT=input[3];

				float pt=output[0];
				float eta=output[1];
				float phi=output[2];
				float charge=output[3];


				if (ptT>dmin && ptT<=dmax) {

                                        float dminmax=dmin+0.5*width;
					float ptIN=((ptT-dminmax)/(0.5*width));
					float etaIN=etaT/EtaMax; // range -1 -1
					float phiIN=phiT/PhiMax; // range -1-1 from -pi - pi
					float chargeIN=chargeT;

                                        // sliced input for NN
                                        float etaINSlice[slices_etaphi-1];
                                        // bin the resolution plots
                                        for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
                                                float d1=-EtaMax+jjj*slicesEta;
                                                float d2=d1+slicesEta;
                                                //float dmm=d1+0.5*slicesEta;
                                                if (etaT>d1  && etaT<=d2) etaINSlice[jjj]=1.0f; //    etaINSlice[jjj]=(etaT-dmm)/(0.5*dmm);
                                                else etaINSlice[jjj]=0;
                                        }

                                        float phiINSlice[slices_etaphi-1];
                                        for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
                                                float d1=-PhiMax+jjj*slicesPhi;
                                                float d2=d1+slicesPhi;
                                                //float dmm=d1+0.5*slicesPhi;
                                                if (phiT>d1  && phiT<=d2) phiINSlice[jjj]=1.0f; //    phiINSlice[jjj]=(phiT-dmm)/(0.5*dmm);
                                                else phiINSlice[jjj]=0;
                                        }


	                   		// assume -1 - 1 range.
					// scale to avoid sharp behaviour near 0
					// for pileup we make the distribution narrower and shift to -0.5 to fit to [-1,1]
					float ptOUT= em_escale*((pt/ptT)-1) + em_eshift;
					//float etaOUT=norm*(eta-etaT);
					float etaOUT=0;
					// make broader rescale for a better efficiency for -1 and 1 range
					if (etaT !=0) etaOUT=em_etascale*(abs(eta/etaT)-1) + em_etashift;
					//float phiOUT=norm*(phi-phiT);
					float  phiOUT=0;
					if (phiT !=0) phiOUT=em_etascale*(abs(phi/phiT)-1) +  em_etashift; 
					//float eeOUT=abs(ee-eeT)/dmax;
					// shrink pT (important for pileup!)
					//float chargeOUT=charge;
					//float shiftOUT=0.0f;
					//if (ee-eeT>0) shiftOUT=1.0f; // gain or positive shift

					fann_type uinput[num_input];
                                        // outputs
					fann_type uoutput1[num_output];
                                        fann_type uoutput2[num_output];
                                        fann_type uoutput3[num_output];
                                        fann_type uoutput4[num_output];
 
					uinput[0] = ptIN;
                                        uinput[1] = 0.0f; // not used 
                                        uinput[2] = etaIN;
                                        uinput[3] = phiIN;
                                        uinput[4] = 0.0f; // not used 

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

                                        // debug input
                                        //cout <<  " " << endl;
                                        //for (int jjj=0; jjj<num_input; jjj++) {cout << uinput[jjj] << " ";}
                                        //cout <<  " " << endl;


					for (int jjj=0; jjj<num_output; jjj++) {uoutput1[jjj]=0; 
                                                                                uoutput2[jjj]=0; 
                                                                                uoutput3[jjj]=0; 
                                                                                uoutput4[jjj]=0; } 

					// bin the resolution plots
					for (int jjj=0; jjj<nBinsNN-1; jjj++) {
						float d1=-1.0+jjj*delta;
						float d2=d1+delta;
						if (ptOUT>d1  && ptOUT<=d2)    uoutput1[jjj]=1.0f;
						if (etaOUT>d1 && etaOUT<=d2)   uoutput2[jjj]=1.0f;
						if (phiOUT>d1 && phiOUT<=d2)   uoutput3[jjj]=1.0f;
						uoutput4[jjj]=0;
					}

					for (int kk=0; kk<num_input; kk++)  {
                                                       dataset1->input[nn][kk] =uinput[kk];
                                                       dataset2->input[nn][kk] =uinput[kk];
                                                       dataset3->input[nn][kk] =uinput[kk];
                                                       dataset4->input[nn][kk] =uinput[kk];
                                        }

 
					for (int kk=0; kk<num_output; kk++) {
                                                        dataset1->output[nn][kk] =uoutput1[kk];
                                                        dataset2->output[nn][kk] =uoutput2[kk];
                                                        dataset3->output[nn][kk] =uoutput3[kk];
                                                        dataset4->output[nn][kk] =uoutput4[kk];
                                        }

                                        // debug output 
                                        //cout <<  " " << endl;
                                        //for (int jjj=0; jjj<nBinsNN-1; jjj++) {cout << uoutput1[jjj] << " ";}
                                        //cout <<  " " << endl;

					h_mu_in1->Fill(ptIN);
					h_mu_in2->Fill(etaIN);
					h_mu_in3->Fill(phiIN);
					h_mu_in4->Fill(chargeIN);

					h_mu_out1->Fill(ptOUT);
					h_mu_out2->Fill(etaOUT);
					h_mu_out3->Fill(phiOUT);
					h_mu_out4->Fill(charge);

					nn++;
				} // end energy range


			} // end over data

                        double mse=0;
			double mse1=0;
                        double mse2=0; 
                        double mse3=0;
                        double mse4=0;
			for (int e=0; e<nEpoch; e++) {
                          mse1 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann1_muons[m], dataset1, num_threads) : fann_train_epoch(ann1_muons[m], dataset1);
                          mse2 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann2_muons[m], dataset2, num_threads) : fann_train_epoch(ann2_muons[m], dataset2);
                          mse3 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann3_muons[m], dataset3, num_threads) : fann_train_epoch(ann3_muons[m], dataset3);
                          mse4 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann4_muons[m], dataset4, num_threads) : fann_train_epoch(ann4_muons[m], dataset4);

                          mse=mse1+mse2+mse3+mse4;  
    			  if (e%100==0 ||  (e<10) || (e%20==0)) cout << " muon epoch=" << e << " MSE=" << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << mse1 << ", " <<  mse2 << ", " << mse3 << ", " << mse4 << " tot=" << mse << endl;
			  if (mse<MSESTOP) break;
			}
			cout << "  Bin=" << m << " with " << eventsMuonBins[m] << " events has total MSE=" <<  mse  << " after epoch Nr=" << nEpoch << endl;


			fann_destroy_train(dataset1) ; // clear
                        fann_destroy_train(dataset2) ; // clear
                        fann_destroy_train(dataset3) ; // clear
                        fann_destroy_train(dataset4) ; // clear

			// empty data for feature NN 
			fann_train_data *  dataset_eff=  fann_create_train(eventsMuonBins[m], num_input_eff, num_output_eff);
			nn=0;
    			cout << "\n  Muons: -> Training for feature/efficiency  = " << m << " sample size=" << finput_muons_eff.size() << endl;
			for (unsigned int i=0; i<finput_muons_eff.size(); i++){
				//cout << i << endl;
				vector<float> input2 = finput_muons_eff[i];
				vector<float> output2 =  foutput_muons_eff[i];
				float ptT=input2[0];
				float etaT=input2[1];
				float phiT=input2[2];
                                float chargeT=input2[3]; // some feature 
                                float isolationT=input2[4]; // some feature 

                                // outputs
				float match=output2[0];
                                float charge=output2[1];
                                h_mu_out5->Fill(match);

                          
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
					uinput[3]=isolationT;
                                        uinput[4]=chargeT; // some feature

                                       // sliced input for NN
                                        float etaINSlice[slices_etaphi-1];
                                        // bin the resolution plots
                                        for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
                                                float d1=-EtaMax+jjj*slicesEta;
                                                float d2=d1+slicesEta;
                                                //float dmm=d1+0.5*slicesEta;
                                                if (etaT>d1  && etaT<=d2) etaINSlice[jjj]=1.0f; //    etaINSlice[jjj]=(etaT-dmm)/(0.5*dmm);
                                                else etaINSlice[jjj]=0;
                                        }

                                        float phiINSlice[slices_etaphi-1];
                                        for (int jjj=0; jjj<slices_etaphi-1; jjj++) {
                                                float d1=-PhiMax+jjj*slicesPhi;
                                                float d2=d1+slicesPhi;
                                                //float dmm=d1+0.5*slicesPhi;
                                                if (phiT>d1  && phiT<=d2) phiINSlice[jjj]=1.0f; //    phiINSlice[jjj]=(phiT-dmm)/(0.5*dmm);
                                                else phiINSlice[jjj]=0;
                                        }


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


                                        //debug input 
                                        //cout <<  " " << endl;
                                        //for (int jjj=0; jjj<num_input_eff; jjj++) {cout << uinput[jjj] << " ";}
                                        //cout <<  " " << endl;


                                        // outputs
					uoutput[0]=(float)match;
                                        uoutput[1]=charge;

					for (unsigned int kk=0; kk<num_input_eff; kk++)  dataset_eff->input[nn][kk] =uinput[kk];
					for (unsigned int kk=0; kk<num_output_eff; kk++)  dataset_eff->output[nn][kk] =uoutput[kk];

					nn++;
				} // end energy range
			}// end of dataset

			for (int e=0; e<nEpoch*2; e++) {
                            float mmse = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann5_muons[m], dataset_eff, num_threads) : fann_train_epoch(ann5_muons[m], dataset_eff);

			     if (e%100==0 || (e<10)) cout << "bin=" << m << " epoch=" << e << " MSE=" << mmse << endl;
		             if (mse<MSESTOP) break;
			}
			fann_destroy_train(dataset_eff) ; // clear


		} // end run over energy bins



		cout << "Muons: ### MSE for this batch after " << nevv << endl;
		for (int m=0; m<nBins-1; m++) eventsMuonBins[m] =0;


		double sumMSE=0;
		for (int m=0; m<nBins-1; m++){
			double mse=fann_get_MSE(ann1_muons[m]);
			cout << "Event=" << nevv << " Bin=" << m << "  MSE=" << mse << endl;
			if (mse<1000) sumMSE=sumMSE+mse;
		};
		cout << "Muons: ## Total MSE for all bins for this batch " << sumMSE << endl;
		// batch if
		finput_muons.clear();
		foutput_muons.clear();
		finput_muons_eff.clear();
		foutput_muons_eff.clear();

	}

	return 0;

}


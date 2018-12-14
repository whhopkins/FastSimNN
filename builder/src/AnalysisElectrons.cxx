// S.Chekanov (ANL)

#include "Ana.h"


Int_t Ana::AnalysisElectrons(vector<LParticle> ElectronsTrue, vector<LParticle> ElectronsReco) {

	const double EtaMax=maxEta;
	const double PhiMax=PI;
	const double delta=2.0/(nBinsNN);                // slices for resolution 
        const double slicesEta=(2*EtaMax)/slices_etaphi; // slices in eta  
        const double slicesPhi=(2*PhiMax)/slices_etaphi; // slices in phi  

	for(unsigned int j = 0; j<ElectronsTrue.size(); j++){
		LParticle tele = (LParticle)ElectronsTrue.at(j);
		TLorentzVector L2 = tele.GetP();
		double phiT = L2.Phi();
		double ptT =  L2.Perp();
		double etaT = L2.PseudoRapidity();
		double massT =  L2.M();

		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			if (ptT>dmin && ptT<=dmax)  eventsElectronBins[m]++;
		}

		// reco 
		int indexMatch=-1;
		for(unsigned int i = 0; i<ElectronsReco.size(); i++){
			LParticle rele = (LParticle)ElectronsReco.at(i);
			TLorentzVector L1 = rele.GetP();
			double phi = L1.Phi();
			double eta = L1.PseudoRapidity();
			double dEta=etaT-eta;
			double dPhi=phiT-phi;
			if (abs(dPhi)>PI) dPhi=PI2-abs(dPhi);
			double dR=sqrt(dEta*dEta+dPhi*dPhi);
			h_dR->Fill(dR);
			if (dR<DeltaR) indexMatch=i;
		}

		if (indexMatch>-1){
			LParticle rmuon = (LParticle)ElectronsReco.at(indexMatch);
			TLorentzVector L1 = rmuon.GetP();
			float phi = L1.Phi();
			float pt  = L1.Perp();
			float eta = L1.PseudoRapidity();
			float mass  = L1.M();

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
			finput_electrons.push_back(input);
			foutput_electrons.push_back(output);
		}

		// statistics limitted for efficiency. No matching 
		vector<float> input2;
		vector<float> output2;
		input2.push_back((float)ptT);
		input2.push_back((float)etaT);
		input2.push_back((float)phiT);
                input2.push_back((float)0); //  some feature 

                float ibb=-1.0f;
		float iout=-1.0f; // no match
		if (indexMatch>-1)  iout=1.0f;

                output2.push_back(iout);
                output2.push_back(ibb);

		finput_electrons_eff.push_back(input2);
		foutput_electrons_eff.push_back(output2);


	}


	if (nevv%nBatch == 0  && nevv>0) {
                cout << "\033[1;31m Training of electrons \033[0m\n";
		cout << "Electrons: ### Reset MSE after =" << nevv << " training " << nEpoch << " epoches " << endl;
		for (int m=0; m<nBins-1; m++){
			double dmin=eBins[m];
			double dmax=eBins[m+1];
			double width=dmax-dmin;
			if (eventsElectronBins[m]<1) continue;
			cout << "-> Electrons: Training for energy bin = " << m << endl;

			// empty data
			fann_train_data *  dataset1=  fann_create_train(eventsElectronBins[m], num_input, num_output);
                        fann_train_data *  dataset2=  fann_create_train(eventsElectronBins[m], num_input, num_output);
                        fann_train_data *  dataset3=  fann_create_train(eventsElectronBins[m], num_input, num_output);
                        fann_train_data *  dataset4=  fann_create_train(eventsElectronBins[m], num_input, num_output);

			// create a dataset for a given bin
			int nn=0;
			for (unsigned int i=0; i<finput_electrons.size(); i++){
				vector<float> input = finput_electrons[i];
				vector<float> output = foutput_electrons[i];
				float ptT=input[0];
				float etaT=input[1];
				float phiT=input[2];
				float massT=input[3];

				float pt=output[0];
				float eta=output[1];
				float phi=output[2];
				float mass=output[3];


				if (ptT>dmin && ptT<=dmax) {

                                        float dminmax=dmin+0.5*width;
					float ptIN=((ptT-dminmax)/(0.5*width));
					float etaIN=etaT/EtaMax; // range -1 -1
					float phiIN=phiT/PhiMax; // range -1-1 from -pi - pi
					float massIN=-1+(massT/(0.5*dmin));

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
					float mOUT=em_mscale*((mass/massT) -1) + em_mshift;
					//float shiftOUT=0.0f;
					//if (ee-eeT>0) shiftOUT=1.0f; // gain or positive shift

					fann_type uinput[num_input];
                                        // outputs
					fann_type uoutput1[num_output];
                                        fann_type uoutput2[num_output];
                                        fann_type uoutput3[num_output];
                                        fann_type uoutput4[num_output];
 
					uinput[0] = ptIN;
                                        uinput[1] = massIN;
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
						if (mOUT>d1   && mOUT<=d2)     uoutput4[jjj]=1.0f;
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

					h_in1->Fill(ptIN);
					h_in2->Fill(etaIN);
					h_in3->Fill(phiIN);
					h_in4->Fill(massIN);

					h_out1->Fill(ptOUT);
					h_out2->Fill(etaOUT);
					h_out3->Fill(phiOUT);
					h_out4->Fill(mOUT);

					nn++;
				} // end energy range


			} // end over data

                        double mse=0;
			double mse1=0;
                        double mse2=0; 
                        double mse3=0;
                        double mse4=0;
			for (int e=0; e<nEpoch; e++) {
                          mse1 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann1_electrons[m], dataset1, num_threads) : fann_train_epoch(ann1_electrons[m], dataset1);
                          mse2 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann2_electrons[m], dataset2, num_threads) : fann_train_epoch(ann2_electrons[m], dataset2);
                          mse3 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann3_electrons[m], dataset3, num_threads) : fann_train_epoch(ann3_electrons[m], dataset3);
                          mse4 = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann4_electrons[m], dataset4, num_threads) : fann_train_epoch(ann4_electrons[m], dataset4);

                          mse=mse1+mse2+mse3+mse4;  
    			  if (e%100==0 ||  (e<10) || (e%20==0)) cout << " electron epoch=" << e << " MSE=" << std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << mse1 << ", " <<  mse2 << ", " << mse3 << ", " << mse4 << " tot=" << mse << endl;
			  if (mse<MSESTOP) break;
			}
			cout << "  Bin=" << m << " with " << eventsElectronBins[m] << " events has total MSE=" <<  mse  << " after epoch Nr=" << nEpoch << endl;


			fann_destroy_train(dataset1) ; // clear
                        fann_destroy_train(dataset2) ; // clear
                        fann_destroy_train(dataset3) ; // clear
                        fann_destroy_train(dataset4) ; // clear

			// empty data for feature NN 
			fann_train_data *  dataset_eff=  fann_create_train(eventsElectronBins[m], num_input_eff, num_output_eff);
			nn=0;
    			cout << "\n  Electrons: -> Training for efficiency  = " << m << " sample size=" << finput_electrons_eff.size() << endl;
			for (unsigned int i=0; i<finput_electrons_eff.size(); i++){
				//cout << i << endl;
				vector<float> input2 = finput_electrons_eff[i];
				vector<float> output2 =  foutput_electrons_eff[i];
				float ptT=input2[0];
				float etaT=input2[1];
				float phiT=input2[2];
                                float btagT=input2[3]; // some feature 
                                // outputs
				float match=output2[0];
                                float btag=output2[1];

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
                                        uinput[3]=0; // some feature

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
                                        uoutput[1]=0;

					for (unsigned int kk=0; kk<num_input_eff; kk++)  dataset_eff->input[nn][kk] =uinput[kk];
					for (unsigned int kk=0; kk<num_output_eff; kk++)  dataset_eff->output[nn][kk] =uoutput[kk];

					nn++;
				} // end energy range
			}// end of dataset

			for (int e=0; e<nEpoch*2; e++) {
                            float mmse = num_threads > 1 ? fann_train_epoch_irpropm_parallel(ann5_electrons[m], dataset_eff, num_threads) : fann_train_epoch(ann5_electrons[m], dataset_eff);

			     if (e%100==0 || (e<10)) cout << "bin=" << m << " epoch=" << e << " MSE=" << mmse << endl;
		             if (mse<MSESTOP) break;
			}
			fann_destroy_train(dataset_eff) ; // clear


		} // end run over energy bins



		cout << "Electrons: ### MSE for this batch after " << nevv << endl;
		for (int m=0; m<nBins-1; m++) eventsElectronBins[m] =0;


		double sumMSE=0;
		for (int m=0; m<nBins-1; m++){
			double mse=fann_get_MSE(ann1_electrons[m]);
			cout << "Event=" << nevv << " Bin=" << m << "  MSE=" << mse << endl;
			if (mse<1000) sumMSE=sumMSE+mse;
		};
		cout << "Electrons: ## Total MSE for all bins for this batch " << sumMSE << endl;
		// batch if
		finput_electrons.clear();
		foutput_electrons.clear();
		finput_electrons_eff.clear();
		foutput_electrons_eff.clear();

	}

	return 0;

}


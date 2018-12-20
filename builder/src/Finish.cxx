#include "Ana.h"

// Finish all calculations at the end   
Int_t Ana::Finish() { 


    cout << "\n\n-- Write output file=" << ffile << endl;
    cout << "\n\n";

    RootFile->Write();
    RootFile->Print();
    RootFile->Close();


/*
   // compare with random RMM
   const unsigned int num_input = NrOfInputs;
   const unsigned int num_output = NrOfInputs;
   struct fann *ann_random   = fann_create_standard(num_layers, num_input, num_neurons_hidden_1, num_output);

    fann_set_activation_function_hidden(ann_random, FANN_SIGMOID);
    fann_set_activation_function_output(ann_random, FANN_SIGMOID);
    fann_randomize_weights(ann_random,0,0.999);
    double mse_random=fann_get_MSE(ann_random);
    fann_destroy(ann_random);
*/


   //config_setting_t *root, *setting, *group, *array;
   static const char *output_file = "out_ann/config.cfg";
   Config cfg;
   Setting &root = cfg.getRoot();
   // general settings
   Setting &nn = root.add("DelphesML", Setting::TypeGroup);
   nn.add("name", Setting::TypeString) = "ML"; // "Backpropogation:"+FANN_SIGMOID_SYMMETRIC;
   nn.add("events_for_training", Setting::TypeInt) = nBatch;
   nn.add("pileup_mu", Setting::TypeInt) = MuPileup;
   nn.add("number_of_kinematics", Setting::TypeInt) = num_kin;
   nn.add("number_of_epoch", Setting::TypeInt) = nEpoch;
   nn.add("number_of_layers", Setting::TypeInt) = num_layers;
   nn.add("number_of_inputs", Setting::TypeInt) = num_input;
   nn.add("number_of_outputs", Setting::TypeInt) = num_output;
   nn.add("number_of_hidden_nodes", Setting::TypeInt) = num_neurons_hidden_1;
   nn.add("number_of_threads", Setting::TypeInt) = num_threads;
   nn.add("number_of_inputs_eff", Setting::TypeInt) = num_input_eff;
   nn.add("number_of_outputs_eff", Setting::TypeInt) = num_output_eff;
   nn.add("dR_btag", Setting::TypeFloat) = dRbtag;
   nn.add("dR_isolation", Setting::TypeFloat) = dRisolation;
   nn.add("bquark_frac", Setting::TypeFloat) =  btag_frac;
 
   // configuration for jets
   Setting &njet = root.add("Jets", Setting::TypeGroup);
   njet.add("MinPT", Setting::TypeFloat) = minPT;
   njet.add("MaxEta", Setting::TypeFloat) = maxEta;
   njet.add("SlicesEtaPhi4Input", Setting::TypeInt) = slices_etaphi;
   njet.add("PtBinsForResolution", Setting::TypeInt) = nBinsNN;
   njet.add("EnergyBinsNr", Setting::TypeInt) = nBins;
   njet.add("EnergyScale", Setting::TypeFloat) = jet_escale;
   njet.add("EnergyShift", Setting::TypeFloat) = jet_eshift;
   njet.add("EtaScale", Setting::TypeFloat) = jet_etascale;
   njet.add("EtaShift", Setting::TypeFloat) = jet_etashift;
   njet.add("MassScale", Setting::TypeFloat) = jet_mscale;
   njet.add("MassShift", Setting::TypeFloat) = jet_mshift;
   // bins
   Setting &array = njet.add("PtBins", Setting::TypeArray);
   for(int i = 0; i <  nBins; ++i) array.add(Setting::TypeFloat) = eBins[i]; 


   // configuration for jets
   Setting &mu = root.add("EMObjects", Setting::TypeGroup);
   mu.add("MinPT", Setting::TypeFloat) = minPT;
   mu.add("MaxEta", Setting::TypeFloat) = maxEta;
   mu.add("SlicesEtaPhi4Input", Setting::TypeInt) = slices_etaphi;
   mu.add("PtBinsForResolution", Setting::TypeInt) = nBinsNN;
   mu.add("EnergyBinsNr", Setting::TypeInt) = nBins;
   mu.add("EnergyScale", Setting::TypeFloat) = em_escale;
   mu.add("EnergyShift", Setting::TypeFloat) = em_eshift;
   mu.add("EtaScale", Setting::TypeFloat) = em_etascale;
   mu.add("EtaShift", Setting::TypeFloat) = em_etashift;
   // bins
   Setting &array1 = mu.add("PtBins", Setting::TypeArray);
   for(int i = 0; i <  nBins; ++i) array1.add(Setting::TypeFloat) = eBins[i];





    
   try
  {
    cfg.writeFile(output_file);
    cerr << "New configuration successfully written to: " << output_file
         << endl;
  }
  catch(const FileIOException &fioex)
  {
    cerr << "I/O error while writing file: " << output_file << endl;
    return(EXIT_FAILURE);
  }




   double sumMSE=0;
   cout << "Saving ANN results " << endl;
   for (int m=0; m<nBins-1; m++){
   double mse1=fann_get_MSE(ann1_jets[m]);
   double mse2=fann_get_MSE(ann2_jets[m]);
   double mse3=fann_get_MSE(ann3_jets[m]);
   double mse4=fann_get_MSE(ann4_jets[m]);
   double mse=mse1+mse2+mse3+mse4;

   finalRME[m]=mse1+mse2+mse3+mse4;
   // cout << m << " Events=" << eventsBins[m]<< " MSE =" << mse << "  Improved =" << (initialRME[m]/mse) << "x. save to " << ann_name[m] << endl;
   if (mse<1000000 ) sumMSE=sumMSE+mse;

   // always save first sample
   if (firstTime[m] == true)  {
              fann_save(ann1_jets[m],  ann1_jets_name[m].c_str());
              fann_save(ann2_jets[m],  ann2_jets_name[m].c_str());
              fann_save(ann3_jets[m],  ann3_jets_name[m].c_str());
              fann_save(ann4_jets[m],  ann4_jets_name[m].c_str());
              fann_save(ann5_jets[m],  ann5_jets_name[m].c_str());
              // muons
              fann_save(ann1_muons[m],  ann1_muons_name[m].c_str());
              fann_save(ann2_muons[m],  ann2_muons_name[m].c_str());
              fann_save(ann3_muons[m],  ann3_muons_name[m].c_str());
              fann_save(ann4_muons[m],  ann4_muons_name[m].c_str());
              fann_save(ann5_muons[m],  ann5_muons_name[m].c_str());
              // electrons 
              fann_save(ann1_electrons[m],  ann1_electrons_name[m].c_str());
              fann_save(ann2_electrons[m],  ann2_electrons_name[m].c_str());
              fann_save(ann3_electrons[m],  ann3_electrons_name[m].c_str());
              fann_save(ann4_electrons[m],  ann4_electrons_name[m].c_str());
              fann_save(ann5_electrons[m],  ann5_electrons_name[m].c_str());
              // photons 
              fann_save(ann1_photons[m],  ann1_photons_name[m].c_str());
              fann_save(ann2_photons[m],  ann2_photons_name[m].c_str());
              fann_save(ann3_photons[m],  ann3_photons_name[m].c_str());
              fann_save(ann4_photons[m],  ann4_photons_name[m].c_str());
              fann_save(ann5_photons[m],  ann5_photons_name[m].c_str());


   };

   // for other samples save only small RMS NN
   if (firstTime[m] == false) { 
       if (mse>0 && mse < 0.1) {
              fann_save(ann1_jets[m],  ann1_jets_name[m].c_str());
              fann_save(ann2_jets[m],  ann2_jets_name[m].c_str());
              fann_save(ann3_jets[m],  ann3_jets_name[m].c_str());
              fann_save(ann4_jets[m],  ann4_jets_name[m].c_str());
              fann_save(ann5_jets[m],  ann5_jets_name[m].c_str());

              fann_save(ann1_muons[m],  ann1_muons_name[m].c_str());
              fann_save(ann2_muons[m],  ann2_muons_name[m].c_str());
              fann_save(ann3_muons[m],  ann3_muons_name[m].c_str());
              fann_save(ann4_muons[m],  ann4_muons_name[m].c_str());
              fann_save(ann5_muons[m],  ann5_muons_name[m].c_str());

              fann_save(ann1_electrons[m],  ann1_electrons_name[m].c_str());
              fann_save(ann2_electrons[m],  ann2_electrons_name[m].c_str());
              fann_save(ann3_electrons[m],  ann3_electrons_name[m].c_str());
              fann_save(ann4_electrons[m],  ann4_electrons_name[m].c_str());
              fann_save(ann5_electrons[m],  ann5_electrons_name[m].c_str());

              fann_save(ann1_photons[m],  ann1_photons_name[m].c_str());
              fann_save(ann2_photons[m],  ann2_photons_name[m].c_str());
              fann_save(ann3_photons[m],  ann3_photons_name[m].c_str());
              fann_save(ann4_photons[m],  ann4_photons_name[m].c_str());
              fann_save(ann5_photons[m],  ann5_photons_name[m].c_str());
 
        } else { 
             cout << m << " Do not save   with RME=" << mse << endl;
        }
    };

   fann_destroy(ann1_jets[m]);
   fann_destroy(ann2_jets[m]);
   fann_destroy(ann3_jets[m]);
   fann_destroy(ann4_jets[m]);
   fann_destroy(ann5_jets[m]);

   fann_destroy(ann1_muons[m]);
   fann_destroy(ann2_muons[m]);
   fann_destroy(ann3_muons[m]);
   fann_destroy(ann4_muons[m]);
   fann_destroy(ann5_muons[m]);

   fann_destroy(ann1_electrons[m]);
   fann_destroy(ann2_electrons[m]);
   fann_destroy(ann3_electrons[m]);
   fann_destroy(ann4_electrons[m]);
   fann_destroy(ann5_electrons[m]);

   fann_destroy(ann1_photons[m]);
   fann_destroy(ann2_photons[m]);
   fann_destroy(ann3_photons[m]);
   fann_destroy(ann4_photons[m]);
   fann_destroy(ann5_photons[m]);

    };

   cout << "Total MSE for all bins for this run is " << sumMSE << endl;



   return 0;

}


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
   nn.add("number_of_epoch", Setting::TypeInt) = nEpoch;
   nn.add("number_of_layers", Setting::TypeInt) = num_layers;
   nn.add("number_of_inputs", Setting::TypeInt) = num_input;
   nn.add("number_of_outputs", Setting::TypeInt) = num_output;
   nn.add("number_of_hidden_nodes", Setting::TypeInt) = num_neurons_hidden_1;
   nn.add("number_of_threads", Setting::TypeInt) = num_threads;
   nn.add("number_of_inputs_eff", Setting::TypeInt) = num_input_eff;
   nn.add("number_of_outputs_eff", Setting::TypeInt) = num_output_eff;

   // configuration for jets
   Setting &njet = root.add("Jets", Setting::TypeGroup);
   njet.add("MinPT", Setting::TypeFloat) = minPT;
   njet.add("MaxEta", Setting::TypeFloat) = maxEta;
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
   };

   // for other samples save only small RMS NN
   if (firstTime[m] == false) { 
       if (mse>0 && mse < 1) {
              fann_save(ann1_jets[m],  ann1_jets_name[m].c_str());
              fann_save(ann2_jets[m],  ann2_jets_name[m].c_str());
              fann_save(ann3_jets[m],  ann3_jets_name[m].c_str());
              fann_save(ann4_jets[m],  ann4_jets_name[m].c_str());
              fann_save(ann5_jets[m],  ann5_jets_name[m].c_str());
        } else { 
             cout << m << " Do not save   with RME=" << mse << endl;
        }
    };

   fann_destroy(ann1_jets[m]);
   fann_destroy(ann2_jets[m]);
   fann_destroy(ann3_jets[m]);
   fann_destroy(ann4_jets[m]);
   fann_destroy(ann5_jets[m]);

    };

   cout << "Total MSE for all bins for this run is " << sumMSE << endl;



   return 0;

}


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

/*
    cout << "Saving all ANN results " << endl;
    for (int m=0; m<nBins-1; m++){
   double mse=fann_get_MSE(ann[m]);
   finalRME[m]=mse;
   cout << m << ") MSE error=" << mse << "Improved =" << int(initialRME[m]/mse) << "x. ANN save to " << ann_name[m] << endl;
   fann_save(ann[m],  ann_name[m].c_str());
   fann_destroy(ann[m]);


    };
*/




   return 0;

}


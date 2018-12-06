// S.Chekanov

#include "Ana.h"

string getEnvVar( std::string const & key ) {
	char * val = getenv( key.c_str() );
	return val == NULL ? std::string("") : std::string(val);
}

// initialize your calculations (histograms)
Int_t Ana::Init() {


	// jet cuts
	minPT=25;
	maxEta=2.5;
	nevv=0;
	DeltaR=0.2; // cone to match with jets
	MSESTOP=0.0001; // error to stop training
	MuPileup=0;

	mcEventWeight=1.0;
	double xBins[] = {minPT,50,70,100,125,150,175,200,225,250,275,300,325,350,375,400};
	//double xBins[] = {minPT,50,100,200,400,600,800,1000,1200,1400,1600,1800,2000,2400,2800,3200};
	const int nxBins=sizeof(xBins)/sizeof(double);


	for (int m=0; m<nBins; m++){
		eBins[m]=xBins[m];
	};

	cout << " ###  Initialize ANN in Nr bins=" <<nBins-1<<  endl;


	for (int m=0; m<nBins-1; m++){
		firstTime[m]=true;
		initialRME[m]=0;
		finalRME[m]=0;
		eventsBins[m]=0; // Nr of events
		ann_jets_name[m]="out_ann/ann_jet_"+std::to_string(m)+".net";
		ann_jets_eff_name[m]="out_ann/ann_jet_eff_"+std::to_string(m)+".net";



		// if 1st ANN is found, we will read previous NN
		std::ifstream ifile(ann_jets_name[m].c_str());
		if ((bool)ifile) {
			firstTime[m]=false;
			cout << "We found ANN file for bin=" << m << " So we improve the ANN training.." << endl;
		}


		if (firstTime[m]) {
			cout << "   - Create new ANN Nr " << ann_jets_name[m] << endl;
			ann_jets[m] = fann_create_standard(num_layers, num_input, num_neurons_hidden_1, num_output);
			fann_set_activation_function_hidden(ann_jets[m], FANN_SIGMOID_SYMMETRIC);
			fann_set_activation_function_output(ann_jets[m], FANN_SIGMOID_SYMMETRIC);
			fann_randomize_weights(ann_jets[m],-1.0,1.0);
			//fann_set_activation_function_hidden(ann[m], FANN_SIGMOID);
			//fann_set_activation_function_output(ann[m], FANN_SIGMOID);
			//fann_randomize_weights(ann[m],0.0,0.9999);


			cout << "   - Create new ANN for efficiency.. " << ann_jets_eff_name[m] << endl;
			int num_layers_eff=3;
			int num_input_eff=3;
			int num_neurons_hidden_eff=3;
			int num_output_eff=1;
			ann_jets_eff[m] = fann_create_standard(num_layers_eff, num_input_eff, num_neurons_hidden_eff, num_output_eff);
			fann_set_activation_function_hidden(ann_jets_eff[m], FANN_SIGMOID_SYMMETRIC);
			fann_set_activation_function_output(ann_jets_eff[m], FANN_SIGMOID_SYMMETRIC);
			fann_randomize_weights(ann_jets_eff[m],-1.0,1.0);

		}


		if (firstTime[m]==false){
			cout << "Found previous ANN. Continue training.. Read from " << ann_jets_name[m] << endl;
			ann_jets[m] = fann_create_from_file(ann_jets_name[m].c_str());
			//fann_set_activation_function_hidden(ann[m], FANN_SIGMOID);
			//fann_set_activation_function_output(ann[m], FANN_SIGMOID);
			//for (int m=0; m<nBins-1; m++){
			//    double mse=fann_get_MSE(ann[m]);
			//    cout << m << ") initial MSE error=" << mse << endl;
			//}
			ann_jets_eff[m] = fann_create_from_file(ann_jets_eff_name[m].c_str());

		}


		if (m==nBins-2) {
			cout << "# Structure this ANN" << endl;
			fann_print_connections(ann_jets[m]);
			fann_print_parameters(ann_jets[m]);
		};



	}



	ffile="Analysis.root";
	cout << "\n -> Output file is =" << ffile << endl;
	RootFile = new TFile(ffile.c_str(), "RECREATE", "Histogram file");
	// define histograms
	h_debug = new TH1D("debug", "debug", 10, 0, 10);
	h_jetpt = new TH1D("jet_PT", "jet_PT reconstructed", 100, 0, 500);  h_jetpt->Sumw2();
	h_jetpt_truth = new TH1D("jet_truth_PT", "jet_PT truth", 100, 0, 500);  h_jetpt_truth->Sumw2();

	h_dR = new TH1D("jet_dR", "truth-jet distance", 500, 0, 7);


	h_in1 = new TH1D("in1", "in1", nBinsNN, -1.2, 1.2);
	h_in2 = new TH1D("in2", "in2", nBinsNN, -1.2, 1.2);
	h_in3 = new TH1D("in3", "in3", nBinsNN, -1.2, 1.2);
	h_in4 = new TH1D("in4", "in4", nBinsNN, -1.2, 1.2);

	h_out1 = new TH1D("out1", "out1", nBinsNN, -1.2, 1.2);
	h_out2 = new TH1D("out2", "out2", nBinsNN, -1.2, 1.2);
	h_out3 = new TH1D("out3", "out3", nBinsNN, -1.2, 1.2);
	h_out4 = new TH1D("out4", "out4", nBinsNN, -1.2, 1.2);
	h_out5 = new TH1D("out5_eff", "out5_eff", nBinsNN, -1.2, 1.2);



	// create ntuple
	m_ntuple  = new TTree("Ntuple","Ntuple");
	m_ntuple->Branch("RunNumber",   &RunNumber,  "RunNumber/I");
	m_ntuple->Branch("EventNumber", &EventNumber,   "EventNumber/I");
	//a_jets = new TClonesArray("LParticle",10);
	//m_ntuple->Branch("jets", "TClonesArray", &a_jets, 256000, 0);

	// read files from data.in file and put to a vector
	string name="data.in";
	ifstream myfile;
	myfile.open(name.c_str(), ios::in);
	if (!myfile) {
		cerr << "Can't open input file:  " << name << endl;
		exit(1);
	} else {
		cout << "-> Read data file=" << name << endl;
	}
	string temp;
	while (myfile >> temp) {
		//the following line trims white space from the beginning of the string
		temp.erase(temp.begin(), std::find_if(temp.begin(), temp.end(), not1(ptr_fun<int, int>(isspace))));
		if (temp.find("#") == 0) continue;

		// detect files with pileup (automatically)
		std::size_t found = temp.find("rfast006");
		if (found!=std::string::npos) MuPileup=140; // detect pileup events from HepSim

		ntup.push_back(temp);
	}
	cout << "-> Number of files=" << ntup.size()  << endl;
	myfile.close();
	for (unsigned int i=0; i<ntup.size(); i++) {
		cout << i << " file to analyse="+ntup[i] << endl;
	}




      // energy and eta scaling factors 
      jet_escale=1.5;
      jet_eshift=0.0;
      jet_etascale=2.0;
      jet_etashift=0.0;
      if (MuPileup>100) { // shrink to fit to -1 - 1 
                          jet_escale=0.5;
                          jet_eshift=-0.5;
                          jet_etascale=2.0;
                          jet_etashift=0.0;
                         }





	return 0;
}


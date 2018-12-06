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
	nev=0;
	DeltaR=0.2; // cone to match with jets
	MuPileup=0;


	static const char *conf_file = "out_ann/config.cfg";
	Config cfg;
	try
	{
		cout << "Read configuration=" << conf_file << endl;
		cfg.readFile(conf_file);
	}
	catch(const FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file: " <<  conf_file << std::endl;
		return(EXIT_FAILURE);
	}
	catch(const ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
		<< " - " << pex.getError() << std::endl;
		return(EXIT_FAILURE);
	}


	const Setting& root = cfg.getRoot();

	// ML
	Setting &mL = root["DelphesML"];
	if (!(mL.lookupValue("pileup_mu", MuPileup) ))  {
		cout << "Error: some values are not in DelphesML configuration!" << endl;
		exit(0);
	}

	Setting &njets = root["Jets"];
	if (!(njets.lookupValue("JetEnergyBinsNr", nBins)
	                && njets.lookupValue("JetMinPT", minPT)
	                && njets.lookupValue("JetMaxEta", maxEta)
	                && njets.lookupValue("JetEnergyBinsForResolution", nBinsNN))

	   ) {
		cout << "Error: some values are not in Jet configuration!" << endl;
		exit(0);
	}


	eBins = new double[nBinsNN];
	const Setting &bin_settings = njets.lookup("JetEnergyBins");
	for (int n = 0; n < bin_settings.getLength(); ++n) eBins[n]=bin_settings[n];
	if (bin_settings.getLength() != nBins) {
		cout << "getLength() =" << bin_settings.getLength() << "  nBinsNN=" << nBins << endl;
		cout << "Error: Bin sizes for jets cannot be read!" << endl;
		exit(0);
	}


	num_input=4;
	num_output= num_input*nBinsNN;

	mcEventWeight=1.0;


	// vector with bins for reco/true
	ann_jets_name = new string[nBinsNN-1];
	ann_jets_eff_name = new string[nBinsNN-1];

	initialRME = new double[nBinsNN-1];
	finalRME   = new double[nBinsNN-1];
	BinOverTrue = new int[nBinsNN];
	ann_jets = new fann*[nBinsNN-1];
	ann_jets_eff = new fann*[nBinsNN-1];
	for (unsigned int jjj=0; jjj<nBinsNN; jjj++) BinOverTrue[jjj]=jjj;

	//eBins = new double[nBins];
	for (int m=0; m<nBins; m++){
		cout << "Used jet bins=" << eBins[m] << endl;
	};

	cout << " ###  Initialize ANN in Nr bins=" <<nBins-1<<  endl;

	for (int m=0; m<nBins-1; m++){

		initialRME[m]=0;
		finalRME[m]=0;

		// jets and efficiencis
		ann_jets_name[m]="out_ann/ann_jet_"+std::to_string(m)+".net";
		ann_jets_eff_name[m]="out_ann/ann_jet_eff_"+std::to_string(m)+".net";

		//ann[m] = fann_create_standard(num_layers, num_input, num_neurons_hidden_1, num_output);
		cout << "Read ANN from " << ann_jets_name[m] << endl;
		ann_jets[m] = fann_create_from_file(ann_jets_name[m].c_str());
		ann_jets_eff[m] = fann_create_from_file(ann_jets_eff_name[m].c_str());
		//fann_set_activation_function_hidden(ann[m], FANN_SIGMOID);
		//fann_set_activation_function_output(ann[m], FANN_SIGMOID);
		//fann_set_activation_function_hidden(ann_jets[m], FANN_SIGMOID_SYMMETRIC);
		//fann_set_activation_function_output(ann_jets[m], FANN_SIGMOID_SYMMETRIC);

		if (m==nBins-2) {
			cout << "# Structure for last ANN" << endl;
			fann_print_connections(ann_jets[m]);
			fann_print_parameters(ann_jets[m]);
		};



	}



	ffile="Analysis.root";
	cout << "\n -> Output file is =" << ffile << endl;
	RootFile = new TFile(ffile.c_str(), "RECREATE", "Histogram file");
	// define histograms
	h_debug = new TH1D("debug", "debug", 10, 0, 10);
	h_jetpt = new TH1D("jet_PT", "jet_PT reconstructed", 100, 0, 1000);  h_jetpt->Sumw2();
	h_jetpt_truth = new TH1D("jet_truth_PT", "jet_PT truth", 100, 0, 1000);  h_jetpt_truth->Sumw2();
	h_jetpt_nn = new TH1D("jet_nn_PT", "jet_PT from NN", 100, 0, 1000);  h_jetpt_nn->Sumw2();
	h_jetcutflow = new TH1D("jet_cutflow", "jet cutflow", 10, 0, 10);

	h_jetres100 = new TH1D("jet_res50", "pT/pTtrue ratio", 200, 0.0, 3.0);
	h_jetres1000 = new TH1D("jet_res200", "pT/pTtrue ratio", 200, 0.0, 3.0);

	h_ptcor = new TH1D("jet_ptcor", "PT correction", 200, 0.0, 3);
	h_etacor = new TH1D("jet_etacor", "Eta correction", 200, 0.0, 3);
	h_phicor = new TH1D("jet_phicor", "Phi correction", 200, 0.0, 3);


	h_in1 = new TH1D("in1", "in1", 140, -1.2, 1.2);
	h_in2 = new TH1D("in2", "in2", 140, -1.2, 1.2);
	h_in3 = new TH1D("in3", "in3", 140, -1.2, 1.2);
	h_in4 = new TH1D("in4", "in4", 140, -1.2, 1.2);

	h_out1 = new TH1D("out1", "out1", nBinsNN, -1., 1.);
	h_out2 = new TH1D("out2", "out2", nBinsNN, -1., 1.);
	h_out3 = new TH1D("out3", "out3", nBinsNN, -1., 1.);
	h_out4 = new TH1D("out4", "out4", nBinsNN, -1., 1.);
	h_out5_eff = new TH1D("out5_eff", "out5 efficiency", 140, -1.2, 1.2);


	// create ntuple
	m_ntuple  = new TTree("Ntuple","Ntuple");
	m_ntuple->Branch("RunNumber",   &RunNumber,  "RunNumber/I");
	m_ntuple->Branch("EventNumber", &EventNumber,   "EventNumber/I");
	m_ntuple->Branch("AntiKt4JetPt",    &m_jetpt);
	m_ntuple->Branch("AntiKt4JetEta",   &m_jeteta);
	m_ntuple->Branch("AntiKt4JetPhi",   &m_jetphi);


	m_ntuple->Branch("AntiKt4TruthJetPT",    &m_gjetpt);
	m_ntuple->Branch("AntiKt4TruthJetEta",   &m_gjeteta);
	m_ntuple->Branch("AntiKt4TruthJetPhi",   &m_gjetphi);

	m_ntuple->Branch("AntiKt4NNJetPt",    &m_nnjetpt);
	m_ntuple->Branch("AntiKt4NNJetEta",   &m_nnjeteta);
	m_ntuple->Branch("AntiKt4NNJetPhi",   &m_nnjetphi);



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
		ntup.push_back(temp);
	}
	cout << "-> Number of files=" << ntup.size()  << endl;
	myfile.close();
	for (unsigned int i=0; i<ntup.size(); i++) {
		cout << i << " file to analyse="+ntup[i] << endl;
	}


	return 0;
}


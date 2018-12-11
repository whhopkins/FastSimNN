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
        jet_escale=0;
        jet_mscale=0;
        jet_etascale=0;

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
	if (!(mL.lookupValue("pileup_mu", MuPileup)
                 && mL.lookupValue("number_of_inputs",   num_input)
                 && mL.lookupValue("number_of_inputs_eff",  num_input_eff) 
                 && mL.lookupValue("number_of_outputs_eff",  num_output_eff) 
                 ))  {
		cout << "Error: some values are not in DelphesML configuration!" << endl;
		exit(0);
	}

	Setting &njets = root["Jets"];
	if (!(njets.lookupValue("EnergyBinsNr", nBins)
	                && njets.lookupValue("MinPT", minPT)
	                && njets.lookupValue("MaxEta", maxEta)
	                && njets.lookupValue("PtBinsForResolution", nBinsNN)
                        && njets.lookupValue("EnergyScale",jet_escale)
                        && njets.lookupValue("EnergyShift",jet_eshift)
                        && njets.lookupValue("MassScale",jet_mscale)
                        && njets.lookupValue("MassShift",jet_mshift)
                        && njets.lookupValue("EtaScale",jet_etascale)
                        && njets.lookupValue("EtaShift",jet_etashift))

	   ) {
		cout << "Error: some values are not in Jet configuration!" << endl;
		exit(0);
	}


        if (jet_escale == 0 || jet_mscale == 0 || jet_etascale == 0){
         cout << "Error: failed to read scale factors from configuration file !" << endl;
         exit(0);
        }

	eBins = new double[nBinsNN];
	const Setting &bin_settings = njets.lookup("PtBins");
	for (unsigned int n = 0; n < bin_settings.getLength(); ++n) eBins[n]=bin_settings[n];
	if (bin_settings.getLength() != (unsigned int)nBins) {
		cout << "getLength() =" << bin_settings.getLength() << "  nBinsNN=" << nBins << endl;
		cout << "Error: Bin sizes for jets cannot be read!" << endl;
		exit(0);
	}


        // number of bins -1 
	num_output= nBinsNN-1;

	mcEventWeight=1.0;

	// vector with bins for reco/true
	ann1_jets_name = new string[nBinsNN-1];
        ann2_jets_name = new string[nBinsNN-1];
        ann3_jets_name = new string[nBinsNN-1];
        ann4_jets_name = new string[nBinsNN-1];
        ann5_jets_name = new string[nBinsNN-1];

        ann1_jets = new fann*[nBinsNN-1];
        ann2_jets = new fann*[nBinsNN-1];
        ann3_jets = new fann*[nBinsNN-1];
        ann4_jets = new fann*[nBinsNN-1];
        ann5_jets = new fann*[nBinsNN-1];

	initialRME = new double[nBinsNN-1];
	finalRME   = new double[nBinsNN-1];
	BinOverTrue = new int[nBinsNN-1];

	for (int jjj=0; jjj<nBinsNN-1; jjj++) BinOverTrue[jjj]=jjj;

	//eBins = new double[nBins];
	for (int m=0; m<nBins; m++){
		cout << "Used jet bins=" << eBins[m] << endl;
	};

	cout << " ###  Initialize ANN in Nr bins=" <<nBins-1<<  endl;

	for (int m=0; m<nBins-1; m++){

		initialRME[m]=0;
		finalRME[m]=0;

		// jets and efficiencis
		ann1_jets_name[m]="out_ann/ann1_jet_"+std::to_string(m)+".net";
                ann2_jets_name[m]="out_ann/ann2_jet_"+std::to_string(m)+".net";
                ann3_jets_name[m]="out_ann/ann3_jet_"+std::to_string(m)+".net";
                ann4_jets_name[m]="out_ann/ann4_jet_"+std::to_string(m)+".net";
                ann5_jets_name[m]="out_ann/ann5_jet_"+std::to_string(m)+".net";

		cout << "Read ANN from files.."  << endl;
		ann1_jets[m] = fann_create_from_file(ann1_jets_name[m].c_str());
                ann2_jets[m] = fann_create_from_file(ann2_jets_name[m].c_str());
                ann3_jets[m] = fann_create_from_file(ann3_jets_name[m].c_str());
                ann4_jets[m] = fann_create_from_file(ann4_jets_name[m].c_str());
                ann5_jets[m] = fann_create_from_file(ann5_jets_name[m].c_str());
  

		if (m==nBins-2) {
			cout << "# Structure for last ANN" << endl;
			fann_print_connections(ann1_jets[m]);
			fann_print_parameters(ann1_jets[m]);
		};



	}



	ffile="Analysis.root";
	cout << "\n -> Output file is =" << ffile << endl;
	RootFile = new TFile(ffile.c_str(), "RECREATE", "Histogram file");
	// define histograms
	h_debug = new TH1D("debug", "debug", 10, 0, 10);
        h_dR = new TH1D("jet_dR", "truth-jet distance", 500, 0, 7);

	h_jetpt = new TH1D("jet_PT", "jet_PT reconstructed", 100, 0, 1000);  h_jetpt->Sumw2();
	h_jetpt_truth = new TH1D("jet_truth_PT", "jet_PT truth", 100, 0, 1000);  h_jetpt_truth->Sumw2();
	h_jetpt_nn = new TH1D("jet_nn_PT", "jet_PT from NN", 100, 0, 1000);  h_jetpt_nn->Sumw2();
	h_jetcutflow = new TH1D("jet_cutflow", "jet cutflow", 10, 0, 10);

	h_jetres100 = new TH1D("jet_res50", "pT/pTtrue ratio", 200, 0.0, 3.0);
	h_jetres1000 = new TH1D("jet_res200", "pT/pTtrue ratio", 200, 0.0, 3.0);

	h_ptcor = new TH1D("jet_ptcor", "PT correction", 200, 0.0, 3);
	h_etacor = new TH1D("jet_etacor", "Eta correction", 200, 0.0, 3);
	h_phicor = new TH1D("jet_phicor", "Phi correction", 200, 0.0, 3);


	h_in1 = new TH1D("in1", "in1", nBinsNN, -1., 1.);
	h_in2 = new TH1D("in2", "in2", nBinsNN, -1., 1.);
	h_in3 = new TH1D("in3", "in3", nBinsNN, -1., 1.);
	h_in4 = new TH1D("in4", "in4", nBinsNN, -1., 1.);

	h_out1 = new TH1D("out1", "out1", nBinsNN, -1., 1.);
	h_out2 = new TH1D("out2", "out2", nBinsNN, -1., 1.);
	h_out3 = new TH1D("out3", "out3", nBinsNN, -1., 1.);
	h_out4 = new TH1D("out4", "out4", nBinsNN, -1., 1.);
	h_out5_eff = new TH1D("out5_eff", "out5 efficiency", nBinsNN, -1., 1.);
        h_out6_btag = new TH1D("out6_btag", "out6 b-tagging", nBinsNN, -1., 1.);

        h_rout1 = new TH1D("rout1", "out1 random bin", nBinsNN, 0, nBinsNN);
        h_rout2 = new TH1D("rout2", "out2 random bin", nBinsNN, 0, nBinsNN);
        h_rout3 = new TH1D("rout3", "out3 random bin", nBinsNN, 0, nBinsNN);
        h_rout4 = new TH1D("rout4", "out4 random bin", nBinsNN, 0, nBinsNN);


	// create ntuple
	m_ntuple  = new TTree("Ntuple","Ntuple");
	m_ntuple->Branch("RunNumber",   &RunNumber,  "RunNumber/I");
	m_ntuple->Branch("EventNumber", &EventNumber,   "EventNumber/I");
	m_ntuple->Branch("AntiKt4JetPt",    &m_jetpt);
	m_ntuple->Branch("AntiKt4JetEta",   &m_jeteta);
	m_ntuple->Branch("AntiKt4JetPhi",   &m_jetphi);
        m_ntuple->Branch("AntiKt4JetM",      &m_jetm);
        m_ntuple->Branch("AntiKt4JetBtag",   &m_jetbtag);
        // matched
        m_ntuple->Branch("AntiKt4MatchedJetPt",    &m_matchedjetpt);
        m_ntuple->Branch("AntiKt4MatchedJetEta",   &m_matchedjeteta);
        m_ntuple->Branch("AntiKt4MatchedJetPhi",   &m_matchedjetphi);
        m_ntuple->Branch("AntiKt4MatchedJetM",     &m_matchedjetm);
        m_ntuple->Branch("AntiKt4MatchedJetBtag",  &m_matchedjetbtag);
        // truth
	m_ntuple->Branch("AntiKt4TruthJetPt",    &m_gjetpt);
	m_ntuple->Branch("AntiKt4TruthJetEta",   &m_gjeteta);
	m_ntuple->Branch("AntiKt4TruthJetPhi",   &m_gjetphi);
        m_ntuple->Branch("AntiKt4TruthJetM",   &m_gjetm);
        m_ntuple->Branch("AntiKt4TruthJetBtag",   &m_gjetbtag);
        // NN jets
	m_ntuple->Branch("AntiKt4NNJetPt",    &m_nnjetpt);
	m_ntuple->Branch("AntiKt4NNJetEta",   &m_nnjeteta);
	m_ntuple->Branch("AntiKt4NNJetPhi",   &m_nnjetphi);
        m_ntuple->Branch("AntiKt4NNJetM",   &m_nnjetm);
        m_ntuple->Branch("AntiKt4NNJetBtag",   &m_nnjetbtag);


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


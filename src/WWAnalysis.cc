#include "WWAnalysis.h"

WWAnalysis aWWAnalysis;


WWAnalysis::WWAnalysis() : Processor("WWAnalysis") {


  // register steering parameters: name, description, class-variable, default value

	registerProcessorParameter( "Printing" ,
	                            "Print certain messages"  ,
	                             _printing,
	                             (int)5 ) ;

	std::string inputMcParticleCollectionName = "x";
	registerInputCollection( LCIO::MCPARTICLE,
				"McParticleCollectionName" ,
				"Name of the MCParticle input collection" ,
				_inputMcParticleCollectionName,
				inputMcParticleCollectionName);

	//input collection parameters
	std::string inputJetCollectionName = "x";
  	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     	"InputJetCollectionName" , 
			     	"Input Jet Collection Name "  ,
			     	_inputJetCollectionName,
			      	inputJetCollectionName);

}

void WWAnalysis::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl;
  // usually a good idea to
  printParameters() ;
  nEvt = 0;
	
	file = new TFile("file.root","RECREATE");

	/* init histograms */
	WmassMuon = new TH1D("Wmassmuon","W mass, mass of dijet or lepton jet from muon event",100, 50.0, 120.0 );
	WmassTau = new TH1D("Wmasstau","W mass, mass of dijet or lepton jet from tau event",100, 50.0, 120.0 );
	WEMuon = new TH1D("WEmuon","W Energy, energy of dijet or lepton jet from muon event",100, 25.0, 300.0);
	WETau = new TH1D("WEtau","W Energy,energy of dijet or lepton jet from tau event",100, 50.0, 250.0 );
	//TH1D* Wm_cosTheta;
	LjetMassMuon=new TH1D("Ljetmassmuon","ljet mass, mass of lepton jet from muon event",100, 0.0 ,3.0 );
	LjetMassTau=new TH1D("Ljetmasstau","ljet mass, mass of lepton jet from tau event",100, 0.0, 5.0 );

	costhetawMuon = new TH1D("costhetawMuon", "production angle of W- in Muon event",100,-1.0,1.0);
	thetaLMuon = new TH1D("thetaLMuon", "polar angle of CM lepton in Muon event",100, 0.0, 3.14);
	phiLMuon = new TH1D("phiLMuon", "azimuthal angle of CM Lepton in Muon event", 100,-3.14,3.14);
	thetaHMuon = new TH1D("thetaHMuon", "polar angle of CM quark in Muon event",100,0.0,3.14);
	phiHMuon = new TH1D("phiHMuon","azimuthal angle of CM quark in Muon event", 100,-3.14,3.14);

	costhetawTau = new TH1D("costhetawTau", "production angle of W- in Tau event",100,-1.0,1.0);
	thetaLTau = new TH1D("thetaLTau", "polar angle of CM lepton in Tau event",100, 0.0, 3.14);
	phiLTau = new TH1D("phiLTau", "azimuthal angle of CM Lepton in Tau event", 100,-3.14,3.14);
	thetaHTau = new TH1D("thetaHTau", "polar angle of CM quark in Tau event",100,0.0,3.14);
	phiHTau = new TH1D("phiHTau","azimuthal angle of CM quark in Tau event", 100,-3.14,3.14);
}

void WWAnalysis::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ;
}

bool WWAnalysis::FindMCParticles( LCEvent* evt ){
   
	bool collectionFound = false;

  	// clear old global MCParticle vector
  	_mcpartvec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;

	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    
		//if found print name and number of elements 
		if(*itname==_inputMcParticleCollectionName){
      			LCCollection* collection = evt->getCollection(*itname);
     			std::cout<< "Located MC Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
      			collectionFound = true;
      
			//add the collection elements to the global vector
			for(int i=0;i<collection->getNumberOfElements();i++){
				MCParticle* mcPart = dynamic_cast<MCParticle*>(collection->getElementAt(i));
				_mcpartvec.push_back(mcPart);

       
      			}
    		}
  	}

  	if(!collectionFound){
		std::cout<<"MC Collection "<< _inputMcParticleCollectionName << "not found"<<std::endl;
	}
  
  	return collectionFound;
}

bool WWAnalysis::FindJets( LCEvent* evt ) {

	bool collectionFound = false;

  	// clear old global pfovector
	_jets.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;
	
	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
     
		//if found print name and number of elements
    		if(*itname==_inputJetCollectionName){ 
			LCCollection* collection = evt->getCollection(*itname);
			std::cout<< "Located Jets Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
			collectionFound = true;

 			//add the collection elements to the global vector
      			for(int i=0; i<collection->getNumberOfElements(); i++){
				ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(collection->getElementAt(i));
				_jets.push_back(recoPart);
      			}
    		}
  	}
	
	if(!collectionFound){
		std::cout<<"Jet Collection "<< _inputJetCollectionName << "not found"<<std::endl;
	}

   
	return collectionFound;
}
int WWAnalysis::identifyLeptonJet( std::vector<ReconstructedParticle*> jets){

	//maybe the lepton is the jet with the least particles
	int indexofminjet = -1;
	int nparticles = 999999;
	for(int i=0; i<_jets.size(); i++){
		std::cout<<"jet "<<i<<" has "<< _jets.at(i)->getParticles().size() << " particles "<<std::endl;
		if( _jets.at(i)->getParticles().size() < nparticles ){
			indexofminjet = i;
			nparticles = _jets.at(i)->getParticles().size();
		}
	}
	
	return indexofminjet;

}
int WWAnalysis::getLeptonJetCharge( ReconstructedParticle* ljet ){
	//assign by leading track charge? or charge sum of reco parts?
	std::vector<ReconstructedParticle*> jetparts = ljet->getParticles();
	int totalcharge = 0;
	int leadingcharge = 0;
	double maxP = 0;
	for(int i=0; i<jetparts.size(); i++){
		//method 1
		totalcharge += jetparts.at(i)->getCharge();

		//method 2
		const double* p = jetparts.at(i)->getMomentum();
		double P = std::sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
		if(P > maxP){
			maxP = P;
			leadingcharge = jetparts.at(i)->getCharge(); 
		}
	}

	//return totalcharge;
	return leadingcharge;

}

void WWAnalysis::processEvent( LCEvent * evt ) {
 FindMCParticles(evt);
 FindJets(evt);
 std::cout << "======================================== event " << nEvt << std::endl ;

	//bools to characterize the lepton decay for this event
	bool isTau = false;
	bool isMuon = false;
	//the true lepton charge
	int trueq;

	for(int i=0; i<_mcpartvec.size(); i++){
		std::vector<int> parentpdgs{};
		std::vector<int> daughterpdgs{};
		std::vector<MCParticle*> mcparents{};
		std::vector<MCParticle*> daughters{};
		daughters = _mcpartvec.at(i)->getDaughters();
		for(int j = 0; j<daughters.size(); j++){
			daughterpdgs.push_back(daughters.at(j)->getPDG());
			
		}
		
	

		//allowed quarks
		std::vector<int> quarks{ 1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6};
		std::vector<int> leptons{11, 12, 13, 14, 15, 16, -11, -12, -13, -14, -15, -16};
		//we require exactly 2 elements from leptons and 2 from quarks
		int lep=0;
		int qrk=0;

		//categorize the event for plotting
		
	
			for(int k=0; k<quarks.size(); k++){
				qrk += std::count(daughterpdgs.begin(),daughterpdgs.end(),quarks.at(k));
			}
			for(int k=0; k<leptons.size(); k++){
				lep += std::count(daughterpdgs.begin(),daughterpdgs.end(),leptons.at(k));
			}

		if( qrk == 2 && lep == 2){
		//found the proper set 
		for(int j=0; j<parentpdgs.size(); j++){
			std::cout<<parentpdgs.at(j)<<" ";
		}
		std::cout<< " -> "<<_mcpartvec.at(i)->getPDG()<<" -> ";
		for(int j=0; j<daughters.size(); j++){
			std::cout<<daughters.at(j)->getPDG()<<" ";
		}
		std::cout<<std::endl;

	

 		 if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 11) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -11) != daughterpdgs.end() ){
			nelec++;
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 13) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -13) != daughterpdgs.end() ){
			nmuon++;
			//identify event containing muon
			isMuon = true;
			//get true charge of the muon
			if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 13) != daughterpdgs.end() ){
				trueq = -1;
			}
			else{
				trueq = 1;
			}
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 15) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -15) != daughterpdgs.end() ){
			ntau++;
			//identify event containing a tau
			isTau = true;
			//identify the true charge of the lepton
			if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 15) != daughterpdgs.end() ){
				trueq = -1;
			}
			else{
				trueq = 1;
			}
				
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 1) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -1) != daughterpdgs.end() ){
			ndwn++;
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 2) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -2) != daughterpdgs.end() ){
			nup++;
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 3) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -3) != daughterpdgs.end() ){
			nstr++;
		}
		if (std::find(daughterpdgs.begin(),daughterpdgs.end(), 4) != daughterpdgs.end() ||
			std::find(daughterpdgs.begin(),daughterpdgs.end(), -4) != daughterpdgs.end() ){
			nchm++;
		}
		 break;

		}//end if 2 qrks and 2 lep

	}//end mcpartvec loop

	//now asses jets
	int ljet_index = identifyLeptonJet( _jets );
	int lq = getLeptonJetCharge( _jets.at(ljet_index) );

	if( trueq == lq){
		std::cout<<" got correct lepton charge "<<std::endl;
		if(isTau) tauqmatch++;
		if(isMuon) muonqmatch++;
	}
	else{ 
		std::cout<<" charge wrong "<<std::endl;
	}

	//put jets into tlvs
	std::vector<TLorentzVector> jets{};
	for(int i=0; i<_jets.size(); i++){
		TLorentzVector j;
		jets.push_back(j);
		jets.at(i).SetXYZM(_jets.at(i)->getMomentum()[0], _jets.at(i)->getMomentum()[1], _jets.at(i)->getMomentum()[2], _jets.at(i)->getMass() );
	}

	TLorentzVector dijet;
	TLorentzVector ljet;

	for(int i=0; i<jets.size(); i++){
		if( i == ljet_index ){
			ljet = jets.at(i);
		}
		else{
			dijet += jets.at(i);
		}
	}
	

	//figure out the muon 
	double missingPx= ljet.Px() - dijet.Px();
	double missingPy= ljet.Py() - dijet.Py();
	double missingPz= ljet.Pz() - dijet.Pz();

	TLorentzVector neutrino;
	neutrino.SetXYZM(missingPx, missingPy, missingPz, 0.0);

	TLorentzVector Wh= dijet;
	TLorentzVector Wl = neutrino + ljet;

	TVector3 Whboost(Wh.Px(),Wh.Py(),Wh.Pz());
	TVector3 Wlboost(Wl.Px(),Wl.Py(),Wl.Pz());

	Whboost = -Whboost;
	Wlboost = -Wlboost;

	std::vector<TLorentzVector> CMjets{};
	for(int i=0; i<jets.size(); i++){
		CMjets.push_back(jets.at(i));
		if(i == ljet_index){
			CMjets.at(i).Boost(wlboost);
		}
		else{
			CMjets.at(i).Boost(whboost);
		}
	}

	if( isTau ){
		WmassTau->Fill( dijet.M() );
		WmassTau->Fill(ljet.M() );
		WETau->Fill(dijet.E() );
		WETau->Fill(ljet.E() );
		LjetMassTau->Fill( ljet.M() );

		//+ events, take prod angle from qq
		if(lq> 0){
			
		}else{

		}
	}
	if( isMuon) {
		WmassMuon->Fill( dijet.M() );
		WmassMuon->Fill(ljet.M() );
		WEMuon->Fill(dijet.E() );
		WEMuon->Fill(ljet.E() );

		LjetMassMuon->Fill( ljet.M() );
	
	}



 nEvt++;
}
void WWAnalysis::end(){

	std::cout<<" nelec "<<nelec<<" nmuon "<< nmuon <<" ntau "<< ntau << std::endl;
	std::cout<<" ndwn "<<ndwn<<" nup "<<nup<<" nstr "<<nstr<<" nchm "<<nchm<<std::endl;

	std::cout<<" nevents "<< nEvt << " mu q match "<< muonqmatch <<  " tau q match "<< tauqmatch <<std::endl;

	file->Write();
}

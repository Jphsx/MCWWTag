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
 
	for(int i=0; i<= ncuts; i++){
		 char cuts[100];
         sprintf(cuts, "_%d", i);
         std::string cutnum(cuts);

		/* init histograms */
		WmassMuon[i] = new TH1D(("Wmassmuon"+cutnum).c_str(),"W^{#pm} Mass, with true #mu",100, 50.0, 120.0 );
		WmassTau[i] = new TH1D(("Wmasstau"+cutnum).c_str(),"W^{#pm} Mass, with true #tau",100, 50.0, 120.0 );
		WEMuon[i] = new TH1D(("WEmuon"+cutnum).c_str(),"W^{#pm} Energy, with true #mu",100, 25.0, 300.0);
		WETau[i] = new TH1D(("WEtau"+cutnum).c_str(),"W^{#pm} Energy, with true #tau ",100, 50.0, 250.0 );
		EtotalMuon[i] = new TH1D(("EtotalMuon"+cutnum).c_str(),"WW Total Energy, with true #mu",100,10,300); 
		EtotalTau[i] = new TH1D(("EtotalTau"+cutnum).c_str(),"WW Total Energy, with true #tau",100,10,300);
		//TH1D* Wm_cosTheta;
		LjetMassMuon[i]=new TH1D(("Ljetmassmuon"+cutnum).c_str(),"Mass of Lepton Jet with true #mu",100, 0.0 ,3.0 );
		LjetMassTau[i]=new TH1D(("Ljetmasstau"+cutnum).c_str(),"Mass of Lepton Jet with ture #tau",100, 0.0, 5.0 );

		costhetawMuon[i] = new TH1D(("costhetawMuon"+cutnum).c_str(), "production angle of W^- in Muon event",100,-1.0,1.0);
		thetaLMuon[i] = new TH1D(("thetaLMuon"+cutnum).c_str(), "polar angle of CM lepton in Muon event",100, 0.0, 3.14);
		phiLMuon[i] = new TH1D(("phiLMuon"+cutnum).c_str(), "azimuthal angle of CM Lepton in Muon event", 100,-3.14,3.14);
		thetaHMuon[i] = new TH1D(("thetaHMuon"+cutnum).c_str(), "polar angle of CM quark in Muon event",100,0.0,3.14);
		phiHMuon[i] = new TH1D(("phiHMuon"+cutnum).c_str(),"azimuthal angle of CM quark in Muon event", 100,-3.14,3.14);

		costhetawTau[i] = new TH1D(("costhetawTau"+cutnum).c_str(), "production angle of W^- in Tau event",100,-1.0,1.0);
		thetaLTau[i] = new TH1D(("thetaLTau"+cutnum).c_str(), "polar angle of CM lepton in Tau event",100, 0.0, 3.14);
		phiLTau[i] = new TH1D(("phiLTau"+cutnum).c_str(), "azimuthal angle of CM Lepton in Tau event", 100,-3.14,3.14);
		thetaHTau[i] = new TH1D(("thetaHTau"+cutnum).c_str(), "polar angle of CM quark in Tau event",100,0.0,3.14);
		phiHTau[i] = new TH1D(("phiHTau"+cutnum).c_str(),"azimuthal angle of CM quark in Tau event", 100,-3.14,3.14);
	
		/* end init histograms */
	}
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
/* identifies the lepton jet with the minimum particle multiplicity */
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
/* identifies Lepton jet charge by taking the charge of the leading track from the jet */
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
/* classify the the event based on the type of lepton in MCParticle info, also set the true charge for that lepton */
/* also tallies the number of muon/electron/tau events */
/* also retrieves the mcparticle which has daughters qqlnu */
MCParticle* WWAnalysis::classifyEvent(bool& isTau, bool& isMuon, int& trueq){
	
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
		std::vector<int> quarks{ 1, 2, 3, 4, -1, -2, -3, -4};
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
			//if we have found the true decay set break out of the mcpart vec loop
			break;
		}//end if 2qrk & 2 lep

	}//end mcpartvec loop

}
/* populate the tlvs based on the identified lepton jet */
void WWAnalysis::populateTLVs(int lindex){

	for(int i=0; i<_jets.size(); i++){

		TLorentzVector* j = new TLorentzVector();
		j->SetXYZM(_jets.at(i)->getMomentum()[0], _jets.at(i)->getMomentum()[1], _jets.at(i)->getMomentum()[2], _jets.at(i)->getMass() );
		jets.push_back(j);

		std::cout<<_jets.at(i)->getMomentum()[0]<<" "<< _jets.at(i)->getMomentum()[1]<<" "<<_jets.at(i)->getMomentum()[2]<< " "<< _jets.at(i)->getMass()<<std::endl;
	}
	
	
	//Wl = new TLorentzVector();
	Wqq = new TLorentzVector();
	TLorentzVector temp1;
	//loop over the new tlv jets and make wl and wqq
	for(int i=0; i<jets.size(); i++){
		if( i == lindex ){
			//right not Wl will be missing its neutrino
			Wl = jets.at(i);
		}
		else{
			temp1 += *jets.at(i);
		std::cout<<temp1.Px()<<" "<<temp1.Py()<<" "<<temp1.Pz()<<" "<<temp1.M()<<std::endl;
		}
	}
	Wqq->SetXYZM(temp1.Px(), temp1.Py(), temp1.Pz(), temp1.M() );
	
	

	//figure out the muon 
	double missingPx= Wl->Px() - Wqq->Px();
	double missingPy= Wl->Py() - Wqq->Py();
	double missingPz= Wl->Pz() - Wqq->Pz();

	std::cout<<"missing P "<< missingPx<<" "<<missingPy<<" "<<missingPz<<std::endl;
	//create the tlv neutrino
	nu = new TLorentzVector();
	nu->SetXYZM(missingPx, missingPy, missingPz, 0.0);

	//add the neutrino to complete the leptonic W
	TLorentzVector temp2;
	temp2 = *Wl + *nu;
	Wl->SetXYZM(temp2.Px(), temp2.Py(), temp2.Pz(), temp2.M());

	std::cout<<"WL and wqq at fn scope ";
	std::cout<<Wqq->Px()<<" "<<Wqq->Py()<<" "<<Wqq->Pz()<<" "<<Wqq->M()<<std::endl;

}
//populate W rest fram versions of the jets to access TGC observables
void WWAnalysis::populateCMTLVs(){
	TVector3 Wqqboost(Wqq->Px(),Wqq->Py(),Wqq->Pz());
	TVector3 Wlboost(Wl->Px(),Wl->Py(),Wl->Pz());

	Wqqboost = -Wqqboost;
	Wlboost = -Wlboost;

	for(int i=0; i<jets.size(); i++){
		CMJets.push_back(jets.at(i));
		if(i == ljet_index){
			CMJets.at(i)->Boost(Wlboost);
		}
		else{
			CMJets.at(i)->Boost(Wqqboost);
		}
	}
	//boost the neutrino into CM
	CMnu = nu;
	CMnu->Boost(Wlboost);
	
}
//get the production angle for W-  (W- . z)
double WWAnalysis::getCosThetaW(){
	
	//our unit z vector along the beam axis
	TVector3 z(0.0,0.0,1.0);
	if(lq <0 ){
		//W- is the lepton
		TVector3 Wm(Wl->Px(),Wl->Py(),Wl->Pz());
		Wm = Wm * (1/Wm.Mag());
		return Wm.Dot(z);
		
	}
	else{
		//infer qq charge to be W-
		TVector3 Wm(Wqq->Px(),Wqq->Py(),Wqq->Pz());
		Wm = Wm * (1/Wm.Mag());
		return Wm.Dot(z);
	}

}
void WWAnalysis::FillHistos(int histNumber){
	if(isTau){
		FillTauHistos(histNumber);
	}
	if(isMuon){
		FillMuonHistos(histNumber);
	}
}
void WWAnalysis::FillMuonHistos(int histNumber){

	WmassMuon[histNumber]->Fill( Wqq->M() );
	//WmassMuon[histNumber]->Fill(Wl.M() );
	WEMuon[histNumber]->Fill(Wqq->E() );
	//WEMuon[histNumber]->Fill(Wl.E() );

	LjetMassMuon[histNumber]->Fill( jets.at(ljet_index )->M() );

	//TGC stuff
/*	//costhetawMuon[histNumber]->Fill(getCosThetaW());
	thetaLMuon[histNumber]->Fill( CMJets.at(ljet_index).Theta());
	phiLMuon[histNumber]->Fill( CMJets.at(ljet_index).Phi());
	for(int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHMuon[histNumber]->Fill( CMJets.at(i).Theta()); 
			phiHMuon[histNumber]->Fill( CMJets.at(i).Phi());
		}
	} */
		
}
void WWAnalysis::FillTauHistos(int histNumber){

	WmassTau[histNumber]->Fill( Wqq->M() );
	//WmassTau[histNumber]->Fill( Wl.M() );
	WETau[histNumber]->Fill( Wqq->E() );
	//WETau[histNumber]->Fill(Wl.E() );

	LjetMassTau[histNumber]->Fill( jets.at(ljet_index)->M() );

	//TGC stuff
/*	//costhetawTau[histNumber]->Fill(getCosThetaW());
	thetaLTau[histNumber]->Fill( CMJets.at(ljet_index).Theta());
	phiLTau[histNumber]->Fill( CMJets.at(ljet_index).Phi());
	for(int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHTau[histNumber]->Fill( CMJets.at(i).Theta()); 
			phiHTau[histNumber]->Fill( CMJets.at(i).Phi());
		}
	} */
}

void WWAnalysis::processEvent( LCEvent * evt ) {

 FindMCParticles(evt);
 FindJets(evt);
 std::cout << "======================================== event " << nEvt << std::endl ;


	//bools to characterize the true lepton decay for this event
	isTau = false;
	isMuon = false;

	//from the mcparticles find what the type of semileptonic decay is present
	parent = classifyEvent(isTau, isMuon, trueq);

	//now asses jets
	//keep the index on _jets of the jet we consider to be the lepton
	ljet_index = identifyLeptonJet( _jets );

	//get the charge of the lepton jet
	lq = getLeptonJetCharge( _jets.at(ljet_index) );


	//check if the assessed charge matches the true charge of the lepton
	if( trueq == lq){
		std::cout<<" got correct lepton charge "<<std::endl;
		//count the number of times we get it right
		if(isTau) tauqmatch++;
		if(isMuon) muonqmatch++;
	}
	else{ 
		std::cout<<" charge wrong "<<std::endl;
	}

	//build up all the different tlvs for calculation
  	populateTLVs(ljet_index);
	std::cout<<"Wl at process scope"<<std::endl;
	std::cout<<Wqq->Px()<<" "<<Wqq->Py()<<" "<<Wqq->Pz()<<" "<<Wqq->M()<<std::endl;

    //boost jets to cm for TGC observables
	populateCMTLVs();

	//fill base histograms and produce histos with sequential cuts hist0 is always no cuts
	FillHistos(0);
	//cut #1 trueq == lq, lepton sign is correctly assessed
	if(trueq == lq){
		FillHistos(1);
	}



 nEvt++;
}
void WWAnalysis::end(){

	/* print out stuff */
	std::cout<<" nelec "<<nelec<<" nmuon "<< nmuon <<" ntau "<< ntau << std::endl;
	std::cout<<" nevents "<< nEvt << " mu q match "<< muonqmatch <<  " tau q match "<< tauqmatch <<std::endl;

	file->Write();
}


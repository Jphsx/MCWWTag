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
	double pi = 3.1416;

 
	for(int i=0; i<= ncuts; i++){
		 char cuts[100];
         sprintf(cuts, "_%d", i);
         std::string cutnum(cuts);

		/* init histograms */
		WmassMuon[i] = new TH1D(("Wmassmuon"+cutnum).c_str(),"W^{#pm} Mass, with true #mu",100, 50.0, 120.0 );
		WmassTau[i] = new TH1D(("Wmasstau"+cutnum).c_str(),"W^{#pm} Mass, with true #tau",100, 50.0, 120.0 );
        qqmassMuon[i] = new TH1D(("qqmassmuon"+cutnum).c_str(),"qq Mass, with true #mu",100,50.0,120.0);
		qqmassTau[i] = new TH1D(("qqmasstau"+cutnum).c_str(),"qq Mass, with true #tau",100,50.0,120.0);
		WEMuon[i] = new TH1D(("WEmuon"+cutnum).c_str(),"W^{#pm} Energy, with true #mu",100, 25.0, 300.0);
		WETau[i] = new TH1D(("WEtau"+cutnum).c_str(),"W^{#pm} Energy, with true #tau ",100, 25.0, 300.0 );
		EtotalMuon[i] = new TH1D(("EtotalMuon"+cutnum).c_str(),"WW Total Energy, with true #mu",100,100,550); 
		EtotalTau[i] = new TH1D(("EtotalTau"+cutnum).c_str(),"WW Total Energy, with true #tau",100,100,550);
		//TH1D* Wm_cosTheta;
		LjetMassMuon[i]=new TH1D(("Ljetmassmuon"+cutnum).c_str(),"Mass of Lepton Jet with true #mu",100, 0.0 ,20.0 );
		LjetMassTau[i]=new TH1D(("Ljetmasstau"+cutnum).c_str(),"Mass of Lepton Jet with true #tau",100, 0.0, 20.0 );

		costhetawMuon[i] = new TH1D(("costhetawMuon"+cutnum).c_str(), "production angle of W^- in Muon event",100,-1.0,1.0);
		thetaLMuon[i] = new TH1D(("thetaLMuon"+cutnum).c_str(), "polar angle of CM lepton in Muon event",100, 0.0,  pi);
		phiLMuon[i] = new TH1D(("phiLMuon"+cutnum).c_str(), "azimuthal angle of CM Lepton in Muon event", 100,-pi,pi);
		thetaHMuon[i] = new TH1D(("thetaHMuon"+cutnum).c_str(), "polar angle of CM quark in Muon event",100,0.0,pi);
		phiHMuon[i] = new TH1D(("phiHMuon"+cutnum).c_str(),"azimuthal angle of CM quark in Muon event", 100,-pi,pi);

		costhetawTau[i] = new TH1D(("costhetawTau"+cutnum).c_str(), "production angle of W^- in Tau event",100,-1.0,1.0);
		thetaLTau[i] = new TH1D(("thetaLTau"+cutnum).c_str(), "polar angle of CM lepton in Tau event",100, 0.0, pi);
		phiLTau[i] = new TH1D(("phiLTau"+cutnum).c_str(), "azimuthal angle of CM Lepton in Tau event", 100,-pi,pi);
		thetaHTau[i] = new TH1D(("thetaHTau"+cutnum).c_str(), "polar angle of CM quark in Tau event",100,0.0,pi);
		phiHTau[i] = new TH1D(("phiHTau"+cutnum).c_str(),"azimuthal angle of CM quark in Tau event", 100,-pi,pi);
	
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
//recursive function to go through and look at the decay chain of a particle
//look at specifically charged particles
void WWAnalysis::exploreDaughterTracks(MCParticle* p , std::vector<MCParticle*>& FSP){
	if(p->isCreatedInSimulation()) return;

	std::cout<<p->id()<<" ";
	std::cout<<p->getPDG()<<" -> ";
	std::vector<MCParticle*> d = p->getDaughters();
	for(int i=0; i< d.size(); i++){
		if( ( d.at(i)->isCreatedInSimulation() ) ){
		//this is an initial final state particle
			FSP.push_back(d.at(i));
		}
		if( (! d.at(i)->isCreatedInSimulation()) ){//&& (d.at(i)->getCharge() != 0) ){
			std::cout<< "( "<< d.at(i)->id()<<" "<<d.at(i)->getPDG() <<" "<< d.at(i)->isDecayedInTracker()<<" "<< d.at(i)->isDecayedInCalorimeter()<<" ) ";
		}
	}
	std::cout<<std::endl;
	for(int i=0; i<d.size(); i++){
		exploreDaughterTracks(d.at(i));
	}
		
}
//look at all particles
/*void WWAnalysis::exploreDaughterParticles(){

}*/
//looks at the number of charged particles/ total particles in the lepton identified jet
//looks at the number of charged particles/ total particles produced directly from the true lepton
void WWAnalysis::getJetMultiplicities(){



  //get the number of particles/tracks for the jet identified as a lepton
  lnparts = _jets.at(ljet_index)->getParticles().size();

  std::vector<ReconstructedParticle*> lparts = _jets.at(ljet_index)->getParticles();
  lntracks = 0;

  for(int i=0; i< lparts.size(); i++){
	if( lparts.at(i)->getCharge() != 0 ){
		lntracks++;
	}
  }
	

  //use the globally stored parent particle, our true lepton is a daughter of the mcparent
  std::vector<MCParticle*> daughters = parent->getDaughters();
  //find the lepton and look at what it directly produces

    std::vector<MCParticle*> lmcFSP{};

  
  for(int i=0; i<daughters.size(); i++){
	if(daughters.at(i)->getPDG() == lpdg){
		//found the lepton
		exploreDaughterTracks(daughters.at(i), lmcFSP );
		std::vector<MCParticle*> ldaughters = daughters.at(i)->getDaughters();
		for(int j=0; j<ldaughters.size(); j++){
			std::cout<< " l daughters "<< ldaughters.at(j)->getPDG();
		}
	}
	
  }
	std::cout<<std::endl;
	std::cout<<"lfsp ";
  for(int i=0; i<lmcFSP.size(); i++){
	std::cout<<lmcFSP.at(i)->getPDG()<<" ";
  }

  std::cout<<std::endl;
  


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
		std::vector<int> quarks{ 1, 2, 3, 4, 5, -1, -2, -3, -4, -5};
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
					lpdg = 13;
				}
				else{
					trueq = 1;
					lpdg = -13;
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
					lpdg = 15;
				}
				else{
					trueq = 1;
					lpdg = -15;
				}	
			}
			//if we have found the true decay set break out of the mcpart vec loop
			return _mcpartvec.at(i);
			//break;
		}//end if 2qrk & 2 lep

	}//end mcpartvec loop

	//if nothing is found return nulll
	return NULL;

}
/* populate the tlvs based on the identified lepton jet */
void WWAnalysis::populateTLVs(int lindex){

	std::vector<TLorentzVector*> tempjets(_jets.size());
	for(int i=0; i<_jets.size(); i++){
	
		TLorentzVector* j = new TLorentzVector();
		j->SetXYZM(_jets.at(i)->getMomentum()[0], _jets.at(i)->getMomentum()[1], _jets.at(i)->getMomentum()[2], _jets.at(i)->getMass() );
		tempjets.at(i) = j;

		std::cout<<_jets.at(i)->getMomentum()[0]<<" "<< _jets.at(i)->getMomentum()[1]<<" "<<_jets.at(i)->getMomentum()[2]<< " "<< _jets.at(i)->getMass()<<std::endl;
	}
	
	//save the tlv vector globally
	jets = tempjets;
	
	//Wl = new TLorentzVector();
	Wqq = new TLorentzVector();
	TLorentzVector temp1;

	//loop over the new tlv jets and make wl and wqq
	for(int i=0; i<tempjets.size(); i++){
		if( i == lindex ){
			//right now Wl will be missing its neutrino
			Wl = new TLorentzVector();
			Wl->SetXYZM( tempjets.at(i)->Px(), tempjets.at(i)->Py(), tempjets.at(i)->Pz(), tempjets.at(i)->M());
		}
		else{
			temp1 += *tempjets.at(i);
		}
	}
	Wqq->SetXYZM(temp1.Px(), temp1.Py(), temp1.Pz(), temp1.M() );
	
	

	//figure out the muon 
	double missingPx= -(Wl->Px() + Wqq->Px());
	double missingPy= -(Wl->Py() + Wqq->Py());
	double missingPz= -(Wl->Pz() + Wqq->Pz());

	std::cout<<"missing P "<< missingPx<<" "<<missingPy<<" "<<missingPz<<std::endl;
	//create the tlv neutrino
	nu = new TLorentzVector();
	nu->SetXYZM(missingPx, missingPy, missingPz, 0.0);

	//add the neutrino to complete the leptonic W
	TLorentzVector temp2;
	temp2 = *Wl + *nu;
	Wl->SetXYZM(temp2.Px(), temp2.Py(), temp2.Pz(), temp2.M());


}
//populate W rest fram versions of the jets to access TGC observables
void WWAnalysis::populateCMTLVs(){
	//TVector3 Wqqboost(Wqq->Px(),Wqq->Py(),Wqq->Pz());
	//TVector3 Wlboost(Wl->Px(),Wl->Py(),Wl->Pz());

	TVector3 Wqqboost = Wqq->BoostVector();
	TVector3 Wlboost = Wl->BoostVector();

	Wqqboost = -Wqqboost;
	Wlboost = -Wlboost;

	Wqqboost.Print();
	Wlboost.Print();

	std::vector<TLorentzVector*> cmtempvec(_jets.size());
	TLorentzVector cmtemp;
	for(int i=0; i<cmtempvec.size(); i++){
		cmtemp = *(jets.at(i));
		std::cout<<cmtemp.Px()<<" "<<cmtemp.Py()<<" "<<cmtemp.Pz()<<" "<<cmtemp.M()<<std::endl;
		if(i == ljet_index){
			cmtemp.Boost(Wlboost.X(),Wlboost.Y(),Wlboost.Z());
		}
		else{
			cmtemp.Boost(Wqqboost.X(),Wqqboost.Y(),Wqqboost.Z());
		}
		cmtempvec.at(i) = new TLorentzVector();
		cmtempvec.at(i)->SetXYZM(cmtemp.Px(),cmtemp.Py(),cmtemp.Pz(),cmtemp.M());
	}
	std::cout<<"cm"<<std::endl;
	for(int i=0; i<cmtempvec.size(); i++){
		std::cout<<cmtempvec.at(i)->Px()<<" "<<cmtempvec.at(i)->Py()<<" "<<cmtempvec.at(i)->Pz()<<" "<<cmtempvec.at(i)->M()<<std::endl;
	}
	CMJets = cmtempvec;
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
	WmassMuon[histNumber]->Fill(Wl->M() );
	qqmassMuon[histNumber]->Fill(Wqq->M() );
	WEMuon[histNumber]->Fill(Wqq->E() );
	WEMuon[histNumber]->Fill(Wl->E() );
	EtotalMuon[histNumber]->Fill(Wqq->E() + Wl->E() );

	LjetMassMuon[histNumber]->Fill( _jets.at( ljet_index )->getMass() );

	//TGC stuff
	costhetawMuon[histNumber]->Fill(getCosThetaW());
	thetaLMuon[histNumber]->Fill( CMJets.at(ljet_index)->Theta());
	phiLMuon[histNumber]->Fill( CMJets.at(ljet_index)->Phi());
	for(int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHMuon[histNumber]->Fill( CMJets.at(i)->Theta()); 
			phiHMuon[histNumber]->Fill( CMJets.at(i)->Phi());
		}
	} 
		
}
void WWAnalysis::FillTauHistos(int histNumber){

	WmassTau[histNumber]->Fill( Wqq->M() );
	WmassTau[histNumber]->Fill( Wl->M() );
	qqmassTau[histNumber]->Fill( Wqq->M() );
	WETau[histNumber]->Fill( Wqq->E() );
	WETau[histNumber]->Fill(Wl->E() );
	EtotalTau[histNumber]->Fill(Wqq->E() + Wl->E());

	LjetMassTau[histNumber]->Fill( _jets.at(ljet_index)->getMass() );

	//TGC stuff
	costhetawTau[histNumber]->Fill(getCosThetaW());
	thetaLTau[histNumber]->Fill( CMJets.at(ljet_index)->Theta());
	phiLTau[histNumber]->Fill( CMJets.at(ljet_index)->Phi());
	for(int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHTau[histNumber]->Fill( CMJets.at(i)->Theta()); 
			phiHTau[histNumber]->Fill( CMJets.at(i)->Phi());
		}
	} 
}

void WWAnalysis::processEvent( LCEvent * evt ) {

 FindMCParticles(evt);
 FindJets(evt);
 std::cout << "======================================== event " << nEvt << std::endl ;


	//bools to characterize the true lepton decay for this event
	isTau = false;
	isMuon = false;

	//from the MCParticles find what type of semileptonic decay is present
    //return the parent mcparticle that has the qqlnu decay
	parent = classifyEvent(isTau, isMuon, trueq);


	//now assess jets
	//keep the index on _jets of the jet we consider to be the lepton
	ljet_index = identifyLeptonJet( _jets );

	//get the charge of the lepton jet
	lq = getLeptonJetCharge( _jets.at(ljet_index) );

	//assess jet multiplicity
	//fill variables pertaining to leptonic jet numbers of particles
	getJetMultiplicities(); 


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


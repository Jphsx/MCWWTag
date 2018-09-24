#include "MCWWTag.h"

MCWWTag aMCWWTag;


MCWWTag::MCWWTag() : Processor("MCWWTag") {


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

void MCWWTag::init() {

  streamlog_out(DEBUG) << "   init called  " << std::endl;
  // usually a good idea to
  printParameters() ;
  nEvt = 0;

}

void MCWWTag::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ;
}

bool MCWWTag::FindMCParticles( LCEvent* evt ){
   
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

bool MCWWTag::FindJets( LCEvent* evt ) {

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
int MCWWTag::identifyLeptonJet( std::vector<ReconstructedParticle*> jets){

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
int MCWWTag::getLeptonJetCharge( ReconstructedParticle* ljet ){
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

void MCWWTag::processEvent( LCEvent * evt ) {
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
	

 nEvt++;
}
void MCWWTag::end(){

	std::cout<<" nelec "<<nelec<<" nmuon "<< nmuon <<" ntau "<< ntau << std::endl;
	std::cout<<" ndwn "<<ndwn<<" nup "<<nup<<" nstr "<<nstr<<" nchm "<<nchm<<std::endl;

	std::cout<<" nevents "<< nEvt << " mu q match "<< muonqmatch <<  " tau q match "<< tauqmatch <<std::endl;
}


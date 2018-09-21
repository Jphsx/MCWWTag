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

void MCWWTag::processEvent( LCEvent * evt ) {
 FindMCParticles(evt);
 std::cout << "======================================== event " << nEvt << std::endl ;


	for(int i=0; i<_mcpartvec.size(); i++){
		std::vector<int> parentpdgs{};
		std::vector<int> daughterpdgs{};
		std::vector<MCParticle*> mcparents{};
		std::vector<MCParticle*> daughters{};
		daughters = _mcpartvec.at(i)->getDaughters();
		for(int j = 0; j<daughters.size(); j++){
			daughterpdgs.push_back(daughters.at(j)->getPDG());
			
		}
		
		/*//if( _mcpartvec.at(i)->getParents().size() == 0 ){
			std::vector<MCParticle*> daughters{};
			daughters = _mcpartvec.at(i)->getDaughters();
				std::cout<<_mcpartvec.at(i)->isCreatedInSimulation()<<" "<<_mcpartvec.at(i)->getPDG()<<" -> ";
			for(int j=0; j<daughters.size(); j++){
				 std::cout<<daughters.at(j)->getPDG() << " ";
			}
			std::cout<<std::endl;
		//}*/
	/*	if(_mcpartvec.at(i)->getParents().size() == 0 ){
			parentpdgs.push_back(-0);
		}
		else{
			mcparents = _mcpartvec.at(i)->getParents();
			for(int j=0; j< mcparents.size(); j++){
				parentpdgs.push_back( mcparents.at(j)->getPDG() );
			}
		}
		
		//print junk
		for(int j=0; j<parentpdgs.size(); j++){
			std::cout<<parentpdgs.at(j)<<" ";
		}
		std::cout<< " -> "<<_mcpartvec.at(i)->getPDG()<<" -> ";
		for(int j=0; j<daughters.size(); j++){
			std::cout<<daughters.at(j)->getPDG()<<" ";
		}
		std::cout<<std::endl;
	*/

		//allowed quarks
		std::vector<int> quarks{ 1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6};
,		std::vector<int> leptons{11, 12, 13, 14, 15, 16, 17, 18, -11, -12, -13, -14, -15, -16, -17, -18};
		//we require exactly 2 elements from leptons and 2 from quarks
		int lep=0;
		int qrk=0;

		
		for(int j=0; j< daughters.size(); j++){
			for(int k=0; k<quarks.size(); k++){
				qrk += std::count(daughterpdgs.begin(),daughterpdgs.end(),quarks.at(k));
			}
			for(int k=0; k<quarks.size(); k++){
				lep += std::count(daughterpdgs.begin(),daughterpdgs.end(),leptons.at(k));
			}
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
		 break;
		}

	}

 nEvt++;
}
void MCWWTag::end(){
	
}


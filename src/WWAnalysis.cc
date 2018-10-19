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

  _nRun = 0;
  _nEvt = 0;

  nEvt = 0;

	
	file = new TFile("file.root","RECREATE");
	double pi = 3.1416;

 
	for(int i=0; i<= ncuts; i++){
		 char cuts[100];
         sprintf(cuts, "_%d", i);
         std::string cutnum(cuts);

		/* init histograms */
		//some general plots
		WmassMuon[i] = new TH1D(("Wmassmuon"+cutnum).c_str(),"W^{#pm} Mass from l #nu, with true #mu",100, 50.0, 120.0 );
		WmassTau[i] = new TH1D(("Wmasstau"+cutnum).c_str(),"W^{#pm} Mass from l #nu, with true #tau",100, 50.0, 120.0 );
                qqmassMuon[i] = new TH1D(("qqmassmuon"+cutnum).c_str(),"W #rightarrow qq Mass, #mu;Mass [GeV];Events Per Bin",100,50.0,120.0);
		qqmassTau[i] = new TH1D(("qqmasstau"+cutnum).c_str(),"W #rightarrow qq Mass;Mass [GeV];Events Per Bin",100,50.0,120.0);
		WEMuon[i] = new TH1D(("WEmuon"+cutnum).c_str(),"W^{#pm} Energy, with true #mu",100, 25.0, 300.0);
		WETau[i] = new TH1D(("WEtau"+cutnum).c_str(),"W^{#pm} Energy, with true #tau ",100, 25.0, 300.0 );
		EtotalMuon[i] = new TH1D(("EtotalMuon"+cutnum).c_str(),"WW Total Energy, with true #mu",100,100,550); 
		EtotalTau[i] = new TH1D(("EtotalTau"+cutnum).c_str(),"WW Total Energy, with true #tau",100,100,550);
		
		//TGC stuff
		LjetMassMuon[i]=new TH1D(("Ljetmassmuon"+cutnum).c_str(),"Mass of Lepton Jet with true #mu",100, 0.0 ,20.0 );
		LjetMassTau[i]=new TH1D(("Ljetmasstau"+cutnum).c_str(),"Mass of Lepton Jet with true #tau",100, 0.0, 20.0 );

		costhetawMuon[i] = new TH1D(("costhetawMuon"+cutnum).c_str(), "Production Angle of W^- #mu;cos #theta_W ; Events Per Bin",100,-1.0,1.0);
		thetaLMuon[i] = new TH1D(("thetaLMuon"+cutnum).c_str(), "polar angle of CM lepton in Muon event",100, 0.0,  pi);
		phiLMuon[i] = new TH1D(("phiLMuon"+cutnum).c_str(), "azimuthal angle of CM Lepton in Muon event", 100,-pi,pi);
		thetaHMuon[i] = new TH1D(("thetaHMuon"+cutnum).c_str(), "polar angle of CM quark in Muon event",100,0.0,pi);
		phiHMuon[i] = new TH1D(("phiHMuon"+cutnum).c_str(),"azimuthal angle of CM quark in Muon event", 100,-pi,pi);

		costhetawTau[i] = new TH1D(("costhetawTau"+cutnum).c_str(), "Production Angle of W^- #tau; cos #theta_W ; Events Per Bin",100,-1.0,1.0);
		thetaLTau[i] = new TH1D(("thetaLTau"+cutnum).c_str(), "polar angle of CM lepton in Tau event",100, 0.0, pi);
		phiLTau[i] = new TH1D(("phiLTau"+cutnum).c_str(), "azimuthal angle of CM Lepton in Tau event", 100,-pi,pi);
		thetaHTau[i] = new TH1D(("thetaHTau"+cutnum).c_str(), "polar angle of CM quark in Tau event",100,0.0,pi);
		phiHTau[i] = new TH1D(("phiHTau"+cutnum).c_str(),"azimuthal angle of CM quark in Tau event", 100,-pi,pi);


        
	 //jet information histograms
   		leptonMCNPartsMuon[i] = new TH1D(("leptonMCNPartsMuon"+cutnum).c_str(), "True Visible Particles in #mu Event;Number of Particles;Entries ", 11,-0.5,10.5);
		leptonMCNTracksMuon[i] = new TH1D(("leptonMCNTracksMuon"+cutnum).c_str(), "True Visible Tracks in #mu  Event;Number of Tracks;Entries ",11,-0.5,10.5);
		jetNpartsMuon[i] = new TH1D(("jetNpartsMuon"+cutnum).c_str(), "Visible Particles per Jet in Muon Event",51,-0.5,50.5);
        jetNtracksMuon[i] = new TH1D(("jetNtracksMuon"+cutnum).c_str(), "Visible Tracks per Jet in Muon Event", 51,-0.5,50.5);		
minjetNpartsMuon[i] = new TH1D(("minjetNpartsMuon"+cutnum).c_str(), "Visible Particle of Jet with least Particles in Muon Event",26,-0.5,25.5);
		minjetNtracksMuon[i] = new TH1D(("minjetNtracksMuon"+cutnum).c_str(), "Visible Tracks of Jet with least Particlesin Muon Event",26,-0.5,25.5);

		leptonMCNPartsTau[i] = new TH1D(("leptonMCNPartsTau"+cutnum).c_str(), "True Visible Particles in #tau Event;Number of Particles;Entries", 11,-0.5,10.5);
		leptonMCNTracksTau[i] = new TH1D(("leptonMCNTracksTau"+cutnum).c_str(), "True Visible Tracks in #tau Event;Number of Particles;Entries",11,-0.5,10.5);
		jetNpartsTau[i] = new TH1D(("jetNpartsTau"+cutnum).c_str(), "Visible Particles per Jet in Tau Event",21,-0.5,20.5);
        jetNtracksTau[i] = new TH1D(("jetNtracksTau"+cutnum).c_str(), "Visible Tracks per Jet in Tau Event", 51,-0.5,50.5);
		minjetNpartsTau[i] = new TH1D(("minjetNpartsTau"+cutnum).c_str(), "Visible Particle of Jet with least Particles in Tau Event",26,-0.5,25.5);
		minjetNtracksTau[i] = new TH1D(("minjetNtracksTau"+cutnum).c_str(), "Visible Tracks of Jet with least Particlesin Tau Event",26,-0.5,25.5);

		//some jet leading track parameters
		ljetleadingd0Muon[i] = new TH1D(("ljetleadingd0Muon"+cutnum).c_str(), "d0 of the leading track of lepton jet in muon event",50,-0.1,0.1);
		ljetleadingptMuon[i] = new TH1D(("ljetleadingptMuon"+cutnum).c_str(), "p of the leading track of lepton jet in muon event", 100,0.0, 200.0);
		ljetd0relerrMuon[i] = new TH1D(("ljetd0relerrMuon"+cutnum).c_str(),"d0/ #delta d0  of leading track of lepton jet in muon event",20,0.0,100);
		qjetleadingd0Muon[i] = new TH1D(("qjetleadingd0Muon"+cutnum).c_str(), "d0 of the leading track of q jet in muon event",50,-0.1,0.1);
		qjetleadingptMuon[i] = new TH1D(("qjetleadingptMuon"+cutnum).c_str(), "p of the leading track of q jet in muon event", 100,0.0, 200.0);
		qjetd0relerrMuon[i] = new TH1D(("qjetd0relerrMuon"+cutnum).c_str(),"d0/ #delta d0  of leading track of q jet in muon event",100,0.0,100);
		ljetleadingd0Tau[i] = new TH1D(("ljetleadingd0Tau"+cutnum).c_str(), "d0 of the leading track of lepton jet in Tau event",50,-.1,.1);
		ljetleadingptTau[i] = new TH1D(("ljetleadingptTau"+cutnum).c_str(), "p of the leading track of lepton jet in Tau event", 100,0.0, 200.0);
		ljetd0relerrTau[i] = new TH1D(("ljetd0relerrTau"+cutnum).c_str(),"d0/ #delta d0  of leading track of lepton jet in Tau event",100,0.0,100);
		qjetleadingd0Tau[i] = new TH1D(("qjetleadingd0Tau"+cutnum).c_str(), "d0 of the leading track of q jet in Tau event",50,-.1,.1);
		qjetleadingptTau[i] = new TH1D(("qjetleadingptTau"+cutnum).c_str(), "p of the leading track of q jet in Tau event", 100,0.0, 200.0);
		qjetd0relerrTau[i] = new TH1D(("qjetd0relerrTau"+cutnum).c_str(),"d0/ #delta d0  of leading track of q jet in Tau event",100,0.0,100);
	
		psiljetmclMuon[i] = new TH1D(("psiljetmclMuon"+cutnum).c_str(),"angle between the lepton jet and the true lepton, muon;cos #theta; Entries per 0.02",100,-1.0,1.0);	
		psiljetmclTau[i]  = new TH1D( ("psiljetmclTau"+cutnum).c_str(),"angle between the lepton jet and the true lepton, tau;cos #theta; Entries per 0.02 bin",100,-1.0,1.0);
		/* end init histograms */
		}
// TTree similar to ttbar.cc to start with

     _tree = new TTree("tree", "tree");
     _tree->Branch("runNumber", &_nRun, "runNumber/I");
     _tree->Branch("eventNumber", &_nEvt, "eventNumber/I");
//     _tree->Branch("crossSection", &_xsec, "crossSection/F");
//     _Process = new TString();
//     _tree->Branch("Process","TString",&_Process,16000,0);
     _tree->Branch("isMuon", &isMuon, "isMuon/O");
     _tree->Branch("isTau", &isTau, "isTau/O");
     _tree->Branch("leptonCharge", &trueq,"leptonCharge/I");

     for(int i = 0 ; i < 4 ; i++)
       {
	 _MCf[i] = new TLorentzVector();
	 std::stringstream name;
	 name << "MCf" << i;
	 _tree->Branch(name.str().c_str(),"TLorentzVector",&_MCf[i],16000,0);
	 name << "_PDG";
	 _tree->Branch(name.str().c_str(), &_MCfpdg[i], (name.str()+"/I").c_str());
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
			for(unsigned int i=0;i<collection->getNumberOfElements();i++){
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
      			for(unsigned int i=0; i<collection->getNumberOfElements(); i++){
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
	for(unsigned int i=0; i<_jets.size(); i++){
		std::cout<<"jet "<<i<<" has "<< _jets.at(i)->getParticles().size() << " particles "<<std::endl;
		if( _jets.at(i)->getParticles().size() < nparticles ){
			indexofminjet = i;
			nparticles = _jets.at(i)->getParticles().size();
		}
	}
	
	return indexofminjet;

}
int WWAnalysis::identifyLeptonJet_bySeparation(std::vector<ReconstructedParticle*> jets){

	double pi = 3.14159;
	TVector3 j1(jets.at(0)->getMomentum()[0],jets.at(0)->getMomentum()[1], jets.at(0)->getMomentum()[2]);
	TVector3 j2(jets.at(1)->getMomentum()[0],jets.at(1)->getMomentum()[1], jets.at(1)->getMomentum()[2]);
	TVector3 j3(jets.at(2)->getMomentum()[0],jets.at(2)->getMomentum()[1], jets.at(2)->getMomentum()[2]);

        double psi12 = acos(j1.Dot(j2)/ ( j1.Mag() * j2.Mag() ));
	double psi13 = acos(j1.Dot(j3)/ ( j1.Mag() * j3.Mag() ));
	double psi23 = acos(j2.Dot(j3)/ ( j2.Mag() * j3.Mag() ));

	if(psi12<0) psi12=psi12+2*pi;
	if(psi13<0) psi13=psi13+2*pi;
	if(psi23<0) psi23=psi23+2*pi;
	
	std::vector<double> angles{ psi12, psi13, psi23};
	std::vector<double> excludedjet{ 3, 2, 1 };
	//what angle is the smallest?
	double min = 999;
	int minindex= -1;
	for(unsigned int i=0; i<angles.size(); i++){
		if(angles.at(i) < min){
			min = angles.at(i);
			minindex = i;
		}
	}
	
	return (excludedjet.at(minindex) -1);
}
//TODO vertex fitting for lepton id with tau
//int WWAnalysis::idljet_vertexfit( std::vec
/* identifies Lepton jet charge by taking the charge of the leading track from the jet */
int WWAnalysis::getLeptonJetCharge( ReconstructedParticle* ljet ){
	//assign by leading track charge? or charge sum of reco parts?
	std::vector<ReconstructedParticle*> jetparts = ljet->getParticles();
	int totalcharge = 0;
	int leadingcharge = 0;
	double maxP = 0;
	for(unsigned int i=0; i<jetparts.size(); i++){
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
void WWAnalysis::getAngleOfljetandMCLepton(){
	TVector3 ljet( _jets.at(ljet_index)->getMomentum()[0], _jets.at(ljet_index)->getMomentum()[1], _jets.at(ljet_index)->getMomentum()[2] ); 
	
	int mclindex;
	for(int i=0; i<4; i++){
		if( abs(_MCfpdg[i]) == 13 || abs(_MCfpdg[i])== 15){
			mclindex = i;
		}
	}
	TVector3 mcl( _MCf[mclindex]->Px(), _MCf[mclindex]->Py(), _MCf[mclindex]->Pz() );

	psi_mcl_ljet = ljet.Dot(mcl)/( ljet.Mag() * mcl.Mag());
	
}
bool WWAnalysis::allChildrenAreSimulation(MCParticle* p){
	std::vector<MCParticle*> d = p->getDaughters();
	bool flag = true;
	for(unsigned int i=0; i<d.size(); i++){
		if( ! d.at(i)->isCreatedInSimulation() ){
			flag= false;
		}
	}
	return flag;
}
//recursive function to go through and look at the decay chain of a particle
//look at specifically charged particles
void WWAnalysis::exploreDaughterParticles(MCParticle* p , std::vector<MCParticle*>& FSP){
	if(p->isCreatedInSimulation()) return;

	//std::cout<<p->id()<<" ";
	//std::cout<<p->getPDG()<<" -> ";
	std::vector<MCParticle*> d = p->getDaughters();
	
	for(unsigned int i=0; i< d.size(); i++){
		if( (! d.at(i)->isCreatedInSimulation() ) && ( allChildrenAreSimulation(d.at(i)) || (d.at(i)->getDaughters().size()==0)  ) ){
		//this is an initial final state particle
			FSP.push_back(d.at(i));
		}
		//if( (! d.at(i)->isCreatedInSimulation()) ){//&& (d.at(i)->getCharge() != 0) ){
		//	std::cout<< "( "<< d.at(i)->id()<<" "<<d.at(i)->getPDG() <<" "<< d.at(i)->isDecayedInTracker()<<" "<< d.at(i)->isDecayedInCalorimeter()<<" ) ";
		//}
	}
	//std::cout<<std::endl;
	for(unsigned int i=0; i<d.size(); i++){
		exploreDaughterParticles(d.at(i), FSP);
	}
		
}
//looks at the number of charged particles/ total particles in the lepton identified jet
//looks at the number of charged particles/ total particles produced directly from the true lepton
void WWAnalysis::getJetMultiplicities(){



  //get the number of particles/tracks for the jet identified as a lepton
  lnparts = _jets.at(ljet_index)->getParticles().size();

  std::vector<ReconstructedParticle*> lparts = _jets.at(ljet_index)->getParticles();
  lntracks = 0;

  for(unsigned int i=0; i< lparts.size(); i++){
	if( lparts.at(i)->getCharge() != 0 ){
		lntracks++;
	}
  }
	

  //use the globally stored parent particle, our true lepton is a daughter of the mcparent
  std::vector<MCParticle*> daughters = parent->getDaughters();
  //find the lepton and look at what it directly produces

    std::vector<MCParticle*> lmcFSP{};
	//quarks have the same children???
 //   std::vector<MCParticle*> q1FSP{};
 //   std::vector<MCParticle*> q2FSP{};

  
  for(unsigned int i=0; i<daughters.size(); i++){
	int pdg = daughters.at(i)->getPDG();
	if(pdg == lpdg){
		//found the lepton
		exploreDaughterParticles(daughters.at(i), lmcFSP );
	}
/*	else if( abs(pdg) != 12 || abs(pdg) != 14 || abs(pdg) != 16 ){
		//found a quark
		if(q1FSP.size()==0){
			
			exploreDaughterParticles( daughters.at(i), q1FSP );
		}
		else{
			//exploreDaughterParticles( daughters.at(i), q2FSP );
		}
	}
*/
	
  }
	std::cout<<std::endl;
	std::cout<<"lfsp ";
  for(unsigned int i=0; i<lmcFSP.size(); i++){
	std::cout<<lmcFSP.at(i)->getPDG()<<" ";
  }

  std::cout<<std::endl;

  //loop over visible particles  in lfsp and save jet details
  int countparts=0;
  int counttracks=0;
  for(unsigned int i=0; i<lmcFSP.size(); i++){
		if( (abs(lmcFSP.at(i)->getPDG()) == 12) || (abs(lmcFSP.at(i)->getPDG()) == 14) || (abs(lmcFSP.at(i)->getPDG()) == 16)){
		//	std::cout<<"neut skipped"<<std::endl;
			continue;
		}
		else{
			
			//not a neutrino
			if(lmcFSP.at(i)->getCharge() != 0){
				counttracks++;
			}
			//count tracks and neutrals
			countparts++;
		//f	std::cout<<"part counted"<<std::endl;
		}
  }

  lnmcparts = countparts;
  lnmctracks = counttracks;

  //count all the tracks and particles from fast jets
  int countjetparts=0;
  int countjettracks=0;
  //temp vectors
   std::vector<int> nparts(_jets.size());
   std::vector<int> ntrks(_jets.size());
  //loop over all jets
  for(unsigned int i=0; i<_jets.size(); i++){
	std::vector<ReconstructedParticle*> p = _jets.at(i)->getParticles();
	countjetparts = p.size();
	for(unsigned int j=0; j<p.size(); j++){
		if(p.at(j)->getCharge() != 0){
			countjettracks++;
		}
	}
	nparts.at(i) = countjetparts;
    ntrks.at(i) = countjettracks;
	countjetparts = 0;
	countjettracks = 0;
  }

  jetNparts = nparts;
  jetNtracks = ntrks;


}
void WWAnalysis::analyzeLeadingTracks(){
	ReconstructedParticle* leader;
	std::vector<ReconstructedParticle*> d;
	std::vector<Track*> dt;
	double maxP = -9999;
    double maxtP = -9999;
	int maxindex= -1;
	int maxtindex = -1;
	double tp;
	const double* mom;
	double p;

	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
    const double BField = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
  	const double eB = BField*c*mm2m*eV2GeV; 
 double p_t=0;
                        double p_z=0 ;
	for(unsigned int i=0; i< _jets.size(); i++){
		d = _jets.at(i)->getParticles();
		for(unsigned int j=0; j< d.size(); j++){
			if(d.at(j)->getCharge() !=0){
				mom = d.at(j)->getMomentum();
				p = std::sqrt( mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
			//	std::cout<<"om "<<om<<" "<<fabs(om)<<" "<<minOm<<std::endl;
			//	std::cout<<" p "<<p<<std::endl;
				if( p > maxP){
					maxP = p;
					maxindex=j;
				}//end max reset
			}//end charge condition
		}//end jet particles

		if(maxindex == -1) { std::cout<< "continuing" <<std::endl; continue;} //no tracks in this jet

		//the maxpt track may have more than 1 associated tracks so pick out the one with highest pt
		dt = d.at(maxindex)->getTracks();
		for(unsigned int j=0; j< dt.size();j++){
			p_t = abs(dt.at(j)->getOmega())*eB;
			p_z = p_t*dt.at(j)->getTanLambda();
			tp = std::sqrt(p_t*p_t + p_z*p_z);
			if(tp > maxtP){
				maxtP = tp;
				maxtindex=j;
			}

		}		
		//look at track of this particle
	    Track* t = dt.at(maxtindex);
		std::cout<<"jet "<<i<< std::endl;
		
			std::cout<<t->getD0()<<" ";
			std::cout<<t->getPhi()<<" ";
			std::cout<<t->getOmega()<<" ";
			std::cout<<t->getZ0()<<" ";
			std::cout<<t->getTanLambda()<<" ";
		 
		std::cout<<std::endl;
    
			 p_t = abs(t->getOmega())*eB;
                         p_z = p_t*t->getTanLambda();
                        

	
		if(i == ljet_index){
			leadingptljet = std::sqrt(p_t*p_t + p_z*p_z);
  
			leadingd0ljet = t->getD0();
		//	leadingd0relerrljet = sqrt(t->getCovMatrix()[0])/fabs(t->getD0());
			leadingd0relerrljet = fabs(t->getD0())/sqrt(t->getCovMatrix()[0]);
		}
		else{
			//we lazily just examing 1 q jet for now
			leadingptqjet = std::sqrt(p_t*p_t + p_z*p_z);

			leadingd0qjet = t->getD0();
		//	leadingd0relerrqjet = sqrt(t->getCovMatrix()[0])/fabs(t->getD0());
			leadingd0relerrqjet = fabs(t->getD0())/sqrt(t->getCovMatrix()[0]);
		}
		//reset maxindex and max p new jet
		//
		p_t = 0; 
		p_z = 0;
		maxtP = -9999;
		maxindex= -1;
		maxtindex = -1;
		maxP = -9999;
	}//end jet loop
 
}
/* classify the the event based on the type of lepton in MCParticle info, also set the true charge for that lepton */
/* also tallies the number of muon/electron/tau events */
/* also retrieves the mcparticle which has daughters qqlnu */
MCParticle* WWAnalysis::classifyEvent(bool& isTau, bool& isMuon, int& trueq, TLorentzVector* (&_MCf)[4], int (&_MCfpdg)[4]){

//MCParticle* WWAnalysis::classifyEvent(bool& isTau, bool& isMuon, int& trueq, int (&_MCfpdg)[4]){
	
	for(unsigned int i=0; i<_mcpartvec.size(); i++){
		std::vector<int> parentpdgs{};
		std::vector<int> daughterpdgs{};
		std::vector<MCParticle*> mcparents{};
		std::vector<MCParticle*> daughters{};
		daughters = _mcpartvec.at(i)->getDaughters();
		for(unsigned int j = 0; j<daughters.size(); j++){
			daughterpdgs.push_back(daughters.at(j)->getPDG());
			
		}
		//allowed quarks
		std::vector<int> quarks{ 1, 2, 3, 4, 5, -1, -2, -3, -4, -5};
		std::vector<int> leptons{11, 12, 13, 14, 15, 16, -11, -12, -13, -14, -15, -16};
		//we require exactly 2 elements from leptons and 2 from quarks
		int lep=0;
		int qrk=0;

		//categorize the event for plotting
		for(unsigned int k=0; k<quarks.size(); k++){
			qrk += std::count(daughterpdgs.begin(),daughterpdgs.end(),quarks.at(k));
		}
		for(unsigned int k=0; k<leptons.size(); k++){
			lep += std::count(daughterpdgs.begin(),daughterpdgs.end(),leptons.at(k));
		}

		if( qrk == 2 && lep == 2){
			//found the proper set 
			for(unsigned int j=0; j<parentpdgs.size(); j++){
				std::cout<<parentpdgs.at(j)<<" ";
			}
			std::cout<< " -> "<<_mcpartvec.at(i)->getPDG()<<" -> ";
			for(unsigned int j=0; j<daughters.size(); j++){
				std::cout<<daughters.at(j)->getPDG()<<" ";
			}
			std::cout<<std::endl;

			for(unsigned int j=0; j<daughters.size(); j++){
				std::cout<<daughters.at(j)->getPDG()<<" " 
                                                                    << daughters.at(j)->getMomentum()[0] << " "
                                                                    << daughters.at(j)->getMomentum()[1] << " "
                                                                    << daughters.at(j)->getMomentum()[2] << " "
                                                                    << daughters.at(j)->getEnergy() << " " << std::endl;
                            TLorentzVector mcVec(TVector3(daughters.at(j)->getMomentum()),daughters.at(j)->getEnergy());
                            *_MCf[j] = mcVec;
                            _MCfpdg[j] = daughters.at(j)->getPDG();
			}



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
	for(unsigned int i=0; i<_jets.size(); i++){
	
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
	for(unsigned int i=0; i<tempjets.size(); i++){
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
	for(unsigned int i=0; i<cmtempvec.size(); i++){
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
	for(unsigned int i=0; i<cmtempvec.size(); i++){
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

//	WmassMuon[histNumber]->Fill( Wqq->M() );
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
	for(unsigned int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHMuon[histNumber]->Fill( CMJets.at(i)->Theta()); 
			phiHMuon[histNumber]->Fill( CMJets.at(i)->Phi());
		}
	} 

    //jet details
    leptonMCNPartsMuon[histNumber]->Fill(lnmcparts);
	leptonMCNTracksMuon[histNumber]->Fill(lnmctracks); 
    for(unsigned int i=0; i<jetNparts.size(); i++){
		jetNpartsMuon[histNumber]->Fill(jetNparts.at(i));
		jetNtracksMuon[histNumber]->Fill(jetNtracks.at(i));
	}
	minjetNpartsMuon[histNumber]->Fill(lnparts);
	minjetNtracksMuon[histNumber]->Fill(lntracks);

	//more jet details
	ljetleadingd0Muon[histNumber]->Fill(leadingd0ljet);
	ljetleadingptMuon[histNumber]->Fill(leadingptljet);
	ljetd0relerrMuon[histNumber]->Fill(leadingd0relerrljet);
	qjetleadingd0Muon[histNumber]->Fill(leadingd0qjet);
	qjetleadingptMuon[histNumber]->Fill(leadingptqjet);
	qjetd0relerrMuon[histNumber]->Fill(leadingd0relerrqjet);


	 psiljetmclMuon[histNumber]->Fill(psi_mcl_ljet);
}
void WWAnalysis::FillTauHistos(int histNumber){

//	WmassTau[histNumber]->Fill( Wqq->M() );
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
	for(unsigned int i=0; i<CMJets.size(); i++){
		if( i != ljet_index ){
			thetaHTau[histNumber]->Fill( CMJets.at(i)->Theta()); 
			phiHTau[histNumber]->Fill( CMJets.at(i)->Phi());
		}
	} 
	
	//jet details
    leptonMCNPartsTau[histNumber]->Fill(lnmcparts);
	leptonMCNTracksTau[histNumber]->Fill(lnmctracks);
	for(unsigned int i=0; i<jetNparts.size(); i++){ 
		jetNpartsTau[histNumber]->Fill(jetNparts.at(i));
	    jetNtracksTau[histNumber]->Fill(jetNtracks.at(i));
	}
	minjetNpartsTau[histNumber]->Fill(lnparts);
	minjetNtracksTau[histNumber]->Fill(lntracks);

	ljetleadingd0Tau[histNumber]->Fill(leadingd0ljet);
	ljetleadingptTau[histNumber]->Fill(leadingptljet);
	ljetd0relerrTau[histNumber]->Fill(leadingd0relerrljet);
	qjetleadingd0Tau[histNumber]->Fill(leadingd0qjet);
	qjetleadingptTau[histNumber]->Fill(leadingptqjet);
	qjetd0relerrTau[histNumber]->Fill(leadingd0relerrqjet);

	 psiljetmclTau[histNumber]->Fill(psi_mcl_ljet);

}

void WWAnalysis::processEvent( LCEvent * evt ) {

  _nEvt++;
 // Get Process name and cross section
 //  *_Process = evt->getParameters().getStringVal("Process");
 // _xsec = evt->getParameters().getFloatVal("CrossSection_fb");

 FindMCParticles(evt);
 FindJets(evt);
 std::cout << "======================================== event " << nEvt << std::endl ;


	//bools to characterize the true lepton decay for this event
	isTau = false;
	isMuon = false;

	//from the MCParticles find what type of semileptonic decay is present
        //return the parent mcparticle that has the qqlnu decay
//	parent = classifyEvent(isTau, isMuon, trueq, _MCf[0], _MCfpdg[0]);
	parent = classifyEvent(isTau, isMuon, trueq, _MCf, _MCfpdg);

	//now assess jets
	//keep the index on _jets of the jet we consider to be the lepton
//	ljet_index = identifyLeptonJet( _jets );
	//try using the most separated jet as the lepton jet
	ljet_index = identifyLeptonJet_bySeparation(_jets);

	getAngleOfljetandMCLepton();

	//get the charge of the lepton jet
	lq = getLeptonJetCharge( _jets.at(ljet_index) );

	//assess jet multiplicity
	//fill variables pertaining to leptonic jet numbers of particles
	getJetMultiplicities(); 

	analyzeLeadingTracks();

	//check if the assessed charge matches the true charge of the lepton
	if( trueq == lq){
		std::cout<<" got correct lepton charge "<<std::endl;
		//count the number of times we get it right
//		if(isTau) tauqmatch++;
//		if(isMuon) muonqmatch++;
	}
	else{ 
//		std::cout<<" charge wrong "<<std::endl;
	}

	//build up all the different tlvs for calculation
  	populateTLVs(ljet_index);

    //boost jets to cm for TGC observables
	populateCMTLVs();

	//fill base histograms and produce histos with sequential cuts hist0 is always no cuts
	FillHistos(0);
	//cut #1 trueq == lq, lepton sign is correctly assessed
//	if(trueq == lq){
//		FillHistos(1);
//	}
	//cut #1 require polar angle of q's or lepton to be > cos(theta)=0.995
	bool PolarAngleRequirementMet = true;
	for(unsigned int i=0; i<4; i++){
		if( abs(_MCfpdg[i]) == 12 || abs(_MCfpdg[i]) == 14 || abs(_MCfpdg[i]) == 16 ){
		//dont worry about neutrino angle
			continue;	
		} 
		if( fabs(_MCf[i]->CosTheta()) > 0.995 ){
			//a particle has failed break out
			PolarAngleRequirementMet = false;
			break;
		}
	}
	if(PolarAngleRequirementMet){
		FillHistos(1);
		if(trueq == lq){
			std::cout<<"got correct charge"<<std::endl;
			if(isTau) tauqmatch++;
			if(isMuon) muonqmatch++;
		}
		else{
			std::cout<<"charge wrong "<<std::endl;
		}
		if(psi_mcl_ljet > 0.7){
			std::cout<<"jet loose match"<<std::endl;
			if(isTau) ljetmatchmctau++;
			if(isMuon) ljetmatchmcmuon++;

		}
	}

  _tree->Fill();

 nEvt++;
}
void WWAnalysis::end(){

	/* print out stuff */
	std::cout<<" nelec "<<nelec<<" nmuon "<< nmuon <<" ntau "<< ntau << std::endl;
	std::cout<<" nevents "<< nEvt << " mu q match "<< muonqmatch <<  " tau q match "<< tauqmatch <<std::endl;
	std::cout<<" jet match tau "<<ljetmatchmctau<<" "<<ljetmatchmcmuon<<std::endl;

	file->Write();
}


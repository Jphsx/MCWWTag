#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include <marlin/Global.h>
#include "gear/BField.h"
#include "lcio.h"
#include "TFile.h"
#include "TH1D.h"
#include <vector>
#include <algorithm>
#include "TLorentzVector.h"
#include "TVector3.h"

#include <iostream>
#include <fstream>

#define ncuts 1

using namespace lcio;

	/** WWAnalysis:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class WWAnalysis : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new WWAnalysis ; }

  WWAnalysis(const WWAnalysis&) = delete ;
  WWAnalysis& operator=(const WWAnalysis&) = delete ;

  WWAnalysis() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

 //collection gathering
  bool FindMCParticles( LCEvent* evt );
  bool FindJets( LCEvent* evt ) ;


  // lepton jet functions
  int identifyLeptonJet( std::vector<ReconstructedParticle*> jets);
  int getLeptonJetCharge( ReconstructedParticle* ljet );

  //classify the type of lepton decay and retrieve the
  //mcparticles for qqlnu
  MCParticle* classifyEvent(bool& isTau, bool& isMuon, int& trueq);

  //populate local datastructures (TLVS)
  void populateTLVs(int lindex);
  void populateCMTLVs();

  //helper function to get production angle of W-
  double getCosThetaW();

  //functions to populate histograms
  void FillHistos(int histNumber);
  void FillMuonHistos(int histNumber);
  void FillTauHistos(int histNumber);

  protected:
//event number
  int nEvt{};

//MC information
 //the true parent that contains qqlnu
  MCParticle* parent;
 //bools to characterize the true lepton decay for this event
  bool isTau;
  bool isMuon;
 //the true lepton charge
  int trueq;

//Lepton Jet variables
 //index of the identified lepton on jet vector
  int ljet_index;
 //the assigned charge for identifed lepton jet
  int lq;

//tallies for the number of each type of true lepton per event
  int ntau=0;
  int nmuon=0;
  int nelec=0;

//the total number of unique cuts applied (for histogram indexing)
//  int ncuts = 1;

  //how many times do we get the proper lepton charge?
  //for muons and for leptons separately
  int muonqmatch=0;
  int tauqmatch=0;
  

  //vector to hold the particles for the event
  std::vector<MCParticle*> _mcpartvec{};
  std::vector<ReconstructedParticle*> _jets{};
  
  //useful structures for calculation/ readability
  std::vector<TLorentzVector*> jets{};
  TLorentzVector* Wl; //l+nu
  TLorentzVector* Wqq; //q+q
  TLorentzVector* nu; //made from missing p with m=0
  std::vector<TLorentzVector*> CMJets{}; //q,q,l boosted into W rest frame
  TLorentzVector* CMnu;//nu boosted into W restframe

	int   _printing{};

  //input collections
  std::string _inputMcParticleCollectionName{};
  std::string _inputJetCollectionName{};


  /* histograms split between muon/tau true events */
	TFile* file;

	TH1D *WmassMuon[ncuts+1], *WmassTau[ncuts+1], qqmassMuon[ncuts+1], qqmassTau[ncuts+1];
	TH1D *WEMuon[ncuts+1], *WETau[ncuts+1], *EtotalMuon[ncuts+1], *EtotalTau[ncuts+1];
	TH1D *Wm_cosTheta[ncuts+1];

	TH1D *LjetMassMuon[ncuts+1], *LjetMassTau[ncuts+1];

	//tgc hists
	TH1D *costhetawMuon[ncuts+1] , *costhetawTau[ncuts+1];
	TH1D *thetaLMuon[ncuts+1], *thetaLTau[ncuts+1];
	TH1D *phiLMuon[ncuts+1], *phiLTau[ncuts+1];
	TH1D *thetaHMuon[ncuts+1], *thetaHTau[ncuts+1];
	TH1D *phiHMuon[ncuts+1], *phiHTau[ncuts+1];
	
 	/* end histograms */

};

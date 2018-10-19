#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/LCCollectionVec.h"
#include <marlin/Global.h>
#include "gear/BField.h"
#include "lcio.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
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
  int identifyLeptonJet_bySeparation(std::vector<ReconstructedParticle*> jets);
  void getAngleOfljetandMCLepton();
  void classifyTauDecay(MCParticle* mctau);

  int getLeptonJetCharge( ReconstructedParticle* ljet );

  //jet analysis helpers
  void getJetMultiplicities();
  void exploreDaughterParticles(MCParticle* p, std::vector<MCParticle*>& FSP);
  bool allChildrenAreSimulation(MCParticle* p);
  void analyzeLeadingTracks();

  //classify the type of lepton decay and retrieve the
  //mcparticles for qqlnu
//  MCParticle* classifyEvent(bool& isTau, bool& isMuon, int& trueq);
  MCParticle* classifyEvent(bool& isTau, bool& isMuon, int& trueq, TLorentzVector* (&_MCf)[4], int (&_MCfpdg)[4]);
//  MCParticle* classifyEvent(bool& isTau, bool& isMuon, int& trueq, int (&_MCfpdg)[4]);

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

//TTree
  TTree* _tree;
  int _nRun;
  int _nEvt;
//  float _xsec;
//  TString *_Process;

  TLorentzVector* _MCf[4];
  int _MCfpdg[4];

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


  //jet matching and jet multiplicity variables
  int lpdg; // true pdg code for the lepton
  int lnparts; // number of particles in lepton jet
  int lntracks; // number of tracks in lepton jet
  int lnmcparts; //true n daughters of lepton
  int lnmctracks; //true n daughter tracks of lepton
  std::vector<int> jetNparts; //number of particles per any jet
  std::vector<int> jetNtracks; //number of tracks per any jet
  double leadingptljet; //pt of the leading track in the lepton jet
  double leadingd0ljet; //d0 of the leading track in the lepton jet
  double leadingd0relerrljet; //relative error of d0 of leading track in lepton jet

  double leadingptqjet; //pt of the leading track in a quark jet
  double leadingd0qjet; //d0 of the leading track in a quark jet
  double leadingd0relerrqjet; //relative error of d0 of leading track in lepton jet

  //opening angle between the lepton jet and mc lepton
  double psi_mcl_ljet;
  int ljetmatchmctau=0;
  int ljetmatchmcmuon=0;

  int qnparts;
  int qntracks;
  int qmcparts;
  int qmctracks;


	int   _printing{};

  //input collections
  std::string _inputMcParticleCollectionName{};
  std::string _inputJetCollectionName{};


  /* histograms split between muon/tau true events */
	TFile* file;

	TH1D *WmassMuon[ncuts+1], *WmassTau[ncuts+1], *qqmassMuon[ncuts+1], *qqmassTau[ncuts+1];
	TH1D *WEMuon[ncuts+1], *WETau[ncuts+1], *EtotalMuon[ncuts+1], *EtotalTau[ncuts+1];
	TH1D *Wm_cosTheta[ncuts+1];

	TH1D *LjetMassMuon[ncuts+1], *LjetMassTau[ncuts+1];

	//tgc hists
	TH1D *costhetawMuon[ncuts+1] , *costhetawTau[ncuts+1];
	TH1D *thetaLMuon[ncuts+1], *thetaLTau[ncuts+1];
	TH1D *phiLMuon[ncuts+1], *phiLTau[ncuts+1];
	TH1D *thetaHMuon[ncuts+1], *thetaHTau[ncuts+1];
	TH1D *phiHMuon[ncuts+1], *phiHTau[ncuts+1];

    //jet information histograms
    TH1D  *leptonMCNPartsMuon[ncuts+1], *leptonMCNTracksMuon[ncuts+1], *leptonMCNPartsTau[ncuts+1], *leptonMCNTracksTau[ncuts+1];
	TH1D  *jetNpartsMuon[ncuts+1], *minjetNpartsMuon[ncuts+1], *jetNpartsTau[ncuts+1], *minjetNpartsTau[ncuts+1];
    TH1D  *jetNtracksMuon[ncuts+1], *minjetNtracksMuon[ncuts+1],  *jetNtracksTau[ncuts+1], *minjetNtracksTau[ncuts+1];

    //lepton jet info
	TH1D *ljetleadingd0Muon[ncuts+1], *ljetleadingd0Tau[ncuts+1], *ljetleadingptMuon[ncuts+1], *ljetleadingptTau[ncuts+1];
    TH1D *ljetd0relerrMuon[ncuts+1], *ljetd0relerrTau[ncuts+1]; 
    TH1D *qjetleadingd0Muon[ncuts+1], *qjetleadingd0Tau[ncuts+1], *qjetleadingptMuon[ncuts+1], *qjetleadingptTau[ncuts+1];
    TH1D *qjetd0relerrMuon[ncuts+1], *qjetd0relerrTau[ncuts+1];

    TH1D *psiljetmclMuon[ncuts+1], *psiljetmclTau[ncuts+1];	
 	/* end histograms */

};

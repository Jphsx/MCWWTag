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

#include <iostream>
#include <fstream>

using namespace lcio;

	/** MCWWTag:<br>
 *
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */

 class MCWWTag : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new MCWWTag ; }

  MCWWTag(const MCWWTag&) = delete ;
  MCWWTag& operator=(const MCWWTag&) = delete ;

  MCWWTag() ;

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

  bool FindMCParticles( LCEvent* evt );
  bool FindJets( LCEvent* evt ) ;

  int identifyLeptonJet( std::vector<ReconstructedParticle*> jets);
  int getLeptonJetCharge( ReconstructedParticle* ljet );

  protected:
  int nEvt{};
  int ntau=0;
  int nmuon=0;
  int nelec=0;

   int nup=0;
   int ndwn=0;
   int nstr=0;
   int nchm=0;

  //how many times do we get the proper lepton charge?
  //for muons and for leptons separately
  int muonqmatch=0;
  int tauqmatch=0;
  
  //vector to hold the particles for the event
  std::vector<MCParticle*> _mcpartvec{};
  std::vector<ReconstructedParticle*> _jets{};
  int   _printing{};

  std::string _inputMcParticleCollectionName{};
  std::string _inputJetCollectionName{};


  /* histograms */
	TFile* file;

	TH1D* WmassMuon;
	TH1D* WmassTau;
	TH1D* WEMuon;
	TH1D* WETau;
	TH1D* Wm_cosTheta;

	TH1D* LjetMassMuon;
	TH1D* LjetMassTau;







 /* end histograms */

};

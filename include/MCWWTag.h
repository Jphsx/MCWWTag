#include "marlin/Processor.h"
#include "EVENT/MCParticle.h"
#include "IMPL/LCCollectionVec.h"
#include <marlin/Global.h>
#include "gear/BField.h"
#include "lcio.h"
#include "TFile.h"
#include <vector>

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

  protected:
  int nEvt{};
  
  //vector to hold the tracks for the event
  std::vector<MCParticle*> _mcpartvec{};
  int   _printing{};

  std::string _inputMcParticleCollectionName{};

};

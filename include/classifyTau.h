#ifndef _TAUDECAY_
#define _TAUDECAY_

#include "EVENT/MCParticle.h"
#include "EVENT/ReconstructedParticle.h"
#include "TString.h"
#include <iostream>
using namespace lcio;
class classifyTau {

	public:
		
		enum { decayChPi=0, decayRho, decayA1_1p, decayA1_3p , decayEl, decayMu , decayW , decayK , decayMultiprong , decayOthersingleprong, decayUnrecognised, NDECAYS};

	 static std::vector < EVENT::MCParticle* > getstablemctauDaughters( EVENT::MCParticle* mctau );
  static std::vector < EVENT::MCParticle* > getmctauDaughters( EVENT::MCParticle* mctau, bool onlyStable=false );

static int getMCdecayMode( std::vector < EVENT::MCParticle* > mcps );
  static int getMCdecayMode( EVENT::MCParticle* mcp ) {
    std::vector < EVENT::MCParticle* > sd = getstablemctauDaughters( mcp );
    return getMCdecayMode ( sd );
  }

  static TString getTauDecLab( int code );
	private:

} 

#endif

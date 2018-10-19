

#include "classifyTau.h"

std::vector < MCParticle* > classifyTau::getstablemctauDaughters( MCParticle* mctau ) {
  return getmctauDaughters ( mctau, true );
}


std::vector < MCParticle* > classifyTau::getmctauDaughters( MCParticle* mctau, bool onlystable ) {
  assert( abs( mctau->getPDG())==15 );
  std::vector < MCParticle* > daughters;
  std::vector < MCParticle* > pp;
  pp.push_back(mctau);
  while ( pp.size()>0 ) {
    std::vector < MCParticle* > pp2;
    for ( size_t i=0; i< pp.size(); i++) {
      if ( pp[i]->getGeneratorStatus() == 1 ) {
	daughters.push_back( pp[i] );
      } else {
	if ( !onlystable ) daughters.push_back( pp[i] );
	for ( size_t k=0; k< pp[i]->getDaughters().size(); k++) {
	  pp2.push_back( pp[i]->getDaughters()[k] );
	}
      }
    }
    pp=pp2;
  }
  return daughters;
}
int classifyTau::getMCdecayMode( std::vector < MCParticle* > mcps ) {

  int npichg(0);
  int nel(0);
  int nmu(0);
  int nnu(0);
  int ngamma(0);
  int nK0L(0);
  int nKchg(0);

  int nChg(0);
  int nNeu(0);

  for ( size_t i=0; i<mcps.size(); i++) {
    switch ( abs( mcps[i]->getPDG() ) ) {
    case 11:
      nel++; 
      nChg++;
      break;
    case 13:
      nmu++; 
      nChg++;
      break;
    case 12:
    case 14:
    case 16:
      nnu++; 
      break;
    case 211:
      npichg++; 
      nChg++;
      break;
    case 22:
      ngamma++; 
      nNeu++;
      break;
    case 130:
      nK0L++;
      nNeu++;
      break;
    case 321:
      nKchg++; 
      nChg++;
      break;
    default:
      cout << "unknown pdg " << mcps[i]->getPDG() << endl;
      assert(0);
    }
  }



  int decayMode(-1);
  if ( nel==1 ) decayMode = decayEl;
  else if ( nmu==1 ) decayMode = decayMu;
  else if ( npichg==nChg && npichg==1 && nNeu==0 ) decayMode = decayChPi;
  else if ( npichg==nChg && npichg==1 && ngamma==2 && nNeu==2 ) decayMode = decayRho;
  else if ( npichg==nChg && npichg==3 && nNeu==0 ) decayMode = decayA1_3p;
  else if ( npichg==nChg && npichg==1 && ngamma==nNeu && ngamma==4 ) decayMode = decayA1_1p;
  else if ( nK0L>0 || nKchg>0 ) decayMode = decayK;
  else if ( npichg+nKchg>0 ) decayMode = decayW; // everything else...
  else {
    cout << "unrecognised tau decay!" << endl;
    cout << " pi " <<  npichg;
    cout << " e " <<  nel;
    cout << " m " <<  nmu;
    cout << " nu " <<  nnu;
    cout << " g " <<  ngamma;
    cout << " kl " <<  nK0L;
    cout << " k " <<  nKchg;
    cout << " chg " <<  nChg;
    cout << " neu " <<  nNeu;
    cout << endl;    
    assert(0);
  }

  return decayMode;
}

TString classifyTau::getTauDecLab( int code ) {
  //  enum { decayChPi=0, decayRho, decayA1_3p, decayA1_1p , decayEl, decayMu , decayW , decayK , decayMultiprong , decayOthersingleprong, decayUnrecognised};
  switch (code) {
  case decayChPi:
    return "PI";
  case decayRho:
    return "RHO";
  case decayA1_3p:
    return "A3P";
  case decayA1_1p:
    return "A1P";
  case decayEl:
    return "E";
  case decayMu:
    return "M";
  case decayW:
    return "HAD";
  case decayK:
    return "K";
  case decayMultiprong:
    return "MP";
  case decayOthersingleprong:
    return "SP";
  case decayUnrecognised:
  default:
    return "OTH";
  }

}

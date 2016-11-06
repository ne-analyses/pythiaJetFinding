// testing jetfinding algorithms
// on pythia
// Nick Elsey

// ROOT is used for histograms
// ROOT Headers
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TObjArray.h"
#include "TString.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TBranch.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TSystem.h"

// My standard includes
// Make use of std::vector,
// std::string, IO and algorithm
// STL Headers
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <random>
#include <time.h>
#include <limits.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

// The analysis is run on FastJet::PseudoJets
// We make use of the jetfinding tools
// And the convenient FastJet::Selectors
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/ClusterSequencePassiveArea.hh"
#include "fastjet/ClusterSequenceActiveArea.hh"
#include "fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh"
#include "fastjet/Selector.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Filter.hh"

// Pythia generator
#include "Pythia8/Pythia.h"

void convertToPseudoJet( Pythia8::Pythia& p, double max_rap, std::vector<fastjet::PseudoJet>& all, std::vector<fastjet::PseudoJet>& charged, std::vector<fastjet::PseudoJet>& part ) {
  
  // clear the event containers
  all.clear();
  charged.clear();
  part.clear();
  
  // get partons first
  // the initial protons are ids 1 & 2,
  // the hard scattering products are ids 5 & 6
  if (p.event[5].status() != -23)
    std::cerr<<"Error: assumption that id 5 is the outgoing parton is not valid."<<std::endl;
  if (p.event[6].status() != -23)
    std::cerr<<"Error: assumption that id 6 is the outgoing parton is not valid."<<std::endl;
  fastjet::PseudoJet part1( p.event[5].px(), p.event[5].py(), p.event[5].pz(), p.event[5].e() );
  part1.set_user_index( 3 * p.event[5].charge() );
  fastjet::PseudoJet part2( p.event[6].px(), p.event[6].py(), p.event[6].pz(), p.event[6].e() );
  part2.set_user_index( 3 * p.event[6].charge() );
  part.push_back( part1 );
  part.push_back( part2 );
  
  // now loop over all particles, and fill the vectors
  for ( int i = 0; i < p.event.size(); ++i ) {
    if ( p.event[i].isFinal() && p.event[i].isVisible() ) {
      fastjet::PseudoJet tmp( p.event[i].px(), p.event[i].py(), p.event[i].pz(), p.event[i].e() );
      tmp.set_user_index( p.event[i].charge() );
      
      // check to make sure its in our rapidity range
      if ( fabs( tmp.rap() ) > max_rap  )
        continue;
      
      all.push_back( tmp );
      if ( p.event[i].charge() )
        charged.push_back( tmp );
      
      
    }
  }
  
}

// Arguments
// 0: xml directory for pythia
// 1: exponent base 10 for number of events
// 2: output location


int main( int argc, const char** argv ) {
  
  // set parameters
  unsigned exponent;
  std::string outFile;
  std::string xmldir;

  switch ( argc ) {
    case 1: {
      exponent = 4;
      outFile = "out/test.root";
      xmldir = "/Users/nick/physics/software/pythia8/share/Pythia8/xmldoc";
      break;
    }
    case 4: {
      xmldir = argv[1];
      exponent = atoi( argv[2] );
      outFile = argv[3];
      break;
    }
    default: {
      std::cerr<<"Error: unexpected number of inputs."<<std::endl;
      return -1;
    }
  }
  
  // set the total number of events as
  // 10^exponent
  unsigned maxEvent = pow( 10, exponent );
  std::cout<<"set for "<<maxEvent<<" events"<<std::endl;
  
  // setup pythia
  // ------------
  
  // create the pythia generator and initialize it with the xmldoc in my pythia8 directory
  Pythia8::Pythia pythia( xmldir );
  
  // settings for LHC pp at 13 TeV
  pythia.readString("Beams:eCM = 13000");
  pythia.readString("HardQCD:all = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0");
  pythia.readString("PhaseSpace:pTHatMin = 5.0");
  
  // initialize the pythia generator
  pythia.init();
  
  // set jet finding parameters
  // --------------------------
  
  // set a hard cut on rapidity for all tracks
  const double max_track_rap = 10;
  const double max_rap = max_track_rap;
  
  // first some base jetfinding definitions
  double baseRadius = 0.8;
  fastjet::JetDefinition antiKtBase( fastjet::antikt_algorithm, baseRadius );
  fastjet::JetDefinition KtBase( fastjet::kt_algorithm, baseRadius );
  fastjet::JetDefinition CaBase( fastjet::cambridge_algorithm, baseRadius );
  
  // but we will also be testing these with different radii, so we'll initialize that here
  // there will be nRadii different radii, in increments of deltaRad;
  int nRadii = 10;
  double deltaRad = 0.1;
  double radii[nRadii];
  fastjet::JetDefinition antiKtDefs[10];
  fastjet::JetDefinition KtDefs[10];
  fastjet::JetDefinition CaDefs[10];
  
  for ( int i = 1; i <= nRadii; ++i ) {
    radii[i-1] = deltaRad * i;
    antiKtDefs[i-1] = fastjet::JetDefinition( fastjet::antikt_algorithm, radii[i-1] );
    KtDefs[i-1] = fastjet::JetDefinition( fastjet::antikt_algorithm, radii[i-1] );
    CaDefs[i-1] = fastjet::JetDefinition( fastjet::cambridge_algorithm, radii[i-1] );
  }
  
  
  // set up our fastjet environment
  // ------------------------------
  
  // We'll be using fastjet:PseudoJet for all our work
  // so we'll make a few helpers
  
  // this will include all final state particles
  // ( minus neutrinos )
  std::vector<fastjet::PseudoJet> allFinal;
  
  // this will include all charged particles in the final state
  std::vector<fastjet::PseudoJet> chargedFinal;
  
  // this will be the two partons from the scattering
  std::vector<fastjet::PseudoJet> partons;
  
  // create an area definition for the clustering
  //----------------------------------------------------------
  // ghosts should go up to the acceptance of the detector or
  // (with infinite acceptance) at least 2R beyond the region
  // where you plan to investigate jets.
  const int ghost_repeat = 1;
  const double ghost_area = 0.01;
  // we'll set ghost_max_rap to the largest applicable based
  // on our radius settings
  const double ghost_max_rap = max_rap + 2.0 * radii[nRadii];
  
  fastjet::GhostedAreaSpec area_spec = fastjet::GhostedAreaSpec( ghost_max_rap, ghost_repeat, ghost_area );
  fastjet::AreaDefinition  area_def = fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, area_spec);
  
  
  // create output histograms using root
  TH1D* multiplicity = new TH1D("mult", "Visible Multiplicity", 300, -0.5, 899.5 );
  TH1D* chargedMultiplicity = new TH1D("chargemult", "Charged Multiplicity", 300, -0.5, 899.5 );
  TH1D* partonPt = new TH1D("partonpt", "Parton Pt", 100, 0, 1000 );
  TH1D* partonE = new TH1D( "parton_e", "Parton Energy", 100, 0, 1000 );
  TH2D* partonEtaPhi = new TH2D("partonetaphi", "Parton Eta x Phi", 100, -12, 12, 100, -TMath::Pi(), TMath::Pi() );
  
  // associated particle information
  TH1D* visiblePt = new TH1D( "finalstatept", "Detected Pt", 200, 0, 100 );
  TH1D* visibleE = new TH1D( "finalstateE", "Detected E", 200, 0, 100 );
  TH2D* visibleEtaPhi = new TH2D( "finaletaphi", "Detected Eta x Phi",  100, -12, 12, 100, -TMath::Pi(), TMath::Pi() );
  TH1D* chargedPt = new TH1D("chargedfstatept", "Detected Charged Pt", 200, 0, 100);
  TH1D* chargedE = new TH1D( "chargedfstateE", "Detected Charged E", 200, 0, 100 );
  TH2D* chargedEtaPhi = new TH2D( "chargedetaphi", "Detected Charged Eta x Phi",  100, -12, 12, 100, -TMath::Pi(), TMath::Pi() );
  
  // jet information
  TH1D* nJetsAntiKtBaseAll = new TH1D("njetsantiktbase", "Jet Multiplicity Anti-Kt Base", 50, 49.5, 299.5);
  TH1D* nJetsKtBaseAll = new TH1D("njetsktbase", "Jet Multiplicity Kt Base", 50, 49.5, 299.5);
  TH1D* nJetsCaBaseAll = new TH1D("njetsCabase", "Jet Multiplicity CA Base", 50, 49.5, 299.5);
  TH1D* nJetsAntiKtBaseCharged = new TH1D("njetsantiktbasecharged", "Jet Multiplicity Anti-Kt Base Charged", 50, 49.5, 299.5);
  TH1D* nJetsKtBaseCharged = new TH1D("njetsktbasecharged", "Jet Multiplicity Kt Base Charged", 50, 49.5, 299.5);
  TH1D* nJetsCaBaseCharged = new TH1D("njetsCabasecharged", "Jet Multiplicity CA Base Charged", 50, 49.5, 299.5);
  
  // start the event loop from event 0
  unsigned currentEvent = 0;
  try{
    while ( currentEvent < maxEvent ) {
      
      // try to generate a new event
      // if it fails, iterate without incrementing
      // current event number
      if ( !pythia.next() )
        continue;
      
      // pythia succeeded, so increment the event
      currentEvent++;
      
      // convert pythia particles into useable pseudojets,
      // only take those in our eta range && that are visible
      // in conventional detectors
      // note: particles user_index() is the charge
      convertToPseudoJet( pythia, max_track_rap, allFinal, chargedFinal, partons );
      
      // event information
      multiplicity->Fill( allFinal.size() );
      chargedMultiplicity->Fill( chargedFinal.size() );
      
      // fill parton information
      for ( int i = 0; i < 2; ++i ) {
        partonEtaPhi->Fill( partons[i].eta(), partons[i].phi_std() );
        partonPt->Fill( partons[i].pt() );
        partonE->Fill( partons[i].E() );
      }
      
      // now fill track information
      for ( int i = 0; i < allFinal.size(); ++i ) {
        visiblePt->Fill( allFinal[i].pt() );
        visibleE->Fill( allFinal[i].E() );
        visibleEtaPhi->Fill( allFinal[i].eta(), allFinal[i].phi_std() );
      }
      
      for ( int i = 0; i < chargedFinal.size(); ++i ) {
        chargedPt->Fill( chargedFinal[i].pt() );
        chargedE->Fill( chargedFinal[i].E() );
        chargedEtaPhi->Fill( chargedFinal[i].eta(), chargedFinal[i].phi_std() );
      }
      
      // now set up the clustering
      fastjet::ClusterSequenceArea clusterAntiKtAll ( allFinal, antiKtBase, area_def );
      std::vector<fastjet::PseudoJet> antiKtBaseJets = clusterAntiKtAll.inclusive_jets();
      fastjet::ClusterSequenceArea clusterKtAll ( allFinal, KtBase, area_def );
      std::vector<fastjet::PseudoJet> KtBaseJets = clusterKtAll.inclusive_jets();
      fastjet::ClusterSequenceArea clusterCaAll ( allFinal, CaBase, area_def );
      std::vector<fastjet::PseudoJet> CaBaseJets = clusterCaAll.inclusive_jets();
      fastjet::ClusterSequenceArea clusterAntiKtCharged ( chargedFinal, antiKtBase, area_def );
      std::vector<fastjet::PseudoJet> antiKtChargedJets = clusterAntiKtCharged.inclusive_jets();
      fastjet::ClusterSequenceArea clusterKtCharged ( chargedFinal, KtBase, area_def );
      std::vector<fastjet::PseudoJet> KtChargedJets = clusterKtCharged.inclusive_jets();
      fastjet::ClusterSequenceArea clusterCaCharged ( chargedFinal, CaBase, area_def );
      std::vector<fastjet::PseudoJet> CaChargedJets = clusterCaCharged.inclusive_jets();
      
      
      // plot the number of jets by the different algorithms
      nJetsAntiKtBaseAll->Fill( antiKtBaseJets.size() );
      nJetsKtBaseAll->Fill( KtBaseJets.size() );
      nJetsCaBaseAll->Fill( CaBaseJets.size() );
      nJetsAntiKtBaseCharged->Fill( antiKtChargedJets.size() );
      nJetsKtBaseCharged->Fill( KtChargedJets.size() );
      nJetsCaBaseCharged->Fill( CaChargedJets.size() );
      
    }
  } catch ( std::exception& e) {
    std::cerr << "Caught " << e.what() << std::endl;
    return -1;
  }
  std::cout<<"processed "<<currentEvent<<" events"<<std::endl;
  
  // print out pythia statistics
  pythia.stat();
  
  // write out to a root file all histograms
  TFile out( outFile.c_str(), "RECREATE" );
  
  // event information
  multiplicity->Write();
  chargedMultiplicity->Write();
  
  // parton information
  partonEtaPhi->Write();
  partonPt->Write();
  partonE->Write();
  
  // track information
  visiblePt->Write();
  visibleE->Write();
  visibleEtaPhi->Write();
  chargedPt->Write();
  chargedE->Write();
  chargedEtaPhi->Write();
  
  // jet multiplicity
  nJetsAntiKtBaseAll->Write();
  nJetsKtBaseAll->Write();
  nJetsCaBaseAll->Write();
  nJetsAntiKtBaseCharged->Write();
  nJetsKtBaseCharged->Write();
  nJetsCaBaseCharged->Write();
  
  
  // close the output file
  out.Close();
  
  return 0;
}




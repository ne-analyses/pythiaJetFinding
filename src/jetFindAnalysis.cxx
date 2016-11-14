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

// the grid does not have std::to_string() for some ungodly reason
// replacing it here. Simply ostringstream
namespace patch {
  template < typename T > std::string to_string( const T& n )
  {
    std::ostringstream stm ;
    stm << n ;
    return stm.str() ;
  }
}

// used to convert pythia events to vectors of pseudojets
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
      //xmldir = "/Users/nick/physics/software/pythia8/share/Pythia8/xmldoc";
      xmldir = "/wsu/home/dx/dx54/dx5412/software/pythia8219/share/Pythia8/xmldoc";
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
  //Pythia8::Pythia pythia( xmldir );
  Pythia8::Pythia pythia;
  
  // settings for LHC pp at 13 TeV
  //pythia.readString("Beams:eCM = 13000");
  pythia.readString("HardQCD:all = on");
  //pythia.readString("Random:setSeed = on");
  //pythia.readString("Random:seed = 0");
  //pythia.readString("PhaseSpace:pTHatMin = 5.0");
  
  // initialize the pythia generator
  pythia.init();
  pythia.next()
  
  /*
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
  
  for ( int i = 0; i <= nRadii; ++i ) {
    radii[i] = deltaRad * (i+1);
    antiKtDefs[i] = fastjet::JetDefinition( fastjet::antikt_algorithm, radii[i] );
    KtDefs[i] = fastjet::JetDefinition( fastjet::antikt_algorithm, radii[i] );
    CaDefs[i] = fastjet::JetDefinition( fastjet::cambridge_algorithm, radii[i] );
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
  
  // make a histogram for all of the differing radii
  TH2D* nJetsAntiKt = new TH2D( "njetsantikt", "Number of Jets - Anti-Kt", nRadii, -0.5, nRadii-0.5, 200, -0.5, 399.5 );
  TH2D* deltaEAntiKt = new TH2D( "deltaEantikt", "#Delta E - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, -100, 100 );
  TH2D* deltaRAntiKt = new TH2D( "deltaRantikt", "#Delta R Leading - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 1.0 );
  TH2D* nPartAntiKt = new TH2D( "npartantikt", "Number of Particles per Jet - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* nPartLeadAntiKt = new TH2D( "npartleadantikt", "Number of Particles per Leading Jet - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* timeAntiKt = new TH2D("clustertimeantikt", "Time Required to cluster - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, -0.5, 99.5);
  TH2D* areaAntiKt = new TH2D("areaantikt", "Jet Area - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );
  TH2D* areaLeadAntiKt = new TH2D("arealeadantikt", "Lead Jet Area - Anti-Kt", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );
  
  TH2D* nJetsKt = new TH2D( "njetskt", "Number of Jets - Kt", nRadii, -0.5, nRadii-0.5, 200, -0.5, 399.5 );;
  TH2D* deltaEKt = new TH2D( "deltaEkt", "#Delta E - Kt", nRadii, -0.5, nRadii-0.5, 100, -100, 100 );
  TH2D* deltaRKt = new TH2D( "deltaRkt", "#Delta R - Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 1.0 );
  TH2D* nPartKt = new TH2D( "npartkt", "Number of Particles per Jet - Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* nPartLeadKt = new TH2D( "npartleadKt", "Number of Particles per Leading Jet - Kt", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* timeKt = new TH2D("clustertimekt", "Time Required to cluster - Kt", nRadii, -0.5, nRadii-0.5, 100, -0.5, 99.5);
  TH2D* areaKt = new TH2D("areakt", "Jet Area - Kt", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );
  TH2D* areaLeadKt = new TH2D("arealeadkt", "Lead Jet Area - Kt", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );
  
  TH2D* nJetsCa = new TH2D( "njetsca", "Number of Jets - CA", nRadii, -0.5, nRadii-0.5, 200, -0.5, 399.5 );
  TH2D* deltaECa = new TH2D( "deltaEca", "#Delta E - CA", nRadii, -0.5, nRadii-0.5, 100, -100, 100 );
  TH2D* deltaRCa = new TH2D( "deltaRca", "#Delta R Leading - CA", nRadii, -0.5, nRadii-0.5, 100, 0, 1.0 );
  TH2D* nPartCa = new TH2D( "npartca", "Number of Particles per Jet - CA", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* nPartLeadCa = new TH2D( "npartleadca", "Number of Particles per Leading Jet - CA", nRadii, -0.5, nRadii-0.5, 100, 0, 100 );
  TH2D* timeCa = new TH2D("clustertimeca", "Time Required to cluster - CA", nRadii, -0.5, nRadii-0.5, 100, -0.5, 99.5);
  TH2D* areaCa = new TH2D("areaca", "Jet Area - CA", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );
  TH2D* areaLeadCa = new TH2D("arealeadca", "Lead Jet Area - CA", nRadii, -0.5, nRadii-0.5, 100, 0, TMath::Pi() );

  // set bin labels to radii
  for ( int i = 1; i <= nRadii; ++i ) {

    nJetsAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaEAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaRAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartLeadAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    timeAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1]).c_str() );
    areaAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1]).c_str() );
    areaLeadAntiKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1]).c_str() );
    
    nJetsKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaEKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaRKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartLeadKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    timeKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    areaKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    areaLeadKt->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    
    nJetsCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaECa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    deltaRCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    nPartLeadCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    timeCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    areaCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
    areaLeadCa->GetXaxis()->SetBinLabel( i, patch::to_string( radii[i-1] ).c_str() );
  }
  
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
      std::vector<fastjet::PseudoJet> antiKtBaseJets = fastjet::sorted_by_pt(clusterAntiKtAll.inclusive_jets());
      fastjet::ClusterSequenceArea clusterKtAll ( allFinal, KtBase, area_def );
      std::vector<fastjet::PseudoJet> KtBaseJets = fastjet::sorted_by_pt(clusterKtAll.inclusive_jets());
      fastjet::ClusterSequenceArea clusterCaAll ( allFinal, CaBase, area_def );
      std::vector<fastjet::PseudoJet> CaBaseJets = fastjet::sorted_by_pt(clusterCaAll.inclusive_jets());
      fastjet::ClusterSequenceArea clusterAntiKtCharged ( chargedFinal, antiKtBase, area_def );
      std::vector<fastjet::PseudoJet> antiKtChargedJets = fastjet::sorted_by_pt(clusterAntiKtCharged.inclusive_jets());
      fastjet::ClusterSequenceArea clusterKtCharged ( chargedFinal, KtBase, area_def );
      std::vector<fastjet::PseudoJet> KtChargedJets = fastjet::sorted_by_pt(clusterKtCharged.inclusive_jets());
      fastjet::ClusterSequenceArea clusterCaCharged ( chargedFinal, CaBase, area_def );
      std::vector<fastjet::PseudoJet> CaChargedJets = fastjet::sorted_by_pt(clusterCaCharged.inclusive_jets());

      // plot the number of jets by the different algorithms
      nJetsAntiKtBaseAll->Fill( antiKtBaseJets.size() );
      nJetsKtBaseAll->Fill( KtBaseJets.size() );
      nJetsCaBaseAll->Fill( CaBaseJets.size() );
      nJetsAntiKtBaseCharged->Fill( antiKtChargedJets.size() );
      nJetsKtBaseCharged->Fill( KtChargedJets.size() );
      nJetsCaBaseCharged->Fill( CaChargedJets.size() );

      // now we'll do the loop over differing radii
      for ( int i = 0; i < nRadii; ++i ) {

        std::string radBin = patch::to_string( radii[i] );
        
        // first perform the clustering
        
        // time the clustering as well
        clock_t t = clock();
        fastjet::ClusterSequenceArea clusterAntiKt( allFinal, antiKtDefs[i], area_def );
        clock_t antiKtTime = clock() - t;
        t = clock();
        fastjet::ClusterSequenceArea clusterKt( allFinal, KtDefs[i], area_def );
        clock_t ktTime = clock() - t;
        t = clock();
        fastjet::ClusterSequenceArea clusterCa( allFinal, CaDefs[i], area_def );
        clock_t caTime = clock() - t;
        
        // fill timing measurements
        timeAntiKt->Fill( radBin.c_str(), antiKtTime, 1 );
        timeKt->Fill( radBin.c_str(), ktTime, 1 );
        timeCa->Fill( radBin.c_str(), caTime, 1 );
        
        std::vector<fastjet::PseudoJet> antiKtJets = fastjet::sorted_by_pt( clusterAntiKt.inclusive_jets() );
        std::vector<fastjet::PseudoJet> KtJets = fastjet::sorted_by_pt( clusterKt.inclusive_jets() );
        std::vector<fastjet::PseudoJet> CaJets = fastjet::sorted_by_pt( clusterCa.inclusive_jets() );
        
        // now start to fill histograms
        // first, number of jets in the event
        nJetsAntiKt->Fill ( radBin.c_str(), antiKtJets.size(), 1 );
        nJetsKt->Fill ( radBin.c_str(), KtJets.size(), 1 );
        nJetsCa->Fill ( radBin.c_str(), CaJets.size(), 1 );
        
        // now, we'll do number of particles, and area, for both both leading jets and inclusive jets
        nPartLeadAntiKt->Fill ( radBin.c_str(), antiKtJets[0].constituents().size(), 1 );
        nPartLeadKt->Fill ( radBin.c_str(), KtJets[0].constituents().size(), 1 );
        nPartLeadCa->Fill ( radBin.c_str(), CaJets[0].constituents().size(), 1 );
        areaLeadAntiKt->Fill ( radBin.c_str(), antiKtJets[0].area(), 1 );
        areaLeadKt->Fill ( radBin.c_str(), KtJets[0].area(), 1 );
        areaLeadCa->Fill ( radBin.c_str(), CaJets[0].area(), 1 );
        
        for ( int j = 0; j < antiKtJets.size(); ++j ) {
          nPartAntiKt->Fill ( radBin.c_str(), antiKtJets[j].constituents().size(), 1 );
          areaAntiKt->Fill ( radBin.c_str(), antiKtJets[j].area(), 1 );
        }
        for ( int j = 0; j < KtJets.size(); ++j ) {
          nPartKt->Fill ( radBin.c_str(), KtJets[j].constituents().size(), 1 );
          areaKt->Fill ( radBin.c_str(), KtJets[j].area(), 1 );
        }
        for ( int j = 0; j < CaJets.size(); ++j ) {
          nPartCa->Fill ( radBin.c_str(), CaJets[j].constituents().size(), 1 );
          areaCa->Fill ( radBin.c_str(), CaJets[j].area(), 1 );
        }
        
        // and compare to the initial partons for delta E and delta R
        // we find the minimum of the delta R between leading jet and parton1 and parton2
        // and use that as the base for both delta R and delta E
        
        // first antikt
        double distToPart1 = partons[0].delta_R(antiKtJets[0]);
        double distToPart2 = partons[1].delta_R(antiKtJets[0]);
        int partonIdx = 0;
        if ( distToPart2 < distToPart1 )
          partonIdx = 1;
        deltaRAntiKt->Fill ( radBin.c_str(), partons[partonIdx].delta_R( antiKtJets[0] ), 1 );
        deltaEAntiKt->Fill ( radBin.c_str(), partons[partonIdx].E() - antiKtJets[0].E(), 1 );
        
        // repeat for Kt and Ca
        distToPart1 = partons[0].delta_R( KtJets[0] );
        distToPart2 = partons[1].delta_R( KtJets[0] );
        partonIdx = 0;
        if ( distToPart2 < distToPart1 )
          partonIdx = 1;
        deltaRKt->Fill ( radBin.c_str(), partons[partonIdx].delta_R( KtJets[0] ), 1 );
        deltaEKt->Fill ( radBin.c_str(), partons[partonIdx].E() - KtJets[0].E(), 1 );
        
        distToPart1 = partons[0].delta_R( CaJets[0] );
        distToPart2 = partons[1].delta_R( CaJets[0] );
        partonIdx = 0;
        if ( distToPart2 < distToPart1 )
          partonIdx = 1;
        deltaRCa->Fill ( radBin.c_str(), partons[partonIdx].delta_R( CaJets[0] ), 1 );
        deltaECa->Fill ( radBin.c_str(), partons[partonIdx].E() - CaJets[0].E(), 1 );
        
      }
      
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
  
  // histograms for differing radii
  nJetsAntiKt->Write();
  nPartAntiKt->Write();
  nPartLeadAntiKt->Write();
  deltaEAntiKt->Write();
  deltaRAntiKt->Write();
  timeAntiKt->Write();
  areaAntiKt->Write();
  areaLeadAntiKt->Write();
  
  nJetsKt->Write();
  nPartKt->Write();
  nPartLeadKt->Write();
  deltaEKt->Write();
  deltaRKt->Write();
  timeKt->Write();
  areaKt->Write();
  areaLeadKt->Write();
  
  nJetsCa->Write();
  nPartCa->Write();
  nPartLeadCa->Write();
  deltaECa->Write();
  deltaRCa->Write();
  timeCa->Write();
  areaCa->Write();
  areaLeadCa->Write();
  
  
  // close the output file
  out.Close();
  */
  return 0;
}




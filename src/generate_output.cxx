// generates 1D plots from the collected data
// from jetFindAnalysis.cxx

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
#include "TStyle.h"
#include "TLegend.h"

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
#include <chrono>
#include <thread>

// only one argument: the input file
// [1]: root file for input


int main ( int argc, const char** argv ) {
  
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  
  // get the input file name
  unsigned exponent;
  std::string inFile;
  
  switch ( argc ) {
    case 1: {
      inFile = "out/addedpythia.root";
      break;
    }
    case 2: {
      inFile = argv[1];
      break;
    }
    default: {
      std::cerr<<"Error: unexpected number of inputs."<<std::endl;
      return -1;
    }
  }

  // load the root file where the histograms are stored
  TFile rootFile( inFile.c_str(), "READ" );
  
  // Current histograms
  // ------------------
  
  // jetfinder names will be combined with what is plotted to
  // give the histogram name
  
  // first, jetfinder names
  const unsigned nJetFinders = 4;
  std::string jfNames[nJetFinders] = { "antikt", "kt", "ca", "sis" };
  std::string jfString[nJetFinders] = { "Anti-Kt", "Kt", "Cambridge-Aachen", "SISCone" };
  
  // now, the histogram names
  const unsigned nHistograms = 14;
  std::string histNames[nHistograms] = { "njets", "deltaE", "deltaR", "npart", "npartlead",
    "clustertime", "area", "arealead", "ptlead", "elead", "eta", "phi", "etalead", "philead" };
  
  // and the relevant jetfinding radii
  const unsigned nRadii = 10;
  std::string radii[nRadii] = { "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0" };
  unsigned baseRad = 7;
  
  // store the histograms in arrays of TH2Ds
  TH2D* histograms[nJetFinders][nHistograms];
  
  // load histograms from file
  for ( int i = 0; i < nHistograms; ++i ) {
    for ( int j = 0; j < nJetFinders; ++j ) {
      histograms[j][i] = (TH2D*) rootFile.Get( (jfNames[j]+histNames[i]).c_str() );
    }
  }
  
  // Now convert 2D histograms into 1D histograms
  TH1D* hist1D[nJetFinders][nHistograms][nRadii];

  
  for ( int i = 0; i < nJetFinders; ++i ) {
    for ( int j = 0; j < nHistograms; ++j ) {
      for ( int k = 0; k < nRadii; ++k ) {
        
        histograms[i][j]->GetXaxis()->SetRange( (k+1),(k+1) );
        hist1D[i][j][k] = (TH1D*) histograms[i][j]->ProjectionY();
        std::string name = jfNames[i] + histNames[j] + radii[k];
        hist1D[i][j][k]->SetName( name.c_str() );
        
      }
    }
  }
  

  
  // first produce measures of number of jets
  TCanvas* c1 = new TCanvas();
  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  for ( int i = 0; i < nJetFinders; ++ i ) {
    hist1D[i][0][baseRad]->SetTitle("Number of Jets");
    hist1D[i][0][baseRad]->GetXaxis()->SetTitle("Jets per Event");
    hist1D[i][0][baseRad]->GetYaxis()->SetTitle("Count");
    hist1D[i][0][baseRad]->SetLineColor(i);
    hist1D[i][0][baseRad]->SetLineWidth(2);
    hist1D[i][0][baseRad]->SetMarkerStyle(20+i);
    hist1D[i][0][baseRad]->SetMarkerColor(i);
    
    leg->AddEntry( hist1D[i][0][baseRad], jfString[i].c_str(), "lep"  );
    if ( i == 0 ) {
      hist1D[i][0][baseRad]->Draw();
    }
    else {
      hist1D[i][0][baseRad]->Draw("SAME");
    }
    
  }
  leg->Draw();
  
  c1->SaveAs("test.pdf");
  
  return 0;
}




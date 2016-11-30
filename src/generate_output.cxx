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
  const unsigned nJetFinders = 3;
  std::string jfNames[nJetFinders] = { "antikt", "kt", "ca" };
  
  // now, the histogram names
  const unsigned nHistograms = 14;
  std::string histNames[nHistograms] = { "njets", "deltaE", "deltaR", "npart", "npartlead",
    "clustertime", "area", "arealead", "ptlead", "elead", "eta", "phi", "etalead", "philead" };
  
  // store the histograms in arrays of TH2Ds
  TH2D* histograms[nJetFinders][nHistograms];
  
  // load histograms from file
  for ( int i = 0; i < nHistograms; ++i ) {
    for ( int j = 0; j < nJetFinders; ++j ) {
      histograms[j][i] = (TH2D*) rootFile.Get( (jfNames[j]+histNames[i]).c_str() );
    }
  }
  
  for ( int i = 0; i < nJetFinders; ++i ) {
    for ( int j = 0; j < nHistograms; ++j ) {
      std::cout<<histograms[i][j]<<std::endl;
    }
  }
  
  return 0;
}




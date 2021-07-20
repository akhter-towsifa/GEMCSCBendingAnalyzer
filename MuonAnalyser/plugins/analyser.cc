#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h>
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"


using namespace std;
using namespace edm;


struct MuonData
{
  void init();
  TTree* book(TTree *t);

  int muon_charge;
  float muon_pt;
  float muon_eta;
  float muon_momentum;

  //Shared prop
  int prop_location[5];
 
  //Prop from tracker
  float prop_inner_GP_GE11[3];
  float prop_inner_LP_GE11[3];
  float prop_inner_GP_startingPoint[3];
  float prop_inner_y_roll_GE11;
  float prop_inner_localphi_rad_GE11;
  float prop_inner_localphi_deg_GE11;
  float prop_inner_chi2_GE11;
  float prop_inner_ndof_GE11;

  //Prop from CSC
  float prop_CSC_GP_GE11[3];
  float prop_CSC_LP_GE11[3];
  float prop_CSC_GP_startingPoint[3];
  float prop_CSC_y_roll_GE11;
  float prop_CSC_localphi_rad_GE11;
  float prop_CSC_localphi_deg_GE11;
  float prop_CSC_chi2_GE11;
  float prop_CSC_ndof_GE11;

  //Prop Inner ME11
  float prop_innerSeg_GP_GE11[3];
  float prop_innerSeg_LP_GE11[3];
  float prop_innerSeg_GP_startingPoint[3];
  float prop_innerSeg_y_roll_GE11;
  float prop_innerSeg_localphi_rad_GE11;
  float prop_innerSeg_localphi_deg_GE11;

  //Prop Outer ME11
  float prop_outerSeg_GP_GE11[3];
  float prop_outerSeg_LP_GE11[3];
  float prop_outerSeg_GP_startingPoint[3];
  float prop_outerSeg_localx_GE11;
  float prop_outerSeg_localy_GE11;
  float prop_outerSeg_y_roll_GE11;
  float prop_outerSeg_localphi_rad_GE11;
  float prop_outerSeg_localphi_deg_GE11;

  //Prop Seg Only
  float prop_Seg_GP_GE11[3];
  float prop_Seg_LP_GE11[3];
  float prop_Seg_GP_startingPoint[3];
  float prop_Seg_localx_GE11;
  float prop_Seg_localy_GE11;
  float prop_Seg_y_roll_GE11;
  float prop_Seg_localphi_rad_GE11;
  float prop_Seg_localphi_deg_GE11;

  bool has_rechit_CSC_GE11;
  int rechit_location_CSC[5];
  float rechit_CSC_GP[3];
  float rechit_CSC_LP[3];
  int rechit_first_strip_CSC;
  int rechit_CLS_CSC;
  int rechit_BunchX_CSC;
  float rechit_y_roll_CSC_GE11;
  float rechit_localphi_rad_CSC_GE11;
  float rechit_localphi_deg_CSC_GE11;

  bool has_rechit_inner_GE11;
  int rechit_location_inner[5];
  float rechit_inner_GP[3];
  float rechit_inner_LP[3];
  int rechit_first_strip_inner;
  int rechit_CLS_inner;
  int rechit_BunchX_inner;
  float rechit_y_roll_inner_GE11;
  float rechit_localphi_rad_inner_GE11;
  float rechit_localphi_deg_inner_GE11;

  bool has_rechit_Seg_GE11;
  int rechit_location_Seg[5];
  float rechit_Seg_GP[3];
  float rechit_Seg_LP[3];
  int rechit_first_strip_Seg;
  int rechit_CLS_Seg;
  int rechit_BunchX_Seg;
  float rechit_y_roll_Seg_GE11;
  float rechit_localphi_rad_Seg_GE11;
  float rechit_localphi_deg_Seg_GE11;

  float RdPhi_inner_GE11;
  float RdPhi_inner_Corrected; // Short/Long chambers are mounted in opposite directions - local coordinates are flipped
  float RdPhi_CSC_GE11;
  float RdPhi_CSC_Corrected;
  float RdPhi_Seg_GE11;
  float RdPhi_Seg_Corrected;
  int det_id;


  bool has_fidcut_inner_GE11;
  bool has_fidcut_CSC_GE11;
  int isGEMmuon;

  int which_track_CSC_GE11;
  int which_track_inner_GE11;
  int hasME11;
  int hasME11RecHit;  
  int hasME11A;
  int hasME11ARecHit;  

  unsigned long long  evtNum;
  unsigned long long  lumiBlock;
  int muonIdx;

  int nCSCSeg;
  int nDTSeg;
  int CSCSeg_region;
  int num_props;

  float simDy;
  int nRecHitsTot;
  int nRecHits5;
  int nRecHits2;

  int closest;

  //ME11 track based prop

  float RdPhi_outerSeg_GE11;
  float RdPhi_innerSeg_GE11;

  int closest_ME11;
  float ME11_startingPoint[3];

  bool has_prop_inner;
  bool has_prop_CSC;
  bool has_prop_innerSeg;
  bool has_prop_outerSeg;
  bool has_prop_Seg;

  //Sim info for MC
  float sim_GP[3];
  float sim_LP[3];
  float sim_localy_roll;
  int nSim;

};

void MuonData::init()
{
  muon_charge = 9999;
  muon_pt = 9999;
  muon_eta = 9999;
  muon_momentum = 9999;

  for (int i=0; i<4; ++i){
    prop_location[i] = 99999;
  }

  //Prop Tracker
  for (int i=0; i<2; ++i){
    prop_inner_GP_GE11[i] = 99999;
    prop_inner_LP_GE11[i] = 99999;
    prop_inner_GP_startingPoint[i] = 99999;
  }
  prop_inner_y_roll_GE11 = 99999;
  prop_inner_localphi_rad_GE11 = 99999;
  prop_inner_localphi_deg_GE11 = 99999;
  prop_inner_chi2_GE11 = 99999;
  prop_inner_ndof_GE11 = 99999;

  //Prop CSC
  for (int i=0; i<2; ++i){
    prop_CSC_GP_GE11[i] = 99999;
    prop_CSC_LP_GE11[i] = 99999;
    prop_CSC_GP_startingPoint[i] = 99999;
  }
  prop_CSC_y_roll_GE11 = 99999;
  prop_CSC_localphi_rad_GE11 = 99999;
  prop_CSC_localphi_deg_GE11 = 99999;
  prop_CSC_chi2_GE11 = 99999;
  prop_CSC_ndof_GE11 = 99999;

  //Prop Outer ME11
  for (int i=0; i<2; ++i){
    prop_outerSeg_GP_GE11[i] = 9999999;
    prop_outerSeg_LP_GE11[i] = 9999999;
    prop_outerSeg_GP_startingPoint[i] = 9999999;
  }
  prop_outerSeg_y_roll_GE11 = 9999999;
  prop_outerSeg_localphi_rad_GE11 = 9999999;
  prop_outerSeg_localphi_deg_GE11 = 9999999;

  //Prop Inner ME11
  for (int i=0; i<2; ++i){
    prop_innerSeg_GP_GE11[i] = 9999999;
    prop_innerSeg_LP_GE11[i] = 9999999;
    prop_innerSeg_GP_startingPoint[i] = 9999999;
  }
  prop_innerSeg_y_roll_GE11 = 9999999;
  prop_innerSeg_localphi_rad_GE11 = 9999999;
  prop_innerSeg_localphi_deg_GE11 = 9999999;

  //Prop Seg Only
  for (int i=0; i<2; ++i){
    prop_Seg_GP_GE11[i] = 9999999;
    prop_Seg_LP_GE11[i] = 9999999;
    prop_Seg_GP_startingPoint[i] = 9999999;
  }
  prop_Seg_y_roll_GE11 = 9999999;
  prop_Seg_localphi_rad_GE11 = 9999999;
  prop_Seg_localphi_deg_GE11 = 9999999;


  has_rechit_CSC_GE11 = false;
  for (int i=0; i<4; ++i){
    rechit_location_CSC[i] = 999999;
  }
  for (int i=0; i<2; ++i){
    rechit_CSC_GP[i] = 999999;
    rechit_CSC_LP[i] = 999999;
  }
  rechit_first_strip_CSC = 999999;
  rechit_CLS_CSC = 999999;
  rechit_BunchX_CSC = 999999;
  rechit_y_roll_CSC_GE11 = 999999;
  rechit_localphi_rad_CSC_GE11 = 999999;
  rechit_localphi_deg_CSC_GE11 = 999999;

  has_rechit_inner_GE11 = false;
  for (int i=0; i<4; ++i){
    rechit_location_inner[i] = 999999;
  }
  for (int i=0; i<2; ++i){
    rechit_inner_GP[i] = 999999;
    rechit_inner_LP[i] = 999999;
  }
  rechit_first_strip_inner = 999999;
  rechit_CLS_inner = 999999;
  rechit_BunchX_inner = 999999;
  rechit_y_roll_inner_GE11 = 999999;
  rechit_localphi_rad_inner_GE11 = 999999;
  rechit_localphi_deg_inner_GE11 = 999999;

  has_rechit_Seg_GE11 = false;
  for (int i=0; i<4; ++i){
    rechit_location_Seg[i] = 999999;
  }
  for (int i=0; i<2; ++i){
    rechit_Seg_GP[i] = 999999;
    rechit_Seg_LP[i] = 999999;
  }
  rechit_first_strip_Seg = 999999;
  rechit_CLS_Seg = 999999;
  rechit_BunchX_Seg = 999999;
  rechit_y_roll_Seg_GE11 = 999999;
  rechit_localphi_rad_Seg_GE11 = 999999;
  rechit_localphi_deg_Seg_GE11 = 999999;

  RdPhi_inner_GE11 = 999999;
  RdPhi_inner_Corrected = 999999;
  RdPhi_CSC_GE11 = 999999;
  RdPhi_CSC_Corrected = 999999;
  RdPhi_Seg_GE11 = 999999;
  RdPhi_Seg_Corrected = 999999;
  det_id = 999999;

  has_fidcut_inner_GE11 = false;
  has_fidcut_CSC_GE11 = false;
  isGEMmuon = 0;

  which_track_CSC_GE11 = 999;
  which_track_inner_GE11 = 999;
  hasME11 = 0;
  hasME11RecHit = 0; 
  hasME11A = 0;
  hasME11ARecHit = 0; 

  evtNum = 0;
  lumiBlock = 0;
  muonIdx = 0;

  nCSCSeg = 0;
  nDTSeg = 0;
  CSCSeg_region = 0;
  num_props = 0;

  simDy = 999;
  nRecHits5 = 0;
  nRecHits2 = 0;
  nRecHitsTot = 0;

  closest = 999;

  RdPhi_outerSeg_GE11 = 9999999;
  RdPhi_innerSeg_GE11 = 9999999;

  closest_ME11 = 999;
  for (int i=0; i<2; ++i){
    ME11_startingPoint[i] = 9999999;
  }

  has_prop_inner = false;
  has_prop_CSC = false;
  has_prop_innerSeg = false;
  has_prop_outerSeg = false;
  has_prop_Seg = false;

  //Sim info for MC
  for (int i=0; i<2; ++i){
    sim_GP[i] = 99999999;
    sim_LP[i] = 99999999;
  }
  sim_localy_roll = 99999999;
  nSim = 99999999;
}

TTree* MuonData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("MuonData", "MuonData");

  t->Branch("muon_charge", &muon_charge);
  t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta);
  t->Branch("muon_momentum", &muon_momentum);

//Shared prop info
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, cha, lay, roll)/I");


//Propagated Inner
  t->Branch("prop_inner_GP_GE11", &prop_inner_GP_GE11, "prop_inner_GP_GE11[3] (x,y,z)/F");
  t->Branch("prop_inner_LP_GE11", &prop_inner_LP_GE11, "prop_inner_LP_GE11[3] (x,y,z)/F");
  t->Branch("prop_inner_GP_startingPoint", &prop_inner_GP_startingPoint, "prop_inner_GP_startingPoint[3] (x,y,z)/F");
  t->Branch("prop_inner_y_roll_GE11", &prop_inner_y_roll_GE11);
  t->Branch("prop_inner_localphi_rad_GE11", &prop_inner_localphi_rad_GE11);
  t->Branch("prop_inner_localphi_deg_GE11", &prop_inner_localphi_deg_GE11);
  t->Branch("prop_inner_chi2_GE11", &prop_inner_chi2_GE11);
  t->Branch("prop_inner_ndof_GE11", &prop_inner_ndof_GE11);

//Propagated CSC
  t->Branch("prop_CSC_GP_GE11", &prop_CSC_GP_GE11, "prop_CSC_GP_GE11[3] (x,y,z)/F");
  t->Branch("prop_CSC_LP_GE11", &prop_CSC_LP_GE11, "prop_CSC_LP_GE11[3] (x,y,z)/F");
  t->Branch("prop_CSC_GP_startingPoint", &prop_CSC_GP_startingPoint, "prop_CSC_GP_startingPoint[3] (x,y,z)/F");
  t->Branch("prop_CSC_y_roll_GE11", &prop_CSC_y_roll_GE11);
  t->Branch("prop_CSC_localphi_rad_GE11", &prop_CSC_localphi_rad_GE11);
  t->Branch("prop_CSC_localphi_deg_GE11", &prop_CSC_localphi_deg_GE11);
  t->Branch("prop_CSC_chi2_GE11", &prop_CSC_chi2_GE11);
  t->Branch("prop_CSC_ndof_GE11", &prop_CSC_ndof_GE11);

//Propagated Inner ME11
  t->Branch("prop_innerSeg_GP_GE11", &prop_innerSeg_GP_GE11, "prop_innerSeg_GP_GE11[3] (x,y,z)/F");
  t->Branch("prop_innerSeg_LP_GE11", &prop_innerSeg_LP_GE11, "prop_innerSeg_LP_GE11[3] (x,y,z)/F");
  t->Branch("prop_innerSeg_GP_startingPoint", &prop_innerSeg_GP_startingPoint, "prop_innerSeg_GP_startingPoint[3] (x,y,z)/F");
  t->Branch("prop_innerSeg_y_roll_GE11", &prop_innerSeg_y_roll_GE11);
  t->Branch("prop_innerSeg_localphi_rad_GE11", &prop_innerSeg_localphi_rad_GE11);
  t->Branch("prop_innerSeg_localphi_deg_GE11", &prop_innerSeg_localphi_deg_GE11);

//Propagated Outer ME11
  t->Branch("prop_outerSeg_GP_GE11", &prop_outerSeg_GP_GE11, "prop_outerSeg_GP_GE11[3] (x,y,z)/F");
  t->Branch("prop_outerSeg_LP_GE11", &prop_outerSeg_LP_GE11, "prop_outerSeg_LP_GE11[3] (x,y,z)/F");
  t->Branch("prop_outerSeg_GP_startingPoint", &prop_outerSeg_GP_startingPoint, "prop_outerSeg_GP_startingPoint[3] (x,y,z)/F");
  t->Branch("prop_outerSeg_y_roll_GE11", &prop_outerSeg_y_roll_GE11);
  t->Branch("prop_outerSeg_localphi_rad_GE11", &prop_outerSeg_localphi_rad_GE11);
  t->Branch("prop_outerSeg_localphi_deg_GE11", &prop_outerSeg_localphi_deg_GE11);

//Propagated Seg Only
  t->Branch("prop_Seg_GP_GE11", &prop_Seg_GP_GE11, "prop_Seg_GP_GE11[3] (x,y,z)/F");
  t->Branch("prop_Seg_LP_GE11", &prop_Seg_LP_GE11, "prop_Seg_LP_GE11[3] (x,y,z)/F");
  t->Branch("prop_Seg_GP_startingPoint", &prop_Seg_GP_startingPoint, "prop_Seg_GP_startingPoint[3] (x,y,z)/F");
  t->Branch("prop_Seg_y_roll_GE11", &prop_Seg_y_roll_GE11);
  t->Branch("prop_Seg_localphi_rad_GE11", &prop_Seg_localphi_rad_GE11);
  t->Branch("prop_Seg_localphi_deg_GE11", &prop_Seg_localphi_deg_GE11);

//Reconstructed
  t->Branch("has_rechit_CSC_GE11", &has_rechit_CSC_GE11);
  t->Branch("rechit_location_CSC", &rechit_location_CSC, "rechit_location_CSC[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("rechit_CSC_GP", &rechit_CSC_GP, "rechit_CSC_GP[3] (x,y,z)/F");
  t->Branch("rechit_CSC_LP", &rechit_CSC_LP, "rechit_CSC_LP[3] (x,y,z)/F");
  t->Branch("rechit_first_strip_CSC", &rechit_first_strip_CSC);
  t->Branch("rechit_CLS_CSC", &rechit_CLS_CSC);
  t->Branch("rechit_BunchX_CSC", &rechit_BunchX_CSC);
  t->Branch("rechit_y_roll_CSC_GE11", &rechit_y_roll_CSC_GE11);
  t->Branch("rechit_localphi_rad_CSC_GE11", &rechit_localphi_rad_CSC_GE11);
  t->Branch("rechit_localphi_deg_CSC_GE11", &rechit_localphi_deg_CSC_GE11);

  t->Branch("has_rechit_inner_GE11", &has_rechit_inner_GE11);
  t->Branch("rechit_location_inner", &rechit_location_inner, "rechit_location_inner[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("rechit_inner_GP", &rechit_inner_GP, "rechit_inner_GP[3] (x,y,z)/F");
  t->Branch("rechit_inner_LP", &rechit_inner_LP, "rechit_inner_LP[3] (x,y,z)/F");
  t->Branch("rechit_first_strip_inner", &rechit_first_strip_inner);
  t->Branch("rechit_CLS_inner", &rechit_CLS_inner);
  t->Branch("rechit_BunchX_inner", &rechit_BunchX_inner);
  t->Branch("rechit_y_roll_inner_GE11", &rechit_y_roll_inner_GE11);
  t->Branch("rechit_localphi_rad_inner_GE11", &rechit_localphi_rad_inner_GE11);
  t->Branch("rechit_localphi_deg_inner_GE11", &rechit_localphi_deg_inner_GE11);

  t->Branch("has_rechit_Seg_GE11", &has_rechit_Seg_GE11);
  t->Branch("rechit_location_Seg", &rechit_location_Seg, "rechit_location_Seg[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("rechit_Seg_GP", &rechit_Seg_GP, "rechit_Seg_GP[3] (x,y,z)/F");
  t->Branch("rechit_Seg_LP", &rechit_Seg_LP, "rechit_Seg_LP[3] (x,y,z)/F");
  t->Branch("rechit_first_strip_Seg", &rechit_first_strip_Seg);
  t->Branch("rechit_CLS_Seg", &rechit_CLS_Seg);
  t->Branch("rechit_BunchX_Seg", &rechit_BunchX_Seg);
  t->Branch("rechit_y_roll_Seg_GE11", &rechit_y_roll_Seg_GE11);
  t->Branch("rechit_localphi_rad_Seg_GE11", &rechit_localphi_rad_Seg_GE11);
  t->Branch("rechit_localphi_deg_Seg_GE11", &rechit_localphi_deg_Seg_GE11);

//Residual
  t->Branch("RdPhi_inner_GE11", &RdPhi_inner_GE11);
  t->Branch("RdPhi_inner_Corrected", &RdPhi_inner_Corrected);
  t->Branch("RdPhi_CSC_GE11", &RdPhi_CSC_GE11);
  t->Branch("RdPhi_CSC_Corrected", &RdPhi_CSC_Corrected);
  t->Branch("RdPhi_Seg_GE11", &RdPhi_Seg_GE11);
  t->Branch("RdPhi_Seg_Corrected", &RdPhi_Seg_Corrected);
  t->Branch("det_id", &det_id);
//Cut
  t->Branch("has_fidcut_inner_GE11", &has_fidcut_inner_GE11);
  t->Branch("has_fidcut_CSC_GE11", &has_fidcut_CSC_GE11);
  t->Branch("isGEMmuon", &isGEMmuon);

  t->Branch("which_track_CSC_GE11", &which_track_CSC_GE11);
  t->Branch("which_track_inner_GE11", &which_track_inner_GE11);
 
  t->Branch("hasME11", &hasME11);
  t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A);
  t->Branch("hasME11ARecHit", &hasME11ARecHit);

  t->Branch("evtNum", &evtNum);
  t->Branch("lumiBlock", &lumiBlock);
  t->Branch("muonIdx", &muonIdx);

  t->Branch("nCSCSeg", &nCSCSeg); 
  t->Branch("nDTSeg", &nDTSeg); 
  t->Branch("CSCSeg_region", &CSCSeg_region);
  t->Branch("num_props", &num_props);

  t->Branch("simDy", &simDy);
  t->Branch("nRecHits5", &nRecHits5);
  t->Branch("nRecHits2", &nRecHits2);
  t->Branch("nRecHitsTot", &nRecHitsTot);

  t->Branch("closest", &closest);

  //ME11 trackbased prop

  t->Branch("RdPhi_outerSeg_GE11", &RdPhi_outerSeg_GE11);
  t->Branch("RdPhi_innerSeg_GE11", &RdPhi_innerSeg_GE11);

  t->Branch("closest_ME11", &closest_ME11);
  t->Branch("ME11_startingPoint", &ME11_startingPoint, "ME11_startingPoint[3] (x,y,z)/F");

  t->Branch("has_prop_inner", &has_prop_inner);
  t->Branch("has_prop_CSC", &has_prop_CSC);
  t->Branch("has_prop_innerSeg", &has_prop_innerSeg);
  t->Branch("has_prop_outerSeg", &has_prop_outerSeg);
  t->Branch("has_prop_Seg", &has_prop_Seg);

  //Sim info for MC
  t->Branch("sim_GP", &sim_GP, "sim_GP[3] (x,y,z)/F");
  t->Branch("sim_LP", &sim_LP, "sim_LP[3] (x,y,z)/F");
  t->Branch("sim_localy_roll", &sim_localy_roll);
  t->Branch("nSim", &nSim);

  return t;
}

float RdPhi(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, edm::ESHandle<GEMGeometry> GEMGeometry_, const GEMEtaPartition* ch);

void propagate_segment(const CSCSegment* ME11_segment, const reco::Track* innerTrack, const GEMEtaPartition* ch, const reco::Muon* mu, MuonServiceProxy* theService_, const GeomDet* segDet, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool debug);

void propagate_track_based_segment(const CSCSegment* ME11_segment, reco::TransientTrack track, const GEMEtaPartition* ch, const reco::Muon* mu, MuonServiceProxy* theService_, const GeomDet* segDet, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool debug);

void propagate_track(reco::TransientTrack track, const GEMEtaPartition* ch, const string inner_or_CSC, MuonServiceProxy* theService_, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool &which_track, bool debug);

bool fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch);

class analyser : public edm::EDAnalyzer {
public:
  explicit analyser(const edm::ParameterSet&);
  ~analyser(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;

  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  edm::ESHandle<GEMGeometry> GEMGeometry_;

  bool tracker_prop;
  bool CSC_prop;
  bool debug;

  TTree * tree_data_;
  MuonData data_;
  const CSCSegment *tmp_ME11_seg;
};

analyser::analyser(const edm::ParameterSet& iConfig)
{
  cout << "Begin analyser" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gemSimHits_ = consumes<vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("gemSimHits"));

  tracker_prop = iConfig.getParameter<bool>("tracker_prop");
  CSC_prop = iConfig.getParameter<bool>("CSC_prop");
  debug = iConfig.getParameter<bool>("debug");
  cout << "tracker_prop " << tracker_prop << " CSC_prop " << CSC_prop << " debug " << debug << std::endl;

  tree_data_ = data_.book(tree_data_);
}


void
analyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);
 
  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  bool isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);
  edm::Handle<vector<PSimHit> > gemSimHits;
  if (isMC) {
    iEvent.getByToken(gemSimHits_, gemSimHits); 
  }
  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;

  //isMC = false;
  int num_props_ME11 = 0;
  int num_props_noME11 = 0;
  int num_props = 0;
  int muons_with_cscSeg = 0;
  cout << "new evt numb is " << iEvent.eventAuxiliary().event() << " and new lumiblock is " << iEvent.eventAuxiliary().luminosityBlock() << endl;
  for (size_t i = 0; i < muons->size(); ++i){
    data_.init();
    cout << "new muon" << endl;
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    cout << "is standalone" << endl;
    data_.isGEMmuon = mu->isGEMMuon();
    data_.nCSCSeg = mu->numberOfSegments(1,2) + mu->numberOfSegments(2,2) + mu->numberOfSegments(3,2) + mu->numberOfSegments(4,2);
    data_.nDTSeg = mu->numberOfSegments(1,1) + mu->numberOfSegments(2,1) + mu->numberOfSegments(3,1) + mu->numberOfSegments(4,1);
    if (data_.nCSCSeg > 0){
      muons_with_cscSeg++;
    }
    auto matches = mu->matches();
    int tmp_station = 99;
    int tmp_ring = 99;
    for ( auto MCM : matches){
      if (MCM.detector() != 2) continue;
      for( auto MSM : MCM.segmentMatches){
        if (debug){std::cout << "Looping over segments" << std::endl;}
        auto cscSegRef = MSM.cscSegmentRef;
        auto cscDetID = cscSegRef->cscDetId();
        data_.CSCSeg_region = cscDetID.endcap();
        if (cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
          if (debug && data_.hasME11 == 1) {std::cout << "Already has an ME11 seg" << std::endl;}
          data_.hasME11 = 1;
          if (debug){std::cout << "has ME11 segment" << std::endl;}
          tmp_ME11_seg = cscSegRef.get();
          for ( auto rh : cscSegRef->specificRecHits()){
            if (rh.cscDetId().ring() == 1) data_.hasME11RecHit = 1;
            if (rh.cscDetId().ring() == 4) data_.hasME11ARecHit = 1;
          }
        }

        //Closest segment finder! Start from station 4 and work inwards
	if (cscDetID.station() < tmp_station){
	  tmp_station = cscDetID.station();
	  tmp_ring = cscDetID.ring();
	}

	if (cscDetID.station() == tmp_station){
	  if (cscDetID.ring() < tmp_ring){
	    tmp_ring = cscDetID.ring();
	  }
	}

	data_.closest = std::stoi(std::to_string(tmp_station)+std::to_string(tmp_ring));
	if (debug){std::cout << "Closest now says " << data_.closest << std::endl;}


        if (cscDetID.station() == 1 and cscDetID.ring() == 4) data_.hasME11A = 1;
      }
    }
    
    data_.evtNum = iEvent.eventAuxiliary().event();
    data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock();
    data_.muonIdx = data_.evtNum*100 + i;

    //Get the tracks
    if (debug){std::cout << "Getting tracks" << std::endl;}
    const reco::Track* globalTrack = 0;
    const reco::Track* innerTrack = 0;
    const reco::Track* outerTrack = 0;
    reco::TransientTrack ttTrack_global;
    reco::TransientTrack ttTrack_tracker;
    reco::TransientTrack ttTrack_CSC;
    if ( mu->globalTrack().isNonnull() ){ globalTrack = mu->globalTrack().get(); ttTrack_global = ttrackBuilder_->build(globalTrack);}
    if ( mu->track().isNonnull() ){ innerTrack = mu->track().get(); ttTrack_tracker = ttrackBuilder_->build(innerTrack);}
    if ( mu->outerTrack().isNonnull() ){ outerTrack = mu->outerTrack().get(); ttTrack_CSC = ttrackBuilder_->build(outerTrack);}

    //The plan
    //I want to loop over all rechits in the track to count the number of SEGMENTS and save the ME1/1 segment
    
    //for (auto CSChit_iter = outerTrack->recHitsBegin(); CSChit_iter != outerTrack->recHitsEnd(); CSChit_iter++){
    int my_CSCseg_counter = 0;
    int my_DTseg_counter = 0;
    if (debug){std::cout << "Total number of hits to loop over = " << outerTrack->recHitsSize() << std::endl;}
    for (size_t RecHit_iter = 0; RecHit_iter != outerTrack->recHitsSize(); RecHit_iter++){
      //std::cout << "Looping the track hits?" << std::endl;
      //TrackingRecHitRef CSChit = outerTrack->recHit(CSChit_iter).Get();
      const TrackingRecHit* RecHit = (outerTrack->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId();
      uint16_t RecHitDetId = RecHitId.det();
      if (RecHitDetId == DetId::Muon){
        uint16_t RecHitSubDet = RecHitId.subdetId();
        if (RecHitSubDet == (uint16_t)MuonSubdetId::CSC){
          if (debug){std::cout << "CSC hit found: Dimensions = " << RecHit->dimension() << " and DetID " << CSCDetId(RecHitId) << std::endl;}
          if (RecHit->dimension() == 4){
            if(debug){std::cout << "CSC Segment! DetID " << CSCDetId(RecHitId) << std::endl;}
            my_CSCseg_counter++;
          }
        }
        if (RecHitSubDet == (uint16_t)MuonSubdetId::DT){
          if (debug){std::cout << "DT hit found: Dimensions = " << RecHit->dimension() << " and DetID " << DTLayerId(RecHitId) << std::endl;}
          if (RecHit->dimension() > 1){
            if (debug){std::cout << "DT Segment! DetID " << DTLayerId(RecHitId) << std::endl;}
            my_DTseg_counter++;
          }
        }
      }
    }
    if (debug){
      std::cout << "Hyunyong's CSC counter = " << data_.nCSCSeg << std::endl;
      std::cout << "My CSC counter         = " << my_CSCseg_counter << std::endl;
      std::cout << "Hyunyong's DT counter  = " << data_.nDTSeg << std::endl;
      std::cout << "My DT counter          = " << my_DTseg_counter << std::endl;
    }
    data_.nCSCSeg = my_CSCseg_counter;
    data_.nDTSeg = my_DTseg_counter;

    //Start the propagations
    float count = 0;
    if (debug){std::cout << "Starting chamber loop" << std::endl;}
    for (const auto& ch : GEMGeometry_->etaPartitions()) {
      if (ch->id().station() != 1) continue; //Only takes GE1/1
      const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
      const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp(); //y distance to the current eta part
      const float prop_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2))).y();

      //Props from ME11 segment
      data_.has_prop_innerSeg = 0;
      data_.has_prop_outerSeg = 0;
      data_.has_prop_Seg = 0;
      GlobalPoint pos_global_ch_outerSeg;
      GlobalPoint pos_global_start_outerSeg;
      GlobalPoint pos_global_ch_innerSeg;
      GlobalPoint pos_global_start_innerSeg;
      GlobalPoint pos_global_ch_Seg;
      GlobalPoint pos_global_start_Seg;

      if (data_.hasME11 == 1 and ch->id().station() == 1 and ch->id().ring() == 1){
        if (debug){std::cout << "Doing segment propagations" << std::endl;}
        DetId segDetId = tmp_ME11_seg->geographicalId();
        const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
        data_.ME11_startingPoint[0] = segDet->toGlobal(tmp_ME11_seg->localPosition()).x(); data_.ME11_startingPoint[1] = segDet->toGlobal(tmp_ME11_seg->localPosition()).y(); data_.ME11_startingPoint[2] = segDet->toGlobal(tmp_ME11_seg->localPosition()).z();
        if (debug){std::cout << "Got segment starting position" << std::endl;}

        //Segment only propagation
        bool has_prop_Seg = 0;
        propagate_segment(tmp_ME11_seg, innerTrack, ch, mu, theService_, segDet, pos_global_ch_Seg, pos_global_start_Seg, has_prop_Seg, debug);
        if (has_prop_Seg){
          data_.has_prop_Seg = 1;
          LocalPoint pos_local_ch_Seg = ch->toLocal(pos_global_ch_Seg);
          if (debug){std::cout << "Seg success, local = " << pos_local_ch_Seg << std::endl;}
          LocalPoint local_to_center_Seg(pos_local_ch_Seg.x(), prop_y_to_center + pos_local_ch_Seg.y(), 0);
          const float prop_Seg_ch_localphi_rad = (3.14159265/2.) - local_to_center_Seg.phi();
          const float prop_Seg_ch_localphi_deg = prop_Seg_ch_localphi_rad*180/3.14159265;
          data_.prop_Seg_GP_GE11[0] = pos_global_ch_Seg.x(); data_.prop_Seg_GP_GE11[1] = pos_global_ch_Seg.y(); data_.prop_Seg_GP_GE11[2] = pos_global_ch_Seg.z();
          data_.prop_Seg_LP_GE11[0] = pos_local_ch_Seg.x(); data_.prop_Seg_LP_GE11[1] = prop_y_to_chamber + pos_local_ch_Seg.y(); data_.prop_Seg_LP_GE11[2] = pos_local_ch_Seg.z();
          data_.prop_Seg_GP_startingPoint[0] = pos_global_start_Seg.x(); data_.prop_Seg_GP_startingPoint[1] = pos_global_start_Seg.y(); data_.prop_Seg_GP_startingPoint[2] = pos_global_start_Seg.z();
          data_.prop_Seg_y_roll_GE11 = pos_local_ch_Seg.y();
          data_.prop_Seg_localphi_rad_GE11 = prop_Seg_ch_localphi_rad;
          data_.prop_Seg_localphi_deg_GE11 = prop_Seg_ch_localphi_deg;
        }


        //Track propagation starting at ME11 segment location
        if (tracker_prop and ttTrack_tracker.isValid()){
          bool has_prop_innerSeg = 0;
          propagate_track_based_segment(tmp_ME11_seg, ttTrack_tracker, ch, mu, theService_, segDet, pos_global_ch_innerSeg, pos_global_start_innerSeg, has_prop_innerSeg, debug);
          if (has_prop_innerSeg){
            data_.has_prop_innerSeg = 1;
            LocalPoint pos_local_ch_innerSeg = ch->toLocal(pos_global_ch_innerSeg);
            if (debug){std::cout << "innerSeg success, local = " << pos_local_ch_innerSeg << std::endl;}

            LocalPoint local_to_center_innerSeg(pos_local_ch_innerSeg.x(), prop_y_to_center + pos_local_ch_innerSeg.y(), 0);
            const float prop_innerSeg_ch_localphi_rad = (3.14159265/2.) - local_to_center_innerSeg.phi();
            const float prop_innerSeg_ch_localphi_deg = prop_innerSeg_ch_localphi_rad*180/3.14159265;

            data_.prop_innerSeg_GP_GE11[0] = pos_global_ch_innerSeg.x(); data_.prop_innerSeg_GP_GE11[1] = pos_global_ch_innerSeg.y(); data_.prop_innerSeg_GP_GE11[2] = pos_global_ch_innerSeg.z();
            data_.prop_innerSeg_LP_GE11[0] = pos_local_ch_innerSeg.x(); data_.prop_innerSeg_LP_GE11[1] = prop_y_to_chamber + pos_local_ch_innerSeg.y(); data_.prop_innerSeg_LP_GE11[2] = pos_local_ch_innerSeg.z();
            data_.prop_innerSeg_GP_startingPoint[0] = pos_global_start_innerSeg.x(); data_.prop_innerSeg_GP_startingPoint[1] = pos_global_start_innerSeg.y(); data_.prop_innerSeg_GP_startingPoint[2] = pos_global_start_innerSeg.z();
            data_.prop_innerSeg_y_roll_GE11 = pos_local_ch_innerSeg.y();
            data_.prop_innerSeg_localphi_rad_GE11 = prop_innerSeg_ch_localphi_rad;
            data_.prop_innerSeg_localphi_deg_GE11 = prop_innerSeg_ch_localphi_deg;
          }
        }

        //CSC propagation starting at ME11 segment location
        if (CSC_prop and ttTrack_CSC.isValid()){
          bool has_prop_outerSeg = 0;
          propagate_track_based_segment(tmp_ME11_seg, ttTrack_CSC, ch, mu, theService_, segDet, pos_global_ch_outerSeg, pos_global_start_outerSeg, has_prop_outerSeg, debug);
          if (has_prop_outerSeg){
            data_.has_prop_outerSeg = 1;
            LocalPoint pos_local_ch_outerSeg = ch->toLocal(pos_global_ch_outerSeg);
            if (debug){std::cout << "outerSeg success, local = " << pos_local_ch_outerSeg << std::endl;}

            LocalPoint local_to_center_outerSeg(pos_local_ch_outerSeg.x(), prop_y_to_center + pos_local_ch_outerSeg.y(), 0);
            const float prop_outerSeg_ch_localphi_rad = (3.14159265/2.) - local_to_center_outerSeg.phi();
            const float prop_outerSeg_ch_localphi_deg = prop_outerSeg_ch_localphi_rad*180/3.14159265;

            data_.prop_outerSeg_GP_GE11[0] = pos_global_ch_outerSeg.x(); data_.prop_outerSeg_GP_GE11[1] = pos_global_ch_outerSeg.y(); data_.prop_outerSeg_GP_GE11[2] = pos_global_ch_outerSeg.z();
            data_.prop_outerSeg_LP_GE11[0] = pos_local_ch_outerSeg.x(); data_.prop_outerSeg_LP_GE11[1] = prop_y_to_chamber + pos_local_ch_outerSeg.y(); data_.prop_outerSeg_LP_GE11[2] = pos_local_ch_outerSeg.z();
            data_.prop_outerSeg_GP_startingPoint[0] = pos_global_start_outerSeg.x(); data_.prop_outerSeg_GP_startingPoint[1] = pos_global_start_outerSeg.y(); data_.prop_outerSeg_GP_startingPoint[2] = pos_global_start_outerSeg.z();
            data_.prop_outerSeg_y_roll_GE11 = pos_local_ch_outerSeg.y();
            data_.prop_outerSeg_localphi_rad_GE11 = prop_outerSeg_ch_localphi_rad;
            data_.prop_outerSeg_localphi_deg_GE11 = prop_outerSeg_ch_localphi_deg;
          }
        }
      }


      //Props from tracks
      data_.has_prop_inner = 0;
      data_.has_prop_CSC = 0;
      GlobalPoint pos_global_ch_CSC;      // Outer props ending pos
      GlobalPoint pos_global_start_CSC;   // Outer props starting pos
      GlobalPoint pos_global_ch_inner;    // Inner props ending pos
      GlobalPoint pos_global_start_inner; // Inner props starting pos

 
      // Tracker propagated
      if (tracker_prop and ttTrack_tracker.isValid()){
        bool has_prop_inner = 0;
        bool which_track_inner = 0;
        const string inner_or_CSC = "inner";
        propagate_track(ttTrack_tracker, ch, inner_or_CSC, theService_, pos_global_ch_inner, pos_global_start_inner, has_prop_inner, which_track_inner, debug);
        if (has_prop_inner){
          data_.has_prop_inner = 1;
          data_.which_track_inner_GE11 = which_track_inner;
          LocalPoint pos_local_ch_inner = ch->toLocal(pos_global_ch_inner);
          if (debug){std::cout << "inner success, local = " << pos_local_ch_inner << std::endl;}
          LocalPoint local_to_center_inner(pos_local_ch_inner.x(), prop_y_to_center + pos_local_ch_inner.y(), 0);
          const float prop_inner_localphi_rad = (3.14159265/2.) - local_to_center_inner.phi();
          const float prop_inner_localphi_deg = prop_inner_localphi_rad*180/3.14169265;

          data_.prop_inner_GP_GE11[0] = pos_global_ch_inner.x(); data_.prop_inner_GP_GE11[1] = pos_global_ch_inner.y(); data_.prop_inner_GP_GE11[2] = pos_global_ch_inner.z();
          data_.prop_inner_LP_GE11[0] = pos_local_ch_inner.x(); data_.prop_inner_LP_GE11[1] = prop_y_to_chamber + pos_local_ch_inner.y(); data_.prop_inner_LP_GE11[2] = pos_local_ch_inner.z();
          data_.prop_inner_GP_startingPoint[0] = pos_global_start_inner.x(); data_.prop_inner_GP_startingPoint[1] = pos_global_start_inner.y(); data_.prop_inner_GP_startingPoint[2] = pos_global_start_inner.z();
          data_.prop_inner_y_roll_GE11 = pos_local_ch_inner.y();
          data_.prop_inner_localphi_rad_GE11 = prop_inner_localphi_rad;
          data_.prop_inner_localphi_deg_GE11 = prop_inner_localphi_deg;
          data_.prop_inner_chi2_GE11 = innerTrack->chi2();
          data_.prop_inner_ndof_GE11 = innerTrack->ndof();
          data_.has_fidcut_inner_GE11 = fidcutCheck(pos_local_ch_inner.y(), prop_inner_localphi_deg, ch);
        }
      }
      // CSC propagated
      if (CSC_prop and ttTrack_CSC.isValid()){
        bool has_prop_CSC = 0;
        bool which_track_CSC = 0;
        const string inner_or_CSC = "CSC";
        propagate_track(ttTrack_CSC, ch, inner_or_CSC, theService_, pos_global_ch_CSC, pos_global_start_CSC, has_prop_CSC, which_track_CSC, debug);
        if (has_prop_CSC){
          data_.has_prop_CSC = 1;
          data_.which_track_CSC_GE11 = which_track_CSC;
          LocalPoint pos_local_ch_CSC = ch->toLocal(pos_global_ch_CSC);
          if (debug){std::cout << "CSC success, local = " << pos_local_ch_CSC << std::endl;}
          LocalPoint local_to_center_CSC(pos_local_ch_CSC.x(), prop_y_to_center + pos_local_ch_CSC.y(), 0);
          const float prop_CSC_localphi_rad = (3.14159265/2.) - local_to_center_CSC.phi();
          const float prop_CSC_localphi_deg = prop_CSC_localphi_rad*180/3.14169265;

          data_.prop_CSC_GP_GE11[0] = pos_global_ch_CSC.x(); data_.prop_CSC_GP_GE11[1] = pos_global_ch_CSC.y(); data_.prop_CSC_GP_GE11[2] = pos_global_ch_CSC.z();
          data_.prop_CSC_LP_GE11[0] = pos_local_ch_CSC.x(); data_.prop_CSC_LP_GE11[1] = prop_y_to_chamber + pos_local_ch_CSC.y(); data_.prop_CSC_LP_GE11[2] = pos_local_ch_CSC.z();
          data_.prop_CSC_GP_startingPoint[0] = pos_global_start_CSC.x(); data_.prop_CSC_GP_startingPoint[1] = pos_global_start_CSC.y(); data_.prop_CSC_GP_startingPoint[2] = pos_global_start_CSC.z();
          data_.prop_CSC_y_roll_GE11 = pos_local_ch_CSC.y();
          data_.prop_CSC_localphi_rad_GE11 = prop_CSC_localphi_rad;
          data_.prop_CSC_localphi_deg_GE11 = prop_CSC_localphi_deg;
          data_.prop_CSC_chi2_GE11 = outerTrack->chi2();
          data_.prop_CSC_ndof_GE11 = outerTrack->ndof();
          data_.has_fidcut_CSC_GE11 = fidcutCheck(pos_local_ch_CSC.y(), prop_CSC_localphi_deg, ch);
        }
      }

      //Skip chambers with no props
      if (!(data_.has_prop_CSC or data_.has_prop_inner or data_.has_prop_outerSeg or data_.has_prop_innerSeg)) continue;

      if(debug){cout << "charge is " << mu->charge() << endl;}
      data_.num_props++;
      data_.muon_charge = mu->charge();
      data_.muon_pt = mu->pt();
      data_.muon_eta = mu->eta();
      data_.muon_momentum = mu->momentum().mag2();

      count++;

      data_.prop_location[0] = ch->id().region(); data_.prop_location[1] = ch->id().station(); data_.prop_location[2] = ch->id().chamber(); data_.prop_location[3] = ch->id().layer(); data_.prop_location[4] = ch->id().roll();
        
      data_.has_rechit_CSC_GE11 = false;
      data_.RdPhi_CSC_GE11 = 999999;
      data_.has_rechit_inner_GE11 = false;
      data_.RdPhi_inner_GE11 = 999999;
      data_.nRecHits5 = 0;
      data_.nRecHits2 = 0;
      data_.nRecHitsTot = 0;
      int rechit_counter = 0;
      int rechit_matches = 0;
      int tmpNRH5 = 0;
      int tmpNRH2 = 0;
      int tmpNRHT = 0;

      for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
        rechit_counter++;
        if ( (hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
          GEMDetId gemid((hit)->geographicalId());
          if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
            cout << "starting rechit" << endl;
            tmpNRHT++;
            const auto& etaPart = GEMGeometry_->etaPartition(gemid);
            float strip = etaPart->strip(hit->localPosition());
            float stripAngle = etaPart->specificTopology().stripAngle(strip);

            float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
            float rechit_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y();
            LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
            float rechit_localphi_rad = (3.14159265/2.) - local_to_center.phi();
            float rechit_localphi_deg = rechit_localphi_rad*180/3.14159265;


            if (ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - data_.prop_CSC_LP_GE11[0]) < 999.0){
              if (abs(RdPhi(stripAngle, hit, data_.prop_CSC_LP_GE11[0], data_.prop_CSC_y_roll_GE11, GEMGeometry_, ch)) < 5) tmpNRH5++;
              if (abs(RdPhi(stripAngle, hit, data_.prop_CSC_LP_GE11[0], data_.prop_CSC_y_roll_GE11, GEMGeometry_, ch)) < 2) tmpNRH2++;
              data_.det_id = gemid.region()*(gemid.station()*100 + gemid.chamber());

              //CSC matcher
              if (data_.has_prop_CSC and abs(data_.RdPhi_CSC_GE11) > abs(RdPhi(stripAngle, hit, data_.prop_CSC_LP_GE11[0], data_.prop_CSC_y_roll_GE11, GEMGeometry_, ch))){
                rechit_matches++;
                if (debug){std::cout << "Overwrite CSC rechit" << std::endl;}

                data_.has_rechit_CSC_GE11 = true;
                data_.rechit_location_CSC[0] = gemid.region(); data_.rechit_location_CSC[1] = gemid.station(); data_.rechit_location_CSC[2] = gemid.chamber(); data_.rechit_location_CSC[3] = gemid.layer(); data_.rechit_location_CSC[4] = gemid.roll();
                data_.rechit_CSC_GP[0] = etaPart->toGlobal((hit)->localPosition()).x(); data_.rechit_CSC_GP[1] = etaPart->toGlobal((hit)->localPosition()).y(); data_.rechit_CSC_GP[2] = etaPart->toGlobal((hit)->localPosition()).z();
                data_.rechit_CSC_LP[0] = (hit)->localPosition().x(); data_.rechit_CSC_LP[1] = rechit_y_to_chamber + (hit)->localPosition().y(); data_.rechit_CSC_LP[2] = (hit)->localPosition().z();

                data_.rechit_first_strip_CSC = (hit)->firstClusterStrip();
                data_.rechit_CLS_CSC = (hit)->clusterSize();
                data_.rechit_BunchX_CSC = (hit)->BunchX();
                data_.rechit_y_roll_CSC_GE11 = (hit)->localPosition().y();
                data_.rechit_localphi_rad_CSC_GE11 = rechit_localphi_rad;
                data_.rechit_localphi_deg_CSC_GE11 = rechit_localphi_deg;
                data_.RdPhi_CSC_GE11 = RdPhi(stripAngle, hit, data_.prop_CSC_LP_GE11[0], data_.prop_CSC_y_roll_GE11, GEMGeometry_, ch);
                data_.RdPhi_CSC_Corrected = data_.RdPhi_CSC_GE11;
                if ((gemid.region() == 1 && gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)){
                  data_.RdPhi_CSC_Corrected = -1.0*data_.RdPhi_CSC_Corrected;
                  if (debug){std::cout << "CORRECTING THE RDPHI" << std::endl;}
                }
              }
              //inner matcher
              if (data_.has_prop_inner and abs(data_.RdPhi_inner_GE11) > abs(RdPhi(stripAngle, hit, data_.prop_inner_LP_GE11[0], data_.prop_inner_y_roll_GE11, GEMGeometry_, ch))){
                rechit_matches++;
                if (debug){std::cout << "Overwrite CSC rechit" << std::endl;}

                data_.has_rechit_inner_GE11 = true;
                data_.rechit_location_inner[0] = gemid.region(); data_.rechit_location_inner[1] = gemid.station(); data_.rechit_location_inner[2] = gemid.chamber(); data_.rechit_location_inner[3] = gemid.layer(); data_.rechit_location_inner[4] = gemid.roll();
                data_.rechit_inner_GP[0] = etaPart->toGlobal((hit)->localPosition()).x(); data_.rechit_inner_GP[1] = etaPart->toGlobal((hit)->localPosition()).y(); data_.rechit_inner_GP[2] = etaPart->toGlobal((hit)->localPosition()).z();
                data_.rechit_inner_LP[0] = (hit)->localPosition().x(); data_.rechit_inner_LP[1] = rechit_y_to_chamber + (hit)->localPosition().y(); data_.rechit_inner_LP[2] = (hit)->localPosition().z();

                data_.rechit_first_strip_inner = (hit)->firstClusterStrip();
                data_.rechit_CLS_inner = (hit)->clusterSize();
                data_.rechit_BunchX_inner = (hit)->BunchX();
                data_.rechit_y_roll_inner_GE11 = (hit)->localPosition().y();
                data_.rechit_localphi_rad_inner_GE11 = rechit_localphi_rad;
                data_.rechit_localphi_deg_inner_GE11 = rechit_localphi_deg;
                data_.RdPhi_inner_GE11 = RdPhi(stripAngle, hit, data_.prop_inner_LP_GE11[0], data_.prop_inner_y_roll_GE11, GEMGeometry_, ch);
                data_.RdPhi_inner_Corrected = data_.RdPhi_CSC_GE11;
                if ((gemid.region() == 1 && gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)){
                  if (debug){std::cout << "CORRECTING THE RDPHI" << std::endl;}
                  data_.RdPhi_inner_Corrected = -1.0*data_.RdPhi_inner_Corrected;
                }
              }
              if(debug){cout << "Starting ME11 rechit match" << endl;}
              //ME11 seg matcher
              if (data_.has_prop_Seg and abs(data_.RdPhi_Seg_GE11) > abs(RdPhi(stripAngle, hit, data_.prop_Seg_LP_GE11[0], data_.prop_Seg_y_roll_GE11, GEMGeometry_, ch))){
                if (debug){std::cout << "Overwrite Seg rechit" << std::endl;}

                data_.has_rechit_Seg_GE11 = true;
                data_.rechit_location_Seg[0] = gemid.region(); data_.rechit_location_Seg[1] = gemid.station(); data_.rechit_location_Seg[2] = gemid.chamber(); data_.rechit_location_Seg[3] = gemid.layer(); data_.rechit_location_Seg[4] = gemid.roll();
                data_.rechit_Seg_GP[0] = etaPart->toGlobal((hit)->localPosition()).x(); data_.rechit_Seg_GP[1] = etaPart->toGlobal((hit)->localPosition()).y(); data_.rechit_Seg_GP[2] = etaPart->toGlobal((hit)->localPosition()).z();
                data_.rechit_Seg_LP[0] = (hit)->localPosition().x(); data_.rechit_Seg_LP[1] = rechit_y_to_chamber + (hit)->localPosition().y(); data_.rechit_Seg_LP[2] = (hit)->localPosition().z();

                data_.rechit_first_strip_Seg = (hit)->firstClusterStrip();
                data_.rechit_CLS_Seg = (hit)->clusterSize();
                data_.rechit_BunchX_Seg = (hit)->BunchX();
                data_.rechit_y_roll_Seg_GE11 = (hit)->localPosition().y();
                data_.rechit_localphi_rad_Seg_GE11 = rechit_localphi_rad;
                data_.rechit_localphi_deg_Seg_GE11 = rechit_localphi_deg;
                data_.RdPhi_Seg_GE11 = RdPhi(stripAngle, hit, data_.prop_Seg_LP_GE11[0], data_.prop_Seg_y_roll_GE11, GEMGeometry_, ch);
                data_.RdPhi_Seg_Corrected = data_.RdPhi_Seg_GE11;
                if ((gemid.region() == 1 && gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)){
                  data_.RdPhi_Seg_Corrected = -1.0*data_.RdPhi_Seg_Corrected;
                  if (debug){std::cout << "CORRECTING THE RDPHI" << std::endl;}
                }
              }
              if (data_.has_prop_outerSeg){
                if (abs(data_.RdPhi_outerSeg_GE11) > abs(RdPhi(stripAngle, hit, data_.prop_outerSeg_LP_GE11[0], data_.prop_outerSeg_y_roll_GE11, GEMGeometry_, ch))){
                  data_.RdPhi_outerSeg_GE11 = RdPhi(stripAngle, hit, data_.prop_outerSeg_LP_GE11[0], data_.prop_outerSeg_y_roll_GE11, GEMGeometry_, ch);
                }
              }
              if (data_.has_prop_innerSeg){
                if (abs(data_.RdPhi_innerSeg_GE11) > abs(RdPhi(stripAngle, hit, data_.prop_innerSeg_LP_GE11[0], data_.prop_innerSeg_y_roll_GE11, GEMGeometry_, ch))){
                  data_.RdPhi_innerSeg_GE11 = RdPhi(stripAngle, hit, data_.prop_innerSeg_LP_GE11[0], data_.prop_innerSeg_y_roll_GE11, GEMGeometry_, ch);
                }
              }
            }
          }
        } 
      }
      data_.nRecHits5 = tmpNRH5;
      data_.nRecHits2 = tmpNRH2;
      data_.nRecHitsTot = tmpNRHT;
      if (isMC and data_.has_prop_CSC) {
        if(debug){cout << "Starting sim info" << endl;}
        data_.nSim = 99999999;
        data_.simDy = 999.;
        float tmpDy = 999.;
        float tmpDr = 999.;
        int tmpSimCounter = 0;
        for (const auto& simHit:*gemSimHits.product()) {
          GEMDetId gemid((simHit).detUnitId());
          if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
            if(debug){cout << "Found a match simhit" << endl;}
            tmpSimCounter ++;
            const auto& etaPart = GEMGeometry_->etaPartition(gemid);
            //GlobalPoint pGlobal = pos_global_ch_CSC;
            float dy = pos_global_ch_CSC.y() - etaPart->toGlobal(simHit.localPosition()).y();
            float dx = pos_global_ch_CSC.x() - etaPart->toGlobal(simHit.localPosition()).x();
            if (dy < tmpDy) tmpDy = dy;
            if (pow(pow(dy, 2) + pow(dx, 2), 0.5) < tmpDr){
              data_.sim_GP[0] = etaPart->toGlobal(simHit.localPosition()).x();
              data_.sim_GP[1] = etaPart->toGlobal(simHit.localPosition()).y();
              data_.sim_GP[2] = etaPart->toGlobal(simHit.localPosition()).z();
              data_.sim_LP[0] = simHit.localPosition().x();
              data_.sim_LP[1] = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y() + simHit.localPosition().y();
              data_.sim_LP[2] = simHit.localPosition().z();
              data_.sim_localy_roll = simHit.localPosition().y();
              tmpDr = pow(pow(dy, 2) + pow(dx, 2), 0.5);
            }
          }
        }
        data_.simDy = tmpDy;
        data_.nSim = tmpSimCounter;
      }

      std::cout << "Sim point was " << data_.sim_LP[0] << ", " << data_.sim_LP[1] << std::endl;

      if (debug){std::cout << "Num of rechits = " << rechit_counter << std::endl;}
      if (debug){std::cout << "Num of matches = " << rechit_matches << std::endl;}
      cout << "Filling!" << endl;
      if (data_.hasME11 == 1){num_props_ME11++;}
      if (data_.hasME11 != 1){num_props_noME11 ++;}
      num_props++;
      tree_data_->Fill();
    }
    //cout << "Filling!" << endl;
    //tree_data_->Fill();
  }
  std::cout << "Muons with cscSegs = " << muons_with_cscSeg << std::endl;
  std::cout << "Muons with prop to GE11 = " << num_props << std::endl;
  std::cout << "Muons with prop to GE11 and ME11 = " << num_props_ME11 << std::endl;
  std::cout << "Muons with prop to GE11 and no ME11 = " << num_props_noME11 << std::endl;
}

float RdPhi(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, edm::ESHandle<GEMGeometry> GEMGeometry_, const GEMEtaPartition* ch){
  GEMDetId gemid((rechit)->geographicalId());
  const auto& etaPart = GEMGeometry_->etaPartition(gemid);
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
  return cos(stripAngle) * (prop_localx - (rechit)->localPosition().x()) + sin(stripAngle) * (prop_localy + deltay_roll);
}

void propagate_segment(const CSCSegment* ME11_segment, const reco::Track* innerTrack, const GEMEtaPartition* ch, const reco::Muon* mu, MuonServiceProxy* theService_, const GeomDet* segDet, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool debug){
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  LocalVector momentum_at_surface = ME11_segment->localDirection(); //No momentum for segments
  if (innerTrack != 0){
    std::cout << "Segment using inner momentum! Wow success" << std::endl;
    momentum_at_surface = momentum_at_surface*(innerTrack->outerP()); //If innerTrack exists, use momentum
  }
  LocalTrajectoryParameters param(ME11_segment->localPosition(), momentum_at_surface, mu->charge());
  AlgebraicSymMatrix mat(5,0);
  mat = ME11_segment->parametersError().similarityT( ME11_segment->projectionMatrix() );
  LocalTrajectoryError error(asSMatrix<5>(mat));
  TrajectoryStateOnSurface tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
  TrajectoryStateOnSurface tsos_ch = propagator->propagate(tsos_seg, ch->surface());
  if (tsos_ch.isValid()){
    const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
    const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
    if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
      has_prop = true;
      pos_global_ch = tsos_ch.globalPosition();
      pos_global_seg = tsos_seg.globalPosition();
    }
  }
}

void propagate_track_based_segment(const CSCSegment* ME11_segment, reco::TransientTrack track, const GEMEtaPartition* ch, const reco::Muon* mu, MuonServiceProxy* theService_, const GeomDet* segDet, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool debug){
  if (track.stateOnSurface(segDet->toGlobal(ME11_segment->localPosition())).isValid()){
    const BoundPlane& bps(ch->surface());
    auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
    GlobalVector tracker_momentum_at_surface = track.trajectoryStateClosestToPoint(segDet->toGlobal(ME11_segment->localPosition())).momentum();
    LocalTrajectoryParameters param(ME11_segment->localPosition(), segDet->toLocal(tracker_momentum_at_surface), mu->charge());
    AlgebraicSymMatrix  mat(5,0);
    mat = ME11_segment->parametersError().similarityT( ME11_segment->projectionMatrix() );
    LocalTrajectoryError error(asSMatrix<5>(mat)); //This is not handled correctly yet
    TrajectoryStateOnSurface tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
    TrajectoryStateOnSurface tsos_ch = propagator->propagate(tsos_seg, ch->surface());
    if (tsos_ch.isValid()){
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
        has_prop = true;
        pos_global_ch = tsos_ch.globalPosition();
        pos_global_seg = tsos_seg.globalPosition();
      }
    }
  }
}



void propagate_track(reco::TransientTrack track, const GEMEtaPartition* ch, const string inner_or_CSC, MuonServiceProxy* theService_, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool &which_track, bool debug){
  if ((inner_or_CSC != "CSC") && (inner_or_CSC != "inner")){std::cout << "Propagation failed bad source" << std::endl; return;}
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  TrajectoryStateOnSurface tsos_ch;
  TrajectoryStateOnSurface tsos_seg;

  if (track.outermostMeasurementState().globalPosition().mag2() - track.innermostMeasurementState().globalPosition().mag2() > 0){
    if (inner_or_CSC == "CSC"){tsos_seg = track.innermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());}
    else{tsos_seg = track.outermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());}
    which_track = 1;
  }
  else{
    if (inner_or_CSC == "CSC"){tsos_seg = track.outermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());}
    else{tsos_seg = track.innermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());}
    which_track = 0;
  }
  if (tsos_ch.isValid()){
    const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
    const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
    if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
      has_prop = true;
      pos_global_ch = tsos_ch.globalPosition();
      pos_global_seg = tsos_seg.globalPosition();
    }
  }
}


bool fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch){
  const float fidcut_angle = 1.0;
  const float cut_chamber = 5.0;
  const float cut_angle = 5.0 - fidcut_angle;
  auto& parameters(ch->specs()->parameters());
  float height(parameters[2]);
  if ((abs(localphi_deg) < cut_angle) && ((local_y < (height - cut_chamber) && ch->id().roll() == 1) || (local_y > -1.0*(height - cut_chamber) && ch->id().roll() == 8) || (ch->id().roll() != 1 && ch->id().roll() != 8))){return 1;}
  else{return 0;}
}



void analyser::beginJob(){}
void analyser::endJob(){}

DEFINE_FWK_MODULE(analyser);

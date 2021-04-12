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

  //Prop from tracker
  int prop_region_GE11;
  int prop_station_GE11;
  int prop_layer_GE11;
  int prop_chamber_GE11;
  int prop_roll_GE11;


 
  float prop_inner_x_GE11;
  float prop_inner_y_GE11;
  float prop_inner_r_GE11;
  float prop_inner_localx_GE11;
  float prop_inner_localy_GE11;
  float prop_inner_y_adjusted_GE11;
  float prop_inner_localphi_rad_GE11;
  float prop_inner_localphi_deg_GE11;
  float prop_inner_chi2_GE11;
  float prop_inner_ndof_GE11;
  float prop_inner_chi2ndof_GE11;

  //Prop from CSC
  float prop_CSC_x_GE11;
  float prop_CSC_y_GE11;
  float prop_CSC_r_GE11;
  float prop_CSC_localx_GE11;
  float prop_CSC_localy_GE11;
  float prop_CSC_y_adjusted_GE11;
  float prop_CSC_localphi_rad_GE11;
  float prop_CSC_localphi_deg_GE11;
  float prop_CSC_chi2_GE11;
  float prop_CSC_ndof_GE11;
  float prop_CSC_chi2ndof_GE11;


  bool has_rechit_CSC_GE11;
  int rechit_region_CSC_GE11;
  int rechit_station_CSC_GE11;
  int rechit_layer_CSC_GE11;
  int rechit_chamber_CSC_GE11;
  int rechit_roll_CSC_GE11;
  int recHit_first_strip_CSC;
  int recHit_CLS_CSC;
  float rechit_x_CSC_GE11;
  float rechit_y_CSC_GE11;
  float rechit_r_CSC_GE11;
  float rechit_localx_CSC_GE11;
  float rechit_localy_CSC_GE11;
  float rechit_y_adjusted_CSC_GE11;
  float rechit_localphi_rad_CSC_GE11;
  float rechit_localphi_deg_CSC_GE11;

  bool has_rechit_inner_GE11;
  int rechit_region_inner_GE11;
  int rechit_station_inner_GE11;
  int rechit_layer_inner_GE11;
  int rechit_chamber_inner_GE11;
  int rechit_roll_inner_GE11;
  int recHit_first_strip_inner;
  int recHit_CLS_inner;
  float rechit_x_inner_GE11;
  float rechit_y_inner_GE11;
  float rechit_r_inner_GE11;
  float rechit_localx_inner_GE11;
  float rechit_localy_inner_GE11;
  float rechit_y_adjusted_inner_GE11;
  float rechit_localphi_rad_inner_GE11;
  float rechit_localphi_deg_inner_GE11;

  float RdPhi_inner_GE11;
  float RdPhi_CSC_GE11;
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

  int nCSCSeg;
  int nDTSeg;
  int CSCSeg_region;
  int num_props;

  float simDy;
  int nRecHitsTot;
  int nRecHits5;
  int nRecHits2;

  int closest;
  float startingPoint_x_CSC;
  float startingPoint_y_CSC;
  float startingPoint_z_CSC;
  float startingPoint_r_CSC;

  float startingPoint_x_inner;
  float startingPoint_y_inner;
  float startingPoint_z_inner;
  float startingPoint_r_inner;

  //ME11 track based prop
  float prop_outerSeg_x_GE11;
  float prop_outerSeg_y_GE11;
  float prop_outerSeg_r_GE11;
  float prop_outerSeg_localx_GE11;
  float prop_outerSeg_localy_GE11;
  float prop_outerSeg_y_adjusted_GE11;
  float prop_outerSeg_localphi_rad_GE11;
  float prop_outerSeg_localphi_deg_GE11;

  float prop_innerSeg_x_GE11;
  float prop_innerSeg_y_GE11;
  float prop_innerSeg_r_GE11;
  float prop_innerSeg_localx_GE11;
  float prop_innerSeg_localy_GE11;
  float prop_innerSeg_y_adjusted_GE11;
  float prop_innerSeg_localphi_rad_GE11;
  float prop_innerSeg_localphi_deg_GE11;

  float RdPhi_outerSeg_GE11;
  float RdPhi_innerSeg_GE11;

  int closest_ME11;
  float startingPoint_r_ME11;
  float startingPoint_z_ME11;
  float startingPoint_x_ME11;
  float startingPoint_y_ME11;

  float startingPoint_r_outerSeg;
  float startingPoint_z_outerSeg;
  float startingPoint_x_outerSeg;
  float startingPoint_y_outerSeg;

  float startingPoint_r_innerSeg;
  float startingPoint_z_innerSeg;
  float startingPoint_x_innerSeg;
  float startingPoint_y_innerSeg;

  bool has_prop_inner;
  bool has_prop_CSC;
  bool has_prop_innerSeg;
  bool has_prop_outerSeg;


  //Sim info for MC
  float sim_r;
  float sim_x;
  float sim_y;
  float sim_z;
  float sim_localx;
  float sim_localy;
  float sim_localy_adjusted;
  int nSim;

};

void MuonData::init()
{
  muon_charge = 9999;
  muon_pt = 9999;

  prop_region_GE11 = 99999;
  prop_station_GE11 = 99999;
  prop_layer_GE11 = 99999;
  prop_chamber_GE11 = 99999;
  prop_roll_GE11 = 99999;



  prop_inner_x_GE11 = 99999;
  prop_inner_y_GE11 = 99999;
  prop_inner_r_GE11 = 99999;
  prop_inner_localx_GE11 = 99999;
  prop_inner_localy_GE11 = 99999;
  prop_inner_y_adjusted_GE11 = 99999;
  prop_inner_localphi_rad_GE11 = 99999;
  prop_inner_localphi_deg_GE11 = 99999;
  prop_inner_chi2_GE11 = 99999;
  prop_inner_ndof_GE11 = 99999;
  prop_inner_chi2ndof_GE11 = 99999;

  prop_CSC_x_GE11 = 99999;
  prop_CSC_y_GE11 = 99999;
  prop_CSC_r_GE11 = 99999;
  prop_CSC_localx_GE11 = 99999;
  prop_CSC_localy_GE11 = 99999;
  prop_CSC_y_adjusted_GE11 = 99999;
  prop_CSC_localphi_rad_GE11 = 99999;
  prop_CSC_localphi_deg_GE11 = 99999;
  prop_CSC_chi2_GE11 = 99999;
  prop_CSC_ndof_GE11 = 99999;
  prop_CSC_chi2ndof_GE11 = 99999;


  has_rechit_CSC_GE11 = false;
  rechit_region_CSC_GE11 = 999999;
  rechit_station_CSC_GE11 = 999999;
  rechit_layer_CSC_GE11 = 999999;
  rechit_chamber_CSC_GE11 = 999999;
  rechit_roll_CSC_GE11 = 999999;
  recHit_first_strip_CSC = 999999;
  recHit_CLS_CSC = 999999;
  rechit_x_CSC_GE11 = 999999;
  rechit_y_CSC_GE11 = 999999;
  rechit_r_CSC_GE11 = 999999;
  rechit_localx_CSC_GE11 = 999999;
  rechit_localy_CSC_GE11 = 999999;
  rechit_y_adjusted_CSC_GE11 = 999999;
  rechit_localphi_rad_CSC_GE11 = 999999;
  rechit_localphi_deg_CSC_GE11 = 999999;

  has_rechit_inner_GE11 = false;
  rechit_region_inner_GE11 = 999999;
  rechit_station_inner_GE11 = 999999;
  rechit_layer_inner_GE11 = 999999;
  rechit_chamber_inner_GE11 = 999999;
  rechit_roll_inner_GE11 = 999999;
  recHit_first_strip_inner = 999999;
  recHit_CLS_inner = 999999;
  rechit_x_inner_GE11 = 999999;
  rechit_y_inner_GE11 = 999999;
  rechit_r_inner_GE11 = 999999;
  rechit_localx_inner_GE11 = 999999;
  rechit_localy_inner_GE11 = 999999;
  rechit_y_adjusted_inner_GE11 = 999999;
  rechit_localphi_rad_inner_GE11 = 999999;
  rechit_localphi_deg_inner_GE11 = 999999;


  RdPhi_inner_GE11 = 999999;
  RdPhi_CSC_GE11 = 999999;
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

  nCSCSeg = 0;
  nDTSeg = 0;
  CSCSeg_region = 0;
  num_props = 0;

  simDy = 999;
  nRecHits5 = 0;
  nRecHits2 = 0;
  nRecHitsTot = 0;

  closest = 999;
  startingPoint_x_CSC = 999;
  startingPoint_y_CSC = 999;
  startingPoint_z_CSC = 999;
  startingPoint_r_CSC = 999;

  startingPoint_x_inner = 999;
  startingPoint_y_inner = 999;
  startingPoint_z_inner = 999;
  startingPoint_r_inner = 999;

  //ME11 track based prop
  prop_outerSeg_x_GE11 = 9999999;
  prop_outerSeg_y_GE11 = 9999999;
  prop_outerSeg_r_GE11 = 9999999;
  prop_outerSeg_localx_GE11 = 9999999;
  prop_outerSeg_localy_GE11 = 9999999;
  prop_outerSeg_y_adjusted_GE11 = 9999999;
  prop_outerSeg_localphi_rad_GE11 = 9999999;
  prop_outerSeg_localphi_deg_GE11 = 9999999;

  prop_innerSeg_x_GE11 = 9999999;
  prop_innerSeg_y_GE11 = 9999999;
  prop_innerSeg_r_GE11 = 9999999;
  prop_innerSeg_localx_GE11 = 9999999;
  prop_innerSeg_localy_GE11 = 9999999;
  prop_innerSeg_y_adjusted_GE11 = 9999999;
  prop_innerSeg_localphi_rad_GE11 = 9999999;
  prop_innerSeg_localphi_deg_GE11 = 9999999;

  RdPhi_outerSeg_GE11 = 9999999;
  RdPhi_innerSeg_GE11 = 9999999;

  closest_ME11 = 999;
  startingPoint_r_ME11 = 999;
  startingPoint_z_ME11 = 999;
  startingPoint_x_ME11 = 999;
  startingPoint_y_ME11 = 999;

  startingPoint_r_outerSeg = 999;
  startingPoint_z_outerSeg = 999;
  startingPoint_x_outerSeg = 999;
  startingPoint_y_outerSeg = 999;

  startingPoint_r_innerSeg = 999;
  startingPoint_z_innerSeg = 999;
  startingPoint_x_innerSeg = 999;
  startingPoint_y_innerSeg = 999;

  has_prop_inner = false;
  has_prop_CSC = false;
  has_prop_innerSeg = false;
  has_prop_outerSeg = false;


  //Sim info for MC
  sim_r = 99999999;
  sim_x = 99999999;
  sim_y = 99999999;
  sim_z = 99999999;
  sim_localx = 99999999;
  sim_localy = 99999999;
  sim_localy_adjusted = 99999999;
  nSim = 99999999;
}

TTree* MuonData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("MuonData", "MuonData");

  t->Branch("muon_charge", &muon_charge);
  t->Branch("muon_pt", &muon_pt);
  t->Branch("prop_region_GE11", &prop_region_GE11);
  t->Branch("prop_station_GE11", &prop_station_GE11);
  t->Branch("prop_layer_GE11", &prop_layer_GE11);
  t->Branch("prop_chamber_GE11", &prop_chamber_GE11);
  t->Branch("prop_roll_GE11", &prop_roll_GE11);




//Propagated Inner
  t->Branch("prop_inner_x_GE11", &prop_inner_x_GE11);
  t->Branch("prop_inner_y_GE11", &prop_inner_y_GE11);
  t->Branch("prop_inner_r_GE11", &prop_inner_r_GE11);
  t->Branch("prop_inner_localx_GE11", &prop_inner_localx_GE11);
  t->Branch("prop_inner_localy_GE11", &prop_inner_localy_GE11);
  t->Branch("prop_inner_y_adjusted_GE11", &prop_inner_y_adjusted_GE11);
  t->Branch("prop_inner_localphi_rad_GE11", &prop_inner_localphi_rad_GE11);
  t->Branch("prop_inner_localphi_deg_GE11", &prop_inner_localphi_deg_GE11);
  t->Branch("prop_inner_chi2_GE11", &prop_inner_chi2_GE11);
  t->Branch("prop_inner_ndof_GE11", &prop_inner_ndof_GE11);
  t->Branch("prop_inner_chi2ndof_GE11", &prop_inner_chi2ndof_GE11);

//Propogated CSC
  t->Branch("prop_CSC_x_GE11", &prop_CSC_x_GE11);
  t->Branch("prop_CSC_y_GE11", &prop_CSC_y_GE11);
  t->Branch("prop_CSC_r_GE11", &prop_CSC_r_GE11);
  t->Branch("prop_CSC_localx_GE11", &prop_CSC_localx_GE11);
  t->Branch("prop_CSC_localy_GE11", &prop_CSC_localy_GE11);
  t->Branch("prop_CSC_y_adjusted_GE11", &prop_CSC_y_adjusted_GE11);
  t->Branch("prop_CSC_localphi_rad_GE11", &prop_CSC_localphi_rad_GE11);
  t->Branch("prop_CSC_localphi_deg_GE11", &prop_CSC_localphi_deg_GE11);
  t->Branch("prop_CSC_chi2_GE11", &prop_CSC_chi2_GE11);
  t->Branch("prop_CSC_ndof_GE11", &prop_CSC_ndof_GE11);
  t->Branch("prop_CSC_chi2ndof_GE11", &prop_CSC_chi2ndof_GE11);


//Reconstructed
  t->Branch("has_rechit_CSC_GE11", &has_rechit_CSC_GE11);
  t->Branch("rechit_region_CSC_GE11", &rechit_region_CSC_GE11);
  t->Branch("rechit_station_CSC_GE11", &rechit_station_CSC_GE11);
  t->Branch("rechit_layer_CSC_GE11", &rechit_layer_CSC_GE11);
  t->Branch("rechit_chamber_CSC_GE11", &rechit_chamber_CSC_GE11);
  t->Branch("rechit_roll_CSC_GE11", &rechit_roll_CSC_GE11);
  t->Branch("recHit_first_strip_CSC", &recHit_first_strip_CSC);
  t->Branch("recHit_CLS_CSC", &recHit_CLS_CSC);
  t->Branch("rechit_x_CSC_GE11", &rechit_x_CSC_GE11);
  t->Branch("rechit_y_CSC_GE11", &rechit_y_CSC_GE11);
  t->Branch("rechit_r_CSC_GE11", &rechit_r_CSC_GE11);
  t->Branch("rechit_localx_CSC_GE11", &rechit_localx_CSC_GE11);
  t->Branch("rechit_localy_CSC_GE11", &rechit_localy_CSC_GE11);
  t->Branch("rechit_y_adjusted_CSC_GE11", &rechit_y_adjusted_CSC_GE11);
  t->Branch("rechit_localphi_rad_CSC_GE11", &rechit_localphi_rad_CSC_GE11);
  t->Branch("rechit_localphi_deg_CSC_GE11", &rechit_localphi_deg_CSC_GE11);

  t->Branch("has_rechit_inner_GE11", &has_rechit_inner_GE11);
  t->Branch("rechit_region_inner_GE11", &rechit_region_inner_GE11);
  t->Branch("rechit_station_inner_GE11", &rechit_station_inner_GE11);
  t->Branch("rechit_layer_inner_GE11", &rechit_layer_inner_GE11);
  t->Branch("rechit_chamber_inner_GE11", &rechit_chamber_inner_GE11);
  t->Branch("rechit_roll_inner_GE11", &rechit_roll_inner_GE11);
  t->Branch("recHit_first_strip_inner", &recHit_first_strip_inner);
  t->Branch("recHit_CLS_inner", &recHit_CLS_inner);
  t->Branch("rechit_x_inner_GE11", &rechit_x_inner_GE11);
  t->Branch("rechit_y_inner_GE11", &rechit_y_inner_GE11);
  t->Branch("rechit_r_inner_GE11", &rechit_r_inner_GE11);
  t->Branch("rechit_localx_inner_GE11", &rechit_localx_inner_GE11);
  t->Branch("rechit_localy_inner_GE11", &rechit_localy_inner_GE11);
  t->Branch("rechit_y_adjusted_inner_GE11", &rechit_y_adjusted_inner_GE11);
  t->Branch("rechit_localphi_rad_inner_GE11", &rechit_localphi_rad_inner_GE11);
  t->Branch("rechit_localphi_deg_inner_GE11", &rechit_localphi_deg_inner_GE11);

//Residual
  t->Branch("RdPhi_inner_GE11", &RdPhi_inner_GE11);
  t->Branch("RdPhi_CSC_GE11", &RdPhi_CSC_GE11);
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

  t->Branch("nCSCSeg", &nCSCSeg); 
  t->Branch("nDTSeg", &nDTSeg); 
  t->Branch("CSCSeg_region", &CSCSeg_region);
  t->Branch("num_props", &num_props);

  t->Branch("simDy", &simDy);
  t->Branch("nRecHits5", &nRecHits5);
  t->Branch("nRecHits2", &nRecHits2);
  t->Branch("nRecHitsTot", &nRecHitsTot);

  t->Branch("closest", &closest);
  t->Branch("startingPoint_x_CSC", &startingPoint_x_CSC);
  t->Branch("startingPoint_y_CSC", &startingPoint_y_CSC);
  t->Branch("startingPoint_z_CSC", &startingPoint_z_CSC);
  t->Branch("startingPoint_r_CSC", &startingPoint_r_CSC);

  t->Branch("startingPoint_x_inner", &startingPoint_x_inner);
  t->Branch("startingPoint_y_inner", &startingPoint_y_inner);
  t->Branch("startingPoint_z_inner", &startingPoint_z_inner);
  t->Branch("startingPoint_r_inner", &startingPoint_r_inner);

  //ME11 trackbased prop
  t->Branch("prop_outerSeg_x_GE11", &prop_outerSeg_x_GE11);
  t->Branch("prop_outerSeg_y_GE11", &prop_outerSeg_y_GE11);
  t->Branch("prop_outerSeg_r_GE11", &prop_outerSeg_r_GE11);
  t->Branch("prop_outerSeg_localx_GE11", &prop_outerSeg_localx_GE11);
  t->Branch("prop_outerSeg_localy_GE11", &prop_outerSeg_localy_GE11);
  t->Branch("prop_outerSeg_y_adjusted_GE11", &prop_outerSeg_y_adjusted_GE11);
  t->Branch("prop_outerSeg_localphi_rad_GE11", &prop_outerSeg_localphi_rad_GE11);
  t->Branch("prop_outerSeg_localphi_deg_GE11", &prop_outerSeg_localphi_deg_GE11);

  t->Branch("prop_innerSeg_x_GE11", &prop_innerSeg_x_GE11);
  t->Branch("prop_innerSeg_y_GE11", &prop_innerSeg_y_GE11);
  t->Branch("prop_innerSeg_r_GE11", &prop_innerSeg_r_GE11);
  t->Branch("prop_innerSeg_localx_GE11", &prop_innerSeg_localx_GE11);
  t->Branch("prop_innerSeg_localy_GE11", &prop_innerSeg_localy_GE11);
  t->Branch("prop_innerSeg_y_adjusted_GE11", &prop_innerSeg_y_adjusted_GE11);
  t->Branch("prop_innerSeg_localphi_rad_GE11", &prop_innerSeg_localphi_rad_GE11);
  t->Branch("prop_innerSeg_localphi_deg_GE11", &prop_innerSeg_localphi_deg_GE11);

  t->Branch("RdPhi_outerSeg_GE11", &RdPhi_outerSeg_GE11);
  t->Branch("RdPhi_innerSeg_GE11", &RdPhi_innerSeg_GE11);

  t->Branch("closest_ME11", &closest_ME11);
  t->Branch("startingPoint_r_ME11", &startingPoint_r_ME11);
  t->Branch("startingPoint_z_ME11", &startingPoint_z_ME11);
  t->Branch("startingPoint_x_ME11", &startingPoint_x_ME11);
  t->Branch("startingPoint_y_ME11", &startingPoint_y_ME11);

  t->Branch("startingPoint_r_outerSeg", &startingPoint_r_outerSeg);
  t->Branch("startingPoint_z_outerSeg", &startingPoint_z_outerSeg);
  t->Branch("startingPoint_x_outerSeg", &startingPoint_x_outerSeg);
  t->Branch("startingPoint_y_outerSeg", &startingPoint_y_outerSeg);

  t->Branch("startingPoint_r_innerSeg", &startingPoint_r_innerSeg);
  t->Branch("startingPoint_z_innerSeg", &startingPoint_z_innerSeg);
  t->Branch("startingPoint_x_innerSeg", &startingPoint_x_innerSeg);
  t->Branch("startingPoint_y_innerSeg", &startingPoint_y_innerSeg);

  t->Branch("has_prop_inner", &has_prop_inner);
  t->Branch("has_prop_CSC", &has_prop_CSC);
  t->Branch("has_prop_innerSeg", &has_prop_innerSeg);
  t->Branch("has_prop_outerSeg", &has_prop_outerSeg);


  //Sim info for MC
  t->Branch("sim_r", &sim_r);
  t->Branch("sim_x", &sim_x);
  t->Branch("sim_y", &sim_y);
  t->Branch("sim_z", &sim_z);
  t->Branch("sim_localx", &sim_localx);
  t->Branch("sim_localy", &sim_localy);
  t->Branch("sim_localy_adjusted", &sim_localy_adjusted);
  t->Branch("nSim", &nSim);

  return t;
}

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
  cout << "new evt numb is " << iEvent.eventAuxiliary().event() << endl;
  for (size_t i = 0; i < muons->size(); ++i){
    cout << "new muon" << endl;
    //cout << "evt number is " << iEvent.eventAuxiliary().event() << endl;
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    cout << "is standalone" << endl;
    data_.init();
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
        if (debug)std::cout << "Looping over segments" << std::endl;
        auto cscSegRef = MSM.cscSegmentRef;
        auto cscDetID = cscSegRef->cscDetId();
        data_.CSCSeg_region = cscDetID.endcap();
        if (cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
          if (debug && data_.hasME11 == 1) std::cout << "Already has an ME11 seg" << std::endl;
          data_.hasME11 = 1;
          if (debug)std::cout << "has ME11 segment" << std::endl;
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
	if (debug)std::cout << "Closest now says " << data_.closest << std::endl;


        if (cscDetID.station() == 1 and cscDetID.ring() == 4) data_.hasME11A = 1;
      }
    }
    
    data_.evtNum = iEvent.eventAuxiliary().event();
    data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock();

    //Get det info
    //int chamber_count = 0;
    //for (const auto& ch : GEMGeometry_->chambers()){
    //  chamber_count ++;
    //}
    //cout << "chamber count = " << chamber_count << endl;


    //Get the tracks
    if (debug)std::cout << "Getting tracks" << std::endl;
    const reco::Track* globalTrack = 0;
    const reco::Track* innerTrack = 0;
    const reco::Track* outerTrack = 0;
    reco::TransientTrack ttTrack_global;
    reco::TransientTrack ttTrack_tracker;
    reco::TransientTrack ttTrack_CSC;
    if ( mu->globalTrack().isNonnull() ){ globalTrack = mu->globalTrack().get(); ttTrack_global = ttrackBuilder_->build(globalTrack);}
    if ( mu->track().isNonnull() ){ innerTrack = mu->track().get(); ttTrack_tracker = ttrackBuilder_->build(innerTrack);}
    if ( mu->outerTrack().isNonnull() ){ outerTrack = mu->outerTrack().get(); ttTrack_CSC = ttrackBuilder_->build(outerTrack);}


    //Start the propagations
    float count = 0;
    if (debug)std::cout << "Starting chamber loop" << std::endl;
    for (const auto& ch : GEMGeometry_->etaPartitions()) {
      if (ch->id().station() != 1) continue; //Only takes GE1/1
      const BoundPlane& bps(ch->surface());

      //Props from ME11 segment
      data_.has_prop_innerSeg = 0;
      data_.has_prop_outerSeg = 0;
      TrajectoryStateOnSurface tsos_ME11trk_outer;
      TrajectoryStateOnSurface tsos_ME11trk_inner;
      LocalPoint pos_local_outerSeg;
      LocalPoint pos_local_innerSeg;
      GlobalPoint pos_global_outerSeg;
      GlobalPoint pos_global_innerSeg;
      if (data_.hasME11 == 1 and ch->id().station() == 1 and ch->id().ring() == 1){
        if (debug)std::cout << "Doing segment propagations" << std::endl; 
        DetId segDetId = tmp_ME11_seg->geographicalId();
        const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
        data_.startingPoint_z_ME11 = segDet->toGlobal(tmp_ME11_seg->localPosition()).z();
        data_.startingPoint_x_ME11 = segDet->toGlobal(tmp_ME11_seg->localPosition()).x();
        data_.startingPoint_y_ME11 = segDet->toGlobal(tmp_ME11_seg->localPosition()).y();
        data_.startingPoint_r_ME11 = pow(pow(data_.startingPoint_x_ME11,2) + pow(data_.startingPoint_y_ME11,2), 0.5);
        if (debug)std::cout << "Got segment starting position" << std::endl;
        

        //We noticed a problem with the segment position and the stateOnSurface position not being equal (apr 8 2021)
        //Attempting to move back to generating our own TSOS by hand with track momentum
        //Using commit from Feb3 to remember how to create our own TSOS
        //
        //LocalTrajectoryParameters param(tmp_ME11_seg->localPosition(), tmp_ME11_seg->localDirection(), mu->charge());
        //AlgebraicSymMatrix mat(5,0);
        //mat = tmp_ME11_seg->parametersError().similarityT( tmp_ME11_seg->projectionMatrix() );
        //LocalTrajectoryError error(asSMatrix<5>(mat));
        //DetId segDetId = tmp_ME11_seg->geographicalId(); //We Have This Already
        //const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId); //We Have This Already
        //TrajectoryStateOnSurface tsos_ME11seg_GE11(param, error, segDet->surface(), &*theService_->magneticField());
        //tsos_CSC_ME11seg = propagator->propagate(tsos_ME11seg_GE11, ch->surface());
        //
        //Issue is now using the segment position but track direction and momentum -- what is the projectionMatrix?



        //Track propagation starting at ME11 segment location
        if (tracker_prop and ttTrack_tracker.isValid()){
          if (ttTrack_tracker.stateOnSurface(segDet->toGlobal(tmp_ME11_seg->localPosition())).isValid()){  //Had segfaults for not-valid trajectory states
            //Working on updating how we make TSOS
            //TrajectoryStateOnSurface innerSeg_TSOS = ttTrack_tracker.stateOnSurface(segDet->toGlobal(tmp_ME11_seg->localPosition()));

            //Begin new test
            //To Do
            //Switch from local to global coordinates and build through that (for compatability with tracks)
            //
            //LocalTrajectoryParameters - previously we had (position, direction, charge) but this seems wrong https://github.com/cms-sw/cmssw/blob/ba6e8604a35283e39e89bc031766843d0afc3240/DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h#L70
            //Here it shows that the second piece should be the momentum, but I don't think segments have momentum info
            //
            //This seems to be working, but I need to ask about the error matrix still
            if (debug)std::cout << "Testing momentums! Wow fun stuff! " << tmp_ME11_seg->localDirection() << std::endl;
            GlobalVector tracker_momentum_at_surface = ttTrack_tracker.trajectoryStateClosestToPoint(segDet->toGlobal(tmp_ME11_seg->localPosition())).momentum();
            LocalTrajectoryParameters param(tmp_ME11_seg->localPosition(), segDet->toLocal(tracker_momentum_at_surface), mu->charge());
            AlgebraicSymMatrix  mat(5,0);
            mat = tmp_ME11_seg->parametersError().similarityT( tmp_ME11_seg->projectionMatrix() );
            LocalTrajectoryError error(asSMatrix<5>(mat));
            TrajectoryStateOnSurface innerSeg_TSOS(param, error, segDet->surface(), &*theService_->magneticField());

            if (debug)std::cout << "Momentum from the tracks closest measurement" << tracker_momentum_at_surface << std::endl;
            //End new test

            tsos_ME11trk_inner = propagator->propagate(innerSeg_TSOS, ch->surface());
            if (debug)std::cout << "inner trkprop check" << std::endl;
            if (tsos_ME11trk_inner.isValid()){
              pos_global_innerSeg = tsos_ME11trk_inner.globalPosition();
              pos_local_innerSeg = ch->toLocal(tsos_ME11trk_inner.globalPosition());
              const GlobalPoint pos2D_global_innerSeg(pos_global_innerSeg.x(), pos_global_innerSeg.y(), 0);
              const LocalPoint pos2D_local_innerSeg(pos_local_innerSeg.x(), pos_local_innerSeg.y(), 0);
              if (!(pos_global_innerSeg.eta() * mu->eta() < 0.0) and bps.bounds().inside(pos2D_local_innerSeg) and ch->id().station() == 1 and ch->id().ring() == 1){
                data_.has_prop_innerSeg = true;
                float startingPoint_x_innerSeg;
                float startingPoint_y_innerSeg;
                float startingPoint_z_innerSeg;
                float startingPoint_r_innerSeg;
                GlobalPoint startingPoint_GP_innerSeg = innerSeg_TSOS.globalPosition();
                startingPoint_x_innerSeg = startingPoint_GP_innerSeg.x();
                startingPoint_y_innerSeg = startingPoint_GP_innerSeg.y();
                startingPoint_z_innerSeg = startingPoint_GP_innerSeg.z();
                startingPoint_r_innerSeg = pow(pow(startingPoint_x_innerSeg,2) + pow(startingPoint_y_innerSeg,2), 0.5);
                data_.startingPoint_x_innerSeg = startingPoint_x_innerSeg;
                data_.startingPoint_y_innerSeg = startingPoint_y_innerSeg;
                data_.startingPoint_z_innerSeg = startingPoint_z_innerSeg;
                data_.startingPoint_r_innerSeg = startingPoint_r_innerSeg;

                const MagneticField* innerSeg_magField = innerSeg_TSOS.magneticField();
                std::cout << "Mag field looks like " << innerSeg_magField << std::endl;
              }
            }
          }
        }
        //CSC propagation starting at ME11 segment location
        if (CSC_prop and ttTrack_CSC.isValid()){
          if (ttTrack_CSC.stateOnSurface(segDet->toGlobal(tmp_ME11_seg->localPosition())).isValid()){  //Had segfaults for not-valid trajectory states
            //Testing new method for TSOS
            //TrajectoryStateOnSurface outerSeg_TSOS = ttTrack_CSC.stateOnSurface(segDet->toGlobal(tmp_ME11_seg->localPosition()));

            if (debug)std::cout << "Testing momentums! Wow fun stuff! " << tmp_ME11_seg->localDirection() << std::endl;
            GlobalVector tracker_momentum_at_surface = ttTrack_CSC.trajectoryStateClosestToPoint(segDet->toGlobal(tmp_ME11_seg->localPosition())).momentum();
            LocalTrajectoryParameters param(tmp_ME11_seg->localPosition(), segDet->toLocal(tracker_momentum_at_surface), mu->charge());
            AlgebraicSymMatrix  mat(5,0);
            mat = tmp_ME11_seg->parametersError().similarityT( tmp_ME11_seg->projectionMatrix() );
            LocalTrajectoryError error(asSMatrix<5>(mat));
            TrajectoryStateOnSurface outerSeg_TSOS(param, error, segDet->surface(), &*theService_->magneticField());

            if (debug)std::cout << "Momentum from the tracks closest measurement" << tracker_momentum_at_surface << std::endl;
            //End new test


            tsos_ME11trk_outer = propagator->propagate(outerSeg_TSOS, ch->surface());
            if (debug)std::cout << "outer trkprop check" << std::endl;
            if (tsos_ME11trk_outer.isValid()){
              pos_global_outerSeg = tsos_ME11trk_outer.globalPosition();
              pos_local_outerSeg = ch->toLocal(tsos_ME11trk_outer.globalPosition());
              const GlobalPoint pos2D_global_outerSeg(pos_global_outerSeg.x(), pos_global_outerSeg.y(), 0);
              const LocalPoint pos2D_local_outerSeg(pos_local_outerSeg.x(), pos_local_outerSeg.y(), 0);
              if (!(pos_global_outerSeg.eta() * mu->eta() < 0.0) and bps.bounds().inside(pos2D_local_outerSeg) and ch->id().station() == 1 and ch->id().ring() == 1){
                data_.has_prop_outerSeg = true;
                float startingPoint_x_outerSeg;
                float startingPoint_y_outerSeg;
                float startingPoint_z_outerSeg;
                float startingPoint_r_outerSeg;
                GlobalPoint startingPoint_GP_outerSeg = outerSeg_TSOS.globalPosition();
                startingPoint_x_outerSeg = startingPoint_GP_outerSeg.x();
                startingPoint_y_outerSeg = startingPoint_GP_outerSeg.y();
                startingPoint_z_outerSeg = startingPoint_GP_outerSeg.z();
                startingPoint_r_outerSeg = pow(pow(startingPoint_x_outerSeg,2) + pow(startingPoint_y_outerSeg,2), 0.5);
                data_.startingPoint_x_outerSeg = startingPoint_x_outerSeg;
                data_.startingPoint_y_outerSeg = startingPoint_y_outerSeg;
                data_.startingPoint_z_outerSeg = startingPoint_z_outerSeg;
                data_.startingPoint_r_outerSeg = startingPoint_r_outerSeg;
              }
            }
          }
        }
      }
 



      // Tracker propagated
      TrajectoryStateOnSurface tsos_inner;
      GlobalPoint pos_global_inner;
      LocalPoint pos_local_inner;
      data_.has_prop_inner = false;
      if(tracker_prop and ttTrack_tracker.isValid()){
        if(debug)std::cout << "Starting tracker prop" << std::endl;
        float startingPoint_x_inner;
        float startingPoint_y_inner;
        float startingPoint_z_inner;
        float startingPoint_r_inner;
        GlobalPoint startingPoint_GP_inner;
        if ( innerTrack->outerPosition().Mag2() - innerTrack->innerPosition().Mag2() > 0){
          tsos_inner = propagator->propagate(ttTrack_tracker.outermostMeasurementState(),ch->surface());
          data_.which_track_inner_GE11 = 1;
          startingPoint_GP_inner = ttTrack_tracker.outermostMeasurementState().globalPosition();
        }
        else{
          tsos_inner = propagator->propagate(ttTrack_tracker.innermostMeasurementState(),ch->surface());
          data_.which_track_inner_GE11 = 0;
          startingPoint_GP_inner = ttTrack_tracker.innermostMeasurementState().globalPosition();
        }
        startingPoint_x_inner = startingPoint_GP_inner.x();
        startingPoint_y_inner = startingPoint_GP_inner.y();
        startingPoint_z_inner = startingPoint_GP_inner.z();
        startingPoint_r_inner = pow(pow(startingPoint_x_inner,2) + pow(startingPoint_y_inner,2), 0.5);

        data_.startingPoint_x_inner = startingPoint_x_inner;
        data_.startingPoint_y_inner = startingPoint_y_inner;
        data_.startingPoint_z_inner = startingPoint_z_inner;
        data_.startingPoint_r_inner = startingPoint_r_inner;
        if (tsos_inner.isValid()){
          pos_global_inner = tsos_inner.globalPosition();
          pos_local_inner = ch->toLocal(tsos_inner.globalPosition());
          const GlobalPoint pos2D_global_inner(pos_global_inner.x(), pos_global_inner.y(), 0);
          const LocalPoint pos2D_local_inner(pos_local_inner.x(), pos_local_inner.y(), 0);
          if (!(pos_global_inner.eta() * mu->eta() < 0.0) and bps.bounds().inside(pos2D_local_inner) and ch->id().station() == 1 and ch->id().ring() == 1){
            data_.has_prop_inner = true;
          }
        }
      }

      // CSC propagated
      TrajectoryStateOnSurface tsos_CSC;
      GlobalPoint pos_global_CSC;
      LocalPoint pos_local_CSC;
      data_.has_prop_CSC = false;
      if(CSC_prop and ttTrack_CSC.isValid()){
        if(debug)std::cout << "Starting CSC prop" << std::endl;
        float startingPoint_x_CSC;
        float startingPoint_y_CSC;
        float startingPoint_z_CSC;
        float startingPoint_r_CSC;
        GlobalPoint startingPoint_GP_CSC;

        //if ( outerTrack->outerPosition().Mag2() - outerTrack->innerPosition().Mag2() > 0){
        if ( outerTrack->outerPosition().Mag2() - outerTrack->innerPosition().Mag2() < 0){ //Flipped test!!! Now uses farthest most position
          tsos_CSC = propagator->propagate(ttTrack_CSC.innermostMeasurementState(),ch->surface());
          //tsos_CSC = propagator_opposite->propagate(ttTrack_CSC.innermostMeasurementState(),ch->surface());
          data_.which_track_CSC_GE11 = 1;
          startingPoint_GP_CSC = ttTrack_CSC.innermostMeasurementState().globalPosition();
        }
        else{
          tsos_CSC = propagator->propagate(ttTrack_CSC.outermostMeasurementState(),ch->surface());
          //tsos_CSC = propagator_opposite->propagate(ttTrack_CSC.outermostMeasurementState(),ch->surface());
          data_.which_track_CSC_GE11 = 0;
          startingPoint_GP_CSC = ttTrack_CSC.outermostMeasurementState().globalPosition();
        }
        startingPoint_x_CSC = startingPoint_GP_CSC.x();
        startingPoint_y_CSC = startingPoint_GP_CSC.y();
        startingPoint_z_CSC = startingPoint_GP_CSC.z();
        startingPoint_r_CSC = pow(pow(startingPoint_x_CSC,2) + pow(startingPoint_y_CSC,2), 0.5);

        data_.startingPoint_x_CSC = startingPoint_x_CSC;
        data_.startingPoint_y_CSC = startingPoint_y_CSC;
        data_.startingPoint_z_CSC = startingPoint_z_CSC;
        data_.startingPoint_r_CSC = startingPoint_r_CSC;
        if (tsos_CSC.isValid()){
          pos_global_CSC = tsos_CSC.globalPosition();
          pos_local_CSC = ch->toLocal(tsos_CSC.globalPosition());
          const GlobalPoint pos2D_global_CSC(pos_global_CSC.x(), pos_global_CSC.y(), 0);
          const LocalPoint pos2D_local_CSC(pos_local_CSC.x(), pos_local_CSC.y(), 0);
          if (!(pos_global_CSC.eta() * mu->eta() < 0.0) and bps.bounds().inside(pos2D_local_CSC) and ch->id().station() == 1 and ch->id().ring() == 1){
            data_.has_prop_CSC = true;
          }
        }
      }
      //Skip chambers with no props
      if (!(data_.has_prop_CSC or data_.has_prop_inner or data_.has_prop_outerSeg or data_.has_prop_innerSeg)) continue;

      if(debug)cout << "charge is " << mu->charge() << endl;
      data_.num_props++;
      data_.muon_charge = mu->charge();
      data_.muon_pt = mu->pt();

      const float fidcut_angle = 1.0;
      const float cut_ang = 5.0 - fidcut_angle;
      const float cut_chamber = 5.0;

      const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());

      const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp();
      count++;

      data_.prop_region_GE11 = ch->id().region();
      data_.prop_station_GE11 = ch->id().station();
      data_.prop_layer_GE11 = ch->id().layer();
      data_.prop_chamber_GE11 = ch->id().chamber();
      data_.prop_roll_GE11 = ch->id().roll();


      // CSC prop
      if (data_.has_prop_CSC){
        cout << "pos_local_CSC = " << pos_local_CSC << endl;
        data_.prop_CSC_chi2_GE11 = outerTrack->chi2();
        data_.prop_CSC_ndof_GE11 = outerTrack->ndof();
        data_.prop_CSC_chi2ndof_GE11 = data_.prop_CSC_chi2_GE11/data_.prop_CSC_ndof_GE11;
        LocalPoint local_to_center_CSC(pos_local_CSC.x(), prop_y_to_center + pos_local_CSC.y(), 0);
        const float prop_CSC_localphi_rad = (3.14159265/2.) - local_to_center_CSC.phi();
        const float prop_CSC_localphi_deg = prop_CSC_localphi_rad*180/3.14169265;


        data_.prop_CSC_x_GE11 = pos_global_CSC.x();
        data_.prop_CSC_y_GE11 = pos_global_CSC.y();
        data_.prop_CSC_r_GE11 = pos_global_CSC.mag();
        data_.prop_CSC_localx_GE11 = pos_local_CSC.x();
        data_.prop_CSC_localy_GE11 = pos_local_CSC.y();
        data_.prop_CSC_y_adjusted_GE11 = prop_y_to_center + pos_local_CSC.y();
        data_.prop_CSC_localphi_rad_GE11 = prop_CSC_localphi_rad;
        data_.prop_CSC_localphi_deg_GE11 = prop_CSC_localphi_deg;

        auto& parameters(ch->specs()->parameters());
        float height(parameters[2]);

        if ((abs(prop_CSC_localphi_deg) < cut_ang && (pos_local_CSC.y()) < (height - cut_chamber) && ch->id().roll() == 1) || (abs(prop_CSC_localphi_deg) < cut_ang && (pos_local_CSC.y()) > (height - cut_chamber) && ch->id().roll() == 8)){
          data_.has_fidcut_CSC_GE11 = true;
        }
        else{
          data_.has_fidcut_CSC_GE11 = false;
        }
      }

      // inner prop
      if (data_.has_prop_inner){
        cout << "pos_local_inner = " << pos_local_inner << endl;
        data_.prop_inner_chi2_GE11 = innerTrack->chi2();
        data_.prop_inner_ndof_GE11 = innerTrack->ndof();
        data_.prop_inner_chi2ndof_GE11 = data_.prop_inner_chi2_GE11/data_.prop_inner_ndof_GE11;
        LocalPoint local_to_center_inner(pos_local_inner.x(), prop_y_to_center + pos_local_inner.y(), 0);
        const float prop_inner_localphi_rad = (3.14159265/2.) - local_to_center_inner.phi();
        const float prop_inner_localphi_deg = prop_inner_localphi_rad*180/3.14169265;


        data_.prop_inner_x_GE11 = pos_global_inner.x();
        data_.prop_inner_y_GE11 = pos_global_inner.y();
        data_.prop_inner_r_GE11 = pos_global_inner.mag();
        data_.prop_inner_localx_GE11 = pos_local_inner.x();
        data_.prop_inner_localy_GE11 = pos_local_inner.y();
        data_.prop_inner_y_adjusted_GE11 = prop_y_to_center + pos_local_inner.y();
        data_.prop_inner_localphi_rad_GE11 = prop_inner_localphi_rad;
        data_.prop_inner_localphi_deg_GE11 = prop_inner_localphi_deg;

        auto& parameters(ch->specs()->parameters());
        float height(parameters[2]);

        if ((abs(prop_inner_localphi_deg) < cut_ang && (pos_local_inner.y()) < (height - cut_chamber) && ch->id().roll() == 1) || (abs(prop_inner_localphi_deg) < cut_ang && (pos_local_inner.y()) > (height - cut_chamber) && ch->id().roll() == 8)){
          data_.has_fidcut_inner_GE11 = true;
        }
        else{
          data_.has_fidcut_inner_GE11 = false;
        }
      }



      //ME11 prop part
      if (data_.has_prop_outerSeg){
        cout << "pos_local_outerSeg = " << pos_local_outerSeg << endl;
        GlobalPoint pos_global_outerSeg = tsos_ME11trk_outer.globalPosition();

        LocalPoint local_to_center_outerSeg(pos_local_outerSeg.x(), prop_y_to_center + pos_local_outerSeg.y(), 0);
        const float prop_outerSeg_localphi_rad = (3.14159265/2.) - local_to_center_outerSeg.phi();
        const float prop_outerSeg_localphi_deg = prop_outerSeg_localphi_rad*180/3.14169265;


        data_.prop_outerSeg_x_GE11 = pos_global_outerSeg.x();
        data_.prop_outerSeg_y_GE11 = pos_global_outerSeg.y();
        data_.prop_outerSeg_r_GE11 = pos_global_outerSeg.mag();
        data_.prop_outerSeg_localx_GE11 = pos_local_outerSeg.x();
        data_.prop_outerSeg_localy_GE11 = pos_local_outerSeg.y();
        data_.prop_outerSeg_y_adjusted_GE11 = prop_y_to_center + pos_local_outerSeg.y();
        data_.prop_outerSeg_localphi_rad_GE11 = prop_outerSeg_localphi_rad;
        data_.prop_outerSeg_localphi_deg_GE11 = prop_outerSeg_localphi_deg;
      }

      if (data_.has_prop_innerSeg){
        cout << "pos_local_innerSeg = " << pos_local_innerSeg << endl;
        GlobalPoint pos_global_innerSeg = tsos_ME11trk_inner.globalPosition();

        LocalPoint local_to_center_innerSeg(pos_local_innerSeg.x(), prop_y_to_center + pos_local_innerSeg.y(), 0);
        const float prop_innerSeg_localphi_rad = (3.14159265/2.) - local_to_center_innerSeg.phi();
        const float prop_innerSeg_localphi_deg = prop_innerSeg_localphi_rad*180/3.14169265;


        data_.prop_innerSeg_x_GE11 = pos_global_innerSeg.x();
        data_.prop_innerSeg_y_GE11 = pos_global_innerSeg.y();
        data_.prop_innerSeg_r_GE11 = pos_global_innerSeg.mag();
        data_.prop_innerSeg_localx_GE11 = pos_local_innerSeg.x();
        data_.prop_innerSeg_localy_GE11 = pos_local_innerSeg.y();
        data_.prop_innerSeg_y_adjusted_GE11 = prop_y_to_center + pos_local_innerSeg.y();
        data_.prop_innerSeg_localphi_rad_GE11 = prop_innerSeg_localphi_rad;
        data_.prop_innerSeg_localphi_deg_GE11 = prop_innerSeg_localphi_deg;
      }
      //////////////////




        
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
          //if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and gemid.region() == ch->id().region()){
          //if (gemid.layer() == ch->id().layer() and gemid.region() == ch->id().region()){
            cout << "starting rechit" << endl;
            tmpNRHT++;
            const auto& etaPart = GEMGeometry_->etaPartition(gemid);
            float strip = etaPart->strip(hit->localPosition());
            float stripAngle = etaPart->specificTopology().stripAngle(strip);
            float cosAngle = cos(stripAngle);
            float sinAngle = sin(stripAngle);

            float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
            LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
            float rechit_localphi_rad = (3.14159265/2.) - local_to_center.phi();
            float rechit_localphi_deg = rechit_localphi_rad*180/3.14159265;
            float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();


            if (ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - pos_local_CSC.x()) < 999.0){
              if (abs(cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll)) < 5) tmpNRH5++;
              if (abs(cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll)) < 2) tmpNRH2++;
              data_.det_id = gemid.region()*(gemid.station()*100 + gemid.chamber());

              //CSC matcher
              if (data_.has_prop_CSC and abs(data_.RdPhi_CSC_GE11) > abs(cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll))){
                rechit_matches++;
                std::cout << "Overwrite" << std::endl;

                data_.has_rechit_CSC_GE11 = true;
                data_.rechit_region_CSC_GE11 = gemid.region();
                data_.rechit_station_CSC_GE11 = gemid.station();
                data_.rechit_layer_CSC_GE11 = gemid.layer();
                data_.rechit_chamber_CSC_GE11 = gemid.chamber();
                data_.rechit_roll_CSC_GE11 = gemid.roll();
                data_.recHit_first_strip_CSC = (hit)->firstClusterStrip();
                data_.recHit_CLS_CSC = (hit)->clusterSize();
                data_.rechit_x_CSC_GE11 = etaPart->toGlobal((hit)->localPosition()).x();
                data_.rechit_y_CSC_GE11 = etaPart->toGlobal((hit)->localPosition()).y();
                data_.rechit_r_CSC_GE11 = etaPart->toGlobal((hit)->localPosition()).mag();
                data_.rechit_localx_CSC_GE11 = (hit)->localPosition().x();
                data_.rechit_localy_CSC_GE11 = (hit)->localPosition().y();
                data_.rechit_y_adjusted_CSC_GE11 = rechit_y_to_center + (hit)->localPosition().y();
                data_.rechit_localphi_rad_CSC_GE11 = rechit_localphi_rad;
                data_.rechit_localphi_deg_CSC_GE11 = rechit_localphi_deg;
                data_.RdPhi_CSC_GE11 = cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll);
              }
              //inner matcher
              if (data_.has_prop_inner and abs(data_.RdPhi_inner_GE11) > abs(cosAngle * (pos_local_inner.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_inner.y() + deltay_roll))){
                std::cout << "Overwrite" << std::endl;

                data_.has_rechit_inner_GE11 = true;
                data_.rechit_region_inner_GE11 = gemid.region();
                data_.rechit_station_inner_GE11 = gemid.station();
                data_.rechit_layer_inner_GE11 = gemid.layer();
                data_.rechit_chamber_inner_GE11 = gemid.chamber();
                data_.rechit_roll_inner_GE11 = gemid.roll();
                data_.recHit_first_strip_inner = (hit)->firstClusterStrip();
                data_.recHit_CLS_inner = (hit)->clusterSize();
                data_.rechit_x_inner_GE11 = etaPart->toGlobal((hit)->localPosition()).x();
                data_.rechit_y_inner_GE11 = etaPart->toGlobal((hit)->localPosition()).y();
                data_.rechit_r_inner_GE11 = etaPart->toGlobal((hit)->localPosition()).mag();
                data_.rechit_localx_inner_GE11 = (hit)->localPosition().x();
                data_.rechit_localy_inner_GE11 = (hit)->localPosition().y();
                data_.rechit_y_adjusted_inner_GE11 = rechit_y_to_center + (hit)->localPosition().y();
                data_.rechit_localphi_rad_inner_GE11 = rechit_localphi_rad;
                data_.rechit_localphi_deg_inner_GE11 = rechit_localphi_deg;
                data_.RdPhi_inner_GE11 = cosAngle * (pos_local_inner.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_inner.y() + deltay_roll);
              }
              if(debug)cout << "Starting ME11 rechit match" << endl;
              //ME11 seg matcher
              if (data_.has_prop_outerSeg){
                if (abs(data_.RdPhi_outerSeg_GE11) > abs(cosAngle * (pos_local_outerSeg.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_outerSeg.y() + deltay_roll))){
                  data_.RdPhi_outerSeg_GE11 = cosAngle * (pos_local_outerSeg.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_outerSeg.y() + deltay_roll);
                }
              }
              if (data_.has_prop_innerSeg){
                if (abs(data_.RdPhi_innerSeg_GE11) > abs(cosAngle * (pos_local_innerSeg.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_innerSeg.y() + deltay_roll))){
                  data_.RdPhi_innerSeg_GE11 = cosAngle * (pos_local_innerSeg.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_innerSeg.y() + deltay_roll);
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
        if(debug)cout << "Starting sim info" << endl;
        data_.nSim = 99999999;
        data_.simDy = 999.;
        float tmpDy = 999.;
        float tmpDr = 999.;
        int tmpSimCounter = 0;
        for (const auto& simHit:*gemSimHits.product()) {
          GEMDetId gemid((simHit).detUnitId());
          if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
            if(debug)cout << "Found a match simhit" << endl;
            tmpSimCounter ++;
            const auto& etaPart = GEMGeometry_->etaPartition(gemid);
            GlobalPoint pGlobal = tsos_CSC.globalPosition();
            float dy = pGlobal.y() - etaPart->toGlobal(simHit.localPosition()).y();
            float dx = pGlobal.x() - etaPart->toGlobal(simHit.localPosition()).x();
            if (dy < tmpDy) tmpDy = dy;
            if (pow(pow(dy, 2) + pow(dx, 2), 0.5) < tmpDr){
              data_.sim_x = etaPart->toGlobal(simHit.localPosition()).x();
              data_.sim_y = etaPart->toGlobal(simHit.localPosition()).y();
              data_.sim_z = etaPart->toGlobal(simHit.localPosition()).z();
              data_.sim_r = pow(pow(data_.sim_x, 2) + pow(data_.sim_y, 2), 0.5);
              data_.sim_localx = simHit.localPosition().x();
              data_.sim_localy = simHit.localPosition().y();
              data_.sim_localy_adjusted = data_.sim_localy + etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
              tmpDr = pow(pow(dy, 2) + pow(dx, 2), 0.5);
            }
          }
        }
        data_.simDy = tmpDy;
        data_.nSim = tmpSimCounter;
      }

      std::cout << "Sim point was " << data_.sim_localx << ", " << data_.sim_localy << std::endl;

      if (debug)std::cout << "Num of rechits = " << rechit_counter << std::endl;
      if (debug)std::cout << "Num of matches = " << rechit_matches << std::endl;
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

void analyser::beginJob(){}
void analyser::endJob(){}

DEFINE_FWK_MODULE(analyser);

#include <assert.h>
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

//user include files below
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectoryFitter.h"
#include "TrackingTools/TrackFitters/interface/KFTrajectorySmoother.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "Alignment/CommonAlignmentAlgorithm/interface/AlignmentAlgorithmBase.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "RecoLocalMuon/CSCSegment/src/CSCSegmentProducer.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;


struct MuonData
{
  void init();
  TTree* book(TTree *t, int prop_type);
  //============ Muon Info ================//
  int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
  unsigned long long evtNum; unsigned long long lumiBlock; int muonIdx;
  int runNum;
  bool has_TightID; bool isPFIsoTightMu;
  //============ Propagation Info =========//
  float prop_GP[3]; float prop_LP[3]; float prop_startingPoint_GP[3];
  float prop_dxdz;  float prop_yroll; float prop_localphi_rad;
  float prop_localphi_deg;            float prop_globalphi_rad;
  bool has_prop;    bool has_fidcut;
  int prop_location[5];
  //============ Track Info ===============//
  float track_chi2; float track_ndof; int n_ME11_segment;
  int which_track;  int hasME11;      int hasME11RecHit;
  int hasME11A;     int hasME11ARecHit;
  int nCSCSeg;      int nDTSeg;       int nME11RecHits;
  float ME11_BunchX;                  int ME11_strip;
  int ME11_location[5];               int inner_or_outer_mom;
  float ME11_Segment_Direction[3];
  float ME11_Segment_slope_dxdz;      float ME11_Segment_slope_dydz;
  int eighthStripDiff;
  int n_ME21_segment; int hasME21;    int hasME21RecHit;
  int nME21RecHits; float ME21_BunchX; int ME21_strip;
  int ME21_location[5];
  float ME21_Segment_Direction[3];
  float ME21_Segment_slope_dxdz;      float ME21_Segment_slope_dydz;
  //============ Rechit Info GE11 =============//
  float rechit_GP[3]; float rechit_LP[3];        bool has_rechit;
  float rechit_yroll; float rechit_localphi_rad; float rechit_localphi_deg;
  int rechit_first_strip;     int rechit_CLS;    int rechit_BunchX;
  float RdPhi;        float RdPhi_Corrected;     int rechit_detId;
  float dPhi;   float dPhi_Corrected;
  float bending_angle;
  int nRecHitsTot;    int nRecHits5;             int nRecHits2;
  int rechit_location[5];
  int nRecHitsRpos1L1; int nRecHitsRpos1L2;
  int nRecHitsRneg1L1; int nRecHitsRneg1L2;
  //============ Rechit Info GE21 =============//
  float rechit_GP_GE21[3]; float rechit_LP_GE21[3];        bool has_rechit_GE21;
  float rechit_yroll_GE21; float rechit_localphi_rad_GE21; float rechit_localphi_deg_GE21;
  int rechit_first_strip_GE21;     int rechit_CLS_GE21;    int rechit_BunchX_GE21;
  float RdPhi_GE21;        float RdPhi_Corrected_GE21;     int rechit_detId_GE21;
  float dPhi_GE21;   float dPhi_Corrected_GE21;
  float bending_angle_GE21;
  int nRecHitsTot_GE21;    int nRecHits5_GE21;             int nRecHits2_GE21;
  int rechit_location_GE21[5];
  int nRecHitsRpos1L1_GE21; int nRecHitsRpos1L2_GE21;
  int nRecHitsRneg1L1_GE21; int nRecHitsRneg1L2_GE21;
  //=========== Sim info for MC ==========//
  float sim_GP[3];   float sim_LP[3];
  float simDy;       float sim_yroll;            int nSim;
};

void MuonData::init()
{
  //=========== Muon Info ===============//
  muon_charge = 9999; muon_pt = 99999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
  has_TightID = 0;    isPFIsoTightMu = 0;
  //=========== Propagation Info =======//
  for(int i=0; i<3; ++i){
    prop_GP[i] = 99999; prop_LP[i] = 99999; prop_startingPoint_GP[i] = 99999;
  }
  prop_dxdz = 99999; prop_yroll = 99999;
  prop_localphi_rad = 99999; prop_localphi_deg = 99999;
  prop_globalphi_rad = 99999;
  has_prop = false; has_fidcut = false;
  for(int i=0 ; i<5; ++i){
    prop_location[i] = 99999;
  }
  //=========== Track Info =============//
  track_chi2 = 999999; track_ndof = 999999; n_ME11_segment = 999999; which_track = 999999;
  hasME11 = 0; hasME11RecHit = 0; hasME11A = 0; hasME11ARecHit = 0;
  nCSCSeg = 999999; nDTSeg = 999999; nME11RecHits = 999999; ME11_BunchX = 999999; ME11_strip = 999999;
  n_ME21_segment = 999999; hasME21 = 0; hasME21RecHit = 0;
  nME21RecHits = 999999; ME21_BunchX = 999999; ME21_strip = 999999;
  for(int i=0; i<5; ++i){
    ME11_location[i] = 999999;
    ME21_location[i] = 999999;
  }
  inner_or_outer_mom = 99999;
  for (int i=0; i<3; ++i){
    ME11_Segment_Direction[i] = 999999;
    ME21_Segment_Direction[i] = 999999;
  }
  ME11_Segment_slope_dxdz = 999999;    ME11_Segment_slope_dydz = 999999;
  ME21_Segment_slope_dxdz = 999999;    ME21_Segment_slope_dydz = 999999;
  eighthStripDiff = 99999;
  //=========== Rechit Info GE11 ===========//
  for (int i=0; i<3; ++i){
    rechit_GP[i] = 999999; rechit_LP[i] = 999999;
  }
  rechit_yroll = 999999; rechit_localphi_rad = 999999; rechit_localphi_deg = 999999;
  has_rechit = false;
  rechit_first_strip = 999999; rechit_CLS = 999999; rechit_BunchX = 999999;
  RdPhi = 999999; RdPhi_Corrected = 999999; rechit_detId = 999999;
  dPhi = 999999; dPhi_Corrected = 999999;
  bending_angle = 999999;
  nRecHitsTot = 999999; nRecHits5 = 999999; nRecHits2 = 999999;
  for (int i=0; i<5; ++i){
    rechit_location[i] = 999999;
  }
  nRecHitsRpos1L1 = 999999; nRecHitsRpos1L2 = 999999;
  nRecHitsRneg1L1 = 999999; nRecHitsRneg1L2 = 999999;
  //=========== Rechit Info GE21 ===========//
  for (int i=0; i<3; ++i){
    rechit_GP_GE21[i] = 999999; rechit_LP_GE21[i] = 999999;
  }
  rechit_yroll_GE21 = 999999; rechit_localphi_rad_GE21 = 999999; rechit_localphi_deg_GE21 = 999999;
  has_rechit_GE21 = false;
  rechit_first_strip_GE21 = 999999; rechit_CLS_GE21 = 999999; rechit_BunchX_GE21 = 999999;
  RdPhi_GE21 = 999999; RdPhi_Corrected_GE21 = 999999; rechit_detId_GE21 = 999999;
  dPhi_GE21 = 999999; dPhi_Corrected_GE21 = 999999;
  bending_angle_GE21 = 999999;
  nRecHitsTot_GE21 = 999999; nRecHits5_GE21 = 999999; nRecHits2_GE21 = 999999;
  for (int i=0; i<5; ++i){
    rechit_location_GE21[i] = 999999;
  }
  nRecHitsRpos1L1_GE21 = 999999; nRecHitsRpos1L2_GE21 = 999999;
  nRecHitsRneg1L1_GE21 = 999999; nRecHitsRneg1L2_GE21 = 999999;

  //Sim info for MC
  for (int i=0; i<3; ++i){
    sim_GP[i] = 9999999; sim_LP[i] = 9999999;
  }
  simDy = 9999999; sim_yroll = 9999999; nSim = 9999999;
}

TTree* MuonData::book(TTree *t, int prop_type){
  edm::Service< TFileService > fs;
  if(prop_type == 1){
    t = fs->make<TTree>("CSC_Prop", "CSC_Prop");
  }
  else if(prop_type == 2){
    t = fs->make<TTree>("Inner_Prop", "Inner_Prop");
  }
  else if(prop_type == 3){
    t = fs->make<TTree>("ME11Seg_Prop", "ME11Seg_Prop");
  }
  else if(prop_type == 4){
    t = fs->make<TTree>("InnerRefit_Prop", "InnerRefit_Prop");
  }
  else if(prop_type == 5){
    t = fs->make<TTree>("ME11SegReco_Prop", "ME11SegReco_Prop");
  }
  else{
    std::cout << "Bad prop type, failure, doesnt fall under the 5 prop_type listed" << std::endl;
  }
  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta);       t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum);           t->Branch("lumiBlock", &lumiBlock);
  t->Branch("runNum", &runNum);           t->Branch("muonIdx", &muonIdx);
  t->Branch("has_TightID", &has_TightID); t->Branch("isPFIsoTightMu", &isPFIsoTightMu);
  //========== Propagation Info =======//
  t->Branch("prop_GP", &prop_GP, "prop_GP[3] (x,y,z)/F");
  t->Branch("prop_LP", &prop_LP, "prop_LP[3] (x,y,z)/F");
  t->Branch("prop_dxdz", &prop_dxdz);
  t->Branch("prop_startingPoint_GP", &prop_startingPoint_GP, "prop_startingPoint_GP[3] (x,y,z)/F");
  t->Branch("prop_yroll", &prop_yroll);
  t->Branch("prop_localphi_rad", &prop_localphi_rad);
  t->Branch("prop_localphi_deg", &prop_localphi_deg);
  t->Branch("prop_globalphi_rad", &prop_globalphi_rad);
  t->Branch("has_prop", &has_prop);
  t->Branch("has_fidcut", &has_fidcut);
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, cha, lay, rol)/I");
  //========== Track Info =============//
  t->Branch("track_chi2", &track_chi2);         t->Branch("track_ndof", &track_ndof);
  t->Branch("n_ME11_segment", &n_ME11_segment); t->Branch("which_track", &which_track);
  t->Branch("hasME11", &hasME11);               t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A);             t->Branch("hasME11ARecHit", &hasME11ARecHit);
  t->Branch("eighthStripDiff", &eighthStripDiff);
  t->Branch("nCSCSeg", &nCSCSeg);               t->Branch("nDTSeg", &nDTSeg);
  t->Branch("nME11RecHits", &nME11RecHits);     t->Branch("ME11_BunchX", &ME11_BunchX);
  t->Branch("ME11_strip", &ME11_strip);
  t->Branch("ME11_location", &ME11_location, "ME11_location[5] (end, sta, ring, cha, lay)/I");
  t->Branch("ME11_Segment_Direction", &ME11_Segment_Direction, "ME11_Segment_Direction[3] (x,y,z)/F");
  t->Branch("ME11_Segment_slope_dxdz", &ME11_Segment_slope_dxdz);
  t->Branch("ME11_Segment_slope_dydz", &ME11_Segment_slope_dydz);
  t->Branch("inner_or_outer_mom", &inner_or_outer_mom, "inner_or_outer_mom (0 = inner, 1 = outer)/I");

  t->Branch("n_ME21_segment", &n_ME21_segment);
  t->Branch("hasME21", &hasME21);               t->Branch("hasME21RecHit", &hasME21RecHit);
  t->Branch("nME21RecHits", &nME21RecHits);     t->Branch("ME21_BunchX", &ME21_BunchX);
  t->Branch("ME21_strip", &ME21_strip);
  t->Branch("ME21_location", &ME21_location, "ME21_location[5] (end, sta, ring, cha, lay)/I");
  t->Branch("ME21_Segment_Direction", &ME21_Segment_Direction, "ME21_Segment_Direction[3] (x,y,z)/F");
  t->Branch("ME21_Segment_slope_dxdz", &ME21_Segment_slope_dxdz);
  t->Branch("ME21_Segment_slope_dydz", &ME21_Segment_slope_dydz);
  //========== Rechit Info GE11 ============//
  t->Branch("rechit_GP", &rechit_GP, "rechit_GP[3] (x,y,z)/F");
  t->Branch("rechit_LP", &rechit_LP, "rechit_LP[3] (x,y,z)/F");
  t->Branch("rechit_yroll", &rechit_yroll);
  t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
  t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
  t->Branch("has_rechit", &has_rechit);
  t->Branch("rechit_first_strip", &rechit_first_strip);
  t->Branch("rechit_CLS", &rechit_CLS);
  t->Branch("rechit_BunchX", &rechit_BunchX);
  t->Branch("RdPhi", &RdPhi);
  t->Branch("RdPhi_Corrected", &RdPhi_Corrected);  
  t->Branch("dPhi", &dPhi);
  t->Branch("dPhi_Corrected", &dPhi_Corrected);
  t->Branch("rechit_detId", &rechit_detId);
  t->Branch("bending_angle", &bending_angle);
  t->Branch("nRecHitsTot", &nRecHitsTot);
  t->Branch("nRecHits2", &nRecHits2);
  t->Branch("nRecHits5", &nRecHits5);
  t->Branch("rechit_location", &rechit_location, "rechit_location[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("nRecHitsRpos1L1", &nRecHitsRpos1L1);
  t->Branch("nRecHitsRpos1L2", &nRecHitsRpos1L2);
  t->Branch("nRecHitsRneg1L1", &nRecHitsRneg1L1);
  t->Branch("nRecHitsRneg1L2", &nRecHitsRneg1L2);
  //========== Rechit Info GE21 ============//
  t->Branch("rechit_GP_GE21", &rechit_GP_GE21, "rechit_GP_GE21[3] (x,y,z)/F");
  t->Branch("rechit_LP_GE21", &rechit_LP_GE21, "rechit_LP_GE21[3] (x,y,z)/F");
  t->Branch("rechit_yroll_GE21", &rechit_yroll_GE21);
  t->Branch("rechit_localphi_rad_GE21", &rechit_localphi_rad_GE21);
  t->Branch("rechit_localphi_deg_GE21", &rechit_localphi_deg_GE21);
  t->Branch("has_rechit_GE21", &has_rechit_GE21);
  t->Branch("rechit_first_strip_GE21", &rechit_first_strip_GE21);
  t->Branch("rechit_CLS_GE21", &rechit_CLS_GE21);
  t->Branch("rechit_BunchX_GE21", &rechit_BunchX_GE21);
  t->Branch("RdPhi_GE21", &RdPhi_GE21);
  t->Branch("RdPhi_Corrected_GE21", &RdPhi_Corrected_GE21);
  t->Branch("dPhi_GE21", &dPhi_GE21);
  t->Branch("dPhi_Corrected_GE21", &dPhi_Corrected_GE21);
  t->Branch("rechit_detId_GE21", &rechit_detId_GE21);
  t->Branch("bending_angle_GE21", &bending_angle_GE21);
  t->Branch("nRecHitsTot_GE21", &nRecHitsTot_GE21);
  t->Branch("nRecHits2_GE21", &nRecHits2_GE21);
  t->Branch("nRecHits5_GE21", &nRecHits5_GE21);
  t->Branch("rechit_location_GE21", &rechit_location_GE21, "rechit_location_GE21[5] (reg, sta, cha, lay, rol)/I");
  t->Branch("nRecHitsRpos1L1_GE21", &nRecHitsRpos1L1_GE21);
  t->Branch("nRecHitsRpos1L2_GE21", &nRecHitsRpos1L2_GE21);
  t->Branch("nRecHitsRneg1L1_GE21", &nRecHitsRneg1L1_GE21);
  t->Branch("nRecHitsRneg1L2_GE21", &nRecHitsRneg1L2_GE21);
  //========== Sim Info ==============//
  t->Branch("sim_GP", &sim_GP, "sim_GP[3] (x,y,z)/F");
  t->Branch("sim_LP", &sim_LP, "sim_LP[3] (x,y,z)/F");
  t->Branch("simDy", &simDy);
  t->Branch("sim_yroll", &sim_yroll);
  t->Branch("nSim", &nSim);
  return t;
}

class analyzer : public edm::one::EDAnalyzer<> {
public:
  explicit analyzer(const edm::ParameterSet&);
  ~analyzer(){};

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  void propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i, const Trajectory* traj_of_muon);
  void CSCSegmentCounter(const reco::Muon* mu, MuonData& data_, int prop_type, std::vector<const CSCSegment*> &RecoSegment);//const CSCSegment* RecoSegment);//const CSCSegmentCollection::const_iterator* RecoSegment);
  void propagate_to_GEM(const reco::Muon* mu, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_, const Trajectory* traj_of_muon, std::vector<const CSCSegment*> RecoSegment);//const CSCSegment* RecoSegment);//const CSCSegmentCollection::const_iterator* RecoSegment);
  void GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_, int prop_type, std::vector<const CSCSegment*> RecoSegment);//const CSCSegment* RecoSegment);//const CSCSegmentCollection::const_iterator* RecoSegment);
  void GEM_simhit_matcher(const GEMEtaPartition* ch, GlobalPoint prop_GP, MuonData& data_);
  float RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch);
  bool fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch);

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::Handle<vector<PSimHit> > gemSimHits;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  //edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegmentsReco_;
  edm::Handle<CSCSegmentCollection> cscSegmentsReco;
  edm::Handle<TrajTrackAssociationCollection> ref_track;
  edm::EDGetTokenT<TrajTrackAssociationCollection> ref_track_;

  edm::Service<TFileService> fs;
  MuonServiceProxy* theService_;
  
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  
  edm::ESHandle<GEMGeometry> GEMGeometry_;
  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool CSC_prop; bool tracker_prop; bool Segment_prop; bool trackerRefit_prop; bool SegmentReco_prop;
  vector<int> prop_list;
  bool debug;
  bool isCosmic;

  MuonData data_;
  TTree* CSC_tree; TTree* Tracker_tree; TTree* Segment_tree; TTree* TrackerRefit_tree; TTree* SegmentReco_tree;
  TH2D* nME11_col_vs_matches = new TH2D("nME11_test", "nME11_test", 5, 0, 5, 5, 0, 5);

  bool isMC;
  const CSCSegment *ME11_segment;

  const edm::ESGetToken<GEMGeometry, MuonGeometryRecord> gemGeomToken_;
  const edm::ESGetToken<CSCGeometry, MuonGeometryRecord> cscGeomToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geomToken_;
};


analyzer::analyzer(const edm::ParameterSet& iConfig)
  : gemGeomToken_(esConsumes()),
    cscGeomToken_(esConsumes()),
    ttkToken_(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))),
    geomToken_(esConsumes())
{
  cout << "Begin analyzer" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  vertexCollection_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection")); 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gemSimHits_ = consumes<vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("gemSimHits"));
  //cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));
  cscSegmentsReco_ = consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegmentsReco"));
  ref_track_ = consumes<TrajTrackAssociationCollection>(iConfig.getParameter<InputTag>("ref_track"));

  tracker_prop = iConfig.getParameter<bool>("tracker_prop");
  CSC_prop = iConfig.getParameter<bool>("CSC_prop");
  Segment_prop = iConfig.getParameter<bool>("Segment_prop");
  trackerRefit_prop = iConfig.getParameter<bool>("trackerRefit_prop");
  SegmentReco_prop = iConfig.getParameter<bool>("SegmentReco_prop");
  debug = iConfig.getParameter<bool>("debug");
  isCosmic = iConfig.getParameter<bool>("isCosmic");
  std::cout << "tracker_prop:trackerRefit_prop " << tracker_prop << ":" << trackerRefit_prop << "\tCSC_prop " << CSC_prop << "\tSegment_prop:SegmentReco_prop " << Segment_prop << ":" << SegmentReco_prop << "\tdebug " << debug << std::endl;

  if(CSC_prop){CSC_tree = data_.book(CSC_tree, 1); prop_list.push_back(1);}
  if(tracker_prop){Tracker_tree = data_.book(Tracker_tree, 2); prop_list.push_back(2);}
  if(Segment_prop){Segment_tree = data_.book(Segment_tree, 3); prop_list.push_back(3);}
  if(trackerRefit_prop){TrackerRefit_tree = data_.book(TrackerRefit_tree, 4); prop_list.push_back(4);}
  if(SegmentReco_prop){SegmentReco_tree = data_.book(SegmentReco_tree, 5); prop_list.push_back(5);}
}


void
analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  GEMGeometry_ = &iSetup.getData(gemGeomToken_);
  CSCGeometry_ = &iSetup.getData(cscGeomToken_);
  ttrackBuilder_ = &iSetup.getData(ttkToken_);
  theTrackingGeometry = &iSetup.getData(geomToken_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  isMC = false;
  if (! iEvent.eventAuxiliary().isRealData()) isMC = true;
  iEvent.getByToken(gemRecHits_, gemRecHits);
  if (isMC) {
    iEvent.getByToken(gemSimHits_, gemSimHits);
  }

  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;

  //adding refit trajectory

  edm::Handle<TrajTrackAssociationCollection> ref_track;
  iEvent.getByToken(ref_track_, ref_track);

  ConstTrajTrackPairs ref_track_pairs;
  for (auto it = ref_track->begin(); it != ref_track->end(); ++it) {
    ref_track_pairs.push_back(ConstTrajTrackPair(&(*(*it).key), &(*(*it).val)));
  } //the loop goes over tracks and saves the key and value of each track as a pair in ref_track_pairs.
  //end of Refit trajectory

  //edm::Handle<CSCSegmentCollection> cscSegments;
  //if (! iEvent.getByToken(cscSegments_, cscSegments)){std::cout << "Bad segments" << std::endl;}

  iEvent.getByToken(cscSegmentsReco_, cscSegmentsReco);
  if (debug) std::cout << "cscSegmentsReco->size() " << cscSegmentsReco->size() << std::endl;
  if (! iEvent.getByToken(cscSegmentsReco_, cscSegmentsReco)){std::cout << "Bad segments" << std::endl;}
  
  if (debug) cout << "New! EventNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    if (not mu->isGlobalMuon()) continue;
    if (debug) cout << "new muon, i = " << i << endl;
    
    //trajectory muon matching--we want to match muons with the tracks
    const Trajectory* traj_of_muon;
    const Trajectory* traj_of_Track;
    const reco::Track* track_of_Track;

    if (!(mu->passed(reco::Muon::PFIsoTight))) continue;
    for (auto it = std::begin(prop_list); it != std::end(prop_list); ++it){
      if (debug) std::cout << "\tprop " << *it << "about to start propagate" << std::endl;
      int prop_type = *it;
      for (ConstTrajTrackPairs::const_iterator iter = ref_track_pairs.begin(); iter != ref_track_pairs.end(); ++iter) {
        traj_of_Track = (*iter).first;
        track_of_Track = (*iter).second;
        if (track_of_Track == mu->track().get()) {
          traj_of_muon = traj_of_Track;
          if (debug) cout << "mu, prop_type, i, cscSegmentsReco->size(): " << mu << "," << prop_type << "," << i << "," << cscSegmentsReco->size() << endl;
          propagate(mu, prop_type, iEvent, i, traj_of_muon); //taking the matched trajectory of muon to the propagate function
        }
      }
    }
  }
}

void analyzer::propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i, const Trajectory* traj_of_muon){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  //if (debug) cout << "\tGetting tree, Track, ttTrack " << endl;

  //===============Vertex selection
  edm::Handle<reco::VertexCollection> vertexCollection;
  iEvent.getByToken(vertexCollection_, vertexCollection);
  reco::Vertex vertexSelection; //choose type of vertex needed
  for (const auto& vertex : *vertexCollection.product()){
    if (vertexCollection.isValid()) {
      vertexSelection = vertex;
      break; //selecting the first valid vertex
    }
  }
  //================end of Vertex selection
  if (prop_type == 1){
    tree = CSC_tree;
    if (!(mu->isGlobalMuon())) {return;} //if(!(mu->isStandAloneMuon())){return;}
    if (!(mu->globalTrack().isNonnull())) {return;}  //if(!(mu->outerTrack().isNonnull())){return;}
    Track = mu->globalTrack().get(); //Track = mu->outerTrack().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 2 or prop_type == 4){
    if (prop_type == 2) {tree = Tracker_tree;}
    else{tree = TrackerRefit_tree;} 
    if (!(mu->isTrackerMuon())) {return;}
    if (!(mu->track().isNonnull())) {return;}
    Track = mu->track().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 3 or prop_type == 5){
    if (prop_type == 3) {tree = Segment_tree;}
    else if (prop_type == 5) {tree = SegmentReco_tree;}
    if(isCosmic){
      if(!(mu->isStandAloneMuon())){return;}
      if(!(mu->outerTrack().isNonnull())){return;}
      Track = mu->outerTrack().get();
      ttTrack = ttrackBuilder_->build(Track);
    }
    else{
      if (!(mu->isGlobalMuon())) {return;}//if (!(mu->isTrackerMuon())) {return;}
      if (!(mu->globalTrack().isNonnull())) {return;}//if (!(mu->track().isNonnull())) {return;}
      Track = mu->globalTrack().get();//Track = mu->track().get();
      ttTrack = ttrackBuilder_->build(Track);
    }
  }
  else{
    if (debug) cout << "Bad prop type, failure, not one of the 5" << endl;
    return;
  }

  if (!(ttTrack.isValid())) {cout << "Bad event, no track" << endl;}
  if (debug) cout << "Got track, now initiating data" << endl;
  data_.init();
  //======================Muon Info=======================
  data_.muon_charge = mu->charge(); data_.muon_pt = mu->pt();  data_.muon_eta = mu->eta();
  data_.muon_momentum = mu->momentum().mag2();                 data_.evtNum = iEvent.eventAuxiliary().event();
  data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock(); data_.muonIdx = data_.evtNum*100 + i;
  data_.runNum = iEvent.run();
  data_.has_TightID = muon::isTightMuon(*mu, vertexSelection);
  //=====================Track Info=======================
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  //const CSCSegment* RecoSegment; //const CSCSegmentCollection::const_iterator* RecoSegment;
  //CSCSegment* RecoSegment;
  std::vector<const CSCSegment*> RecoSegment;
  CSCSegmentCounter(mu, data_, prop_type, RecoSegment);
  if (prop_type == 3 or prop_type == 5) {
    if (data_.ME11_location[1]==1){
      if (data_.hasME11 != 1) {
        cout << "data_.ME11_location[1]==1 and data_.hasME11 == 0 " << data_.ME11_location[1] << " " << data_.hasME11 << endl;
        return;
      }
    }
    else if (data_.ME21_location[1]==2){
      if (data_.hasME21 != 1) {
        cout << "data_.ME21_location[1]==1 and data_.hasME21 == 0 " << data_.ME21_location[1] << " " << data_.hasME21 << endl;
        return;
      }
    }
    else {
      cout << "data_.ME11_location[1] and data_.ME21_location[1] " << data_.ME11_location[1] << " " << data_.ME21_location[1] << endl;
      return;
    }
  }
  //================Propagation Info===================
  if (debug) cout << "starting chamber loop" << endl;
  for (const auto& ch : GEMGeometry_->etaPartitions()) {
    if (ch->id().station() > 2) continue; //only concerned about GE1/1 and GE2/1
    GlobalPoint tmp_prop_GP;        bool tmp_has_prop = 0;
    propagate_to_GEM(mu, ch, prop_type, tmp_has_prop, tmp_prop_GP, data_, traj_of_muon, RecoSegment);
    if (tmp_has_prop){
      LocalPoint tmp_prop_LP = ch->toLocal(tmp_prop_GP);
      //==============RecHit Info======================
      if (debug) cout << "before rechit_matcher" << endl;
      GEM_rechit_matcher(ch, tmp_prop_LP, data_, prop_type, RecoSegment);
      if (isMC){
        GEM_simhit_matcher(ch, tmp_prop_GP, data_);
      }
      tree->Fill();
    }
  }
}

void analyzer::CSCSegmentCounter(const reco::Muon* mu, MuonData& data_, int prop_type, std::vector<const CSCSegment*> &RecoSegment){//const CSCSegment* RecoSegment){//const CSCSegmentCollection::const_iterator* RecoSegment){
  if (!(mu->isGlobalMuon())) {return;} //!(mu->isStandAloneMuon())
  if (!(mu->globalTrack().isNonnull())) {return;} //!(mu->outerTrack().isNonnull())
  const reco::Track* Track = mu->globalTrack().get();

  int tmp_CSC_counter = 0;   int tmp_DT_counter = 0;   int tmp_ME11_counter = 0;
  int tmp_ME11RecHit_counter = 0; float tmp_ME11_BunchX = 99999;
  int tmp_ME11_strip = 99999; bool tmp_hasME11A = 0;
  float tmp_me11_segment_x; float tmp_me11_segment_y; float tmp_me11_segment_z;
  float tmp_me11_segment_slope_dxdz; float tmp_me11_segment_slope_dydz;

  //ME2/1
  int tmp_ME21_counter = 0; int tmp_ME21RecHit_counter = 0; float tmp_ME21_BunchX = 99999;
  int tmp_ME21_strip = 99999;
  float tmp_me21_segment_x; float tmp_me21_segment_y; float tmp_me21_segment_z;
  float tmp_me21_segment_slope_dxdz; float tmp_me21_segment_slope_dydz;
  
  //Below part is for Cosmic Muon Study
  if(isCosmic){
    tmp_CSC_counter = mu->numberOfSegments(1,2) + mu->numberOfSegments(2,2) + mu->numberOfSegments(3,2) + mu->numberOfSegments(4,2);
    tmp_DT_counter = mu->numberOfSegments(1,1) + mu->numberOfSegments(2,1) + mu->numberOfSegments(3,1) + mu->numberOfSegments(4,1);
    auto matches = mu->matches();
    for (auto MCM : matches){
      if(MCM.detector() != 2) continue;
      for(auto MSM : MCM.segmentMatches){
        auto cscSegRef = MSM.cscSegmentRef;
        auto cscDetID = cscSegRef->cscDetId();
        if(cscDetID.station() == 1 and (cscDetID.ring() == 1 or cscDetID.ring() == 4)){
          if(cscDetID.ring() == 4){tmp_hasME11A = 1;}
          tmp_ME11_counter++;
          if (debug){std::cout << "isCosmic = True! Getting ME11 Segment" << std::endl;}
          ME11_segment = cscSegRef.get();
          tmp_ME11RecHit_counter = (cscSegRef.get())->nRecHits(); // Find the real function for this. Bad if multiple segments.
          tmp_ME11_BunchX = ME11_segment->time();
          auto cscDetID_FAKE = CSCDetId(cscDetID.endcap(), cscDetID.station(), cscDetID.ring(), cscDetID.chamber(), 3);
          const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
          const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
          tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());
          data_.ME11_location[0] = cscDetID.endcap(); data_.ME11_location[1] = cscDetID.station(); data_.ME11_location[2] = cscDetID.ring(); data_.ME11_location[3] = cscDetID.chamber(); data_.ME11_location[4] = cscDetID.layer();
        }
      }
    }
  }
  //end of Cosmics
  
  else{
    //if (debug) cout << "Track->recHitsSize(): " << Track->recHitsSize() << endl;
    for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
      const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
      DetId RecHitId = RecHit->geographicalId(); 
      uint16_t RecHitDetId = RecHitId.det();

      if (RecHitDetId == DetId::Muon){
        uint16_t RecHitSubDet = RecHitId.subdetId();
        if (RecHitSubDet == (uint16_t)MuonSubdetId::CSC){
          if (CSCDetId(RecHitId).station() == 1 and (CSCDetId(RecHitId).ring() == 1 or CSCDetId(RecHitId).ring() == 4) and RecHit->dimension() == 4){
            tmp_ME11_counter++;
            if (CSCDetId(RecHitId).ring() == 4) {tmp_hasME11A = 1;}
            if (debug) cout << "tmp_hasME11A: " << tmp_hasME11A << endl;

            RecSegment* Rec_segment = (RecSegment*)RecHit;

            if (prop_type==5){
              float tmp_mutrack_x = RecHit->localPosition().x(); float tmp_mutrack_y = RecHit->localPosition().y(); float tmp_mutrack_z = RecHit->localPosition().z();
       	      CSCSegmentCollection::const_iterator RecoSeg;
              for (RecoSeg = cscSegmentsReco->begin(); RecoSeg !=cscSegmentsReco->end(); RecoSeg++){
                LocalVector ME11_RecoSeg = RecoSeg->localDirection();
                DetId RecoSegId = RecoSeg->geographicalId(); uint16_t RecoSegDetId = RecoSegId.det();
                float tmp_min_distance_seg = 99999;
		if (RecoSegDetId==(uint16_t)MuonSubdetId::CSC){
                  if (CSCDetId(RecHitId).endcap()==CSCDetId(RecoSegId).endcap() and
                      CSCDetId(RecHitId).station()==CSCDetId(RecoSegId).station() and
                      CSCDetId(RecHitId).ring()==CSCDetId(RecoSegId).ring() and
                      CSCDetId(RecHitId).chamber()==CSCDetId(RecoSegId).chamber() and
		      CSCDetId(RecHitId).layer()==CSCDetId(RecoSegId).layer() ){
                      
		    float reco_segment_x = ME11_RecoSeg.x(); float reco_segment_y = ME11_RecoSeg.y(); float reco_segment_z = ME11_RecoSeg.z();
		    float tmp_distance = sqrt((tmp_mutrack_x-reco_segment_x)*(tmp_mutrack_x-reco_segment_x)
                                              +(tmp_mutrack_y-reco_segment_y)*(tmp_mutrack_y-reco_segment_y)
                                              +(tmp_mutrack_z-reco_segment_z)*(tmp_mutrack_z-reco_segment_z));
                    if (tmp_distance < tmp_min_distance_seg){
		      //const CSCSegmentCollection::const_iterator* RecoSegment = RecoSeg->clone();
                      if (debug) cout << "RecoSeg: " << *RecoSeg << endl;
                      //const CSCSegment* RecoSegment = (CSCSegment*)RecoSeg;
                      //CSCSegment* tmp_segment = *RecoSeg;
                      //RecoSegment = (CSCSegment*)tmp_segment;
                      //const CSCSegment* RecoSegment = ((CSCSegment*)RecoSeg)->clone();
                      //CSCSegment* RecoSegment = *RecoSeg;
                      //const CSCSegment* RecoSegment = ((CSCSegment*)RecoSeg).get();
                      RecoSegment.push_back(&(*RecoSeg));
                      //if (debug) cout << "RecoSegment: " <<RecoSegment[0]<< endl;

                      tmp_me11_segment_x = reco_segment_x;
		      tmp_me11_segment_y = reco_segment_y;
		      tmp_me11_segment_z = reco_segment_z;
                      tmp_me11_segment_slope_dxdz = tmp_me11_segment_x / tmp_me11_segment_z;
                      tmp_me11_segment_slope_dydz = tmp_me11_segment_y / tmp_me11_segment_z;

                      tmp_ME11_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
                      auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
                      const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
		      const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
		      tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(RecoSegment[0]->localPosition());

                      if (debug) cout << "test reco seg x:y:z "<< tmp_me11_segment_x << ":" << tmp_me11_segment_y << ":" << tmp_me11_segment_z << endl;
                    }
		  }
                }
		/*
                if (debug){
                  cout << "mu endcap:station:ring:chamber:layer " << CSCDetId(RecHitId).endcap() << ":" << CSCDetId(RecHitId).station() << ":" << CSCDetId(RecHitId).ring() << ":" << CSCDetId(RecHitId).chamber() << ":" << CSCDetId(RecHitId).layer() << endl;
                  if (RecoSegDetId==(uint16_t)MuonSubdetId::CSC) cout << "RecoSeg endcap:station:ring:chamber:layer" << CSCDetId(RecoSegId).endcap() << ":" << CSCDetId(RecoSegId).station() << ":" << CSCDetId(RecoSegId).ring() << ":" << CSCDetId(RecoSegId).chamber() << ":" << CSCDetId(RecoSegId).layer() << endl;
                }
		*/
       	      } //end of the csc segment collection loop
	      //if (debug) cout << "new reco seg x:y:z "<< tmp_me11_segment_x << ":" << tmp_me11_segment_y << ":" << tmp_me11_segment_z << endl;

            }

	    if (prop_type!=5){
              ME11_segment = (CSCSegment*)Rec_segment;
              //if (debug) cout<< "ME11_segment type after type casting: "<< typeid(ME11_segment).name() << endl;
              tmp_me11_segment_x = ME11_segment->localDirection().x();
              tmp_me11_segment_y = ME11_segment->localDirection().y();
              tmp_me11_segment_z = ME11_segment->localDirection().z();
              //if (debug) cout << "ME11segment direction x:y:z" << tmp_me11_segment_x << ":" << tmp_me11_segment_y << ":" << tmp_me11_segment_z << endl;
	  
              tmp_me11_segment_slope_dxdz = tmp_me11_segment_x / tmp_me11_segment_z;
              tmp_me11_segment_slope_dydz = tmp_me11_segment_y / tmp_me11_segment_z;
              tmp_ME11_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
              auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
              const CSCLayer* tmp_ME11_layer = CSCGeometry_->layer(cscDetID_FAKE);
              const CSCLayerGeometry* tmp_ME11_layer_geo = tmp_ME11_layer->geometry();
              tmp_ME11_strip = tmp_ME11_layer_geo->nearestStrip(ME11_segment->localPosition());

            }
	    
	    data_.ME11_Segment_Direction[0] = tmp_me11_segment_x;   data_.ME11_Segment_Direction[1] = tmp_me11_segment_y; data_.ME11_Segment_Direction[2] = tmp_me11_segment_z;
            //if (debug) cout << "CSC Endcap:Station:Ring:SC:Layer " << CSCDetId(RecHitId).endcap() << ":" << CSCDetId(RecHitId).station() << ":" << CSCDetId(RecHitId).ring() << ":" << CSCDetId(RecHitId).chamber() << ":" << CSCDetId(RecHitId).layer() << endl;
            if (debug) cout << "ME11 segment direction x:y:z " << data_.ME11_Segment_Direction[0] << ":" << data_.ME11_Segment_Direction[1] << ":" << data_.ME11_Segment_Direction[2] << endl;
            data_.ME11_Segment_slope_dxdz = tmp_me11_segment_slope_dxdz;  data_.ME11_Segment_slope_dydz= tmp_me11_segment_slope_dydz;
            //if (debug) cout << "segment slope dx/dz:dy/dz " << data_.ME11_Segment_slope_dxdz << ":" << data_.ME11_Segment_slope_dydz << endl;
            data_.ME11_location[0] = CSCDetId(RecHitId).endcap();
            data_.ME11_location[1] = CSCDetId(RecHitId).station();
            data_.ME11_location[2] = CSCDetId(RecHitId).ring();
            data_.ME11_location[3] = CSCDetId(RecHitId).chamber();
            data_.ME11_location[4] = CSCDetId(RecHitId).layer();
          }
          //ME2/1 below
          if(CSCDetId(RecHitId).station() == 2 and CSCDetId(RecHitId).ring() == 1 and RecHit->dimension() == 4){
            tmp_ME21_counter++;
	    RecSegment* Rec_segment = (RecSegment*)RecHit;

	    if (prop_type==5){
	      float tmp_mutrack_x = RecHit->localPosition().x(); float tmp_mutrack_y = RecHit->localPosition().y(); float tmp_mutrack_z = RecHit->localPosition().z();
       	      CSCSegmentCollection::const_iterator RecoSeg;
	      for (RecoSeg = cscSegmentsReco->begin(); RecoSeg !=cscSegmentsReco->end(); RecoSeg++){
                LocalVector ME21_RecoSeg = RecoSeg->localDirection();
                DetId RecoSegId = RecoSeg->geographicalId(); uint16_t RecoSegDetId = RecoSegId.det();
                float tmp_min_distance_seg = 99999;
		if (RecoSegDetId==(uint16_t)MuonSubdetId::CSC){
                  if (CSCDetId(RecHitId).endcap()==CSCDetId(RecoSegId).endcap() and
                      CSCDetId(RecHitId).station()==CSCDetId(RecoSegId).station() and
                      CSCDetId(RecHitId).ring()==CSCDetId(RecoSegId).ring() and
                      CSCDetId(RecHitId).chamber()==CSCDetId(RecoSegId).chamber() and
		      CSCDetId(RecHitId).layer()==CSCDetId(RecoSegId).layer() ){
                      
		    float reco_segment_x = ME21_RecoSeg.x(); float reco_segment_y = ME21_RecoSeg.y(); float reco_segment_z = ME21_RecoSeg.z();
		    float tmp_distance = sqrt((tmp_mutrack_x-reco_segment_x)*(tmp_mutrack_x-reco_segment_x)
                                              +(tmp_mutrack_y-reco_segment_y)*(tmp_mutrack_y-reco_segment_y)
                                              +(tmp_mutrack_z-reco_segment_z)*(tmp_mutrack_z-reco_segment_z));
                    if (tmp_distance < tmp_min_distance_seg){
                                            RecoSegment.push_back(&(*RecoSeg));
		      //if (debug) cout << "RecoSegment: " <<RecoSegment[0]<< endl;

		      tmp_me21_segment_x = reco_segment_x;
		      tmp_me21_segment_y = reco_segment_y;
		      tmp_me21_segment_z = reco_segment_z;
                      tmp_me21_segment_slope_dxdz = tmp_me21_segment_x / tmp_me21_segment_z;
                      tmp_me21_segment_slope_dydz = tmp_me21_segment_y / tmp_me21_segment_z;

                      tmp_ME21_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
                      auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
                      const CSCLayer* tmp_ME21_layer = CSCGeometry_->layer(cscDetID_FAKE);
		      const CSCLayerGeometry* tmp_ME21_layer_geo = tmp_ME21_layer->geometry();
		      tmp_ME21_strip = tmp_ME21_layer_geo->nearestStrip(RecoSegment[0]->localPosition());

                      if (debug) cout << "test reco seg x:y:z "<< tmp_me21_segment_x << ":" << tmp_me21_segment_y << ":" << tmp_me21_segment_z << endl;
                    }
		  }
                }
		/*
                if (debug){
                  cout << "mu endcap:station:ring:chamber:layer " << CSCDetId(RecHitId).endcap() << ":" << CSCDetId(RecHitId).station() << ":" << CSCDetId(RecHitId).ring() << ":" << CSCDetId(RecHitId).chamber() << ":" << CSCDetId(RecHitId).layer() << endl;
                  if (RecoSegDetId==(uint16_t)MuonSubdetId::CSC) cout << "RecoSeg endcap:station:ring:chamber:layer" << CSCDetId(RecoSegId).endcap() << ":" << CSCDetId(RecoSegId).station() << ":" << CSCDetId(RecoSegId).ring() << ":" << CSCDetId(RecoSegId).chamber() << ":" << CSCDetId(RecoSegId).layer() << endl;
                }
		*/
       	      } //end of the csc segment collection loop
	      //if (debug) cout << "new reco seg x:y:z "<< tmp_me21_segment_x << ":" << tmp_me21_segment_y << ":" << tmp_me21_segment_z << endl;

            }

	    if (prop_type!=5){
              ME11_segment = (CSCSegment*)Rec_segment;
              //if (debug) cout<< "ME21_segment type after type casting: "<< typeid(ME11_segment).name() << endl;
              tmp_me21_segment_x = ME11_segment->localDirection().x();
              tmp_me21_segment_y = ME11_segment->localDirection().y();
              tmp_me21_segment_z = ME11_segment->localDirection().z();
              //if (debug) cout << "ME21segment direction x:y:z" << tmp_me21_segment_x << ":" << tmp_me21_segment_y << ":" << tmp_me21_segment_z << endl;
	  
              tmp_me21_segment_slope_dxdz = tmp_me21_segment_x / tmp_me21_segment_z;
              tmp_me21_segment_slope_dydz = tmp_me21_segment_y / tmp_me21_segment_z;
              tmp_ME21_BunchX = ((CSCRecHit2D*)RecHit)->wgroupsBX();
              auto cscDetID_FAKE = CSCDetId(CSCDetId(RecHitId).endcap(), CSCDetId(RecHitId).station(), CSCDetId(RecHitId).ring(), CSCDetId(RecHitId).chamber(), 3);
              const CSCLayer* tmp_ME21_layer = CSCGeometry_->layer(cscDetID_FAKE);
              const CSCLayerGeometry* tmp_ME21_layer_geo = tmp_ME21_layer->geometry();
              tmp_ME21_strip = tmp_ME21_layer_geo->nearestStrip(ME11_segment->localPosition());

            }
	    
	    data_.ME21_Segment_Direction[0] = tmp_me21_segment_x;   data_.ME21_Segment_Direction[1] = tmp_me21_segment_y; data_.ME21_Segment_Direction[2] = tmp_me21_segment_z;
            //if (debug) cout << "CSC Endcap:Station:Ring:SC:Layer " << CSCDetId(RecHitId).endcap() << ":" << CSCDetId(RecHitId).station() << ":" << CSCDetId(RecHitId).ring() << ":" << CSCDetId(RecHitId).chamber() << ":" << CSCDetId(RecHitId).layer() << endl;
            if (debug) cout << "ME21 segment direction x:y:z " << data_.ME21_Segment_Direction[0] << ":" << data_.ME21_Segment_Direction[1] << ":" << data_.ME21_Segment_Direction[2] << endl;
            data_.ME21_Segment_slope_dxdz = tmp_me21_segment_slope_dxdz;  data_.ME21_Segment_slope_dydz= tmp_me21_segment_slope_dydz;
            //if (debug) cout << "segment slope dx/dz:dy/dz " << data_.ME21_Segment_slope_dxdz << ":" << data_.ME21_Segment_slope_dydz << endl;
            data_.ME21_location[0] = CSCDetId(RecHitId).endcap();
            data_.ME21_location[1] = CSCDetId(RecHitId).station();
            data_.ME21_location[2] = CSCDetId(RecHitId).ring();
            data_.ME21_location[3] = CSCDetId(RecHitId).chamber();
            data_.ME21_location[4] = CSCDetId(RecHitId).layer();
          }
	  
          if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1){tmp_ME11RecHit_counter++;}
	  if (CSCDetId(RecHitId).station() == 2 and CSCDetId(RecHitId).ring() == 1){tmp_ME21RecHit_counter++;}
          if (RecHit->dimension() == 4) {tmp_CSC_counter++;}
          //if (debug) cout << "tmp_CSC_counter: " << tmp_CSC_counter << endl;
        }
        if (RecHitSubDet == (uint16_t)MuonSubdetId::DT){
          if (RecHit->dimension() > 1) {tmp_DT_counter++;}
        }
      }
    }
  }
  data_.nCSCSeg = tmp_CSC_counter; data_.nDTSeg = tmp_DT_counter;
  data_.n_ME11_segment = tmp_ME11_counter;
  data_.nME11RecHits = tmp_ME11RecHit_counter;
  data_.ME11_BunchX = tmp_ME11_BunchX;
  data_.ME11_strip = tmp_ME11_strip;
  data_.hasME11A = tmp_hasME11A;
  if (data_.n_ME11_segment >=1 and data_.n_ME11_segment < 1000) {data_.hasME11 = 1;}
  if (debug) cout << "data_.hasME11: " << data_.hasME11 << "\tdata_.n_ME11_segment: " << data_.n_ME11_segment << endl;

  data_.n_ME21_segment = tmp_ME21_counter;
  data_.nME21RecHits = tmp_ME21RecHit_counter;
  data_.ME21_BunchX = tmp_ME21_BunchX;
  data_.ME21_strip = tmp_ME21_strip;
  if (data_.n_ME21_segment >=1 and data_.n_ME21_segment < 1000) {data_.hasME21 = 1;}
  if (debug) cout << "data_.hasME21: " << data_.hasME21 << "\tdata_.n_ME21_segment: " << data_.n_ME21_segment << endl;
}

void analyzer::propagate_to_GEM(const reco::Muon* mu, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_, const Trajectory* traj_of_muon, std::vector<const CSCSegment*> RecoSegment){//const CSCSegment* RecoSegment){
  //if (debug) cout << "propagate_to_gem ch: station,chamber, layer,roll,region" << ch->id().station() << ch->id().chamber() << ch->id().layer() << ch->id().roll() << ch->id().region() << endl;
  const reco::Track* Track;
  reco::TransientTrack ttrack;
  tmp_has_prop = false;
  int tmp_inner_or_outer_mom = 99999;
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  TrajectoryStateOnSurface tsos;
  TrajectoryStateOnSurface tsos_ch; TrajectoryStateOnSurface tsos_seg;
  GlobalPoint pos_startingPoint_GP;
  float prop_dxdz = 99999;

  TrajectoryStateOnSurface previous_trackTSOS;
  double previous_trackTSOS_globalPositionR = 0.0; 
  std::vector<TrajectoryMeasurement> traj_measurement = traj_of_muon->measurements();

  if (prop_type==1 or prop_type==2 or prop_type==4){
    if (prop_type==1){
      if (!(mu->isGlobalMuon())) return;
      Track = mu->globalTrack().get();
      ttrack = ttrackBuilder_->build(Track);
    }
    if (prop_type==2 or prop_type==4){
      if (!(mu->isTrackerMuon())) return;
      Track = mu->track().get();
      ttrack = ttrackBuilder_->build(Track);
    } 
    if (prop_type==4){
      for (std::vector<TrajectoryMeasurement>::const_iterator it_traj_measurement = traj_measurement.begin(); it_traj_measurement != traj_measurement.end(); ++it_traj_measurement){
        TrajectoryMeasurement iTraj_measurement = *it_traj_measurement;
        TrajectoryStateOnSurface tsos = TrajectoryStateCombiner().combine(iTraj_measurement.forwardPredictedState(), iTraj_measurement.backwardPredictedState());
        if (tsos.isValid()){
          double tsosGlobalPositionR = sqrt(tsos.globalPosition().x() * tsos.globalPosition().x() +
                                            tsos.globalPosition().y() * tsos.globalPosition().y());
       
          if (tsosGlobalPositionR > previous_trackTSOS_globalPositionR){
            previous_trackTSOS = tsos;
            previous_trackTSOS_globalPositionR = tsosGlobalPositionR;
          }          
        }
      }

      tsos_ch = propagator->propagate(previous_trackTSOS, ch->surface());

      if (tsos_ch.isValid()){
        const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
        const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
        if (!(tsos_ch.globalPosition().z() * previous_trackTSOS.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
          tmp_has_prop = true;
          pos_GP = tsos_ch.globalPosition();
          pos_startingPoint_GP = previous_trackTSOS.globalPosition();
        }
      }
    }

    if (prop_type==2 or prop_type==1){
      float inner_delta = abs(ttrack.innermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
      float outer_delta = abs(ttrack.outermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
      float used_delta = 0;

      if (inner_delta < outer_delta){
	tsos_seg = ttrack.innermostMeasurementState();
	tsos_ch = propagator->propagate(tsos_seg, ch->surface());
	used_delta = inner_delta;
	if (prop_type==1){data_.which_track = 1;}
	else{data_.which_track = 0;}
      }
      else{
	tsos_seg = ttrack.outermostMeasurementState();
	tsos_ch = propagator->propagate(tsos_seg, ch->surface());
	used_delta = outer_delta;
	if (prop_type==1){data_.which_track = 0;}
	else{data_.which_track = 1;}
      }
      
      if (tsos_ch.isValid()){
	const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
	const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
	if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
	  tmp_has_prop = true;
	  if (debug) {cout << "Delta to GEM = " << used_delta << "\tprop " << prop_type << std::endl;}
	  pos_GP = tsos_ch.globalPosition();
	  pos_startingPoint_GP = tsos_seg.globalPosition();
	}
      }
    }
  }

  if (prop_type == 3){
    LocalVector momentum_at_surface = ME11_segment->localDirection(); //No momentum for segments;
    DetId segDetId = ME11_segment->geographicalId();
    
    const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
    if (mu->isTrackerMuon()) {
      Track = mu->track().get();
      if (Track != 0) {
        tmp_inner_or_outer_mom = 0;
        momentum_at_surface = momentum_at_surface*(Track->outerP()); //If inner track exists, use momentum
        //if (debug) cout << "Got inner momentum = " << momentum_at_surface << endl;
      }
    }
    else if (mu->isGlobalMuon()){
      Track = mu->globalTrack().get();
      if (Track != 0){
        tmp_inner_or_outer_mom = 1;
        momentum_at_surface = momentum_at_surface*(Track->outerP()); //If no inner 
        //if (debug) cout << "Got outer momentum = " << momentum_at_surface << endl;
      }
      else{
        //if (debug) cout << "No tracks!" << endl;
        return;
      }
    }
    else{
      return;
    }
    //if (debug) cout << "Got momentum" << endl;
    LocalTrajectoryParameters param(ME11_segment->localPosition(), momentum_at_surface, mu->charge());
    AlgebraicSymMatrix mat(5,0);
    mat = ME11_segment->parametersError().similarityT( ME11_segment->projectionMatrix() );
    LocalTrajectoryError error(asSMatrix<5>(mat));
    TrajectoryStateOnSurface tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
    TrajectoryStateOnSurface tsos_ch = propagator->propagate(tsos_seg, ch->surface());
    if (tsos_ch.isValid()){
      //if (debug) cout << "propagate_to_gem tsos_ch ch station,chamber,layer,roll,region " << ch->id().station() << ch->id().chamber() << ch->id().layer() << ch->id().roll() << ch->id().region() << endl;
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      const LocalVector direction_local_ch = ch->toLocal(tsos_ch.globalDirection());
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
        tmp_has_prop = true;
        if (debug) cout << "tmp_has_prop " << tmp_has_prop << endl;
        pos_GP = tsos_ch.globalPosition();
        pos_startingPoint_GP = tsos_seg.globalPosition();
        prop_dxdz = direction_local_ch.x()/direction_local_ch.z();
      }
    }
  }

  //testing region:
  if (prop_type == 5 and RecoSegment.size()>0){
    //if (debug) cout << "prop 5 loop in propagate_to_gem, RecoSegment size:" << RecoSegment.size() << endl;
    //CSCSegmentCollection::const_iterator RecoSeg;
    //for (RecoSeg = cscSegmentsReco->begin(); RecoSeg !=cscSegmentsReco->end(); RecoSeg++){
    //LocalVector momentum_at_surface = RecoSegment->localDirection();
    //DetId segDetId = RecoSegment->geographicalId();
    LocalVector momentum_at_surface = RecoSegment[0]->localDirection();
    DetId segDetId = RecoSegment[0]->geographicalId();
    const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
    //if (debug) cout << "RecoSegment[0]->localDirection().x(): " <<  momentum_at_surface.x() << endl;
    if (mu->isTrackerMuon()){
      Track = mu->track().get();
      if (Track != 0) {
        tmp_inner_or_outer_mom = 0;
        momentum_at_surface = momentum_at_surface*(Track->outerP());
      }
    }
    else if (mu->isGlobalMuon()){
      Track = mu->globalTrack().get();
      if (Track != 0){
        tmp_inner_or_outer_mom = 1;
        momentum_at_surface = momentum_at_surface*(Track->outerP());
      }
      else {return;}
    }
    else {return;}

    //LocalTrajectoryParameters param(RecoSegment->localPosition(), momentum_at_surface, mu->charge());
    LocalTrajectoryParameters param(RecoSegment[0]->localPosition(), momentum_at_surface, mu->charge());
    AlgebraicSymMatrix mat(5,0);
    //mat = RecoSegment->parametersError().similarityT(RecoSegment->projectionMatrix());
    mat = RecoSegment[0]->parametersError().similarityT(RecoSegment[0]->projectionMatrix());
    LocalTrajectoryError error(asSMatrix<5>(mat));
    TrajectoryStateOnSurface tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
    TrajectoryStateOnSurface tsos_ch = propagator->propagate(tsos_seg, ch->surface());
    if(tsos_ch.isValid()){
      const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
      const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
      const LocalVector direction_local_ch = ch->toLocal(tsos_ch.globalDirection());
      if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1) {
	tmp_has_prop = true;
	if (debug) cout << "tmp_has_prop " << tmp_has_prop << endl;
	pos_GP = tsos_ch.globalPosition();
	pos_startingPoint_GP = tsos_seg.globalPosition();
	prop_dxdz = direction_local_ch.x()/direction_local_ch.z();
      }
    }
  }
  //end of testing region

  if (tmp_has_prop){
    //if (debug) cout << "Valid GEM prop" << endl;
    const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
    const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp(); //y distance to the current eta part
    const float prop_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2))).y();
    LocalPoint tmp_prop_LP = ch->toLocal(pos_GP);
    data_.prop_GP[0] = pos_GP.x(); data_.prop_GP[1] = pos_GP.y(); data_.prop_GP[2] = pos_GP.z();
    data_.prop_LP[0] = tmp_prop_LP.x(); data_.prop_LP[1] = tmp_prop_LP.y() + prop_y_to_chamber; data_.prop_LP[2] = tmp_prop_LP.z();
    data_.prop_dxdz = prop_dxdz;
    data_.prop_startingPoint_GP[0] = pos_startingPoint_GP.x(); data_.prop_startingPoint_GP[1] = pos_startingPoint_GP.y(); data_.prop_startingPoint_GP[2] = pos_startingPoint_GP.z();
    data_.prop_yroll = tmp_prop_LP.y();
    LocalPoint local_to_center(tmp_prop_LP.x(), tmp_prop_LP.y() + prop_y_to_center, 0);
    float local_phi = local_to_center.phi();
    data_.prop_localphi_rad = (3.14159265/2.) - local_phi;
    data_.prop_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
    data_.prop_globalphi_rad = pos_GP.phi();
    data_.has_prop = tmp_has_prop;
    data_.has_fidcut = fidcutCheck(tmp_prop_LP.y(), ((3.14159265/2.) - local_phi)*(180./3.14159265), ch);
    data_.prop_location[0] = ch->id().region(); data_.prop_location[1] = ch->id().station(); data_.prop_location[2] = ch->id().chamber(); data_.prop_location[3] = ch->id().layer(); data_.prop_location[4] = ch->id().roll();
    //if (debug) cout << "prop region:station:chamber:layer:roll " << data_.prop_location[0] << ":" << data_.prop_location[1] << ":"<< data_.prop_location[2] << ":" << data_.prop_location[3] << ":" << data_.prop_location[4] << endl; 
    data_.inner_or_outer_mom = tmp_inner_or_outer_mom;
    //if (debug) cout << "bunch of data.branches filled" << endl;
  }
}

void analyzer::GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_, int prop_type, std::vector<const CSCSegment*> RecoSegment){//const CSCSegment* RecoSegment){//const CSCSegmentCollection::const_iterator* RecoSegment){
  float tmp_rechit_GP_x; float tmp_rechit_GP_y; float tmp_rechit_GP_z; float tmp_rechit_GP_x_GE21; float tmp_rechit_GP_y_GE21; float tmp_rechit_GP_z_GE21;
  float tmp_rechit_LP_x; float tmp_rechit_LP_y; float tmp_rechit_LP_z; float tmp_rechit_LP_x_GE21; float tmp_rechit_LP_y_GE21; float tmp_rechit_LP_z_GE21;
  float tmp_rechit_yroll; float tmp_rechit_localphi_rad; float tmp_rechit_localphi_deg; float tmp_rechit_yroll_GE21; float tmp_rechit_localphi_rad_GE21; float tmp_rechit_localphi_deg_GE21;
  bool tmp_has_rechit = false; bool tmp_has_rechit_GE21 = false;
  int tmp_rechit_first_strip; int tmp_rechit_CLS; int tmp_rechit_BunchX; int tmp_rechit_first_strip_GE21; int tmp_rechit_CLS_GE21; int tmp_rechit_BunchX_GE21;
  float tmp_RdPhi = 9999.; float tmp_RdPhi_Corrected; int tmp_rechit_detId; float tmp_RdPhi_GE21 = 9999.; float tmp_RdPhi_Corrected_GE21; int tmp_rechit_detId_GE21;
  float tmp_dPhi = 9999.; float tmp_dPhi_Corrected; float tmp_dPhi_GE21 = 9999.; float tmp_dPhi_Corrected_GE21;
  float tmp_bending_angle = 9999.; float tmp_bending_angle_GE21 = 9999.;
  int tmp_nRecHitsTot = 0; int tmp_nRecHits5 = 0; int tmp_nRecHits2 = 0; int tmp_nRecHitsTot_GE21 = 0; int tmp_nRecHits5_GE21 = 0; int tmp_nRecHits2_GE21 = 0;
  int tmp_rechit_region; int tmp_rechit_station; int tmp_rechit_chamber; int tmp_rechit_layer; int tmp_rechit_roll; int tmp_rechit_region_GE21; int tmp_rechit_station_GE21; int tmp_rechit_chamber_GE21; int tmp_rechit_layer_GE21; int tmp_rechit_roll_GE21;
  int tmp_nRecHitsRpos1L1 = 0; int tmp_nRecHitsRpos1L2 = 0; int tmp_nRecHitsRneg1L1 = 0; int tmp_nRecHitsRneg1L2 = 0; int tmp_nRecHitsRpos1L1_GE21 = 0; int tmp_nRecHitsRpos1L2_GE21 = 0; int tmp_nRecHitsRneg1L1_GE21 = 0; int tmp_nRecHitsRneg1L2_GE21 = 0;

  //if (debug) cout << "gem_rechit_matcher ch: station,chamber, layer,roll,region" << ch->id().station() << ch->id().chamber() << ch->id().layer() << ch->id().roll() << ch->id().region() << endl;
  if (debug) cout << "gemRecHits->size(): " << gemRecHits->size() << endl;
  for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++) {
    //if (debug) cout << "gemRecHits loop within the rechit matcher function" << endl;
    if ((hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM) {
      //if (debug) cout << "getting detector info" << endl;
      GEMDetId gemid((hit)->geographicalId());
      if (gemid.region() == 1) {
        if (gemid.layer() == 1) {tmp_nRecHitsRpos1L1++;}
        if (gemid.layer() == 2) {tmp_nRecHitsRpos1L2++;}
      }
      if (gemid.region() == -1) {
        if (gemid.layer() == 1) {tmp_nRecHitsRneg1L1++;}
        if (gemid.layer() == 2) {tmp_nRecHitsRneg1L2++;}
      }
      int test_st = ch->id().station(); int test_ch = ch->id().chamber(); int test_lay = ch->id().layer(); int test_rol = ch->id().roll(); int test_reg = ch->id().region();
      if (debug) cout << "gemid: station,chamber,layer,roll,region " << gemid.station() << gemid.chamber() << gemid.layer() << gemid.roll() << gemid.region() << " ch " << test_st << test_ch << test_lay << test_rol << test_reg << endl;

      if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()) {
        if (debug) cout << "matching gemid chamber with ch" << endl;
        const auto& etaPart = GEMGeometry_->etaPartition(gemid);
        if (debug) cout << "etaPart" << endl;
        float strip = etaPart->strip(hit->localPosition());
        float stripAngle = etaPart->specificTopology().stripAngle(strip);
        if (debug) cout << "stripAngle: " << stripAngle << endl;
        float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
        float rechit_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y();
        LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
        //below is for GE1/1
        if (ch->id().station()==1 and ch->id().ring()==1 and fabs((hit)->localPosition().x() - prop_LP.x())<999.0) {
          tmp_nRecHitsTot++;
          if (abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 5) {tmp_nRecHits5++;}
          if (abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 2) {tmp_nRecHits2++;}
          if (abs(tmp_RdPhi) > abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch))) {
            tmp_rechit_GP_x = etaPart->toGlobal((hit)->localPosition()).x();
            tmp_rechit_GP_y = etaPart->toGlobal((hit)->localPosition()).y();
            tmp_rechit_GP_z = etaPart->toGlobal((hit)->localPosition()).z();
            tmp_rechit_LP_x = (hit)->localPosition().x();   tmp_rechit_LP_y = rechit_y_to_chamber + (hit)->localPosition().y();    tmp_rechit_LP_z = (hit)->localPosition().z();
	    tmp_rechit_yroll = (hit)->localPosition().y();
	    float local_phi = local_to_center.phi(); 
	    tmp_rechit_localphi_rad = (3.14159265/2.) - local_phi;
	    tmp_rechit_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
	    tmp_has_rechit = true;
            tmp_rechit_first_strip = (hit)->firstClusterStrip();
	    tmp_rechit_CLS = (hit)->clusterSize();
	    tmp_rechit_BunchX = (hit)->BunchX();

            //Calculating the bending angle = CSC segment phi - GEM rechit phi
            if (data_.hasME11) {
              float CSC_segment_phi;
              if (prop_type!=5){
                DetId segDetId = ME11_segment->geographicalId();
                const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
                CSC_segment_phi = (segDet->toGlobal(ME11_segment->localPosition())).phi();
                if (debug) cout << "CSC_segment_phi: " << CSC_segment_phi << endl;
              }

              if (prop_type==5){
		//CSCSegmentCollection::const_iterator RecoSeg;
                //for (RecoSeg = cscSegmentsReco->begin(); RecoSeg !=cscSegmentsReco->end(); RecoSeg++){
                //DetId segDetId = RecoSegment->geographicalId();
                DetId segDetId = RecoSegment[0]->geographicalId();
                const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
                //CSC_segment_phi = (segDet->toGlobal(RecoSegment->localPosition())).phi();
                CSC_segment_phi = (segDet->toGlobal(RecoSegment[0]->localPosition())).phi();
                if (debug) cout<< "CSC_segment_phi: " << CSC_segment_phi << endl;
              }

              float GEM_hit_phi = (etaPart->toGlobal(hit->localPosition())).phi();
              tmp_bending_angle = CSC_segment_phi - GEM_hit_phi;
              if (debug) cout << "Bending Angle = " << tmp_bending_angle << "\tpT = " << data_.muon_pt << endl;
            }
            
            tmp_RdPhi = RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch);
            tmp_RdPhi_Corrected = tmp_RdPhi;
            tmp_dPhi = tmp_rechit_localphi_rad - data_.prop_localphi_rad; //units of radian
            tmp_dPhi_Corrected = tmp_dPhi;
            if ((gemid.region() == 1 and gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)) {
              tmp_RdPhi_Corrected = -1.0*tmp_RdPhi_Corrected;
            }
            if ((gemid.region() == -1 and gemid.chamber()%2 == 1) || (gemid.region() == 1 && gemid.chamber()%2 == 0)) {
              tmp_dPhi_Corrected = -1.0*tmp_dPhi_Corrected;
            }
            /*
            Comment on SIgn Convention: notice that for the residual calculation in phi we have dphi = rechit_phi - prop_phi.
            Also dphi_corrected is flipped from what rdphi_corrected was. these are because in firmware, the -endcap's phi matches the global phi.
	     */
            tmp_rechit_detId = gemid.region()*(gemid.station()*100 + gemid.chamber());
            tmp_rechit_region = gemid.region();  tmp_rechit_station = gemid.station();
            tmp_rechit_chamber = gemid.chamber();  tmp_rechit_layer = gemid.layer();
            tmp_rechit_roll = gemid.roll();
            if (debug) cout << "rechit_detId:RdPhi:RdPhi_Corrected:dPhi:dPhi_Corrected\t" << tmp_rechit_detId << ":" << tmp_RdPhi << ":" << tmp_RdPhi_Corrected << ":" << tmp_dPhi << ":" << tmp_dPhi_Corrected << endl;
          }
        }
        //below is for GE2/1:(incomplete)
        if (ch->id().station()==2 and ch->id().ring()==1 and fabs((hit)->localPosition().x() - prop_LP.x())<999.0) {
          tmp_nRecHitsTot_GE21++;
          if (abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 5) {tmp_nRecHits5_GE21++;}
          if (abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 2) {tmp_nRecHits2_GE21++;}
          if (abs(tmp_RdPhi) > abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch))) {
            tmp_rechit_GP_x_GE21 = etaPart->toGlobal((hit)->localPosition()).x();
            tmp_rechit_GP_y_GE21 = etaPart->toGlobal((hit)->localPosition()).y();
            tmp_rechit_GP_z_GE21 = etaPart->toGlobal((hit)->localPosition()).z();
            tmp_rechit_LP_x_GE21 = (hit)->localPosition().x();   tmp_rechit_LP_y_GE21 = rechit_y_to_chamber + (hit)->localPosition().y();    tmp_rechit_LP_z_GE21 = (hit)->localPosition().z();
            tmp_rechit_yroll_GE21 = (hit)->localPosition().y();
            float local_phi_GE21 = local_to_center.phi();
            tmp_rechit_localphi_rad_GE21 = (3.14159265/2.) - local_phi_GE21;
            tmp_rechit_localphi_deg_GE21 = ((3.14159265/2.) - local_phi_GE21)*(180./3.14159265);
            tmp_has_rechit_GE21 = true;
            tmp_rechit_first_strip_GE21 = (hit)->firstClusterStrip();
            tmp_rechit_CLS_GE21 = (hit)->clusterSize();
            tmp_rechit_BunchX_GE21 = (hit)->BunchX();

            //Calculating the bending angle = CSC segment phi - GEM rechit phi
            if (data_.hasME21) {
              float CSC_segment_phi;
              if (prop_type!=5){
                DetId segDetId = ME11_segment->geographicalId();
                const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
                CSC_segment_phi = (segDet->toGlobal(ME11_segment->localPosition())).phi();
                if (debug) cout << "CSC_segment_phi: " << CSC_segment_phi << endl;
              }

              if (prop_type==5){
		//CSCSegmentCollection::const_iterator RecoSeg;
                //for (RecoSeg = cscSegmentsReco->begin(); RecoSeg !=cscSegmentsReco->end(); RecoSeg++){
                //DetId segDetId = RecoSegment->geographicalId();
                DetId segDetId = RecoSegment[0]->geographicalId();
                const GeomDet* segDet = theTrackingGeometry->idToDet(segDetId);
                //CSC_segment_phi = (segDet->toGlobal(RecoSegment->localPosition())).phi();
                CSC_segment_phi = (segDet->toGlobal(RecoSegment[0]->localPosition())).phi();
                //if (debug) cout<< "CSC_segment_phi: " << CSC_segment_phi << endl;
              }

              float GEM_hit_phi = (etaPart->toGlobal(hit->localPosition())).phi();
              tmp_bending_angle = CSC_segment_phi - GEM_hit_phi;
              if (debug) cout << "Bending Angle = " << tmp_bending_angle << "\tpT = " << data_.muon_pt << endl;
            }
	    
            tmp_RdPhi_GE21 = RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch);
            tmp_RdPhi_Corrected_GE21 = tmp_RdPhi_GE21;
            tmp_dPhi_GE21 = tmp_rechit_localphi_rad_GE21 - data_.prop_localphi_rad; //units of radian
            tmp_dPhi_Corrected_GE21 = tmp_dPhi_GE21;
            //check the sign convention below with GE2/1 Firmware
            if ((gemid.region() == 1 and gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)) {
              tmp_RdPhi_Corrected_GE21 = -1.0*tmp_RdPhi_Corrected_GE21;
            }
            if ((gemid.region() == -1 and gemid.chamber()%2 == 1) || (gemid.region() == 1 && gemid.chamber()%2 == 0)) {
              tmp_dPhi_Corrected_GE21 = -1.0*tmp_dPhi_Corrected_GE21;
            }
            tmp_rechit_detId_GE21 = gemid.region()*(gemid.station()*100 + gemid.chamber());
            tmp_rechit_region_GE21 = gemid.region();  tmp_rechit_station_GE21 = gemid.station();
            tmp_rechit_chamber_GE21 = gemid.chamber();  tmp_rechit_layer_GE21 = gemid.layer();
            tmp_rechit_roll_GE21 = gemid.roll();
            if (debug) cout << "GE21 rechit_detId:RdPhi:RdPhi_Corrected:dPhi:dPhi_Corrected\t" << tmp_rechit_detId_GE21 << ":" << tmp_RdPhi_GE21 << ":" << tmp_RdPhi_Corrected_GE21 << ":" << tmp_dPhi_GE21 << ":" << tmp_dPhi_Corrected_GE21 << endl;
          }
        }
      }
    }
  }
  if (tmp_has_rechit){
    //below is for GE11
    data_.rechit_GP[0] = tmp_rechit_GP_x; data_.rechit_GP[1] = tmp_rechit_GP_y; data_.rechit_GP[2] = tmp_rechit_GP_z;
    data_.rechit_LP[0] = tmp_rechit_LP_x; data_.rechit_LP[1] = tmp_rechit_LP_y; data_.rechit_LP[2] = tmp_rechit_LP_z;
    data_.rechit_yroll = tmp_rechit_yroll;
    data_.rechit_localphi_rad = tmp_rechit_localphi_rad;
    data_.rechit_localphi_deg = tmp_rechit_localphi_deg;
    data_.has_rechit = tmp_has_rechit;
    data_.rechit_first_strip = tmp_rechit_first_strip;
    data_.rechit_CLS = tmp_rechit_CLS;
    data_.rechit_BunchX = tmp_rechit_BunchX;
    data_.RdPhi = tmp_RdPhi;                      data_.dPhi = tmp_dPhi;
    data_.RdPhi_Corrected = tmp_RdPhi_Corrected;  data_.dPhi_Corrected = tmp_dPhi_Corrected;
    data_.rechit_detId = tmp_rechit_detId;
    data_.bending_angle = tmp_bending_angle;
    data_.nRecHitsTot = tmp_nRecHitsTot; data_.nRecHits5 = tmp_nRecHits5; data_.nRecHits2 = tmp_nRecHits2;
    data_.rechit_location[0] = tmp_rechit_region; data_.rechit_location[1] = tmp_rechit_station; data_.rechit_location[2] = tmp_rechit_chamber; data_.rechit_location[3] = tmp_rechit_layer; data_.rechit_location[4] = tmp_rechit_roll;
    data_.nRecHitsRpos1L1 = tmp_nRecHitsRpos1L1; data_.nRecHitsRpos1L2 = tmp_nRecHitsRpos1L2; data_.nRecHitsRneg1L1 = tmp_nRecHitsRneg1L1; data_.nRecHitsRneg1L2 = tmp_nRecHitsRneg1L2;
    //below is for GE21
    data_.rechit_GP_GE21[0] = tmp_rechit_GP_x_GE21; data_.rechit_GP_GE21[1] = tmp_rechit_GP_y_GE21; data_.rechit_GP_GE21[2] = tmp_rechit_GP_z_GE21;
    data_.rechit_LP_GE21[0] = tmp_rechit_LP_x_GE21; data_.rechit_LP_GE21[1] = tmp_rechit_LP_y_GE21; data_.rechit_LP_GE21[2] = tmp_rechit_LP_z_GE21;
    data_.rechit_yroll_GE21 = tmp_rechit_yroll_GE21;
    data_.rechit_localphi_rad_GE21 = tmp_rechit_localphi_rad_GE21;
    data_.rechit_localphi_deg_GE21 = tmp_rechit_localphi_deg_GE21;
    data_.has_rechit_GE21 = tmp_has_rechit_GE21;
    data_.rechit_first_strip_GE21 = tmp_rechit_first_strip_GE21;
    data_.rechit_CLS_GE21 = tmp_rechit_CLS_GE21;
    data_.rechit_BunchX_GE21 = tmp_rechit_BunchX_GE21;
    data_.RdPhi_GE21 = tmp_RdPhi_GE21;                      data_.dPhi_GE21 = tmp_dPhi_GE21;
    data_.RdPhi_Corrected_GE21 = tmp_RdPhi_Corrected_GE21;  data_.dPhi_Corrected_GE21 = tmp_dPhi_Corrected_GE21;
    data_.rechit_detId_GE21 = tmp_rechit_detId_GE21;
    data_.bending_angle_GE21 = tmp_bending_angle_GE21;
    data_.nRecHitsTot_GE21 = tmp_nRecHitsTot_GE21; data_.nRecHits5_GE21 = tmp_nRecHits5_GE21; data_.nRecHits2_GE21 = tmp_nRecHits2_GE21;
    data_.rechit_location_GE21[0] = tmp_rechit_region_GE21; data_.rechit_location_GE21[1] = tmp_rechit_station_GE21; data_.rechit_location_GE21[2] = tmp_rechit_chamber_GE21; data_.rechit_location_GE21[3] = tmp_rechit_layer_GE21; data_.rechit_location_GE21[4] = tmp_rechit_roll_GE21;
    data_.nRecHitsRpos1L1_GE21 = tmp_nRecHitsRpos1L1_GE21; data_.nRecHitsRpos1L2_GE21 = tmp_nRecHitsRpos1L2_GE21; data_.nRecHitsRneg1L1_GE21 = tmp_nRecHitsRneg1L1_GE21; data_.nRecHitsRneg1L2_GE21 = tmp_nRecHitsRneg1L2_GE21;
  }
}

void analyzer::GEM_simhit_matcher(const GEMEtaPartition* ch, GlobalPoint prop_GP, MuonData& data_){
  float tmpDy = 999.; float tmpDr = 999.; int tmpSimCounter = 0;
  float tmp_sim_GP_x; float tmp_sim_GP_y; float tmp_sim_GP_z;
  float tmp_sim_LP_x; float tmp_sim_LP_y; float tmp_sim_LP_z;
  bool has_tmp = false;
  for (const auto& simHit:*gemSimHits.product()){
    GEMDetId gemid((simHit).detUnitId());
    if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
      tmpSimCounter++;
      const auto& etaPart = GEMGeometry_->etaPartition(gemid);
      float dy = prop_GP.y() - etaPart->toGlobal(simHit.localPosition()).y();
      float dx = prop_GP.x() - etaPart->toGlobal(simHit.localPosition()).x();
      if (dy < tmpDy) tmpDy = dy;
      if (pow(pow(dy, 2) + pow(dx, 2), 0.5) < tmpDr){
	tmp_sim_GP_x = etaPart->toGlobal(simHit.localPosition()).x();
	tmp_sim_GP_y = etaPart->toGlobal(simHit.localPosition()).y();
	tmp_sim_GP_z = etaPart->toGlobal(simHit.localPosition()).z();
	tmp_sim_LP_x = simHit.localPosition().x();
	tmp_sim_LP_y = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y() + simHit.localPosition().y();
        tmp_sim_LP_z = simHit.localPosition().z();
	tmpDr = pow(pow(dy, 2) + pow(dx, 2), 0.5);
        has_tmp = true;
      }
    }
  }
  if (has_tmp){
    data_.sim_GP[0] = tmp_sim_GP_x; data_.sim_GP[1] = tmp_sim_GP_y; data_.sim_GP[2] = tmp_sim_GP_z;
    data_.sim_LP[0] = tmp_sim_LP_x; data_.sim_LP[1] = tmp_sim_LP_y; data_.sim_LP[2] = tmp_sim_LP_z;
    data_.simDy = tmpDy;
    data_.nSim = tmpSimCounter;
  }
}

float analyzer::RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch){
  GEMDetId gemid((rechit)->geographicalId());
  const auto& etaPart = GEMGeometry_->etaPartition(gemid); //eta partition of the reconstructed hit location
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id()); //eta partition of the propagated hit location
  float deltay_roll = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp(); //global position of the center of the propagated y and subtract the rechit chamber center eta
  return cos(stripAngle) * (prop_localx - (rechit)->localPosition().x()) - sin(stripAngle) * (prop_localy + deltay_roll);
  /*
  RdPhi clarification: The residual calculation is RdPhi = cos(Angle) * delta_x + sin(Angle) * delta_y.
  Mapping of the angle in GEM is clockwise, whereas this equation considers counterclockwise strip angle to be positive.
  Since cosine is an even function, no sign change is needed for the first part of residual calculation.
  But to the second part, for sine being an odd function, we add the minus sign to account for the CW angles. 
   */
}

bool analyzer::fidcutCheck(float local_y, float localphi_deg, const GEMEtaPartition* ch){
  const float fidcut_angle = 1.0;
  const float cut_chamber = 5.0;
  const float cut_angle = cut_chamber - fidcut_angle;
  auto& parameters(ch->specs()->parameters());
  float height(parameters[2]);
  if ((abs(localphi_deg) < cut_angle) &&
      ((local_y < (height - cut_chamber) && ch->id().roll() == 1) || (local_y > -1.0*(height - cut_chamber) && ch->id().roll() == 8) || (ch->id().roll() != 1 && ch->id().roll() != 8))
     )
    {return 1;}
  else {return 0;}
}

void analyzer::beginJob(){}
void analyzer::endJob(){
  if (debug) nME11_col_vs_matches->Write();
}

DEFINE_FWK_MODULE(analyzer);

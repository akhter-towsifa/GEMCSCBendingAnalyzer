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
  TTree* book(TTree *t, int prop_type);
  //Muon Info//////////////////////////////////////////////////////
  int muon_charge; float muon_pt; float muon_eta; float muon_momentum;
  unsigned long long  evtNum; unsigned long long  lumiBlock; int muonIdx;
  //Propagation Info//////////////////////////////////////////////////////
  std::vector<float> prop_GP_x; std::vector<float> prop_GP_y; std::vector<float> prop_GP_z;
  std::vector<float> prop_LP_x; std::vector<float> prop_LP_y; std::vector<float> prop_LP_z;
  std::vector<float> prop_startingPoint_GP_x; std::vector<float> prop_startingPoint_GP_y; std::vector<float> prop_startingPoint_GP_z;
  std::vector<float> prop_yroll; std::vector<float> prop_localphi_rad; std::vector<float> prop_localphi_deg;
  std::vector<bool> has_prop; std::vector<bool> has_fidcut;
  std::vector<int> prop_region; std::vector<int> prop_station; std::vector<int> prop_chamber; std::vector<int> prop_layer; std::vector<int> prop_roll;
  //Track Info//////////////////////////////////////////////////////
  float track_chi2; float track_ndof; int n_ME11_segment; int which_track;
  int hasME11; int hasME11RecHit; int hasME11A; int hasME11ARecHit;
  int nCSCSeg; int nDTSeg;
  //Rechit Info//////////////////////////////////////////////////////
  std::vector<float> rechit_GP_x; std::vector<float> rechit_GP_y; std::vector<float> rechit_GP_z;
  std::vector<float> rechit_LP_x; std::vector<float> rechit_LP_y; std::vector<float> rechit_LP_z;
  std::vector<float> rechit_yroll; std::vector<float> rechit_localphi_rad; std::vector<float> rechit_localphi_deg;
  std::vector<bool> has_rechit;
  std::vector<int> rechit_first_strip; std::vector<int> rechit_CLS; std::vector<int> rechit_BunchX;
  std::vector<float> RdPhi; std::vector<float> RdPhi_Corrected; std::vector<int> rechit_detId;
  std::vector<int> nRecHitsTot; std::vector<int> nRecHits5; std::vector<int> nRecHits2;
  std::vector<int> rechit_region; std::vector<int> rechit_station; std::vector<int> rechit_chamber; std::vector<int> rechit_layer; std::vector<int> rechit_roll;
  //Sim info for MC
  float sim_GP[3]; float sim_LP[3]; float sim_localy_roll; int nSim;
};

void MuonData::init()
{
  //Muon Info//////////////////////////////////////////////////////
  muon_charge = 9999; muon_pt = 9999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999;
  //Propagation Info//////////////////////////////////////////////////////
  prop_GP_x.clear(); prop_GP_y.clear(); prop_GP_z.clear();
  prop_LP_x.clear(); prop_LP_y.clear(); prop_LP_z.clear();
  prop_startingPoint_GP_x.clear(); prop_startingPoint_GP_y.clear(); prop_startingPoint_GP_z.clear();
  prop_yroll.clear(); prop_localphi_rad.clear(); prop_localphi_deg.clear();
  has_prop.clear(); has_fidcut.clear();
  prop_region.clear(); prop_station.clear(); prop_chamber.clear(); prop_layer.clear(); prop_roll.clear();
  //Track Info//////////////////////////////////////////////////////
  track_chi2 = 999999; track_ndof = 999999; n_ME11_segment = 999999; which_track = 999999;
  hasME11 = 999999; hasME11RecHit = 999999; hasME11A = 999999; hasME11ARecHit = 999999;
  nCSCSeg = 999999; nDTSeg = 999999;
  //Rechit Info//////////////////////////////////////////////////////
  rechit_GP_x.clear(); rechit_GP_y.clear(); rechit_GP_z.clear();
  rechit_LP_x.clear(); rechit_LP_y.clear(); rechit_LP_z.clear();
  rechit_yroll.clear(); rechit_localphi_rad.clear(); rechit_localphi_deg.clear();
  has_rechit.clear();
  rechit_first_strip.clear(); rechit_CLS.clear(); rechit_BunchX.clear();
  RdPhi.clear(); RdPhi_Corrected.clear(); rechit_detId.clear();
  nRecHitsTot.clear(); nRecHits5.clear(); nRecHits2.clear();
  rechit_region.clear(); rechit_station.clear(); rechit_chamber.clear(); rechit_layer.clear(); rechit_roll.clear();
  //Sim info for MC
  for (int j=0; j<3; ++j){
    sim_GP[j] = 9999999; sim_LP[j] = 9999999;
  }
  sim_localy_roll = 9999999; nSim = 9999999;
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
  else{
    std::cout << "Bad prop type, failure" << std::endl;
  }
  //Muon Info//////////////////////////////////////////////////////
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta); t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum); t->Branch("lumiBlock", &lumiBlock); t->Branch("muonIdx", &muonIdx);
  //Propagation Info//////////////////////////////////////////////////////
  t->Branch("prop_GP_x", &prop_GP_x);
  t->Branch("prop_GP_y", &prop_GP_y);
  t->Branch("prop_GP_z", &prop_GP_z);
  t->Branch("prop_LP_x", &prop_LP_x);
  t->Branch("prop_LP_y", &prop_LP_y);
  t->Branch("prop_LP_z", &prop_LP_z);
  t->Branch("prop_startingPoint_GP_x", &prop_startingPoint_GP_x);
  t->Branch("prop_startingPoint_GP_y", &prop_startingPoint_GP_y);
  t->Branch("prop_startingPoint_GP_z", &prop_startingPoint_GP_z);
  t->Branch("prop_yroll", &prop_yroll);
  t->Branch("prop_localphi_rad", &prop_localphi_rad);
  t->Branch("prop_localphi_deg", &prop_localphi_deg);
  t->Branch("has_prop", &has_prop);
  t->Branch("has_fidcut", &has_fidcut);
  t->Branch("prop_region", &prop_region);
  t->Branch("prop_station", &prop_station);
  t->Branch("prop_chamber", &prop_chamber);
  t->Branch("prop_layer", &prop_layer);
  t->Branch("prop_roll", &prop_roll);
  //Track Info//////////////////////////////////////////////////////
  t->Branch("track_chi2", &track_chi2); t->Branch("track_ndof", &track_ndof);
  t->Branch("n_ME11_segment", &n_ME11_segment); t->Branch("which_track", &which_track);
  t->Branch("hasME11", &hasME11); t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A); t->Branch("hasME11ARecHit", &hasME11ARecHit);
  t->Branch("nCSCSeg", &nCSCSeg); t->Branch("nDTSeg", &nDTSeg);
  //Rechit Info//////////////////////////////////////////////////////
  t->Branch("rechit_GP_x", &rechit_GP_x);
  t->Branch("rechit_GP_y", &rechit_GP_y);
  t->Branch("rechit_GP_z", &rechit_GP_z);
  t->Branch("rechit_LP_x", &rechit_LP_x);
  t->Branch("rechit_LP_y", &rechit_LP_y);
  t->Branch("rechit_LP_z", &rechit_LP_z);
  t->Branch("rechit_yroll", &rechit_yroll);
  t->Branch("rechit_localphi_rad", &rechit_localphi_rad);
  t->Branch("rechit_localphi_deg", &rechit_localphi_deg);
  t->Branch("has_rechit", &has_rechit);
  t->Branch("rechit_first_strip", &rechit_first_strip);
  t->Branch("rechit_CLS", &rechit_CLS);
  t->Branch("rechit_BunchX", &rechit_BunchX);
  t->Branch("RdPhi", &RdPhi);
  t->Branch("RdPhi_Corrected", &RdPhi_Corrected);
  t->Branch("rechit_detId", &rechit_detId);
  t->Branch("nRecHitsTot", &nRecHitsTot);
  t->Branch("nRecHits2", &nRecHits2);
  t->Branch("nRecHits5", &nRecHits5);
  t->Branch("rechit_region", &rechit_region);
  t->Branch("rechit_station", &rechit_station);
  t->Branch("rechit_chamber", &rechit_chamber);
  t->Branch("rechit_layer", &rechit_layer);
  t->Branch("rechit_roll", &rechit_roll);
  //Sim info for MC
  float sim_GP[3]; float sim_LP[3]; float sim_localy_roll; int nSim;
  t->Branch("sim_GP", &sim_GP, "sim_GP[3] (x,y,z)/F");
  t->Branch("sim_LP", &sim_LP, "sim_LP[3] (x,y,z)/F");
  t->Branch("sim_localy_roll", &sim_localy_roll);
  t->Branch("nSim", &nSim);
  return t;
}


void propagate(const reco::Muon* mu, int prop_type);

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
  void propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i);
  void CSCSegmentCounter(const reco::Muon* mu, int& n_ME11_segment, int& nCSCSeg, int& nDTSeg);
  //void propagate_to_GEM(reco::TransientTrack track, const GEMEtaPartition* ch, int prop_type, GlobalPoint &pos_GP, GlobalPoint &pos_startingPoint_GP, bool &has_prop, bool &in_or_out);
  void propagate_to_GEM(reco::TransientTrack track, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_);
  void GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_);
  float RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch);

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::Handle<vector<PSimHit> > gemSimHits;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;

  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  edm::ESHandle<GEMGeometry> GEMGeometry_;

  bool tracker_prop;
  bool CSC_prop;
  bool Segment_prop;
  bool debug;

  MuonData data_;
  TTree* CSC_tree; TTree* Tracker_tree; TTree* Segment_tree;

  const CSCSegment *ME11_segment;
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
  Segment_prop = iConfig.getParameter<bool>("Segment_prop");
  debug = iConfig.getParameter<bool>("debug");
  cout << "tracker_prop " << tracker_prop << " CSC_prop " << CSC_prop << " debug " << debug << std::endl;

  bool CSC_prop = 1;
  if(CSC_prop){
    CSC_tree = data_.book(CSC_tree, 1);
  }
  if(tracker_prop){
    Tracker_tree = data_.book(Tracker_tree, 2);
  }
  if(Segment_prop){
    Segment_tree = data_.book(Segment_tree, 3);
  }
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
  iEvent.getByToken(gemRecHits_, gemRecHits);
  if (isMC) {
    iEvent.getByToken(gemSimHits_, gemSimHits); 
  }
  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;
  if (muons->size() == 0) return;

  cout << "new evt numb is " << iEvent.eventAuxiliary().event() << " and new lumiblock is " << iEvent.eventAuxiliary().luminosityBlock() << endl;
  for (size_t i = 0; i < muons->size(); ++i){
    //cout << "new muon" << endl;
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    //if (mu->pt() < 2.0) continue;  //can apply a pt cut later
    if (not mu->standAloneMuon()) continue;
    cout << "new standalone" << endl;
    int prop_type = 1;
    propagate(mu, prop_type, iEvent, i);

  }
}


void propagate_segment(const CSCSegment* ME11_segment, const reco::Track* innerTrack, const GEMEtaPartition* ch, const reco::Muon* mu, MuonServiceProxy* theService_, const GeomDet* segDet, GlobalPoint &pos_global_ch, GlobalPoint &pos_global_seg, bool &has_prop, bool debug){
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  LocalVector momentum_at_surface = ME11_segment->localDirection(); //No momentum for segments
  if (innerTrack != 0){
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


float analyser::RdPhi_func(float stripAngle, const edm::OwnVector<GEMRecHit, edm::ClonePolicy<GEMRecHit> >::const_iterator rechit, float prop_localx, float prop_localy, const GEMEtaPartition* ch){
  GEMDetId gemid((rechit)->geographicalId());
  const auto& etaPart = GEMGeometry_->etaPartition(gemid);
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  float deltay_roll =  etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp() - etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
  return cos(stripAngle) * (prop_localx - (rechit)->localPosition().x()) + sin(stripAngle) * (prop_localy + deltay_roll);
}
void analyser::CSCSegmentCounter(const reco::Muon* mu, int& n_ME11_segment, int& nCSCSeg, int& nDTSeg){
  const reco::Track* Track = mu->outerTrack().get();
  int tmp_CSC_counter = 0; int tmp_DT_counter = 0; int tmp_ME11_counter = 0;
  for (size_t RecHit_iter = 0; RecHit_iter != Track->recHitsSize(); RecHit_iter++){
    const TrackingRecHit* RecHit = (Track->recHit(RecHit_iter)).get();
    DetId RecHitId = RecHit->geographicalId();
    uint16_t RecHitDetId = RecHitId.det();
    if (RecHitDetId == DetId::Muon){
      uint16_t RecHitSubDet = RecHitId.subdetId();
      if (RecHitSubDet == (uint16_t)MuonSubdetId::CSC){
        if (CSCDetId(RecHitId).station() == 1 and CSCDetId(RecHitId).ring() == 1 and RecHit->dimension() == 4){tmp_ME11_counter++;}
        if (RecHit->dimension() == 4){tmp_CSC_counter++;}
      }
      if (RecHitSubDet == (uint16_t)MuonSubdetId::DT){
        if (RecHit->dimension() > 1){tmp_DT_counter++;}
      }
    }
  }
  nCSCSeg = tmp_CSC_counter; nDTSeg = tmp_DT_counter;
  n_ME11_segment = tmp_ME11_counter;
  std::cout << "Counted nSeg = " << nCSCSeg << std::endl;
}
//void analyser::propagate_to_GEM(reco::TransientTrack track, const GEMEtaPartition* ch, int prop_type, GlobalPoint &pos_GP, GlobalPoint &pos_startingPoint_GP, bool &has_prop, bool &in_or_out){
void analyser::propagate_to_GEM(reco::TransientTrack track, const GEMEtaPartition* ch, int prop_type, bool &tmp_has_prop, GlobalPoint &pos_GP, MuonData& data_){
  tmp_has_prop = false; bool tmp_has_fidcut = 0;
  const BoundPlane& bps(ch->surface());
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");
  const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
  TrajectoryStateOnSurface tsos_ch; TrajectoryStateOnSurface tsos_seg;
  GlobalPoint pos_startingPoint_GP;
  if(prop_type == 1 or prop_type == 2){
    float inner_delta = abs(track.innermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
    float outer_delta = abs(track.outermostMeasurementState().globalPosition().z() - GEMGeometry_->etaPartition(ch->id())->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).z());
    if (inner_delta < outer_delta){
      tsos_seg = track.innermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());
      if(prop_type == 1){data_.which_track = 1;}
      else{data_.which_track = 0;}
    }
    else{
      tsos_seg = track.outermostMeasurementState(); tsos_ch = propagator->propagate(tsos_seg, ch->surface());
      if(prop_type == 1){data_.which_track = 0;}
      else{data_.which_track = 1;}
    }
  }
  if(prop_type == 3){
    track = mu->track().get();
    LocalVector momentum_at_surface = ME11_segment->localDirection(); //No momentum for segments
    if (track != 0){
      momentum_at_surface = momentum_at_surface*(track->outerP()); //If innerTrack exists, use momentum
    }
    LocalTrajectoryParameters param(ME11_segment->localPosition(), momentum_at_surface, mu->charge());
    AlgebraicSymMatrix mat(5,0);
    mat = ME11_segment->parametersError().similarityT( ME11_segment->projectionMatrix() );
    LocalTrajectoryError error(asSMatrix<5>(mat));
    tsos_seg(param, error, segDet->surface(), &*theService_->magneticField());
    tsos_ch = propagator->propagate(tsos_seg, ch->surface());
  }
}




  }
  if (tsos_ch.isValid()){
    const LocalPoint pos_local_ch = ch->toLocal(tsos_ch.globalPosition());
    const LocalPoint pos2D_local_ch(pos_local_ch.x(), pos_local_ch.y(), 0);
    if (!(tsos_ch.globalPosition().z() * tsos_seg.globalPosition().z() < 0) and bps.bounds().inside(pos2D_local_ch) and ch->id().station() == 1 and ch->id().ring() == 1){
      tmp_has_prop = true;
      pos_GP = tsos_ch.globalPosition();
      pos_startingPoint_GP = tsos_seg.globalPosition();
    }
  }
  if(tmp_has_prop){
    const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
    const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp(); //y distance to the current eta part
    const float prop_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2))).y();
    data_.prop_GP_x.push_back(pos_GP.x()); data_.prop_GP_y.push_back(pos_GP.y()); data_.prop_GP_z.push_back(pos_GP.z());
    LocalPoint tmp_prop_LP = ch->toLocal(pos_GP);
    data_.prop_LP_x.push_back(tmp_prop_LP.x()); data_.prop_LP_y.push_back(tmp_prop_LP.y() + prop_y_to_chamber); data_.prop_LP_z.push_back(tmp_prop_LP.z());
    data_.prop_startingPoint_GP_x.push_back(pos_startingPoint_GP.x()); data_.prop_startingPoint_GP_y.push_back(pos_startingPoint_GP.y()); data_.prop_startingPoint_GP_z.push_back(pos_startingPoint_GP.z());
    data_.prop_yroll.push_back(tmp_prop_LP.y());
    LocalPoint local_to_center(tmp_prop_LP.x(), tmp_prop_LP.y() + prop_y_to_center, 0);
    float local_phi = local_to_center.phi();
    data_.prop_localphi_rad.push_back((3.14159265/2.) - local_phi);
    data_.prop_localphi_deg.push_back(((3.14159265/2.) - local_phi)*(180./3.14159265));
    data_.has_prop.push_back(tmp_has_prop);
    data_.has_fidcut.push_back(tmp_has_fidcut);
    data_.prop_region.push_back(ch->id().region()); data_.prop_station.push_back(ch->id().station()); data_.prop_chamber.push_back(ch->id().chamber()); data_.prop_layer.push_back(ch->id().layer()); data_.prop_roll.push_back(ch->id().roll());
  }
}
void analyser::GEM_rechit_matcher(const GEMEtaPartition* ch, LocalPoint prop_LP, MuonData& data_){
  float tmp_rechit_GP_x; float tmp_rechit_GP_y; float tmp_rechit_GP_z;
  float tmp_rechit_LP_x; float tmp_rechit_LP_y; float tmp_rechit_LP_z;
  float tmp_rechit_yroll; float tmp_rechit_localphi_rad; float tmp_rechit_localphi_deg;
  bool tmp_has_rechit = false;
  int tmp_rechit_first_strip; int tmp_rechit_CLS; int tmp_rechit_BunchX;
  float tmp_RdPhi = 9999.; float tmp_RdPhi_Corrected; int tmp_rechit_detId;
  int tmp_nRecHitsTot = 0; int tmp_nRecHits5 = 0; int tmp_nRecHits2 = 0;
  int tmp_rechit_region; int tmp_rechit_station; int tmp_rechit_chamber; int tmp_rechit_layer; int tmp_rechit_roll;
  for(auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
    if((hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
      GEMDetId gemid((hit)->geographicalId());
      if(gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
        const auto& etaPart = GEMGeometry_->etaPartition(gemid);
        float strip = etaPart->strip(hit->localPosition());
        float stripAngle = etaPart->specificTopology().stripAngle(strip);
        float rechit_y_to_center = etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2)).perp();
        float rechit_y_to_chamber = (GEMGeometry_->chamber(ch->id()))->toLocal(etaPart->toGlobal(etaPart->centreOfStrip(etaPart->nstrips()/2))).y();
        LocalPoint local_to_center((hit)->localPosition().x(), rechit_y_to_center + (hit)->localPosition().y(), 0);
        if (ch->id().station() == 1 and ch->id().ring() == 1 and fabs((hit)->localPosition().x() - prop_LP.x()) < 999.0){
          tmp_nRecHitsTot++;
          if(abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 5){tmp_nRecHits5++;}
          if(abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch)) < 2){tmp_nRecHits2++;}
          if(abs(tmp_RdPhi) > abs(RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch))){
            tmp_rechit_GP_x = etaPart->toGlobal((hit)->localPosition()).x(); tmp_rechit_GP_y = etaPart->toGlobal((hit)->localPosition()).y(); tmp_rechit_GP_z = etaPart->toGlobal((hit)->localPosition()).z();
            tmp_rechit_LP_x = (hit)->localPosition().x(); tmp_rechit_LP_y = rechit_y_to_chamber + (hit)->localPosition().y(); tmp_rechit_LP_z = (hit)->localPosition().z();
            tmp_rechit_yroll = (hit)->localPosition().y();
            float local_phi = local_to_center.phi();
            tmp_rechit_localphi_rad = (3.14159265/2.) - local_phi;
            tmp_rechit_localphi_deg = ((3.14159265/2.) - local_phi)*(180./3.14159265);
            tmp_has_rechit = true;
            tmp_rechit_first_strip = (hit)->firstClusterStrip();
            tmp_rechit_CLS = (hit)->clusterSize();
            tmp_rechit_BunchX = (hit)->BunchX();
            tmp_RdPhi = RdPhi_func(stripAngle, hit, prop_LP.x(), prop_LP.y(), ch);
            tmp_RdPhi_Corrected = tmp_RdPhi;
            if((gemid.region() == 1 and gemid.chamber()%2 == 1) || (gemid.region() == -1 && gemid.chamber()%2 == 0)){
              tmp_RdPhi_Corrected = -1.0*tmp_RdPhi_Corrected;
            }
            tmp_rechit_detId = gemid.region()*(gemid.station()*100 + gemid.chamber());
            tmp_rechit_region = gemid.region(); tmp_rechit_station = gemid.station(); tmp_rechit_chamber = gemid.chamber(); tmp_rechit_layer = gemid.layer(); tmp_rechit_roll = gemid.roll();
          }
        }
      }
    }
  }
  if(tmp_has_rechit){
    data_.rechit_GP_x.push_back(tmp_rechit_GP_x); data_.rechit_GP_y.push_back(tmp_rechit_GP_y); data_.rechit_GP_z.push_back(tmp_rechit_GP_z);
    data_.rechit_LP_x.push_back(tmp_rechit_LP_x); data_.rechit_LP_y.push_back(tmp_rechit_LP_y); data_.rechit_LP_z.push_back(tmp_rechit_LP_z);
    data_.rechit_yroll.push_back(tmp_rechit_yroll);
    data_.rechit_localphi_rad.push_back(tmp_rechit_localphi_rad);
    data_.rechit_localphi_deg.push_back(tmp_rechit_localphi_deg);
    data_.has_rechit.push_back(tmp_has_rechit);
    data_.rechit_first_strip.push_back(tmp_rechit_first_strip);
    data_.rechit_CLS.push_back(tmp_rechit_CLS);
    data_.rechit_BunchX.push_back(tmp_rechit_BunchX);
    data_.RdPhi.push_back(tmp_RdPhi);
    data_.RdPhi_Corrected.push_back(tmp_RdPhi_Corrected);
    data_.rechit_detId.push_back(tmp_rechit_detId);
    data_.nRecHitsTot.push_back(tmp_nRecHitsTot); data_.nRecHits5.push_back(tmp_nRecHits5); data_.nRecHits2.push_back(tmp_nRecHits2);
    data_.rechit_region.push_back(tmp_rechit_region); data_.rechit_station.push_back(tmp_rechit_station); data_.rechit_chamber.push_back(tmp_rechit_chamber); data_.rechit_layer.push_back(tmp_rechit_layer); data_.rechit_roll.push_back(tmp_rechit_roll);
  }
}
void analyser::propagate(const reco::Muon* mu, int prop_type, const edm::Event& iEvent, int i){
  const reco::Track* Track;
  reco::TransientTrack ttTrack;
  TTree* tree;
  if(prop_type == 1){ //If want to swith to global, use mu->globalTrack().get()
    tree = CSC_tree;
    Track = mu->outerTrack().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 2){
    tree = Tracker_tree;
    Track = mu->track().get();
    ttTrack = ttrackBuilder_->build(Track);
  }
  else if (prop_type == 3){
    tree = Segment_tree;
  }
  else{
    std::cout << "Bad prop type, failure." << std::endl; return;
  }
  data_.init();
  //Muon Info//////////////////////////////////////////////////////
  data_.muon_charge = mu->charge(); data_.muon_pt = mu->pt(); data_.muon_eta = mu->eta(); data_.muon_momentum = mu->momentum().mag2();
  data_.evtNum = iEvent.eventAuxiliary().event(); data_.lumiBlock = iEvent.eventAuxiliary().luminosityBlock(); data_.muonIdx = data_.evtNum*100 + i;
  //Track Info//////////////////////////////////////////////////////
  data_.track_chi2 = Track->chi2(); data_.track_ndof = Track->ndof();
  CSCSegmentCounter(mu, data_.n_ME11_segment, data_.nCSCSeg, data_.nDTSeg);
  if(data_.n_ME11_segment >= 1){data_.hasME11 = 1;}
  //which_track
  //Propagation Info//////////////////////////////////////////////////////
  for (const auto& ch : GEMGeometry_->etaPartitions()) {
    if (ch->id().station() != 1) continue; //Only takes GE1/1
    GlobalPoint tmp_pos_GP; bool tmp_has_prop = 0;
    propagate_to_GEM(ttTrack, ch, prop_type, tmp_has_prop, tmp_pos_GP, data_);
    if(tmp_has_prop){
      LocalPoint tmp_prop_LP = ch->toLocal(tmp_pos_GP);
      //Rechit Info//////////////////////////////////////////////////////
      GEM_rechit_matcher(ch, tmp_prop_LP, data_);


     
    }
  }
  tree->Fill();
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

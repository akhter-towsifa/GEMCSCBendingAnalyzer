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
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"

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
  //============ Rechit Info =============//
  float rechit_GP[3]; float rechit_LP[3];        bool has_rechit;
  float rechit_yroll; float rechit_localphi_rad; float rechit_localphi_deg;
  int rechit_first_strip;     int rechit_CLS;    int rechit_BunchX;
  float RdPhi;        float RdPhi_Corrected;     int rechit_detId;
  float bending_angle;
  int nRecHitsTot;    int nRecHits5;             int nRecHits2;
  int rechit_location[5];
  int nRecHitsRpos1L1; int nRecHitsRpos1L2;
  int nRecHitsRneg1L1; int nRecHitsRneg1L2;
  //=========== Sim info for MC ==========//
  float sim_GP[3];   float sim_LP[3];
  float simDy;       float sim_yroll;            int nSim;
};

void MuonData::init()
{
  //=========== Muon Info ===============//
  muon_charge = 9999; muon_pt = 99999; muon_eta = 9999; muon_momentum = 9999;
  evtNum = 99999999; lumiBlock = 99999999; muonIdx = 99999999; runNum = 99999999;
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
  for(int i=0; i<5; ++i){
    ME11_location[i] = 999999;
  }
  inner_or_outer_mom = 99999;
  //=========== Rechit Info ===========//
  for (int i=0; i<3; ++i){
    rechit_GP[i] = 999999; rechit_LP[i] = 999999;
  }
  rechit_yroll = 999999; rechit_localphi_rad = 999999; rechit_localphi_deg = 999999;
  has_rechit = false;
  rechit_first_strip = 999999; rechit_CLS = 999999; rechit_BunchX = 999999;
  RdPhi = 999999; RdPhi_Corrected = 999999; rechit_detId = 999999;
  bending_angle = 999999;
  nRecHitsTot = 999999; nRecHits5 = 999999; nRecHits2 = 999999;
  for (int i=0; i<5; ++i){
    rechit_location[i] = 999999;
  }
  nRecHitsRpos1L1 = 999999; nRecHitsRpos1L2 = 999999;
  nRecHitsRneg1L1 = 999999; nRecHitsRneg1L2 = 999999;
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
  else{
    std::cout << "Bad prop type, failure, doesnt fall under the 3 prop_type listed" << std::endl;
  }
  //=========== Muon Info =============//
  t->Branch("muon_charge", &muon_charge); t->Branch("muon_pt", &muon_pt);
  t->Branch("muon_eta", &muon_eta);       t->Branch("muon_momentum", &muon_momentum);
  t->Branch("evtNum", &evtNum);           t->Branch("lumiBlock", &lumiBlock);
  t->Branch("runNum", &runNum);           t->Branch("muonIdx", &muonIdx);
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
  t->Branch("prop_location", &prop_location, "prop_location[5] (reg, sta, cha, lay, rol)I");
  //========== Track Info =============//
  t->Branch("track_chi2", &track_chi2);         t->Branch("track_ndof", &track_ndof);
  t->Branch("n_ME11_segment", &n_ME11_segment); t->Branch("which_track", &which_track);
  t->Branch("hasME11", &hasME11);               t->Branch("hasME11RecHit", &hasME11RecHit);
  t->Branch("hasME11A", &hasME11A);             t->Branch("hasME11ARecHit", &hasME11ARecHit);
  t->Branch("nCSCSeg", &nCSCSeg);               t->Branch("nDTSeg", &nDTSeg);
  t->Branch("nME11RecHits", &nME11RecHits);     t->Branch("ME11_BunchX", &ME11_BunchX);
  t->Branch("ME11_strip", &ME11_strip);
  t->Branch("ME11_location", &ME11_location, "ME11_location[5] (end, sta, ring, cha, lay)/I");
  t->Branch("inner_or_outer_mom", &inner_or_outer_mom, "inner_or_outer_mom (0 = inner, 1 = outer)/I");
  //========== Rechit Info ============//
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

  //add lines 230 - 236 from Devin's code afterwards here

  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::Handle<GEMRecHitCollection> gemRecHits;
  edm::EDGetTokenT<vector<PSimHit> > gemSimHits_;
  edm::Handle<vector<PSimHit> > gemSimHits;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<CSCSegmentCollection> cscSegments_;

  edm::Service<TFileService> fs;
  MuonServiceProxy* theService_;
  
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  
  edm::ESHandle<GEMGeometry> GEMGeometry_;
  edm::ESHandle<CSCGeometry> CSCGeometry_;

  bool CSC_prop; bool tracker_prop; bool Segment_prop;
  vector<int> prop_list;
  bool debug;
  bool isCosmic;

  MuonData data_;
  TTree* CSC_tree; TTree* Tracker_tree; TTree* Segment_tree;
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
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  gemSimHits_ = consumes<vector<PSimHit> >(iConfig.getParameter<edm::InputTag>("gemSimHits"));
  cscSegments_ = consumes<CSCSegmentCollection>(edm::InputTag("cscSegments"));

  tracker_prop = iConfig.getParameter<bool>("tracker_prop");
  CSC_prop = iConfig.getParameter<bool>("CSC_prop");
  Segment_prop = iConfig.getParameter<bool>("Segment_prop");
  debug = iConfig.getParameter<bool>("debug");
  isCosmic = iConfig.getParameter<bool>("isCosmic");
  std::cout << "tracker_prop " << tracker_prop << "\tCSC_prop " << CSC_prop << "\tSegment_prop " << Segment_prop << "\tdebug " << debug << std::endl;

  if(CSC_prop){CSC_tree = data_.book(CSC_tree, 1); prop_list.push_back(1);}
  if(tracker_prop){Tracker_tree = data_.book(Tracker_tree, 2); prop_list.push_back(2);}
  if(Segment_prop){Segment_tree = data_.book(Segment_tree, 3); prop_list.push_back(3);}
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

  edm::Handle<CSCSegmentCollection> cscSegments;
  if (! iEvent.getByToken(cscSegments_, cscSegments)){std::cout << "Bad segments" << std::endl;}

  if (debug) cout << "New! EventNumber = " << iEvent.eventAuxiliary().event() << " LumiBlock = " << iEvent.eventAuxiliary().luminosityBlock() << " RunNumber = " << iEvent.run() << endl;

  for (size_t i = 0; i < muons->size(); ++i){
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();
    if (not mu->isGlobalMuon()) continue;
    if (debug) cout << "new muon" << endl;
    for (auto it = std::begin(prop_list); it != std::end(prop_list); ++it){
      if (debug) std::cout << "prop " << *it << "about to start propagate" << std::endl;
      int prop_type = *it;
      //propagate(mu, prop_type, iEvent, i);
    }
  }
}

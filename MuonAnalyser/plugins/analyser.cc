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


using namespace std;
using namespace edm;


struct MuonData
{
  void init();
  TTree* book(TTree *t);

  int muon_charge;
  float muon_pt;

  //Prop from tracker
  bool has_prop_GE11;
  int prop_region_GE11;
  int prop_station_GE11;
  int prop_layer_GE11;
  int prop_chamber_GE11;
  int prop_roll_GE11;

/*
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
*/

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


  bool has_rechit_GE11;
  int rechit_region_GE11;
  int rechit_station_GE11;
  int rechit_layer_GE11;
  int rechit_chamber_GE11;
  int rechit_roll_GE11;
  float rechit_x_GE11;
  float rechit_y_GE11;
  float rechit_r_GE11;
  float rechit_localx_GE11;
  float rechit_localy_GE11;
  float rechit_y_adjusted_GE11;
  float rechit_localphi_rad_GE11;
  float rechit_localphi_deg_GE11;


  bool has_rechit_GE11_looseEta;
  int rechit_region_GE11_looseEta;
  int rechit_station_GE11_looseEta;
  int rechit_layer_GE11_looseEta;
  int rechit_chamber_GE11_looseEta;
  int rechit_roll_GE11_looseEta;
  float rechit_x_GE11_looseEta;
  float rechit_y_GE11_looseEta;
  float rechit_r_GE11_looseEta;
  float rechit_localx_GE11_looseEta;
  float rechit_localy_GE11_looseEta;
  float rechit_y_adjusted_GE11_looseEta;
  float rechit_localphi_rad_GE11_looseEta;
  float rechit_localphi_deg_GE11_looseEta;


  //float RdPhi_inner_GE11;
  float RdPhi_CSC_GE11;
  int det_id;

  float RdPhi_CSC_GE11_looseEta;
  int det_id_looseEta;

  //bool has_fidcut_inner_GE11;
  bool has_fidcut_CSC_GE11;
  int isGEMmuon;

  int which_track_CSC_GE11;
  //int which_track_inner_GE11;
  int hasME11;

  unsigned long long  evtNum;
};

void MuonData::init()
{
  muon_charge = 9999;
  muon_pt = 9999;

  has_prop_GE11 = false;
  prop_region_GE11 = 99999;
  prop_station_GE11 = 99999;
  prop_layer_GE11 = 99999;
  prop_chamber_GE11 = 99999;
  prop_roll_GE11 = 99999;

/*
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
*/

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

  has_rechit_GE11 = false;
  rechit_region_GE11 = 999999;
  rechit_station_GE11 = 999999;
  rechit_layer_GE11 = 999999;
  rechit_chamber_GE11 = 999999;
  rechit_roll_GE11 = 999999;
  rechit_x_GE11 = 999999;
  rechit_y_GE11 = 999999;
  rechit_r_GE11 = 999999;
  rechit_localx_GE11 = 999999;
  rechit_localy_GE11 = 999999;
  rechit_y_adjusted_GE11 = 999999;
  rechit_localphi_rad_GE11 = 999999;
  rechit_localphi_deg_GE11 = 999999;


  has_rechit_GE11_looseEta = false;
  rechit_region_GE11_looseEta = 999999;
  rechit_station_GE11_looseEta = 999999;
  rechit_layer_GE11_looseEta = 999999;
  rechit_chamber_GE11_looseEta = 999999;
  rechit_roll_GE11_looseEta = 999999;
  rechit_x_GE11_looseEta = 999999;
  rechit_y_GE11_looseEta = 999999;
  rechit_r_GE11_looseEta = 999999;
  rechit_localx_GE11_looseEta = 999999;
  rechit_localy_GE11_looseEta = 999999;
  rechit_y_adjusted_GE11_looseEta = 999999;
  rechit_localphi_rad_GE11_looseEta = 999999;
  rechit_localphi_deg_GE11_looseEta = 999999;


  //RdPhi_inner_GE11 = 999999;
  RdPhi_CSC_GE11 = 999999;
  det_id = 999999;

  RdPhi_CSC_GE11_looseEta = 999999;
  det_id_looseEta = 999999;

  //has_fidcut_inner_GE11 = false;
  has_fidcut_CSC_GE11 = false;
  isGEMmuon = 0;

  which_track_CSC_GE11 = 999;
  //which_track_inner_GE11 = 999;
  hasME11 = 0;
 
  evtNum = 0;
}

TTree* MuonData::book(TTree *t){
  edm::Service< TFileService > fs;
  t = fs->make<TTree>("MuonData", "MuonData");

  t->Branch("muon_charge", &muon_charge);
  t->Branch("muon_pt", &muon_pt);
//Propogated Inner
  t->Branch("has_prop_GE11", &has_prop_GE11);
  t->Branch("prop_region_GE11", &prop_region_GE11);
  t->Branch("prop_station_GE11", &prop_station_GE11);
  t->Branch("prop_layer_GE11", &prop_layer_GE11);
  t->Branch("prop_chamber_GE11", &prop_chamber_GE11);
  t->Branch("prop_roll_GE11", &prop_roll_GE11);

/*
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
*/

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
//Reconstructed
  t->Branch("has_rechit_GE11", &has_rechit_GE11);
  t->Branch("rechit_region_GE11", &rechit_region_GE11);
  t->Branch("rechit_station_GE11", &rechit_station_GE11);
  t->Branch("rechit_layer_GE11", &rechit_layer_GE11);
  t->Branch("rechit_chamber_GE11", &rechit_chamber_GE11);
  t->Branch("rechit_roll_GE11", &rechit_roll_GE11);
  t->Branch("rechit_x_GE11", &rechit_x_GE11);
  t->Branch("rechit_y_GE11", &rechit_y_GE11);
  t->Branch("rechit_r_GE11", &rechit_r_GE11);
  t->Branch("rechit_localx_GE11", &rechit_localx_GE11);
  t->Branch("rechit_localy_GE11", &rechit_localy_GE11);
  t->Branch("rechit_y_adjusted_GE11", &rechit_y_adjusted_GE11);
  t->Branch("rechit_localphi_rad_GE11", &rechit_localphi_rad_GE11);
  t->Branch("rechit_localphi_deg_GE11", &rechit_localphi_deg_GE11);
//Reconstructed
  t->Branch("has_rechit_GE11_looseEta", &has_rechit_GE11_looseEta);
  t->Branch("rechit_region_GE11_looseEta", &rechit_region_GE11_looseEta);
  t->Branch("rechit_station_GE11_looseEta", &rechit_station_GE11_looseEta);
  t->Branch("rechit_layer_GE11_looseEta", &rechit_layer_GE11_looseEta);
  t->Branch("rechit_chamber_GE11_looseEta", &rechit_chamber_GE11_looseEta);
  t->Branch("rechit_roll_GE11_looseEta", &rechit_roll_GE11_looseEta);
  t->Branch("rechit_x_GE11_looseEta", &rechit_x_GE11_looseEta);
  t->Branch("rechit_y_GE11_looseEta", &rechit_y_GE11_looseEta);
  t->Branch("rechit_r_GE11_looseEta", &rechit_r_GE11_looseEta);
  t->Branch("rechit_localx_GE11_looseEta", &rechit_localx_GE11_looseEta);
  t->Branch("rechit_localy_GE11_looseEta", &rechit_localy_GE11_looseEta);
  t->Branch("rechit_y_adjusted_GE11_looseEta", &rechit_y_adjusted_GE11_looseEta);
  t->Branch("rechit_localphi_rad_GE11_looseEta", &rechit_localphi_rad_GE11_looseEta);
  t->Branch("rechit_localphi_deg_GE11_looseEta", &rechit_localphi_deg_GE11_looseEta);
//Residual
  //t->Branch("RdPhi_inner_GE11", &RdPhi_inner_GE11);
  t->Branch("RdPhi_CSC_GE11", &RdPhi_CSC_GE11);
  t->Branch("det_id", &det_id);
//Residual loose eta
  t->Branch("RdPhi_CSC_GE11_looseEta", &RdPhi_CSC_GE11_looseEta);
  t->Branch("det_id_looseEta", &det_id_looseEta);
//Cut
  //t->Branch("has_fidcut_inner_GE11", &has_fidcut_inner_GE11);
  t->Branch("has_fidcut_CSC_GE11", &has_fidcut_CSC_GE11);
  t->Branch("isGEMmuon", &isGEMmuon);

  t->Branch("which_track_CSC_GE11", &which_track_CSC_GE11);
  //t->Branch("which_track_inner_GE11", &which_track_inner_GE11);
 
  t->Branch("hasME11", &hasME11);

  t->Branch("evtNum", &evtNum);
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
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;

  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  edm::ESHandle<GEMGeometry> GEMGeometry_;

  TTree * tree_data_;
  MuonData data_;

};

analyser::analyser(const edm::ParameterSet& iConfig)
{
  cout << "Begin analyser" << endl;
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector());

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));

  tree_data_ = data_.book(tree_data_);
}


void
analyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<View<reco::Muon> > muons;
  if (! iEvent.getByToken(muons_, muons)) return;

  if (muons->size() == 0) return;


  for (size_t i = 0; i < muons->size(); ++i){
    cout << "new muon" << endl;
    edm::RefToBase<reco::Muon> muRef = muons->refAt(i);
    const reco::Muon* mu = muRef.get();

    if (mu->pt() < 2.0) continue;
    if (not mu->standAloneMuon()) continue;
    cout << "is standalone" << endl;
    data_.init();
    data_.isGEMmuon = mu->isGEMMuon();

    auto matches = mu->matches();
    for ( auto MCM : matches){
      if (MCM.detector() != 2) continue;
      for( auto MSM : MCM.segmentMatches){
        auto cscSeg = MSM.cscSegmentRef;
        auto cscDetID = cscSeg->cscDetId();
        if (cscDetID.station() == 1 and cscDetID.ring() == 1) data_.hasME11 = 1;
      }
    }
    data_.evtNum = iEvent.eventAuxiliary().event();
/*
    if (not mu->innerTrack()) continue;
    const reco::Track* innerTrack = mu->track().get();
    const reco::Track* muonTrack = 0;
    const reco::Track* muonOuterTrack = 0;
    if ( mu->globalTrack().isNonnull() ){
      muonTrack = mu->globalTrack().get(); 
      muonOuterTrack = mu->outerTrack().get();
      }
    else 
      continue;
*/

    //Get det info
    //int chamber_count = 0;
    //for (const auto& ch : GEMGeometry_->chambers()){
    //  chamber_count ++;
    //}
    //cout << "chamber count = " << chamber_count << endl;

    const reco::Track* muonTrack = 0;
    if ( mu->outerTrack().isNonnull() ){
      muonTrack = mu->outerTrack().get();
      }
    //reco::TransientTrack ttTrack_tracker = ttrackBuilder_->build(innerTrack); //tracker to GEM

    reco::TransientTrack ttTrack_CSC = ttrackBuilder_->build(muonTrack); //CSC to GEM

    float count = 0;

    for (const auto& ch : GEMGeometry_->etaPartitions()) {
      //eta_count ++;
      //cout << "eta count = " << eta_count << endl;
      if (ch->id().station() != 1) continue; //Only takes GE1/1
      const BoundPlane& bps(ch->surface());

/*
      // Tracker propagated
      TrajectoryStateOnSurface tsos_inner;
      if ( innerTrack->outerPosition().Mag2() - innerTrack->innerPosition().Mag2() > 0){
        tsos_inner = propagator->propagate(ttTrack_tracker.outermostMeasurementState(),ch->surface());
        data_.which_track_inner_GE11 = "outer most";
      }
      else{
        tsos_inner = propagator->propagate(ttTrack_tracker.innermostMeasurementState(),ch->surface());
        data_.which_track_inner_GE11 = "inner most";
      }
      if (!tsos_inner.isValid()) continue;
      GlobalPoint pos_global_inner = tsos_inner.globalPosition();
      LocalPoint pos_local_inner = ch->toLocal(tsos_inner.globalPosition());
      const GlobalPoint pos2D_global_inner(pos_global_inner.x(), pos_global_inner.y(), 0);
      const LocalPoint pos2D_local_inner(pos_local_inner.x(), pos_local_inner.y(), 0);
      if (pos_global_inner.eta() * mu->eta() < 0.0) continue;
*/

      // CSC propagated
      TrajectoryStateOnSurface tsos_CSC;
      if ( muonTrack->outerPosition().Mag2() - muonTrack->innerPosition().Mag2() > 0){
        tsos_CSC = propagator->propagate(ttTrack_CSC.outermostMeasurementState(),ch->surface());
        data_.which_track_CSC_GE11 = 0;
      }
      else{
        tsos_CSC = propagator->propagate(ttTrack_CSC.innermostMeasurementState(),ch->surface());
        data_.which_track_CSC_GE11 = 1;
      }
      if (!tsos_CSC.isValid()) continue;
      GlobalPoint pos_global_CSC = tsos_CSC.globalPosition();
      LocalPoint pos_local_CSC = ch->toLocal(tsos_CSC.globalPosition());
      const GlobalPoint pos2D_global_CSC(pos_global_CSC.x(), pos_global_CSC.y(), 0);
      const LocalPoint pos2D_local_CSC(pos_local_CSC.x(), pos_local_CSC.y(), 0);
      //if (pos_global_CSC.eta() * mu->eta() < 0.0) continue;
     

      //if (bps.bounds().inside(pos2D_local_inner) and bps.bounds().inside(pos2D_local_CSC) and ch->id().station() == 1 and ch->id().ring() == 1){
      if (bps.bounds().inside(pos2D_local_CSC) and ch->id().station() == 1 and ch->id().ring() == 1){
      //if (ch->id().station() == 1 and ch->id().ring() == 1){
        //cout << "Track used was " << data_.which_track_CSC_GE11 << endl;
        cout << "pos_local_CSC = " << pos_local_CSC << endl;
        //cout << "bps.bounds().inside(pos2D_local_CSC) = " << bps.bounds().inside(pos2D_local_CSC) << endl;

        //cout << "charge is " << mu->charge() << endl;
        data_.muon_charge = mu->charge();
        data_.muon_pt = mu->pt();

        //data_.prop_inner_chi2_GE11 = innerTrack->chi2();
        //data_.prop_inner_ndof_GE11 = innerTrack->ndof();

        data_.prop_CSC_chi2_GE11 = muonTrack->chi2();
        data_.prop_CSC_ndof_GE11 = muonTrack->ndof();


        const float fidcut_angle = 1.0;
        const float cut_ang = 5.0 - fidcut_angle;
        const float cut_chamber = 5.0;
        //const float fidcut_y = 5.0;
        //const float cut_even_high = 250.0 - fidcut_y;
        //const float cut_odd_high = 250.0 - fidcut_y;
        //const float cut_even_low = 130.0 + fidcut_y;
        //const float cut_odd_low = 130.0 + fidcut_y;

        const auto& etaPart_ch = GEMGeometry_->etaPartition(ch->id());
        //float strip = etaPart_ch->strip(pos_local);

        // Tracker prop
        const float prop_y_to_center = etaPart_ch->toGlobal(etaPart_ch->centreOfStrip(etaPart_ch->nstrips()/2)).perp();
        //LocalPoint local_to_center_inner(pos_local_inner.x(), prop_y_to_center + pos_local_inner.y(), 0);
        //const float prop_inner_localphi_rad = (3.14159265/2.) - local_to_center_inner.phi();
        //const float prop_inner_localphi_deg = prop_inner_localphi_rad*180/3.14169265;
        count++;


        data_.has_prop_GE11 = true;
        data_.prop_region_GE11 = ch->id().region();
        data_.prop_station_GE11 = ch->id().station();
        data_.prop_layer_GE11 = ch->id().layer();
        data_.prop_chamber_GE11 = ch->id().chamber();
        data_.prop_roll_GE11 = ch->id().roll();
/*
        data_.prop_inner_x_GE11 = pos_global_inner.x();
        data_.prop_inner_y_GE11 = pos_global_inner.y();
        data_.prop_inner_r_GE11 = pos_global_inner.mag();
        data_.prop_inner_localx_GE11 = pos_local_inner.x();
        data_.prop_inner_localy_GE11 = pos_local_inner.y();
        data_.prop_inner_y_adjusted_GE11 = prop_y_to_center + pos_local_inner.y();
        data_.prop_inner_localphi_rad_GE11 = prop_inner_localphi_rad;
        data_.prop_inner_localphi_deg_GE11 = prop_inner_localphi_deg;

        if (ch->id().chamber()%2 == 0){
          if (prop_inner_localphi_deg > -cut_ang && prop_inner_localphi_deg < cut_ang && prop_y_to_center + pos_local_inner.y() > cut_low && prop_y_to_center + pos_local_inner.y() < cut_even_high){
            data_.has_fidcut_inner_GE11 = true;
          }
        }
        if (ch->id().chamber()%2 == 1){
          if (prop_inner_localphi_deg > -cut_ang && prop_inner_localphi_deg < cut_ang && prop_y_to_center + pos_local_inner.y() > cut_low && prop_y_to_center + pos_local_inner.y() < cut_odd_high){
            data_.has_fidcut_inner_GE11 = true;
          }
        }
*/
        // CSC prop
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






        //if (ch->id().chamber()%2 == 0){
        //  if (prop_CSC_localphi_deg > -cut_ang && prop_CSC_localphi_deg < cut_ang && prop_y_to_center + pos_local_CSC.y() > cut_even_low && prop_y_to_center + pos_local_CSC.y() < cut_even_high){
        //    data_.has_fidcut_CSC_GE11 = true;
        //  }
        //}
        //if (ch->id().chamber()%2 == 1){
        //  if (ch->id().roll() == 1 && ){
        //    if (prop_CSC_localphi_deg > -cut_ang && prop_CSC_localphi_deg < cut_ang && prop_y_to_center + pos_local_CSC.y() > cut_odd_low && prop_y_to_center + pos_local_CSC.y() < cut_odd_high){
        //      data_.has_fidcut_CSC_GE11 = true;
        //    }
        //  }
        //}




        
        for (auto hit = gemRecHits->begin(); hit != gemRecHits->end(); hit++){
          if ( (hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::GEM){
            GEMDetId gemid((hit)->geographicalId());
            if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
              cout << "starting rechit" << endl;
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


                if (abs(data_.RdPhi_CSC_GE11) > abs(cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll))){
                  std::cout << "Overwrite" << std::endl;


                  data_.has_rechit_GE11 = true;
                  data_.rechit_region_GE11 = gemid.region();
                  data_.rechit_station_GE11 = gemid.station();
                  data_.rechit_layer_GE11 = gemid.layer();
                  data_.rechit_chamber_GE11 = gemid.chamber();
                  data_.rechit_roll_GE11 = gemid.roll();
                  data_.rechit_x_GE11 = etaPart->toGlobal((hit)->localPosition()).x();
                  data_.rechit_y_GE11 = etaPart->toGlobal((hit)->localPosition()).y();
                  data_.rechit_r_GE11 = etaPart->toGlobal((hit)->localPosition()).mag();
                  data_.rechit_localx_GE11 = (hit)->localPosition().x();
                  data_.rechit_localy_GE11 = (hit)->localPosition().y();
                  data_.rechit_y_adjusted_GE11 = rechit_y_to_center + (hit)->localPosition().y();
                  data_.rechit_localphi_rad_GE11 = rechit_localphi_rad;
                  data_.rechit_localphi_deg_GE11 = rechit_localphi_deg;
                  //data_.RdPhi_inner_GE11 = cosAngle * (pos_local_inner.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_inner.y() + deltay_roll);
                  data_.RdPhi_CSC_GE11 = cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll);
                  data_.det_id = gemid.region()*(gemid.station()*100 + gemid.chamber());
                }
              }
            }





            if (gemid.station() == ch->id().station() and gemid.chamber() == ch->id().chamber() and gemid.layer() == ch->id().layer() and abs(gemid.roll() - ch->id().roll()) <= 1 and gemid.region() == ch->id().region()){
              cout << "starting rechit" << endl;
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


                if (abs(data_.RdPhi_CSC_GE11_looseEta) > abs(cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll))){
                  std::cout << "Overwrite" << std::endl;


                  data_.has_rechit_GE11_looseEta = true;
                  data_.rechit_region_GE11_looseEta = gemid.region();
                  data_.rechit_station_GE11_looseEta = gemid.station();
                  data_.rechit_layer_GE11_looseEta = gemid.layer();
                  data_.rechit_chamber_GE11_looseEta = gemid.chamber();
                  data_.rechit_roll_GE11_looseEta = gemid.roll();
                  data_.rechit_x_GE11_looseEta = etaPart->toGlobal((hit)->localPosition()).x();
                  data_.rechit_y_GE11_looseEta = etaPart->toGlobal((hit)->localPosition()).y();
                  data_.rechit_r_GE11_looseEta = etaPart->toGlobal((hit)->localPosition()).mag();
                  data_.rechit_localx_GE11_looseEta = (hit)->localPosition().x();
                  data_.rechit_localy_GE11_looseEta = (hit)->localPosition().y();
                  data_.rechit_y_adjusted_GE11_looseEta = rechit_y_to_center + (hit)->localPosition().y();
                  data_.rechit_localphi_rad_GE11_looseEta = rechit_localphi_rad;
                  data_.rechit_localphi_deg_GE11_looseEta = rechit_localphi_deg;
                  data_.RdPhi_CSC_GE11_looseEta = cosAngle * (pos_local_CSC.x() - (hit)->localPosition().x()) + sinAngle * (pos_local_CSC.y() + deltay_roll);
                  data_.det_id_looseEta = gemid.region()*(gemid.station()*100 + gemid.chamber());
                }
              }
            }
          }
        } 
      }
    }
    cout << "Filling!" << endl;
    tree_data_->Fill();
  }
}

void analyser::beginJob(){}
void analyser::endJob(){}

DEFINE_FWK_MODULE(analyser);

import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('analyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('TrackingTools.TrackRefitter.globalMuonTrajectories_cff')
process.load('TrackingTools.TrackFitters.TrackFitters_cff')

from Configuration.AlCa.GlobalTag import GlobalTag


### This is the misalignment part
misalign = True
do_GEM = True
do_CSC = False
if misalign:
  #db_file = 'sqlite_file:dummy_dx1.db'
  gem_db_file = 'sqlite_file:2022D_backPropModiefiedRefitTracker_alcareco_v0.db' #for GEM
  #csc_db_file = 'sqlite_file:Run3v1.db' #for csc alignment only in this case
  #gpr_db_file = 'sqlite_file:Run3v1.db' #for gpr only in this case
  process.GlobalTag.toGet = cms.VPSet(
    #GE11 rec/tag
    cms.PSet(
        connect = cms.string(gem_db_file),
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEMAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string(gem_db_file),
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('GEMAlignmentErrorExtendedRcd')
    ),
    #ME11 rec/tag
    cms.PSet(
        connect = cms.string(csc_db_file),
        record = cms.string('CSCAlignmentRcd'),
        tag = cms.string('CSCAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string(csc_db_file),
        record = cms.string('CSCAlignmentErrorExtendedRcd'),
        tag = cms.string('CSCAlignmentErrorExtendedRcd')
    ),
    cms.PSet(
        connect = cms.string(gpr_db_file), 
        record = cms.string('GlobalPositionRcd'), 
        tag = cms.string('GlobalPositionRcd') #cms.string('IdealGeometry')
    )
  )

  process.GEMGeometryESModule.applyAlignment = cms.bool(do_GEM)
  #process.CSCGeometryESModule.applyAlignment = cms.bool(do_CSC)
################################




#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_design', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v14', '')

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(options.nEvents)
)
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
				fileNames = cms.untracked.vstring(options.inputFiles), 
				inputCommands = cms.untracked.vstring(
			"keep *", 
			"drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			"drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			)
				)

#testfile = "/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_83.root"
outfile = "out_ge11.root"
#process.source.fileNames.append('file:'+testfile)
process.source.fileNames.append('root://cms-xrd-global.cern.ch//store/data/Run2022D/Muon/ALCARECO/MuAlCalIsolatedMu-PromptReco-v2/000/357/734/00000/20e6e175-9a53-4d4a-b233-8cb4fae82b0b.root')

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) 

from TrackingTools.TrackRefitter.globalMuonTrajectories_cff import *
process.MuonAlignmentFromReferenceGlobalMuonRefit = globalMuons.clone()
process.MuonAlignmentFromReferenceGlobalMuonRefit.Tracks = cms.InputTag("ALCARECOMuAlCalIsolatedMu:TrackerOnly")
process.MuonAlignmentFromReferenceGlobalMuonRefit.TrackTransformer.RefitRPCHits = cms.bool(False)



process.analyzer = cms.EDAnalyzer('analyzer', 
	process.MuonServiceProxy,
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("ALCARECOMuAlCalIsolatedMu:SelectedMuons"),
        ref_track = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted"),
	vertexCollection = cms.InputTag("offlinePrimaryVertices"),
        tracker_prop = cms.bool(True),
        CSC_prop = cms.bool(False),
        Segment_prop = cms.bool(True),
        trackerRefit_prop = cms.bool(True),
        debug = cms.bool(False),
        isCosmic = cms.bool(False)
)

process.p = cms.Path(process.MuonAlignmentFromReferenceGlobalMuonRefit + process.analyzer)

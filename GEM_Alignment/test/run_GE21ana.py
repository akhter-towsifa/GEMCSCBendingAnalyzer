import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('analyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('TrackingTools.TrackRefitter.globalMuonTrajectories_cff')
process.load('TrackingTools.TrackFitters.TrackFitters_cff')
process.load('RecoLocalMuon.CSCSegment.cscSegments_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag

### This is the misalignment part
misalign = True
do_GEM = False
do_CSC = False
if misalign:
  #db_file = 'sqlite_file:dummy_dx1.db'
  # gem_db_file = 'sqlite_file:gem.db' #for GEM
  # csc_db_file = 'sqlite_file:csc.db' #for csc alignment only in this case
  # gpr_db_file = 'sqlite_file:gpr.db' #for gpr only in this case
  process.GlobalTag.toGet = cms.VPSet(
    #GE11 rec/tag
    # cms.PSet(
    #     connect = cms.string(gem_db_file),
    #     record = cms.string('GEMAlignmentRcd'),
    #     tag = cms.string('GEMAlignmentRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(gem_db_file),
    #     record = cms.string('GEMAlignmentErrorExtendedRcd'),
    #     tag = cms.string('GEMAlignmentErrorExtendedRcd')
    # ),
    #ME11 rec/tag
    # cms.PSet(
    #     connect = cms.string(csc_db_file),
    #     record = cms.string('CSCAlignmentRcd'),
    #     tag = cms.string('CSCAlignmentRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(csc_db_file),
    #     record = cms.string('CSCAlignmentErrorExtendedRcd'),
    #     tag = cms.string('CSCAlignmentErrorExtendedRcd')
    # ),
    # cms.PSet(
    #     connect = cms.string(gpr_db_file), 
    #     record = cms.string('GlobalPositionRcd'), 
    #     tag = cms.string('GlobalPositionRcd') #cms.string('IdealGeometry')
    # )
  )

  process.GEMGeometryESModule.applyAlignment = cms.bool(do_GEM)
  process.CSCGeometryESModule.applyAlignment = cms.bool(do_CSC)
################################


#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_prompt', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_Candidate_2024_10_08_09_43_55', '')


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
#process.source.fileNames.append('file:'+testfile)
outfile = "out_ge21.root"

#process.source.fileNames.append('root://cms-xrd-global.cern.ch/')
process.source.fileNames.append('root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/ALCARECO/MuAlCalIsolatedMu-PromptReco-v1/000/385/836/00000/26c36d10-e695-4cec-a930-67303655004d.root')

process.options = cms.untracked.PSet(
                        #SkipEvent = cms.untracked.vstring('ProductNotFound')
                        TryToContinue = cms.untracked.vstring('ProductNotFound') #SkipEvent parameter does not work for CMSSW_13_3_X and above
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) 

from TrackingTools.TrackRefitter.globalMuonTrajectories_cff import *
process.MuonAlignmentFromReferenceGlobalMuonRefit = globalMuons.clone()
process.MuonAlignmentFromReferenceGlobalMuonRefit.Tracks = cms.InputTag("ALCARECOMuAlCalIsolatedMu:TrackerOnly")
process.MuonAlignmentFromReferenceGlobalMuonRefit.TrackTransformer.RefitRPCHits = cms.bool(False)

from RecoLocalMuon.CSCSegment.cscSegments_cfi import *
process.cscSegments = cscSegments.clone()

process.analyzer = cms.EDAnalyzer('ge21analyzer', 
	process.MuonServiceProxy,
        cscSegmentsReco = cms.InputTag("cscSegments"),
	      gemRecHits = cms.InputTag("gemRecHits"), 
	      gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("ALCARECOMuAlCalIsolatedMu:SelectedMuons"),
        ref_track = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted"),
	      vertexCollection = cms.InputTag("offlinePrimaryVertices"),
        tracker_prop = cms.bool(False),
        CSC_prop = cms.bool(False),
        Segment_prop = cms.bool(False),
        trackerRefit_prop = cms.bool(False),
        SegmentReco_prop = cms.bool(True),
        debug = cms.bool(True),
        isCosmic = cms.bool(False)
)

process.p = cms.Path(process.MuonAlignmentFromReferenceGlobalMuonRefit + process.cscSegments + process.analyzer)

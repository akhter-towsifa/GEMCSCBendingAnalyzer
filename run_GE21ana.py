import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('analyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
#process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
# process.load('TrackingTools.TrackRefitter.globalMuonTrajectories_cff')
process.load('RecoMuon.GlobalMuonProducer.globalMuons_cfi')
process.load('TrackingTools.TrackFitters.TrackFitters_cff')
process.load('RecoLocalMuon.CSCSegment.cscSegments_cfi')
process.load('RecoLocalMuon.GEMRecHit.gemRecHits_cfi')
process.load('Geometry.GEMGeometryBuilder.gemGeometryDB_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag

### This is the misalignment part
misalign = True
do_GEM = False
do_CSC = False
if misalign:
  #db_file = 'sqlite_file:dummy_dx1.db'
  # gem_db_file = 'sqlite_file:/afs/cern.ch/work/t/toakhter/public/GEM_Alignment/2024H_prompt_reco_v0.db' #for GEM
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
    # )
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
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_dataRun3_v5', '')

### Try to do RecoLocalMuon on all muon detectors ###
#####################################################
from RecoLocalMuon.Configuration.RecoLocalMuon_cff import *
process.localreco = cms.Sequence(muonlocalreco)

# Skip Digi2Raw and Raw2Digi steps for All Muon detectors (skipping DT, RPC)#
##########################################################
process.gemRecHits.gemDigiLabel = cms.InputTag('muonGEMDigis', 'GEMDigi')
process.csc2DRecHits.wireDigiTag = cms.InputTag("muonCSCDigis","MuonCSCWireDigi")
process.csc2DRecHits.stripDigiTag = cms.InputTag("muonCSCDigis","MuonCSCStripDigi")

process.gemRecHits = cms.EDProducer("GEMRecHitProducer",
    recAlgoConfig = cms.PSet(),
    recAlgo = cms.string('GEMRecHitStandardAlgo'),
    gemDigiLabel = cms.InputTag("muonGEMDigis"),
    ge21Off = cms.bool(False),
)
# process.rechit_step  = cms.Path(process.gemRecHits)

process.MessageLogger.cerr.FwkReport.reportEvery = 5000

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register ('nEvents',
			-1, #Max number of events 
			VarParsing.multiplicity.singleton, 
			VarParsing.varType.int, 
			"Number of events")
options.parseArguments()

# process.maxEvents = cms.untracked.PSet(
#   input = cms.untracked.int32(options.nEvents)
# )
process.maxEvents.input = cms.untracked.int32(-1)


process.source = cms.Source("PoolSource", 
			fileNames = cms.untracked.vstring(options.inputFiles), 
			inputCommands = cms.untracked.vstring(
			  "keep *", 
			  "drop TotemTimingDigiedmDetSetVector_totemTimingRawToDigi_TotemTiming_reRECO", 
			  "drop TotemTimingRecHitedmDetSetVector_totemTimingRecHits__reRECO"
			)
      # SelectEvents = cms.untracked.PSet(
      #   SelectEvents = cms.vstring('rechit_step')
      # )
		)

#testfile = "/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_83.root"
#process.source.fileNames.append('file:'+testfile)
outfile = "out_ge21.root"

#process.source.fileNames.append('root://cms-xrd-global.cern.ch/')
#process.source.fileNames.append('root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/385/836/00000/1501c861-c37a-42e5-a1a0-48ee5dc1be01.root')
process.source.fileNames.append('file:/eos/user/t/toakhter/tamu_mual/2024/2024H/GE21_rechit_out_local_reco_test.root')

process.options = cms.untracked.PSet(
                        TryToContinue = cms.untracked.vstring('ProductNotFound') #SkipEvent parameter does not work for CMSSW_13_3_X and above
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) 

# from TrackingTools.TrackRefitter.globalMuonTrajectories_cff import *
# process.MuonAlignmentFromReferenceGlobalMuonRefit = globalMuons.clone()
# process.MuonAlignmentFromReferenceGlobalMuonRefit.Tracks = cms.InputTag('muons')#("ALCARECOMuAlCalIsolatedMu:TrackerOnly")
# process.MuonAlignmentFromReferenceGlobalMuonRefit.TrackTransformer.RefitRPCHits = cms.bool(False)

from RecoLocalMuon.CSCSegment.cscSegments_cfi import *
process.cscSegments = cscSegments.clone()

process.analyzer = cms.EDAnalyzer('ge21analyzer', 
	      process.MuonServiceProxy,
        cscSegmentsReco = cms.InputTag("cscSegments"),
	      gemRecHits = cms.InputTag("gemRecHits", "", "GEMLocalRECO"), 
	      gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("muons"),#("ALCARECOMuAlCalIsolatedMu:SelectedMuons"),
        # ref_track = cms.InputTag("MuonAlignmentFromReferenceGlobalMuonRefit:Refitted"),
	      vertexCollection = cms.InputTag("offlinePrimaryVertices"),
        tracker_prop = cms.bool(False),
        CSC_prop = cms.bool(False),
        Segment_prop = cms.bool(False),
        trackerRefit_prop = cms.bool(False),
        SegmentReco_prop = cms.bool(True),
        debug = cms.bool(False),
        isCosmic = cms.bool(False)
)

# process.p = cms.Path(process.MuonAlignmentFromReferenceGlobalMuonRefit + process.cscSegments + process.gemRecHits + process.analyzer)
process.p = cms.Path( process.gemRecHits +
                      process.cscSegments +
                      process.analyzer
                    )

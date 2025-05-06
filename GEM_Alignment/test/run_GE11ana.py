import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('analyzer',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

from Configuration.AlCa.GlobalTag import GlobalTag


### This is the misalignment part

misalign = True
do_GEM = True
do_CSC = True
if misalign:
#   #db_file = 'sqlite_file:dummy_dx1.db'
  gem_db_file = 'sqlite_file:myDB.db' #for GEM
  csc_db_file = 'sqlite_file:myDB.db' #for csc alignment only in this case
#   #gpr_db_file = 'sqlite_file:Run3v1.db' #for gpr only in this case
  process.GlobalTag.toGet = cms.VPSet(
    #GE11 rec/tag
    cms.PSet(
        connect = cms.string(gem_db_file),
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEMAlignment_prompt_v2')
    ),
    cms.PSet(
        connect = cms.string(gem_db_file),
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('GEMAlignmentErrorExtended_6x6_prompt_v2')
    ),
    #ME11 rec/tag
    cms.PSet(
        connect = cms.string(csc_db_file),
        record = cms.string('CSCAlignmentRcd'),
        tag = cms.string('CSCAlignment_2009_v2_express')
    ),
    cms.PSet(
        connect = cms.string(csc_db_file),
        record = cms.string('CSCAlignmentErrorExtendedRcd'),
        tag = cms.string('CSCAlignmentErrorExtended_6x6_express')
    )
#   #  cms.PSet(
#   #      connect = cms.string(gpr_db_file), 
#   #      record = cms.string('GlobalPositionRcd'), 
#   #      tag = cms.string('GlobalPositionRcd') #cms.string('IdealGeometry')
#   #  )
  )


  process.GEMGeometryESModule.applyAlignment = cms.bool(do_GEM)
  process.CSCGeometryESModule.applyAlignment = cms.bool(do_CSC)

################################

#process.GEMGeometryESModule.applyAlignment = cms.bool(True)
#process.CSCGeometryESModule.applyAlignment = cms.bool(True)


#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_design', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '150X_dataRun3_Prompt_v1', '')


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


outfile = "out_GE11ana.root"
process.source.fileNames.append('root://cms-xrd-global.cern.ch//store/data/Run2024H/Muon0/RAW-RECO/ZMu-PromptReco-v1/000/385/836/00000/21e6e5ae-42c4-4a0d-a06a-a27fc605fdf3.root')
# process.source.fileNames.append('file:/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_109.root')


process.options = cms.untracked.PSet(
                        TryToContinue = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) #variable name set above

process.analyzer = cms.EDAnalyzer('analyzer', 
	process.MuonServiceProxy,
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("muons"),
	      vertexCollection = cms.InputTag("offlinePrimaryVertices"),
        tracker_prop = cms.bool(True),
        CSC_prop = cms.bool(False),
        Segment_prop = cms.bool(True),
        debug = cms.bool(False), #set to False before submitting a crab job
        isCosmic = cms.bool(False)
)

process.p = cms.Path(process.analyzer)

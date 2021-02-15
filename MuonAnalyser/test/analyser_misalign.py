import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyser',Phase2C9)
process = cms.Process('analyser',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

from Configuration.AlCa.GlobalTag import GlobalTag

#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')



### This is the misalignment part
process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        connect = cms.string('sqlite_file:Alignments.db'),
        record = cms.string('GEMAlignmentRcd'),
        tag = cms.string('GEM')
    ),
    cms.PSet(
        connect = cms.string('sqlite_file:Alignments.db'),
        record = cms.string('GEMAlignmentErrorExtendedRcd'),
        tag = cms.string('test')
    ),
    cms.PSet(record=cms.string('GlobalPositionRcd'), tag = cms.string('IdealGeometry'))
)
################################



process.GEMGeometryESModule.applyAlignment = cms.bool(True)





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

#process.source.fileNames.append('file:A6D8F0A8-462B-0344-91D7-35EC85B5DB6A.root')
process.source.fileNames.append('file:MCsegfault_test3.root')

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string("out_MCsegfault_test3.root"))

process.analyser = cms.EDAnalyzer('analyser', 
	process.MuonServiceProxy, 
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("muons"),
	vertexCollection = cms.InputTag("offlinePrimaryVerticies")
)

#process.p = cms.EndPath(process.analyser)
process.p = cms.Path(process.analyser)

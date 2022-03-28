import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
from Configuration.Eras.Era_Run3_cff import Run3

#process = cms.Process('analyzer',Phase2C9)
process = cms.Process('ME11ana',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_0T_cff') #0T for cruzet runs

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('RecoMuon.TrackingTools.MuonServiceProxy_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')

from Configuration.AlCa.GlobalTag import GlobalTag


### This is the misalignment part
misalign = False
if misalign:
  process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(
        connect = cms.string('sqlite_file:ME11_UL2018IOV_iter1_fixedFromXML.csv.db'),
        record = cms.string('CSCAlignmentRcd'),
        tag = cms.string('CSCAlignmentRcd')
    ),
    cms.PSet(
        connect = cms.string('sqlite_file:ME11_UL2018IOV_iter1_fixedFromXML.csv.db'),
        record = cms.string('CSCAlignmentErrorExtendedRcd'),
        tag = cms.string('CSCAlignmentErrorExtendedRcd')
    ),
    cms.PSet(record=cms.string('GlobalPositionRcd'), tag = cms.string('IdealGeometry'))
  )

  #process.GEMGeometryESModule.applyAlignment = cms.bool(True)
  process.CSCGeometryESModule.applyAlignment = cms.bool(True)
################################




#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_design', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data_prompt', '') #Antonello Comparison

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

testfile = "/eos/cms/store/group/alca_muonalign/singleMuonGun_11_3_4_2021_design/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3/210816_170519/0000/step2_83.root"
outfile = "out_ME11ana_test.root"
process.source.fileNames.append('file:'+testfile)
#process.source.fileNames.append('file:step2.root')
#process.source.fileNames.append('root://cms-xrd-global.cern.ch//store/data/Commissioning2021/Cosmics/AOD/PromptReco-v1/000/339/479/00000/72dc8a4a-a3cd-4eba-9f8b-7b38d3c3ec38.root')
#process.source.fileNames.append('file:CRUZET_344064_testfile.root')
#process.source.fileNames.append('file:2018runCtest.root')

process.options = cms.untracked.PSet(
                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                        )

process.TFileService = cms.Service("TFileService", fileName = cms.string(outfile)) #variable name set above

process.ME11ana = cms.EDAnalyzer('ME11ana', 
	process.MuonServiceProxy,
        csc2DRecHits = cms.InputTag("csc2DRecHits"), 
	gemRecHits = cms.InputTag("gemRecHits"), 
	gemSimHits = cms.InputTag("g4SimHits", "MuonGEMHits"), 
        muons = cms.InputTag("muons"),
	vertexCollection = cms.InputTag("offlinePrimaryVerticies"),
        debug = cms.bool(False),
        isCosmic = cms.bool(False)
)

#process.p = cms.EndPath(process.analyzer)
process.p = cms.Path(process.ME11ana)
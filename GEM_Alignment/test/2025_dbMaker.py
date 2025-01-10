import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process("TEST", Run3)
# Message logger service
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '142X_mcRun3_2025_realistic_Candidate_2025_01_08_09_08_01', '')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.source = cms.Source("EmptySource")

import Geometry.DTGeometryBuilder.dtGeometryDB_cfi
process.DTGeometryMuonMisalignedProducer = Geometry.DTGeometryBuilder.dtGeometryDB_cfi.DTGeometryESModule.clone()
process.DTGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.DTGeometryMuonMisalignedProducer.applyAlignment = cms.bool(False)
import Geometry.CSCGeometryBuilder.cscGeometryDB_cfi
process.CSCGeometryMuonMisalignedProducer = Geometry.CSCGeometryBuilder.cscGeometryDB_cfi.CSCGeometryESModule.clone()
process.CSCGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.CSCGeometryMuonMisalignedProducer.applyAlignment = cms.bool(False) #by default, in the GT, this is True. if this is set to False, need to comment it out. set to True when doing alignment for CSC specificlly 
import Geometry.GEMGeometryBuilder.gemGeometryDB_cfi
process.GEMGeometryMuonMisalignedProducer = Geometry.GEMGeometryBuilder.gemGeometryDB_cfi.GEMGeometryESModule.clone()
process.GEMGeometryMuonMisalignedProducer.appendToDataLabel = 'idealForMuonMisalignedProducer'
process.GEMGeometryMuonMisalignedProducer.applyAlignment = cms.bool(False) #set to true when GEM alignment is needed. otherwise False during CSC alignment

process.GEMAlDBWriter = cms.EDAnalyzer("GEMAlDBWriter",
                                       doChamber = cms.untracked.bool(True),
                                       doEndcap = cms.untracked.bool(False),
                                       doME11Chamber = cms.untracked.bool(False),
                                       doCSCEndcap = cms.untracked.bool(False),
                                       chamberFile = cms.untracked.string('/afs/cern.ch/work/t/toakhter/public/GEM_Alignment/2024H_prompt_reco_v0.csv'),      # GEM Chamber Alignment csv
                                       endcapFile = cms.untracked.string(''),       # GEM Endcap Alignment csv
                                       ME11ChamberFile = cms.untracked.string(''),      # ME1/1 Chamber Alignment csv
                                       CSCEndcapFile = cms.untracked.string('')     # ME1/1 Endcap Alignment csv
                                      )

from CondCore.DBCommon.CondDBSetup_cfi import CondDBSetup

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    CondDBSetup,
    toPut = cms.VPSet(
        cms.PSet(
        record = cms.string('DTAlignmentRcd'),
        tag = cms.string('DTAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('DTAlignmentErrorExtendedRcd'),
            tag = cms.string('DTAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('CSCAlignmentRcd'),
            tag = cms.string('CSCAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('CSCAlignmentErrorExtendedRcd'),
            tag = cms.string('CSCAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('GEMAlignmentRcd'),
            tag = cms.string('GEMAlignmentRcd')
        ),
        cms.PSet(
            record = cms.string('GEMAlignmentErrorExtendedRcd'),
            tag = cms.string('GEMAlignmentErrorExtendedRcd')
        ),
        cms.PSet(
            record = cms.string('GeometryFileExtended2024'),
            tag = cms.string('GeometryFileExtended2024')
        ),
        cms.PSet(
            record = cms.string('GEMRecoGeometryRcd'),
            tag = cms.string('GEMRecoGeometryRcd')
	),
        cms.PSet(
            record = cms.string('GlobalPositionRcd'),
            tag = cms.string('GlobalPositionRcd')
        )
    ),

    connect = cms.string('sqlite_file:2025_GEM_data_v1.db')
)
process.p1 = cms.Path(process.GEMAlDBWriter)
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string('INFO'),
    default = cms.untracked.PSet(
        limit = cms.untracked.int32(10000000)
    )
)

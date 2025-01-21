#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

#section General
config.General.requestName = 'Run2023C_v0' 
config.General.workArea = 'Run3_crabLogs'#working directory
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_GE11ana.py'
config.JobType.numCores = 1

misalign = False  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  config.JobType.inputFiles =  ['./db_file.db']

#section Data
#config.Data.runRange = '362695'
config.Data.inputDataset = '/Muon0/Run2023C-ZMu-PromptReco-v4/RAW-RECO'
#config.Data.userInputFiles = open('singleMuonGun_11_3_4_2021_design.list').readlines()
config.Data.lumiMask = '/eos/user/c/cmsdqm/www/CAF/certification/Collisions23/Cert_Collisions2023_366442_370790_Muon.json'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/toakhter/tamu_mual/2023/2023C'
config.Data.publication = False

config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_CH_CERNBOX' #'T3_US_FNALLPC'

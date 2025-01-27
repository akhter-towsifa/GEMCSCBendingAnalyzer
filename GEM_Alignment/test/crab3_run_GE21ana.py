#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
#section general
config.General.requestName = 'Run2024H_muon0_v0' #Run2023D_muon0_alignedreco_v1
config.General.workArea = 'crabLogs'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_GE21ana.py'
config.JobType.numCores = 1

misalign = False  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  config.JobType.inputFiles =  ['./filename.db']

#section Data
#config.Data.runRange = '348776,348773,349073'
config.Data.inputDataset = '/Muon0/Run2024H-ZMu-PromptReco-v1/RAW-RECO'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/toakhter/tamu_mual/2024/2024H'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'

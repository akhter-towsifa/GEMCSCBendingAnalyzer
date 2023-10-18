#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
#section general
config.General.requestName = 'Run2023B_ALCARECO_v7'
config.General.workArea = 'crabLogs'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_GE11ana.py'
config.JobType.numCores = 1

misalign = False  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  config.JobType.inputFiles =  ['./2022D_backPropModiefiedRefitTracker_alcareco_v0.db']

#section Data
#config.Data.runRange = '348776,348773,349073'
config.Data.inputDataset = '/Muon/Run2022D-MuAlCalIsolatedMu-PromptReco-v2/ALCARECO'

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'

config.Data.outLFNDirBase = '/store/user/your_username/Run3_ALCARECO'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName

config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.storageSite = 'T3_CH_CERNBOX'
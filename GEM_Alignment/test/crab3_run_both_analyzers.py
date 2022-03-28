#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'RdPhiAna_CRAFT2022_250322_ExpressCosmics'
config.General.workArea = 'RdPhiAna_CRAFT2022_ExpressCosmics'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_both_analyzers.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

misalign = True  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  config.JobType.inputFiles = ['./CRAFT2022_ME11Only_Iter2.db']

#section Data
config.Data.inputDataset = '/ExpressCosmics/Commissioning2022-Express-v1/FEVT' #ExpressCosmics is bigger than Cosmics
config.Data.runRange = '348776,348773,349073'


#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/daebi/CRAFT2022/'
config.Data.publication = False
config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
config.Site.ignoreGlobalBlacklist = True

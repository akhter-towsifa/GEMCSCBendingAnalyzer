#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'Run2023B_ALCARECO_v7'
config.General.workArea = 'crabLogs'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_GE11ana.py'
#config.JobType.maxMemoryMB = 3500
#config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
#config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

misalign = False  #Make sure to change the run_GE11ana.py too!!!
if misalign:
  config.JobType.inputFiles =  ['./2022D_backPropModiefiedRefitTracker_alcareco_v0.db']

#section Data
#config.Data.runRange = '348776,348773,349073'
#config.Data.inputDataset = '/Muon/Run2022D-MuAlCalIsolatedMu-PromptReco-v2/ALCARECO'
#config.Data.inputDataset = '/Muon/Run2022D-ZMu-PromptReco-v2/RAW-RECO'
config.Data.inputDataset = '/Muon0/Run2023B-MuAlCalIsolatedMu-PromptReco-v1/ALCARECO'

#config.Data.userInputFiles = open('Run2023B-MuAlCalIsolatedMu-PromptReco-v1-exp-2023-08-05.list').readlines()
#config.Data.runRange = '342810,342966,343034,343082,343171,343266,343387,344134,344186,344266,344366'

#config.Data.inputDBS = 'phys03'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/toakhter/tamu_mual/2023/2023B'
#config.Data.outLFNDirBase = '/store/group/lpcgem/'
config.Data.publication = False
#import FWCore.PythonUtilities.LumiList as LumiList
##lumiList = LumiList(filename='my_original_lumi_mask.json')
#lumiList = LumiList(filename='320887_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt')
#lumiList.selectRuns(runs = [321475, 321461,  321457,  321434,  321433,  321432,  321431,  321415,  321414,  321396,  321393,  321313,  321312,  321311, 321310,  321305,  321218,  321178,  321177,  321167,  321166,  321165,  321164,  321162,  321149,  321140,  321138,  321134,  321126, 321123,  321122,  321121,  321119,  321069,  321068,  321067,  321055,  321051,  320996,  320995])
#lumiList.writeJSON('my_lumi_mask.json')
#config.Data.lumiMask = 'my_lumi_mask.json'
#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()
#config.Data.runRange = '%d-%d'%(runstart, runend)#'315257-315270'#'278820-278820' # '193093-194075'
config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.ignoreGlobalBlacklist = True
#config.Site.whitelist = ["T2_KR_KISTI"]
#config.Site.whitelist = ["T0_CH_CERN_MSS"]

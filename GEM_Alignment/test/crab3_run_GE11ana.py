#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()

#section General
config.General.requestName = 'Run362695_ZMu_PromptReco_RAWRECO_v0' 
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
config.Data.runRange = '362695' #'348776,348773,349073'
config.Data.inputDataset = '/Muon/Run2022G-ZMu-PromptReco-v1/RAW-RECO'
#config.Data.userInputFiles = open('singleMuonGun_11_3_4_2021_design.list').readlines()

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/toakhter/Run3'
config.Data.publication = False
#import FWCore.PythonUtilities.LumiList as LumiList
##lumiList = LumiList(filename='my_original_lumi_mask.json')
#lumiList = LumiList(filename='320887_13TeV_PromptReco_Collisions18_JSON_MuonPhys.txt')
#lumiList.selectRuns(runs = [321475, 321461,  321457,  321434,  321433,  321432,  321431,  321415,  321414,  321396,  321393,  321313,  321312,  321311, 321310,  321305,  321218,  321178,  321177,  321167,  321166,  321165,  321164,  321162,  321149,  321140,  321138,  321134,  321126, 321123,  321122,  321121,  321119,  321069,  321068,  321067,  321055,  321051,  320996,  320995])
#lumiList.writeJSON('my_lumi_mask.json')
#config.Data.lumiMask = 'my_lumi_mask.json'
#process.source.lumisToProcess = LumiList.LumiList(filename = 'goodList.json').getVLuminosityBlockRange()

config.Data.outputDatasetTag = config.General.requestName
config.Site.storageSite = 'T3_US_FNALLPC'

#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
import os
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'analyser_misalign'
config.General.workArea = '3dof_1134GT2021design_bigger_initial/3dof_first_10keach_bigger_1134GT2021design_1203_feb21'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyser_misalign.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

pwd = os.getcwd()
config.JobType.inputFiles = ['{pwd}/GEMposition_from_ME11_3times_firstfit_10k_each.db'.format(pwd = pwd)]

#section Data
#config.Data.inputDataset = '/Cosmics/Commissioning2021-CosmicTP-PromptReco-v1/RAW-RECO'
config.Data.inputDataset = '/singleMuonGun_pT_20_200_CMSSW_11_3_4_GT_2021_design/hyunyong-crab_singleMuonGun_11_3_4_2021_design_RAW2DIGI_RECO_v3-ce33467258ba6d6e1b97c4b94c6a3a02/USER'
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/daebi/3dof_1134GT2021design_bigger_initial_feb17'
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
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.ignoreGlobalBlacklist = True
#config.Site.whitelist = ["T2_KR_KISTI"]
#config.Site.whitelist = ["T0_CH_CERN_MSS"]

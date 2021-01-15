#from CRABClient.UserUtilities import config, getUsernameFromSiteDB
from CRABClient.UserUtilities import config
config = config()
###2018runA  314472-318876
#section general
config.General.requestName = 'analyser'
config.General.workArea = 'cosmics_MC_0115'#working dir 
config.General.transferOutputs = True
config.General.transferLogs = True

#section JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'analyser.py'
config.JobType.maxMemoryMB = 2000
config.JobType.maxJobRuntimeMin = 1440 # 1440min = 24hours
config.JobType.numCores = 1
config.JobType.allowUndistributedCMSSW = True
#config.JobType.generator
#config.JobType.pyCfgParams
#config.JobType.inputFiles

#config.JobType.inputFiles = ['/uscms/home/daebi/nobackup/analyser/CMSSW_11_0_0/src/GEMCSCBendingAnalyzer/MuonAnalyser/test/test.db']

#section Data
#config.Data.inputDataset = '/SLHC23_patch1_2023Muon_gen_sim_Pt2_50_1M/tahuang-SLHC25_patch1_2023Muon_1M_L1_PU0_Pt2_50_updategemeta-1bf93df4dfbb43dc918bd6e47dedbf79/USER'
#config.Data.inputDataset = '/SingleMuon/Run2016G-v1/RAW'
#config.Data.inputDataset = '/SingleMuon/Run2016H-v1/RAW'
#config.Data.inputDataset = '/SingleMuon/tahuang-RERECO_Run2018D_singlemuon_GEMon_320995-321475_20180917-8a1254d3d0a422ad143b7aead0544ce7/USER'
#config.Data.inputDataset = '/SingleMuon/tahuang-RERECO_Run2018D_singlemuon_GEMon_323470-324200_20181005-4201715c0f9d22f1f7baffdbca473c7b/USER'
#config.Data.inputDataset = '/RelValSingleMuPt15Eta1p7_2p7/CMSSW_10_3_0_pre4-103X_upgrade2023_realistic_v2_2023D17noPUEA1000-v1/GEN-SIM-RECO'
#config.Data.inputDataset = '/SingleMuon/Run2018D-ZMu-PromptReco-v2/RAW-RECO'
#config.Data.inputDataset = '/singleMuonGun_pT-30to200_1102_phase1_2021_realistic/hyunyong-crab_singleMuonGun_pT-30to200_1102_phase1_2021_realistic_step2-1b4eba2dcd577d6bb642bb3e45609e5f/USER'
#config.Data.inputDataset = '/Cosmics/Commissioning2020-MuAlGlobalCosmics-PromptReco-v1/ALCARECO'
config.Data.inputDataset = '/UndergroundCosmicMu_cfi_1102_phase1_2021/hyunyong-crab_UndergroundCosmicMu_1102_phase1_2021_step3-c9a68e0f82896adf960bb4af1c4e9584/USER'
#config.Data.inputDataset = '/singleMuonGun_30pT_1102_phase1_2021_realistic/hyunyong-crab_singleMuonGun_30pT_1102_phase1_2021_realistic_step2-1b4eba2dcd577d6bb642bb3e45609e5f/USER'
#config.Data.inputDataset = '/SingleMu_Pt50_500K_01_22_20/daebi-SingleMu_Pt30_to_200_500K_RAW2DIGI_RECO_phase2-26430e2fb821d1f1477c362253785707/USER'
config.Data.inputDBS = 'phys03'
#config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/daebi/'
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

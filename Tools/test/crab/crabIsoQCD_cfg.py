from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
#config.General.requestName = 'QCD_Pt-15to20_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-20to30_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-30to50_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-50to80_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-80to120_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-120to170_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-170to300_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
#config.General.requestName = 'QCD_Pt-300to470_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'
config.General.requestName = 'QCD_Pt-470to600_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26_Nmu1'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['globalTag=81X_upgrade2017_realistic_v26',
                              'ntupleName=muonPOGNtuple_IsolationStudies.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=10.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
#config.Data.inputDataset = '/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'
config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'

config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 10
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
# config.Data.ignoreLocality  = True
# config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.blacklist = ['T2_US_*']
#config.Site.whitelist = ['T2_IT_Bari']

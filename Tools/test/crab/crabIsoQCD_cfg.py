from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
#config.General.requestName = 'QCD_Pt-15to20_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-20to30_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-30to50_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-50to80_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-80to120_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-120to170_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-170to300_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-300to470_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-470to600_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-600to800_DRPremix-92X_upgrade2017_realistic_v10'
config.General.requestName = 'QCD_Pt-800to1000_DRPremix-92X_upgrade2017_realistic_v10'
#config.General.requestName = 'QCD_Pt-1000toInf_DRPremix-92X_upgrade2017_realistic_v10'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['globalTag=92X_upgrade2017_realistic_v10',
                              'ntupleName=muonPOGNtuple_IsolationStudies_QCD.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=5.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
#config.Data.inputDataset = '/QCD_Pt-15to20_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'
#config.Data.inputDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'

config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 10
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
# config.Data.ignoreLocality  = True
# config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.blacklist = ['T2_US_*']
#config.Site.whitelist = ['T2_IT_Bari']

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'DY_FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['globalTag=81X_upgrade2017_realistic_v26',
                              'ntupleName=muonPOGNtuple_IsolationStudies.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=10.0',
                              'minNMu=2'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/PhaseIFall16DR-FlatPU28to62HcalNZSRAW_81X_upgrade2017_realistic_v26-v1/AODSIM'


config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 2
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.ignoreLocality  = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'


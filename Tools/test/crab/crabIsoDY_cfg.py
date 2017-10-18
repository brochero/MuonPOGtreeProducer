from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferLogs = True
config.General.requestName = 'DY_RunIISummer17DRPremix-92X_upgrade2017_realistic_v10'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['globalTag=92X_upgrade2017_realistic_v10',
                              'ntupleName=muonPOGNtuple_IsolationStudies_DY.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=5.0',
                              'minNMu=1'
               ]
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17DRPremix-92X_upgrade2017_realistic_v10-v2/AODSIM'


config.Data.splitting    = 'FileBased'
config.Data.unitsPerJob  = 20
config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
config.Data.ignoreLocality  = True
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'


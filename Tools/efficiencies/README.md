# Estimation of ISOLATION efficiencies using the MuonPOG ntuples

# Ntuples production  
Procedure to produce the ntuples is well explained in MuonPOGtreeProducer/README.md.

# Ntuples location
Ntuple produced are storage ```/eos/cms/store/group/phys_muon/Commissioning/Ntuples/Commissioning2017``` in directory.

# Efficiency estimation
1. ISOEfficiency.C contains definitions and definitions of all histograms.
   a. Histograms are produced, so far, for all different ID (trk, glb, loose, medium, tight, soft and high pT). 
2. To compile:
```
./ISOEfficiency
```
3. To run ONE job (it runs locally):
```
./ISOEfficiency.run config/CONFIG.ini OutputDir OutputFileName
```
MuonResults is the default output directory.
4. In order to run the full set of samples (can be DY, QCD, ttbar, WJets, etc) the ```SubmitJobs.py``` macro has been created. To run it, first, you must check the samples in each file config included into ```Configs``` variable (line 9). After it, just
```
python SubmitJobs.py OutputFileName
```
5. To check the jobs
```
bjobs
```
Outputs are saved into ```MuonResults```
6. Macros (readme in preparation ;) )

#!/Usr/bin/env python

import sys
import ROOT 
from array import array
import math
import re
from ROOT import *
import numpy as np
import argparse
from mt2 import mt2

try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Example1.py input_file")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass      

inputFile = sys.argv[1]
outputFile = sys.argv[2]
stopqmass = sys.argv[3]
LSPmass = sys.argv[4]
parser = argparse.ArgumentParser()
#parser.add_argument("-o",
#                    help="out put file name,  example -o output.root ")
#parser.add_argument("-mstop",
#                    help="mass of STOP quark,  example -mstop 1000   ", action='store_true')
#parser.add_argument("-mlsp",
#                   help="mass of LSP,  example -mlsp 1   ", action='store_true')
#args =  parser.parse_args()

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
#branchJet      = treeReader.UseBranch("JetPUPPITight")
branchElectron = treeReader.UseBranch("ElectronMedium")
branchWeight   = treeReader.UseBranch("Weight")
branchEvent    = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch('JetPUPPI')
branchPuppiMissingET  = treeReader.UseBranch('PuppiMissingET')
branchPuppiCandidate  = treeReader.UseBranch('ParticleFlowCandidate')
branchRho             = treeReader.UseBranch('Rho')
branchParticle        = treeReader.UseBranch('Particle') 
branchGenJet          = treeReader.UseBranch('GenJet')
branchPuppiJetLoose   = treeReader.UseBranch('JetPUPPILoose')
#branchPuppiJetTight   = treeReader.UseBranch('JetPUPPITight')
branchFatJet          = treeReader.UseBranch('JetPUPPIAK8')


# Book histograms
outputfile = ROOT.TFile(outputFile, 'RECREATE')
outputfile.cd()
stop_lsp = 'stop_'+stopqmass+'_lsp_'+LSPmass
print("stop_lsp is ",stop_lsp)
outputfile.mkdir(stop_lsp)
outputfile.cd(stop_lsp)

nEvents = ROOT.TH1F("nEvents", "total events", 10, 0.0, 10)
nEvents_1 = ROOT.TH1F("nEvents_1", "total events", 10, 0.0, 10)
nEvents_2 = ROOT.TH1F("nEvents_2", "total events", 10, 0.0, 10)
nEvents_3 = ROOT.TH1F("nEvents_3", "total events", 10, 0.0, 10)
nEvents_4 = ROOT.TH1F("nEvents_4", "total events", 10, 0.0, 10)
nEvents_13 = ROOT.TH1F("nEvents_13", "total events", 10, 0.0, 10)
nEntries_1 = ROOT.TH1F("nEntries_1", "total events", 10, 0.0, 10)
nEntries_2 = ROOT.TH1F("nEntries_2", "total events", 10, 0.0, 10)
nEntries_3 = ROOT.TH1F("nEntries_3", "total events", 10, 0.0, 10)
nEntries_4 = ROOT.TH1F("nEntries_4", "total events", 10, 0.0, 10)
nEntries_5 = ROOT.TH1F("nEntries_5", "total events", 10, 0.0, 10)
nEntries_6 = ROOT.TH1F("nEntries_6", "total events", 10, 0.0, 10)
nEntries_7 = ROOT.TH1F("nEntries_7", "total events", 10, 0.0, 10)
nEntries_8 = ROOT.TH1F("nEntries_8", "total events", 10, 0.0, 10)
nEntries_9 = ROOT.TH1F("nEntries_9", "total events", 10, 0.0, 10)
nEntries_10 = ROOT.TH1F("nEntries_10", "total events", 10, 0.0, 10)
nEntries_11 = ROOT.TH1F("nEntries_11", "total events", 10, 0.0, 10)


nEvts_bin1 = ROOT.TH1F("nEvts_bin1", "total events", 10, 0.0, 10)
nEvts_bin2 = ROOT.TH1F("nEvts_bin2", "total events", 10, 0.0, 10)
nEvts_bin3 = ROOT.TH1F("nEvts_bin3", "total events", 10, 0.0, 10)
nEvts_bin4 = ROOT.TH1F("nEvts_bin4", "total events", 10, 0.0, 10)
nEvts_bin5 = ROOT.TH1F("nEvts_bin5", "total events", 10, 0.0, 10)
nEvts_bin6 = ROOT.TH1F("nEvts_bin6", "total events", 10, 0.0, 10)
nEvts_bin7 = ROOT.TH1F("nEvts_bin7", "total events", 10, 0.0, 10)
nEvts_bin8 = ROOT.TH1F("nEvts_bin8", "total events", 10, 0.0, 10)
nEvts_bin9 = ROOT.TH1F("nEvts_bin9", "total events", 10, 0.0, 10)
nEvts_bin10 = ROOT.TH1F("nEvts_bin10", "total events", 10, 0.0, 10)

nEvts_bin11 = ROOT.TH1F("nEvts_bin11", "total events", 10, 0.0, 10)
nEvts_bin12 = ROOT.TH1F("nEvts_bin12", "total events", 10, 0.0, 10)
nEvts_bin13 = ROOT.TH1F("nEvts_bin13", "total events", 10, 0.0, 10)
nEvts_bin14 = ROOT.TH1F("nEvts_bin14", "total events", 10, 0.0, 10)
nEvts_bin15 = ROOT.TH1F("nEvts_bin15", "total events", 10, 0.0, 10)
nEvts_bin16 = ROOT.TH1F("nEvts_bin16", "total events", 10, 0.0, 10)
nEvts_bin17 = ROOT.TH1F("nEvts_bin17", "total events", 10, 0.0, 10)
nEvts_bin18 = ROOT.TH1F("nEvts_bin18", "total events", 10, 0.0, 10)
nEvts_bin19 = ROOT.TH1F("nEvts_bin19", "total events", 10, 0.0, 10)
nEvts_bin20 = ROOT.TH1F("nEvts_bin20", "total events", 10, 0.0, 10)

nEvts_bin21 = ROOT.TH1F("nEvts_bin21", "total events", 10, 0.0, 10)
nEvts_bin22 = ROOT.TH1F("nEvts_bin22", "total events", 10, 0.0, 10)
nEvts_bin23 = ROOT.TH1F("nEvts_bin23", "total events", 10, 0.0, 10)
nEvts_bin24 = ROOT.TH1F("nEvts_bin24", "total events", 10, 0.0, 10)
nEvts_bin25 = ROOT.TH1F("nEvts_bin25", "total events", 10, 0.0, 10)
nEvts_bin26 = ROOT.TH1F("nEvts_bin26", "total events", 10, 0.0, 10)
nEvts_bin27 = ROOT.TH1F("nEvts_bin27", "total events", 10, 0.0, 10)
nEvts_bin28 = ROOT.TH1F("nEvts_bin28", "total events", 10, 0.0, 10)
nEvts_bin29 = ROOT.TH1F("nEvts_bin29", "total events", 10, 0.0, 10)
nEvts_bin30 = ROOT.TH1F("nEvts_bin30", "total events", 10, 0.0, 10)

tauPT_1 = ROOT.TH1F("tau1_pt", "tau P_{T}", 100, 30.0, 5000.0)
tauPT_2 = ROOT.TH1F("tau2_pt", "tau P_{T}", 100, 30.0, 5000.0)
metPT = ROOT.TH1F("MET", "MET", 50, 0.0, 5000.0)
HT_Tot = ROOT.TH1F("HT", "Sum P_{T}", 100, 0.0, 5000.0)
HT_Tot_1 = ROOT.TH1F("HT_1", "Sum P_{T}", 100, 0.0, 5000.0)
ptratio_tau1 = ROOT.TH1F("ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
ptratio_tau2 = ROOT.TH1F("ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)

genmatch_ptratio_tau1 = ROOT.TH1F("genmatch_ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_1030 = ROOT.TH1F("genmatchedR_PT1_tau0", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_3050 = ROOT.TH1F("genmatchedR_PT1_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_5080 = ROOT.TH1F("genmatchedR_PT2_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_80130 = ROOT.TH1F("genmatchedR_PT3_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_130180 = ROOT.TH1F("genmatchedR_PT4_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_180230 = ROOT.TH1F("genmatchedR_PT5_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_230280 = ROOT.TH1F("genmatchedR_PT6_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_280330 = ROOT.TH1F("genmatchedR_PT7_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_330380 = ROOT.TH1F("genmatchedR_PT8_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_R_380500 = ROOT.TH1F("genmatchedR_PT9_tau1", "P_{T} Ratio", 50, 0.0, 2.0)

genmatch_ptratio_tau2 = ROOT.TH1F("genmatch_ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_1030 = ROOT.TH1F("genmatchedR_PT0_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_3050 = ROOT.TH1F("genmatchedR_PT1_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_5080 = ROOT.TH1F("genmatchedR_PT2_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_80130 = ROOT.TH1F("genmatchedR_PT3_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_130180 = ROOT.TH1F("genmatchedR_PT4_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_180230 = ROOT.TH1F("genmatchedR_PT5_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_230280 = ROOT.TH1F("genmatchedR_PT6_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_280330 = ROOT.TH1F("genmatchedR_PT7_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_330380 = ROOT.TH1F("genmatchedR_PT8_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_R_380500 = ROOT.TH1F("genmatchedR_PT9_tau2", "P_{T} Ratio", 50, 0.0, 2.0)

notgenmatch_ptratio_tau1 = ROOT.TH1F("notgenmatch_ptratio_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_1030 = ROOT.TH1F("notgenmatched_R_PT0_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_3050 = ROOT.TH1F("notgenmatched_R_PT1_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_5080 = ROOT.TH1F("notgenmatched_R_PT2_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_80130 = ROOT.TH1F("notgenmatched_R_PT3_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_130180 = ROOT.TH1F("notgenmatched_R_PT4_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_180230 = ROOT.TH1F("notgenmatched_R_PT5_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_230280 = ROOT.TH1F("notgenmatched_R_PT6_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_280330 = ROOT.TH1F("notgenmatched_R_PT7_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_330380 = ROOT.TH1F("notgenmatched_R_PT8_tau1", "P_{T} Ratio", 50, 0.0, 2.0)
t1_notgR_380500 = ROOT.TH1F("notgenmatched_R_PT9_tau1", "P_{T} Ratio", 50, 0.0, 2.0)

notgenmatch_ptratio_tau2 = ROOT.TH1F("notgenmatch_ptratio_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_1030 = ROOT.TH1F("notgenmatched_R_PT0_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_3050 = ROOT.TH1F("notgenmatched_R_PT1_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_5080 = ROOT.TH1F("notgenmatched_R_PT2_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_80130 = ROOT.TH1F("notgenmatched_R_PT3_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_130180 = ROOT.TH1F("notgenmatched_R_PT4_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_180230 = ROOT.TH1F("notgenmatched_R_PT5_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_230280 = ROOT.TH1F("notgenmatched_R_PT6_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_280330 = ROOT.TH1F("notgenmatched_R_PT7_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_330380 = ROOT.TH1F("notgenmatched_R_PT8_tau2", "P_{T} Ratio", 50, 0.0, 2.0)
t2_notgR_380500 = ROOT.TH1F("notgenmatched_R_PT9_tau2", "P_{T} Ratio", 50, 0.0, 2.0)

MT = ROOT.TH1F("MT", "MT", 50, 0.0, 5000.0)
DR_daughter = ROOT.TH1F("DeltaR","Delta R ",2000,0.0,3)
DR_nr_genreco1 = ROOT.TH1F("DeltaR_1","Delta R ",2000,0.0,3)
DR_nr_genreco2 = ROOT.TH1F("DeltaR_2","Delta R ",2000,0.0,3)

#histElectronPT = ROOT.TH1F("Electron_pt", "electron P_{T}", 100, 0.0, 1000.0)


#------------ functions for gen index ----------
def get_gentau(igen):
  tauindex = igen
  return tauindex



# Loop over all events
count_1 = 0
count_2 = 0
count_3 = 0
count_4 = 0
count_5 = 0
count_6 = 0
count_7 = 0
count_8 = 0
count_9 = 0
count_10 = 0
count_11 = 0
count_13 = 0
Bin_1 = 0
Bin_2 = 0
Bin_3 = 0
Bin_4 = 0
Bin_5 = 0
Bin_6 = 0
Bin_7 = 0
Bin_8 = 0
Bin_9 = 0
Bin_10 = 0
Bin_11 = 0
Bin_12 = 0
Bin_13 = 0
Bin_14 = 0
Bin_15 = 0
Bin_16 = 0
Bin_17 = 0
Bin_18 = 0
Bin_19 = 0
Bin_20 = 0
Bin_21 = 0
Bin_22 = 0
Bin_23 = 0
Bin_24 = 0
Bin_25 = 0
Bin_26 = 0
Bin_27 = 0
Bin_28 = 0
Bin_29 = 0
Bin_30 = 0

for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)
  nEvents.Fill(1)
  nEvents.Fill(2,numberOfEntries)
  # Choose STOP and LSP mass
  STOP_idx = -1
  LSP_idx = -1

  for igen,gen in enumerate(branchParticle):
    if(abs(gen.PID) == 1000006):
      STOP_idx = igen
    if(abs(gen.PID) == 1000022 and igen != STOP_idx):
      LSP_idx = igen                                                                                             
  if(not(STOP_idx >= 0 and LSP_idx>= 0)): continue
  STOP = branchParticle.At(STOP_idx)
  LSP = branchParticle.At(LSP_idx)
  if(not(STOP.Mass == stopqmass and LSP.Mass == LSPmass)): continue #[STOP,LSP] = [1000,1], [300,1], [500,350]
  count_1 += 1
  nEvents_1.Fill(1)
  nEvents_1.Fill(2,count_1)
 
  # Take first jet and get tau1 and tau1 
  tau1_idx = -1
  tau2_idx = -1
  tau1_tau2_HT = -1
  Tltau1_p4 = TLorentzVector()
  Tltau2_p4 = TLorentzVector()
  i = 0
  #print("coming in entries,", entry)
  for iTau1, tau1 in enumerate(branchJet) :
    i +=1
    #print("coming inside tau1")
    tautagOk1 = ( tau1.TauTag & (1 << 2) )
    if (not (tau1.PT >30 and abs(tau1.Eta) < 3. and tautagOk1 and abs(tau1.Charge) == 1)): continue
    for iTau2, tau2 in enumerate(branchJet):
      #print("coming inside tau2")
      if (iTau1 == iTau2) : continue
      tautagOk2 = ( tau2.TauTag & (1 << 2) )
      if (not (tau2.PT >30 and abs(tau2.Eta) < 3. and tautagOk2 and abs(tau2.Charge) == 1)): continue
      if (tau1.Charge*tau2.Charge < 0):
        HT = tau1.PT + tau2.PT
        if (HT > tau1_tau2_HT) :
          tau1_idx = iTau1
          tau2_idx = iTau2
          tau1_tau2_HT = HT
  
  if (not(tau1_idx >= 0 and tau2_idx >= 0)): continue
  tau1 = branchJet.At(tau1_idx)
  tau2 = branchJet.At(tau2_idx)
  if(tau1.PT < tau2.PT):
    tau = tau1
    tau1 = tau2
    tau2 = tau
  Tltau1_p4.SetPtEtaPhiM(tau1.PT, tau1.Eta, tau1.Phi, tau1.Mass)
  Tltau2_p4.SetPtEtaPhiM(tau2.PT, tau2.Eta, tau2.Phi, tau2.Mass)
  tau1tau2_m = (Tltau1_p4 + Tltau2_p4).M()

  btag_idx = -1    
  for ibjet, bjet in enumerate(branchJet) :
    if (ibjet == tau1_idx or ibjet == tau2_idx): continue
    
    btagok = (bjet.BTag & (1 << 1) )
    if (bjet.PT > 30 and abs(bjet.Eta) < 5. and btagok):
      btag_idx = ibjet
  #print( "btag index" , btag_idx)

  HT_Total = -1
  for ijet, jet in enumerate(branchJet) :
    HT_Total += jet.PT
  #print("Total HT", HT_Total)
    
  Met_PT = -1
  imet_idx = -1
  Met_Phi = 0
  for imet, met in enumerate(branchPuppiMissingET):
    Met_PT = met.MET
    Met_Phi = met.Phi
    imet_idx = imet
  #print(" imet_idx", imet_idx)
 
  if (not (tau1_idx >= 0 and tau2_idx >= 0 and btag_idx >= 0 and HT_Total > 100 and Met_PT > 50 and tau1tau2_m > 100)): continue

  #tauPT_1.Fill(tau1pt)
  #tauPT_2.Fill(tau2pt)
  #  metP.Fill(metpt)
  #HT_Tot_1.Fill(HT_Total)
  count_2 += 1
  nEvents_2.Fill(1)
  nEvents_2.Fill(2,count_2)

  tau_1 = branchJet.At(tau1_idx)
  tau_2 = branchJet.At(tau2_idx)
  met_pt = branchPuppiMissingET.At(imet_idx)
  tau1pt = tau_1.PT
  tau2pt = tau_2.PT
  metpt = met_pt.MET
  
  leadchtau1 = -1
  leadchtau2 = -1

  ######----------leading chtau2 --------------#########
  isLeptonic = False  
  tau1_leadCH = None
  all_consti1_p4 = TLorentzVector()
  for consti in tau1.Constituents:
    const_p4 =  TLorentzVector()
    const_p4.SetPtEtaPhiM(consti.PT, consti.Eta, consti.Phi, consti.Mass)
    all_consti1_p4 += const_p4
    ids = consti.PID
    if (abs(ids) in [11, 12, 13, 14]):
      isLeptonic = True
    if (isLeptonic == True or consti.Charge == 0) :
      continue
    const_p4 = TLorentzVector()
    const_p4.SetPtEtaPhiM(consti.PT, consti.Eta, consti.Phi, consti.Mass)
    dR = const_p4.DeltaR(Tltau1_p4)
    if (dR < 0.1):
      #print(dR)
      chpt = consti.PT
      if (tau1_leadCH is None or consti.PT > tau1_leadCH.PT):
        tau1_leadCH = consti
  if (tau1_leadCH is not None):
    leadchtau1 =  tau1_leadCH.PT/all_consti1_p4.Pt()
    print("leadchtau1 in function before gen match: ", leadchtau1)
    ptratio_tau1.Fill(leadchtau1)


  ######----------leading chtau2 --------------#########
  isLeptonic2 = False  
  tau2_leadCH = None
  all_consti2_p4= TLorentzVector()
  for consti2 in tau2.Constituents:
    const2_p4 = TLorentzVector()
    const2_p4.SetPtEtaPhiM(consti2.PT, consti2.Eta, consti2.Phi, consti2.Mass)
    all_consti2_p4 += const2_p4
    ids2 = consti2.PID
    if (abs(ids2) in [11, 12, 13, 14]):
      isLeptonic2 = True
    if (isLeptonic2 == True or consti2.Charge == 0) :
      continue
    const2_p4 = TLorentzVector()
    const2_p4.SetPtEtaPhiM(consti2.PT, consti2.Eta, consti2.Phi, consti2.Mass)
    dR2 = const2_p4.DeltaR(Tltau2_p4)
    
    if (dR2 < 0.1):
      #print("coming inside 3 dr2 :",dR2)  
      chpt = consti.PT
      if (tau2_leadCH is None or consti2.PT > tau2_leadCH.PT):
        tau2_leadCH = consti2
  if (tau2_leadCH is not None):
    leadchtau2 =  tau2_leadCH.PT/all_consti2_p4.Pt()
    #print("leadchtau2 is ",leadchtau2)
    ptratio_tau2.Fill(leadchtau2)

    
  ##########-----MT-----------#######
  tau1_px = Tltau1_p4.Px()
  tau1_py = Tltau1_p4.Py()
  tau1_m = tau_1.Mass
  tau2_px = Tltau2_p4.Px()
  tau2_py = Tltau2_p4.Py()
  tau2_m = tau_2.Mass
  MET_px = Met_PT*(math.cos(Met_Phi))
  MET_py = Met_PT*(math.sin(Met_Phi))
  MT_ = mt2(
    tau1_m,tau1_px,tau1_py,
    tau2_m,tau2_px,tau2_py,
    MET_px,MET_py,
    0,0)
  MT.Fill(MT_)
  

  if(leadchtau1 <= 0.5 and leadchtau2 <= 0.5):
    count_3 += 1
    nEvents_3.Fill(1)
    nEvents_3.Fill(2,count_3)
    if(metpt >= 50 and  metpt < 200):
      if(MT < 40):
        if(HT >= 100 and HT < 300):
          Bin_1 += 1
          nEvts_bin1.Fill(1)
          nEvts_bin1.Fill(2,Bin_1)
        elif(HT >= 300 and HT < 700): 
          Bin_2 += 1
          nEvts_bin2.Fill(1)
          nEvts_bin2.Fill(2,Bin_2)
        elif(HT >= 700):
          Bin_3 += 1
          nEvts_bin3.Fill(1)
          nEvts_bin3.Fill(2,Bin_3)
      
      elif(MT >= 40 and MT < 80 ):
        if(HT >= 100 and HT < 300):
          Bin_6 += 1
          nEvts_bin6.Fill(1)
          nEvts_bin6.Fill(2,Bin_6)
        elif(HT >= 300 and HT < 700): 
          Bin_7 += 1
          nEvts_bin7.Fill(1)
          nEvts_bin7.Fill(2,Bin_7)
        elif(HT >= 700):
          Bin_8 += 1
          nEvts_bin8.Fill(1)
          nEvts_bin8.Fill(2,Bin_8)
        
      elif(MT >= 80):
        if(HT >= 100 and HT < 300):
          Bin_11 += 1
          nEvts_bin11.Fill(1)
          nEvts_bin11.Fill(2,Bin_11)
        elif(HT >= 300 and HT < 700): 
          Bin_12 += 1
          nEvts_bin12.Fill(1)
          nEvts_bin12.Fill(2,Bin_12)
        elif(HT >= 700):
          Bin_13 += 1
          nEvts_bin13.Fill(1)
          nEvts_bin13.Fill(2,Bin_13)
     
    elif(metpt >= 200):   ############         if(metpt >= 50 and metpt < 200):
      if(MT < 40):
        if(HT >= 100 and HT < 700):
          Bin_4 += 1
          nEvts_bin4.Fill(1)
          nEvts_bin4.Fill(2,Bin_4)
        elif(HT >= 700):
          Bin_5 += 1
          nEvts_bin5.Fill(1)
          nEvts_bin5.Fill(2,Bin_5)
      
      elif(MT >= 40 and MT < 80 ):
        if(HT >= 100 and HT < 700):
          Bin_9 += 1
          nEvts_bin9.Fill(1)
          nEvts_bin9.Fill(2,Bin_9)
        elif(HT >= 700):
          Bin_10 += 1
          nEvts_bin10.Fill(1)
          nEvts_bin10.Fill(2,Bin_10)
          
      elif(MT >= 80):
        if(HT >= 100 and HT < 700):
          Bin_14 += 1
          nEvts_bin14.Fill(1)
          nEvts_bin14.Fill(2,Bin_14)
        elif(HT >= 700):
          Bin_15 += 1
          nEvts_bin15.Fill(1)
          nEvts_bin15.Fill(2,Bin_15)
          
          ######################       if(leadchtau1 <= 0.5 and leadchtau2 <= 0.5)):
  else:
    count_4 += 1
    nEvents_4.Fill(1)
    nEvents_4.Fill(2,count_4)
    if(metpt >= 50 and  metpt < 200):
      if(MT < 40):
        if(HT >= 100 and HT < 300):
          Bin_16 += 1
          nEvts_bin16.Fill(1)
          nEvts_bin16.Fill(2,Bin_16)
        elif(HT >= 300 and HT < 700): 
          Bin_17 += 1
          nEvts_bin17.Fill(1)
          nEvts_bin17.Fill(2,Bin_17)
        elif(HT >= 700):
          Bin_18 += 1
          nEvts_bin18.Fill(1)
          nEvts_bin18.Fill(2,Bin_18)
      
      elif(MT >= 40 and MT < 80 ):
        if(HT >= 100 and HT < 300):
          Bin_21 += 1
          nEvts_bin21.Fill(1)
          nEvts_bin21.Fill(2,Bin_21)
        elif(HT >= 300 and HT < 700): 
          Bin_22 += 1
          nEvts_bin22.Fill(1)
          nEvts_bin22.Fill(2,Bin_22)
        elif(HT >= 700):
          Bin_23 += 1
          nEvts_bin23.Fill(1)
          nEvts_bin23.Fill(2,Bin_23)
        
      elif(MT >= 80):
        if(HT >= 100 and HT < 300):
          Bin_26 += 1
          nEvts_bin26.Fill(1)
          nEvts_bin26.Fill(2,Bin_26)
        elif(HT >= 300 and HT < 700): 
          Bin_27 += 1
          nEvts_bin27.Fill(1)
          nEvts_bin27.Fill(2,Bin_27)
        elif(HT >= 700):
          Bin_28 += 1
          nEvts_bin28.Fill(1)
          nEvts_bin28.Fill(2,Bin_28)

    elif(metpt >= 200):   ############         if(metpt >= 50 and metpt < 200):
      if(MT < 40):
        if(HT >= 100 and HT < 700):
          Bin_19 += 1
          nEvts_bin19.Fill(1)
          nEvts_bin19.Fill(2,Bin_19)
        elif(HT >= 700):
          Bin_20 += 1
          nEvts_bin20.Fill(1)
          nEvts_bin20.Fill(2,Bin_20)
      
      elif(MT >= 40 and MT < 80 ):
        if(HT >= 100 and HT < 700):
          Bin_24 += 1
          nEvts_bin24.Fill(1)
          nEvts_bin24.Fill(2,Bin_24)
        elif(HT >= 700):
          Bin_25 += 1
          nEvts_bin25.Fill(1)
          nEvts_bin25.Fill(2,Bin_25)
          
      elif(MT >= 80):
        if(HT >= 100 and HT < 700):
          Bin_29 += 1
          nEvts_bin29.Fill(1)
          nEvts_bin29.Fill(2,Bin_29)
        elif(HT >= 700):
          Bin_30 += 1
          nEvts_bin30.Fill(1)
          nEvts_bin30.Fill(2,Bin_30)



  tauPT_1.Fill(tau1pt)
  tauPT_2.Fill(tau2pt)
  metPT.Fill(metpt)
  HT_Tot.Fill(HT_Total)

  gen1_idx = gen2_idx = -1           
  taucand1 = taucand2 = -1    
  gen_1 = gen_2 = None 
  gen_1PT = -1
  gen_2PT = -1
  nele = 0
  min_dr_1 = 999.9
  min_dr_2 = 999.9
  for igen,gen in enumerate(branchParticle):

    #print("gen index at starting of gen loop",igen)
    gen_tau = None        
    #
    #dr_dau = -1
    gen1_p4 = gen2_p4 = TLorentzVector()
    gen_p4 = TLorentzVector()
    gen_tau_p4 = TLorentzVector()
    if(not(abs(gen.PID) == 15)): continue
    count_7 += 1
    gen_tau = igen
    gen_tau_p4.SetPtEtaPhiM(gen.PT, gen.Eta, gen.Phi, gen.Mass)
      
    #here we are checking for the hadronically decay of the gentau. 
    #As daughter info was not stored in delphys, so checking by dR < 0.1  b/w gen taus and the gen particles
    #if tau decayed leptonically, then there would be a lepton(e/mu) in dR < 0.1

    for jgen,genlep in enumerate(branchParticle):
      gen_p4.SetPtEtaPhiM(genlep.PT, genlep.Eta, genlep.Phi, genlep.Mass)
      dr_gentau = gen_p4.DeltaR(gen_tau_p4)
      #if(jgen == igen or dr_gentau > 0.1   ): continue
      if(dr_gentau > 0.1   ): continue
      if((abs(genlep.PID) ==  11) or (abs(genlep.PID) ==  13) ): continue
      dr_1 = gen_tau_p4.DeltaR(Tltau1_p4)
      dr_2 = gen_tau_p4.DeltaR(Tltau2_p4)
      if(dr_1 < 0.3 and dr_1 < min_dr_1 ):
        min_dr_1 = dr_1
        taucand1 = get_gentau(igen)
            
      if(dr_2 < 0.3 and dr_2 < min_dr_2):
        min_dr_2 = dr_2  
        taucand2 = get_gentau(igen)
      
    print("taucand1 in loop ", taucand1)
    print("dr_1 in loop ", dr_1)
    print("min_dr_1 in loop ", min_dr_1)

    print("taucand2 in loop ", taucand2)
    print("dr_2 in loop ", dr_2)
    print("min_dr_2 in loop ", min_dr_2)
        
          
  DR_nr_genreco1.Fill(min_dr_1)
  DR_nr_genreco2.Fill(min_dr_2)
  print("gen1_index outside gen loop ", taucand1)
  print("gen2_index outside gen loop ", taucand2)
  
  if (taucand1 >= 0):
    gen_1 = branchParticle.At(taucand1)
    gen1_p4.SetPtEtaPhiM(gen_1.PT, gen_1.Eta, gen_1.Phi, gen_1.Mass)     
  
  if (taucand2 >= 0):
    gen_2 = branchParticle.At(taucand2)
    gen2_p4.SetPtEtaPhiM(gen_2.PT, gen_2.Eta, gen_2.Phi, gen_2.Mass)     


  if(gen_1 is not None):
    print("gen1_index ", taucand1)
    #gen_1PT = gen_1.PT
    
    #nEntries_3.Fill(1)
    nEntries_3.Fill(2,count_3)
    genmatch_ptratio_tau1.Fill(leadchtau1)
    if (gen_1PT > 10  and gen_1PT < 30):
      t1_R_1030.Fill(leadchtau1)

    if (gen_1PT > 30  and gen_1PT < 50):
      t1_R_3050.Fill(leadchtau1)

    if (gen_1PT > 50 and gen_1PT < 80):
      t1_R_5080.Fill(leadchtau1)

    if (gen_1PT > 80 and gen_1PT < 130):
      t1_R_80130.Fill(leadchtau1)

    if (gen_1PT > 130 and gen_1PT < 180):
      t1_R_130180.Fill(leadchtau1)

    if (gen_1PT > 180 and gen_1PT < 230):
      t1_R_180230.Fill(leadchtau1)

    if (gen_1PT > 230 and gen_1PT < 280):
      t1_R_230280.Fill(leadchtau1)

    if (gen_1PT > 280 and gen_1PT < 330):
      t1_R_280330.Fill(leadchtau1)

    if (gen_1PT > 330 and gen_1PT < 380):
      t1_R_330380.Fill(leadchtau1)

    if (gen_1PT > 380 and gen_1PT < 500):
      t1_R_380500.Fill(leadchtau1)

  else: 
    #print (leadchtau1)
    notgenmatch_ptratio_tau1.Fill(leadchtau1)
    if (gen_1PT > 10  and gen_1PT < 30):
      t1_notgR_1030.Fill(leadchtau1)

    if (gen_1PT > 30  and gen_1PT < 50):
      t1_notgR_3050.Fill(leadchtau1)

    if (gen_1PT > 50 and gen_1PT < 80):
      t1_notgR_5080.Fill(leadchtau1)

    if (gen_1PT > 80 and gen_1PT < 130):
      t1_notgR_80130.Fill(leadchtau1)

    if (gen_1PT > 130 and gen_1PT < 180):
      t1_notgR_130180.Fill(leadchtau1)

    if (gen_1PT > 180 and gen_1PT < 230):
      t1_notgR_180230.Fill(leadchtau1)

    if (gen_1PT > 230 and gen_1PT < 280):
      t1_notgR_230280.Fill(leadchtau1)

    if (gen_1PT > 280 and gen_1PT < 330):
      t1_notgR_280330.Fill(leadchtau1)

    if (gen_1PT > 330 and gen_1PT < 380):
      t1_notgR_330380.Fill(leadchtau1)

    if (gen_1PT > 380 and gen_1PT < 500):
      t1_notgR_380500.Fill(leadchtau1)


  if (gen_2 is not None):
    gen_2PT = gen_2.PT  
    genmatch_ptratio_tau2.Fill(leadchtau2)
    if (gen_2PT > 30  and gen_2PT < 50):
      #print("gen_2 pt3050 is : ",gen_2PT)
      t2_R_3050.Fill(leadchtau2)

    if (gen_2PT > 10  and gen_2PT < 30):
      #print("gen_2 pt3050 is : ",gen_2PT)
      t2_R_1030.Fill(leadchtau2)

    if (gen_2PT > 50 and gen_2PT < 80):
      #print("gen_2 pt5080 is : ",gen_2PT)
      t2_R_5080.Fill(leadchtau2)

    if (gen_2PT > 80 and gen_2PT < 130):
      t2_R_80130.Fill(leadchtau2)

    if (gen_2PT > 130 and gen_2PT < 180):
      t2_R_130180.Fill(leadchtau2)

    if (gen_2PT > 180 and gen_2PT < 230):
      t2_R_180230.Fill(leadchtau2)

    if (gen_2PT > 230 and gen_2PT < 280):
      t2_R_230280.Fill(leadchtau2)

    if (gen_2PT > 280 and gen_2PT < 330):
      t2_R_280330.Fill(leadchtau2)

    if (gen_2PT > 330 and gen_2PT < 380):
      t2_R_330380.Fill(leadchtau2)

    if (gen_2PT > 380 and gen_2PT < 500):
      t2_R_380500.Fill(leadchtau2)

  

  else:
    if (gen_2PT > 10  and gen_2PT < 30):
      t2_notgR_1030.Fill(leadchtau2)

    if (gen_2PT > 30  and gen_2PT < 50):
      t2_notgR_3050.Fill(leadchtau2)

    if (gen_2PT > 50 and gen_2PT < 80):
      t2_notgR_5080.Fill(leadchtau2)

    if (gen_2PT > 80 and gen_2PT < 130):
      t2_notgR_80130.Fill(leadchtau2)

    if (gen_2PT > 130 and gen_2PT < 180):
      t2_notgR_130180.Fill(leadchtau2)

    if (gen_2PT > 180 and gen_2PT < 230):
      t2_notgR_180230.Fill(leadchtau2)

    if (gen_2PT > 230 and gen_2PT < 280):
      t2_notgR_230280.Fill(leadchtau2)

    if (gen_2PT > 280 and gen_2PT < 330):
      t2_notgR_280330.Fill(leadchtau2)

    if (gen_2PT > 330 and gen_2PT < 380):
      t2_notgR_330380.Fill(leadchtau2)

    if (gen_2PT > 380 and gen_2PT < 500):
      t2_notgR_380500.Fill(leadchtau2)



  lep_p4 = TLorentzVector()
  lep1_dr = lep2_dr = -1
  min_dr = 999.9 
  for jgen,gen in enumerate(branchParticle):

    if(abs(gen.PID) in [11,13]):
      lep_p4.SetPtEtaPhiM(gen.PT, gen.Eta, gen.Phi, gen.Mass)
      if (gen_1 is not None): 
        lep1_dr = lep_p4.DeltaR(gen1_p4)       
      elif(gen_2 is not None):   
        lep2_dr = lep_p4.DeltaR(gen2_p4)
      if (lep1_dr < min_dr):
        min_dr = lep1_dr
      elif (lep2_dr < min_dr):
        min_dr = lep2_dr
  DR_daughter.Fill(min_dr)
  #print("min dr is : ",min_dr)



  
  

#cnv.cd(2)
#outputfile.cd()
nEvents.Write()
nEvents_1.Write()
nEvents_2.Write()
nEvents_3.Write()
nEvents_4.Write()
nEvents_13.Write()

nEntries_1.Write()
nEntries_2.Write()
nEntries_7.Write()
nEntries_8.Write()
nEntries_9.Write()
nEntries_3.Write()
nEntries_4.Write()
nEntries_5.Write()
nEntries_6.Write()
nEntries_10.Write()
nEntries_11.Write()

nEvts_bin1.Write()
nEvts_bin2.Write()
nEvts_bin3.Write()
nEvts_bin4.Write()
nEvts_bin5.Write()
nEvts_bin6.Write()
nEvts_bin7.Write()
nEvts_bin8.Write()
nEvts_bin9.Write()
nEvts_bin10.Write()

nEvts_bin11.Write()
nEvts_bin12.Write()
nEvts_bin13.Write()
nEvts_bin14.Write()
nEvts_bin15.Write()
nEvts_bin16.Write()
nEvts_bin17.Write()
nEvts_bin18.Write()
nEvts_bin19.Write()
nEvts_bin20.Write()

nEvts_bin21.Write()
nEvts_bin22.Write()
nEvts_bin23.Write()
nEvts_bin24.Write()
nEvts_bin25.Write()
nEvts_bin26.Write()
nEvts_bin27.Write()
nEvts_bin28.Write()
nEvts_bin29.Write()
nEvts_bin30.Write()


tauPT_1.Write()
tauPT_2.Write()
metPT.Write()
HT_Tot.Write()
HT_Tot_1.Write()
ptratio_tau1.Write()
ptratio_tau2.Write()
genmatch_ptratio_tau1.Write()
genmatch_ptratio_tau2.Write()
notgenmatch_ptratio_tau1.Write()
notgenmatch_ptratio_tau2.Write()
MT.Write()
DR_daughter.Write()
DR_nr_genreco1.Write()
DR_nr_genreco2.Write()
t1_R_1030.Write()
t1_R_3050.Write()
t1_R_5080.Write()
t1_R_80130.Write()
t1_R_130180.Write()
t1_R_180230.Write()
t1_R_230280.Write()
t1_R_280330.Write()
t1_R_330380.Write()
t1_R_380500.Write()

t1_notgR_1030.Write()
t1_notgR_3050.Write()
t1_notgR_5080.Write()
t1_notgR_80130.Write()
t1_notgR_130180.Write()
t1_notgR_180230.Write()
t1_notgR_230280.Write()
t1_notgR_280330.Write()
t1_notgR_330380.Write()
t1_notgR_380500.Write()

t2_R_1030.Write()
t2_R_3050.Write()
t2_R_5080.Write()
t2_R_80130.Write()
t2_R_130180.Write()
t2_R_180230.Write()
t2_R_230280.Write()
t2_R_280330.Write()
t2_R_330380.Write()
t2_R_380500.Write()

t2_notgR_1030.Write()
t2_notgR_3050.Write()
t2_notgR_5080.Write()
t2_notgR_80130.Write()
t2_notgR_130180.Write()
t2_notgR_180230.Write()
t2_notgR_230280.Write()
t2_notgR_280330.Write()
t2_notgR_330380.Write()
t2_notgR_380500.Write()


print( nEvents.GetEntries())
print( nEvents_1.GetEntries())
print( nEvents_2.GetEntries())
print( nEntries_1.GetEntries())
print( nEntries_2.GetEntries())
print( nEntries_7.GetEntries())
print( nEntries_8.GetEntries())
print( nEntries_9.GetEntries())
print( nEntries_3.GetEntries())
print( nEntries_4.GetEntries())
print( nEntries_5.GetEntries())
print( nEntries_6.GetEntries())
print( nEntries_10.GetEntries())
print( nEntries_11.GetEntries())
print( tauPT_1.GetEntries())
print( tauPT_2.GetEntries())
print( metPT.GetEntries())
print( HT_Tot.GetEntries())
print( ptratio_tau1.GetEntries())
print( ptratio_tau2.GetEntries())
print( genmatch_ptratio_tau1.GetEntries())
print( genmatch_ptratio_tau2.GetEntries())
print( notgenmatch_ptratio_tau1.GetEntries())
print( notgenmatch_ptratio_tau2.GetEntries())
print( MT.GetEntries())
print( DR_daughter.GetEntries())
print( DR_nr_genreco1.GetEntries())
print( DR_nr_genreco2.GetEntries())

print( t1_R_3050.GetEntries())
print( t1_R_5080.GetEntries())
print( t1_R_80130.GetEntries())
print( t1_R_130180.GetEntries())
print( t1_R_180230.GetEntries())
print( t1_R_230280.GetEntries())
print( t1_R_280330.GetEntries())
print( t1_R_330380.GetEntries())
print( t1_R_380500.GetEntries())


print( t1_notgR_3050.GetEntries())
print( t1_notgR_5080.GetEntries())
print( t1_notgR_80130.GetEntries())
print( t1_notgR_130180.GetEntries())
print( t1_notgR_180230.GetEntries())
print( t1_notgR_230280.GetEntries())
print( t1_notgR_280330.GetEntries())
print( t1_notgR_330380.GetEntries())
print( t1_notgR_380500.GetEntries())

print( t2_R_3050.GetEntries())
print( t2_R_5080.GetEntries())
print( t2_R_80130.GetEntries())
print( t2_R_130180.GetEntries())
print( t2_R_180230.GetEntries())
print( t2_R_230280.GetEntries())
print( t2_R_280330.GetEntries())
print( t2_R_330380.GetEntries())
print( t2_R_380500.GetEntries())

print( t2_notgR_3050.GetEntries())
print( t2_notgR_5080.GetEntries())
print( t2_notgR_80130.GetEntries())
print( t2_notgR_130180.GetEntries())
print( t2_notgR_180230.GetEntries())
print( t2_notgR_230280.GetEntries())
print( t2_notgR_280330.GetEntries())
print( t2_notgR_330380.GetEntries())
print( t2_notgR_380500.GetEntries())

print( t1_R_1030.GetEntries())
print( t1_notgR_1030.GetEntries())
print( t2_R_1030.GetEntries())
print( t2_notgR_1030.GetEntries())



outputfile.Close()
input("Press Enter to continue...")


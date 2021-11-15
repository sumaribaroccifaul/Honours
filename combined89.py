import ROOT
import math
import array
from ROOT import TMath
from ROOT import Math
from ROOT import TLorentzVector
import numpy as np

if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_89GeV.root')
  nEntries = ntuple.GetEntries()

  m_mu = 0.1057 #Muon mass in GeV


  # create histograms 
  m_hist_mu = ROOT.TH1F('m_hist_mu','Muon invariant mass @ ecm 89 GeV;m_{mumu} [GeV];entries',250,85,92)
  doublehist_mu = ROOT.TH2F('doublehist_mu','jet invariant mass vs muon invariant mass after muon selection;m_{mumu} [GeV];m_{jets} [GeV]',250,85,92,250,85,92)

  m_hist_jet = ROOT.TH1F('m_hist_jet','Jet invariant mass @ ecm 89 GeV;m_{jet} [GeV];entries',250,85,92)
  doublehist_jet = ROOT.TH2F('doublehist_jet','jet invariant mass vs muon invariant mass after jet selection;m_{mumu} [GeV];m_{jets} [GeV]',250,85,92,250,85,92)

  m_hist_jetmu = ROOT.TH1F('m_hist_jetmu','Jet invariant mass after muon selection @ ecm 89 GeV;m_{jet} [GeV];entries',250,85,92)
  m_hist_mujet = ROOT.TH1F('m_hist_mujet','Muon invariant mass after jet selection @ ecm 89 GeV;m_{mumu} [GeV];entries',250,85,92)
  doublehist_both = ROOT.TH2F('doublehist_both','jet invariant mass vs muon invariant mass after muon and jet selection;m_{mumu} [GeV];m_{jets} [GeV]',250,85,92,250,85,92)

  deltaR_hist1 = ROOT.TH1F('deltaR_hist1',' deltaR between muon1 and jet1 ; deltaR ; entries ', 50, -0.001,0.005)
  deltaR_hist2 = ROOT.TH1F('deltaR_hist2',' detaR between muon2 and jet2 ; deltaR ; entries ', 50, -0.001,0.005)

  # loop over entries
  for entry in range(nEntries):
    entryCheck = ntuple.GetEntry(entry)
    if entryCheck <= 0:  continue

    if not ntuple.jets_px.size() == 2: continue
    if not ntuple.muons_px.size() == 2: continue
    
    m_z_mu = 0
    m_z_jet = 0 

    # loop over each muon in the event and create muon pairs
    for i in range(0,ntuple.muons_px.size()):

      for j in range(i+1,ntuple.muons_px.size()):

        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.muons_px[i]**2 + ntuple.muons_py[i]**2 + ntuple.muons_pz[i]**2)
        E_i = math.sqrt(m_mu**2 + P_i**2)

        P_j = math.sqrt(ntuple.muons_px[j]**2 + ntuple.muons_py[j]**2 + ntuple.muons_pz[j]**2)
        E_j = math.sqrt(m_mu**2 + P_j**2)

        #Reconstruct the momemtum 4-vector 
        # Vecmu1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.muons_px[i],ntuple.muons_py[i],ntuple.muons_pz[i],E_i)
        # Vecmu2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.muons_px[j],ntuple.muons_py[j],ntuple.muons_pz[j],E_j)
        Vecmu1 = ROOT.TLorentzVector(ntuple.muons_px[i],ntuple.muons_py[i],ntuple.muons_pz[i],E_i)
        Vecmu2 = ROOT.TLorentzVector(ntuple.muons_px[j],ntuple.muons_py[j],ntuple.muons_pz[j],E_j)
        sumVecmu = Vecmu1+Vecmu2

        m_z_mu = sumVecmu.M() 

        eta_mu1 = Vecmu1.Eta()
        eta_mu2 = Vecmu2.Eta()
        phi_mu1 = Vecmu1.Phi()
        phi_mu2 = Vecmu2.Phi()
        
    # loop over each jet in the event and create jet pairs
    for i in range(0,ntuple.jets_px.size()):

      for j in range(i+1,ntuple.jets_px.size()):

        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.jets_px[i]**2 + ntuple.jets_py[i]**2 + ntuple.jets_pz[i]**2)
        E_i = math.sqrt(ntuple.jets_M[i]**2 + P_i**2)

        P_j = math.sqrt(ntuple.jets_px[j]**2 + ntuple.jets_py[j]**2 + ntuple.jets_pz[j]**2)
        E_j = math.sqrt(ntuple.jets_M[j]**2 + P_j**2)

        #Reconstruct the momemtum 4-vector 
        # Vecjet1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[i],ntuple.jets_py[i],ntuple.jets_pz[i],E_i)
        # Vecjet2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[j],ntuple.jets_py[j],ntuple.jets_pz[j],E_j)
        Vecjet1 = ROOT.TLorentzVector(ntuple.jets_px[i],ntuple.jets_py[i],ntuple.jets_pz[i],E_i)
        Vecjet2 = ROOT.TLorentzVector(ntuple.jets_px[j],ntuple.jets_py[j],ntuple.jets_pz[j],E_j)

        sumVecjet = Vecjet1+Vecjet2

        m_z_jet = sumVecjet.M()

        if m_z_mu >= 85 and m_z_mu <= 92: 
          m_hist_mu.Fill(m_z_mu) # mass cut 
          doublehist_mu.Fill(m_z_mu, m_z_jet)

          if m_z_jet >= 85 and m_z_jet <= 92:
            m_hist_jetmu.Fill(m_z_jet)
            m_hist_mujet.Fill(m_z_mu)
            doublehist_both.Fill(m_z_mu, m_z_jet)

        if m_z_jet >= 85 and m_z_jet <= 92: 
          m_hist_jet.Fill(m_z_jet) # mass cut
          doublehist_jet.Fill(m_z_mu, m_z_jet)

        eta_jet1 = Vecjet1.Eta()
        phi_jet1 = Vecjet1.Phi()
        eta_jet2 = Vecjet2.Eta()
        phi_jet2 = Vecjet2.Phi()

        deltaR_1 = Vecmu1.DeltaR(Vecjet1)
        deltaR_2 = Vecmu2.DeltaR(Vecjet2)
        
        if ( Vecmu1.Phi()*Vecjet1.Phi() ) > 0:
          deltaR_hist1.Fill(deltaR_1)
        
        if ( Vecmu2.Phi()*Vecjet2.Phi() ) > 0:
          deltaR_hist2.Fill(deltaR_2)
        
  outfile = ROOT.TFile('output_combined89.root','recreate')

  m_hist_mu.Write()
  doublehist_mu.Write()

  m_hist_jet.Write()
  doublehist_jet.Write()

  m_hist_jetmu.Write()
  m_hist_mujet.Write()
  doublehist_both.Write()

  deltaR_hist1.Write()
  deltaR_hist2.Write()

  outfile.Close()

  can = ROOT.TCanvas("can", "mu", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist_mu.Draw()
  can.SaveAs("89plots/mu89.pdf")

  can2 = ROOT.TCanvas("can2", "jet", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist_jet.Draw()
  can2.SaveAs("89plots/jet89.pdf")

  can3 = ROOT.TCanvas("can3", "jetmu", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist_jetmu.Draw()
  can3.SaveAs("89plots/jetmu89.pdf")

  can4 = ROOT.TCanvas("can4", "mujet", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist_mujet.Draw()
  can4.SaveAs("89plots/mujet89.pdf")

  can5 = ROOT.TCanvas("can5", "doublemu", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  doublehist_mu.Draw()
  can5.SaveAs("89plots/doublehistmu89.pdf")

  can6 = ROOT.TCanvas("can6", "doublejet", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  doublehist_jet.Draw()
  can6.SaveAs("89plots/doublehistjet89.pdf")

  can7 = ROOT.TCanvas("can7", "doubleboth", 400, 400, 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  doublehist_both.Draw()
  can7.SaveAs("89plots/doublehistboth89.pdf")

  can8 = ROOT.TCanvas("can8", "deltaR_hist1", 200, 300, 200, 300)
  ROOT.gStyle.SetOptStat("eM")
  deltaR_hist1.Draw()
  can8.SaveAs("89plots/overlap/deltaR_hist1.pdf")

  can9 = ROOT.TCanvas("can9", "deltaR_hist2", 200, 300, 200, 300)
  ROOT.gStyle.SetOptStat("eM")
  deltaR_hist2.Draw()
  can9.SaveAs("89plots/overlap/deltaR_hist2.pdf")




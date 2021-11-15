import ROOT
import math
import array

#########Main############
if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_97GeV.root')
  nEntries = ntuple.GetEntries()

  m_mu = 0.1057 #Muon mass in GeV

  #Create some ROOT histograms: MUONS
  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [GeV];entries/bin',100,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [GeV];entries/bin',100,0,200)#Momentum

  m_hist = ROOT.TH1F('m_hist','Invariant mass for ee -> Z -> mumu events generated at 96 GeV;m_{mumu} [GeV];entries',180,88,100)
  etaplot = ROOT.TH2F ('eta', '97 GeV: eta_{jet_1} vs eta_{jet_2}; eta_{jet_1}; eta_{jet_2}',50,-2,2,50,-2,2)
  phiplot = ROOT.TH2F ('phi', '97 GeV:phi_{jet_1} vs phi_{jet_2}; phi_{jet_1}; phi_{jet_2}',50,-4,4,50,-4,4)
  etaphi1plot = ROOT.TH2F ('etaphi1', '97 GeV:eta_{jet_1} vs phi_{jet_1}; eta_{jet_1}; phi_{jet_1}',50,-4,4,50,-4,4)
  etaphi2plot = ROOT.TH2F ('etaphi2', '97 GeV:eta_{jet_2} vs phi_{jet_2}; eta_{jet_2}; phi_{jet_2}',50,-4,4,50,-4,4)
  detadphi = ROOT.TH2F ('detadphi', '97 GeV: dEta vs dPhi (all entries)',50,-4,4,50,-4,4)

  N=0

  #Loop over each entry in the ntuple
  for entry in range( nEntries ):
   # give the user something to look at...
    #if entry%1000 == 0:
    #  print(entry)

    # check that the event is read properly
    entryCheck = ntuple.GetEntry( entry )
    if entryCheck <= 0:  continue

    #require exactly 2 muons
    if not ntuple.muons_px.size() == 2 : continue

    # loop over each muon in the event and create muon pairs
    for i in range(0,ntuple.muons_px.size()):

      for j in range(i+1,ntuple.muons_px.size()):
        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.muons_px[i]**2 + ntuple.muons_py[i]**2 + ntuple.muons_pz[i]**2)
        E_i = math.sqrt(m_mu**2 + P_i**2)

        P_j = math.sqrt(ntuple.muons_px[j]**2 + ntuple.muons_py[j]**2 + ntuple.muons_pz[j]**2)
        E_j = math.sqrt(m_mu**2 + P_j**2)

        #Reconstruct the momemtum 4-vector 
        Vecmu1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.muons_px[i],ntuple.muons_py[i],ntuple.muons_pz[i],E_i)
        Vecmu2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.muons_px[j],ntuple.muons_py[j],ntuple.muons_pz[j],E_j)
        sumVecmu = Vecmu1+Vecmu2

        m_z = sumVecmu.M()

        if abs(m_z-97.0)<1.0: N=N+1

        eta1 = Vecmu1.Eta()
        phi1 = Vecmu1.Phi()
        if phi1 < 0: phi1 = phi1 + 2*math.pi

        eta2 = Vecmu2.Eta()
        phi2 = Vecmu2.Phi()
        if phi2 < 0: phi2 = phi2 + 2*math.pi

        dEta = eta1 - eta2
        dPhi = phi1 - phi2

        if m_z >= 88 and m_z <= 100: 
          m_hist.Fill(m_z)
        
        etaplot.Fill(eta1,eta2)
        phiplot.Fill(phi1,phi2)
        etaphi1plot.Fill(eta1,phi1)
        etaphi2plot.Fill(eta2,phi2)
        detadphi.Fill(dEta,dPhi)
#        Ei_hist.Fill(E_i)
  print("NUMBER OF ENTRIES IN Z-PEAK:",N)

  gaus1 = " [2]*exp(-0.5*((x-[1])/[0])*((x-[1])/[0]))  "
  gaus2 = " [2]*exp(-0.5*((x-[1])/[0])*((x-[1])/[0])) "

  gausFit1 = ROOT.TF1("gausFit1",gaus1,89.25,92.75)
  gausFit1.SetParameters(1,91,760)
  gausFit2 = ROOT.TF1("gausFit2",gaus2,96.3,99)
  gausFit2.SetParameters(0.5,97,2000)
# Histograms to ROOT file
  outfile = ROOT.TFile('output_2muons97.root','recreate')

  Ei_hist.Write()
  Pi_hist.Write()
  m_hist.Write()

  outfile.Close()

  #Save output
  can = ROOT.TCanvas("can", "histograms", 600, 600)
  ROOT.gStyle.SetOptStat("eM")
  m_hist.Draw("pe")
  m_hist.Fit(gausFit1,"R")
  m_hist.Fit(gausFit2,"+","",96.3,99)
  chi21 = gausFit1.GetChisquare()
  ndf1 = gausFit1.GetNDF()
  chiperdof1 = chi21/ndf1
  print("Chi square per degrees of freedom: ",chiperdof1)
  chi22 = gausFit2.GetChisquare()
  ndf2 = gausFit2.GetNDF()
  chiperdof2 = chi22/ndf2
  print("Chi square per degrees of freedom: ",chiperdof2)
  can.SaveAs("97plots/2muons97.pdf")

  can2 = ROOT.TCanvas("can2","eta",400,400)
  etaplot.Draw()
  can2.SaveAs("97plots/eta_mm97.pdf")

  can3 = ROOT.TCanvas("can3","phi",400,400)
  phiplot.Draw()
  can3.SaveAs("97plots/phi_mm97.pdf")

  can4 = ROOT.TCanvas("can4","etaphi1plot",400,400)
  etaphi1plot.Draw()
  can4.SaveAs("97plots/etaphi1_mm97.pdf")

  can5 = ROOT.TCanvas("can5","etaphi2plot",400,400)
  etaphi2plot.Draw()
  can5.SaveAs("97plots/etaphi2_mm97.pdf")

  # can6 = ROOT.TCanvas("can6","detadphilow",400,400)
  # detadphiLow.Draw()
  # can6.SaveAs("97plots/detadphiLOW_mm97.pdf")

  # can7 = ROOT.TCanvas("can7","detadphihigh",400,400)
  # detadphiHigh.Draw()
  # can7.SaveAs("97plots/detadphiHIGH_mm97.pdf")

  can6 = ROOT.TCanvas("can6","detadphi",400,400)
  detadphi.Draw()
  can6.SaveAs("97plots/detadphi_mm97.pdf")

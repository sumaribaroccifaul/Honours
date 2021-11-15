import ROOT
import math
import array
import numpy as np

#########Main############
if __name__=="__main__":

  # load the ROOT ntuple  
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_97GeV.root')
  nEntries = ntuple.GetEntries()
  #print 'There are',nEntries,'entries in your ntuple'

  #Create some ROOT histograms 
  jets_hist = ROOT.TH1F('jets_hist', 'Number of jets per entry; number of jets @ 97 GeV/c^2 ; entries',5,0,5)
  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [GeV];entries/bin',100,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [GeV];entries/bin',100,0,200)#Momentum

  m_hist = ROOT.TH1F('m_hist','Invariant mass for ee->Z->jets events generated at 97 GeV;m_{jj} [GeV];entries/bin',180,83,100)

  etaplot = ROOT.TH2F ('eta', '97 GeV: eta_{jet_1} vs eta_{jet_2}; eta_{jet_1}; eta_{jet_2}',50,-2,2,50,-2,2)
  phiplot = ROOT.TH2F ('phi', '97 GeV:phi_{jet_1} vs phi_{jet_2}; phi_{jet_1}; phi_{jet_2}',50,-7,7,50,-7,7)
  etaphi1plot = ROOT.TH2F ('etaphi1', '97 GeV:eta_{jet_1} vs phi_{jet_1}; eta_{jet_1}; phi_{jet_1}',50,-2,2,50,-1,7)
  etaphi2plot = ROOT.TH2F ('etaphi2', '97 GeV:eta_{jet_2} vs phi_{jet_2}; eta_{jet_2}; phi_{jet_2}',50,-2,2,50,-1,7)

  detadphiLow = ROOT.TH2F ('detadphilow', '97 GeV: dEta vs dPhi (low peak)',50,-4,4,50,-4,4)
  detadphiHigh = ROOT.TH2F ('detadphihigh', '97 GeV: dEta vs dPhi (high peak)',50,-4,4,50,-4,4)
  detadphi = ROOT.TH2F ('detadphi', '97 GeV: dEta vs dPhi (all entries)',50,-4,4,50,-4,4)

  #m_hist_manual = ROOT.TH1F('m_hist_manual','m;m_{mm} [GeV];entries/bin',100,0,200)
  N=0

  #Loop over each entry in the ntuple
  for entry in range( nEntries ):
   # give the user something to look at...
  #  if entry%1000 == 0:
   #   print(entry)

    # check that the event is read properly
    entryCheck = ntuple.GetEntry( entry )
    if entryCheck <= 0:  continue
    
    jets = ntuple.jets_px.size()
    jets_hist.Fill(jets)

    #require exactly 2 jets
    if not ntuple.jets_px.size() == 2: continue
    # loop over each muon in the event and create muon pairs
    for i in range(0,ntuple.jets_px.size()):
      for j in range(i+1,ntuple.jets_px.size()):
        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.jets_px[i]**2 + ntuple.jets_py[i]**2 + ntuple.jets_pz[i]**2)
        E_i = math.sqrt(ntuple.jets_M[i]**2 + P_i**2)

        P_j = math.sqrt(ntuple.jets_px[j]**2 + ntuple.jets_py[j]**2 + ntuple.jets_pz[j]**2)
        E_j = math.sqrt(ntuple.jets_M[j]**2 + P_j**2)

        Vecjet1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[i],ntuple.jets_py[i],ntuple.jets_pz[i],E_i)
        Vecjet2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[j],ntuple.jets_py[j],ntuple.jets_pz[j],E_j)
        sumVecjet = Vecjet1+Vecjet2
        m_z = sumVecjet.M() #This works perfectly

        eta1 = Vecjet1.Eta()
        phi1 = Vecjet1.Phi()
        #if phi1 > 0: phi1 = phi1 - 2*math.pi

        eta2 = Vecjet2.Eta()
        phi2 = Vecjet2.Phi()
        
        #if phi2 < 0: phi1 = phi1 - 2*math.pi # this gives an empty phi vs phi plot

        dEta = eta1 - eta2
        dPhi = phi1 - phi2

        if abs(m_z-97.0)<1.0: N=N+1 #N is the number of entries to be used in the cross section calculation 

        if m_z >= 83 and m_z <= 100:
          m_hist.Fill(m_z)

        if m_z >= 90.2 and m_z <= 92.2:
          detadphiLow.Fill(dEta,dPhi)
        
        if m_z >= 96.5 and m_z <= 97.5:
         detadphiHigh.Fill(dEta,dPhi)

        
        Ei_hist.Fill(E_i)
        Pi_hist.Fill(P_i)

        etaplot.Fill(eta1,eta2)
        phiplot.Fill(phi1,phi2)
        etaphi1plot.Fill(eta1,phi1)
        etaphi2plot.Fill(eta2,phi2)
        detadphi.Fill(dEta,dPhi)

        # # calculate dEta and dPhi for entries for a narrow range around the two peaks
        # lowpeakentries = np.where(m_z >= 90.2 and m_z <= 92.2)

        # for l in range(0,len(lowpeakentries) - 1):
        #   for k in range(l+1, len(lowpeakentries) - 1):
        #     P_l = math.sqrt(ntuple.jets_px[l]**2 + ntuple.jets_py[l]**2 + ntuple.jets_pz[l]**2)
        #     E_l = math.sqrt(ntuple.jets_M[l]**2 + P_l**2)

        #     P_k = math.sqrt(ntuple.jets_px[k]**2 + ntuple.jets_py[k]**2 + ntuple.jets_pz[k]**2)
        #     E_k = math.sqrt(ntuple.jets_M[k]**2 + P_k**2)

        #     fourvecJet1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[l],ntuple.jets_py[l],ntuple.jets_pz[l],E_l)
        #     fourvecJet2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.jets_px[k],ntuple.jets_py[k],ntuple.jets_pz[k],E_k)
        #     totalvec = fourvecJet1+fourvecJet2

        #     eta1 = fourvecJet1.Eta()
        #     phi1 = fourvecJet1.Phi()
        #     eta2 = fourvecJet2.Eta()
        #     phi2 = fourvecJet2.Phi()

        #     dEta = eta1 - eta2
        #     dPhi = phi1 - phi2

        #     detadphiLow.Fill(dEta,dPhi)

  print("number of entries:",N)

  gaus1 = " [2]*exp(-0.5*((x-[1])/[0])*((x-[1])/[0]))  "
  gaus2 = " [2]*exp(-0.5*((x-[1])/[0])*((x-[1])/[0])) "
  #gausSum = "[4]*(%s) + [5]*(%s)"%(gaus1,gaus2)
  #gausFit = ROOT.TF1("gausFit",gausSum,88,97)
  #gausFit.SetParameters(1,91,0.5,95,700,2000)
  #gausFit.SetParameters(1.58221e+00, 9.17909e+01,4.58308e-01,9.48713e+01,6.67719e+02,2.21456e+03)
  #gausFit.SetParLimits(1, 91.6, 91.82)
  #gausFit.SetParLimits(3, 94,96) 

  gausFit1 = ROOT.TF1("gausFit1",gaus1,89,92.5)
  gausFit1.SetParameters(1,91,760)
  gausFit2 = ROOT.TF1("gausFit2",gaus2,96.2,99)
  gausFit2.SetParameters(0.5,95,2000)

  # BW fit: function and setting param initial guesses 
  bw_A = "2*sqrt(2)*[0]*[1]*sqrt([0]*[0]*([0]*[0]+[1]*[1]))"
  bw_B = "3.14159*sqrt([0]*[0] + sqrt([0]*[0]*([0]*[0]+[1]*[1])))" 
  bw_C = "(x*x-[0]*[0])*(x*x-[0]*[0]) + [0]*[0]*[1]*[1]"
  bw = "[2]*((%s)/(%s))/(%s)"%(bw_A, bw_B, bw_C)

  bwFit = ROOT.TF1("bwfit",bw,85,93) 

  bwFit.SetParameter(0,91) 
  bwFit.SetParameter(1,1)
  bwFit.SetParameter(2,3000)
  bwFit.SetLineColor(6)

#Write the histograms to a ROOT file
  outfile = ROOT.TFile('outputjj97.root','recreate')
  Ei_hist.Write()
  Pi_hist.Write()
  #m_hist_manual.Write()
  m_hist.Write()
  jets_hist.Write()
  outfile.Close()

  can = ROOT.TCanvas("can", "histograms", 800, 800)
  ROOT.gStyle.SetOptStat("eM")

  m_hist.Draw("pe")

  m_hist.Fit(gausFit1,"R")
  m_hist.Fit(gausFit2,"+","",96.2,99)
  m_hist.Fit(bwFit,"+","",85,93)

  chi21 = gausFit1.GetChisquare()
  ndf1 = gausFit1.GetNDF()
  chiperdof1 = chi21/ndf1
  print("Chi square per degrees of freedom: ",chiperdof1)
  
  chi22 = gausFit2.GetChisquare()
  ndf2 = gausFit2.GetNDF()
  chiperdof2 = chi22/ndf2
  print("Chi square per degrees of freedom: ",chiperdof2)

  chi23 = bwFit.GetChisquare()
  ndf3 = bwFit.GetNDF()
  chiperdof3 = chi23/ndf3
  print("Chi square per degrees of freedom: ",chiperdof3)

  can.SaveAs("97plots/2jets97.pdf")

  can2 = ROOT.TCanvas("can2","eta",400,400)
  etaplot.Draw()
  can2.SaveAs("97plots/eta_jj97.pdf")

  can3 = ROOT.TCanvas("can3","phi",400,400)
  phiplot.Draw()
  can3.SaveAs("97plots/phi_jj97.pdf")

  can4 = ROOT.TCanvas("can4","etaphi1plot",400,400)
  etaphi1plot.Draw()
  can4.SaveAs("97plots/etaphi1_jj97.pdf")

  can5 = ROOT.TCanvas("can5","etaphi2plot",400,400)
  etaphi2plot.Draw()
  can5.SaveAs("97plots/etaphi2_jj97.pdf")

  can6 = ROOT.TCanvas("can6","detadphilow",400,400)
  detadphiLow.Draw()
  can6.SaveAs("97plots/detadphiLOW_jj97.pdf")

  can7 = ROOT.TCanvas("can7","detadphihigh",400,400)
  detadphiHigh.Draw()
  can7.SaveAs("97plots/detadphiHIGH_jj97.pdf")

  can8 = ROOT.TCanvas("can8","detadphi",400,400)
  detadphi.Draw()
  can8.SaveAs("97plots/detadphi_jj97.pdf")
#  can2 = ROOT.TCanvas("can2", "jets histogram", 800,800)
#  ROOT.gStyle.SetOptStat("eM")
#  jets_hist.Draw()
#  can2.SaveAs("jetshist97.pdf") 

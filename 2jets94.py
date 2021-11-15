import ROOT
import math
import array

if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_94GeV.root')
  nEntries = ntuple.GetEntries()
  #print 'There are',nEntries,'entries in your ntuple'

  #Create some ROOT histograms 
  jets_hist = ROOT.TH1F('jets_hist', 'Number of jets per entry @ 94 GeV/c^2 ; number of jets; entries',5,0,5) # number of jets 

  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [GeV];entries/bin',100,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [GeV];entries/bin',100,0,200)#Momentum
 
  m_hist = ROOT.TH1F('m_hist','Invariant mass for ee->Z->jets events generated at 94 GeV;m_{jj} [GeV];entries/bin',180,83,97)
  etaplot = ROOT.TH2F ('eta', 'Jet 1 eta vs jet 2 eta; eta1; eta2',50,-2,2,50,-2,2)
  phiplot = ROOT.TH2F ('phi', 'Jet 1 phi vs jet 2 phi; phi1; phi2',50,-4,4,50,-4,4)
  etaphi1plot = ROOT.TH2F ('etaphi1', 'Jet 1: Eta vs Phi; eta1; phi1',50,-4,4,50,-4,4)
  etaphi2plot = ROOT.TH2F ('etaphi2', 'Jet 2: Eta vs Phi; eta2; phi2',50,-4,4,50,-4,4)
  detadphi = ROOT.TH2F ('detadphi', 'dEta vs dPhi',50,-4,4,50,-4,4)
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
        eta2 = Vecjet2.Eta()
        phi2 = Vecjet2.Phi()

        dEta = eta1 - eta2
        dPhi = phi1 - phi2

        if abs(m_z-94.0)<1.0: N=N+1 #N is the number of entries to be used in the crosss section calculation 

        if m_z >= 83 and m_z <= 97: m_hist.Fill(m_z)
        #m_hist.Fill(m_z)
        Ei_hist.Fill(E_i)
        Pi_hist.Fill(P_i)
        etaplot.Fill(eta1,eta2)
        phiplot.Fill(phi1,phi2)
        etaphi1plot.Fill(eta1,phi1)
        etaphi2plot.Fill(eta2,phi2)
        detadphi.Fill(dEta,dPhi)

  print("NUMBER OF ENTRIES IN Z-PEAK:",N)

  # Do the Double Gauss fit 
  gaus1 = " exp(-0.5*((x-[1])/[0])*((x-[1])/[0]))  "
  gaus2 = " exp(-0.5*((x-[3])/[2])*((x-[3])/[2])) "
  gausSum = "[4]*(%s) + [5]*(%s)"%(gaus1,gaus2)
  gausFit = ROOT.TF1("gausFit",gausSum,88.5,96)
  gausFit.SetParameters(1,91,0.5,93,600,2000)


 #Write the histograms to a ROOT file
  outfile = ROOT.TFile('output_2jets94.root','recreate')
  jets_hist.Write()
  Ei_hist.Write()
  Pi_hist.Write()
  #m_hist_manual.Write()
  m_hist.Write()
  outfile.Close()

  can = ROOT.TCanvas("can", "histograms", 800, 800)
  ROOT.gStyle.SetOptStat("eMou")

  m_hist.Draw("pe")
  #m_hist.Fit("fDoubleGauss","M", 83,96)
 # m_hist.Fit(fDoubleGauss)
  #m_hist.Fit("gaus")
  m_hist.Fit(gausFit,"R")
  chi2 = gausFit.GetChisquare()
  ndf = gausFit.GetNDF()
  chiperdof = chi2/ndf
  print("Chi square per degrees of freedom: ",chiperdof)
  can.SaveAs("94plots/2jets94.pdf")

  can2 = ROOT.TCanvas("can2","eta",400,400)
  etaplot.Draw()
  can2.SaveAs("94plots/eta_jj94.pdf")

  can3 = ROOT.TCanvas("can3","phi",400,400)
  phiplot.Draw()
  can3.SaveAs("94plots/phi_jj94.pdf")

  can4 = ROOT.TCanvas("can4","etaphi1plot",400,400)
  etaphi1plot.Draw()
  can4.SaveAs("94plots/etaphi1_jj94.pdf")

  can5 = ROOT.TCanvas("can5","etaphi2plot",400,400)
  etaphi2plot.Draw()
  can5.SaveAs("94plots/etaphi2_jj94.pdf")

  can6 = ROOT.TCanvas("can6","detadphi",400,400)
  detadphi.Draw()
  can6.SaveAs("94plots/detadphi_jj94.pdf")

  #can2 = ROOT.TCanvas("can2", "jets histogram", 800,800)
  #ROOT.gStyle.SetOptStat("eM")
  #jets_hist.Draw()
  #can2.SaveAs("jetshist94.pdf") 

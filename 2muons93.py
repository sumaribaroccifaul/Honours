import ROOT
import math
import array

#########Main############
if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_93GeV.root')
  nEntries = ntuple.GetEntries()

  m_mu = 0.1057 #Muon mass in GeV

  #Create some ROOT histograms: MUONS
  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [GeV];entries/bin',100,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [GeV];entries/bin',100,0,200)#Momentum

  m_hist = ROOT.TH1F('m_hist','Invariant mass for ee -> Z -> mumu events generated at 93 GeV;m_{mumu} [GeV];entries',180,87,95)

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

        if m_z >= 87 and m_z <= 95: m_hist.Fill(m_z)

        if abs(m_z-93.0)<1.0: N=N+1

  print("NUMBER OF ENTRIES IN Z-PEAK:",N)

# Fitting: 
  gaus1 = " exp(-0.5*((x-[1])/[0])*((x-[1])/[0]))  "
  gaus2 = " exp(-0.5*((x-[3])/[2])*((x-[3])/[2])) "
  gausSum = "[4]*(%s) + [5]*(%s)"%(gaus1,gaus2)
  gausFit = ROOT.TF1("gausFit",gausSum,90,95)
  gausFit.SetParameters(1,91,0.5,93,600,2000)

# Histograms to ROOT file
  outfile = ROOT.TFile('output_2muons93.root','recreate')

 # Ei_hist.Write()
 # Pi_hist.Write()
  m_hist.Write()

  outfile.Close()

  #Save output
  can = ROOT.TCanvas("can", "histograms", 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist.Draw("pe")
  m_hist.Fit(gausFit,"R")

  chi2 = gausFit.GetChisquare()
  ndf = gausFit.GetNDF()
  chiperdof = chi2/ndf
  print("Chi square per degrees of freedom: ",chiperdof)
  can.SaveAs("2muons93.pdf")

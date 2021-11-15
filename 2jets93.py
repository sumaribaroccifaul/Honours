import ROOT
import math
import array
import numpy as np

#########Main############
if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZX_93GeV.root')
  nEntries = ntuple.GetEntries()
  
  #Create ROOT histograms 
  jets_hist = ROOT.TH1F('jets_hist', 'Number of jets per entry @ 93 GeV ; number of jets; entries',10,0,10)

  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [GeV];entries/bin',100,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [GeV];entries/bin',100,0,200)#Momentum

  m_hist = ROOT.TH1F('m_hist','Invariant mass for ee->Z->jets events generated at 93 GeV;m_{jj} [GeV];entries/bin',100,75,107)
  m_hist_gen = ROOT.TH1F('m_hist_gen','gen_particle invariant mass at \sqrt{s} = 93 GeV ;m_{jj} [GeV];entries/bin',100,75,107)

  etaplot = ROOT.TH2F ('eta', 'Jet 1 eta vs jet 2 eta; eta1; eta2',50,-2,2,50,-2,2)
  phiplot = ROOT.TH2F ('phi', 'Jet 1 phi vs jet 2 phi; phi1; phi2',50,-4,4,50,-4,4)
  etaphi1plot = ROOT.TH2F ('etaphi1', 'Jet 1: Eta vs Phi; eta1; phi1',50,-4,4,50,-4,4)
  etaphi2plot = ROOT.TH2F ('etaphi2', 'Jet 2: Eta vs Phi; eta2; phi2',50,-4,4,50,-4,4)
  detadphi = ROOT.TH2F ('detadphi', 'dEta vs dPhi',50,-4,4,50,-4,4)
  
  N=0
  #pt= np.zeros(nEntries)
  #pt_values = []

  #Loop over each entry in the ntuple
  for entry in range( nEntries ):

    # check that the event is read properly
    entryCheck = ntuple.GetEntry( entry )
    if entryCheck <= 0:  continue

    # distribution of nr jets in entries
    jets = ntuple.jets_px.size()
    jets_hist.Fill(jets)

    # calculate pt for each jet in event
    # for jet in range(0,ntuple.jets_px.size()):
    #   pt[jet] = math.sqrt(ntuple.jets_px[jet]**2 + ntuple.jets_py[jet]**2 + ntuple.jets_pz[jet]**2)

    # # find index (which jet) at which pt is highest and second highest
    # maximum_pt_1 = np.where(pt==max(pt))
    # maximum_pt_2 = np.where(pt==max(pt, key = lambda x: min(pt) - 1 if (x == maximum_pt_1) else x))

    #print("Index of jet with highest pt:",maximum_pt_1)
    #print("Index of jet with second highest pt:",maximum_pt_2)

    #require exactly 2 jets
    #if not ntuple.jets_px.size() == 2: continue

    #require >= 2 jets
    if ntuple.jets_px.size() < 2: continue

    # loop over each jet in the event and create jet pairs
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

        if abs(m_z-93.0)<1.0: N=N+1 #N is the number of entries to be used in the crosss section calculation 

        if m_z >= 75 and m_z <= 107: m_hist.Fill(m_z)
        #m_hist.Fill(m_z)
        Ei_hist.Fill(E_i)
        Pi_hist.Fill(P_i)
        etaplot.Fill(eta1,eta2)
        phiplot.Fill(phi1,phi2)
        etaphi1plot.Fill(eta1,phi1)
        etaphi2plot.Fill(eta2,phi2)
        detadphi.Fill(dEta,dPhi)

    # loop over genparticles:
    for i in range(0,ntuple.gen_particle_px.size()):
      for j in range(i+1,ntuple.gen_particle_px.size()):
        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.gen_particle_px[i]**2 + ntuple.gen_particle_py[i]**2 + ntuple.gen_particle_pz[i]**2)
        E_i = math.sqrt(ntuple.gen_particle_m[i]**2 + P_i**2)

        P_j = math.sqrt(ntuple.gen_particle_px[j]**2 + ntuple.gen_particle_py[j]**2 + ntuple.gen_particle_pz[j]**2)
        E_j = math.sqrt(ntuple.gen_particle_m[j]**2 + P_j**2)

        Vec1 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.gen_particle_px[i],ntuple.gen_particle_py[i],ntuple.gen_particle_pz[i],E_i)
        Vec2 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(ntuple.gen_particle_px[j],ntuple.gen_particle_py[j],ntuple.gen_particle_pz[j],E_j)
        sumVec = Vec1+Vec2

        #if (abs(ntuple.gen_particle_pdgid) < 10) and (abs(ntuple.gen_particle_pdgid) > 30):
        m_z_gen = sumVec.M()
        m_hist_gen.Fill(m_z_gen)

  #print("NUMBER OF ENTRIES IN Z-PEAK:",N)

  gaus1 = " exp(-0.5*((x-[1])/[0])*((x-[1])/[0]))  "
  gaus2 = " exp(-0.5*((x-[3])/[2])*((x-[3])/[2])) "
  gausSum = "[4]*(%s) + [5]*(%s)"%(gaus1,gaus2)
  gausFit = ROOT.TF1("gausFit",gausSum,87,107)
  gausFit.SetParameters(1,91,0.5,93,600,2000)
  #gausFit.SetParLimits(1, 91, 91.5)
  gausOnepeak = " [2]*(exp(-0.5*((x-[1])/[0])*((x-[1])/[0])))  "
  gausOnepeakfit = ROOT.TF1("gausOnepeakfit",gausOnepeak,88,107)
  gausOnepeakfit.SetParameters(1,93,1000)


  # fitting a crystal ball function to the mass hist:
  # def FitFunction(m_z, par):

  #   x = m_z[0] 
  #   alpha = par[0]
  #   n = par[1]
  #   sigma = par[2]
  #   mean = par[3]
    
  #   y = ROOT.Math.crystalball_function(x, alpha, n, sigma, mean)

  #   return y
  
  # cbFit = ROOT.TF1("cbFit",FitFunction,75,107,4)
  # cbFit.SetParameters(1000,1,1,95)
  # cbFit.SetLineColor(6)


 #Write the histograms to a ROOT file
  outfile = ROOT.TFile('outputjj93.root','recreate')
  Ei_hist.Write()
  Pi_hist.Write()
  m_hist.Write()
  jets_hist.Write()
  m_hist_gen.Write()
  outfile.Close()

  can = ROOT.TCanvas("can", "histograms", 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist.Draw("pe")
  m_hist.Fit(gausOnepeakfit,"R")
  chi2 = gausOnepeakfit.GetChisquare()
  ndf = gausOnepeakfit.GetNDF()
  chiperdof = chi2/ndf
  print("Chi square per degrees of freedom: ",chiperdof)

  # m_hist.Fit(cbFit,"+","",75,107)
 
  # chi2_cbFit = cbFit.GetChisquare() 
  # ndof_cbFit = cbFit.GetNDF()
  # print(chi2_cbFit,ndof_cbFit,chi2_cbFit/ndof_cbFit)
  can.SaveAs("93plots/2jets93.pdf")

  can_gen = ROOT.TCanvas("can_gen", "histograms", 400, 400)
  ROOT.gStyle.SetOptStat("eM")
  m_hist_gen.Draw("pe")
  m_hist_gen.Fit(gausOnepeakfit,"R")
  chi2gen = gausOnepeakfit.GetChisquare()
  ndfgen = gausOnepeakfit.GetNDF()
  chiperdofgen = chi2gen/ndfgen
  print("Chi square per degrees of freedom (gen particles): ",chiperdofgen)

  can_gen.SaveAs("93plots/2jets93_gen.pdf")

  can2 = ROOT.TCanvas("can2","eta",400,400)
  etaplot.Draw()
  can2.SaveAs("93plots/eta_jj93.pdf")

  can3 = ROOT.TCanvas("can3","phi",400,400)
  phiplot.Draw()
  can3.SaveAs("93plots/phi_jj93.pdf")

  can4 = ROOT.TCanvas("can4","etaphi1plot",400,400)
  etaphi1plot.Draw()
  can4.SaveAs("93plots/etaphi1_jj93.pdf")

  can5 = ROOT.TCanvas("can5","etaphi2plot",400,400)
  etaphi2plot.Draw()
  can5.SaveAs("93plots/etaphi2_jj93.pdf")

  can6 = ROOT.TCanvas("can6","detadphi",400,400)
  detadphi.Draw()
  can6.SaveAs("93plots/detadphi_jj93.pdf")

  can7 = ROOT.TCanvas("can7", "jets histogram", 400,400)
  ROOT.gStyle.SetOptStat("eM")
  jets_hist.Draw()
  can7.SaveAs("Njets93.pdf")

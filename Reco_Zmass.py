import ROOT
import math
import array



#########Main############
if __name__=="__main__":

  # load the ROOT ntuple
  ntuple = ROOT.TChain('events')
  ntuple.Add('FCCee_ZH.root')
  nEntries = ntuple.GetEntries()
  print 'There are',nEntries,'entries in your ntuple'

  m_mu = 105.7 #Muon mass in MeV

  #Create some ROOT histograms 
  Ei_hist = ROOT.TH1F('ei_hist','m;E_{i} [MeV];entries/bin',50,0,200)#Energy
  Pi_hist = ROOT.TH1F('pi_hist','m;P_{i} [MeV];entries/bin',50,0,200)#Momentum

  #Loop over each entry in the ntuple
  for entry in range( nEntries ):
   # give the user something to look at...
    if entry%1000 == 0:
      print 'Entry:',entry

    # check that the event is read properly
    entryCheck = ntuple.GetEntry( entry )
    if entryCheck <= 0:  continue


    # loop over each muon in the event and create muon pairs
    for i in range(0,ntuple.n_muons):

      for j in range(i+1,ntuple.n_muons):
        #Get the magnitude of the momentum and energy
        P_i = math.sqrt(ntuple.muon_px[i]**2 + ntuple.muon_py[i]**2 + ntuple.muon_pz[i]**2)
        E_i = math.sqrt(m_mu**2 + P_i**2)

        P_j = math.sqrt(ntuple.muon_px[j]**2 + ntuple.muon_py[j]**2 + ntuple.muon_pz[j]**2)
        E_j = math.sqrt(m_mu**2 + P_j**2)



        Ei_hist.Fill(E_i)
        Pi_hist.Fill(P_i)

    #Write the histograms to a ROOT file
    outfile = ROOT.TFile('output.root','recreate')
    Ei_hist.Write()
    Pi_hist.Write()
    outfile.Close()

  

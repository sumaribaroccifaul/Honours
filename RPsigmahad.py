import ROOT
import math
import array
import numpy as np

# DATA # 
# array of beam energies
xdata = array.array('d',[89, 91, 92, 93, 94, 95, 96, 97])
print(xdata)

# number of entries in the Z-boson peak 
entries = array.array('d',[53473, 57371, 53743, 45128, 37472, 31874, 28234, 25749])
print(entries)

# number of generated events in the sample
sample_entries = array.array('d',[100000, 100000, 100000, 100000, 100000, 100000, 100000, 100000])

# hadronic cross section values form the MC samples in [fb^{-1}]
sample_had = array.array('d', [2.114e6, 8.953e6, 6.496e6, 2.965e6, 1.525e6, 9.059e5, 5.971e5, 4.227e5])

#calculate luminosities
L = np.zeros(8)

for l in range(0, len(xdata)):
    L[l] = sample_entries[l]/sample_had[l]

print("Luminosity:",L)

efficiency = np.zeros(8)

for i in range(0, len(xdata)):
    efficiency[i] = entries[i]/sample_entries[i]

print("efficiency",efficiency)

# effective luminosities from MC samples in [fb^{-1}]

#L = array.array('d',[4.7873e-07, 2.0358e-07, 8.9986e-07, 3.1561e-07, 2.6179e-07, 1.3240e-08, 6.9108e-08, 3.7728e-08]) # real values
#L = array.array('d',[4.7873e-07, 2.0358e-07, 4.9986e-07, 5.1561e-07, 6.6179e-07, 7.3240e-07, 8.9108e-07, 9.7728e-07]) # illustrative values
#print(L)

# create empty array for sigmahad values
sigmahad = np.zeros(8)

# CALCULATION #
# calculate the cross section 
for s in range(0, len(xdata)):
    sigmahad[s] = entries[s]/L[s]

print("sigmahad:",sigmahad)

# PLOT # 
n = 8

plot = ROOT.TGraph(n,xdata,sigmahad)
plot.SetTitle("Hadronic cross section calculated from di-jet events as a function of c.m. energy ;c.m. energy [GeV]; cross section [fb]")

# FIT # 
def sigmaFitFunction(xdata,par):

    s = xdata[0] 
    sigma0 = par[0]
    M_Z = par[1]
    Gamma_Z = par[2]

    y = sigma0*( s / ( ( ((s - M_Z)/Gamma_Z)*((s - M_Z)/Gamma_Z) ) + (s/(M_Z*M_Z)) ))

    return y
  
sigmaFit = ROOT.TF1("sigmaFit",sigmaFitFunction,89,97,3)
sigmaFit.SetParameters(1, 91, 1)
#sigmaFit.SetParLimits(1, 90.8, 91.3)
sigmaFit.SetLineColor(6)

bw_A = "2*sqrt(2)*[0]*[1]*sqrt([0]*[0]*([0]*[0]+[1]*[1]))"
bw_B = "3.14159*sqrt([0]*[0] + sqrt([0]*[0]*([0]*[0]+[1]*[1])))" 
bw_C = "(x*x-[0]*[0])*(x*x-[0]*[0]) + [0]*[0]*[1]*[1]"
bw = "[2]*((%s)/(%s))/(%s)"%(bw_A, bw_B, bw_C)
bwFit = ROOT.TF1("bwfit",bw,88,92) 
bwFit.SetParameter(0,91) 
bwFit.SetParameter(1,4)
bwFit.SetParameter(2,30)

# OUTPUT # 
outfile = ROOT.TFile('output_sigmahad.root','recreate')
plot.Write()
outfile.Close()

can = ROOT.TCanvas("can", "sigmahad",500,300)
plot.Draw("*pa")
plot.Fit("sigmaFit","","",89,97)
#plot.Fit("bwFit","+","",89,97)
can.SaveAs("sigmahad.pdf")

# # GOODNESS-OF-FIT: CHI2 PER DOF #
# print("Sigma Function Fit:")
# sigma_chi2 = sigmaFit.GetChisquare()
# sigma_ndof = sigmaFit.GetNDF()
# sigma_chi2perdof = sigma_chi2/sigma_ndof
# print("sigma fit: chi2 per dof =", sigma_chi2perdof)

# print("Breit Wigner Function Fit:")
# bw_chi2 = bwFit.GetChisquare()
# bw_ndof = bwFit.GetNDF()
# bw_chi2perdof = bw_chi2/bw_ndof
# print("bw fit: chi2 per dof =", bw_chi2perdof)



# # Neutirno fit 
# # initial guesses 
# gamZ = 2.4952 # Z decay width PDG 
# gamL = (3.3658/100)*gamZ # lepton decay width PDG
# gamNu = (20/100)*gamZ # invisible (nu) decay width PDG
# gamHad = (69.911/100)*gamZ # hadronic decay width PDG
# R_l = gamHad/gamL

# # fit function 
# num = "12*pi*[0]"
# den = "( ( ( [1]*( [2]/[3] ) )+[0]+3 )**2 )*( [4]*[4] )"
# sig = "((%s)/(%s))"%(num,den)
# neutrinoFit = ROOT.TF1("neutrinoFit",sig,93,97,5) 
# neutrinoFit.SetLineColor(9)

# # make inital guesses of fitting params 
# neutrinoFit.SetParameter(0,R_l) # R_l 
# neutrinoFit.SetParameter(1,3) # N_nu
# neutrinoFit.SetParameter(2,gamNu) # Gamma_nu 
# neutrinoFit.SetParameter(3,gamL) # Gamma_l 
# neutrinoFit.SetParameter(4,91.2) # m_Z 
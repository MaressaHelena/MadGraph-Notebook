#!/usr/bin/env python
# coding: utf-8

# In[5]:


import ROOT

# onde entra os parametros?

###

#Double_t RBWGaus(Double_t *x, Double_t *par) { #preservar o que vem depois do Double_t

      #Fit parameters:
      #par[0]=Width (scale) Breit-Wigner
      #par[1]=Most Probable (MP, location) Breit mean
      #par[2]=Total area (integral -inf to inf, normalization constant)
      #par[3]=Width (sigma) of convoluted Gaussian function

def RBWGaus(x, par, xx, fbw):

    np = 100 # #number of convolution steps
    sc = 4 #convolution extends to +-sc Gaussian sigmas

          #Range of convolution integral
    xmin = x[0] - sc * par[3]
    xmax = x[0] + sc * par[3]
    step = (xmax-xmin) / np
    
      # Convolution integral of Breit-Wigner and Gaussian by sum
    sum = 0.0

###

    for i in range (1, np/2):
        xx = xmin + (i-.5) * step
        fbw = ROOT.BreitWigner(xx,par[1],par[0])
        sum += fbw * ROOT.Gaus(x[0],xx,par[3])

        xx = xmax - (i-.5) * step
        fbw = ROOT.BreitWigner(xx,par[1],par[0])
        sum += fbw * ROOT.Gaus(x[0],xx,par[3])
    return (par[2] * step * sum * (1./sqrt(2*ROOT.Pi())) / par[3])


# In[28]:


####
def Zmass_Zcharge(Zmass, Zcharge):
    c = ROOT.TCanvas()
    f = ROOT.TFile.Open("Zbosons.root")
    tr = f.Get("ztree")
    tr.SetBranchAddress("Zmass", Zmass)
    tr.SetBranchAddress("Zcharge", Zcharge)
    hist = ROOT.TH1F("hmass","", 60,60,120)
    hist.Draw()
    c.Draw()

    for j in range(0, tr.GetEntries()):

        tr.GetEntry(j)

        if(Zcharge != 0):
            continue

        hist.Fill(Zmass)
        hist.Draw()
        c.Draw()

####

c1 = ROOT.TCanvas("c1","Dimuon mass", 600, 600)

c1.SetTopMargin(0.05)
c1.SetRightMargin(0.05)
c1.SetBottomMargin(0.12)
c1.SetLeftMargin(0.13)
c1.SetTickx(1)
c1.SetTicky(1)



hist.SetMarkerStyle(20)
hist.SetMarkerColor(ROOT.kRed)
hist.SetLineColor(ROOT.kRed)
hist.GetXaxis().SetTitle("Mass (GeV)")
hist.GetYaxis().SetTitle("Counts")
hist.GetXaxis().SetTitleSize(0.05)
hist.GetYaxis().SetTitleSize(0.05)  
hist.GetYaxis().SetTitleOffset(1.2)
   #gStyle.SetOptStat(0) ---->>>> o que era pra ter aqui??
#c1.Update()

hist.Draw()
c.Draw()

###########
#f = ROOT.TF1("f",RBWGaus,60,120,4) #ROBLEMA AQUI


f.SetParameters(2.495, 91.0, 2000.0, 2.0)
f.SetParNames("BW width","BW mean","Area","Sigma")
f.FixParameter(0, 2.495); #/PDG value
f.SetParLimits(1, 86, 96)
f.SetLineColor(ROOT.kBlue)  

#FitResultPtr 
ff= hist.FitResultPtr(f,"RNS","")   #itResultPtr
hist.Draw("ep")
f.Draw("same")


# In[2]:


l = ROOT.TLegend(0.18,0.78,0.34,0.90)
l.SetTextSize(0.04)
l.AddEntry(hist,"Z#rightarrow#mu#mu","lp") # ???
l.AddEntry(f,"Fit","l")
l.Draw()

tx = ROOT.TLatex()
tx.SetTextSize(0.03)
tx.SetTextAlign(12)
tx.SetTextFont(42)
tx.SetNDC(ROOT.kTRUE)

tx.DrawLatex(0.63,0.87,ROOT.Form("#chi^{2}/ndf = %g/%d",ff.Chi2(),ff.Ndf()))
tx.DrawLatex(0.63,0.83,ROOT.Form("Yield = %g#pm%g",ff.Parameter(2),ff.ParError(2)))

c1.SaveAs("Zpeak.png")


# In[ ]:





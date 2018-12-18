# global constants and functions
# S.Chekanov (ANL)

import sys
sys.path.append("modules/")

from array import *
from math import *
from AtlasStyle import *
from AtlasUtils import *
from initialize  import *
from module_functions import *
from ROOT import TH1D,TF1,TRandom3,TFile,TLatex,TLegend,TPaveText,TGraphErrors,kRed,kGreen,kYellow,kTRUE
import math
import ROOT


PTlab="|#eta|<2.5"
dlab="Delphes ATLAS-like"
tlab="MG5 t#bar{t}+jets"
nnlab="ANN ATLAS-like"

# max bins considered
MaxJetBins=21

# make a table with systematics 
def genLatexTable(file,ISYS,NTOT,Value,Error,Chi2NDF):
  v=str(int(Value))
  e=str(int(Error))
  chi2="%.2f" % (Chi2NDF)
  s=str(ISYS)+"  &  $"+v+"  \pm  "+e+ "$  &  "+chi2+ " \\\\ \\hline\n" 

  if ISYS==0: 
     f = open (file,"w")
     f.write ('\\begin{tabular}{|c|cc|}\n')
     f.write ('\\hline \n')
     f.write (" Nr systematics & Fitted candidates  & $\chi^2/ndf$  \\\\ \\hline \n"  )
     f.write(s)
     f.close()
  else:  
     f = open (file,"a")
     f.write(s) 
     f.close()
  if (ISYS==NTOT-1):
     f = open (file,"a")
     f.write ('\\end{tabular} \n')
     f.close()


# write data for the current extraction using Fit 
def writeCurrent(ISYS,Value,Error,Chi2NDF):
  v=str(int(Value))
  e=str(int(Error))
  chi2="%.2f" % (Chi2NDF)
  s=str(ISYS)+" & "+v+" &  "+e+ "  &  "+chi2+ "\n"
  print "Write figs/current_sys.txt with the values=",s
  f = open ("figs/current_sys.txt","w")
  f.write(s)
  f.close()

# write data for the current extraction using MiniIsolation 
def writeCurrentMuonIso(ISYS,Value,Error):
  v=str(int(Value))
  e=str(int(Error))
  s=str(ISYS)+" & "+v+" &  "+e+ " \n"
  print "Write figs/current_sys_muoniso.txt with the values=",s
  f = open ("figs/current_sys_muoniso.txt","w")
  f.write(s)
  f.close()





# make a table with systematics 
def genTable(file,ISYS,NTOT,Value,Error,Chi2NDF):
  v=str(int(Value))
  e=str(int(Error))
  chi2="%.2f" % (Chi2NDF)
  s=str(ISYS)+"  &  "+v+" &  "+e+ "  &  "+chi2+ "\n"

  if ISYS==0:
     f = open (file,"w")
     f.write ('# table with systematic uncertainties \n')
     f.write(s)
     f.close()
  else:
     f = open (file,"a")
     f.write(s)
     f.close()

def drawXAxis(sf,gPad,XMIN,YMIN,XMAX,YMAX,nameX,nameY):
 h=gPad.DrawFrame(XMIN,YMIN,XMAX,YMAX);
 ay=h.GetYaxis();
 ay.SetLabelFont(42)

 if (sf==1): 
             ay.SetLabelSize(0.05)
             ay.SetTitleSize(0.05)

 if (sf==2 or sf==3): 
             ay.SetLabelSize(0.10)
 if (sf==20):
             ay.SetLabelSize(0.18)
 if (sf==30):
             ay.SetLabelSize(0.12)
             
# ay.SetTitleSize(0.1)
 ay.SetNdivisions(505);
 if (sf==1): ay.SetTitle( nameY )
 # ay.Draw("same")
 ax=h.GetXaxis(); 
 if (sf==1 or sf==2): ax.SetTitle( nameX );
 if (sf==30): ax.SetTitle( nameX );
 ax.SetTitleOffset(1.2)
 ay.SetTitleOffset(1.2)

 ax.SetLabelFont(42)
 # ax.SetTitleFont(62)
 ay.SetTitleFont(62)
 ay.SetLabelFont(42)
 # ay.SetTitleFont(42)
 ax.SetLabelSize(0.10)
 ax.SetTitleSize(0.12)
 if (sf==30): 
          ax.SetLabelSize(0.12)
          ax.SetTitleSize(0.14)
 ax.Draw("same");


def getStyle(hh,marker_type=20, color=1):
  hh.SetTitle("")
  hh.SetLineStyle(1)
  hh.SetLineWidth(2)
  hh.SetMarkerStyle(marker_type)
  hh.SetMarkerColor(color)
  hh.SetMarkerSize(0.6)
  hh.SetStats(0)
  return hh
 
# get resolution
def getResolutionGraph(file):
 print "Read=",file
 i=0
 g1=TGraphErrors()
 for line in open(file,"r"):
   columns = line.split(',')
   pt=float(columns[0])
   pt_err=float( columns[1])
   sigma=float( columns[4])
   err=float( columns[5])
   g1.SetPoint(i,pt,sigma)
   ex=pt_err;
   ey=err;
   g1.SetPointError(i,ex,ey)
   i=i+1;

 g1.SetMarkerColor( 1 )
 g1.SetMarkerStyle( 20 )
 g1.SetMarkerSize( 0.5 )

 # g1->Print();
 return g1
 
# get response 
def getResponseGraph(file):
 print "Read=",file
 i=0
 g1=TGraphErrors()
 for line in open(file,"r"):
   columns = line.split(',')
   pt=float(columns[0])
   pt_err=float( columns[1])
   sigma=float( columns[2])
   err=float( columns[3])
   g1.SetPoint(i,pt,sigma)
   ex=pt_err;
   ey=err;
   g1.SetPointError(i,ex,ey)
   i=i+1;

 g1.SetMarkerColor( 1 )
 g1.SetMarkerStyle( 20 )
 g1.SetMarkerSize( 0.5 )

 # g1->Print();
 return g1


 
# input:
# h1- pT(rec)/pT(true)
# h2- pT(true) - for mean
# isFit=False - use RMS instead of fits
def plot(c1,h1,h2, peak=1.0, sigma=0.1, MyMin=0, MyMax=2, mu=140,Xmin=0.0,Xmax=2.0,Ymin=0,Ymax=500, filename="root/resolution.txt",isFit=True):
  nameX="p_{T}^{reco, jet} / p_{T}^{true, jet}"
  nameY="Entries"
  mean=h2.GetMean()
  # if (mean<1): return

  # mean_err=h2.GetMeanError()
  mean_err=h2.GetRMS()
  xmean=h1.GetMean()
  xrms=h1.GetRMS()
  txt='#LT p_{T} #GT =%.0f GeV'%( mean )
  print txt
  h1.Sumw2()
  h1=getStyle(h1,20,1)
  h1.SetAxisRange(Ymin, Ymax,"y");
  h1.SetAxisRange(Xmin, Xmax,"x");
  h1.Draw("pe same")
  xmean=h1.GetMean()
  xrms=h1.GetRMS()
  MyMin=xmean-4*xrms;
  MyMax=xmean+4*xrms;
  signal=TF1("signal",Gauss(),MyMin,MyMax,3);
  signal.SetNpx(200); signal.SetLineColor(2); signal.SetLineStyle(1)
  signal.SetLineWidth(2)
  binmax = h1.GetMaximumBin();
  yyy = 0.9*h1.GetBinContent(binmax);
  signal.SetParameter(0,yyy) # amplitude  
  signal.SetParameter(1,xmean)
  signal.SetParameter(2,xrms)
  signal.SetParLimits(2,0.00001,1.0)
  h1.Fit(signal,"MRE+","",MyMin,MyMax)
  par = signal.GetParameters()
  # signal.Draw("same")
  chi2= signal.GetChisquare()
  ndf=signal.GetNDF()
  chi2ndf=0
  if (ndf>0): 
           chi2ndf=chi2/ndf
           print "Chi2=", chi2," ndf=",ndf, " chi2/ndf=",chi2ndf
  par = signal.GetParameters()
  err=signal.GetParErrors()
  c='Peak=%.3f'%( par[1] ) + " #pm " + '%.3f GeV'%(err[1])
  b='#sigma=%.3f'%( par[2] ) + " #pm " + '%.3f GeV'%(err[2])
  d='#chi^{2}/ndf=%.2f'%( chi2ndf )
  ax=h1.GetXaxis(); ax.SetTitleOffset(1.0)
  ax.SetTitle( nameX );
  ay=h1.GetYaxis();
  ay.SetTitle( nameY );
  ay.SetTitleSize(0.072);
  ax.SetTitleSize(0.072);
  ay.SetLabelSize(0.06)
  ax.SetTitleOffset(1.1);
  ay.SetTitleOffset(0.6)
  ay.SetLabelFont(42)
  # ay.SetTitleFont(42)
  ax.SetLabelFont(42)
  ax.SetTitleFont(62)
  ay.SetTitleFont(62)

  ax.SetLabelSize(0.06)
  myText( 0.16,0.90,1,0.065,txt)
  myText( 0.16,0.80,1,0.065,c)
  myText( 0.16,0.70,1,0.065,b)
  myText( 0.16,0.60,1,0.065,d)
  ss=par[2]*TMath.Sqrt(mean)*100
  # S='%.0f / #sqrt{ E_{T} }'%( ss )
  # myText( 0.16,0.5,4,0.08,S)
  # myText( 0.78,0.8,2,0.08,"#mu="+str(mu))

  # if fit fail, switch to RMS
  if (par[2] !=0): 
           print "Relative error on width=",err[2]/par[2] 
           if (err[2]/par[2]>0.1):
                        isFit=False;
                        print "Fit is unreliable! use RMS 90%"
  if (par[2] ==0): isFit=False;

  # Use RMS of histogram, no fit
  if (isFit == False):
           par[1]=h1.GetMean()
           err[1]=h1.GetMeanError()
           par[2]=rms90(h1); 
           # par[2]=h1.GetRMS()
           err[2]=h1.GetRMSError()
           print "Replace fits with histogram RMS!"

  print "Write="+filename
  file=open(filename,"a")
  file.write(str(mean)+","+str(mean_err)+","+str(par[1])+","+str(err[1])+","+str(par[2])+","+str(err[2])+"\n")
  file.close()

# h1- pT(rec)/pT(true)
# h2- pT(true) - for mean
# isFit=False - use RMS instead of fits
def ploteta(c1,h1,h2, peak=1.0, sigma=0.1, MyMin=0, MyMax=2, mu=140,Xmin=0.0,Xmax=2.0,Ymin=0,Ymax=500, filename="root/resolution_eta.txt",isFit=True):
  nameX="#eta_^{reco, jet} / #eta^{true, jet}"
  nameY="Entries"
  mean=h2.GetMean()
  # if (mean<1): return
  # mean_err=h2.GetMeanError()
  mean_err=h2.GetRMS()
  xmean=h1.GetMean()
  xrms=h1.GetRMS()
  txt='#LT p_{T} #GT =%.0f GeV'%( mean )
  print txt
 
  par=[]
  par.append(0); 
  par.append(0);
  par.append(0);
  err=[]
  err.append(0);
  err.append(0);
  err.append(0);

  # Use RMS of histogram, no fit
  par[1]=h1.GetMean()
  err[1]=h1.GetMeanError()
  par[2]=rms90(h1); 
  # par[2]=h1.GetRMS()
  err[2]=h1.GetRMSError()
  print "Replace fits with histogram RMS!"

  print "Write="+filename
  file=open(filename,"a")
  file.write(str(mean)+","+str(mean_err)+","+str(par[1])+","+str(err[1])+","+str(par[2])+","+str(err[2])+"\n")
  file.close()



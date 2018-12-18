import sys
sys.path.append("modules/")

	
from global_module import *

# is it batch?
# if batch, set input to batch
myinput="no any input"
if (len(sys.argv) > 1):
   myinput = sys.argv[1]
print "Mode=",myinput


# import atlas styles
from AtlasStyle import *
from AtlasUtils import *
from initialize  import *
from module_functions import *
import ROOT


gROOT.Reset()
figdir="figs/"
epsfig=figdir+__file__.replace("py","eps")
Ymin=0.0 
Ymax=4000.0 
Xmin=-2  
Xmax=2

myinput="interactive"
xdir=""

tag="rfast008"
if (len(sys.argv) ==2):
   myinput = sys.argv[1]
print "mode =",myinput

if (len(sys.argv) ==3):
   tag = sys.argv[1]
   myinput = sys.argv[2]
print "TAG =",tag
print "mode =",myinput


filename="root/resolution_jet1_"+tag+".txt"
import os
if os.path.isfile(filename):
     os.remove(filename)


######################################################
gROOT.SetStyle("Plain");


c1=TCanvas("c","BPRE",10,10,700,700);
ps1 = TPostScript( epsfig,113)
c1.Divide(4,5,0.005,0.005);
c1.SetTickx()
c1.SetTicky()
c1.SetTitle("")
c1.SetLineWidth(3)
c1.SetBottomMargin(0.12)
c1.SetTopMargin(0.05)
c1.SetLeftMargin(0.13)
c1.SetFillColor(0)

sig=TFile("out/"+tag+"_histo.root")
sig.ls()

for i in range(MaxJetBins):
  xx="%02d" % (i,)
  name="jet1_resolution_"+xx
  name1="jet1_resolution_pt_"+xx
  c1.cd(i)
  gPad.SetBottomMargin(0.18)
  h1=sig.Get(name)
  mean=sig.Get(name1)
  binmax = h1.GetMaximumBin();
  yyy = h1.GetBinContent(binmax);
  Ymax=1.2*yyy
  mu=int(mean.GetMean())
  print "Processing pT =", mu, " GeV"
  plot(c1,h1,mean,0.8,0.1,0.0,4.0,mu,0.0,4.0,Ymin,Ymax,filename,False)
  gPad.RedrawAxis()
 


ps1.Close()
if (myinput != "-b"):
              if (raw_input("Press any key to exit") != "-9999"): 
                         c1.Close(); sys.exit(1);



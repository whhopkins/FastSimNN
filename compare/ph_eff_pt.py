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

Ymin=0.0 
Ymax=2.0 
Xmin=0  
Xmax=3000

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

gROOT.Reset()
figdir="figs/"
epsfig=figdir+__file__.replace("py","eps")
epsfig=epsfig.replace(".eps","_"+tag+".eps")

######################################################
gROOT.SetStyle("Plain");


c1=TCanvas("c","BPRE",10,10,500,400);
ps1 = TPostScript( epsfig,113)
c1.Divide(1,3,0.005,0.005);
c1.SetTickx()
c1.SetTicky()
c1.SetTitle("")
c1.SetLineWidth(3)
c1.SetBottomMargin(0.13)
c1.SetTopMargin(0.05)
c1.SetLeftMargin(0.14)
c1.SetRightMargin(0.05)
c1.SetFillColor(0)


h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax)
gPad.SetLogy(0)
gPad.SetLogx(0)
ax=h.GetXaxis(); ax.SetTitleOffset(1.0)
ax.SetTitle( "p_{T}^{#gamma} [GeV]"  );
ay=h.GetYaxis();
ay.SetTitle( "#gamma efficiency" );
ay.SetTitleSize(0.05);
ax.SetTitleSize(0.05);
ay.SetLabelSize(0.04)
ax.SetTitleOffset(1.1);
ay.SetTitleOffset(1.25)
ay.SetLabelFont(42)
ax.SetLabelFont(42)
ax.SetLabelSize(0.04)

ax.SetTitleFont(62)
ay.SetTitleFont(62)


sig=TFile("out/"+tag+"_histo.root")
sig.ls()

h1=sig.Get("ph_reco_pt")
h2=sig.Get("ph_true_pt")
h1.Divide(h2)
 
h11=sig.Get("ph_nn_pt")
h22=sig.Get("ph_true_pt")
h11.Divide(h22)

p1=TH1toTGraph(h1)
p1.SetMarkerSize(1.0)
p1.SetMarkerStyle(20)


p2=TH1toTGraph(h11)
p2.SetMarkerColor(2)
p2.SetMarkerStyle(24)
p2.SetMarkerSize(1.0)
p2.Draw("pe same")
p1.Draw("pe same")

leg2=TLegend(0.6, 0.2, 0.85, 0.4);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.AddEntry(p1,dlab, "pl")
leg2.AddEntry(p2,nnlab, "pl")
leg2.Draw("same");

myText( 0.2,0.84,4,0.04,"Direct #gamma")

txtt="<#mu>=0"
if (tag == "rfast007"): txtt="<#mu>=40" 
if (tag == "rfast006"): txtt="<#mu>=140"
myText( 0.2,0.9,2,0.05,txtt)


ps1.Close()
if (myinput != "-b"):
              if (raw_input("Press any key to exit") != "-9999"): 
                         c1.Close(); sys.exit(1);



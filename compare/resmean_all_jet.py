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

tag="rfast004"
if (len(sys.argv) ==2):
   tag = sys.argv[1]
print "TAG Mode=",myinput

gROOT.Reset()
figdir="figs/"
epsfig=figdir+__file__.replace("py","eps")
epsfig=epsfig.replace(".eps","_"+tag+".eps")
Ymin=0.5 
Ymax=1.5  
Xmin=20   
Xmax=1000 

if (tag == "rfast007"):
        Ymax=2.5;
if (tag == "rfast006"):
        Ymax=3.5;


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
c1.SetLeftMargin(0.13)
c1.SetRightMargin(0.05)
c1.SetFillColor(0)

h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax)
gPad.SetLogy(0)
gPad.SetLogx(1)
ax=h.GetXaxis(); ax.SetTitleOffset(1.0)
ax.SetTitle( "p_{T}^{jet} [GeV]"  );
ay=h.GetYaxis();
ay.SetTitle( "Jet response" );
ay.SetTitleSize(0.05);
ax.SetTitleSize(0.05);
ay.SetLabelSize(0.04)
ax.SetTitleOffset(1.1);
ay.SetTitleOffset(1.1)
ay.SetLabelFont(42)
ax.SetLabelFont(42)
ax.SetTitleFont(62)
ay.SetTitleFont(62)
ax.SetLabelSize(0.04)

ptmin=20
ptmax=600
s1="sqrt(([0]/sqrt(x))*([0]/sqrt(x))+([1]/x)*([1]/x)+[2]*[2])"

"""
https://twiki.cern.ch/twiki/bin/view/Main/JetResolutionPlots
a) sampling and statistical fluctuations (a^2/E), which is dominating over most of the useful range of calorimeters
b) noise, which dominates the low-energy performance of calorimeters (b^2/E^2)
c) effects coming from calibration errors, non-uniformities and non-linearities 
"""

f1=TF1("f1",s1,ptmin,ptmax);
f1.SetLineStyle(1)
f1.SetLineWidth(2)
f1.SetParameter(0,0.5)
f1.SetParameter(1,4.58662e-01)
f1.SetParameter(2,2.02282e-02)

h1=getResponseGraph( "root/resolution_jet1_"+tag+".txt" )
h2=getResponseGraph( "root/resolution_jet2_"+tag+".txt" )


#h1.SetMarkerColor(ROOT.kBlack)
h1.SetMarkerSize(1.0)
h1.SetMarkerStyle(20)
h1.Draw("P same")

h2.SetMarkerColor(2)
h2.SetMarkerStyle(24)
h2.SetMarkerSize(1.0)
h2.Draw("pe same")

# h1.Fit(f1,"RE+","",ptmin,ptmax);
par1 = f1.GetParameters()
a1='%.2f'%( par1[0] )
b1='%.2f'%( par1[1] )
k1='%.2f'%( par1[2] )
s1="#sqrt{"+a1+" ^{2} /p_{T} +"+b1+" ^{2}/p_{T}^{2} +"+k1+"^{2}}" 


# myText( 0.2,0.85,4,0.05,dlab)
myText( 0.2,0.8,4,0.04,"antiKT R=0.4 jets")

leg2=TLegend(0.6, 0.7, 0.85, 0.9);
leg2.SetBorderSize(0);
leg2.SetTextFont(62);
leg2.SetFillColor(10);
leg2.SetTextSize(0.04);
leg2.AddEntry(h1,dlab, "pl")
leg2.AddEntry(h2,nnlab, "pl")
leg2.Draw("same");


if (tag == "rfast004"):
      myText( 0.2,0.9,2,0.05,"<#mu>=0")
if (tag == "rfast007"):
      myText( 0.2,0.9,2,0.05,"<#mu>=40")
if (tag == "rfast006"):
      myText( 0.2,0.9,2,0.05,"<#mu>=140")


#myText( 0.2,0.9,4,0.04,lab1)
#myText( 0.2,0.85,4,0.04,lab2)
# myText( 0.2,0.85,4,0.03,lab3)

line =TLine(Xmin,1,Xmax,1);
line.Draw()



gPad.RedrawAxis()



ps1.Close()
if (myinput != "-b"):
              if (raw_input("Press any key to exit") != "-9999"): 
                         c1.Close(); sys.exit(1);



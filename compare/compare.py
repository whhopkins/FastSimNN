# S.Chekanov
from ROOT import gROOT,gPad,gStyle,TCanvas,TSpline3,TFile,TLine,TLatex,TAxis,TLegend,TPostScript
from ROOT import TH2D,TF2,TArrow,TCut,TPad,TH1D,TF1,TObject,TPaveText,TGraph,TGraphErrors,TGraphAsymmErrors
from ROOT import TGraph2D,TTree,TMultiGraph,TBranch,gSystem,gDirectory
from ROOT import TPaveStats,TProfile
import ROOT
import os,sys,math
#sys.path.append("nnet/")
#import bpnn

print ('Number of arguments:', len(sys.argv), 'arguments.') 
print ('Argument List:', str(sys.argv))
print ('Use as: script.py -b 0 (or 1,2)') 
myinput="interactive"

myinput="../validate/out_root/rfast004.root" 

if (len(sys.argv) ==2):
   myinput = sys.argv[1]
print ("Mode=",myinput) 

outfile=os.path.basename(myinput)


gROOT.Reset()
figdir="figs/"
epsfig=figdir+__file__.replace(".py",".eps")



######################################################
gROOT.SetStyle("Plain");
gStyle.SetLabelSize(0.035,"xyz");


rfile=ROOT.TFile.Open(myinput)

etaMax=2.5
ptMin=25

ROOT.gROOT.ProcessLine('.L Loader.C+')

nbins=10
ymax=2.5
ptmax=400

prof_eta=TProfile("prof_eta", "p-ptrue/ptrue Delphes", nbins, -1*etaMax, etaMax, 0, 10) 
prof_pt=TProfile("prof_pt", "p-ptrue/ptrue Delphes", nbins, ptMin, ptmax, 0, 10)
prof_eta_nn=TProfile("prof_eta_nn", "p-ptrue/ptrue NN", nbins, -1*etaMax, etaMax, 0, 10)
prof_pt_nn=TProfile("prof_pt_nn", "p-ptrue/ptrue NN", nbins, ptMin, ptmax, 0, 10)

nn=0
for e in rfile.Ntuple:
     if ((nn%1000==0) or nn<100): print "Event=",nn," Jet size=",e.AntiKt4TruthJetPt.size() 
     nn +=1
     for j1 in range(e.AntiKt4TruthJetPt.size()):
        ptT=e.AntiKt4TruthJetPt[j1]
        etaT=e.AntiKt4TruthJetEta[j1]
        phiT=e.AntiKt4TruthJetPhi[j1]
        if (abs(etaT)>etaMax): continue;
        if (ptT<ptMin): continue;
        for j2 in range(e.AntiKt4JetPt.size()):
           pt=e.AntiKt4JetPt[j2]
           eta=e.AntiKt4JetEta[j2]
           phi=e.AntiKt4JetPhi[j2]
           deta=etaT-eta
           dphi=phiT-phi
           if (abs(dphi)>math.pi): dphi=2*math.pi-abs(dphi) 
           dr=math.sqrt(deta*deta + dphi*dphi)
           if (dr<0.2):
               ratio=(pt/ptT) -1.0 
               prof_eta.Fill(etaT,ratio);
               prof_pt.Fill(ptT,ratio);
        for j2 in range(e.AntiKt4NNJetPt.size()):
           pt=e.AntiKt4NNJetPt[j2]
           eta=e.AntiKt4NNJetEta[j2]
           phi=e.AntiKt4NNJetPhi[j2]
           deta=etaT-eta
           dphi=phiT-phi
           if (abs(dphi)>math.pi): dphi=2*math.pi-abs(dphi)
           dr=math.sqrt(deta*deta + dphi*dphi)
           if (dr<0.2):
               ratio=(pt/ptT) -1.0
               prof_eta_nn.Fill(etaT,ratio);
               prof_pt_nn.Fill(ptT,ratio);
 
c1=TCanvas("c","FastNN",10,10,600,600);

Xmin=-2.5
Xmax=2.5
Ymin=0.0
Ymax=0.4
nameX="eta"
nameY="pT-pT(true)/pT(true)"
h=gPad.DrawFrame(Xmin,Ymin,Xmax,Ymax)
ax=h.GetXaxis(); ax.SetTitleOffset(0.8)
ax.SetTitle( nameX );
ay=h.GetYaxis(); ay.SetTitleOffset(0.8)
ay.SetTitle( nameY );
ax.SetTitleOffset(1.1); ay.SetTitleOffset(1.4)
ax.Draw("same")
ay.Draw("same")

prof_eta.SetMarkerSize(1.1)
prof_eta.SetMarkerColor(1)
prof_eta.Draw("pe same")
prof_eta_nn.SetMarkerSize(1)
prof_eta_nn.SetMarkerColor(2)
prof_eta_nn.Draw("pe same")

out="histo_"+outfile
hfile=TFile(out,"RECREATE","DijetsMC")
c1.Write()
prof_eta.Write()
prof_pt.Write()

prof_eta_nn.Write()
prof_pt_nn.Write()
hfile.ls()
hfile.Close()
print  "Write with cross sections written=",out
     
c1.Close()
rfile.Close()




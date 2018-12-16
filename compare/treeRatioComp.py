#!/usr/bin/env python

from ROOT import TH1F, TFile, TCanvas, gROOT, TPad, gPad, TLine, TChain
import argparse, os, math
from setStyle import setStyle, CompLegend

compDict = {
    '(AntiKt4JetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':['(AntiKt4NNJetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]'],
    '(AntiKt4JetPt[0])/AntiKt4TruthJetPt[0]':['(AntiKt4NNJetPt[0])/AntiKt4TruthJetPt[0]'],
}
compDictNames = {
    '(AntiKt4JetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':'Reco',
    '(AntiKt4NNJetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':'NN',
    '(AntiKt4JetPt[0])/AntiKt4TruthJetPt[0]':'Reco',
    '(AntiKt4NNJetPt[0])/AntiKt4TruthJetPt[0]':'NN',

    }
varNames = [
    '(AntiKt4JetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]',
    '(AntiKt4NNJetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]',    
    '(AntiKt4JetPt[0])/AntiKt4TruthJetPt[0]',
    '(AntiKt4NNJetPt[0])/AntiKt4TruthJetPt[0]',    
    ]
varNamesUnits = {
    '(AntiKt4JetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':('(p_{T,reco}^{0}-p_{T,truth}^{0})/p_{T,truth}^{0}', ''),
    '(AntiKt4NNJetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':('(p_{T,NN}^{0}-p_{T,truth}^{0})/p_{T,truth}^{0}', ''),    
    '(AntiKt4JetPt[0])/AntiKt4TruthJetPt[0]':('p_{T,reco}^{0}/p_{T,truth}^{0}', ''),
    '(AntiKt4NNJetPt[0])/AntiKt4TruthJetPt[0]':('p_{T,reco}^{0}/p_{T,truth}^{0}', ''),
    }
binning = {
    '(AntiKt4JetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':[50, -1, 1],
    '(AntiKt4NNJetPt[0]-AntiKt4TruthJetPt[0])/AntiKt4TruthJetPt[0]':[50, -1, 1],    
    '(AntiKt4JetPt[0])/AntiKt4TruthJetPt[0]':[50, 0, 2],
    '(AntiKt4NNJetPt[0])/AntiKt4TruthJetPt[0]':[50, 0, 2],
    }

def resComp(inFPath):
    """! 
    Plot two distributions on the same plot and also plot the ratio of the two. 
    The input file must have trees.
    """
    gROOT.SetBatch(True);
    gROOT.ForceStyle();
    setStyle()
    inF = TFile.Open(inFPath)
    tree = inF.Get('Ntuple')
    print tree.GetEntries()
    c = TCanvas('c', '', 800, 600)
    c.cd();
    pad = TPad("plotPad","",0,0,1,1);
    pad.SetPad(0.,0.37,1.,1.)
    pad.SetBottomMargin(0.02)
    pad.Draw();
    c.cd()
    ratioPad = TPad("ratioPad", "",0.,0.,1.,1.);
    ratioPad.SetPad(0.,0.05,1.,0.35)
    ratioPad.SetTopMargin(0.04)
    ratioPad.SetBottomMargin(0.3)
    ratioPad.Draw()
    refHists   = {}
    compHists  = {}
    ratioHists = {}
    
    for varName in compDict:
        pad.Clear();
        pad.SetLogy(False)
        ratioPad.Clear();
        refHistName = 'hist_'+varName.replace('[','_').replace(']','_').replace('-','_').replace('/','_').replace('(','_').replace(')','_')
        refHists[varName] = TH1F(refHistName, '', binning[varName][0], binning[varName][1], binning[varName][2])
        tree.Draw(varName+'>>'+refHistName)
        compHists[varName] = []
        ratioHists[varName] = []
        print varName, refHistName
        print refHists[varName].Integral()
        refHists[varName].Scale(1./refHists[varName].Integral())
        maxY = refHists[varName].GetMaximum()
        for compVarI in range(len(compDict[varName])):
            compHistName = 'hist_'+compDict[varName][compVarI].replace('[','_').replace(']','_').replace('-','_').replace('/','_').replace('(','_').replace(')','_')
            compHist = TH1F(compHistName, '', binning[varName][0], binning[varName][1], binning[varName][2])
            tree.Draw(compDict[varName][compVarI]+'>>'+compHistName);
            compHists[varName].append(compHist)
            compHists[varName][compVarI].Scale(1./compHists[varName][compVarI].Integral())

            if maxY < compHists[varName][compVarI].GetMaximum():
                maxY = compHists[varName][compVarI].GetMaximum()

        pad.Clear();
        refHists[varName].SetTitle('')
        refHists[varName].SetMaximum(1.2*maxY)
        refHists[varName].Draw();
        unitsStr = ''
        if varNamesUnits[varName][1] != '':
            unitsStr = ' ['+varNamesUnits[varName][1]+']'
        refHists[varName].GetXaxis().SetTitle(varNamesUnits[varName][0]+unitsStr)
        binSize = refHists[varName].GetBinLowEdge(2)-refHists[varName].GetBinLowEdge(1)
        roundTo = -int(math.floor(math.log10(abs(binSize))));
        roundedBinSize = round(binSize, roundTo)
        if roundTo < 0:
            binSizeStr = "%d" % roundedBinSize
        else:
            binSizeStr = "{0:.{1}f}".format(roundedBinSize, roundTo)
        yLabel = "1/N dN/(%s) [%s %s]^{-1}" % (varNamesUnits[varName][0], binSizeStr, varNamesUnits[varName][1])

        refHists[varName].GetYaxis().SetTitle(yLabel)
            
        ratioPad.Clear()
        for compVarI in range(len(compDict[varName])):
            pad.cd()
            compHists[varName][compVarI].SetLineColor(compVarI+2);
            compHists[varName][compVarI].SetMarkerColor(compVarI+2);
            compHists[varName][compVarI].Draw('same')
            ratioPad.cd()

            ratioHists[varName].append(refHists[varName].Clone())
            ratioHists[varName][compVarI].GetYaxis().SetTitleOffset(0.4)
            ratioHists[varName][compVarI].GetYaxis().SetDecimals(True);
            ratioHists[varName][compVarI].GetYaxis().CenterTitle(True);
            ratioHists[varName][compVarI].Divide(compHists[varName][compVarI])
            ratioHists[varName][compVarI].SetMaximum(2)
            ratioHists[varName][compVarI].SetMinimum(0.)
            ratioHists[varName][compVarI].GetYaxis().SetTitleSize(compHists[varName][compVarI].GetXaxis().GetTitleSize()*2.1)
            ratioHists[varName][compVarI].GetXaxis().SetTitleSize(compHists[varName][compVarI].GetYaxis().GetTitleSize()*2.1)

            ratioHists[varName][compVarI].GetYaxis().SetLabelSize(compHists[varName][compVarI].GetXaxis().GetLabelSize()*2.1)
            ratioHists[varName][compVarI].GetXaxis().SetLabelSize(compHists[varName][compVarI].GetYaxis().GetLabelSize()*2.1)
            ratioHists[varName][compVarI].GetYaxis().SetTitle('Reco/NN');
            
            ratioHists[varName][compVarI].SetLineColor(compVarI+2);
            ratioHists[varName][compVarI].SetMarkerColor(compVarI+2);
            if compVarI == 0:
                ratioHists[varName][compVarI].Draw()
            else:
                ratioHists[varName][compVarI].Draw('same')
            line = TLine()
            line.DrawLine(refHists[varName].GetXaxis().GetXmin(),1,refHists[varName].GetXaxis().GetXmax(),1)
            line.SetLineWidth(1);
            line.SetLineStyle(2);
            line.SetLineColor(1);
            c.cd();
            leg = CompLegend((0.2, 0.8, 0.4, 0.9), [refHists[varName]]+compHists[varName], [compDictNames[varName]]+[compDictNames[compVarName] for compVarName in compDict[varName]])
            leg.Draw('same')
            c.Print(varName.replace('/','_').replace('[','_').replace(']','_').replace('(', '').replace(')', '').rstrip('_').lstrip('_')+'.pdf')
            pad.SetLogy();
            refHists[varName].SetMaximum(10*maxY)
            c.Print(varName.replace('/','_').replace('[','_').replace(']','_').replace('(', '').replace(')', '').rstrip('_').lstrip('_')+'_log.pdf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare jet kinematics between different trees ')
    parser.add_argument('--inFPath', help='Location of ROOT files.', default='/users/chekanov/public/FastNNdelphes/rfast008.root')

    args = parser.parse_args() 
    resComp(args.inFPath)



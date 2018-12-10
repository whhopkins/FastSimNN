from ROOT import gStyle, TLegend, TColor
from array import array

class CompLegend(TLegend):
    def __init__(self, (x1,x2,y1,y2), hists, titles):
        TLegend.__init__(self, x1, x2, y1, y2)
        for i in range(len(hists)):
            #self.AddEntry(hists[i], titles[i], "f")
            self.AddEntry(hists[i], titles[i], "p")

        self.SetTextAlign(12);
        self.SetFillColor(0);
        self.SetEntrySeparation(0.5);
        self.SetTextSize(0.03);
        #self.SetTextSize(0.03);
        self.SetShadowColor(10);
        self.SetLineColor(10);
        self.SetFillColor(10);
        self.SetFillStyle(0); 
        self.SetBorderSize(0);

    def SetCoords(self, x1, x2, y1, y2):
        self.SetX1NDC(x1); self.SetX2NDC(x2); self.SetY1NDC(y1); self.SetY2NDC(y2); 
        
def setStyle():
    TColor.InitializeColors();
    stops =array('d', [0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000]);
    red   =array('d', [0.9764, 0.9956, 0.8186, 0.5301, 0.1802, 0.0232, 0.0780, 0.0592, 0.2082]);
    green =array('d', [0.9832, 0.7862, 0.7328, 0.7492, 0.7178, 0.6419, 0.5041, 0.3599, 0.1664]);
    blue  =array('d', [0.0539, 0.1968, 0.3499, 0.4662, 0.6425, 0.7914, 0.8385, 0.8684, 0.5293]);
    TColor.CreateGradientColorTable(9, stops, red, green, blue, 255, 1.0);
    gStyle.SetHatchesLineWidth(1) 
    gStyle.SetNumberContours(255)
    gStyle.SetNdivisions(505,"xyz");
    gStyle.SetTitleBorderSize(0)
    gStyle.SetTitleColor(1)
    gStyle.SetTitleStyle(3013)
    gStyle.SetTitleFillColor(0)
    gStyle.SetTitleX(0.05)
    gStyle.SetTitleW(0.4)
    gStyle.SetTitleY(0.965)
    gStyle.SetTitleH(0.065)
    gStyle.SetCanvasBorderMode(0)
    gStyle.SetCanvasBorderSize(3)
    gStyle.SetCanvasColor(0)
    gStyle.SetPadColor(0)
    gStyle.SetPadBorderMode(0)
    gStyle.SetFuncWidth(3)
    gStyle.SetPadGridY(False)
    gStyle.SetPadGridX(False)
    gStyle.SetFrameLineWidth(1)
    gStyle.SetMarkerSize(5)
    gStyle.SetPadTickX(True)
    gStyle.SetPadTickY(True)
    #gStyle.SetPalette(1)
    gStyle.SetHistLineWidth(3)
    gStyle.SetHistLineColor(1)
    gStyle.SetOptStat(0)
    gStyle.SetOptFit(0)
    gStyle.SetStatW(0.25)
    gStyle.SetStatH(0.25)
    gStyle.SetStatX(0.9)
    gStyle.SetStatY(0.9)
    gStyle.SetStatColor(0)
    gStyle.SetStatFormat("6.4g");
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetPadRightMargin(0.05)
    gStyle.SetPadBottomMargin(0.16)
    gStyle.SetPadLeftMargin(0.12)
    font = 42;
    gStyle.SetTextSize(.055)

    gStyle.SetTextFont(font);
    gStyle.SetLabelFont(font,"x");
    gStyle.SetTitleFont(font,"x");
    gStyle.SetLabelFont(font,"y");
    gStyle.SetTitleFont(font,"y");
    gStyle.SetLabelFont(font,"z");
    gStyle.SetTitleFont(font,"z");
    
    gStyle.SetTitleSize(.055, "x")
    gStyle.SetTitleSize(.055, "y")
    gStyle.SetTitleSize(.05, "z")
    gStyle.SetLabelSize(.05, "x")
    gStyle.SetLabelSize(.05, "y")
    gStyle.SetLabelSize(.05, "z")
    gStyle.SetLabelOffset(0.014,"x");
    gStyle.SetLabelOffset(0.006,"y");
    gStyle.SetLabelOffset(0.008,"z");
    gStyle.SetTitleOffset(1, "y");
    gStyle.SetTitleXOffset(1.2);
   
    # use bold lines and markers
    gStyle.SetMarkerStyle(20);
    gStyle.SetMarkerColor(1);
    gStyle.SetMarkerSize(1.2);
    gStyle.SetLineStyleString(2,"[12 12]"); # postscript dashes

    gStyle.SetPadTickX(1);
    gStyle.SetPadTickY(1);



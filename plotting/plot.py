#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time

import ROOT

ROOT.gROOT.ProcessLine(".L ~/tdrstyle.C");
ROOT.setTDRStyle();
ROOT.gStyle.SetPadLeftMargin(0.16);


############################################################

if __name__ == '__main__':

    file = ROOT.TFile("out.root");
    tree = file.Get("tree");
    
    vars = ["mWW","mWLep","mWHad","costheta1","costheta2","costhetastar","phi","phi1","dEtajj","dPhijj","mjj"]
    vars_x = ["m_{WW}","m_{WLep}","m_{WHad}","costheta1","costheta2","costhetastar","phi","phi1","dEtajj","dPhijj","mjj"]
    vars_nbin = [20]*11;
    vars_lo = [   0,   0,   0,  -1,   0,   -1, -ROOT.TMath.Pi(), -ROOT.TMath.Pi(), 0,                0,    0]
    vars_hi = [1000, 200, 200,   1,   1,    1,  ROOT.TMath.Pi(),  ROOT.TMath.Pi(),10,  ROOT.TMath.Pi(), 1000]    
    
    h_sig = [];
    h_bkg = [];
    
    for i in range(len(vars)):
        h_sig.append( ROOT.TH1F("hs_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );
        h_bkg.append( ROOT.TH1F("hb_"+vars[i],";"+vars_x[i]+";n.d.",vars_nbin[i],vars_lo[i],vars_hi[i]) );   

    
    for i in range( tree.GetEntries() ):
        tree.GetEntry(i);
        if i%2000 == 0: print "i = ",i
        
        if tree.isSignal == 1:
            for j in range(len(vars)):
                h_sig[j].Fill( getattr( tree, vars[j] ) );
        else:
            for j in range(len(vars)):
                h_bkg[j].Fill( getattr( tree, vars[j] ) );
                
    ##### plot
    # scale and color
    for j in range(len(vars)):
        h_sig[j].Scale( 1./h_sig[j].Integral() );
        h_bkg[j].Scale( 1./h_bkg[j].Integral() );
        h_bkg[j].SetLineColor( 2 )
        
    # on canvas
    for j in range(len(vars)):
        
        tmpcan = ROOT.TCanvas("can"+str(j),"can"+str(j),800,800);
        h_sig[j].SetMaximum( 1.2*max( h_sig[j].GetMaximum(), h_bkg[j].GetMaximum() ) );
        h_sig[j].SetMinimum( 0 );        
        h_sig[j].Draw();
        h_bkg[j].Draw("sames");
        tmpcan.SaveAs("plots/"+vars[j]+".eps");
        tmpcan.SaveAs("plots/"+vars[j]+".png");        
        
        
        
        
        
        
        
        
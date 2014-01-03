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
ROOT.gStyle.SetPadRightMargin(0.10);


############################################################

def makeCanvas(hists, names, canname, isLog = False):

    max = -999.;
    for hist in hists:
        if max < hist.GetMaximum(): max = hist.GetMaximum();

    leg = ROOT.TLegend(0.7,0.7,0.9,0.9);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    for i in range(len(names)):
        leg.AddEntry(h_mWW[i], names[i], "l");


    can = ROOT.TCanvas("can"+canname,"can"+canname,1000,800);
    hists[0].SetMaximum( 1.2*max );
    if not isLog: hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    if isLog: ROOT.gPad.SetLogy();
    can.SaveAs("figs/"+canname+".eps");
    can.SaveAs("figs/"+canname+".png");

    for hist in hists:
        hist.Scale(1./hist.Integral());

    cann = ROOT.TCanvas("cann"+canname,"cann"+canname,1000,800);
    hists[0].SetMaximum( 1.25*hists[0].GetMaximum() );
    if not isLog: hists[0].SetMinimum( 0 );    
    hists[0].Draw();    
    for i in range(1,len(hists)):
        hists[i].Draw("sames");
    leg.Draw();
    if isLog: ROOT.gPad.SetLogy();
    cann.SaveAs("figs/"+canname+"_norm.eps");
    cann.SaveAs("figs/"+canname+"_norm.png");

if __name__ == '__main__':

    prefix = "/uscms_data/d2/ntran/physics/WWScattering/LHEAnalyzers/trees/";

    ## polarization samples from Patricia - v3
    filenames = [prefix+"madgraph-v3_TOT_wplepwmhad.root",prefix+"madgraph_TT_WPlepWMhad.root",prefix+"madgraph-v3_LT_wplepwmhad.root",prefix+"madgraph-v3_LL_wplepwmhad.root"];
    names = ["TOT","TT","LT","LL"];
    crossSections = [0.0394*1000,0.235e-1*1000,0.0139*1000,0.0022*1000];

#    ## phantom samples from Pietro
#    filenames = [prefix+"phantom_EWK6_noh.root",prefix+"phantom_EWK6_h126.root"];
#    names = ["no higgs","w/ Higgs"];
#    crossSections = [7.299e-02*2*1000, 7.607e-02*2*1000];    
    
    colors = [1,2,4,6];
    #nEvents = [1000000.]*len(filenames);
    nEvents = [];
    
    vars = ["mWW"];
    
    h_mWW = [];
    h_mWW_aftercuts = [];    
    h_costheta1 = [];
    h_costheta2 = [];
    h_phi = [];  
    h_costhetastar = [];  
    h_phi1 = [];            
    h_dEtajj = [];
    h_mjj = [];
    
    for i in range(len(filenames)):
        h_mWW.append( ROOT.TH1F("h_mWW_"+str(i),"; m_{WW}; events/fb^{-1}", 50, 0, 1500 ) );
        h_mWW[i].SetLineColor(colors[i]); h_mWW[i].SetLineWidth(2);

        h_mWW_aftercuts.append( ROOT.TH1F("h_mWW_aftercuts_"+str(i),"; m_{WW}; events/fb^{-1}", 50, 0, 1500 ) );
        h_mWW_aftercuts[i].SetLineColor(colors[i]); h_mWW_aftercuts[i].SetLineWidth(2);

        h_costheta1.append( ROOT.TH1F("h_costheta1_"+str(i),"; costheta1; events/fb^{-1}", 20, -1, 1 ) );
        h_costheta1[i].SetLineColor(colors[i]); h_costheta1[i].SetLineWidth(2);
        
        h_costheta2.append( ROOT.TH1F("h_costheta2_"+str(i),"; costheta2; events/fb^{-1}", 20, -1, 1 ) );
        h_costheta2[i].SetLineColor(colors[i]); h_costheta2[i].SetLineWidth(2);
        
        h_phi.append( ROOT.TH1F("h_phi_"+str(i),"; phi; events/fb^{-1}", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi() ) );
        h_phi[i].SetLineColor(colors[i]); h_phi[i].SetLineWidth(2);

        h_costhetastar.append( ROOT.TH1F("h_costhetastar_"+str(i),"; costhetastar; events/fb^{-1}", 20, -1, 1 ) );
        h_costhetastar[i].SetLineColor(colors[i]); h_costhetastar[i].SetLineWidth(2);

        h_phi1.append( ROOT.TH1F("h_phi1_"+str(i),"; phi_{1}; events/fb^{-1}", 20, -ROOT.TMath.Pi(), ROOT.TMath.Pi() ) );
        h_phi1[i].SetLineColor(colors[i]); h_phi1[i].SetLineWidth(2);
        
        h_dEtajj.append( ROOT.TH1F("h_dEtajj_"+str(i),"; #Delta #eta; events/fb^{-1}", 20, 0, 6) );
        h_dEtajj[i].SetLineColor(colors[i]); h_dEtajj[i].SetLineWidth(2);        

        h_mjj.append( ROOT.TH1F("h_mjj_"+str(i),"; m_{jj}; events/fb^{-1}", 20, 0, 600) );
        h_mjj[i].SetLineColor(colors[i]); h_mjj[i].SetLineWidth(2);        


    # fill
    for i in range(len(filenames)):
        print "file: ", filenames[i]
        tmpf = ROOT.TFile(filenames[i]);
        tmpt = tmpf.Get("tree");
        
        nEvents.append( float(tmpt.GetEntries()) );
        
        for j in range(tmpt.GetEntries()):
            if j%100000 == 0: print "j = ", j
            if j > 200000000: break;
            tmpt.GetEntry(j);
            if tmpt.isSignal <= 1:
                h_mWW[i].Fill( tmpt.mWW, crossSections[i]/nEvents[i] );
                if tmpt.mWW > 700:
                    h_dEtajj[i].Fill( tmpt.dEtajj, crossSections[i]/nEvents[i] );                
                    h_mjj[i].Fill( tmpt.mjj, crossSections[i]/nEvents[i] );                                
                    if tmpt.dEtajj > 3 and tmpt.mjj > 300:
                        h_mWW_aftercuts[i].Fill( tmpt.mWW, crossSections[i]/nEvents[i] );
                        h_costheta1[i].Fill( tmpt.costheta1, crossSections[i]/nEvents[i] );
                        h_costheta2[i].Fill( tmpt.costheta2, crossSections[i]/nEvents[i] );                
                        h_phi[i].Fill( tmpt.phi, crossSections[i]/nEvents[i] );
                        h_costhetastar[i].Fill( tmpt.costhetastar, crossSections[i]/nEvents[i] );
                        h_phi1[i].Fill( tmpt.phi1, crossSections[i]/nEvents[i] );                                                
                
        del tmpf;
        del tmpt;
                
    # plot
    makeCanvas(h_mWW,names,"mWW",1);
    makeCanvas(h_mWW_aftercuts,names,"mWW_aftercuts");
    makeCanvas(h_costheta1,names,"costheta1");
    makeCanvas(h_costheta2,names,"costheta2");
    makeCanvas(h_phi,names,"phi");
    makeCanvas(h_costhetastar,names,"costhetastar");
    makeCanvas(h_phi1,names,"phi1");        
    makeCanvas(h_dEtajj,names,"dEtajj");
    makeCanvas(h_mjj,names,"mjj");
            
        
        
        

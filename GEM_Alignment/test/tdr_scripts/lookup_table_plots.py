import ROOT, tdrstyle, sys, os, array, datetime
now = datetime.datetime.now()

# ARGUMENTS ************ file, shift (-.1cm), direction (X), Zlow, Zhigh

run24 = "2024C"
run23 = "2023BC"
run22 = "2022G"
low_pt = 30
high_pt = 200

#txt_description = "dPhi=rechit local phi - prophit local phi [created in units of mRad]. dPhi corrected calculation takes the firmware long and short superchambers orientation into account on +/- endcaps, i.e.\n"+"        R+1    R-1\n"+"Even   -1       +1\n"+"Odd   +1       -1\n"+"Shift [strip(8)] is calculated by the following conversion rule: 0.37 [mRad] = 1 [strip(8)]\n"+"cuts: 30GeV<muon_pt<200GeV, has_fidcut, n_ME11_segment=1, |RdPhi_Corrected|<2cm, |dPhi_Corrected|<3mRad, has_TightID, endcap, station, superchamber, eta\n"+"[filename for TA's notekeeping]: ntuple: Run2022D_ZMu_PromptReco_RAWRECO_globalMu_pfisotight_idealGEMonCSC_v1.root]\n"+"=================================================================================================================\n"

f24 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2024/Run{r}_LUT_v1.root".format(r=run24))
f23 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2023/{r}_LUT_v0.root".format(r=run23))
f22 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2022/{r}/Run{r}_ZMu_PromptReco_RAWRECO_idealGEMidealCSC_v1.root".format(r=run22))
event24 = f24.Get("analyzer/ME11Seg_Prop")
event23 = f23.Get("analyzer/ME11Seg_Prop")
event22 = f22.Get("analyzer/ME11Seg_Prop")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

nX_row = 9
nY_col = 4

H_ref = 800
W_ref = 800
W = W_ref# * nX_row
H = H_ref# * nY_col

T = 0.12*H
B = 0.16*H
L = 0.16*W
R = 0.08*W

xbins = 100
xlow = -4
xhigh = 4


canvas = ROOT.TCanvas("c1", "c1", 100, 100, W, H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)

def strip_conversion(phi): #0.37 [mRad] = 1 [strip(8)]
  strip = round(phi/0.37, 0)
  return strip

'''#below excerpt is for all plots on one canvas is needed
name_mark = ROOT.TText(0.89, 0.25, "Towsifa")
name_mark.SetTextSize(0.025)
name_mark.SetTextColor(ROOT.kViolet+1)
name_mark.Draw()
canvas.Divide(nX_row, nY_col)
'''
for reg in [-1, 1]:
  if reg>0:
    reg_string="+"+str(abs(reg))
  elif reg<0:
    reg_string="-"+str(abs(reg))


  t = open("Run{run}_LUT_R{r}_st1.txt".format(run=run24, r=reg), "w")
  t.writelines("SuperchamberEta | Shift [strip(8)]\n")

  ### .xml file format
  x0 = open("Run{run}_LUT_R{r}_st1.xml".format(run=run24, r=reg), "w")
  x0.writelines('<?xml version="1.0" encoding="utf-8"?>\n')
  x0.writelines('<GEMAlignment generationTime="'+ now.strftime("%Y-%m-%d %H:%M:%S")+'" dataSet="2024C">\n')

  for st in [1]: #[1, 2] for both stations
    for ch in range(1, 37):
      if ch<10:
        x0.writelines(' <CSC label="ME{r}/{st}/0{ch}">\n'.format(r=reg_string, st=st, ch=ch))
      else:
        x0.writelines(' <CSC label="ME{r}/{st}/{ch}">\n'.format(r=reg_string, st=st, ch=ch))
      x0.writelines('   <TMB\n')
      h_index = -1
      h1_index = -1
      h2_index = -1
      h_list = []
      h1_list = []
      h2_list = []
      for eta in range(1, 9):
        h_index+=1
        h1_index+=1
        h2_index+=1
        if st==2:
          t.writelines("{ch}{e} {strip}\n".format(ch=ch, e=eta, strip=0))
          continue
        
        #canvas.cd(ch)
        #ROOT.gPad.SetGridx()
        #ROOT.gPad.SetGridy()
        h_list.append(ROOT.TH1D("h_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "h_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), xbins, xlow, xhigh))
        #h1_list.append(ROOT.TH1D("h1_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "h1_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), xbins, xlow, xhigh))
        #h2_list.append(ROOT.TH1D("h2_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "h2_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), xbins, xlow, xhigh))

        xAxis = h_list[h_index].GetXaxis()
        xAxis.SetTitleOffset(1)
        xAxis.SetTitleSize(0.05)
        xAxis.SetTitle("#Delta#phi [mRad]")
        #xAxis.CenterTitle()

        yAxis = h_list[h_index].GetYaxis()
        yAxis.SetTitleOffset(0)
        yAxis.SetTitleSize(0.05)
        #yAxis.SetTitle("A.U.")
        yAxis.SetTitle("Entries")
        #yAxis.CenterTitle()

        event24.Project("h_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "1000*dPhi_Corrected", "abs(RdPhi_Corrected) < 2 && 1000*abs(dPhi_Corrected)<4 && n_ME11_segment == 1 && has_fidcut && muon_pt > {low} && muon_pt<{high} && has_TightID==1 && prop_location[0]=={reg} && prop_location[1]=={st} && prop_location[2]=={ch} && prop_location[4]=={e}".format(reg=reg, st=st, ch=ch, e=eta, low=low_pt, high=high_pt))
        #event1.Project("h1_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "1000*dPhi_Corrected", "abs(RdPhi_Corrected)< 2 && 1000*abs(dPhi_Corrected)<4 && n_ME11_segment == 1 && has_fidcut && muon_pt > {low} && muon_pt<{high} && has_TightID==1 && prop_location[0]=={reg} && prop_location[1]=={st} && prop_location[2]=={ch} && prop_location[4]=={e}".format(reg=reg, st=st, ch=ch, e=eta, low=low_pt, high=high_pt))

        #event.Project("h1_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta), "RdPhi_Corrected_mrad", "abs(RdPhi_Corrected) < 2 && abs(RdPhi_Corrected_mrad)<4 && n_ME11_segment == 1 && has_fidcut && muon_pt > {low} && muon_pt<{high} && has_TightID==1 && prop_location[0]=={reg} && prop_location[1]=={st} && prop_location[2]=={ch} && prop_location[3]==1 && prop_location[4]=={e}".format(reg=reg, st=st, ch=ch, e=eta, low=low_pt, high=high_pt))
        #event.Project("h2_{r}{st}{sc}eta{e}".format(r=reg, st=st, sc=ch, e=eta),"RdPhi_Corrected_mrad", "abs(RdPhi_Corrected) < 2 && abs(RdPhi_Corrected_mrad)<4 && n_ME11_segment == 1 && has_fidcut && muon_pt > {low} && muon_pt<{high} && has_TightID==1 && prop_location[0]=={reg} && prop_location[1]=={st} && prop_location[2]=={ch} && prop_location[3]==2 && prop_location[4]=={e}".format(reg=reg, st=st, ch=ch, e=eta, low=low_pt, high=high_pt))

        h=h_list[h_index]
        #h1=h1_list[h1_index]
        #h2=h2_list[h2_index]

        if (h.Integral() >0):# and h1.Integral() >0):
          h.Scale(1/h.Integral())
          #h1.Scale(1/h1.Integral())

        h.SetLineWidth(3)
        #h1.SetLineWidth(3)
        #h2.SetLineWidth(2)

        h.SetLineStyle(1)
        #h1.SetLineStyle(1)
        #h2.SetLineStyle(1)

        h.SetMarkerSize(0)
        #h1.SetMarkerSize(0)
        #h2.SetMarkerSize(0)

        h.SetLineColor(ROOT.kViolet+2)
        #h1.SetLineColor(ROOT.kGreen+2)
        #h2.SetLineColor(ROOT.kRed)

        yAxis.SetRangeUser(0, 1.1*max(h.GetMaximum()))#, h1.GetMaximum()))

        #h.Draw("hist")
        #h1.Draw("hist same")
        #h2.Draw("hist same")

        rootkde = ROOT.TLegend(0.6,0.75,0.9,0.87)
        rootkde.AddEntry(h, "2024C: mean: {mean}".format(mean=round(h.GetMean(), 3)))
        #rootkde.AddEntry(h1, "362695: mean: {mean}".format(mean=round(h1.GetMean(), 3)))
        #rootkde.AddEntry(h2, "Layer 2 mean: {mean}".format(mean=round(h2.GetMean(), 3)))
        rootkde.SetBorderSize(0)
        rootkde.Draw()

        #t.writelines("{ch}{e} {strip0} {strip1} {diff}\n".format(ch=ch, e=eta, strip0=round(strip_conversion(h.GetMean())), strip1=round(strip_conversion(h1.GetMean())), diff=round(strip_conversion(h.GetMean()))-round(strip_conversion(h1.GetMean())) ))
        t.writelines("{ch}{e} {strip}\n".format(ch=ch, e=eta, strip=round(strip_conversion(h.GetMean())) ))
        #t1.writelines("{ch}{e} {strip}\n".format(ch=ch, e=eta, strip=round(strip_conversion(h1.GetMean())) ))

        x0.writelines('   gem_xshift_eta'+str(eta-1)+'="{strip}"\n'.format(strip=round(strip_conversion(h.GetMean())) ))
        #x1.writelines('   gem_xshift_eta'+str(eta-1)+'="{strip}"\n'.format(strip=round(strip_conversion(h1.GetMean())) ))

        #t.writelines("{reg} {st} {ch} {e} {strip} {h} {h1} {h2}\n".format(reg=reg, st=st, ch=ch, e=eta, strip=strip_conversion(h.GetMean()), h=h.GetMean(), h1=h1.GetMean(), h2=h2.GetMean()))
        #t1.writelines("{reg} {st} {ch} {e} {strip} {h}\n".format(reg=reg, st=st, ch=ch, e=eta, strip=strip_conversion(h.GetMean()), h=round(h.GetMean(), 3)))
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.4*canvas.GetTopMargin())

        latex.SetTextAlign(32)
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "(13.6 TeV)")
        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.2*canvas.GetTopMargin(), "Run {run}".format(run=run))

        latex.SetTextSize(0.25*canvas.GetTopMargin())
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.6*canvas.GetTopMargin(), "Region{reg}, Station{st}, SC{sc}, eta{e}".format(reg=reg, st=st, sc=ch, e=eta))

        latex.SetTextSize(0.5*canvas.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.4*canvas.GetTopMargin())
        latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")

        latex.SetTextFont(42)
        latex.SetTextSize(0.4*canvas.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()
        #if os.path.exists("Run{y}/dphi_plots/R{r}/st{st}/SC{sc}".format(y=run, r=reg, st=st, sc=ch)) == False:
        #  os.mkdir("Run{y}/dphi_plots/R{r}/st{st}/SC{sc}".format(y=run, r=reg, st=st, sc=ch))
        #canvas.SaveAs("Run{y}/dphi_plots/R{r}/st{st}/SC{sc}/dPhi_R{r}_st{st}_SC{sc}_eta{e}.png".format(y=run, e=eta, r=reg, st=st, sc=ch))
        #canvas.SaveAs("test{r}{ch}{eta}.png".format(r=reg,ch=ch, eta=eta))
      x0.writelines('   />\n')
      #x1.writelines('   />\n')
      x0.writelines(' </CSC>\n')
      #x1.writelines(' </CSC>\n')

  x0.writelines('</GEMAlignment>')
  #x1.writelines('</GEMAlignment>')

t.close()
#t0.close()
#t1.close()

x0.close()
#x1.close()

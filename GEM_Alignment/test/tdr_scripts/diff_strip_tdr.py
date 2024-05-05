import ROOT, tdrstyle, sys, os, array
#Following tutorial available on https://root.cern/doc/master/ratioplotOld_8C.html


y_var = "strip" #"strip" or "rdphi" or "BA"
low_pt=30
high_pt=200
endcap = -1
layer = 1
folder = ""
run24 = "2024C"
run23 = "2023BC"
run22 = "2022G"

f24 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2024/Run{r}_LUT_v1.root".format(r=run24))
f23 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2023/{r}_LUT_v0.root".format(r=run23))
f22 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2022/{r}/Run{r}_ZMu_PromptReco_RAWRECO_idealGEMidealCSC_v1.root".format(r=run22))
event24 = f24.Get("analyzer/ME11Seg_Prop")
event23 = f23.Get("analyzer/ME11Seg_Prop")
event22 = f22.Get("analyzer/ME11Seg_Prop")


def strip_conversion(phi): #0.346 [mRad] = 1 [strip(8)]
  strip = round(phi/0.346, 0)
  return strip

if y_var == "BA":
  ylow = -5
  yhigh = 5
  y_plot = "1000*bending_angle"
  y_axis = "Bending Angle [mrad]"
elif y_var == "rdphi":
  ylow = -2
  yhigh = 2
  y_plot = "RdPhi_Corrected"
  y_axis = "#DeltaR#phi [cm]"
elif y_var == "strip":
  ylow = strip_conversion(-5.5)
  yhigh = strip_conversion(5.5)
  y_plot = "1000*dPhi_Corrected"
  y_axis = "#Delta 1/8 strip"

if endcap==1:
  reg_string="+"
elif endcap==-1:
  reg_string="-"

ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

H_ref = 800
W_ref = 1200
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.20*W_ref
R = 0.20*W_ref

ybins = 100
xbins = 36
xlow = 0.5
xhigh = 36.5

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
canvas.SetGrid()
#canvas.SetLogy()

#h2023 = ROOT.TH1D("h2023", ";GE1/1 Superchambers; {y_axis}".format(y_axis=y_axis), xbins, xlow, xhigh)
h2024 = ROOT.TH1D("h2024", "h2024", xbins, xlow, xhigh)
h2023 = ROOT.TH1D("h2023", "h2023", xbins, xlow, xhigh)
h2022 = ROOT.TH1D("h2022", "h2022", xbins, xlow, xhigh)
hdiff = ROOT.TH1D("hdiff", "hdiff", xbins, xlow, xhigh)

#layer = prop_location[3]
#chamber = prop_location[2]
#region or endcap = prop_location[0]

h2024_tmp = []
h2023_tmp = []
h2022_tmp = []
for ch in range(1, 37):
  h2024_tmp.append(ROOT.TH1D("24{ch}".format(ch=ch), "24{ch}".format(ch=ch), 300, ylow, yhigh))
  h2023_tmp.append(ROOT.TH1D("23{ch}".format(ch=ch), "23{ch}".format(ch=ch), 300, ylow, yhigh))
  h2022_tmp.append(ROOT.TH1D("22{ch}".format(ch=ch), "22{ch}".format(ch=ch), 300, ylow, yhigh))
  event24.Project("24{ch}".format(ch=ch), "{x}".format(x=y_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[2]=={ch} && prop_location[3]=={lay}".format(low=low_pt, high=high_pt, reg=endcap, ch=ch, lay=layer))
  event23.Project("23{ch}".format(ch=ch), "{x}".format(x=y_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[2]=={ch} && prop_location[3]=={lay}".format(low=low_pt, high=high_pt, reg=endcap, ch=ch, lay=layer))
  event22.Project("22{ch}".format(ch=ch), "{x}".format(x=y_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[2]=={ch} && prop_location[3]=={lay}".format(low=low_pt, high=high_pt, reg=endcap, ch=ch, lay=layer))

  h2024.SetBinContent(ch, strip_conversion(h2024_tmp[ch-1].GetMean()))
  h2024.SetBinError(ch, strip_conversion(h2024_tmp[ch-1].GetMeanError()))
  h2023.SetBinContent(ch, strip_conversion(h2023_tmp[ch-1].GetMean()))
  h2023.SetBinError(ch, strip_conversion(h2023_tmp[ch-1].GetMeanError()))
  h2022.SetBinContent(ch, strip_conversion(h2022_tmp[ch-1].GetMean()))
  h2022.SetBinError(ch, strip_conversion(h2022_tmp[ch-1].GetMeanError()))

upPad = ROOT.TPad("up", "up", 0, 0.3, 1, 1)
upPad.SetBottomMargin(0)
upPad.SetGrid()
upPad.Draw()
upPad.cd()

h2024.Draw("EP")
h2023.Draw("EP SAME")
h2022.Draw("EP SAME")

xAxis = h2024.GetXaxis()
xAxis.SetLabelSize(0.12)
xAxis.SetTitleOffset(1)
xAxis.SetTitleSize(0.15)
xAxis.SetTitle("GE1/1 Superchambers")
#xAxis.CenterTitle()

yAxis = h2024.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
if y_var == "BA":
  yAxis.SetNdivisions(-505)
yAxis.SetTitle("{y_axis}".format(y_axis=y_axis))
yAxis.SetRangeUser(ylow, yhigh)
#yAxis.CenterTitle()

canvas.cd()
lowPad = ROOT.TPad("lp", "lp", 0, 0.05, 1, 0.28)
lowPad.SetTopMargin(0)
lowPad.SetBottomMargin(0.3)
lowPad.SetGrid()
lowPad.Draw()
lowPad.cd()
hdiff = h2024.Clone()
hdiff.SetLineColor(ROOT.kBlack)
#hdiff.Sumw2()
hdiff.Add(h2023, -1)
hdiff.Draw("EP")

diff_y = hdiff.GetYaxis()
diff_y.SetLabelSize(0.09)
diff_y.SetTitleOffset(0)
diff_y.SetTitleSize(0.1)
diff_y.SetNdivisions(505)
diff_y.SetTitle("2024C-2023BC")

h2024.SetLineWidth(3)
h2023.SetLineWidth(3) #3
h2022.SetLineWidth(3)

h2024.SetMarkerStyle(4)
h2023.SetMarkerStyle(2)
h2022.SetMarkerStyle(5)

h2024.SetMarkerSize(3)
h2023.SetMarkerSize(3)
h2022.SetMarkerSize(3)

h2024.SetLineColor(ROOT.kGreen+2)
h2023.SetLineColor(ROOT.kBlue)
h2022.SetLineColor(ROOT.kRed)

h2024.SetMarkerColor(ROOT.kGreen+2)
h2023.SetMarkerColor(ROOT.kBlue)
h2022.SetMarkerColor(ROOT.kRed)

canvas.cd()
upPad.cd()
legend = ROOT.TLegend(0.8, 0.7, 0.9, 0.85)
legend.AddEntry(h2024, "2024 C")
legend.AddEntry(h2023, "2023 BC")
legend.AddEntry(h2022, "2022 G")
legend.SetTextSize(0.)
legend.SetBorderSize(0)
legend.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)

latex.SetTextFont(42)
latex.SetTextSize(0.4*canvas.GetTopMargin())

latex.SetTextAlign(5)
latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "(13.6 TeV)")
latex.SetTextAlign(22)

latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.DrawLatex(1-1.8*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.1*canvas.GetTopMargin(), "{low} GeV".format(low=low_pt)+" < p_{T}^{GLB} < "+"{high} GeV".format(high=high_pt)) 
latex.DrawLatex(1-1.8*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.5*canvas.GetTopMargin(), "{reg}Endcap Layer{lay} chambers".format(reg=reg_string, lay=layer))

latex.SetTextSize(0.5*canvas.GetTopMargin())
latex.SetTextFont(61)
latex.SetTextAlign(42)
latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
#latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.27*canvas.GetTopMargin(), "CMS")
latex.SetTextFont(52)
latex.SetTextAlign(48)
latex.SetTextSize(0.4*canvas.GetTopMargin())
latex.DrawLatex(1.4*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
#latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Preliminary")
#latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Work in Progress")

latex.SetTextFont(42)
latex.SetTextSize(0.4*canvas.GetTopMargin())

###fit function part below
latex.SetTextAlign(12)
latex.SetTextSize(0.03)
latex.SetTextFont(61)
#latex.DrawLatex(0.3*canvas.GetLeftMargin(), -0.04+ 0.5*canvas.GetBottomMargin(), "{zer}* exp(-0.5*((x- {one}) / {two})^2) + {three}* exp(-0.5*((x- {four}) / {five})^2)".format(zer = round(f1.GetParameter(0), 3), one = round(f1.GetParameter(1), 3), two = round(f1.GetParameter(2), 3), three = round(f1.GetParameter(3), 3), four = round(f1.GetParameter(4), 3), five = round(f1.GetParameter(5), 3),))

###


frame = canvas.GetFrame()
frame.Draw()


if os.path.exists("{folder}".format(folder=folder)) == False:
  os.mkdir("{folder}".format(folder=folder))
canvas.SaveAs("{folder}/stripDiff_pt{low}to{high}_R{reg}L{lay}.png".format(folder=folder, low=low_pt, high=high_pt, reg=endcap, lay=layer))

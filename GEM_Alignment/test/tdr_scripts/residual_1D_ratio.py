import ROOT, tdrstyle, sys, os, array

low_pt=20
high_pt=200
supch_cut = "both" #0 even, 1 odd, for both "both"
endcap = 1
layer = 1
#charge = -1
x_var = "rdphi" #"BA" or "rdphi" 

year = "2026B" #change between 2022B and 2022C for example

f = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2026/Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_2025alignment_trackerprop.root")
f1 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2026/Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_trackerprop_aligned.root")
#fMC = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2022/singleMuonGun_11_3_4_2021_design_v0.root")

#event0 = f0.Get("analyzer/ME11SegReco_Prop")
#event1 = f0.Get("analyzer/InnerRefit_Prop")
event = f.Get("analyzer/Inner_Prop")
event1 = f1.Get("analyzer/Inner_Prop")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

if x_var == "BA":
  xlow = -10
  xhigh = 10
  x_plot = "1000*bending_angle"
  x_axis = "Bending Angle [mrad]"
elif x_var == "rdphi":
  xlow = -2
  xhigh = 2
  x_plot = "RdPhi_Corrected"
  x_axis = "#DeltaR#phi [cm]"

if endcap==1:
  reg_string="+"
elif endcap==-1:
  reg_string="-"
'''
if charge==1:
  mu_str = "posmu"
  mu_notation = "#mu^{+}"
elif charge==-1:
  mu_str = "negmu"
  mu_notation ="#mu^{-}"
'''
H_ref = 800
W_ref = 800
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.16*W_ref
R = 0.08*W_ref

xbins = 100
ybins = 100
xlow = -2
xhigh = 2
ylow = -300
yhigh = 300
zlow = -2
zhigh = 2

canvas = ROOT.TCanvas("c1", "c1", 800, 800)
pad_main = ROOT.TPad("pad_main", "pad_main", 0, 0.3, 1, 1)
pad_ratio = ROOT.TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.3)
canvas.SetLeftMargin(0.16)
canvas.SetRightMargin(0.08)
canvas.SetTopMargin(0.12)
canvas.SetBottomMargin(0.16)
pad_main.SetBottomMargin(0)
pad_ratio.SetTopMargin(0)
pad_ratio.SetBottomMargin(0.3)
pad_main.Draw()
pad_ratio.Draw()

pad_main.cd()
pad_main.SetGrid()

h = ROOT.TH1D("h", "h", xbins, xlow, xhigh)
h1 = ROOT.TH1D("h1", "h1", xbins, xlow, xhigh)
#h2 = ROOT.TH1D("h2", "h2", xbins, xlow, xhigh)
xAxis = h.GetXaxis()
xAxis.SetTitleOffset(0)
xAxis.SetTitleSize(0.05)
#xAxis.SetNdivisions(-505)
xAxis.SetTitle("#DeltaR#phi [cm]")
#xAxis.SetTitle("|bending angle| [mrad]")
#xAxis.SetTitle("p_{T} [GeV]")
#xAxis.CenterTitle()

yAxis = h.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
#yAxis.SetTitle("A.U.")
yAxis.SetTitle("Normalized Events")
#yAxis.CenterTitle()

ch_list = []
ch_list_even = []
ch_list_odd = []
even_cut = ""
odd_cut = ""

for i in range(1,37):
  ch_list.append(str(i))
  if i%2 == 0:
    ch_list_even.append(str(i))
    if i==2:
      even_cut += "(prop_location[2]== {i} ||".format(i=i)
    elif i==36:
      even_cut += " prop_location[2] == {i})".format(i=i)
    else:
      even_cut += " prop_location[2] == {i} ||".format(i=i)
  else:
    ch_list_odd.append(str(i))
    if i==1:
      odd_cut += "(prop_location[2]== {i} ||".format(i=i)
    elif i==35:
      odd_cut += " prop_location[2] == {i})".format(i=i)
    else:
      odd_cut += " prop_location[2] == {i} ||".format(i=i)

if supch_cut==0:
  cut = even_cut
  ch_string = "Even"
elif supch_cut==1:
  cut = odd_cut
  ch_string = "Odd"
elif supch_cut=="both":
  ch_string = "All"

#print(even_cut)
#print(odd_cut)
#&& n_ME11_segment==1
event.Project("h", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer))
event1.Project("h1", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer))

#h.ResetStats()
#h.GetSumOfWeights()
h.Scale(1/h.Integral())
h1.Scale(1/h1.Integral())
#h2.Scale(1/h2.Integral())

h.SetLineWidth(3) #3
h1.SetLineWidth(3)
# h2.SetLineWidth(3)

h.SetMarkerSize(0)
h1.SetMarkerSize(0)
# h2.SetMarkerSize(0)

h.SetLineColor(ROOT.kBlue)
h1.SetLineColor(ROOT.kGreen+2)
# h2.SetLineColor(ROOT.kRed)
#h.SetFillColorAlpha(ROOT.kBlue, 0.3)
#h1.SetFillColorAlpha(ROOT.kGreen+2, 0.3)
#h2.SetFillColorAlpha(ROOT.kRed, 0.3)

yAxis.SetRangeUser(0, 1.6*max(h.GetMaximum(), h1.GetMaximum() ))
yAxis.SetMaxDigits(3)


h.Draw("HIST")
h1.Draw("HIST SAME")
# h2.Draw("HIST SAME")

# gaussian fits
ROOT.gStyle.SetOptFit(0)
t = ROOT.TF1("t", "gaus", -0.5, 0.5)
t1 = ROOT.TF1("t1", "gaus", -0.5, 0.5)
h.Fit(t)
h1.Fit(t1)

legend = ROOT.TLegend(0.5, 0.75, 0.95, 0.85)
legend.AddEntry(h, f"GEM alignment from 2025: {t.GetParameter(1):.3f} #pm {t.GetParameter(2):.3f}")
legend.AddEntry(h1, f"GEM alignment from 2026: {t1.GetParameter(1):.3f} #pm {t1.GetParameter(2):.3f}")
legend.SetTextSize(0.)
legend.SetBorderSize(0)
legend.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)

latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())

latex.SetTextAlign(32)
latex.DrawLatex(1-0.5*canvas.GetRightMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "(13.6 TeV)")
latex.SetTextAlign(12)

latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.0*canvas.GetTopMargin(), "Run {run}".format(run=year))

latex.SetTextSize(0.2*canvas.GetTopMargin())
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.5*canvas.GetTopMargin(), "{low} GeV".format(low=low_pt)+" < p_{T}^{GLB} < "+"{high} GeV".format(high=high_pt)) 
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.8*canvas.GetTopMargin(), "{reg}Endcap Layer {lay}".format(reg=reg_string, lay=layer))
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-2.1*canvas.GetTopMargin(), "{ch_string} chambers".format(ch_string=ch_string))

latex.SetTextSize(0.4*canvas.GetTopMargin())
latex.SetTextFont(61)
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
latex.SetTextFont(52)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Preliminary")

frame = canvas.GetFrame()
frame.Draw()


# Ratio plot
pad_ratio.cd()
pad_ratio.SetGrid()
h_ratio = h.Clone("h_ratio")
h_ratio.Divide(h1)
h_ratio.SetLineWidth(1)
h_ratio.SetLineColor(ROOT.kBlack)
h_ratio.SetMarkerSize(0)
h_ratio.GetXaxis().SetTitle("#DeltaR#phi [cm]")
h_ratio.GetXaxis().SetTitleSize(0.09)
h_ratio.GetXaxis().SetLabelSize(0.09)
h_ratio.GetYaxis().SetTitleOffset(0)
h_ratio.GetYaxis().SetTitle("Ratio")
h_ratio.GetYaxis().SetTitleSize(0.09)
h_ratio.GetYaxis().SetRangeUser(0.5, 1.5)
h_ratio.GetYaxis().SetNdivisions(505)
h_ratio.GetYaxis().SetLabelSize(0.09)
h_ratio.GetYaxis().SetTitleOffset(0)
h_ratio.Draw("E")
# end of ratio plot

canvas.SaveAs("Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_trackerprop/Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_trackerprop_{x_var}_1D_R{reg}L{lay}_{ch_string}chambers.png".format(reg=endcap, ch_string=ch_string, lay=layer, x_var=x_var))
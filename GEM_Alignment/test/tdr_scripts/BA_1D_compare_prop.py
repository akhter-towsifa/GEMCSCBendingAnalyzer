import ROOT, tdrstyle, sys, os, array

low_pt=30
high_pt=40
supch_cut = 1 #0 even, 1 odd, for both "both"
endcap = 1
layer = 1
charge = -1
x_var = "BA" #"BA" or "rdphi" 

year = "2024H" #change between 2022B and 2022C for example

fdata = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2024/2024H/Run2024H_muon0_150X_dataRun3_Prompt_v1.root")
fMC = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2022/singleMuonGun_11_3_4_2021_design_v0.root")

eventme11 = fdata.Get("analyzer/ME11SegReco_Prop")
eventtrack = fdata.Get("analyzer/InnerRefit_Prop")
eventMCme11 = fMC.Get("analyzer/ME11Seg_Prop")
eventMCtrack = fMC.Get("analyzer/Inner_Prop")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

if x_var == "BA":
  xlow = -5
  xhigh = 5
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

if charge==1:
  mu_str = "posmu"
  mu_notation = "#mu^{+}"
elif charge==-1:
  mu_str = "negmu"
  mu_notation ="#mu^{-}"

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
ylow = -300
yhigh = 300
zlow = -2
zhigh = 2

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

#h0 = ROOT.TH1D("h0", "h0", xbins, xlow, xhigh)
hMCme11 = ROOT.TH1D("hMCme11", "hMCme11", xbins, xlow, xhigh)
hMCtrack = ROOT.TH1D("hMCtrack", "hMCtrack", xbins, xlow, xhigh)
hme11 = ROOT.TH1D("hme11", "hme11", xbins, xlow, xhigh)
htrack = ROOT.TH1D("htrack", "htrack", xbins, xlow, xhigh)
xAxis = hMCme11.GetXaxis()
xAxis.SetTitleOffset(0)
xAxis.SetTitleSize(0.05)
#xAxis.SetNdivisions(-505)
# xAxis.SetTitle("#DeltaR#phi [cm]")
xAxis.SetTitle(f"{x_axis}")
#xAxis.SetTitle("p_{T} [GeV]")
#xAxis.CenterTitle()

yAxis = hMCme11.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitleSize(0.05)
yAxis.SetTitle("A.U.")
# yAxis.SetTitle("Entries")
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

# eventMCme11.Project("hMCme11", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, cut=cut))
# eventMCtrack.Project("hMCtrak", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, cut=cut))
# eventme11.Project("hme11", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, cut=cut))
# eventtrack.Project("htrack", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, cut=cut))
eventMCme11.Project("hMCme11", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && muon_charge=={charge} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, charge=charge, cut=cut))
eventMCtrack.Project("hMCtrak", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && muon_charge=={charge} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, charge=charge, cut=cut))
eventme11.Project("hme11", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && muon_charge=={charge} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, charge=charge, cut=cut))
eventtrack.Project("htrack", "{x}".format(x=x_plot), "muon_pt>{low} && muon_pt<{high} && n_ME11_segment==1 && has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[0]=={reg} && prop_location[3]=={lay} && muon_charge=={charge} && {cut}".format(low=low_pt, high=high_pt, reg=endcap, lay=layer, charge=charge, cut=cut))


#h.ResetStats()
#h.GetSumOfWeights()
# hMCme11.Scale(1/hMCme11.Integral())
# hMCtrack.Scale(1/hMCtrack.Integral())
# hme11.Scale(1/hme11.Integral())
# htrack.Scale(1/htrack.Integral())

hMCme11.SetLineWidth(3) #3
hMCtrack.SetLineWidth(3)
hme11.SetLineWidth(3)
htrack.SetLineWidth(3)

hMCme11.SetLineStyle(1) 
hMCtrack.SetLineStyle(2)
hme11.SetLineStyle(1)
htrack.SetLineStyle(2)

hMCme11.SetMarkerSize(0)
hMCtrack.SetMarkerSize(0)
hme11.SetMarkerSize(0)
htrack.SetMarkerSize(0)


hMCme11.SetLineColor(ROOT.kRed)
hMCtrack.SetLineColor(ROOT.kRed)
hme11.SetLineColor(ROOT.kBlue)
htrack.SetLineColor(ROOT.kGreen+2)
#h.SetFillColorAlpha(ROOT.kBlue, 0.3)
#h1.SetFillColorAlpha(ROOT.kGreen+2, 0.3)
#h2.SetFillColorAlpha(ROOT.kRed, 0.3)

#yAxis.SetRangeUser(0, 1.1*h.GetMaximum())
yAxis.SetRangeUser(0, 1.6*max(hMCme11.GetMaximum(), hMCtrack.GetMaximum(), hme11.GetMaximum(), htrack.GetMaximum() ))
yAxis.SetMaxDigits(3)

hMCme11.Draw("HIST")
hMCtrack.Draw("HIST SAME")
hme11.Draw("HIST SAME")
htrack.Draw("HIST SAME")

#f1 = ROOT.TF1("f1", "[0]* exp([1]*x**2 + [2]) + [3]* exp([4]*x**2 + [5])", -1.5, 1.5)
#f1 = ROOT.TF1("f1", "[0]* exp(-0.5*((x-[1])/[2])**2) + [3]* exp(-0.5*((x-[4])/[5])**2)", gaus_low, gaus_high)
#f1 = ROOT.TF1("f1", "gaus(0)+gaus(3)", gaus_low, gaus_high)
#f1.SetParameters(h.GetMaximum(), h.GetMean(), h.GetStdDev(), h.GetMaximum(), h.GetMean(), h.GetStdDev())
#f1.SetParameters(.1,.1,.1,.1,.1,.1)
#f1.SetLineColor(ROOT.kRed)
#f1.SetMarkerSize(0)
#h.Fit("f1")
#f1.Draw("same")

#legend.AddEntry(h, "mean: {m}".format(m=round(h.GetMean(), 3)))
#legend.AddEntry(h, "std dev: {s}".format(s=round(h.GetStdDev(), 3)))
#legend.AddEntry(f1, "mean: {m}".format(m=round(f1.GetParameter(1), 3)))
#legend.AddEntry(f1, "std dev: {s}".format(s=round(f1.GetParameter(2), 3)))

legend = ROOT.TLegend(0.5, 0.7, 0.9, 0.85)
legend.AddEntry(hMCme11, f"MC back prop")
legend.AddEntry(hMCtrack, f"MC track prop")
legend.AddEntry(hme11, f"data: back prop")
legend.AddEntry(htrack, f"data: track prop")
# legend.AddEntry(hMCme11, f"MC back prop: {hMCme11.GetMean():.3f} #pm {hMCme11.GetStdDev():.3f}")
# legend.AddEntry(hMCtrack, f"MC track prop: {hMCtrack.GetMean():.3f} #pm {hMCtrack.GetStdDev():.3f}")
# legend.AddEntry(hme11, f"data: back prop: {hme11.GetMean():.3f} #pm {hme11.GetStdDev():.3f}")
# legend.AddEntry(htrack, f"data: track prop: {htrack.GetMean():.3f} #pm {htrack.GetStdDev():.3f}")
# legend.AddEntry(h2, "75 GeV < p_{T} < 200 GeV")
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
#latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h.GetEntries())))
#latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h.GetMean(),3)))
#latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h.GetStdDev(),3)))
latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "(13.6 TeV)")
latex.SetTextAlign(12)

latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.9*canvas.GetTopMargin(), "Run {run}".format(run=year))

latex.SetTextSize(0.25*canvas.GetTopMargin())
#latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.2*canvas.GetTopMargin(), "{low} GeV".format(low=low_pt)+" < p_{T}^{GLB} < "+"{high} GeV".format(high=high_pt)) 
latex.DrawLatex(0.65-0.3*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-1.6*canvas.GetTopMargin(), "{reg}Endcap Layer {lay}".format(reg=reg_string, lay=layer))
latex.DrawLatex(0.65-0.3*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-1.9*canvas.GetTopMargin(), "{ch_string} chambers".format(ch_string=ch_string))
latex.DrawLatex(0.55-0.3*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-2.1*canvas.GetTopMargin(), f"{mu_notation}")

latex.SetTextSize(0.5*canvas.GetTopMargin())
latex.SetTextFont(61)
#latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.27*canvas.GetTopMargin(), "CMS")
latex.SetTextFont(52)
latex.SetTextSize(0.3*canvas.GetTopMargin())
#latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Preliminary")
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


#if os.path.exists("Run{run}/{version}".format(run=year, version=version)) == False:
#  os.mkdir("Run{run}/{version}".format(run=year, version=version))
#canvas.SaveAs("Run{run}/{version}/1D_RdPhi_pt{low}to{high}_{ch_string}chambers_R{reg}_L{lay}.png".format(run=year, version=version, low=low_pt, high=high_pt, reg=endcap, lay=layer, ch_string=ch_string))
# canvas.SaveAs("{x_var}_1D_R{reg}L{lay}_{ch_string}chambers.png".format(reg=endcap, ch_string=ch_string, lay=layer, x_var=x_var))
canvas.SaveAs("{x_var}_1D_R{reg}L{lay}_{ch_string}chambers_{mu_str}.png".format(reg=endcap, ch_string=ch_string, lay=layer, x_var=x_var, mu_str=mu_str))

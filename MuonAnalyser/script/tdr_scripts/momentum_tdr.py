import ROOT, tdrstyle, sys, os, array

# ARGUMENTS ************ file, shift (-.1cm), direction (X), Zlow, Zhigh

f = ROOT.TFile("{}".format(sys.argv[1]))
event = f.Get("analyser/ME11Seg_Prop")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

H_ref = 600
W_ref = 800
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.16*W_ref
R = 0.08*W_ref

xbins = 100
ybins = 100
xlow = 0
xhigh = 50
ylow = -300
yhigh = 300

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

h = ROOT.TH1D("h", "h", xbins, xlow, xhigh)

xAxis = h.GetXaxis()
xAxis.SetTitle("pT [GeV]")
#xAxis.CenterTitle()

yAxis = h.GetYaxis()
yAxis.SetTitleOffset(0)
yAxis.SetTitle("Entries")
#yAxis.CenterTitle()

event.Project("h", "muon_pt", "has_prop && n_ME11_segment == 1")

h.SetLineWidth(3)
h.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextAngle(0)
latex.SetTextColor(ROOT.kBlack)

latex.SetTextFont(42)
#latex.SetTextAlign(31)
latex.SetTextSize(0.4*canvas.GetTopMargin())
#latex.DrawLatex(1-1.2*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = h.GetEntries()))
#latex.DrawLatex(1-1.2*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.6*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h.GetMean(),3)))
#latex.DrawLatex(1-1.2*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.9*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h.GetStdDev(),3)))

latex.DrawLatex(0+4.4*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h.GetEntries())))
latex.DrawLatex(0+4.4*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h.GetMean(),3)))
latex.DrawLatex(0+4.4*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h.GetStdDev(),3)))


latex.SetTextSize(0.6*canvas.GetTopMargin())
latex.SetTextFont(61)
latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
latex.SetTextFont(52)
latex.SetTextSize(0.5*canvas.GetTopMargin())
#latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Work in Progress")

latex.SetTextFont(42)
latex.SetTextSize(0.5*canvas.GetTopMargin())
#latex.DrawLatex(4*canvas.GetLeftMargin(), 2*canvas.GetBottomMargin(), "{} {} Shift".format(sys.argv[3], sys.argv[2]))

frame = canvas.GetFrame()
frame.Draw()

if os.path.exists("plots/") == False:
  os.mkdir("plots/")
canvas.SaveAs("muonpt_1seg.png")





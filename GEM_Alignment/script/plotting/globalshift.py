import ROOT, os, sys
from make_plot import *

#fname = "out_RdPhiAna_Craft2022_TerukiRuns_ME11Iter0_GEMIter0_080522.root"
fname = "{}".format(sys.argv[1])
trees = ["ME11ana/Inner_Prop", "analyzer/ME11Seg_Prop"]
RdPhi_Low = -1
RdPhi_High = 1
cut = "abs(RdPhi) < 5 && muon_pt > 3"

f = ROOT.TFile(fname)
for tree in trees:
  event = f.Get(tree)
  branch = "RdPhi"
  if "ME11ana" in tree: branch = branch+":prop_location[3]"
  if "analyzer" in tree: branch = branch + "_Corrected:prop_location[2]"
  plotname = ""
  if "ME11ana" in tree: plotname = "ME11_CRAFT2022"
  if "analyzer" in tree: plotname = "GE11_CRAFT2022"

  for reg in [-1, 1]:
    h1 = plot1Dprofile(branch, "{plotname} RdPhi Profile R{reg}".format(plotname = plotname, reg = reg), event, [36, 1, 37, 10, -100, 100], "Chamber", "RdPhi", cut + " && prop_location[0] == {reg}".format(reg = reg), "RdPhi", False)
    c1 = ROOT.TCanvas("", "", 800, 600)
    hist = ROOT.TH1D("{plotname} RdPhi Profile R{reg}".format(plotname = plotname, reg = reg), "{plotname} RdPhi Profile R{reg}".format(plotname = plotname, reg = reg), 36, -(3.14159265)/36., 2*3.14159265 - (3.14159265)/36.)
    for i in range(1, 37):
      hist.SetBinContent(i, h1.GetBinContent(i))

    hist.GetYaxis().SetRangeUser(RdPhi_Low, RdPhi_High)
    hist.GetXaxis().SetTitle("Phi")
    hist.GetYaxis().SetTitle("RdPhi [cm]")
    hist.Sumw2()
    hist.Draw("Hist")
    hist.SetStats(False)
    c1.SaveAs("{plotname}_R{reg}_GlobalShift.png".format(plotname = plotname, reg = reg))
    del h1
    del hist

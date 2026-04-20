import ROOT, tdrstyle, os

# Detector and object specific variables
endcap = 1
layer = 1
supch_cut = "all" # 0 for even, 1 for odd, "all" for all
low_pt = 20
high_pt = 200
year = "2026B"

# comparing between 2025 GEM alignment (f) and 2026 GEM alignment (f1) using 2026B data
f = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2026/Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_2025alignment_trackerprop.root")
f1 = ROOT.TFile("/eos/user/t/toakhter/tamu_mual/2026/Run2026B_muon0_ZMu_160X_dataRun3_Prompt_frozen260223_v0_trackerprop_aligned.root")

event = f.Get("analyzer/Inner_Prop")
event1 = f1.Get("analyzer/Inner_Prop")
ROOT.gROOT.SetBatch(1)
tdrstyle.setTDRStyle()

class Variables:
    def __init__(self, endcap, layer, supch_cut, low_pt, high_pt):
        self.endcap = endcap
        self.layer = layer
        self.supch_cut = supch_cut
        self.low_pt = low_pt
        self.high_pt = high_pt
    
    def plotting_cut(self):
        cut = f"prop_location[0] == {self.endcap} && prop_location[3] == {self.layer} && muon_pt > {self.low_pt} && muon_pt < {self.high_pt}"
        if self.supch_cut == 0:
            cut += " && prop_location[2] % 2 == 0"
        elif self.supch_cut == 1:
            cut += " && prop_location[2] % 2 == 1"
        return cut
    
    def plot_axis(self):
        if self.endcap == 1:
            x_low = 0
            x_high = 37
        if self.endcap == -1:
            x_low = -37
            x_high = 0
        
        if self.supch_cut == "all":
            x_bin = (x_high - x_low)
        if self.supch_cut in [0, 1]:
            x_bin = (x_high - x_low) // 2 # 18
        return x_low, x_high, x_bin
    
    def plot_legends(self):
        endcap_str = "+Endcap" if self.endcap ==1 else "-Endcap"
        layer_str = f"Layer {self.layer}"
        pt_str = f"{self.low_pt} < pT < {self.high_pt}"
        if self.supch_cut == 0:
            supch_str = "Even chambers"
        elif self.supch_cut == 1:
            supch_str = "Odd chambers"
        else:
            supch_str = "All chambers"
        return endcap_str, layer_str, pt_str, supch_str

    def plot_save_name(self):
        name = f"R{self.endcap}L{self.layer}"
        if self.supch_cut == 0:
            name += "_even"
        elif self.supch_cut == 1:
            name += "_odd"
        name += f"_pt{self.low_pt}to{self.high_pt}"
        return name

plotting_variables = Variables(endcap, layer, supch_cut, low_pt, high_pt)
# print("Plotting cuts:", plotting_variables.plotting_cut())
x_low, x_high, x_bin = plotting_variables.plot_axis()
# print("X axis range: ", x_low, "to", x_high, "with", x_bin, "bins")
endcap_str, layer_str, pt_str, supch_str = plotting_variables.plot_legends()
# print("Legend info: ", endcap_str, layer_str, pt_str, supch_str)
# print("Plot save name: ", plotting_variables.plot_save_name())

chamber_stats = {}

for i in range(1, 37):
    h = ROOT.TH1D(f"h_{i}", f"h_{i}", 100, -2, 2)
    h1 = ROOT.TH1D(f"h1_{i}", f"h1_{i}", 100, -2, 2)

    event.Project(f"h_{i}", "RdPhi_Corrected", plotting_variables.plotting_cut() + f"&& has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[2] == {i}")
    event1.Project(f"h1_{i}", "RdPhi_Corrected", plotting_variables.plotting_cut() + f"&& has_fidcut && abs(RdPhi_Corrected) < 2 && prop_location[2] == {i}")

    mean = h.GetMean()
    mean1 = h1.GetMean()
    stdev = h.GetStdDev()
    stdev1 = h1.GetStdDev()

    chamber_stats[i] = {"mean_2025Geom": mean, "stdev_2025Geom": stdev, "mean_2026Geom": mean1, "stdev_2026Geom": stdev1}
    h.Delete()
    h1.Delete()

# print("Chamber stats: ", chamber_stats)

H_ref = 800
W_ref = 1200
W = W_ref
H = H_ref

T = 0.12*H_ref
B = 0.16*H_ref
L = 0.16*W_ref
R = 0.08*W_ref

canvas = ROOT.TCanvas("c1", "c1", 1200, 800, W, H)
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


plot = ROOT.TGraphErrors()
plot.SetTitle(f"Comparison of #DeltaR#phi between 2025 and 2026 GEM alignment;Chamber GE{endcap}/1;#DeltaR#phi [cm]")
plot.SetMarkerStyle(20)
plot.SetMarkerSize(1)
plot.SetLineColor(ROOT.kBlue)
plot.SetLineWidth(2)
plot.SetMarkerColor(ROOT.kBlue)
for i in range(1, 37):
    plot.SetPoint(i-1, i, chamber_stats[i]["mean_2025Geom"])
    plot.SetPointError(i-1, 0, chamber_stats[i]["stdev_2025Geom"]) 
plot.Draw("AP")
plot.GetYaxis().SetRangeUser(-1, 1)
plot1 = ROOT.TGraphErrors()
plot1.SetMarkerStyle(21)
plot1.SetMarkerSize(1)
plot1.SetLineColor(ROOT.kRed)
plot1.SetLineWidth(2)
plot1.SetMarkerColor(ROOT.kRed)

for i in range(1, 37):
    plot1.SetPoint(i-1, i, chamber_stats[i]["mean_2026Geom"])
    plot1.SetPointError(i-1, 0, chamber_stats[i]["stdev_2026Geom"])
plot1.Draw("P SAME")
plot1.GetYaxis().SetRangeUser(-1, 1)

legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.85)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.AddEntry(plot, "2025 GEM alignment", "P")
legend.AddEntry(plot1, "2026 GEM alignment", "P")
legend.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.DrawLatex(0.7, 0.35, endcap_str)
latex.DrawLatex(0.7, 0.3, layer_str)
latex.DrawLatex(0.7, 0.25, pt_str)
latex.DrawLatex(0.7, 0.2, supch_str)

latex.SetTextFont(42)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.SetTextAlign(32)
latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "(13.6 TeV)")

latex.SetTextAlign(12)
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-1.0*canvas.GetTopMargin(), "Run {run}".format(run=year))

latex.SetTextSize(0.5*canvas.GetTopMargin())
latex.SetTextFont(61)
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.27*canvas.GetTopMargin(), "CMS")
latex.SetTextFont(52)
latex.SetTextSize(0.3*canvas.GetTopMargin())
latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Preliminary")



if not os.path.exists("compare_stats"):
    os.makedirs("compare_stats")

canvas.SaveAs(f"compare_stats/{plotting_variables.plot_save_name()}.png")
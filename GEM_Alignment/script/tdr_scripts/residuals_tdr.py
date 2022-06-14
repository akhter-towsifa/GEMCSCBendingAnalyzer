import ROOT, tdrstyle, sys, os, array

# ARGUMENTS ************ file, shift (-.1cm), direction (X), Zlow, Zhigh

f = ROOT.TFile("{}".format(sys.argv[1]))

cut = "has_prop && abs(RdPhi) < 5 && muon_pt > 20"
cut_string = "MC Muons, p_{T} > 20 GeV"
top_string = "Simulation Preliminary"
extra_string = "Initial"
ME11dict = {
  "tree": "ME11ana/Inner_Prop",
  "RdPhi_Branch": "RdPhi",
  "cut": cut,
  "name": "ME"
}
GE11dict = {
  "tree": "analyser/ME11Seg_Prop",
  "RdPhi_Branch": "RdPhi_Corrected",
  "cut": cut + " && n_ME11_segment == 1 && has_fidcut",
  "name": "GE"
}

for dict in [GE11dict]:
  for dims in [1, 2, "globalphi"]:
    bins = []
    boarder = []
    H_ref = 800
    W_ref = 800
    if dims == 1:
      boarder = [0.12, 0.16, 0.16, 0.08]
      bins = [100, -2, 2, 0, 0, 0, 0, 0]
    if dims == 2:
      boarder = [0.12, 0.16, 0.16, 0.20]
      bins = [100, -300, 300, 100, -300, 300, -1, 1]
      H_ref = 800
      W_ref = 900
    if dims == "globalphi":
      boarder = [0.12, 0.16, 0.16, 0.08]
      bins = [0, 0, 0, 0, 0, 0, -0.5, 0.5]
      H_ref = 800
      W_ref = 1200



    event = f.Get(dict["tree"])
    ROOT.gROOT.SetBatch(1)
    tdrstyle.setTDRStyle()

    #H_ref = 800
    #W_ref = 800
    W = W_ref
    H = H_ref

    T = boarder[0]*H_ref
    B = boarder[1]*H_ref
    L = boarder[2]*W_ref
    R = boarder[3]*W_ref

    xbins = bins[0]
    xlow  = bins[1]
    xhigh = bins[2]
    ybins = bins[3]
    ylow  = bins[4]
    yhigh = bins[5]
    zlow  = bins[6]
    zhigh = bins[7]

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
    canvas.SetGridx()
    canvas.SetGridy()

    if dims == 1:
      h = ROOT.TH1D("h", "h", xbins, xlow, xhigh)

      xAxis = h.GetXaxis()
      xAxis.SetTitleOffset(2)
      xAxis.SetTitleSize(0.04)
      xAxis.SetTitle("#DeltaR#phi [cm]")
      #xAxis.CenterTitle()

      yAxis = h.GetYaxis()
      yAxis.SetTitleOffset(0)
      yAxis.SetTitleSize(0.04)
      yAxis.SetTitle("Entries")
      #yAxis.CenterTitle()


      for reg in [-1, 0, 1]:
        regname = str(reg)
        if reg == 1:
          regname = "+"+regname
        print("Regname = ", regname)

        if reg == 0:
          event.Project("h", dict["RdPhi_Branch"], dict["cut"])
        else:
          event.Project("h", dict["RdPhi_Branch"], dict["cut"] + " && prop_location[0] == {reg}".format(reg = reg))

        h.SetLineWidth(3)
        h.Draw()

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.4*canvas.GetTopMargin())

        latex.SetTextAlign(32)
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h.GetEntries())))
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), "Mean: {mean}".format(mean = round(h.GetMean(),3)))
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-1.1*canvas.GetTopMargin(), "Std Dev: {stddev}".format(stddev = round(h.GetStdDev(),3)))
        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Region: {reg}".format(reg = regname))
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), dict["name"])


        latex.SetTextSize(0.5*canvas.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.4*canvas.GetTopMargin())
        #latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Preliminary")
        latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "Work in Progress")

        latex.SetTextFont(42)
        latex.SetTextSize(0.5*canvas.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()

        if os.path.exists("plots/") == False:
          os.mkdir("plots/")
        canvas.SaveAs("plots/{name}{reg}1_res{dims}D.png".format(reg = reg, name = dict["name"], dims = dims))

    elif dims == 2:
      h = ROOT.TProfile2D("h", "h", xbins, xlow, xhigh, ybins, ylow, yhigh)

      ablue = array.array("d", [1,1,0])
      ared = array.array("d", [0,1,1])
      agreen = array.array("d", [0,1,0])
      astop = array.array("d", [0,.5,1])
      myPalette = []
      fi = ROOT.TColor.CreateGradientColorTable(3, astop, ared, agreen, ablue, 100)
      for x in range(100):
        myPalette.append(fi+x)
      ROOT.gStyle.SetPalette(100, array.array("i", myPalette))

      xAxis = h.GetXaxis()
      xAxis.SetTitleOffset(2)
      xAxis.SetTitleSize(0.04)
      xAxis.SetTitle("Global X [cm]")
      #xAxis.CenterTitle()

      yAxis = h.GetYaxis()
      yAxis.SetTitleOffset(0)
      yAxis.SetTitleSize(0.04)
      yAxis.SetTitle("Global Y [cm]")
      #yAxis.CenterTitle()

      zAxis = h.GetZaxis()
      zAxis.SetTitle("#DeltaR#phi [cm]")
      #zAxis.CenterTitle()
      zAxis.SetTitleOffset(2)
      zAxis.SetTitleSize(0.04)
      zAxis.SetRangeUser(zlow, zhigh)

      for reg in [-1, 0, 1]:
        if reg == 0:
          event.Project("h", dict["RdPhi_Branch"]+":prop_GP[1]:prop_GP[0]", dict["cut"])
        else:
          event.Project("h", dict["RdPhi_Branch"]+":prop_GP[1]:prop_GP[0]", dict["cut"] + " && prop_location[0] == {reg}".format(reg = reg))

        #h.SetLineWidth(3)
        #h.Draw()
        h.Draw("colz")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.4*canvas.GetTopMargin())

        latex.SetTextAlign(32)
        latex.DrawLatex(1-1.1*canvas.GetRightMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "Entries: {entries}".format(entries = int(h.GetEntries())))
        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "{name}{R}/1".format(R = regname, name = dict["name"]))
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), extra_string)
        latex.DrawLatex(1.1*canvas.GetLeftMargin(), 1.2*canvas.GetBottomMargin(), cut_string)



        latex.SetTextSize(0.5*canvas.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.4*canvas.GetTopMargin())
        latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), top_string)

        latex.SetTextFont(42)
        latex.SetTextSize(0.5*canvas.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()

        if os.path.exists("plots/") == False:
          os.mkdir("plots/")
        canvas.SaveAs("plots/{name}{reg}1_res{dims}D.pdf".format(reg = reg, name = dict["name"], dims = dims))

    elif dims == "globalphi":
      """
      W = W_ref
      H = H_ref

      T = boarder[0]*H
      B = boarder[1]*H
      L = boarder[2]*W
      R = boarder[3]*W

      canvas = ROOT.TCanvas("c1", "c1", W, H)
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
      canvas.SetGridx()
      canvas.SetGridy()
      """

      for i, region in enumerate([-1, 1]):
        regname = str(region)
        if region == 1:
          regname = "+"+regname
        print("Regname = ", regname)

        hlist_tmp = []
        hist = ROOT.TH1D("GEM RdPhi Profile R{reg}".format(reg = regname), "GEM RdPhi Profile R{reg}".format(reg = regname), 36, -(3.14159265)/36., 2*3.14159265 - (3.14159265)/36.)
        for ch in range(1, 37):
          hlist_tmp.append(ROOT.TH1D("ch{ch}".format(ch = ch), "ch{ch}".format(ch = ch), 300, -5, 5))
          event.Project("ch{ch}".format(ch = ch), "RdPhi_Corrected", cut + " && prop_location[0] == {reg} && prop_location[2] == {ch}".format(reg = region, ch = ch))
          hist.SetBinContent(ch, hlist_tmp[ch-1].GetMean())
          hist.SetBinError(ch, hlist_tmp[ch-1].GetMeanError())

        hist.GetYaxis().SetRangeUser(zlow, zhigh)
        hist.GetXaxis().SetTitle("Phi")
        hist.GetYaxis().SetTitle("RdPhi [cm]")
        #hlist[i].Sumw2()
        hist.SetMarkerStyle(22)
        hist.SetMarkerSize(2)
        hist.Draw("E")
        hist.SetStats(False)



        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAngle(0)
        latex.SetTextColor(ROOT.kBlack)

        latex.SetTextFont(42)
        latex.SetTextSize(0.4*canvas.GetTopMargin())

        latex.SetTextAlign(12)
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.3*canvas.GetTopMargin(), "{name}{R}/1".format(R = regname, name = dict["name"]))
        latex.DrawLatex(0+1.1*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()-0.7*canvas.GetTopMargin(), extra_string)
        latex.DrawLatex(1.1*canvas.GetLeftMargin(), 1.2*canvas.GetBottomMargin(), cut_string)

        latex.SetTextSize(0.5*canvas.GetTopMargin())
        latex.SetTextFont(61)
        latex.DrawLatex(canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), "CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.4*canvas.GetTopMargin())
        latex.DrawLatex(1.9*canvas.GetLeftMargin(), 1-canvas.GetTopMargin()+0.2*canvas.GetTopMargin(), top_string)

        latex.SetTextFont(42)
        latex.SetTextSize(0.5*canvas.GetTopMargin())

        frame = canvas.GetFrame()
        frame.Draw()



        canvas.SaveAs("plots/{name}{reg}1_GlobalShift.pdf".format(name = dict["name"], reg = region))

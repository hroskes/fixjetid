#!/usr/bin/env python

import errno, math, os, subprocess
import ROOT

oldmaindir = "/work-zfs/lhc/CJLSTtrees/190821/"
newmaindir = "/work-zfs/lhc/CJLSTtrees/190821_fixjetid/"

eventsperjob = 10000

def mkdir_p(path):
  """http://stackoverflow.com/a/600612/5228524"""
  try:
    os.makedirs(path)
  except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise


for a in "Data_2016", "Data_2017", "Data_2018", "MC_2016", "MC_2016_anomalous", "MC_2017", "MC_2017_anomalous", "MC_2018", "MC_2018_anomalous":
  for b in os.listdir(os.path.join(oldmaindir, a)):
    if b == "AAAOK": continue
    if b.lower() == "chunks": continue
    if b == "Chunks_2018_MC_b10d7cd8": continue
    if "Data" in a and b != "AllData": continue
    olddir = os.path.join(oldmaindir, a, b)
    oldfilename = os.path.join(olddir, "ZZ4lAnalysis.root")
    newdir = os.path.join(newmaindir, a, b)
    mkdir_p(newdir)
    newfilename = os.path.join(newdir, "ZZ4lAnalysis.root")
    if os.path.exists(newfilename): continue

    f = ROOT.TFile(oldfilename)
    t = f.Get("ZZTree/candTree")
    events = t.GetEntries()
    if "AllData" in b:
      t = f.Get("CRZLLTree/candTree")
      events = max(events, t.GetEntries())
    njobs = int(math.ceil(1.0 * events / eventsperjob))

    dohadd = True
    haddcommand = ["hadd", newfilename]

    for i in range(njobs):
      firstevent = i*eventsperjob
      lastevent = (i+1)*eventsperjob - 1
      newsubfilename = os.path.join(newdir, "ZZ4lAnalysis_{}_{}.root".format(firstevent, lastevent))
      cmdline = [
        "sbatch",
        "--job-name="+os.path.basename(newsubfilename),
        "--time=1:0:0",
        "--nodes=1",
        "--mem=3000",
        "--partition=shared",
        "--output="+os.path.join(newdir, newsubfilename.replace(".root", ".out")),
        "--error="+os.path.join(newdir, newsubfilename.replace(".root", ".err")),

        "run.sh", oldfilename, newsubfilename,
        "--first-event", str(firstevent), "--last-event", str(lastevent),
      ]
      if "AllData" in b:
        cmdline.append("--doCRZLL")

      haddcommand.append(newsubfilename)

      if not os.path.exists(newsubfilename):
        subprocess.check_call(cmdline)
        dohadd = False

    if dohadd:
      subprocess.check_call(haddcommand)

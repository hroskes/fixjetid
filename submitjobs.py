#!/usr/bin/env python

import argparse, errno, getpass, math, os, shutil, tempfile, subprocess

p = argparse.ArgumentParser()
p.add_argument("--testing", action="store_true")
p.add_argument("--print-order", action="store_true")
p.add_argument("--remove-JER", action="store_true")
args = p.parse_args()

import ROOT

oldmaindir = "/work-zfs/lhc/CJLSTtrees/190821/"
if args.remove_JER:
  newmaindir = "/work-zfs/lhc/CJLSTtrees/190821_fixjetid_removeJER/"
else:
  newmaindir = "/work-zfs/lhc/CJLSTtrees/190821_fixjetid/"

eventsperjob = 10000

def jobisrunning(name):
  output = subprocess.check_output(["squeue", "--name", name])
  for line in output.strip().split("\n"):
    if "JOBID" not in line: return True
  return False

def mkdir_p(path):
  """http://stackoverflow.com/a/600612/5228524"""
  try:
    os.makedirs(path)
  except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir(path):
      pass
    else:
      raise


def run(a, b, testing=False, print_order=False, remove_JER=False):
  olddir = os.path.join(oldmaindir, a, b)
  oldfilename = os.path.join(olddir, "ZZ4lAnalysis.root")
  newdir = os.path.join(newmaindir, a, b)
  mkdir_p(newdir)
  newfilename = os.path.join(newdir, "ZZ4lAnalysis.root")
  if os.path.exists(newfilename):
    f = ROOT.TFile(newfilename)
    assert f.Get("ZZTree/candTree").GetEntries()
    if b == "AllData": assert f.Get("CRZLLTree/candTree").GetEntries()
    return

  f = ROOT.TFile(oldfilename)
  t = f.Get("ZZTree/candTree")
  events = t.GetEntries()
  if "AllData" in b:
    t = f.Get("CRZLLTree/candTree")
    events = max(events, t.GetEntries())
  njobs = int(math.ceil(1.0 * events / eventsperjob))

  dohadd = True
  haddcommand = ["hadd", None]

  if print_order:
    print "{:20} {:30} {:>3}".format(a, b, njobs)
    return

  for i in range(njobs):
    firstevent = i*eventsperjob
    lastevent = (i+1)*eventsperjob - 1
    newsubfilename = os.path.join(newdir, "ZZ4lAnalysis_{}_{}.root".format(firstevent, lastevent))
    if njobs == 1:
      newsubfilename = newfilename
      dohadd = False
    cmdline = [
      "sbatch",
      "--job-name="+os.path.join(a, b, os.path.basename(newsubfilename)),
      "--time=5:0:0",
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
    if remove_JER:
      cmdline.append("--remove-JER")

    haddcommand.append(newsubfilename)

    if jobisrunning(os.path.join(a, b, os.path.basename(newsubfilename))):
      dohadd = False
      continue

    if os.path.exists(newsubfilename):
      f = ROOT.TFile(newsubfilename)
      if not f.Get("ZZTree/candTree"):
        del f
        os.remove(newsubfilename)
        dohadd = False

    if not os.path.exists(newsubfilename):
      subprocess.check_call(cmdline)
      dohadd = False
      if testing: return "submitted"

  if dohadd:
    with tempfile.NamedTemporaryFile(bufsize=0, suffix=".root") as f:
      os.remove(f.name)
      haddcommand[1] = f.name
      subprocess.check_call(haddcommand)
      shutil.move(f.name, newfilename)
      with open(f.name, "w"): pass

folders = []

for a in "Data_2016", "Data_2017", "Data_2018", "MC_2016", "MC_2016_anomalous", "MC_2017", "MC_2017_anomalous", "MC_2018", "MC_2018_anomalous":
  for b in os.listdir(os.path.join(oldmaindir, a)):
    if b == "AAAOK": continue
    if b.lower() == "chunks": continue
    if b == "Chunks_2018_MC_b10d7cd8": continue
    if "Data" in a and b != "AllData": continue
    if args.remove_JER and "MC_2018" not in a: continue

    folders.append((a, b))

def sortorder(folder):
  a, b = folder
  if "anomalous" in a or "Data" in a: return 0
  if "ZZTo4l" in b or "minlo" in b or b in ("ggH125", "VBFH125", "ZH125", "WplusH125", "WminusH125", "ttH125"): return 0, a, b
  if b.startswith("ggTo") or b == "bbH125" or "125_tune" in b: return 1, a, b
  return 2, a, b

folders.sort(key=sortorder)

for a, b in folders:
  result = run(a, b, **args.__dict__)
  if args.testing and result == "submitted": break

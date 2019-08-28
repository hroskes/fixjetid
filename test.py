#!/usr/bin/env python

import argparse

p = argparse.ArgumentParser()
p.add_argument("jec", choices="nominal up down".split())
p.add_argument("--print-jets", action="store_true")
p.add_argument("--mavjj-true", action="store_true")
p.add_argument("--pt-sort", action="store_true")
args = p.parse_args()

from itertools import izip

def maketlv(pt, eta, phi, m):
  result = ROOT.TLorentzVector()
  result.SetPtEtaPhiM(pt, eta, phi, m)
  return result

jetpt = [62.22626495361328, 61.98612976074219, 61.104454040527344, 49.428287506103516, 40.636104583740234, 28.233963012695312, 25.520647048950195]
jeteta = [3.114255428314209, 2.0353877544403076, -0.7179408669471741, 2.6708858013153076, 0.6869650483131409, 2.6115543842315674, -0.23656509816646576]
jetphi = [-0.8394545316696167, 0.04743506759405136, -1.2854363918304443, 3.066879987716675, -1.7294487953186035, -0.622397780418396, -1.4944148063659668]
jetmass = [10.948182106018066, 8.58421802520752, 10.160293579101562, 12.757057189941406, 6.870137691497803, 6.41053581237793, 5.189918518066406]
jetsigma = [0.08672846853733063, 0.033871352672576904, 0.014012597501277924, 0.17160306870937347, 0.017611980438232422, 0.1811283975839615, 0.025665774941444397]

if args.jec == "up":
  jetpt = [pt * (1+sigma) for pt, sigma in izip(jetpt, jetsigma)]
  jetmass = [mass * (1+sigma) for mass, sigma in izip(jetmass, jetsigma)]
elif args.jec == "down":
  jetpt = [pt * (1-sigma) for pt, sigma in izip(jetpt, jetsigma)]
  jetmass = [mass * (1-sigma) for mass, sigma in izip(jetmass, jetsigma)]

lepid = [-11, 11, -11, 11]
leppt = [45.29861068725586, 110.45686340332031, 25.646299362182617, 21.341201782226562]
lepeta = [0.07161509990692139, 1.263501524925232, -0.24340572953224182, -0.1200866773724556]
lepphi = [2.6503918170928955, 2.617446184158325, 0.9813779592514038, 3.002392530441284]

import ROOT

jets = [(0, maketlv(*_)) for _ in izip(jetpt, jeteta, jetphi, jetmass)]
leptons = [(id, maketlv(pt, eta, phi, 0)) for id, pt, eta, phi in izip(lepid, leppt, lepeta, lepphi)]

if args.pt_sort:
  jets.sort(key=lambda x: x[1].Pt(), reverse=True)

if args.print_jets:
  for id, _ in jets:
    _.Print()
  assert False

from ZZMatrixElement.MELA.mela import Mela, TVar, SimpleParticleCollection_t

jets = SimpleParticleCollection_t(jets)
leptons = SimpleParticleCollection_t(leptons)

m = Mela()


m.setInputEvent(leptons, jets, None, False)
m.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.Had_ZH)
print m.computeDijetConvBW(args.mavjj_true)

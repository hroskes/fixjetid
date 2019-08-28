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

jetpt = [37.2165641784668, 33.699684143066406, 30.714012145996094]
jeteta = [-1.0453099012374878, -2.7907354831695557, 0.07553417980670929]
jetphi = [2.8454203605651855, -0.6791841983795166, -1.0391074419021606]
jetmass = [7.891648769378662, 8.078144073486328, 6.681240081787109]
jetsigma = [0.019978823140263557, 0.1735641360282898, 0.02200276032090187]

if args.jec == "up":
  jetpt = [pt * (1+sigma) for pt, sigma in izip(jetpt, jetsigma)]
  jetmass = [mass * (1+sigma) for mass, sigma in izip(jetmass, jetsigma)]
elif args.jec == "down":
  jetpt = [pt * (1-sigma) for pt, sigma in izip(jetpt, jetsigma)]
  jetmass = [mass * (1-sigma) for mass, sigma in izip(jetmass, jetsigma)]

lepid = [-11, 11, -13, 13]
leppt = [56.49150466918945, 28.424598693847656, 14.310689926147461, 60.03144454956055]
lepeta = [-1.9336107969284058, -0.8147843480110168, -2.110483169555664, -1.4010090827941895]
lepphi = [2.3444762229919434, -1.4932364225387573, 1.7744771242141724, -0.11107052117586136]

import ROOT

jets = [(0, maketlv(*_)) for _ in izip(jetpt, jeteta, jetphi, jetmass)]
leptons = [(id, maketlv(pt, eta, phi, 0)) for id, pt, eta, phi in izip(lepid, leppt, lepeta, lepphi)]

if args.pt_sort:
  jets.sort(key=lambda x: x[1].Pt(), reverse=True)
jets = [j for j in jets if j[1].Pt() > 30]

if args.print_jets:
  for id, _ in jets:
    print _.Pt(); _.Print()
  assert False

from ZZMatrixElement.MELA.mela import Mela, TVar, SimpleParticleCollection_t

jets = SimpleParticleCollection_t(jets)
leptons = SimpleParticleCollection_t(leptons)

m = Mela()


m.setInputEvent(leptons, jets, None, False)
m.setProcess(TVar.HSMHiggs, TVar.JHUGen, TVar.Had_ZH)
print m.computeDijetConvBW(args.mavjj_true)

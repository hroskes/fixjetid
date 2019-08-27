#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  p.add_argument("infile")
  p.add_argument("outfile")
  p.add_argument("--no-id", action="store_false", dest="applyid")
  p.add_argument("--no-pu-id", action="store_false", dest="applypuid")
  args = p.parse_args()

import abc, os
from itertools import izip

import numpy as np
import ROOT

class TFile(object):
  def __init__(self, filename, *args, **kwargs):
    self.__filename = filename
    self.__args = args
    self.__deleteifbad = kwargs.pop("deleteifbad", False)
    self.__entered = False
  def __enter__(self):
    import ROOT
    self.__bkpdirectory = ROOT.gDirectory.GetDirectory(ROOT.gDirectory.GetPath())
    self.__f = ROOT.TFile.Open(self.__filename, *self.__args)
    self.__entered = True
    if not self.__f:
      raise IOError(self.__filename+" is a null pointer, see above for details.")
    if self.__f.IsZombie():
      self.__exit__()
      raise IOError(self.__filename+" is a zombie, see above for details.")

    try:
      openoption = self.__args[0].upper()
    except IndexError:
      openoption = ""

    self.__write = {
      "": False,
      "READ": False,
      "NEW": True,
      "CREATE": True,
      "RECREATE": True,
      "UPDATE": True,
    }[openoption]

    return self.__f

  def __exit__(self, *errorstuff):
    try:
      if self.__write and (not any(errorstuff) or not self.__deleteifbad):
        self.__f.Write()
      self.__f.Close()
      self.__bkpdirectory.cd()
    except:
      if self.__write and self.__deleteifbad:
        os.remove(self.__filename)
      raise

    if self.__write and self.__deleteifbad and any(errorstuff):
      os.remove(self.__filename)

class Branch(object):
  __metaclass__ = abc.ABCMeta
  def __init__(self, name, typ):
    self.__name = name
    self.__type = typ

  @abc.abstractproperty
  def thingforsetbranchaddress(self): pass

  def attachtotree(self, t):
    t.SetBranchAddress(self.__name, self.thingforsetbranchaddress)
    target = type(getattr(t, self.__name))
    if target != self.__type:
      if target == int and self.__type == np.short:
        pass #ok
      else:
        raise ValueError("Wrong type for {}: {}, should be {}".format(self.__name, self.__type, type(getattr(t, self.__name))))

  def convertforxcheck(self, thing): return thing

  @abc.abstractmethod
  def setthingfortree(self, t, applyid, applypuid): pass

  def setbranchvalue(self, t, applyid, applypuid, doxcheck):
    newvalue = self.setthingfortree(t, applyid, applypuid)
    if doxcheck:
      old = self.convertforxcheck(getattr(t, self.__name))
      new = self.convertforxcheck(newvalue)
      if new != old:
        raise ValueError("old value of {} is {}, not the same as the new calculated value {}".format(self.__name, old, new))

  @property
  def name(self): return self.__name


class NormalBranch(Branch):
  def __init__(self, name, function, typ):
    super(NormalBranch, self).__init__(name, typ)
    self.__function = function
    self.__array = np.array([0], typ)
  @property
  def thingforsetbranchaddress(self): return self.__array
  def setthingfortree(self, t, applyid, applypuid):
    self.__array[0] = self.__function(t, applyid, applypuid)
    return self.__array[0]


class JetVectorBranch(Branch):
  def __init__(self, name, typ):
    vectortyp = ROOT.vector(typ)
    super(JetVectorBranch, self).__init__(name, vectortyp)
    self.__vector = vectortyp(0)
  @property
  def thingforsetbranchaddress(self): return self.__vector
  def setthingfortree(self, t, applyid, applypuid):
    self.__vector.clear()
    for value, id, puid in izip(getattr(t, self.name), t.JetID, t.JetPUID):
      if (id or not applyid) and (puid or not applypuid):
        self.__vector.push_back(value)
    return self.__vector

  def convertforxcheck(self, thing): return list(thing)

class NJetsBranch(NormalBranch):
  def __init__(self, name, jetcondition):
    def njetsfunction(t, applyid, applypuid):
      return sum(1 for id, puid, pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn in izip(t.JetID, t.JetPUID, t.JetPt, t.JetSigma, t.JetIsBtagged, t.JetIsBtaggedWithSF, t.JetIsBtaggedWithSFUp, t.JetIsBtaggedWithSFDn) if (id or not applyid) and (puid or not applypuid) and jetcondition(pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn))
    super(NJetsBranch, self).__init__(name, njetsfunction, np.short)

branches = [
  NJetsBranch("nCleanedJets", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: True),
  NJetsBranch("nCleanedJetsPt30", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30),
  NJetsBranch("nCleanedJetsPt30_jecUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1+sigma)>30),
  NJetsBranch("nCleanedJetsPt30_jecDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1-sigma)>30),
  NJetsBranch("nCleanedJetsPt30BTagged", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtagged),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsf),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSFUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsfup),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSFDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsfdn),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF_jecUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1+sigma)>30 and isbtagged),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF_jecDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1-sigma)>30 and isbtagged),

  JetVectorBranch("JetPt", float),
  JetVectorBranch("JetEta", float),
  JetVectorBranch("JetPhi", float),
  JetVectorBranch("JetMass", float),
  JetVectorBranch("JetBTagger", float),
  JetVectorBranch("JetIsBtagged", float),
  JetVectorBranch("JetIsBtaggedWithSF", float),
  JetVectorBranch("JetIsBtaggedWithSFUp", float),
  JetVectorBranch("JetIsBtaggedWithSFDn", float),
  JetVectorBranch("JetQGLikelihood", float),
  JetVectorBranch("JetAxis2", float),
  JetVectorBranch("JetMult", float),
  JetVectorBranch("JetPtD", float),
  JetVectorBranch("JetSigma", float),
  JetVectorBranch("JetHadronFlavour", "short"),
  JetVectorBranch("JetPartonFlavour", "short"),
  JetVectorBranch("JetRawPt", float),
  JetVectorBranch("JetPtJEC_noJER", float),
  JetVectorBranch("JetJERUp", float),
  JetVectorBranch("JetJERDown", float),
  JetVectorBranch("JetID", "short"),
  JetVectorBranch("JetPUID", "short"),
  JetVectorBranch("JetPUValue", float),
]

def fixjetid(infile, outfile, applyid=True, applypuid=True, folders=["ZZTree"]):
  print "Processing", infile, "-->", outfile
  with TFile(infile) as f, TFile(outfile, "CREATE", deleteifbad=True) as newf:
    for foldername in folders:
      folder = f.Get(foldername)
      newfolder = newf.mkdir(foldername)
      Counters = folder.Get("Counters")
      Counters.SetDirectory(newfolder)

      failedt = folder.Get("candTree_failed")
      if failedt:
        failedt.SetDirectory(newfolder)

      t = folder.Get("candTree")
      newt = t.CloneTree(0, "fast")
      newt.SetDirectory(newfolder)

      for branch in branches:
        branch.attachtotree(newt)

      nentries = t.GetEntries()
      nxchecks = 0
      for i, entry in enumerate(t, start=1):
        doxcheck = (all(t.JetID) or not applyid) and (all(t.JetPUID) or not applypuid)
        nxchecks += doxcheck
        for branch in branches:
          branch.setbranchvalue(t, applyid, applypuid, doxcheck)
        newt.Fill()

        if i%10000 == 0 or i == nentries:
          print i, "/", nentries

    print "Done!"
    print "Did xchecks on", nxchecks, "that had all jets passing ID and PUID: all confirmed unchanged"
    print

if __name__ == "__main__":
  fixjetid(**args.__dict__)

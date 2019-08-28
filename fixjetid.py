#!/usr/bin/env python

if __name__ == "__main__":
  import argparse
  p = argparse.ArgumentParser()
  p.add_argument("infile")
  p.add_argument("outfile")
  p.add_argument("--no-id", action="store_false", dest="applyid")
  p.add_argument("--no-pu-id", action="store_false", dest="applypuid")
  p.add_argument("--test-no-mela", action="store_true")
  p.add_argument("--first-event", type=int, default=0, dest="firstevent")
  p.add_argument("--last-event", type=int, default=0, dest="lastevent")
  p.add_argument("--debug", action="store_true")
  p.add_argument("--doCRZLL", action="store_true")
  args = p.parse_args()


class MelaWrapper(object):
  def __init__(self, *melaargs, **kwargs):
    self.__melaargs = melaargs
    self.__mela = kwargs.pop("mela", None)
    self.__melaeventdummybranch = None
    assert not kwargs
  @property
  def mela(self):
    if self.__mela is None:
      assert False
      from ZZMatrixElement.MELA.mela import Mela
      self.__mela = Mela(*self.__melaargs)
    return self.__mela
  def __getattr__(self, attr): return getattr(self.mela, attr)
  def __setattr__(self, attr, value):
    if "__" in attr: return super(MelaWrapper, self).__setattr__(attr, value)
    return setattr(self.mela, attr, value)

  def setInputEvent(self, melaeventdummybranch):
    if self.__melaeventdummybranch is melaeventdummybranch: return
    self.resetInputEvent()
    self.__melaeventdummybranch = melaeventdummybranch
    self.mela.setInputEvent(*melaeventdummybranch.lastsetbranchvalue)

  def resetInputEvent(self):
    self.__melaeventdummybranch = None
    if self.__mela is not None:
      self.__mela.resetInputEvent()

if args.test_no_mela:
  mela = MelaWrapper()
else:
  from ZZMatrixElement.MELA.mela import Mela
  mela = MelaWrapper(mela=Mela())

import abc, collections, os, re
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

class FailedXcheckError(Exception):
  def __init__(self, branch, old, new):
    self.branch = branch
    self.old = old
    self.new = new
    super(FailedXcheckError, self).__init__("old value of {} is {}, not the same as the new calculated value {}".format(branch, old, new))

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
      elif target == float and self.__type == np.float32:
        pass #ok
      else:
        raise ValueError("Wrong type for {}: {}, should be {}".format(self.__name, self.__type, type(getattr(t, self.__name))))

  def convertforxcheck(self, thing): return thing
  @abc.abstractmethod
  def compareforxcheck(self, new, old, t): pass

  @abc.abstractmethod
  def setthingfortree(self, t, applyid, applypuid): pass

  def setbranchvalue(self, t, applyid, applypuid, doxcheck):
    newvalue = self.setthingfortree(t, applyid, applypuid)
    self.lastsetbranchvalue = newvalue
    if doxcheck:
      old = self.convertforxcheck(getattr(t, self.__name))
      new = self.convertforxcheck(newvalue)
      if not self.compareforxcheck(new, old, t):
        raise FailedXcheckError(self.__name, old, new)

  @property
  def name(self): return self.__name

class DummyBranch(Branch):
  """
  Branch that doesn't actually go into the tree but can be used to calculate things
  """
  def __init__(self, name):
    return super(DummyBranch, self).__init__(name, "DUMMY")
  @property
  def thingforsetbranchaddress(self): return None
  def attachtotree(self, t): pass

  def setbranchvalue(self, t, applyid, applypuid, doxcheck):
    super(DummyBranch, self).setbranchvalue(t, applyid, applypuid, False)

  def compareforxcheck(self, new, old, t): assert False

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
  def compareforxcheck(self, new, old, t): return np.isclose(new, old, atol=0, rtol=1e-6)

class NotRecalculatedBranch(NormalBranch):
  def __init__(self, name, dummyvalue, typ):
    super(NotRecalculatedBranch, self).__init__(name, lambda *stuff: dummyvalue, typ)
  def setbranchvalue(self, t, applyid, applypuid, doxcheck):
    super(NotRecalculatedBranch, self).setbranchvalue(t, applyid, applypuid, False)

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
  def compareforxcheck(self, new, old, t): return list(new) == list(old)

class NJetsBranch(NormalBranch):
  def __init__(self, name, jetcondition):
    def njetsfunction(t, applyid, applypuid):
      return sum(1 for id, puid, pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn in izip(t.JetID, t.JetPUID, t.JetPt, t.JetSigma, t.JetIsBtagged, t.JetIsBtaggedWithSF, t.JetIsBtaggedWithSFUp, t.JetIsBtaggedWithSFDn) if (id or not applyid) and (puid or not applypuid) and jetcondition(pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn))
    super(NJetsBranch, self).__init__(name, njetsfunction, np.short)

def maketlv(pt, eta, phi, m):
  result = ROOT.TLorentzVector()
  result.SetPtEtaPhiM(pt, eta, phi, m)
  return result

class FirstNJetMomenta(DummyBranch):
  def __init__(self, n, ptbranch, etabranch, phibranch, massbranch, sigmabranch=None, jecdirection=None, rejectiftoomany=False):
    super(FirstNJetMomenta, self).__init__("first{}jetsjec{}".format(n, jecdirection))
    self.__n = n
    self.__ptbranch = ptbranch
    self.__etabranch = etabranch
    self.__phibranch = phibranch
    self.__massbranch = massbranch
    self.__rejectiftoomany = rejectiftoomany
    self.__sigmabranch = sigmabranch
    self.__jecdirection = jecdirection
  def setthingfortree(self, t, applyid, applypuid):
    pt = self.__ptbranch.lastsetbranchvalue
    mass = self.__massbranch.lastsetbranchvalue
    if self.__sigmabranch is not None:
      sigma = self.__sigmabranch.lastsetbranchvalue
      pt = [_ * (1 + self.__jecdirection * s) for _, s in izip(pt, sigma)]
      mass = [_ * (1 + self.__jecdirection * s) for _, s in izip(mass, sigma)]

    eta = self.__etabranch.lastsetbranchvalue
    phi = self.__phibranch.lastsetbranchvalue

    result = [maketlv(*_) for _ in izip(pt, eta, phi, mass)]
    result = [j for j in result if j.Pt() > 30]
#    result.sort(key=lambda x: x.Pt(), reverse=True)

    if self.__rejectiftoomany and len(result) > self.__n: result = []

    result = result[:self.__n]
    return result
  @property
  def n(self):
    return self.__n

class FirstJetsVariable(NormalBranch):
  def __init__(self, name, functiononjets, typ, firstjetmomentabranch, fallbackvalue):
    def function(t, applyid, applypuid):
      jets = firstjetmomentabranch.lastsetbranchvalue
      if len(jets) < firstjetmomentabranch.n: return fallbackvalue
      return functiononjets(*jets)
    super(FirstJetsVariable, self).__init__(name, function, typ)

class MELAEventDummyBranch(DummyBranch):
  def __init__(self, firstjetmomentabranch):
    super(MELAEventDummyBranch, self).__init__(firstjetmomentabranch.name+"MELA")
    self.__firstjetmomentabranch = firstjetmomentabranch
  def setthingfortree(self, t, applyid, applypuid):
    if args.test_no_mela: assert False
    from ZZMatrixElement.MELA.mela import SimpleParticleCollection_t
    leptons = SimpleParticleCollection_t(
      (lepid, maketlv(pt, eta, phi, 0))
        for lepid, pt, eta, phi in izip(t.LepLepId, t.LepPt, t.LepEta, t.LepPhi)
    )
    for pt, eta, phi, which in izip(t.fsrPt, t.fsrEta, t.fsrPhi, t.fsrLept):
      leptons[which-1].second += maketlv(pt, eta, phi, 0)
    jets = SimpleParticleCollection_t(
      (0, jet) for jet in self.__firstjetmomentabranch.lastsetbranchvalue
    )
    mothers = None
    isGen = False
    return leptons, jets, mothers, isGen
  @property
  def firstjetmomentabranch(self): return self.__firstjetmomentabranch

class MELABranch(FirstJetsVariable):
  def compareforxcheck(self, new, old, t):
    if (
      ("_JVBF_" in self.name or "_JQCD_" in self.name)
      and (new == -1 or "pConst" in self.name and new == 1)
      and old > 0
      and (
        "JECUp" in self.name and t.nCleanedJetsPt30_jecUp > 1
        or "JECNominal" in self.name and t.nCleanedJetsPt30 > 1
      )
    ): return True
    return np.isclose(new, old, atol=1e-15, rtol=1e-2)

class MELAProbability(MELABranch):
  def __init__(self, name, melaeventdummybranch, setprocessargnames, couplings, melafunctionname, melafunctionargs, subtractbranches=[], initialQQ=False):
    setprocessargs = []
    def function(*jets):
      mela.setInputEvent(melaeventdummybranch)
      if not setprocessargs:
        from ZZMatrixElement.MELA.mela import TVar
        for _ in setprocessargnames: setprocessargs.append(getattr(TVar, _))
      mela.setProcess(*setprocessargs)
      for coupling, value in couplings.iteritems():
        setattr(mela, coupling, value)
      result = getattr(mela, melafunctionname)(*melafunctionargs)

      if initialQQ:
        result = 0
        iorcd = mela.getIORecord()
        mearray = iorcd.getWeightedMEArray()
        MEsq = iorcd.getUnweightedMEArray()
        partonWeight = iorcd.getPartonWeights()
        for i in range(-5, 6):
          for j in range(-5, 6):
            if i==0 or j==0: continue
            result += mearray[i+5][j+5]#MEsq[i+5][j+5]*partonWeight[0][i+5]*partonWeight[1][j+5]
        cMEAvg=mela.getConstant()
        result=result*cMEAvg
      return result - sum(_.lastsetbranchvalue for _ in subtractbranches)
    super(MELAProbability, self).__init__(name, function, np.float32, melaeventdummybranch.firstjetmomentabranch, -1)
    self.firstjetmomentabranch = melaeventdummybranch.firstjetmomentabranch
    self.__subtractbranches = subtractbranches

  def compareforxcheck(self, new, old, t):
    if super(MELAProbability, self).compareforxcheck(new, old, t): return True
    if self.__subtractbranches:
      if all(abs(new / _.lastsetbranchvalue) < 1e-4 for _ in self.__subtractbranches): return True
    return False

class MELAPConst(MELABranch):
  def __init__(self, melaprobability):
    assert melaprobability.name.startswith("p_")
    name = melaprobability.name.replace("p_", "pConst_", 1)
    def function(*jets):
      return mela.getConstant()
    super(MELAPConst, self).__init__(name, function, np.float32, melaprobability.firstjetmomentabranch, 1)

class MELAPAux(MELABranch):
  def __init__(self, melaprobability):
    assert melaprobability.name.startswith("p_")
    name = melaprobability.name.replace("p_", "pAux_", 1)
    def function(*jets):
      return mela.getPAux()
    super(MELAPAux, self).__init__(name, function, np.float32, melaprobability.firstjetmomentabranch, 1)

branches = [
  NJetsBranch("nCleanedJets", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: True),
  NJetsBranch("nCleanedJetsPt30", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30),
  NJetsBranch("nCleanedJetsPt30_jecUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1+sigma)>30),
  NJetsBranch("nCleanedJetsPt30_jecDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1-sigma)>30),
  NJetsBranch("nCleanedJetsPt30BTagged", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtagged),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsf),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSFUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsfup),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSFDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt>30 and isbtaggedsfdn),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF_jecUp", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1+sigma)>30 and isbtaggedsf),
  NJetsBranch("nCleanedJetsPt30BTagged_bTagSF_jecDn", lambda pt, sigma, isbtagged, isbtaggedsf, isbtaggedsfup, isbtaggedsfdn: pt*(1-sigma)>30 and isbtaggedsf),
]

jetpt = JetVectorBranch("JetPt", "float")
jeteta = JetVectorBranch("JetEta", "float")
jetphi = JetVectorBranch("JetPhi", "float")
jetmass = JetVectorBranch("JetMass", "float")
jetsigma = JetVectorBranch("JetSigma", "float")

branches += [
  jetpt,
  jeteta,
  jetphi,
  jetmass,
  JetVectorBranch("JetBTagger", "float"),
  JetVectorBranch("JetIsBtagged", "float"),
  JetVectorBranch("JetIsBtaggedWithSF", "float"),
  JetVectorBranch("JetIsBtaggedWithSFUp", "float"),
  JetVectorBranch("JetIsBtaggedWithSFDn", "float"),
  JetVectorBranch("JetQGLikelihood", "float"),
  JetVectorBranch("JetAxis2", "float"),
  JetVectorBranch("JetMult", "float"),
  JetVectorBranch("JetPtD", "float"),
  jetsigma,
  JetVectorBranch("JetHadronFlavour", "short"),
  JetVectorBranch("JetPartonFlavour", "short"),
  JetVectorBranch("JetRawPt", "float"),
  JetVectorBranch("JetPtJEC_noJER", "float"),
  JetVectorBranch("JetJERUp", "float"),
  JetVectorBranch("JetJERDown", "float"),
  JetVectorBranch("JetID", "short"),
  JetVectorBranch("JetPUID", "short"),
  JetVectorBranch("JetPUValue", "float"),
]

first2jetmomenta = FirstNJetMomenta(2, jetpt, jeteta, jetphi, jetmass)
onlyjetmomentum = FirstNJetMomenta(1, jetpt, jeteta, jetphi, jetmass, rejectiftoomany=True)
first2jetmomentajecup = FirstNJetMomenta(2, jetpt, jeteta, jetphi, jetmass, sigmabranch=jetsigma, jecdirection=1)
onlyjetmomentumjecup = FirstNJetMomenta(1, jetpt, jeteta, jetphi, jetmass, rejectiftoomany=True, sigmabranch=jetsigma, jecdirection=1)
first2jetmomentajecdn = FirstNJetMomenta(2, jetpt, jeteta, jetphi, jetmass, sigmabranch=jetsigma, jecdirection=-1)
onlyjetmomentumjecdn = FirstNJetMomenta(1, jetpt, jeteta, jetphi, jetmass, rejectiftoomany=True, sigmabranch=jetsigma, jecdirection=-1)
dijetmass = FirstJetsVariable("DiJetMass", lambda jet1, jet2: (jet1+jet2).M(), np.float32, first2jetmomenta, -99)
dijetdeta = FirstJetsVariable("DiJetDEta", lambda jet1, jet2: jet1.Eta() - jet2.Eta(), np.float32, first2jetmomenta, -99)

branches += [
  first2jetmomenta,
  onlyjetmomentum,
  first2jetmomentajecup,
  onlyjetmomentumjecup,
  first2jetmomentajecdn,
  onlyjetmomentumjecdn,
  dijetmass,
  dijetdeta,
  NotRecalculatedBranch("DiJetFisher", -99, np.float32),
]

MELAJECnominal2jets = MELAEventDummyBranch(first2jetmomenta)
MELAJECnominal1jet = MELAEventDummyBranch(onlyjetmomentum)
MELAJECup2jets = MELAEventDummyBranch(first2jetmomentajecup)
MELAJECup1jet = MELAEventDummyBranch(onlyjetmomentumjecup)
MELAJECdn2jets = MELAEventDummyBranch(first2jetmomentajecdn)
MELAJECdn1jet = MELAEventDummyBranch(onlyjetmomentumjecdn)

branches += [
  MELAJECnominal1jet,
  MELAJECnominal2jets,
  MELAJECup1jet,
  MELAJECup2jets,
  MELAJECdn1jet,
  MELAJECdn2jets,
]

MELAbranchesdict = {}

nominalbranches = []
upbranches = []
downbranches = []

with open(os.path.join(os.environ["CMSSW_BASE"], "src", "ZZAnalysis", "AnalysisStep", "test", "prod", "pyFragments", "RecoProbabilities.py")) as f:
  for line in f:
    line = line.split("#")[0]
    if "Name:" not in line: continue
    options = dict(re.findall(r"(\w+):([^ '\"]+)", line))
    if "BestDJJ" in options["Name"] or "ttHUn" in options["Name"]:
      branches.append(NotRecalculatedBranch("p_"+options["Name"], -999, np.float32))
      continue
    if options["Production"] in ("ZZGG", "Lep_ZH", "Lep_WH", "ZZQQB", "ZZINDEPENDENT"): continue
    if options["Process"] == "bkgZJets": continue

    if "Options" in options:
      options.update(dict(re.findall(r"(\w+)=([^;]+)", options["Options"])))

    couplings = {}
    if "Couplings" in options:
      for couplingsetting in options["Couplings"].split(";"):
        name, value = couplingsetting.split("=")
        value = complex(value.replace(",", "+")+"j")
        couplings[name] = value

    initialQQ = False
    if "ForceIncomingFlavors" in options:
      assert options["ForceIncomingFlavors"] == "-21,-21", options["ForceIncomingFlavors"]
      initialQQ = True

    melafunctionname = {
      "JHUGen": "computeProdP",
      "MCFM": "computeProdDecP",
    }[options["MatrixElement"]]
    melafunctionargs = ()
    if int(options.get("isPMaVJJ", 0)):
      melafunctionname = "computeDijetConvBW"
    if int(options.get("isPMaVJJTrue", 0)):
      melafunctionname = "computeDijetConvBW"
      melafunctionargs = True,

    for jec, lst in ("Nominal", nominalbranches), ("Up", upbranches), ("Dn", downbranches):
      name = options["Name"]
      assert "Nominal" in name
      name = name.replace("Nominal", jec)

      inputevent = {
        ("J1JECNominal", "Nominal"): MELAJECnominal1jet,
        ("J2JECNominal", "Nominal"): MELAJECnominal2jets,
        ("J1JECNominal", "Up"): MELAJECup1jet,
        ("J2JECNominal", "Up"): MELAJECup2jets,
        ("J1JECNominal", "Dn"): MELAJECdn1jet,
        ("J2JECNominal", "Dn"): MELAJECdn2jets,
      }[options["Cluster"], jec]

      subtractbranches = []
      if "SubtractP" in options:
        subtractbranches = [MELAbranchesdict[_.replace("Nominal", jec)] for _ in options["SubtractP"].split(",") if "ttH" not in _]

      lst.append(
        MELAProbability(
          "p_"+name,
          inputevent,
          (options["Process"], options["MatrixElement"], options["Production"]),
          couplings,
          melafunctionname,
          melafunctionargs,
          subtractbranches,
          initialQQ,
        )
      )
      MELAbranchesdict[name] = lst[-1]

      if "AddPConst" in options:
        assert int(options["AddPConst"]) == 1
        lst.append(MELAPConst(MELAbranchesdict[name]))
      if "AddPAux" in options:
        assert int(options["AddPAux"]) == 1
        lst.append(MELAPAux(MELAbranchesdict[name]))

branches += nominalbranches + upbranches + downbranches

def fixjetid(infile, outfile, applyid=True, applypuid=True, doCRZLL=False, test_no_mela=False, firstevent=0, lastevent=float("inf"), debug=False):
  print "Processing", infile, "-->", outfile
  folders = ["ZZTree"]
  if doCRZLL: folders.append("CRZLLTree")

  cache = []
  with TFile(infile) as f, TFile(outfile, "CREATE", deleteifbad=True) as newf:
    for foldername in folders:
      folder = f.Get(foldername)
      newfolder = newf.mkdir(foldername)

      if firstevent <= 0:  #so that we can hadd without duplicating things
        Counters = folder.Get("Counters")
        Counters.SetDirectory(newfolder)
        cache.append(Counters)
        failedt = folder.Get("candTree_failed")
        if failedt:
          newfailedt = failedt.CloneTree(-1, "fast")
          newfailedt.SetDirectory(newfolder)
          cache.append(failedt)
          cache.append(newfailedt)

      t = folder.Get("candTree")
      newt = t.CloneTree(0, "fast")
      newt.SetDirectory(newfolder)
      cache.append(t)
      cache.append(newt)

      for branch in branches:
        branch.attachtotree(newt)

      nentries = t.GetEntries()
      nxchecks = 0
      nbadxchecks = collections.Counter()
      worstbadxcheck = collections.defaultdict(lambda: None)
      for i, entry in enumerate(t):
        if i < firstevent: continue
        if i > lastevent: break
        doxcheck = (all(t.JetID) or not applyid) and (all(t.JetPUID) or not applypuid)
        nxchecks += doxcheck
        for branch in branches:
          try:
            branch.setbranchvalue(t, applyid, applypuid, doxcheck)
          except FailedXcheckError as e:
            if debug:
              t.Show()
              print list(t.JetID)
              print list(t.JetPUID)
              print list(t.JetPt)
              print list(t.JetEta)
              print list(t.JetPhi)
              print list(t.JetMass)
              print list(t.JetSigma)
              print list(t.LepLepId)
              print list(t.LepPt)
              print list(t.LepEta)
              print list(t.LepPhi)
              raise
            nbadxchecks[branch.name] += 1
        newt.Fill()
        mela.resetInputEvent()

        if (i+1)%10000 == 0 or (i+1) == nentries or True:
          print i+1, "/", nentries

    print "Done!"
    print "Did xchecks on", nxchecks, " branches that had all jets passing ID and PUID"
    if nbadxchecks:
      print "The following branches had some disagreements:"
      for k, v in sorted(nbadxchecks.iteritems(), key=lambda x: (x[1], x[0])):
        print "{:50} {}".format(k, v)
    else:
      print "All xchecks pass"
    print

if __name__ == "__main__":
  fixjetid(**args.__dict__)

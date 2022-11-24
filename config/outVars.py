import awkward as ak
def fillbranch_R(events,goodevent,weights,reg="sr"):
    events = events[goodevent]
    doc = {"DiJetMass":events.DiJetMass[:,0],"DiJetPt":events.DiJetPt[:,0],"DiJetEta":events.DiJetEta[:,0],"DiJetPhi":events.DiJetPhi[:,0],
           "Jet1Pt":events.bjetpt[:,0],"Jet1Eta":events.bjeteta[:,0],"Jet1Phi":events.bjetphi[:,0],
           "Jet2Pt":events.bjetpt[:,1],"Jet2Eta":events.bjeteta[:,1],"Jet2Phi":events.bjetphi[:,1],
           "MET":events.st_METXYCorr_Met,
           "CaloMET":events.st_pfpatCaloMETPt,
           "isak4JetBasedHemEvent":events.st_isak4JetBasedHemEvent,
           "isak8JetBasedHemEvent":events.st_isak8JetBasedHemEvent,
           "ismetphiBasedHemEvent1":events.st_ismetphiBasedHemEvent1,
           "ismetphiBasedHemEvent2":events.st_ismetphiBasedHemEvent2,
           "eventId":events.st_eventId,"runId":events.st_runId, "lumiSection":events.st_lumiSection,
           "minDphi":ak.to_numpy(events.minDphi_jetMet),
           "DPhi_trkpfMET":events.Dphi_trkpfMet,
           "weight":weights.weight()[goodevent] 
           #"METSFUp":weights.weight("metSFUp")[goodevent],"METSFDown":weights.weight("metSFDown")[goodevent],
           #"btagSFUp":weights.weight("btagSFUp")[goodevent],"btagSFDown":weights.weight("btagSFDown")[goodevent],
           #"fakebSFUp":weights.weight("fakebSFUp")[goodevent],"fakebSFDown":weights.weight("fakebSFDown")[goodevent],
           #"pileupSFUp":weights.weight("pileupSFUp")[goodevent],"pileupSFDown":weights.weight("pileupSFDown")[goodevent],
           #"l1prefireUp":weights.weight("l1prefireUp")[goodevent],"l1prefireDown":weights.weight("l1prefireDown")[goodevent]
          }
    if reg=="tope" : 
         doc["RECOIL"] = events.werecoilPt[:,0]
    if reg=="topmu": 
         doc["RECOIL"] = events.wmurecoilPt[:,0]
    if reg=="zmumu": 
         doc["RECOIL"] = events.zmumurecoilPt[:,0]
    if reg=="zee"  : 
         doc["RECOIL"] = events.zeerecoilPt[:,0]
    return doc
         


def fillbranch_B(events,goodevent,weights,isCR=False):
    events = events[goodevent]
    doc = {"FJetPt":events.fjetptsel[:,0],
           "FJetEta":events.fjetetasel[:,0],
           "FJetPhi":events.fjetphisel[:,0],
           "FJetMass":events.fjetmasssel[:,0],
           "FJetCSV":events.fjetcsvsel[:,0],
           "MET":events.st_METXYCorr_Met,
           "CaloMET":events.st_pfpatCaloMETPt,
           "isak4JetBasedHemEvent":events.st_isak4JetBasedHemEvent,
           "isak8JetBasedHemEvent":events.st_isak8JetBasedHemEvent,
           "ismetphiBasedHemEvent1":events.st_ismetphiBasedHemEvent1,
           "ismetphiBasedHemEvent2":events.st_ismetphiBasedHemEvent2,
           "eventId":events.st_eventId,"runId":events.st_runId, "lumiSection":events.st_lumiSection,
           "min_dPhi":ak.to_numpy(events.minDphi_jetMet),
           "DPhi_trkpfMET":events.Dphi_trkpfMet,
           "weight":weights.weight()[goodevent]
           #"METSFUp":weights.weight("metSFUp")[goodevent],"METSFDown":weights.weight("metSFDown")[goodevent],
           #"btagSFUp":weights.weight("btagSFUp")[goodevent],"btagSFDown":weights.weight("btagSFDown")[goodevent],
           #"fakebSFUp":weights.weight("fakebSFUp")[goodevent],"fakebSFDown":weights.weight("fakebSFDown")[goodevent],
           #"pileupSFUp":weights.weight("pileupSFUp")[goodevent],"pileupSFDown":weights.weight("pileupSFDown")[goodevent],
           #"l1prefireUp":weights.weight("l1prefireUp")[goodevent],"l1prefireDown":weights.weight("l1prefireDown")[goodevent]
          }
    return doc

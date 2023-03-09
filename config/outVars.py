import awkward as ak
def fillbranch_R(events,goodevent,weights,reg="sr"):
    events = events[goodevent]
    #test = events.st_THINjetUncSources[:,0]
    #print ("fjetcsv ",events.st_fjetProbHbb)
    #print ("Jet1Pt",events.bjetpt)
    if not len(events)==0:
        fjetcsv = events.st_fjetProbHbb
        fjetpt = events.fjetpt
    else:
        fjetcsv=events.st_fjetProbHbb[:,0]
        fjetpt = events.fjetpt[:,0]
    #lenght = [[999]]*len(events)
    #print ('lenght',len(events))
    #print ("testing",ak.concatenate([ak.Array(fjetcsv), ak.Array(lenght)], axis=1)) 
    doc = {"DiJetMass":events.DiJetMass[:,0],"DiJetPt":events.DiJetPt[:,0],"DiJetEta":events.DiJetEta[:,0],"DiJetPhi":events.DiJetPhi[:,0],
           "Jet1Pt":events.bjetpt[:,0],"Jet1Eta":events.bjeteta[:,0],"Jet1Phi":events.bjetphi[:,0],"Jet1CSV":events.bjetcsv[:,0],
           "Jet2Pt":events.bjetpt[:,1],"Jet2Eta":events.bjeteta[:,1],"Jet2Phi":events.bjetphi[:,1],"Jet2CSV":events.bjetcsv[:,1],
           "nJets":ak.num(events.jetpt),#"nfjet":ak.num(events.fjetptsel),"fjetcsv":events.fjetcsvsel,"fjetpt":events.fjetptsel,
           "MET":events.st_METXYCorr_Met,"met_Phi":events.st_METXYCorr_MetPhi,
           "CaloMET":events.st_pfpatCaloMETPt,"CaloMETPhi":events.st_pfpatCaloMETPhi,
           "isak4JetBasedHemEvent":events.st_isak4JetBasedHemEvent,
           "isak8JetBasedHemEvent":events.st_isak8JetBasedHemEvent,
           "ismetphiBasedHemEvent1":events.st_ismetphiBasedHemEvent1,
           "ismetphiBasedHemEvent2":events.st_ismetphiBasedHemEvent2,
           "event":events.st_eventId,"run":events.st_runId, "lumi":events.st_lumiSection,
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
           "FJetCSV":events.fjetcsvsel[:,0],"nJets":ak.num(events.isojetpt),
           "MET":events.st_METXYCorr_Met,"met_Phi":events.st_METXYCorr_MetPhi,
           "CaloMET":events.st_pfpatCaloMETPt,"CaloMETPhi":events.st_pfpatCaloMETPhi,
           "isak4JetBasedHemEvent":events.st_isak4JetBasedHemEvent,
           "isak8JetBasedHemEvent":events.st_isak8JetBasedHemEvent,
           "ismetphiBasedHemEvent1":events.st_ismetphiBasedHemEvent1,
           "ismetphiBasedHemEvent2":events.st_ismetphiBasedHemEvent2,
           "event":events.st_eventId,"run":events.st_runId, "lumi":events.st_lumiSection,
           "minDphi":ak.to_numpy(events.minDphi_jetMet),
           "DPhi_trkpfMET":events.Dphi_trkpfMet,
           "weight":weights.weight()[goodevent]
           #"METSFUp":weights.weight("metSFUp")[goodevent],"METSFDown":weights.weight("metSFDown")[goodevent],
           #"btagSFUp":weights.weight("btagSFUp")[goodevent],"btagSFDown":weights.weight("btagSFDown")[goodevent],
           #"fakebSFUp":weights.weight("fakebSFUp")[goodevent],"fakebSFDown":weights.weight("fakebSFDown")[goodevent],
           #"pileupSFUp":weights.weight("pileupSFUp")[goodevent],"pileupSFDown":weights.weight("pileupSFDown")[goodevent],
           #"l1prefireUp":weights.weight("l1prefireUp")[goodevent],"l1prefireDown":weights.weight("l1prefireDown")[goodevent]
          }
    return doc

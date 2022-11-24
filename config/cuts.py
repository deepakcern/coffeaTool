

regions_R = {"sr": {"trigger":True,"noElectron":True, "noMuon":True, "noPhoton":True,"noTau":True,"metcut":True,"nJets":True,"twobJets":True,"bmass":True,"minDphi":True},
             "qcd":{"trigger":True,"noElectron":True, "noMuon":True, "noPhoton":True,"noTau":True,"metcut":True,"nJets":True,"twobJets":True,"bmass":True,"invrtminDphi":True}}



regions_B = {"sr": {"trigger":True,"noElectron":True, "noMuon":True, "noPhoton":True,"noTau":True,"metcut250":True,"nIsojet":True,"noIsobjet":True,"nfjet":True,"minDphi":True},
            "qcd":{"trigger":True,"noElectron":True, "noMuon":True, "noPhoton":True,"noTau":True,"metcut250":True,"nIsojet":True,"noIsobjet":True,"nfjet":True,"invrtminDphi":True}}


cutflow_R = {"trigger":{"trigger"},"eleVeto":{"trigger","noElectron"},"muonVeo":{"trigger","noElectron","noMuon"},"TauVeto":{"trigger","noElectron","noMuon","noTau"},"phoVeto":{"trigger","noElectron","noMuon","noPhoton"},"metcut":{"trigger","noElectron","noMuon","noPhoton","metcut"},"nJets":{"trigger","noElectron","noMuon","noPhoton","metcut","nJets"},"bJets":{"trigger","noElectron","noMuon","noPhoton","metcut","nJets","twobJets"},"mass":{"trigger","noElectron","noMuon","noPhoton","metcut","nJets","twobJets","bmass"},"minDphi":{"trigger","noElectron","noMuon","noPhoton","metcut","nJets","twobJets","bmass","minDphi"}}


cutflow_B = {"trigger":{"trigger"},"eleVeto":{"trigger","noElectron"},"muonVeo":{"trigger","noElectron","noMuon"},"TauVeto":{"trigger","noElectron","noMuon","noTau"},"phoVeto":{"trigger","noElectron","noMuon","noPhoton"},"metcut":{"trigger","noElectron","noMuon","noPhoton","metcut250"},"nfjet":{"trigger","noElectron","noMuon","noPhoton","metcut250","nfjet"},"nIsojet":{"trigger","noElectron","noMuon","noPhoton","metcut250","nfjet","nIsojet"},"noIsobjet":{"trigger","noElectron","noMuon","noPhoton","metcut250","nfjet","nIsojet","noIsobjet"},"minDphi":{"trigger","noElectron","noMuon","noPhoton","metcut250","nfjet","nIsojet","noIsobjet","minDphi"}}

import awkward as ak
from coffea import processor
from coffea.nanoevents.methods import candidate
import hist
from coffea.analysis_tools import PackedSelection
import numpy as np
from coffea.btag_tools import BTagScaleFactor
from coffea.analysis_tools import Weights
from coffea.lookup_tools.dense_lookup import dense_lookup
import data_2017.SFFactory_2017 as sf17
import data_2017.SFFactorySystUp_2017 as sf17Up
import data_2017.SFFactorySystDown_2017 as sf17Down

import data_2018.SFFactory_2018 as sf18
import data_2018.SFFactorySystUp_2018 as sf18Up
import data_2018.SFFactorySystDown_2018 as sf18Down

from config.functions import *
from config.outVars import *
from coffea.lookup_tools import extractor
import uproot
from config.cuts import *


class monoHbbProcessor(processor.ProcessorABC):
    def __init__(self):
        pass

    def setupYearDependency(self,events):
        self.year = events.metadata["year"]
        if self.year=="2017":
            self.LWP = 0.1522
            self.MWP = 0.4941
            self.ResolvedmetSF = {"center":sf17.R_metTrig_firstmethod,"Up":sf17Up.R_metTrig_firstmethod_SystUp,"Down":sf17Down.R_metTrig_firstmethod_SystDown}
            self.ResolvedmetSFRange = sf17.R_metTrig_firstmethod_X_range
            self.BoostedmetSF  = {"center":sf17.B_metTrig_firstmethod,"Up":sf17Up.B_metTrig_firstmethod_SystUp,"Down":sf17Down.B_metTrig_firstmethod_SystDown}
            self.BoostedmetSFRange = sf17.B_metTrig_firstmethod_X_range
            self.PileReweightHisto  = {"center":sf17.pileup2017histo,"Up":sf17Up.pileup2017histo_SystUp,"Down":sf17Down.pileup2017histo_SystDown}
            self.PileReweightRange  = sf17.pileup2017histo_X_range
            self.btagEffroot   = "data_2017/bTagEffs_2017.root"
            self.btagCSVfile   = "data_2017/DeepCSV_94XSF_V5_B_F.csv.gz"

        if self.year=="2018":
            self.LWP = 0.1241
            self.MWP = 0.4184        

            self.ResolvedmetSF = {"center":sf18.R_metTrig_firstmethod,"Up":sf18Up.R_metTrig_firstmethod_SystUp,"Down":sf18Down.R_metTrig_firstmethod_SystDown}
            self.ResolvedmetSFRange = sf18.R_metTrig_firstmethod_X_range

            self.BoostedmetSF  = {"center":sf18.B_metTrig_firstmethod,"Up":sf18Up.B_metTrig_firstmethod_SystUp,"Down":sf18Down.B_metTrig_firstmethod_SystDown}
            self.BoostedmetSFRange = sf18.B_metTrig_firstmethod_X_range

            self.PileReweightHisto  = {"center":sf18.pileup2018histo,"Up":sf18Up.pileup2018histo_SystUp,"Down":sf18Down.pileup2018histo_SystDown}
            self.PileReweightRange  = sf18.pileup2018histo_X_range
            self.btagEffroot   = "data_2018/bTagEffs_2018.root"
            self.btagCSVfile   = "data_2018/DeepCSV_102XSF_V2.csv"

    
    def addWeights(self,events):
        weights = {"Boosted":Weights(len(events)),"Resolved":Weights(len(events))}
        #weights = Weights(len(events))
        corr_met = dense_lookup(np.array(self.ResolvedmetSF["center"]),np.array(self.ResolvedmetSFRange))
        corr_met_B = dense_lookup(np.array(self.BoostedmetSF["center"]),np.array(self.BoostedmetSFRange) )
        corr_pu  = dense_lookup(np.array(self.PileReweightHisto["center"]),np.array(self.PileReweightRange))
        
        weights["Boosted"].add("metSF",weight=corr_met_B(events.st_METXYCorr_Met),weightUp=corr_met_B(events.st_METXYCorr_Met),weightDown=corr_met_B(events.st_METXYCorr_Met))
        weights["Boosted"].add("pileupSF",weight=corr_pu(events.st_pu_nTrueInt),weightUp=corr_pu(events.st_pu_nTrueInt),weightDown=corr_pu(events.st_pu_nTrueInt))
    

        weights["Resolved"].add("metSF",weight=corr_met(events.st_METXYCorr_Met),weightUp=corr_met(events.st_METXYCorr_Met),weightDown=corr_met(events.st_METXYCorr_Met))
        weights["Resolved"].add("pileupSF",weight=corr_pu(events.st_pu_nTrueInt),weightUp=corr_pu(events.st_pu_nTrueInt),weightDown=corr_pu(events.st_pu_nTrueInt))


        weights["Boosted"].add("l1prefire",weight=events.st_prefiringweight,weightUp=events.st_prefiringweightup,weightDown=events.st_prefiringweightdown)
        weights["Resolved"].add("l1prefire",weight=events.st_prefiringweight,weightUp=events.st_prefiringweightup,weightDown=events.st_prefiringweightdown)
        weights["Boosted"].add("mcweight",weight=events.mcweight)
        weights["Resolved"].add("mcweight",weight=events.mcweight)


        ext = extractor()
        ext.add_weight_sets(["btag_eff_mwp efficiency_btag_mwp "+self.btagEffroot])
        ext.add_weight_sets(["ctag_eff_mwp efficiency_ctag_mwp "+self.btagEffroot])
        ext.add_weight_sets(["lighttag_eff_mwp efficiency_lighttag_mwp "+self.btagEffroot])

        ext.add_weight_sets(["btag_eff_lwp efficiency_btag_lwp "+self.btagEffroot])
        ext.add_weight_sets(["ctag_eff_lwp efficiency_ctag_lwp "+self.btagEffroot])
        ext.add_weight_sets(["lighttag_eff_lwp efficiency_lighttag_lwp "+self.btagEffroot])

        ext.finalize()
        evaluator = ext.make_evaluator()


        '''
        COMPUTE B TAG WEGIHTS FOR RESOLVED CATEGORY [FOR TAG AND NON TAG]
        '''
        MWP = self.MWP
        LWP = self.LWP
        btagIndex  = events.st_THINjetDeepCSV>MWP    
        nonbtagIndex = ~btagIndex

        btag_sf_evaluator = BTagScaleFactor(self.btagCSVfile, "medium")
        
        
        btagCal=getbTagWegiht(btag_sf_evaluator,evaluator,events.st_THINjetHadronFlavor,events.jeteta,events.jetpt,btagIndex,nonbtagIndex,wptype="medium")


        weights["Resolved"].add("btagSF",weight=btagCal["btagSF"],weightUp=btagCal["btagSFUp"],weightDown=btagCal["btagSFDown"])
        weights["Resolved"].add("fakebSF",weight=btagCal["fakebSF"],weightUp=btagCal["fakebSFUp"],weightDown=btagCal["fakebSFDown"])

        '''
        ================= Compute b tag weights for boosted category ====
        '''
        isobtag_sf_evaluator = BTagScaleFactor(self.btagCSVfile, "loose")
        isobtagIndex  = events.isojetDeepCSV>LWP
        isononbtagIndex = ~isobtagIndex

        isobtagCal=getbTagWegiht(isobtag_sf_evaluator,evaluator,events.isojetHadronFlavor,events.isojeteta,events.isojetpt,isobtagIndex,isononbtagIndex,wptype="loose")        
        weights["Boosted"].add("btagSF",weight=isobtagCal["btagSF"],weightUp=isobtagCal["btagSFUp"],weightDown=isobtagCal["btagSFDown"])
        weights["Boosted"].add("fakebSF",weight=isobtagCal["fakebSF"],weightUp=isobtagCal["fakebSFUp"],weightDown=isobtagCal["fakebSFDown"])        

        if events.metadata['dataset']=="TT":
            weights["Boosted"].add("topPtRewieght",weight=getTopPtReWgt(events.st_genParPt))
            weights["Resolved"].add("topPtRewieght",weight=getTopPtReWgt(events.st_genParPt))
        else:
            weights["Boosted"].add("topPtRewieght",np.ones(len(events)))


        return weights


#    def isobJetweight():
        

  
    def updateJetColl(self,events):
            ak4jeteta           =   geteta(events.st_THINjetPx,events.st_THINjetPy,events.st_THINjetPz)
            ak4jeteta_index     =   np.abs(ak4jeteta)<2.5
            events['st_THINjetPx']        =   events.st_THINjetPx[ak4jeteta_index]
            events['st_THINjetPy']        =   events.st_THINjetPy[ak4jeteta_index]
            events['st_THINjetPz']        =   events.st_THINjetPz[ak4jeteta_index]
            events['st_THINjetEnergy']    =   events.st_THINjetEnergy[ak4jeteta_index]
            events['st_THINjetDeepCSV']   =   events.st_THINjetDeepCSV[ak4jeteta_index]
            events['st_THINjetHadronFlavor']= events.st_THINjetHadronFlavor[ak4jeteta_index]
            events['st_THINjetCorrUnc']   =   events.st_THINjetCorrUnc[ak4jeteta_index]
            events['st_THINbRegNNCorr']   =   events.st_THINbRegNNCorr[ak4jeteta_index]

            events['st_THINjetNHadEF']    =   events.st_THINjetNHadEF[ak4jeteta_index]
            events['st_THINjetCHadEF']    =   events.st_THINjetCHadEF[ak4jeteta_index]
            events['st_THINjetCEmEF']     =   events.st_THINjetCEmEF[ak4jeteta_index]
            return events

 
    def addMainColoumns(self,events):

        events['jetpt']     =getpt(events.st_THINjetPx,events.st_THINjetPy)
        events['jeteta']    =geteta(events.st_THINjetPx,events.st_THINjetPy,events.st_THINjetPz)
        events['jetphi']    =getphi(events.st_THINjetPx,events.st_THINjetPy)


        events["minDphi_jetMet"]= np.abs(ak.min(DeltaPhi(events.jetphi, events.st_METXYCorr_MetPhi), axis=-1))
        events["Dphi_trkpfMet"] = np.abs(DeltaPhi(events.st_pfTRKMETPhi,events.st_METXYCorr_MetPhi))

        events['fjetpt']    =getpt(events.st_fjetPx,events.st_fjetPy)
        events['fjeteta']   =geteta(events.st_fjetPx,events.st_fjetPy,events.st_fjetPz)
        events['fjetphi']   =getphi(events.st_fjetPx,events.st_fjetPy)
 
        fjet_sel = (np.abs(events.fjeteta) <2.5) & (events.st_fjetSDMass>70) & (events.st_fjetSDMass < 150) & (events.st_fjetProbHbb > 0.86)
        events['nfjetsel'] = ak.num(events.fjetpt[fjet_sel])
        events['fjetptsel']= events.fjetpt[fjet_sel]
        events['fjetetasel']= events.fjeteta[fjet_sel]
        events['fjetphisel']= events.fjetphi[fjet_sel]
        events['fjetmasssel'] =events.st_fjetSDMass[fjet_sel]
        events['fjetcsvsel'] = events.st_fjetProbHbb[fjet_sel]

        MWP = self.MWP#0.4941
        bjetCond            =(events.st_THINjetDeepCSV>MWP) 
        events['bjetpt']    =events.jetpt[bjetCond]
        events['bjeteta']   =events.jeteta[bjetCond]
        events['bjetphi']   =events.jetphi[bjetCond]
        events['bjetHadronFla'] = events.st_THINjetHadronFlavor[bjetCond]
        events['bjetE']       =events.st_THINjetEnergy[bjetCond]
        
        events['DiJetMass']     = getMassPair(events.st_THINjetPx[bjetCond],events.st_THINjetPy[bjetCond],events.st_THINjetPz[bjetCond],events.st_THINjetEnergy[bjetCond])
        parivars = getPair_ptetaphi(events.st_THINjetPx[bjetCond],events.st_THINjetPy[bjetCond],events.st_THINjetPz[bjetCond],events.st_THINjetEnergy[bjetCond])
        events['DiJetPt'] = parivars["pt"]
        events['DiJetEta'] = parivars["eta"]
        events['DiJetPhi'] = parivars["phi"]


        LWP = self.LWP#0.1522


        cleaned = isclean(events.jeteta,events.fjeteta,events.jetphi,events.fjetphi,cut_=0.8)
        events['isojetpt']      = events.jetpt[cleaned]
        events['isojeteta']     = events.jeteta[cleaned]
        events['isojetphi']     = events.jetphi[cleaned]
        events['isojetDeepCSV'] = events.st_THINjetDeepCSV[cleaned]
        events['isojetHadronFlavor'] =events.st_THINjetHadronFlavor[cleaned]

        bjetCondForBoosted      = events.isojetDeepCSV>LWP

        events['isobjetpt']     =events.isojetpt[bjetCondForBoosted]
        events['isobjeteta']    =events.isojeteta[bjetCondForBoosted]
        events['isobjetphi']    =events.isojetphi[bjetCondForBoosted]

        events['elept']         = getpt(events.st_elePx, events.st_elePy)	
        events['eleeta']        =geteta(events.st_elePx, events.st_elePy, events.st_elePz)
        events['elephi']        =getphi(events.st_elePx, events.st_elePy)

        events['nlooseEle']     = ak.num(events.elept[(events.st_eleIsPassLoose)])
        events['ntightEle']     = ak.num(events.elept[(events.st_eleIsPassTight) & (events.elept > 40)])

        events['mupt']          = getpt(events.st_muPx,events.st_muPy)	
        events['mueta']         = geteta(events.st_muPx,events.st_muPy,events.st_muPz)
        events['muphi']         = getphi(events.st_muPx,events.st_muPy)

        events['nlooseMu']      = ak.num(events.mupt)
        events['ntightMu']      = ak.num(events.mupt[(events.st_isTightMuon) & (events.mupt>30)])

        events['phopt']         = getpt(events.st_phoPx, events.st_phoPy)
        events['phophi']        = getphi(events.st_phoPx, events.st_phoPy)
        events['phoeta']        = geteta(events.st_phoPx, events.st_phoPy,events.st_phoPz)

        cleanedPho = (isclean(events.phoeta,events.jeteta,events.phophi,events.jetphi,cut_=0.4)) & (events.phopt>20)
        events['npho']          = ak.num(events.phopt[cleanedPho])
        
        events['werecoilPt']    = getrecoil1(events.st_elePx,events.st_elePy,events.st_pfMetCorrPt,events.st_pfMetCorrPhi) 
        events['wmurecoilPt']   = getrecoil1(events.st_muPx,events.st_muPy,events.st_pfMetCorrPt,events.st_pfMetCorrPhi) 
        events['zeerecoilPt']   = getrecoil2(events.st_elePx,events.st_elePy,events.st_pfMetCorrPt,events.st_pfMetCorrPhi)
        events['zmumurecoilPt'] = getrecoil2(events.st_muPx,events.st_muPy,events.st_pfMetCorrPt,events.st_pfMetCorrPhi)
        
        events['isJetBasedHem'] = isJetBasedHemEvent(self.year,events.metadata['isData'],events.st_isak4JetBasedHemEvent,events.st_isak8JetBasedHemEvent)
        events['isMetBasedHem'] = isLowmetBasedHemEvent(self.year,events.metadata['isData'],events.st_ismetphiBasedHemEvent1)

                

        return events
 
    def process(self, events):
        print ('\n'+"================ Processor  is running now =================="+'\n')

        dataset = events.metadata['dataset']
        #events = events[(events.st_runId==306154) & (events.st_lumiSection==676) & (events.st_eventId==1161700016)]
        self.setupYearDependency(events)
        events  = self.updateJetColl(events)
        events  = self.addMainColoumns(events)
        '''
        ========================
        SISNAL REGION SELECTIONS
        =======================
        '''
        selection = PackedSelection()

        selection.add("trigger", (events.st_mettrigdecision) & (~events.isJetBasedHem) & (~events.isMetBasedHem))
        selection.add("noElectron", (events.nlooseEle == 0 ) & (events.st_mettrigdecision))
        selection.add("noMuon", events.nlooseMu == 0)
        selection.add("noTau", events.st_nTau_discBased_looseElelooseMuVeto == 0)
        selection.add("noPhoton", (events.npho == 0) & (events.st_nTau_discBased_looseElelooseMuVeto == 0))
        selection.add("metcut", (events.st_METXYCorr_Met > 200))# & events.st_mettrigdecision)
        selection.add("metcut250", events.st_METXYCorr_Met > 250)
        selection.add("nJets", ak.num(events.jetpt)<=4)
        selection.add("twobJets", (ak.num(events.bjetpt)==2) & ak.any(events.bjetpt >= 50.0, axis=1))

        selection.add("leadbJetPt50", ak.any(events.bjetpt >= 50.0, axis=1))
        selection.add("bmass", ak.any(events.DiJetMass < 150,axis=1) & ak.any(events.DiJetMass > 100,axis=1))
        selection.add("minDphi",events.minDphi_jetMet>0.4)
        selection.add("invrtminDphi",events.minDphi_jetMet<0.4)
        selection.add("nfjet",events.nfjetsel==1)
        selection.add("noIsobjet",ak.num(events.isobjetpt)==0)
        selection.add("nIsojet",ak.num(events.isojetpt)<=2)

        selection.add("OneElectron", events.ntightEle==1)
        selection.add("OneMuon",    events.ntightMu==1)
        #selection.add("isJetBasedHem", ~events.isJetBasedHem)
        #selection.add("isMetBasedHem", ~events.isMetBasedHem)
        #selection.add("Recoil200", (ak.any(events.werecoilPt>200, axis=1)) | (ak.any(events.wmurecoilPt>200, axis=1)) | (ak.any(events.zmumurecoilPt>200, axis=1)) | (ak.any(events.zeerecoilPt>200, axis=1)))
        #selection.add("Recoil250", ak.any(events.werecoilPt>250) | ak.any(events.wmurecoilPt>250) | ak.any(events.zmumurecoilPt>250) | ak.any(events.zeerecoilPt>250))

        '''
        --------------------------------------------
        GET EVENT WEIGHTS FOR MC AND SET 1 FOR DATA
        -------------------------------------------
        '''

        if not events.metadata['isData']:
            weights = self.addWeights(events)
        else:
            weights = {"Boosted":Weights(len(events)),"Resolved":Weights(len(events))} 


        '''
        ============ START FILLING TREE ===========
        '''
	
        fout = uproot.recreate(events.metadata['outputpath']+'/'+events.metadata['filename'])

        for region, cuts in regions_B.items():
           goodevent = selection.require(**cuts)
           if region.startswith("sr"):
                fout["monoHbb_SR_boosted"] = fillbranch_B(events,goodevent,weights["Boosted"])
           if region.startswith("qcd"):
                fout["monoHbb_qcd_boosted"] = fillbranch_B(events,goodevent,weights["Boosted"])


        for region, cuts in regions_R.items():
            goodevent = (selection.require(**cuts)) & (~selection.require(**regions_B[region]))
         #   weights = self.addWeights(events)[goodevent]
            if region.startswith("sr"):
                fout["monoHbb_SR_resolved"] = fillbranch_R(events,goodevent,weights["Resolved"])
            if region.startswith("qcd"):
                fout["monoHbb_qcd_resolved"] = fillbranch_R(events,goodevent,weights["Resolved"])

        totalevents      =events.metadata['totalevents']
        totalweightevents=events.metadata['totalweightevents']
        h_total_mcweight = (hist.Hist.new.Reg(2, 0, 2, name="h_total_mcweight", label="").Weight())
        h_total          = (hist.Hist.new.Reg(2, 0, 2, name="h_total", label="").Weight())
        h_cutflow_R      = (hist.Hist.new.Reg(12, 0, 12, name="h_total", label="").Weight())
        h_cutflow_B      = (hist.Hist.new.Reg(12, 0, 12, name="h_total", label="").Weight())
        h_total_mcweight.fill(1)
        h_total.fill(1)
        fout["h_total_mcweight"] = h_total_mcweight * totalweightevents
        fout["h_total"]          = h_total * totalevents
        


        '''
        ============ CUTFLOW ===========
        '''
        values = {"Boosted":[],"Resolved":[]}
        for lebel, cut in cutflow_B.items():
             goodevent = selection.all(*(cut))
             #nev       = np.sum(goodevent)
             nev       = weights["Boosted"].weight()[goodevent].sum()
             values["Boosted"].append(nev)
             print(f"Boosted Events passing for : {lebel}: {nev}")

        print ("")
        for lebel, cut in cutflow_R.items():
             goodevent = (selection.all(*(cut))) & (~selection.require(**regions_B["sr"]))
             #nev       = np.sum(goodevent)
             nev       = weights["Resolved"].weight()[goodevent].sum()
             values["Resolved"].append(nev)
             print(f"Resolved Events passing for {lebel}: {nev}")
             
        fout["h_cutflow_sr_boosted"] = fillcutflow(h_cutflow_B,values["Boosted"])
        fout["h_cutflow_sr_resolved"]=fillcutflow(h_cutflow_R,values["Resolved"])
        fout.close()

        return {
            dataset: {
                "entries":totalweightevents, #len(events),
                #"mass": masshist,
        
              }
              }

    def postprocess(self, accumulator):
        pass


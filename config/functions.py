import numpy as np
import awkward as ak
import numpy
def isClean(obj_A, obj_B, drmin=0.4):
    objB_near, objB_DR = obj_A.nearest(obj_B, return_metric=True)
    mask = ak.fill_none(objB_DR > drmin, True)
    return (mask)

def geteta(mupx, mupy,mupz):
    mup = np.sqrt(mupx**2 + mupy**2 + mupz**2)
    mueta = np.log((mup + mupz)/(mup - mupz))/2
    return (mueta)



def getphi(mupx, mupy):
    muphi = np.arctan2(mupy, mupx)
    return (muphi)



def getpt(mupx, mupy):
    mupt = np.sqrt(mupx**2 + mupy**2)
    return (mupt)


def Delta_R(eta1, eta2, phi1,phi2):
    deltaeta = eta1-eta2
    deltaphi = DeltaPhi(phi1,phi2)
    DR = numpy.sqrt ( deltaeta**2 + deltaphi**2 )
    return DR 

def jetcleaning(ak4eta, lepeta, ak4phi, lepphi, DRCut):
    ## usage: (obj_to_clean, obj_cleaned_against, so on
    dr_ = Delta_R(ak4eta, lepeta, ak4phi, lepphi)
    
    return (dr_ > DRCut)


def Phi_mpi_pi(x):
    y = numpy.add(x, numpy.pi)
    y = numpy.mod(y, 2*numpy.pi)
    y = numpy.subtract(y, numpy.pi)
    return y

def DeltaPhi(phi1,phi2):
    #print (phi1-phi2)
    phi = Phi_mpi_pi(phi1-phi2)
    return abs(phi)


def isclean(phoeta, jeteta, phophi, jetphi, cut_=0.4):
    phoeta_unzip, jeteta_unzip = ak.unzip(ak.cartesian([phoeta,jeteta], nested=True))
    phophi_unzip, jetphi_unzip = ak.unzip(ak.cartesian([phophi,jetphi], nested=True))
    deta_unzip = phoeta_unzip - jeteta_unzip
    dphi_unzip = Phi_mpi_pi(phophi_unzip - jetphi_unzip)
    dr_unzip = numpy.sqrt(deta_unzip**2 + dphi_unzip**2)
    dr_pho_jet_status = ~ak.any(dr_unzip<=cut_,axis=-1)  ## use axis in new version of awkward
    return dr_pho_jet_status


def getrecoil1(elepx_,elepy_,met_,metphi_):
    WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
    WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
    WenuRecoilPt = (numpy.sqrt(WenuRecoilPx**2  +  WenuRecoilPy**2))
    #print ('elepx_',elepx_)
    return WenuRecoilPt


def getRecoilPhi1(elepx_,elepy_,met_,metphi_):
    WenuRecoilPx = -( met_*numpy.cos(metphi_) + elepx_)
    WenuRecoilPy = -( met_*numpy.sin(metphi_) + elepy_)
    WenurecoilPhi = numpy.arctan2(WenuRecoilPx,WenuRecoilPy)
    return WenurecoilPhi

def getPair_pxpypz(bjetpx,bjetpy):
    bjets = ak.zip({"px":bjetpx,"py":bjetpy})
    pairs = ak.combinations(bjets, 2)
    bjet1, bjet2 = ak.unzip(pairs)
    return bjet1, bjet2

def getrecoil2(elepx_,elepy_,met_,metphi_):
    mu1 , mu2  = getPair_pxpypz(elepx_,elepy_)
    return getrecoil1(mu1.px+mu2.px,mu1.py+mu2.py,met_,metphi_)
 
def getN(var_, i):
    return ak.mask(var_, ak.num(var_, axis=1)>i, highlevel=False)[:,i]

def getMassPair(bjetpx,bjetpy,bjetpz,bjetE):
    
    bjets = ak.zip({"px":bjetpx,"py":bjetpy,"pz":bjetpz,"En":bjetE})
    pairs = ak.combinations(bjets, 2)
    bjet1, bjet2 = ak.unzip(pairs)
    mass = numpy.sqrt((bjet1.En+bjet2.En)**2 - (bjet1.px+bjet2.px)**2 - (bjet1.py+bjet2.py)**2-(bjet1.pz+bjet2.pz)**2)
    return mass
def getPair_ptetaphi(bjetpx,bjetpy,bjetpz,bjetE):
    bjets = ak.zip({"px":bjetpx,"py":bjetpy,"pz":bjetpz,"En":bjetE})
    pairs = ak.combinations(bjets, 2)
    bjet1, bjet2 = ak.unzip(pairs)
    Px = (bjet1.px+bjet2.px)
    Py = (bjet1.py+bjet2.py)
    Pz = (bjet1.pz+bjet2.pz)
    pt = getpt(Px,Py)
    eta=geteta(Px,Py,Py)
    phi=getphi(Px,Py)
    pairVars ={"pt":pt,"eta":eta,"phi":phi}
    return pairVars


def getTopPtReWgt(pt):
    w = 0.103*(numpy.exp(-0.0118*pt)) - 0.000134*pt + 0.973
    return numpy.sqrt(ak.prod(w,axis=1))


def getbTagWegiht(sf_evaluator,eff_evaluator,Flavor,jeteta,jetpt,btagIndex,nonbtagIndex,wptype):

    btagweights  = {}

    bjet_sf         = sf_evaluator.eval("central", Flavor, abs(jeteta), jetpt)
    bjet_sfUp         =  sf_evaluator.eval("up", Flavor, abs(jeteta), jetpt)
    bjet_sfDown       =  sf_evaluator.eval("down", Flavor, abs(jeteta), jetpt)   

    #print ("inside function btag sf",bjet_sf)
 
    btag_sf         = bjet_sf[btagIndex]
    btag_sfUp       = bjet_sfUp[btagIndex]
    btag_sfDown     = bjet_sfDown[btagIndex]


    nonbtag_sf      = bjet_sf[nonbtagIndex]
    nonbtag_sfUp    = bjet_sfUp[nonbtagIndex]
    nonbtag_sfDown  = bjet_sfDown[nonbtagIndex]
    arr_bflveff     = (Flavor==5) & (nonbtagIndex)
    arr_cflveff     = (Flavor==4) & (nonbtagIndex)
    arr_lightflveff = (Flavor==0) & (nonbtagIndex)
    
    if wptype=="medium":

        bflveff         = eff_evaluator["btag_eff_mwp"](jeteta,jetpt)[arr_bflveff]
        cflveff         = eff_evaluator["ctag_eff_mwp"](jeteta,jetpt)[arr_cflveff]
        lightflveff     = eff_evaluator["lighttag_eff_mwp"](jeteta,jetpt)[arr_lightflveff]
    else:    
        bflveff         = eff_evaluator["btag_eff_lwp"](jeteta,jetpt)[arr_bflveff]
        cflveff         = eff_evaluator["ctag_eff_lwp"](jeteta,jetpt)[arr_cflveff]
        lightflveff     = eff_evaluator["lighttag_eff_lwp"](jeteta,jetpt)[arr_lightflveff]


    nonbjetweight_b         = (1 - ((bjet_sf[arr_bflveff]) * (bflveff))) / (1 - bflveff)
    nonbjetweight_bUp       = (1 - ((bjet_sfUp[arr_bflveff]) * (bflveff))) / (1 - bflveff)
    nonbjetweight_bDown     = (1 - ((bjet_sfDown[arr_bflveff]) * (bflveff))) / (1 - bflveff)


    nonbjetweight_c         = (1 - ((bjet_sf[arr_cflveff]) * (cflveff))) / (1 - cflveff)
    nonbjetweight_cUp       = (1 - ((bjet_sfUp[arr_cflveff]) * (cflveff))) / (1 - cflveff)
    nonbjetweight_cDown     = (1 - ((bjet_sfDown[arr_cflveff]) * (cflveff))) / (1 - cflveff)


    nonbjetweight_light     = (1 - ((bjet_sf[arr_lightflveff]) * (lightflveff))) / (1 - lightflveff)
    nonbjetweight_lightUp   = (1 - ((bjet_sfUp[arr_lightflveff]) * (lightflveff))) / (1 - lightflveff)
    nonbjetweight_lightDown = (1 - ((bjet_sfDown[arr_lightflveff]) * (lightflveff))) / (1 - lightflveff)

    fake_bSF   = ak.prod(nonbjetweight_b,axis=1) * ak.prod(nonbjetweight_c,axis=1) * ak.prod(nonbjetweight_light,axis=1)
    fake_bSFUp = ak.prod(nonbjetweight_bUp,axis=1) * ak.prod(nonbjetweight_cUp,axis=1) * ak.prod(nonbjetweight_lightUp,axis=1)
    fake_bSFDown = ak.prod(nonbjetweight_bDown,axis=1) * ak.prod(nonbjetweight_cDown,axis=1) * ak.prod(nonbjetweight_lightDown,axis=1)
    #print ("bflveff ",bflveff)
    #print ("cflveff ",cflveff)
    #print ("lightflveff",lightflveff)

    #print ("")
    #print ("nonbjetweight_b",nonbjetweight_b)
    #print ("nonbjetweight_c",nonbjetweight_c)
   # print ("nonbjetweight_light",nonbjetweight_light)

    btagweights["btagSF"] 	= ak.prod(btag_sf,axis=1)
    btagweights["btagSFUp"]	= ak.prod(btag_sfUp,axis=1)
    btagweights["btagSFDown"]	= ak.prod(btag_sfDown,axis=1)
    btagweights["fakebSF"] 	= fake_bSF
    btagweights["fakebSFUp"]	= fake_bSFUp
    btagweights["fakebSFDown"]	= fake_bSFDown

    return btagweights

def fillcutflow(h_hist,values):
    for ii, value in enumerate(values):
        h_hist.fill(ii+1,weight=value)
    return h_hist


def isJetBasedHemMCEvent(st_isak4JetBasedHemEvent,st_isak8JetBasedHemEvent):
    ishem = (st_isak4JetBasedHemEvent) | (st_isak8JetBasedHemEvent)
    np.random.seed(123)
    array_index = np.where(ishem)[0]
    prob =  np.random.rand(1,len(array_index))[0]
    Nothemregion = np.where(prob>0.646)[0]
    array_index = array_index[Nothemregion]
    ak.to_numpy(ishem)[array_index] = [False] * len(array_index)
    return ishem

def isLowmetBasedHemMCEvent(ismetphiBasedHemEvent1):
    np.random.seed(123)
    array_index = np.where(ismetphiBasedHemEvent1)[0]
    prob =  np.random.rand(1,len(array_index))
    Nothemregion = np.where(prob>0.646)[0]
    array_index = array_index[Nothemregion]
    ak.to_numpy(ismetphiBasedHemEvent1)[array_index] = [False] * len(array_index)
    return ismetphiBasedHemEvent1

def isJetBasedHemEvent(year,isData,st_isak4JetBasedHemEvent,st_isak8JetBasedHemEvent):
    if year=="2018":
        if not isData: ishem = isJetBasedHemMCEvent(st_isak4JetBasedHemEvent,st_isak8JetBasedHemEvent)
        if isData:ishem = (st_isak4JetBasedHemEvent) | (st_isak8JetBasedHemEvent) # Note: In root file this is apply for CD Era Only
    else:
        ishem = st_isak4JetBasedHemEvent # Note: For 2017, its False for all events
    return ishem


def isLowmetBasedHemEvent(year,isData,ismetphiBasedHemEvent1):
    if year=="2018":
        if not isData: ishem = isLowmetBasedHemMCEvent(ismetphiBasedHemEvent1)
        if isData:ishem = ismetphiBasedHemEvent1 # Note: In root file this is apply for CD Era Only
    else:
        ishem = ismetphiBasedHemEvent1 # Note: For 2017, its False for all events
    return ishem

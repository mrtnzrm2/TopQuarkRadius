#! /usr/bin/env python
import ROOT
import optparse
import json
import sys
import os
import numpy
from TopLJets2015.TopAnalysis.storeTools import getEOSlslist
from TopLJets2015.TopAnalysis.nuSolutions import *

"""
a dummy converter
"""
def convertToPtEtaPhiM(lVec,xyz,m=0.):    
    en=ROOT.TMath.Sqrt(xyz[0]**2+xyz[1]**2+xyz[2]**2)
    p4=ROOT.TLorentzVector(xyz[0],xyz[1],xyz[2],en)
    return lVec(p4.Pt(),p4.Eta(),p4.Phi(),p4.M())

"""
stop quark/chi0 filter the arguments are the masses to filter out
"""
def stopChiFilter(tree,filtArgs):
    weight=-1.0
    mstop,mchi0=0,0
    for i in xrange(0,tree.nt):
        if abs(tree.t_id[i])==1000006 : mstop=tree.t_m[i]
        if abs(tree.t_id[i])==1000022 : mchi0=tree.t_m[i]
    if mstop==float(filtArgs[0]) and mchi0==float(filtArgs[1]): weight=1.0
    return weight


"""
top radius filter
"""
def topRadiusFilter(tree,filtArgs):
    dphill = -1
    if len(filtArgs)==0 : return 1.0,dphill
    if filtArgs[0] is None : return 1.0,dphill

    #select the b-jets and leptons at generator level
    bjets=[]
    leptons=[]
    lVec = ROOT.Math.LorentzVector(ROOT.Math.PtEtaPhiM4D('double'))
    for i in xrange(0,tree.nt):
        p4=lVec(tree.t_pt[i],tree.t_eta[i],tree.t_phi[i],tree.t_m[i])
        if abs(tree.t_id[i])==5 : bjets.append(p4)
        if abs(tree.t_id[i])==11 or abs(tree.t_id[i])==13: leptons.append(p4)

    #require at least two b-jets and two leptons    

    #angle between bb
    if len(bjets)<2 : return 1.0 ,dphill
    dphibb=ROOT.Math.VectorUtil.DeltaPhi(bjets[0],bjets[1])
    xbin=filtArgs[0].GetXaxis().FindBin(dphibb)

    #angle between ll in laboratory frame
    if len(leptons)<2 : return 1.0,dphill
    dphill=ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(leptons[0],leptons[1]))
    dphibb = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjets[0],bjets[1]))
    if  ROOT.TMath.Pi()/5 < dphibb < 2*ROOT.TMath.Pi()/5: dphill+=ROOT.TMath.Pi()
    elif  2*ROOT.TMath.Pi()/5 < dphibb < 3*ROOT.TMath.Pi()/5: dphill+=(2*ROOT.TMath.Pi())
    elif  3*ROOT.TMath.Pi()/5 < dphibb < 4*ROOT.TMath.Pi()/5: dphill+=(3*ROOT.TMath.Pi())
    elif  4*ROOT.TMath.Pi()/5 < dphibb: dphill+=(4*ROOT.TMath.Pi())
    xbin=filtArgs[0].GetXaxis().FindBin(dphill)

    return filtArgs[0].GetBinContent(xbin),dphill


"""
Analysis loop
"""
def runAnomalousTopProductionAnalysis2(fileName,outFileName,filterName):
     
    print '....analysing',fileName,'with output @',outFileName

    systVarNames=['nominal','btagup','btagdn','jerup','jerdn']
    #book histograms
    observablesH={}
    for k in ['emu','ll']:
        observablesH['ht_'+k]=ROOT.TH1F('ht_'+k,';H_{T} [GeV];Events',20,0,1000)

        observablesH['dphibb_'+k]=ROOT.TH1F('dphibb_'+k,';#Delta#phi(b,#bar{b}) [rad];Events',20,0,3.15)
        observablesH['cosbstar_'+k]=ROOT.TH1F('cosbstar_'+k,';cos(#theta*_{b});Events',20,-1,1)
        observablesH['cosbstarprod_'+k]=ROOT.TH1F('cosbstarprod_'+k,';cos(#theta*_{b_{1}})cos(#theta*_{b_{2}});Events',20,-1,1)

        observablesH['dphill_'+k]=ROOT.TH1F('dphill_'+k,';#Delta#phi(l,l) [rad];Events',20,0,3.15)
        observablesH['coslstar_'+k]=ROOT.TH1F('coslstar_'+k,';cos(#theta*_{l});Events',20,-1,1)
        observablesH['coslstarprod_'+k]=ROOT.TH1F('coslstarprod_'+k,';cos(#theta*_{l_{1}})cos(#theta*_{l_{2}});Events',20,-1,1)

	#Top+Top_ rest frame

        observablesH['etattbar_'+k]=ROOT.TH1F('etattbar_'+k,';#eta_{t#xoverline{t}};Events',20,-10,10)
        observablesH['cll**_'+k]=ROOT.TH1F('cll**_'+k,';cos(#theta_{l_{1}l_{2}}^{**});Events',20,-1,1)
        observablesH['cbb**_'+k]=ROOT.TH1F('cbb**_'+k,';cos(#theta_{b_{1}b_{2}}^{**});Events',20,-1,1)
        observablesH['dphill**_'+k]=ROOT.TH1F('dphill**_'+k,';#phi_{l_{1}l_{2}};Events',20,0,3.1415)
        observablesH['dphibb**_'+k]=ROOT.TH1F('dphibb**_'+k,';#phi_{b_{1}b_{2}};Events',20,0,3.1415)
        observablesH['dphillH_'+k]=ROOT.TH1F('dphillH_'+k,';#Delta#phi_{ll};Events',50,0,5*ROOT.TMath.Pi())
        observablesH['dphijj_'+k]=ROOT.TH1F('dphijj_'+k,';#Delta#phi(j,j) [rad];Events',20,0,3.15)
        observablesH['phiRphiG_'+k]=ROOT.TH2F('phiRphiG_'+k,';#Delta #phi^{Gen}_{ll};#Delta #phi^{Rec}_{ll}',50,0,5*ROOT.TMath.Pi(),50,0,5*ROOT.TMath.Pi())

        observablesH['dphillHshapes_exp_'+k]=ROOT.TH2F('dphillHshapes_exp_'+k,';#Delta#phi_{ll};Systematic variation;Events',50,0,5*ROOT.TMath.Pi(),len(systVarNames),0,len(systVarNames))
        for ybin in xrange(1,observablesH['dphillHshapes_exp_'+k].GetYaxis().GetNbins()+1):
            observablesH['dphillHshapes_exp_'+k].GetYaxis().SetBinLabel(ybin,systVarNames[ybin-1])


	


	#------------------
    for var in observablesH:
        observablesH[var].SetDirectory(0)
        observablesH[var].Sumw2()



    #open file
    puNormSF=1.0
    if 'MC13TeV' in fileName:
        fIn=ROOT.TFile.Open(fileName)
        puCorrH=fIn.Get('puwgtctr')
        nonWgt=puCorrH.GetBinContent(1)
        wgt=puCorrH.GetBinContent(2)
        if wgt>0 : puNormSF=nonWgt/wgt
        fIn.Close()
    tree=ROOT.TChain('twev')
    tree.AddFile(fileName)

    #if some filtering has to be applied, specialize it
    filtFunc,filtArgs=None,None
    if filterName:
        filtFunc,filtArgsList=filterName.split('=')
        if filtFunc=='topRadiusFilter':
            args=filtArgsList.split(',')
            bsmUrl,smUrl,distName=args

            #get the target distribution
            bsmFile=ROOT.TFile.Open(bsmUrl)
            weightH=bsmFile.Get(distName)
            weightH.SetDirectory(0)
            bsmFile.Close()

            #get the SM distribution 
            smFile=ROOT.TFile.Open(smUrl)
            weightH.Divide(smFile.Get(distName))
            smFile.Close()
            filtArgs=[weightH]
        else:
            filtArgs=filtArgsList.split(',')
                                   

    #loop over events in the tree and fill histos
    totalEntries=tree.GetEntries()
    lVec = ROOT.Math.LorentzVector(ROOT.Math.PtEtaPhiM4D('double')) 
    for i in xrange(0,totalEntries):

        tree.GetEntry(i)

        if i%100==0 : sys.stdout.write('\r [ %d/100 ] done' %(int(float(100.*i)/float(totalEntries))) )

        filtWeight=1.0
	dphill_gen=-1
        if filtFunc and filtArgs:
            filtWeight,dphill_gen=globals()[filtFunc](tree,filtArgs)
            if filtWeight<0 : continue
            

        if abs(tree.cat)<100 : continue
        evcat = 'emu' if abs(tree.cat)==11*13 else 'll'

        evWeight=puNormSF*tree.weight[0]*filtWeight

        #leptons
        leptons=[]
        for il in xrange(0,tree.nl):
            leptons.append( lVec(tree.l_pt[il],tree.l_eta[il],tree.l_phi[il],tree.l_m[il]) )
        if len(leptons)<2 : continue
	#JER
	#bjetsJesUp, otherjetsJesUp = [],[]
	#bjetsJesDn, otherjetsJesDn= [],[] 

	#for ik in xrange(0,tree.nj):

		#jesSF1=(tree.j_jes[ik] >> 1)
		#jesSF2=(tree.j_jes[ik] & 0x1)
		#print jesSF1

		#jesUpP4 = jesSF[1]*p4
		#jesDnP4 = jesSF[2]*p4
		#jerUpP4 = jerSF[1]*p4
		#jerDnP4 = jerSF[2]*p4
	#print tree.nj		

		
        #preselect the b-jets (save always the jet and the gen jet)
        bjets,otherjets=[],[]
        bjetsUp,otherjetsUp=[],[]
        bjetsDn,otherjetsDn=[],[]
	bjetsJerUp, otherjetsJerUp= [],[]
	bjetsJerDn, otherjetsJerDn= [],[]
        for ij in xrange(0,tree.nj):
	    jetP4=lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij])
	    jerFactor=ROOT.TMath.Abs(1-tree.j_jer[ij])
	    jetP4Up=jetP4*(1+jerFactor)
	    jetP4Dn=jetP4*(1-jerFactor)
	    
            #nominal b-tag decision
            btagVal=(tree.j_btag[ij] & 0x1)
            if btagVal!=0 and len(bjets)<2: 
                bjets.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
		if jetP4Up.pt()>30  and ROOT.TMath.Abs(jetP4Up.eta())<2.5 :
			bjetsJerUp.append(jetP4Up)
		if jetP4Dn.pt()>30  and ROOT.TMath.Abs(jetP4Dn.eta())<2.5 :
			bjetsJerDn.append(jetP4Dn)
            else:
                otherjets.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
		otherjetsJerUp.append(jetP4Up)
		otherjetsJerDn.append(jetP4Dn)

            #up varied b-tag decision
            btagValUp=((tree.j_btag[ij] >> 1) & 0x1)
            if btagValUp!=0 and len(bjetsUp)<2: 
                bjetsUp.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
            else:
                otherjetsUp.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )

            #down varied b-tag decision
            btagValDn=((tree.j_btag[ij] >> 2) & 0x1)
            if btagValDn!=0 and len(bjetsDn)<2: 
                bjetsDn.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )
            else:
                otherjetsDn.append( lVec(tree.j_pt[ij],tree.j_eta[ij],tree.j_phi[ij],tree.j_m[ij]) )

        #loop over syst variations
        for iSystVar in xrange(0,5):

            dphill = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(leptons[0],leptons[1]))
	    dphillJerUp = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(leptons[0],leptons[1]))
	    dphillJerDn = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(leptons[0],leptons[1]))
            dphibb=None
	    dphibbJerUp=None
	    dphibbJerDn=None
            if iSystVar==0 :
		if len(bjets)!=2 : continue		
               	dphibb = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjets[0],bjets[1]))
            if iSystVar==3 :
		if len(bjetsJerUp)!=2 : continue
		dphibbJerUp = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjetsJerUp[0],bjetsJerUp[1]))
            if iSystVar==4 :
		if len(bjetsJerDn)!=2 : 
			continue
		dphibbJerDn = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjetsJerDn[0],bjetsJerDn[1]))
		print("1")
            if iSystVar==1:
                if len(bjetsUp)!=2 : continue
                dphibb = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjetsUp[0],bjetsUp[1]))
            if iSystVar==2:
                if len(bjetsDn)!=2 : continue
                dphibb = ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjetsDn[0],bjetsDn[1]))
          

            if dphibb < 2*ROOT.TMath.Pi()/5 and dphibb >= ROOT.TMath.Pi()/5:     dphill+=ROOT.TMath.Pi()
            elif dphibb < ROOT.TMath.Pi()*3/5 and dphibb >= ROOT.TMath.Pi()*2/5: dphill+=2*ROOT.TMath.Pi()
            elif dphibb < ROOT.TMath.Pi()*4/5 and dphibb >= ROOT.TMath.Pi()*3/5: dphill+=3*ROOT.TMath.Pi()
            elif dphibb <= ROOT.TMath.Pi() and dphibb >= ROOT.TMath.Pi()*4/5:    dphill+=4*ROOT.TMath.Pi()
	       
            if dphibbJerUp < 2*ROOT.TMath.Pi()/5 and dphibbJerUp >= ROOT.TMath.Pi()/5:     dphillJerUp+=ROOT.TMath.Pi()
            elif dphibbJerUp < ROOT.TMath.Pi()*3/5 and dphibbJerUp >= ROOT.TMath.Pi()*2/5: dphillJerUp+=2*ROOT.TMath.Pi()
            elif dphibbJerUp < ROOT.TMath.Pi()*4/5 and dphibbJerUp >= ROOT.TMath.Pi()*3/5: dphillJerUp+=3*ROOT.TMath.Pi()
            elif dphibbJerUp <= ROOT.TMath.Pi() and dphibbJerUp >= ROOT.TMath.Pi()*4/5:    dphillJerUp+=4*ROOT.TMath.Pi()
	

            if dphibbJerDn < 2*ROOT.TMath.Pi()/5 and dphibbJerDn >= ROOT.TMath.Pi()/5:     dphillJerDn+=ROOT.TMath.Pi()
            elif dphibbJerDn < ROOT.TMath.Pi()*3/5 and dphibbJerDn >= ROOT.TMath.Pi()*2/5: dphillJerDn+=2*ROOT.TMath.Pi()
            elif dphibbJerDn < ROOT.TMath.Pi()*4/5 and dphibbJerDn >= ROOT.TMath.Pi()*3/5: dphillJerDn+=3*ROOT.TMath.Pi()
            elif dphibbJerDn <= ROOT.TMath.Pi() and dphibbJerDn >= ROOT.TMath.Pi()*4/5:    dphillJerDn+=4*ROOT.TMath.Pi()
            
	    if iSystVar == 0 or iSystVar == 1 or iSystVar ==2 :
		   observablesH['dphillHshapes_exp_'+evcat].Fill(dphill,iSystVar,evWeight)
	    if iSystVar == 3 :
		   observablesH['dphillHshapes_exp_'+evcat].Fill(dphillJerUp,iSystVar,evWeight)
	    if iSystVar == 4 :
		   observablesH['dphillHshapes_exp_'+evcat].Fill(dphillJerDn,iSystVar,evWeight)

            if iSystVar==0:
                observablesH['dphillH_'+evcat].Fill(dphill,evWeight)
                observablesH['phiRphiG_'+evcat].Fill(dphill_gen,dphill,evWeight)  

        if len(bjets)<2: continue
        #met
        metx,mety=tree.met_pt*ROOT.TMath.Cos(tree.met_phi),tree.met_pt*ROOT.TMath.Sin(tree.met_phi)

        #try to solve the kinematics (need to swap bl assignments)
        allSols=[]
        try:
            sols=doubleNeutrinoSolutions( (bjets[0],   bjets[1]), 
                                          (leptons[0], leptons[1]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):               
                top=bjets[0]+leptons[0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1]+leptons[1]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][1],0.)
                allSols.append( (0,top,top_) )
        except numpy.linalg.linalg.LinAlgError:
            pass        
        try:
            sols=doubleNeutrinoSolutions( (bjets[0],   bjets[1]), 
                                          (leptons[1], leptons[0]),
                                          (metx,mety) )
            for isol in xrange(0,len(sols.nunu_s)):
                top=bjets[0]+leptons[1]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][0],0.)
                top_=bjets[1]+leptons[0]+convertToPtEtaPhiM(lVec,sols.nunu_s[isol][1],0.)
                allSols.append( (1,top,top_) )
        except numpy.linalg.linalg.LinAlgError :
            pass

        #sort solutions by increasing m(ttbar)
        if len(allSols)==0: continue
        allSols=sorted(allSols, key=lambda sol: (sol[1]+sol[2]).mass() )        
        l1idx=0 if allSols[0][0]==0 else 1
        l2idx=1 if allSols[0][0]==0 else 0

        #setup the Lorentz transformations to the top/anti-top rest frames
        topBoost, top_Boost = allSols[0][1].BoostToCM(), allSols[0][2].BoostToCM()

        #measure b-jet angles
        cosb1 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(bjets[0],topBoost), allSols[0][1] ) )
        cosb2 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(bjets[1],top_Boost), allSols[0][2] ) )
        observablesH['dphibb_'+evcat].Fill(ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(bjets[0],bjets[1])),evWeight)
        observablesH['cosbstar_'+evcat].Fill(cosb1,evWeight)
        observablesH['cosbstar_'+evcat].Fill(cosb2,evWeight)
        observablesH['cosbstarprod_'+evcat].Fill(cosb1*cosb2,evWeight)

        #measure leptonic angles
        cosl1 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(leptons[l1idx],topBoost), allSols[0][1] ) )
        cosl2 = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(leptons[l2idx],top_Boost), allSols[0][2] ) )
        observablesH['dphill_'+evcat].Fill(ROOT.TMath.Abs(ROOT.Math.VectorUtil.DeltaPhi(leptons[l1idx],leptons[l2idx])),evWeight)
        observablesH['coslstar_'+evcat].Fill(cosl1,evWeight)
        observablesH['coslstar_'+evcat].Fill(cosl2,evWeight)
        observablesH['coslstarprod_'+evcat].Fill(cosl1*cosl2,evWeight)

        if len(otherjets)>=2:
            observablesH['dphijj_'+evcat].Fill(ROOT.Math.VectorUtil.DeltaPhi(otherjets[0],otherjets[1]),evWeight)
         
        #other control variables
        observablesH['ht_'+evcat].Fill(bjets[0].pt()+bjets[1].pt()+leptons[0].pt()+leptons[1].pt()+tree.met_pt,evWeight)

	#Top + Top_ reference frame

	ttbarBoost = (allSols[0][1]+allSols[0][2]).BoostToCM()
	cll_tt = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle( ROOT.Math.VectorUtil.boost(leptons[l1idx],ttbarBoost), ROOT.Math.VectorUtil.boost(leptons[l2idx],ttbarBoost)))
	cbb_tt = ROOT.TMath.Cos( ROOT.Math.VectorUtil.Angle(ROOT.Math.VectorUtil.boost(bjets[0],ttbarBoost),ROOT.Math.VectorUtil.boost(bjets[1],ttbarBoost)))
	
	lepton1_tt = ROOT.Math.VectorUtil.boost(leptons[l1idx],ttbarBoost)
	lepton2_tt = ROOT.Math.VectorUtil.boost(leptons[l2idx],ttbarBoost)
	bjets1_tt = ROOT.Math.VectorUtil.boost(bjets[0],ttbarBoost)
	bjets2_tt = ROOT.Math.VectorUtil.boost(bjets[1],ttbarBoost)

	dphill_tt = ROOT.Math.VectorUtil.DeltaPhi(lepton1_tt,lepton2_tt)
	dphibb_tt =  ROOT.Math.VectorUtil.DeltaPhi( bjets1_tt,bjets2_tt)
        observablesH['cll**_'+evcat].Fill(cll_tt,evWeight)
        observablesH['cbb**_'+evcat].Fill(cbb_tt,evWeight)
        observablesH['dphill**_'+evcat].Fill(dphill_tt,evWeight)
        observablesH['dphibb**_'+evcat].Fill(dphibb_tt,evWeight)
	
	etattbar = allSols[0][1].eta()-allSols[0][2].eta()
	observablesH['etattbar_'+evcat].Fill(etattbar,evWeight)	





    #save results

    if filterName:
        postfix=outFileName.split('_')[-1]
        outFileName=outFileName.replace(postfix,'%s_%s'%(filterName,postfix))
        outFileName=outFileName.replace('=','_')
        outFileName=outFileName.replace(',','_')
    fOut=ROOT.TFile.Open(outFileName,'RECREATE')
    for var in observablesH: 
        observablesH[var].Write()
    fOut.Close()

 
"""
Wrapper for when the analysis is run in parallel
"""
def runAnomalousTopProductionAnalysisPacked(args):
    try:
        fileNames,outFileName,filterName=args
        runAnomalousTopProductionAnalysis2(fileNames,outFileName,filterName)
    except : # ReferenceError:
        print 50*'<'
        print "  Problem with", name, "continuing without"
        print 50*'<'
        return False
    
"""
Create analysis tasks
"""
def createAnalysisTasks(opt):

    onlyList=opt.only.split('v')

    ## Local directory
    file_list=[]
    if os.path.isdir(opt.input):
        for file_path in os.listdir(opt.input):
            if file_path.endswith('.root'):
                file_list.append(os.path.join(opt.input,file_path))
    elif opt.input.startswith('/store/'):
        file_list = getEOSlslist(opt.input)
    elif '.root' in opt.input:
        file_list.append(opt.input)

    #list of files to analyse
    tasklist=[]
    for filename in file_list:
        baseFileName=os.path.basename(filename)      
        tag,ext=os.path.splitext(baseFileName)
        if len(onlyList)>0:
            processThis=False
            for filtTag in onlyList:
                if filtTag in tag:
                    processThis=True
            if not processThis : continue
        tasklist.append((filename,'%s/%s'%(opt.output,baseFileName),opt.filter))

    #loop over tasks
    if opt.queue=='local':
        if opt.jobs>1:
            print ' Submitting jobs in %d threads' % opt.jobs
            import multiprocessing as MP
            pool = MP.Pool(opt.jobs)
            pool.map(runAnomalousTopProductionAnalysisPacked,tasklist)
        else:
            for fileName,outFileName,filterName in tasklist:
                runAnomalousTopProductionAnalysis2(fileName,outFileName,filterName)
    else:
        cmsswBase=os.environ['CMSSW_BASE']
        for fileName,_,filterName in tasklist:
            localRun='python %s/src/TopLJets2015/TopAnalysis/scripts/runAnomalousTopProductionAnalysis2.py -i %s -o %s -q local'%(cmsswBase,fileName,opt.output)
            if filterName : localRun+=' --filter %s'%filterName
            cmd='bsub -q %s %s/src/TopLJets2015/TopAnalysis/scripts/wrapLocalAnalysisRun.sh \"%s\"' % (opt.queue,cmsswBase,localRun)
            print cmd
            os.system(cmd)

"""
steer
"""
def main():
	usage = 'usage: %prog [options]'
	parser = optparse.OptionParser(usage)
	parser.add_option('-i', '--input',
                          dest='input',   
                          default='/afs/cern.ch/user/p/psilva/work/Stop',
                          help='input directory with the files [default: %default]')
	parser.add_option('-f', '--filter',
                          dest='filter',   
                          default=None,
                          help='apply this filter function')
	parser.add_option('--jobs',
                          dest='jobs', 
                          default=1,
                          type=int,
                          help='# of jobs to process in parallel the trees [default: %default]')
	parser.add_option('--only',
                          dest='only', 
                          default='',
                          type='string',
                          help='csv list of tags to process')
	parser.add_option('-o', '--output',
                          dest='output', 
                          default='analysis',
                          help='Output directory [default: %default]')
	parser.add_option('-q', '--queue',
                          dest='queue',
                          default='local',
                          help='Batch queue to use [default: %default]')
	(opt, args) = parser.parse_args()

        ROOT.FWLiteEnabler.enable() 
	os.system('mkdir -p %s' % opt.output)

        createAnalysisTasks(opt)
        

if __name__ == "__main__":
	sys.exit(main())

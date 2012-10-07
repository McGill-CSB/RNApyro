#/usr/bin/python
import os, sys, re, os.path, random
import matplotlib
import matplotlib.pyplot as plt;
from mpl_toolkits.axes_grid1 import make_axes_locatable;
import numpy as np;
from pylab import *;

## usage

def usage(softname):
    print "%s [-cvswx] -i <RNApyro infile> -o <RNApyro outfile>" % (softname)
    sys.exit(1)

## compute bp set

def ssa2bp(ssa):
    
    stack=[]
    bpdic={}
    
    for i in range(len(ssa)):
        if ssa[i]=='(':
            stack.append(i)
        elif ssa[i]==')':
            j=stack.pop()
            bpdic[i]=j
            bpdic[j]=i
    
    if len(stack)>0:
        print 'Corrupted secondary structure'
        sys.exit(1)
    
    return bpdic

## compute boundaries

def subseqindex(fullseq):
    
    iopen=0
    while fullseq[iopen]=='-' or fullseq[iopen]=='.':
        iopen+=1
    
    iclose=len(fullseq)-1
    while fullseq[iclose]=='-' or fullseq[iclose]=='.':
        iclose-=1
    
    return iopen,iclose


## match sequence and structure

def fitssa2seq(fullssa,fullseq,includeextbp):
    iopen,iclose=subseqindex(fullseq)
    bpdic=ssa2bp(fullssa)
    newseq=''
    newssa=''
    newconsensus=''
    outbp=[]
    #    for i in range(len(fullseq)):
    for i in range(iopen,iclose+1,1):
        if fullseq[i]!='-':
            if fullssa[i]=='(' or fullssa[i]==')':
                j=bpdic[i]
                if fullssa[j]=='(' or fullssa[j]==')':
                    newseq+=fullseq[i]
                    if i<j:
                        if j>=iopen and j<=iclose:
                            newssa+='('
                        else:
                            if includeextbp:
                                newssa+='['
                            else:
                                newssa+='.'
                    else:
                        if j>=iopen and j<=iclose:
                            newssa+=')'
                        else:
                            if includeextbp:
                                newssa+=']'
                            else:
                                newssa+='.'
                else:
                    print 'Corrupted ssa at index (%d,%d)' % (i,j)
                    sys.exit(1)
            elif fullssa[i]=='.':
                newseq+=fullseq[i]
                newssa+='.'
            else:
                print 'Corrupted ssa at index %d' % (i)
                sys.exit(1)
        else:
            # gap at index i
            if fullssa[i]=='(' or fullssa[i]==')':
                j=bpdic[i]
                if fullseq[j]!='-':
                    if fullssa[j]=='(' or fullssa[j]==')':
                        if i>iopen and i<iclose:
                            newseq+='N'
                            if i<j:
                                newssa+='('
                            else:
                                newssa+=')'
                        else:
                            outbp.append(i)
                    else:
                        print 'Corrupted ssa at index (%d,%d)' % (i,j)
                        sys.exit(1)
            elif fullssa[i]=='.' and i>iopen and i<iclose:
                newseq+='-'
                newssa+='.'
    
    
    #print fullseq
    #print fullssa
    #print newseq
    #print newssa
    
    return newseq,newssa

## compute mutations

def hammingdistance(data):

    read=data['read'].replace('.','').replace('-','').upper()
    original=data['target'].replace('.','').replace('-','').upper()

    if len(read)!=len(original):
        print 'ERROR'
        sys.exit(1)

    mutations=[]
    for i in range(len(read)):
        if read[i]!=original[i]:
            mutations.append((original[i],i,read[i]))
    
    return mutations


## read matrix

def readoutput(filename):

    data_re = re.compile("(?P<index>\d+)(\s+)(?P<probA>[0-9\.e\-\+]+)(\s+)(?P<probC>[0-9\.e\-\+]+)(\s+)(?P<probG>[0-9\.e\-\+]+)(\s+)(?P<probU>[0-9\.e\-\+]+)\n")
    
    f=open(filename,'r')
    data={} 
    for buffer in f.readlines():
        cell = data_re.match(buffer)
        if cell:
            index = int(cell.group('index'))
            if data.has_key(index):
                print 'ERROR'
                sys.exit(1)
            data[index]={}
            data[index]['A'] = float(cell.group('probA'))
            data[index]['C'] = float(cell.group('probC'))
            data[index]['G'] = float(cell.group('probG'))
            data[index]['U'] = float(cell.group('probU'))
        else:
            print 'skip'
            print buffer
            sys.exit(1)
                
    f.close()

    return data

## read infile

def readinput(filename):
    
    f=open(filename,'r')
    
    data={}
    mode='none'
    for buffer in f.readlines():
        if buffer.startswith('>'):
            if mode!='none':
                print 'ERROR: file corrupted'
                sys.exit(1)            
            if buffer.endswith('read\n'):
                mode='read'
            elif buffer.endswith('structure\n'):
                mode='structure'
            elif buffer.endswith('consensus\n'):
                mode='consensus'
            elif buffer.endswith('target\n'):
                mode='target'
            else:   
                print 'ERROR: file corrupted'
                sys.exit(1)            
        else:    
            if mode=='read':
                data['read']=buffer.replace('\n','')
                mode='none'
            elif mode=='structure':
                data['ssa']=buffer.replace('\n','')
                mode='none'
            elif mode=='consensus':
                data['consensus']=buffer.replace('\n','')
                mode='none'
            elif mode=='target':
                data['target']=buffer.replace('\n','')
                mode='none'
            else:   
                print 'ERROR: file corrupted'
                sys.exit(1)            

    f.close()

    return data

## read infile

def randommutants(data,nmut):
    mutationindex=[]
    mutationlist=[]
    original=data['target']
    for k in range(nmut):

        i=random.randint(0,len(original)-1)
        while i in mutationindex:
            i=random.randint(0,len(original)-1)
        mutationindex.append(i)
                
        mymut=original[i]
        while mymut==original[i]:
            mymut=random.choice('ACGU')
            
        mutationlist.append((original[i],i,mymut))

    return mutationlist

    
## correlation 2

def probprofile(indatafull,outdatafull,withread):
    
    scores={}
    counter={}
    
    for alpha in indatafull.keys():

        indata = indatafull[alpha]
        outdata = outdatafull[alpha]
        scores[alpha]={'TP':0.0,'TN':0.0,'FP':0.0,'FN':0.0}
        counter[alpha]={'TP':0,'TN':0,'FP':0,'FN':0}

        for mkey in indatafull[alpha].keys():
            read=indata[mkey]['read'].replace('.','').replace('-','').upper()
            original=indata[mkey]['target'].replace('.','').replace('-','').upper()
            
            if len(read)!=len(original):
                print 'ERROR'
                sys.exit(1)
            
            for i in range(len(read)):
                #print i,(outdata[i]['A']+outdata[i]['C']+outdata[i]['G']+outdata[i]['U'])
                for nt in ['A','C','G','U']:
                    if nt!=read[i] or withread:
                        if nt==original[i]:
                            # nt is the good one (T)
                            if original[i]!=read[i]: # this is a mutation site
                                scores[alpha]['TP']+=outdata[mkey][i][nt]
                                counter[alpha]['TP']+=1
                            else: # nt is NOT a mutation site (N)
                                scores[alpha]['TN']+=outdata[mkey][i][nt]
                                counter[alpha]['TN']+=1
                        else:
                            # nt is not correct (F)
                            if original[i]!=read[i]: # nt is a mutation site(P)
                                scores[alpha]['FP']+=outdata[mkey][i][nt]
                                counter[alpha]['FP']+=1
                            else:
                                # nt is NOT a mutation site (N)
                                scores[alpha]['FN']+=outdata[mkey][i][nt]
                                counter[alpha]['FN']+=1
    
    return scores,counter



## compute ROC curve and AUC

def makeROC(indatafull,outdatafull,onlybpregion,verbose,customsteps,withread,withfigs,rocfile):
    
    output={}
    xcoord={}
    ycoord={}
    
    for alpha in indatafull.keys():
    
        indata = indatafull[alpha]
        outdata = outdatafull[alpha]
    
        sensitivity={}
        specificity={}
        PPV={}
        
        bestFM=0
        bestsens=0
        bestspec=0
        bestT=0
        startrange=0
        endrange=customsteps
        mystep=1
        myrange=float(customsteps)
            
        for myvalue in range(startrange,endrange+1,mystep):
            TP=0
            TN=0
            FP=0
            FN=0
            threshold=float(myvalue)/myrange
            for mkey in indata.keys():
                read=indata[mkey]['read'].replace('.','').replace('-','').upper()
                original=indata[mkey]['target'].replace('.','').replace('-','').upper()
                myseq,myssa=fitssa2seq(indata[mkey]['ssa'],indata[mkey]['target'],False)

                ## FIXME: filter mutations in bp regions
                
                if len(read)!=len(original):
                    print 'ERROR: length do not match.'
                    sys.exit(1)
                
                for i in range(len(read)):
                    if myssa[i]!='.' or onlybpregion: # mutation in a base pair
                        for nt in ['A','C','G','U']:
                            if withread or nt!=read[i]:
                                if outdata[mkey][i][nt]>=threshold: # predicted (P)
                                    if nt==original[i]: # nt is correct (T)
                                        TP+=1
                                    else: # nt is wrong (F)
                                        FP+=1
                                else: # not predicted (N)
                                    if nt!=original[i]: # nt is NOT the good nt (T)
                                        TN+=1
                                    else: # nt is the good one (F)
                                        FN+=1
                        
                        
            #print TP,TN,FP,FN
            if TP+FN>0:
                mysensitivity=float(TP)/(TP+FN)
            else:
                mysensitivity=0
            if TN+FP>0:
                myspecificity=float(TN)/(TN+FP)
            else:
                myspecificity=0
            if TP+FP>0:
                myPPV=float(TP)/(TP+FP)
            else:
                myPPV=1.0
            if mysensitivity+myPPV>0:
                Fmeasure= 2.0 * (mysensitivity*myPPV)/(mysensitivity+myPPV)
            else:
                Fmeasure=0
            if verbose:
                print "%.2f: %.2f (%.2f,%.2f,%.2f)" % (threshold,Fmeasure,mysensitivity,myPPV,myspecificity)
            if Fmeasure>bestFM:
                bestFM=Fmeasure
                bestsens=mysensitivity
                bestspec=myspecificity
                bestT=threshold
            sensitivity[myvalue]=mysensitivity
            specificity[myvalue]=myspecificity
            PPV[myvalue]=myPPV

        # compute ROC
        prevsens=0
        prevspec=0
        myroc=0.0
        xcoord[alpha]=[0.0]
        ycoord[alpha]=[0.0]
        for myvalue in range(endrange,startrange-1,-mystep):
        #for myvalue in range(startrange,endrange+1,mystep):
            mysens=sensitivity[myvalue]
            myspec=1.0 - specificity[myvalue]
            myarea=(myspec-prevspec)*(mysens+prevsens)/2.0
            myroc+=myarea
            prevsens=mysens
            prevspec=myspec
            xcoord[alpha].append(myspec)
            ycoord[alpha].append(mysens)
        myroc+=(1.0-prevspec)*(1.0+prevsens)/2.0
        output[alpha] = {'FM':bestFM,'ROC':myroc}
    
    #print xcoord,ycoord
    # make figure
    if withfigs:
        fig = plt.figure();
        ax = fig.add_subplot(111)
        ax.set_xlabel('Specificity',size='xx-large')
        ax.set_ylabel('Sensitivity',size='xx-large')
        
        ax.plot(xcoord[100],ycoord[100],'r.-',label=r'$\alpha=1.0$, $AUC=%.2f$'% output[100]['ROC'])
        ax.plot(xcoord[80],ycoord[80],'c.-',label=r'$\alpha=0.8$, $AUC=%.2f$'% output[80]['ROC'])
        ax.plot(xcoord[50],ycoord[50],'g.-',label=r'$\alpha=0.5$, $AUC=%.2f$'% output[50]['ROC'])
        ax.plot(xcoord[0],ycoord[0],'m.-',label=r'$\alpha=0.0$, $AUC=%.2f$'% output[0]['ROC'])
        ax.plot([0.0,1.0],[0.0,1.0],'b--')
        ax.legend(loc=4, prop={'size':'x-large'})
        fig.savefig(rocfile,transparent=True)
        print 'ROC curve saved in file',rocfile

    return output


## main

def main(indir,infileprefix,onlybpregion,verbose,mysteps,withread,withfigs):
    
    indata = {}
    outdata = {}
    #nmut_list = [6,12,24]
    nmut_list = [6,12]
    alpha_list = [100,80,50,0]
    
    for nmut in nmut_list:
        indata[nmut]={}
        outdata[nmut]={}
        for alpha in alpha_list:
            indata[nmut][alpha]={}
            outdata[nmut][alpha]={}
            for mkey in range(45):
                infilename = infileprefix + '_' + str(mkey) + '.ref'
                outfilename = indir + infileprefix + '_' + str(mkey) + '_nmut_' + str(nmut) + '_alpha_' + str(alpha) + '.out'
                #print infilename
                if os.path.exists(infilename) and os.path.exists(outfilename):
                    indata[nmut][alpha][mkey]=readinput(infilename)
                    outdata[nmut][alpha][mkey]=readoutput(outfilename)
                    print infilename,'and',outfilename,'stored'


    print ">> Mutational Profile"
    for nmut in nmut_list:
        scores,counter = probprofile(indata[nmut],outdata[nmut],withread)
        for alpha in alpha_list:
            val = {'TP':0.0,'TN':0.0,'FP':0.0,'FN':0.0}
            for mkey in val.keys():
                if counter[alpha][mkey]>0:
                    #                    val[mkey] = scores[alpha][mkey]/counter[alpha][mkey]
                    val[mkey] = scores[alpha][mkey]/119
            print 'Mutations: %d; alpha: %.1f; Good=%.2f; Bad=%.2f' % (nmut,float(alpha)/100,val['TP']+val['TN'],val['FP']+val['FN'])
        
    
    ## Accuracy of prediction using threshold
    print ">> ROC curve"
    for nmut in nmut_list:
        rocfile="ROC_" + str(nmut) + ".eps"
        scores = makeROC(indata[nmut],outdata[nmut],onlybpregion,verbose,mysteps,withread,withfigs,rocfile)
        for alpha in alpha_list:
            print 'Mutations: %d; Alpha: %.1f; Best F-measure: %.2f; ROC: %.2f' %(nmut,float(alpha)/100,scores[alpha]['FM'],scores[alpha]['ROC'])
    

###############################################################################

if __name__ == '__main__':

    onlybpregion=True
    verbose=False
    mysteps=100
    withread=True
    withfigs=False
    indir=''
    
    if len(sys.argv)==1:
        usage(sys.argv[0])

    for i in range(len(sys.argv)):
        myarg=sys.argv[i]
        if myarg=='-h' or myarg=='--help':
            usage(sys.argv[0])
        if myarg=='-c' or myarg=='--bp':
            onlybpregion=False
        if myarg=='-v' or myarg=='--verbose':
            verbose=True
        if myarg=='-i':
            inputfileprefix=sys.argv[i+1]
        if myarg=='-o':
            outputfile=sys.argv[i+1]
        if myarg=='-s':
            mysteps=int(sys.argv[i+1])
        if myarg=='-w':
            withread=False
        if myarg=='-x':
            withfigs=True
        if myarg=='-d':
            indir=sys.argv[i+1] + '/'

    main(indir,inputfileprefix,onlybpregion,verbose,mysteps,withread,withfigs)
    










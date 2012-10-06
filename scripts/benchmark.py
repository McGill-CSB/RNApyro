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

## read infile

def basiccorrelation(mutlist,data):

    score=0
    for nt1,index,nt2 in mutlist:
        score += data[index][nt1]

    if False:
        sumscore=0
        for index,probs in data.iteritems():
            subtotal = probs['A']+probs['C']+probs['G']+probs['U']
            if False and subtotal!=1:
                print index,probs
            sumscore+=subtotal
        print sumscore
    
    return score/len(mutlist)

## read infile

def fullcorrelation(indata,outdata,withread):
    
    read=indata['read'].replace('.','').replace('-','').upper()
    original=indata['target'].replace('.','').replace('-','').upper()
    
    if len(read)!=len(original):
        print 'ERROR'
        sys.exit(1)
    
    goodscore=0.
    badscore=0.
    cmptgood=0
    cmptbad=0
    scores={'TP':0.0,'TN':0.0,'FP':0.0,'FN':0.0}
    counter={'TP':0,'TN':0,'FP':0,'FN':0}
    for i in range(len(read)):
        #print i,(outdata[i]['A']+outdata[i]['C']+outdata[i]['G']+outdata[i]['U'])
        for nt in ['A','C','G','U']:
            if nt!=read[i] or withread:
                if nt==original[i]:
                    # nt is the good one (T)
                    if original[i]!=read[i]: # this is a mutation site
                        scores['TP']+=outdata[i][nt]
                        counter['TP']+=1
                    else: # nt is NOT a mutation site (N)
                        scores['TN']+=outdata[i][nt]
                        counter['TN']+=1
                else:
                    # nt is not correct (F)
                    if original[i]!=read[i]: # nt is a mutation site(P)
                        scores['FP']+=outdata[i][nt]
                        counter['FP']+=1
                    else:
                        # nt is NOT a mutation site (N)
                        scores['FN']+=outdata[i][nt]
                        counter['FN']+=1

    return scores,counter
    

## read infile

def checkconsistency(mutlist,indata,outdata):
    
    read=indata['read'].replace('.','').replace('-','').upper()
    original=indata['target'].replace('.','').replace('-','').upper()
    
    if len(read)!=len(original):
        print 'ERROR'
        sys.exit(1)
    
    goodscore=0.
    badscore=0.
    cmptgood=0
    cmptbad=0
    for i in range(len(read)):
        for nt in ['A','C','G','U']:
            if nt!=original[i]:
                if nt!=read[i]:
                    badscore+=outdata[i][nt]
                    cmptbad+=1
                else:
                    goodscore+=outdata[i][nt]
                    cmptgood+=1
    print goodscore,badscore,goodscore+badscore


## read infile

def fullprediction(mutlist,indata,outdata,onlybpregion,verbose,customsteps,withread,withfigs):
    
    read=indata['read'].replace('.','').replace('-','').upper()
    original=indata['target'].replace('.','').replace('-','').upper()
    myseq,myssa=fitssa2seq(indata['ssa'],indata['target'],False)

    ## FIXME: filter mutations in bp regions
    
    if len(read)!=len(original):
        print 'ERROR'
        sys.exit(1)
    
    
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
        for i in range(len(read)):
            if myssa[i]!='.' or onlybpregion: # mutation in a base pair
                for nt in ['A','C','G','U']:
                    if withread or nt!=read[i]:
                        if outdata[i][nt]>=threshold: # predicted (P)
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
    xcoord=[0.0]
    ycoord=[0.0]
    for myvalue in range(endrange,startrange-1,-mystep):
    #for myvalue in range(startrange,endrange+1,mystep):
        mysens=sensitivity[myvalue]
        myspec=1.0 - specificity[myvalue]
        myarea=(myspec-prevspec)*(mysens+prevsens)/2.0
        myroc+=myarea
        prevsens=mysens
        prevspec=myspec
        xcoord.append(myspec)
        ycoord.append(mysens)
    myroc+=(1.0-prevspec)*(1.0+prevsens)/2.0
    #print xcoord,ycoord
    # make figure
    if withfigs:
        fig = plt.figure();
        ax = fig.add_subplot(111)
        ax.plot(xcoord,ycoord,'ro-')
        ax.plot([0.0,1.0],[0.0,1.0],'b--')
        fig.savefig('roccurve.pdf',transparent=True)
        print 'ROC curve saved in file roccurve.pdf'

    return bestFM,myroc



## main

def main(infilename,outfilename,onlybpregion,verbose,mysteps,withread,withfigs,printfile):
    
    indata=readinput(infilename)
    outdata=readoutput(outfilename)

    if printfile!='':
        f = open(printfile,'a')
    
    mutationlist = hammingdistance(indata)

    print 'number of mutations: ', len(mutationlist)

    ## withreads means that proba of nucleotide in read is considered as well
    
    ## Correlation
    scores,counter = fullcorrelation(indata,outdata,withread)
    val = {'TP':0.0,'TN':0.0,'FP':0.0,'FN':0.0}
    for mkey in val.keys():
        if counter[mkey]>0:
            val[mkey] = scores[mkey]/counter[mkey]
    print 'Correlation: TP=%.2f; TN=%.2f; FP=%.2f; FN=%.2f' % (val['TP'],val['TN'],val['FP'],val['FN'])
    
    ## Accuracy of prediction using threshold
    myFM,myroc = fullprediction(mutationlist,indata,outdata,onlybpregion,verbose,mysteps,withread,withfigs)
    print 'Best F-measure: %.2f; ROC: %.2f' %(myFM,myroc)
    
    if printfile!='':
        print >>f,val['TP'],val['TN'],val['FP'],val['FN'],myFM,myroc
                
    if printfile!='':
        f.close()


###############################################################################

if __name__ == '__main__':

    onlybpregion=True
    verbose=False
    mysteps=100
    withread=True
    withfigs=False
    printfile=''
    
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
            inputfile=sys.argv[i+1]
        if myarg=='-o':
            outputfile=sys.argv[i+1]
        if myarg=='-s':
            mysteps=int(sys.argv[i+1])
        if myarg=='-w':
            withread=False
        if myarg=='-x':
            withfigs=True
        if myarg=='-p':
            printfile=sys.argv[i+1]

    main(inputfile,outputfile,onlybpregion,verbose,mysteps,withread,withfigs,printfile)
    










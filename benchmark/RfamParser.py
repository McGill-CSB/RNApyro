import os,sys,string,re,getopt,operator,math,random;

ntdef={};
ntdef['A']=['A'];
ntdef['C']=['C'];
ntdef['G']=['G'];
ntdef['U']=['U'];
ntdef['R']=['A','G'];
ntdef['Y']=['C','U'];
ntdef['K']=['G','U'];
ntdef['M']=['A','C'];
ntdef['S']=['G','C'];
ntdef['W']=['A','U'];
ntdef['B']=['C','G','U'];
ntdef['D']=['A','G','U'];
ntdef['H']=['A','C','U'];
ntdef['V']=['A','C','G'];
ntdef['N']=['A','C','G','U'];
ntdef['.']=[];

###############################################################################
# dummy
def mutate(sequence,nom):
    mmap={}
    while (len(mmap)<nom):
        imut = random.randint(0,len(sequence)-1)
        if not mmap.has_key(imut):
            inuc = random.randint(0,3)
            if inuc==0:
                nuc='A'
            elif inuc==1:
                nuc='C'
            elif inuc==2:
                nuc='G'
            elif inuc==3:
                nuc='U'
            else:
                print "ERROR"
                sys.exit(1)
            if nuc != sequence[imut]:
                mmap[imut]=nuc

    mutant = ''
    for i in range(len(sequence)):
        if mmap.has_key(i):
            mutant += mmap[i]
        else:
            mutant += sequence[i]

    return mutant

###############################################################################

def computeEntropy(seqdic,target,maxsize):

    minseqs=2;
    lenseq=0;
    outdic={};
    
    outdic['gc']=0.0;
    ntfrac={};
    for i in range(maxsize):
        ntfrac[i]={};
        ntfrac[i]['A']=0.0;
        ntfrac[i]['C']=0.0;
        ntfrac[i]['G']=0.0;
        ntfrac[i]['U']=0.0;

    lenseq=len(target);

    if len(seqdic) > 0:
        for id,seq in seqdic.iteritems():
            if lenseq!=len(seq):
                print "sequencelist corrupted";
                sys.exit(1);
            outdic['gc'] += float(seq.count('G')+seq.count('C'))/len(seq);
            for i in range(len(seq)):
                for newnt in ntdef[seq[i]]:
                    ntfrac[i][newnt]+=1.0/len(ntdef[seq[i]]);
        outdic['gc'] /= len(seqdic);

    # compute entropy

    outdic['nseqs'] = len(seqdic);
    ntlist=['A','C','G','U'];
    outdic['avg']=0.0;
    cmptbp=0;
    if outdic['nseqs'] >= minseqs:
        for i in range(lenseq):
            outdic[i]=0.0;
            for nt in ntlist:
                if ntfrac[i][nt] > 0:
                    myfrac = float(ntfrac[i][nt])/outdic['nseqs'];
                    val = -1.0 * myfrac * math.log(myfrac,4);
                else:
                    myfrac = 0.0;
                    val = 0.0
                outdic[i] += val;
                if True or target[i]=='(' or target[i]==')':
                    outdic['avg'] += val;
            if True or target[i]=='(' or target[i]==')':
                cmptbp+=1;
    outdic['avg'] /= cmptbp;

    return outdic;

###############################################################################

def dic2string(bpdic):
  bpstring='';
  lip=bpdic.keys();
  lip.sort();
  for ip in lip:
    if ip<bpdic[ip]:
      bpstring+='[';
    else:
      bpstring+=']';
  return bpstring;

###############################################################################

def fitstruct2seq(seq,struct,removeNonWC):

  if len(seq)!=len(struct):
    print "sequence and structure lengths do not match";
    sys.exit(1);

  # build bpdic consistent with seq
  bpdic={};
  stack=[];
  for k in range(len(struct)):
    carac=struct[k];
    if carac=='(':
      stack.append(k);
    elif carac==')':
      closei=stack.pop();
      openi=k;
      if seq[openi]!='.' and seq[closei]!='.':
        bpdic[openi]=closei;
        bpdic[closei]=openi;
  # build structure
  newseq='';
  newstruct='';
  for k in range(len(seq)):
    nt=seq[k];
    if nt!='.':
      newseq+=nt;
      if bpdic.has_key(k):
        if bpdic[k]>k:
          newstruct+='(';
        else:
          newstruct+=')';
      else:
        newstruct+='.';

  return newseq,newstruct;

#################################################################################

def hashcode(sequenceWithGaps):
    scode=''
    for i in range(len(sequenceWithGaps)):
        if sequenceWithGaps[i]=='.':
            scode += '1'
        else:
            scode += '0'
    return scode

def clusterUnGapped(data):
    clusters={}
    iddic={'num2code':{},'code2num':{}}
    cmpt=0
    for myid,myseq in data['seqs'].iteritems():
        mycode = hashcode(myseq)
        if not clusters.has_key(mycode):
            clusters[mycode]={}
            iddic['num2code'][cmpt]=mycode
            iddic['code2num'][mycode]=cmpt
            cmpt+=1
        clusters[mycode][myid]=myseq
    return clusters,iddic

#################################################################################

def main(filename,nmut,idc,listflag):

    seqline_re = re.compile("(\S+)(\s+)(\S+)");

    # read file store lines

    f = open(filename,"r")
    lines = f.readlines()
    f.close()

    # enum
    info={};
    info['seqs']={};
    lindex = range(len(lines))
    for index in lindex:

      cline = lines[index].replace('\n','');
      
      if len(cline)>0:
        if cline[0]=='#':
          if cline.startswith("#=GC SS_cons"):
            if not info.has_key('consensus'):
              info['consensus']='';
            info['consensus'] += cline[13:].lstrip().replace('<','(').replace('>',')')
          if cline.startswith("#=GF AC"):
            id = cline[8:].replace(' ','');
            info['id']=id;
          if cline.startswith("#=GF ID"):
            info['name']=cline[8:].replace(' ','');
          if cline.startswith("#=GF SQ"):
            info['nseq']=int(cline[8:].replace(' ',''));
        elif cline.startswith("//"):
          if not info.has_key('name'):
            print "WARNING: %s missing name." % (id);
          if not info.has_key('nseq'):
            print "WARNING: %s missing number of sequences." % (id);
          if not info.has_key('consensus') or len(info['consensus'])==0:
            print "WARNING: %s missing consensus structure" % (id);
        else:
          # store sequences
          fields = seqline_re.match(cline);
          if fields:
            seqid = fields.group(1);
            seqcontent = fields.group(3);
          else:
            print "Cannot parse \"%s\"" % cline;
            sys.exit(1);
          if not info['seqs'].has_key(seqid):
            info['seqs'][seqid] = '';
          info['seqs'][seqid] += seqcontent;

    cleandata,idMap = clusterUnGapped(info)

    if listflag:
        for hcode in cleandata.keys():
            if len(cleandata[hcode])>1:
                print '>>>',idMap['code2num'][hcode],', size=' + str(len(cleandata[hcode]))
        sys.exit(1)
            
    hcode=idMap['num2code'][idc]
    listout = []
    listout.append(random.randint(0,len(cleandata[hcode])-1))
    for iseqout in listout:
        ssamatch=True
        prevssa=''
        idout=''
        cmpt=0
        print '>>> cluster',idc,'size =',len(cleandata[hcode])
        for myid,myseq in cleandata[hcode].iteritems():
            cleanseq,cleanssa = fitstruct2seq(myseq,info['consensus'],False)
            if cmpt!=iseqout:
                print '>',myid
                print cleanseq
            else:
                idout=myid
            prevssa=cleanssa
            cmpt+=1
            if prevssa != '' and prevssa != cleanssa:
                ssamatch=False
        print '> consensus'
        print prevssa
        if not ssamatch:
            print 'WARNING: Structure do not match'

        print '>',idout,'(' + str(nmut) + '-mutant)'
        cleanseq = cleandata[hcode][idout].replace('.','')
        mutant = mutate(cleanseq,nmut)
        print mutant
    
    
    sys.exit(1)

##########################################################################################

if __name__ == '__main__':

  try:
      opts, args = getopt.getopt(sys.argv[1:], "hf:m:n:l", ["help", "file=", "mut=", "id=", "list"]);
  except getopt.GetoptError:
      usage(sys.argv[0]);

  filename='';
  nmut=0;
  id=0;
  listflag=False;

  argStart=len(sys.argv);
  for o,a in opts:
    if o in ("-h", "--help"):
      usage(sys.argv[0]);
    if o in ("-f", "--file"):
      filename = a;
      argStart-=2;  
    if o in ("-m", "--file"):
      nmut = int(a);
      argStart -=2;  
    if o in ("-n", "--file"):
      id = int(a);
      argStart-=2;  
    if o in ("-l", "--list"):
      listflag = True;
      argStart-=1;  

  main(filename,nmut,id,listflag);



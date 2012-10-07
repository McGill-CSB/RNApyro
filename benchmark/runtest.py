import sys,os,os.path

outdir=""

def main(nmut,alpha,prefix,penalty):

    for numex in range(100):
        infile = "RF00001_%d.in" % numex
        if os.path.exists(infile):
            alphanorm = int (alpha * 100)
            outfile = "%sRF00001_%d_nmut_%d_alpha_%d.out" % (outdir,numex,nmut,alphanorm)
            commandline1="python ../src/RNAPyroEx.py %s %d %.1f %.2f > %s" % (infile,nmut,alpha,penalty,outfile)
            print commandline1
            if not os.path.exists(outfile):
                os.system(commandline1)
            else:
                print "WARNING: %s already present" % (outfile)
            commandline2="python ../scripts/benchmark.py -i RF00001_%d.ref -o RF00001_%d_nmut_%d_alpha_%d.out" % (numex,numex,nmut,alphanorm)
            print commandline2
            os.system(commandline2)


###############################################################################

if __name__ == '__main__':
    
    nmut=10
    alpha=0.5
    prefix="RF00001"
    penalty=1.5
    
    for i in range(len(sys.argv)):
        myarg=sys.argv[i]
        if myarg=='-h' or myarg=='--help':
            usage(sys.argv[0])
        if myarg=='-m':
            nmut=int(sys.argv[i+1])
        if myarg=='-a':
            alpha=float(sys.argv[i+1])
        if myarg=='-p':
            prefix=sys.argv[i+1]
        if myarg=='-e':
            penalty=float(sys.argv[i+1])
    
    
    main(nmut,alpha,prefix,penalty)





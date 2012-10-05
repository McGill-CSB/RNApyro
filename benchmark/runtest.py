import sys,os,os.path

def main(nmut,alpha,prefix):

    for numex in range(10):
        infile = "RF00001_%d.in" % numex
        if os.path.exists(infile):
            alphanorm = int (alpha * 100)
            commandline1="python ../src/RNAPyroEx.py RF00001_%d.in %d %.1f > RF00001_%d_alpha_%d.out" % (numex,nmut,alpha,numex,alphanorm)
            os.system(commandline1)
            commandline1="python ../scripts/benchmark.py -i RF00001_%d.ref -o RF00001_%d_alpha_%d.out" % (numex,numex,alphanorm)
            os.system(commandline1)



###############################################################################

if __name__ == '__main__':
    
    nmut=10
    alpha=0.5
    prefix="RF00001"
    
    if len(sys.argv)==1:
        usage(sys.argv[0])
    
    for i in range(len(sys.argv)):
        myarg=sys.argv[i]
        if myarg=='-h' or myarg=='--help':
            usage(sys.argv[0])
        if myarg=='-n':
            nmut=int(sys.argv[i+1])
        if myarg=='-a':
            alpha=float(sys.argv[i+1])
        if myarg=='-p':
            prefix=sys.argv[i+1]
    
    
    main(nmut,alpha,prefix)





#!/usr/bin/env python
# encoding: utf-8


from pyroc import *
import pylab

def loadMutationalProfiles(path):
	fileHandler = open(path,'r')
	lines = [line[:-1] for line in fileHandler.readlines()]
	seq = lines[0]
	struct = lines[1]
	nbMut = int(lines[2])
	profiles = eval(lines[3])
	fileHandler.close()
	return seq,struct,nbMut,profiles

BASES = ["A","C","G","U"]

def classifyProfile(profile,seq,struct,mode):
	result = []
	for i in range(len(profile)):
		if mode == "a":
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				if b.lower() == seq[i].lower():
					result.append((1,prob))
				else:
					result.append((0,prob))
		elif mode == "b" and struct[i]!=".":
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				if b.lower() == seq[i].lower():
					result.append((1,prob))
				else:
					result.append((0,prob))
		elif mode == "c" and struct[i]!=".":
			bestB,probB = "",0.
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				if prob > probB:
					bestB,probB = b,prob
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				if (b.lower() == seq[i].lower()):
					cl=1
				else:
					cl=0
				if b.lower() == bestB.lower():
					result.append((cl,prob))
				else:
					result.append((cl,2.0))
	print result
	return result

def computePlotROCs(path,mode):
	seq,struc,nbMut,profiles = loadMutationalProfiles(path)
	rocs = []
	lbls = []
	pylab.rc('text', usetex=True)
	for alpha in profiles:
		c = classifyProfile(profiles[alpha],seq,struc,mode)
		roc = ROCData(c)
		rocs.append(roc)
		lbls.append("$\\alpha=%.1f,$ {\\sf AUC} $=%.1f\\%%$"%(alpha,(100.*roc.auc(0))))
		print "Alpha: %s, ROC AUC: %s" % (alpha,str(roc.auc(0)))
	plt = plot_multiple_roc(rocs,"", lbls)
	return plt

if __name__ == '__main__':
       	print "createROC - ROC Curve Generator from RNAPyro profiles"

        
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-f', '--file', dest='origFile', help="Path to a file with the mutational profile.")
	parser.add_option("-o",'--output', dest= 'outpath' , default='' , help = 'Destination (must be PDF).')
	parser.add_option("-m",'--mode', dest= 'mode' , default='a' , help = 'Evaluation mode: "a" (default) - Conservative prediction, predicts multiple mutations per position, anywhere in the sequence; "b" - Multiple mutations per positions, only in helices; "c" - Single mutations, only in helices')
	
	(options,args) = parser.parse_args()
       	if (not options.origFile):
		parser.print_help()
		exit()


	if (not options.origFile):
		parser.print_help()
		exit()
	plt = computePlotROCs(options.origFile,options.mode)
       	if options.outpath:
               	plt.savefig(options.outpath, format='pdf',bbox_inches='tight')
        else:
                plt.show()


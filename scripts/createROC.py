#!/usr/bin/env python
# encoding: utf-8


from pyroc import *
import pylab
import os

def loadMutationalProfiles(path):
	fileHandler = open(path,'r')
	lines = [line[:-1] for line in fileHandler.readlines()]
	seq = lines[0]
	struct = lines[1]
	nbMut = int(lines[2])
	profiles = eval(lines[3])
	for alpha in profiles:
		profile = profiles[alpha]
		for i in range(len(profile)):
			for j in range(len(profile[i])):
				profile[i][j] = float(profile[i][j])
	fileHandler.close()
	return seq,struct,nbMut,profiles

BASES = ["a","c","g","u"]

def decide(classes,b,mutchar,seqchar,prob):
        if b!=mutchar:
                if b!=seqchar:
                        classes.append((0,prob))
                else:
                        classes.append((1,prob))

def classifyProfile(profile,seq,struct,mode,mutseq):
	result = []
	mutseq = mutseq.lower()
	seq = seq.lower()
	for i in range(len(profile)):
		if mode == "a":
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j].lower()
				decide(result,b,mutseq[i],seq[i],prob)
		elif mode == "b" and struct[i]!=".":
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				decide(result,b,mutseq[i],seq[i],prob)
		elif mode == "c" and struct[i]!=".":
			bestB,probB = "",0.
			for j in range(len(profile[i])):
				prob = profile[i][j]
				b = BASES[j]
				if prob > probB:
					bestB,probB = b,prob
			for j in range(len(profile[i])):
				if b.lower() == bestB.lower():
               				decide(result,b,mutseq[i],seq[i],prob)
				else:
               				decide(result,b,mutseq[i],seq[i],0.)
	#print result
	return result


def evaluatePrediction(profile,seq,nbMut):
	ok = 0.
	seq = seq.lower()
	for i in range(len(profile)):
		bestB,probB = "",0.
		for j in range(len(profile[i])):
			prob = profile[i][j]
			b = BASES[j]
			if prob > probB:
				bestB,probB = b,prob
		if bestB==seq[i]:
			ok += 1
	return nbMut-(len(seq)-ok)

def getMutatedPositions(seq,mutseq):
        result = set()
        if len(seq)!= len(mutseq):
                raise Exception("Mutated and original sequence lengths do not match!")
        seq = seq.lower()
        mutseq = mutseq.lower()
        for i in range(len(seq)):
                if (seq[i]!=mutseq[i]):
                        result.add(i)
        return result

def computePlotROCs(paths,mode,mutpaths):
	allalphas = set()
	exps = []
	pylab.rc('text', usetex=True)
	for i in range(len(paths)):
		path = paths[i]
		mutpath = mutpaths[i]
		mutseq = open(mutpath,"r").readlines()[0]
		seq,struc,nbMut,profiles = loadMutationalProfiles(path)
		mutpos = getMutatedPositions(seq,mutseq)
		exps.append((path,mutseq,seq,struc,nbMut,profiles,mutpos))
		allalphas.update(profiles.keys())
	rocs = []
	lbls = []
	alphas = list(allalphas)
	alphas.sort()
	for alpha in alphas:
		c = []
		delta = 0
		for (path,mutseq,seq,struc,nbMut,profiles,mutpos) in exps:
			c += classifyProfile(profiles[alpha],seq,struc,mode,mutseq)
			delta += evaluatePrediction(profiles[alpha],seq,len(mutpos))
		roc = ROCData(c)
		rocs.append(roc)
		lbls.append("$\\alpha=%.1f,$ {\\sf AUC} $=%.1f\\%%$"%(alpha,(100.*roc.auc(0))))
		if len(paths)==1:
			cap = paths[0].split(os.sep)[-1]
		else:
			cap = "Multiple"
		print "%s\t%s\t%s\t%.2f\t%s\t%s\t%s" % (cap,mode,alpha,(100.*roc.auc(0)),delta,nbMut,len(mutpos))
	plt = plot_multiple_roc(rocs,"", lbls) 
	return plt


def listConversionCallback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

if __name__ == '__main__':
	#print "createROC - ROC Curve Generator from RNAPyro profiles"

	
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='origDir', help="Path to a directory with mutational profiles (full generation mode).")
	parser.add_option("-o",'--output', dest= 'outpath' , default='out.pdf' , help = 'Destination (must be PDF).')
	parser.add_option("-m",'--mode', dest= 'mode' , default='a' , help = 'Evaluation mode: "a" (default) - Conservative prediction, predicts multiple mutations per position, anywhere in the sequence; "b" - Multiple mutations per positions, only in helices; "c" - Single mutations, only in helices')
	parser.add_option("-f", "--file", dest="origFile",type='string',action='callback',callback=listConversionCallback,help="Path to a (set of) file(s) with the mutational profile.",default=[])
	parser.add_option("-k",'--muts', dest= 'mutpath' ,type='string',action='callback',callback=listConversionCallback, help = 'Path to a (set of) mutated sequence(s).')
  
	(options,args) = parser.parse_args()
	if ((not options.origFile) and (not options.origDir)) or (not options.mutpath) :
		parser.print_help()
		exit()

	#print options.origFile, options.mutpath
	if options.origFile:
		plt = computePlotROCs(options.origFile,options.mode,options.mutpath)
		plt.savefig(options.outpath, format='pdf',bbox_inches='tight')
	elif options.origDir:
		for f in os.listdir(options.origDir):
			if (f.endswith(".txt")):
				basename = f[:-4]
				for m in ["a","b","c"]:
					fname = options.origDir+os.sep+basename+"-"+m+".pdf"
					plt = computePlotROCs(options.origDir+os.sep+f,m,options.mutpath)
					plt.savefig(fname, format='pdf',bbox_inches='tight')


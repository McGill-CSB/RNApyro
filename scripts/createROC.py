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
					result.append((cl,0.0))
	#print result
	return result

def evaluatePrediction(profile,seq,nbMut):
	ok = 0.
	for i in range(len(profile)):
		bestB,probB = "",0.
		for j in range(len(profile[i])):
			prob = profile[i][j]
			b = BASES[j]
			if prob > probB:
				bestB,probB = b,prob
		if bestB==seq[i]:
			ok += 1
	return ok-(len(seq)-nbMut)

def computePlotROCs(path,mode):
	seq,struc,nbMut,profiles = loadMutationalProfiles(path)
	rocs = []
	lbls = []
	pylab.rc('text', usetex=True)
	alphas = profiles.keys()
	alphas.sort()
	for alpha in alphas:
		c = classifyProfile(profiles[alpha],seq,struc,mode)
		delta = evaluatePrediction(profiles[alpha],seq,nbMut)
		roc = ROCData(c)
		rocs.append(roc)
		lbls.append("$\\alpha=%.1f,$ {\\sf AUC} $=%.1f\\%%$"%(alpha,(100.*roc.auc(0))))
		print "File:%s, Mode:%s, Alpha:%s, ROC AUC:%.2f, Delta:%s" % (path.split(os.sep)[-1],mode,alpha,(100.*roc.auc(0)),delta)
	plt = plot_multiple_roc(rocs,"", lbls)
	return plt

if __name__ == '__main__':
	print "createROC - ROC Curve Generator from RNAPyro profiles"

	
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-f', '--file', dest='origFile', help="Path to a file with the mutational profile.")
	parser.add_option('-d', '--dir', dest='origDir', help="Path to a directory with mutational profiles (full generation mode).")
	parser.add_option("-o",'--output', dest= 'outpath' , default='out.pdf' , help = 'Destination (must be PDF).')
	parser.add_option("-m",'--mode', dest= 'mode' , default='a' , help = 'Evaluation mode: "a" (default) - Conservative prediction, predicts multiple mutations per position, anywhere in the sequence; "b" - Multiple mutations per positions, only in helices; "c" - Single mutations, only in helices')
	
	(options,args) = parser.parse_args()
	if (not options.origFile) and (not options.origDir) :
		parser.print_help()
		exit()

	if options.origFile:
		plt = computePlotROCs(options.origFile,options.mode)
		plt.savefig(options.outpath, format='pdf',bbox_inches='tight')
	elif options.origDir:
		for f in os.listdir(options.origDir):
			if (f.endswith(".txt")):
				basename = f[:-4]
				for m in ["a","b","c"]:
					fname = options.origDir+os.sep+basename+"-"+m+".pdf"
					plt = computePlotROCs(options.origDir+os.sep+f,m)
					plt.savefig(fname, format='pdf',bbox_inches='tight')


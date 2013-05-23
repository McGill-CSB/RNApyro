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

def classifyProfile(profile,seq):
	result = []
	for i in range(len(profile)):
		for j in range(len(profile[i])):
			prob = profile[i][j]
			b = BASES[j]
			if b.lower() == seq[i].lower():
				result.append((1,prob))
			else:
				result.append((0,prob))
	return result

def computePlotROCs(path):
	seq,struc,nbMut,profiles = loadMutationalProfiles(path)
	rocs = []
	lbls = []
	pylab.rc('text', usetex=True)
	for alpha in profiles:
		c = classifyProfile(profiles[alpha],seq)
		roc = ROCData(c)
		rocs.append(roc)
		lbls.append("$\\alpha=%.1f,$ {\\sf AUC} $=%.1f\\%%$"%(alpha,(100.*roc.auc(0))))
		print "Alpha: %s, ROC AUC: %s" % (alpha,str(roc.auc(0)))
	plt = plot_multiple_roc(rocs,"", lbls)
	return plt

if __name__ == '__main__':
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option('-f', '--file', dest='origFile', help="Path to a file with the class and decision function. The first column of each row is the class, and the second the decision score.")
	parser.add_option("-o",'--output', dest= 'outpath' , default='' , help = 'Destination.')
	
	(options,args) = parser.parse_args()


	if (not options.origFile):
		parser.print_help()
		exit()
	plt = computePlotROCs(options.origFile)
       	if options.outpath:
               	plt.savefig(options.outpath, format='pdf',bbox_inches='tight')
        else:
                plt.show()


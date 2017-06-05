#! /usr/bin/python

#
#Script fills given amount of base pairs with sequences from provided list.
#The selection can be completely random, or it can select sequences to maximize 
#the number of genes they cover
#
__author__ = "Matus Kempa"
__date__ = "$Feb 23, 2016 9:47:00 AM$"

import sys
import getopt
import random
import re
import operator

class Probes:
    
    def __init__(self):
        self.thesaurus = [] #list to choose from
        self.genesPicked = dict() # genes already picked
    
    """
    Loads sequences to choose from but only those not in the 'excluded' file
    infile - file containing all sequences
    excluded - file containing names of sequences to exclude from
    maxgenes - true if we want to maximize number of genes covered
    """
    def loadSequences(self, infile, excluded, maxgenes = False):
        print "Loading excluded sequences"
        e = open(excluded, 'r')
        excl = set() #probes to exclude from random search
        #load excluded annotations into list
        line = e.readline()
        while line:
            l = line.strip()
            excl.add(l)
            line = e.readline()
        e.close()
        print "Loaded ", len(excl), " sequences for exclusion"
        
        if maxgenes:
            print "Flag to maximize genes is True, initializing list of used genes"
            for s in excl: #load genes into dictionary. each gene name has coverage (number of loci covering it) assigned
                gname = s[9:21]
                
                if gname in self.genesPicked:
                    self.genesPicked[gname] += 1
                else:
                    self.genesPicked[gname] = 1
            print len(self.genesPicked), " genes loaded"
        
        print "Creating thesaurus"
        f = open(infile, 'r')
        #now we create the thesaurus, we don't want sequences in the exclusion list
        line = f.readline()
        while line:
            if line[0] == '>': #annotation
                annot = line.strip()
                line = f.readline() #sequence line
                if annot[1:] in excl: #annotation is in the list of excluded
                    line = f.readline()
                    continue
                if maxgenes and (annot[10:22] in self.genesPicked): #gene the annotation is covering is in the list of used genes
                    line = f.readline()
                    continue
                if re.match('[ACGT]+', line): #sequence
                    seq = line.strip()
                    self.thesaurus.append((annot[1:], seq, len(seq))) #storing a triple (annotation, sequence, length of sequence)
                else:
                    print "no sequence found"
                line = f.readline()
            else:
                break
        f.close()
        print "Thesaurus contains ", len(self.thesaurus), " sequences"

    def selectRandom(self, maxLength, minSeqLength, maxgenes = False):
        print "Picking random baits"
        geneExcl = dict(self.genesPicked)
        totalLength = 0
        final = [] #list of result probes
        thesLength = len(self.thesaurus)
        numgenes = self.genesCount([a for (a, s, l) in self.thesaurus]) + len(geneExcl) #count the total number of genes available
        #used = [] #storing indices of already chosen probes to avoid duplicity
        
        available = [i for i in range(0, thesLength - 1)] #available indices of sequences for random picking
        again = False
        
        while True:
            r = random.choice(available)
            el = self.thesaurus[r] #pick element at that index
            gene = el[0][9:21]
            
            if maxgenes and gene in geneExcl and not again: # check if covers used gene when considering genes
                continue
            
            available.remove(r) #mark the index as used
            if totalLength + el[2] <= maxLength:
                final.append(el)
                totalLength += el[2]
                if maxgenes: # add to the list of used genes only if the sequence is used
                    if gene in geneExcl:
                        geneExcl[gene] += 1
                    else:
                        geneExcl[gene] = 1 #adding gene for the first time
                    #print "Adding to used genes: ", el[0][9:21],". Total ", len(self.geneExcl), " genes."
            
            if maxLength - totalLength < minSeqLength:
                break #there is not enough room for any sequences to be added since we have defined shortest possible sequence
                
            #all genes are covered but there is still space left that must be filled
            if maxgenes and len(geneExcl) == numgenes and not again: #reached maximum number of genes
                #print "Reached maximum number of genes", len(geneExcl), ", going again"
                again = True

        #print "Total lenght of selected sequences: ", totalLength #total length of selected sequences
        #print "Number of selected sequences: ", len(final) #Number of selected sequences
        #print "Number of genes covered: (with genes provided in exclusion list)", len(self.geneExcl) #Number of genes covered: (with genes provided in exclusion list)
        #print len(available)
        
        print totalLength, "#", len(final), "#", len(geneExcl) if maxgenes else self.genesCount([a for (a, s, l) in final])
        return final
            
    def writeFasta(self, probes, outfile):
        f = open(outfile, 'w')
        for (annot, seq, l) in probes:
            f.write('>' + annot + "\n")
            f.write(seq + "\n")
            
        f.close()
        
    def genesCount(self, annotations):
        start = 10 if annotations[1] == '>' else 9;
        end = 22 if annotations[1] == '>' else 21;
        l = [x[start:end] for x in annotations]
        return len(set(l))

def init(argv):
    input_file = ''
    exclude_file = ''
    maxbp = 0
    output_file = ''
    maxgenes = False #if true, algorithm tries to cover as many genes as possible, else genes not considered and selection is totally random
    try:
        opts, args = getopt.getopt(argv, "hi:e:c:o:m", ["ifile=", "efile=", "bc-count=", "ofile=", "max-genes"])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -e <excludefile> -c <maxbp> [-m] -o <outputfile>'
        print 'use test.py -h for help'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'test.py -i <inputfile> -e <excludefile> -c <maxbp> -m -o <outputfile>'
            print '-i, --ifile \t path to the file containing all bait sequences as were output by Sondovac Part_b'
            print '-e, --efile \t path to the file containing names of the sequences to be excluded from the random search'
            print '-c, --maxbp \t number of base pairs that we want to fill'
            print '-m, --max-genes \t flag to indicate maximizing targeted genes'
            print '-o, --ofile \t outfile'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-e", "--efile"):
            exclude_file = arg
        elif opt in ("-c", "--maxbp"):
            maxbp = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg
        elif opt in ("-m", "--max-genes"):
            maxgenes = True
    return [input_file, exclude_file, maxbp, output_file, maxgenes]


if __name__ == "__main__":
    
    cmdln = True #run in commandline?
    
    if cmdln:
        args = main(sys.argv[1:])
        pr = Probes()
        pr.loadSequences(args[0], args[1], args[4])
        final = pr.selectRandom(int(args[2]), args[4])
        pr.writeFasta(final, args[3])
    #else:
    #    pr = Probes()
    #    pr.loadSequences('workfiles/parviflora_probes.fasta', 'workfiles/probes_aligned_to_genes_from_papers.txt', True)
    #    final = pr.selectRandom(729032, 120, True)
    #    pr.writeFasta(final, 'workfiles/random_probes_t.fasta')
    
        #-----------
        #for testing
        #-----------
    #    for i in range(0, 10) :
    #        print i, " - run"
    #        final = pr.selectRandom(729032, 120, True)
    #        pr.writeFasta(final, 'workfiles/random_probes_' + str(i) + '.fasta')
    #        print "----------------"
    
        
        print 'done'

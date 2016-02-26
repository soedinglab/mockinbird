import os
import sys
from stammp.obj import functions

class Genome:
    """
        :class:`~stammp.obj.genome.Genome` represents the sequence data of a complete genome. Sequences are loaded upon object creation. *location* should be multiple fasta or a directory which only contains one fasta file per chromosome.
        
        .. warning:: This is a bit tricky. The final product of the pre-process procedure is a pileup file which stores sequencing information for a chromosome-identifier and a position at that chromosome. The :ref:`ref_binding-site-detection` also stores this chromosome information. In order to use post-processing steps correctly the identifiers that are used in the pileup file have to be identical to the fasta header in the multiple fasta file (everything which follow after '>') or to the fasta filenames of directory of the :class:`~stammp.obj.genome.Genome` object. Otherwise you won't be able to access the correct sequences. Take a look at the :ref:`ref_tutorial` for further explanations.
        
        
        Args:
            location (str): Mutliple fasta file or path to a directory only containing one fasta file for each chromosome of the genome.
            verbose (bool): 
        Example:
            >>> from stammp.obj import genome
            >>> yeast = genome.Genome('/path/to/fastafiles/')
    """
    def __init__(self, location, verbose = True):
        if os.path.isdir(location):
            self.genome = {}
            filenames   = os.listdir(location)
            count = 0
            for f in filenames:
                count += 1
                self.genome.update({f.split('.')[0] : self.readFasta(str(location+f))})
                if verbose: functions.showProgress(count, len(filenames),'Loading genome')
            if verbose: print('')
        elif os.path.isfile(location):
            self.genome = {}
            fc = open(location, 'r')
            line = fc.readline()
            tmp_name = line.split('>')[1].split('\n')[0]
            seq = ''
            line = fc.readline()
            while line:
                if line[0] == '>':
                    self.genome[tmp_name] = seq
                    seq = ''
                    tmp_name = line.split('>')[1].split('\n')[0]
                else:
                    line = line.replace('\n', '')
                    seq  = seq+line.upper()
                line = fc.readline()
            fc.close()
        else:
            sys.stderr.write('No file or directory. Exit')
            #sys.exit(-1)
    
    def readFasta(self, filename):
        """
        Returns a single string, containing the sequence of *filename*.fasta.
        
        Args:
            filename (str): Input \*.fasta
        Returns:
            str : Content of *filename*
        """
        fc = open(filename, 'r')
        line = fc.readline()
        line = fc.readline()
        seq = ''
        while(line):
            line = line.replace('\n', '')
            seq  = seq+line.upper()
            line = fc.readline()
        fc.close()
        return(seq)
    
    def getSequence(self,chrname, start, stop):
        """
        Returns the genomic sequence as specified by *chrname*, *start* and *stop* or -1 if the parameters are not part of the genome.
        
        Args:
            chrname (str): chromosome identifier
            start (int): start index
            stop (int): stop index
        """
        try:
            if chrname in self.genome:
                if start < 0 or stop < start or stop < 0 or start > (len(self.genome[chrname])-1) or stop > (len(self.genome[chrname])-1):
                    raise Exception
                else:
                    return self.genome[chrname][start:stop]
            else:
                raise Exception
        except:
            return (-1)
    
    def __str__(self):
        total = 0
        s = ''
        k = list(self.genome.keys())
        k.sort()
        for key in k:
            total += len(self.genome[key])
            s += str(key)+'\t: '+str(len(self.genome[key]))+' nt\n'
        s += '\n'
        s += 'Total\t: '+str(total)+' nt\n'
        return(s)

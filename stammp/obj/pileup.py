import math
import os
from stammp.obj import functions
class Pileup:
    
    def __init__(self, filename):
        if os.path.isfile(filename):
            self.chr = []
            self.pos = []
            self.coverage_forward = []
            self.coverage_reverse = []
            self.chrPositions = {}
            fc = open(filename, 'r')
            line = fc.readline()
            while line:
                split = line.split('\t')
                self.chr.append(split[0])
                self.pos.append(int(split[1]))
                self.coverage_forward.append(functions.getCounts(split[4], True)[0])
                self.coverage_reverse.append(functions.getCounts(split[4], False)[0])
                line = fc.readline()
            fc.close()
            self.getChromosomePositions()
        else:
            raise Exception
    
    def size(self):
        return(len(self.chr))
    
    def print(self, start, stop):
        if start >= 0 and start < (self.size()-1) and stop >= start:
            for i in range(start, stop):
                print(self.chr[i]+'\t'+str(self.pos[i])+'\t'+str(self.coverage_forward[i])+'\t'+str(self.coverage_reverse[i]))
    
    def getChromosomePositions(self):
        cur_chr = self.chr[0]
        start = 0
        for i in range(self.size()):
            if cur_chr != self.chr[i]:
                self.chrPositions[cur_chr] = [start, (i-1)]
                start = i
                cur_chr = self.chr[i]
        self.chrPositions[cur_chr] = [start, (self.size()-1)]
    
    def exactSearch(self, chrname, position):
        """
        Binary search which returns the index of the pileup entry, if a genomic position is found at the given parameters and -1 otherwise.
        
        Args:
            chrname (str): chromosome identifier
            position (int): chromsome coordinate
        """
        count = 0
        try:
            index_min = self.chrPositions[chrname][0]
            index_max = self.chrPositions[chrname][1]
            if position < self.pos[index_min]:
                return [index_min, False]
            if position > self.pos[index_max]:
                return [index_max, False]
            while index_min <= index_max:
                index_mid = index_min + math.trunc((index_max-index_min)/2)
                if position == self.pos[index_mid]:
                    return [index_mid, True]
                if position < self.pos[index_mid]:
                    index_max = (index_mid-1)
                else:
                    index_min = (index_mid+1)
                    
                count += 1
            return [index_mid, False]
        except:
            return [-1, False]
    
    def getValues(self, chrname, position, strand, upstream, downstream):
        """
        Gets all coverage values in the interval [position-upstream, position+downstream] for the given strand.
        
        Args:
            chrname (str): chromosome identidifier
            position (int): postion within the chromosome
            strand (str): \'+\' or \'-\'
            upstream (int): number of nucleotides upstream of the position
            downstream (int): number of nucleotides downstream of the position
        
        Returns:
            list : list of length *upstream* + *downstream* +1 with the occupancy value at that genomic position
        """
        values = [0]*(upstream+downstream+1)
        
        subtract = position-upstream
        start  = self.exactSearch(chrname, (position-upstream))
        stop   = self.exactSearch(chrname, (position+downstream))
        
        if start[0] >= 0 and stop[0] >= 0 and (stop[0]+1) < self.size():
            for i in range(start[0], (stop[0]+1)):
                index = self.pos[i]-subtract
                if index >= 0 and index < len(values):
                    if strand == '+':
                        values[index] = self.coverage_forward[i]
                    else:
                        values[index] = self.coverage_reverse[i]
            return values
        else:
            return None


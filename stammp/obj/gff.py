class GFF:
    """
    :class:`~stammp.obj.gff.GFF` represents the annotation data of a \*.gff file. The class offers methods to load, save, sort and modify the annotation data. The objects loads the data of a \*.gff file upon creation.
    
    The last column of a \*.gff file often contains additional annotation information. On default, the complete string is stored. However, you can set *splitA* and *splitB* to extract specifc information. If you set *splitA* and *splitB* the string is first splitted by *splitA* and the right hand side of the split is then splitted by *splitB* and the left hand side of that split is finally stored.
    
    Example:
        >>> splitA = ":"
        >>> splitB = "-"
        >>> info = "string:RESULT-toSplit"
        >>> info.split(splitA)[1].split(splitB)[0]
        'RESULT'
    
    Args:
        filename (str): Filename of the \*.gff to read.
        splitA (str): [default = None]
        splitB (str): [default = None]
    
    Example:
        >>> from stammp.obj import gff
        >>> annotation = gff.GFF('/path/to/gff-annotation.gff')
        
    .. warning:: This is a bit tricky. The final product of the pre-process procedure is a pileup file which stores the sequencing information for a chromosome-identifier and a position at that chromosome. The :ref:`ref_binding-site-detection` also stores this chromosome information. In order to use post-processing steps correctly the identifiers that are used in the pileup file have to be **identical** to the chromosome names of the used annotation for :class:`~stammp.obj.gff.GFF` object. Otherwise you won't be able to access anything. Take a look at the :ref:`ref_tutorial` for further explanations.
    
    """
    def __init__(self, filename, splitA=None, splitB=None):
        fc     = open(filename, 'r')
        line   = fc.readline()
        self.chr    = []
        self.start  = []
        self.stop   = []
        self.strand = []
        self.info   = []
        self.typeof = []
        self.chrStartIndices = {}
        index = 0
        tmp_chr = ''
        while(line):
            split = line.split('\t')
            if len(split) == 9:
                self.chr.append(split[0])
                self.typeof.append(split[2])
                self.start.append(int(split[3]))
                self.stop.append(int(split[4]))
                self.strand.append(split[6])
                if splitA != None and splitB != None:
                    self.info.append(split[8].split(splitA)[1].split(splitB)[0]) #for SGD mRNA: splitA='Name=' and splitB='_'
                else:
                    self.info.append(split[8].split('\n')[0])
                if tmp_chr != split[0]:
                    tmp_chr = split[0]
                    self.chrStartIndices[tmp_chr] = index
                index += 1
            line = fc.readline()
        fc.close()
    
    def saveGFF(self, filename, name='.', typeof='.'):
        """
        Saves a :class:`~stammp.obj.gff.GFF` as gff file.
        
        Args:
            filename (str): filename
            name (str): String that is saved in the second column
            typeof (str): String that is saved in the third column
        """
        fc = open(filename, 'w')
        for i in range(self.size()):
            fc.write(self.chr[i]+'\t'+name+'\t'+typeof+'\t'+str(self.start[i])+'\t'+str(self.stop[i])+'\t.\t'+self.strand[i]+'\t.\t'+self.info[i]+'\n')
        fc.close()
    
    def removeEntry(self, index):
        """
        Removes an entry of a :class:`~stammp.obj.gff.GFF` object at the given *index*.
        
        Args:
            index (int): Index of the entry that will be removed
        """
        if index >= 0 and index < len(self.chr):
            del self.chr[index]
            del self.start[index]
            del self.stop[index]
            del self.strand[index]
            del self.info[index]
            del self.typeof[index]
    
    def size(self):
        """
        Returns the number of elements of a :class:`~stammp.obj.gff.GFF` object.
        
        Returns:
            (int) : Number of elements
        """
        return(len(self.chr))
    
    def print(self, start, stop):
        if start < 0 or start >= len(self.chr) or stop < start or stop > len(self.chr):
            raise IndexError
        print('Chr\tStart\tStop\tStrand\tInfo\tType')
        for i in range(start, stop):
            print(self.chr[i]+'\t'+str(self.start[i])+'\t'+str(self.stop[i])+'\t'+self.strand[i]+'\t'+self.info[i]+'\t'+self.typeof[i])
    
    def filterSize(self, minSize, maxSize):
        i = 0
        while i < len(self.chr):
            if self.stop[i] - self.start[i] < minSize or self.stop[i] - self.start[i] > maxSize:
                self.removeEntry(i)
                i -= 1
            i += 1
    
    def getAnnotationLength(self):
        l = [0]*len(self.chr)
        for i in range(len(self.chr)):
            l[i] = self.stop[i] - self.start[i]
        return(l)
    
    def sort2size(self):
        diffs = []
        for i in range(self.size()):
            diffs.append((i, self.stop[i]-self.start[i]))
        sorted_diffs = sorted(diffs, key= lambda d: d[1])
        chrs   = []
        start  = []
        stop   = []
        strand = []
        info   = []
        typeof = []
        for i in range(len(sorted_diffs)):
            chrs.append(self.chr[sorted_diffs[i][0]])
            typeof.append(self.typeof[sorted_diffs[i][0]])
            start.append(self.start[sorted_diffs[i][0]])
            stop.append(self.stop[sorted_diffs[i][0]])
            strand.append(self.strand[sorted_diffs[i][0]])
            info.append(self.info[sorted_diffs[i][0]])
        self.chr    = chrs
        self.typeof = typeof
        self.start  = start
        self.stop   = stop
        self.strand = strand
        self.info   = info




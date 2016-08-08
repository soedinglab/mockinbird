import math
from stammp.obj import functions


class ParclipSites:
    """
        :class:`~stammp.obj.parclipsites.ParclipSites` represents the data
        obtained by :ref:`ref_binding-site-detection`. The class offers methods
        to load, save, sort and modify PAR-CLIP data.

        Args:
            name (str): name of the experiment or protein
        Example:
            >>> from stammp.obj import parclipsites
            >>> sites = parclipsites.ParclipSites('name')
    """
    def __init__(self, name=''):
        self.chrs = []
        self.pos = []
        self.m = []
        self.r = []
        self.result = []
        self.strand = []
        self.occ = []
        self.name = name
        self.sorted = False
        self.chrPositions = {}

    def loadFromFile(self, filename):
        """
        Loads PAR-CLIP binding sites obtained by
        :ref:`ref_binding-site-detection` into a
        :class:`~stammp.obj.parclipsites.ParclipSites` object.

        Args:
            filename (str): filename

        Example:
            >>> from stammp.obj import parclipsites
            >>> sites = parclipsites.ParclipSites('name')
            >>> sites.loadFromFile('/path/to/parclip.table')
        """
        self.chrs = []
        self.pos = []
        self.m = []
        self.r = []
        self.result = []
        self.strand = []
        self.occ = []
        fc = open(filename, 'r')
        line = fc.readline()
        line = fc.readline()
        while(line):
            split = line.split('\t')
            if len(split) == 7:
                self.chrs.append(split[0])
                self.pos.append(int(split[1]))
                self.m.append(int(split[2]))
                self.r.append(int(split[3]))
                self.result.append(float(split[4]))
                self.strand.append(split[5])
                self.occ.append(float(split[6]))
            line = fc.readline()
        fc.close()
        if len(self.chrs) == 0:
            return
        else:
            self.getChromosomePositions()

    def save2File(self, filename):
        """
        Saves a :class:`~stammp.obj.parclipsites.ParclipSites` object as tab
        delimeted text file.

        Args:
            filename (str): Filename

        Example:
            >>> from stammp.obj import parclipsites
            >>> sites = parclipsites.ParclipSites('name')
            >>> sites.loadFromFile('/path/to/parclip.table')
            >>> sites.save2File('/path/to/another.table')
        """
        fc = open(filename, 'w')
        fc.write('chromosome\tposition\tm\tr\tpvalue\tstrand\tocc\n')
        for i in range(self.size()):
            fc.write(self.chrs[i] + '\t' + str(self.pos[i]) + '\t' +
                     str(self.m[i]) + '\t' + str(self.r[i]) + '\t' +
                     str(self.result[i]) + '\t'+self.strand[i] + '\t' +
                     str(self.occ[i])+'\n')
        fc.close()

    def print(self, start, stop):
        """
        Prints PAR-CLIP data between given *start* and *stop* indices of an
        :class:`~stammp.obj.parclipsites.ParclipSites` object to stdout

        Args:
            start (int): Start index
            stop (int): Stop index

        Example:
            >>> from stammp.obj import parclipsites
            >>> sites = parclipsites.ParclipSites('name')
            >>> sites.loadFromFile('/path/to/parclip.table')
            >>> sites.print(0,5)
             	: Chr	Position	m	r	pval	strand	occ
            0	: chrI	5375	3	8	0.001412	-	0.5
            1	: chrI	6212	2	5	0.002277	-	0.333
            2	: chrI	8955	2	5	0.002277	-	0.5
            3	: chrI	9614	2	5	0.002277	+	0.285
            4	: chrI	11454	4	19	0.003157	+	0.333
        """
        if start < 0 or stop < start or start >= len(self.chrs):
            print('start and/or stop indices are out of range!')
            raise IndexError
        print(' \t: Chr\tPosition\tm\tr\tpval\tstrand\tocc')
        for i in range(start, stop):
            print(str(i) + '\t: ' + self.chrs[i] + '\t' + str(self.pos[i]) +
                  '\t' + str(self.m[i]) + '\t' + str(self.r[i]) + '\t' +
                  str(self.result[i]) + '\t' + self.strand[i] + '\t' +
                  str(self.occ[i]))

    def head(self):
        self.print(0, 5)

    def tail(self):
        self.print(self.size()-5, self.size())

    def size(self):
        """
        Returns the number of binding sites stored in the
        :class:`~stammp.obj.parclipsites.ParclipSites` object.

        Returns:
            int
        """
        return(len(self.chrs))

    def sort(self, key='occ'):
        """
        Sorts the PAR-CLIP data of a
        :class:`~stammp.obj.parclipsites.ParclipSites` object in descending
        order according to the given *key*.

        *key* can be one of the following values:
            * *occ* [default] = sites are sorted according to their occupancy
            * *m* = sites are sorted according to the number of mutations
            * *r* = sites are sorted according to the coverage
            * *mr* = sites are sorted according to the ratio of m/r
            * *pvalue* = sites are sorted according to the p-value calculated by :mod:`stammp.scripts.bsfinder`

        Args:
            key (str): see above

        """
        tmp_chr = []
        tmp_pos = []
        tmp_result = []
        tmp_m = []
        tmp_r = []
        tmp_strand = []
        tmp_occ = []
        sortedOccupancies = []
        for i in range(len(self.chrs)):
            if key == 'occ':
                sortedOccupancies.append((i, self.occ[i]))
            if key == 'm':
                sortedOccupancies.append((i, self.m[i]))
            if key == 'r':
                sortedOccupancies.append((i, self.r[i]))
            if key == 'pvalue':
                sortedOccupancies.append((i, self.result[i]))
            if key == 'mr':
                sortedOccupancies.append((i, (self.m[i]/self.r[i])))
        sortedOccupancies = sorted(sortedOccupancies, key=lambda d: d[1])

        if key == 'pvalue':
            start = 0
            stop = len(sortedOccupancies)
            step = 1
        else:
            start = len(sortedOccupancies) - 1
            stop = 0
            step = -1
        for i in range(start, stop, step):
            tmp_chr.append(self.chrs[sortedOccupancies[i][0]])
            tmp_pos.append(self.pos[sortedOccupancies[i][0]])
            tmp_result.append(self.result[sortedOccupancies[i][0]])
            tmp_m.append(self.m[sortedOccupancies[i][0]])
            tmp_r.append(self.r[sortedOccupancies[i][0]])
            tmp_strand.append(self.strand[sortedOccupancies[i][0]])
            tmp_occ.append(self.occ[sortedOccupancies[i][0]])
        self.chrs = tmp_chr
        self.pos = tmp_pos
        self.m = tmp_m
        self.r = tmp_r
        self.result = tmp_result
        self.strand = tmp_strand
        self.occ = tmp_occ
        self.sorted = True

    def removeSite(self, index):
        """
        Removes a single binding site entry of a
        :class:`~stammp.obj.parclipsites.ParclipSites` object at the given
        index

        Args:
            index (int): binding site index that will be removed
        """
        if index >= 0 and index < len(self.chrs):
            del self.chrs[index]
            del self.pos[index]
            del self.m[index]
            del self.r[index]
            del self.result[index]
            del self.strand[index]
            del self.occ[index]
        else:
            print('Index: ' + str(index) + ' out of range. Object has only #' +
                  str(self.size()) + ' elements.')
            raise IndexError

    def addSite(self, chrname, position, m, r, pvalue, strand, occupancy):
        """
        Adds a binding site to a :class:`~stammp.obj.parclipsites.ParclipSites`
        object.

        Args:
            chrname (str): chromosome name
            position (int): position of the binding site within the chromosome
            m (int): number PAR-CLIP induced mutations at that position
            r (int): number of reads of the PAR-CLIP experiment covering that position
            pvalue (float): p-value as obtained by the statistical model
            strand (str): **+** or **-** to indicate on which strand the binding site was measured
            occupancy (float): m/RNAseq signal
        """
        self.chrs.append(chrname)
        self.pos.append(position)
        self.m.append(m)
        self.r.append(r)
        self.result.append(pvalue)
        self.strand.append(strand)
        self.occ.append(occupancy)

    def removeSitesLocatedInGFF(self, gff, width=0):
        """
        Removes binding sites of an
        :class:`~stammp.obj.parclipsites.ParclipSites` object which overlap
        with a given :class:`~stammp.obj.gff.GFF` object.

        Args:
            gff(GFF): :class:`~stammp.obj.gff.GFF` object
            width (int): additional width which is added to the annotation information [default = 0]
        """
        i = 0
        while i < self.size():
            for j in range(gff.size()):
                if gff.strand[j] == self.strand[i] and self.chrs[i] == gff.chr[j] and self.pos[i] >= (gff.start[j]-width) and self.pos[i] <= (gff.stop[j]+width):
                    self.removeSite(i)
                    i -= 1
                    break
            i += 1

    def save2Fasta(self, genome, filename, start, stop, width=15):
        """
        For each :class:`~stammp.obj.parclipsites.ParclipSites` between the
        indices *start* and *stop* genomic sequences +/- *width* nt from a
        :class:`~stammp.obj.genome.Genome` object are saved as a fasta file
        specified by *filename*.

        Args:
            genome (Genome): :class:`~stammp.obj.genome.Genome` object
            filename (str): filename of the resulting fasta file
            start (int): start index
            stop (int): Stop index
            width (int): number of nucleotides that are saved +/- the position of a :class:`~stammp.obj.parclipsites.ParclipSites` entry.
        """
        failedSeqs = 0
        if start < 0 or stop < start or start >= len(self.chrs):
            print('start and/or stop indices are out of range!')
            raise IndexError
        if stop >= len(self.chrs):
            stop = len(self.chrs)
            warning = (
                'WARNING: stop index is higher than the number of available '
                'PAR-CLIP sites. Stop index was set to %s.'
            )
            print(warning_tmpl % stop)
        fc = open(filename, 'w')
        for i in range(start, stop):
            if self.strand[i] == '+':
                seq = genome.getSequence(self.chrs[i], (self.pos[i]-width),
                                         (self.pos[i]+width+1))
                if seq != -1:
                    fc.write('>seq:' + str(i) + ':' + self.chrs[i] + ':' +
                             str(self.pos[i]) + ':' + self.strand[i] + ':' +
                             str(round(self.occ[i], 2)) + ':' + str(self.m[i]) +
                             ':' + str(round(self.m[i]/self.occ[i], 2)) + '\n')
                    fc.write(seq+'\n')
                else:
                    failedSeqs += 1
            else:
                seq = genome.getSequence(self.chrs[i], (self.pos[i]-width-2),
                                         (self.pos[i]+width-1))
                if seq != -1:
                    fc.write('>seq:' + str(i) + ':' + self.chrs[i] + ':' +
                             str(self.pos[i]) + ':' + self.strand[i] + ':' +
                             str(round(self.occ[i], 2)) + ':' + str(self.m[i]) +
                             ':' + str(round(self.m[i]/self.occ[i], 2)) + '\n')
                    fc.write(functions.makeReverseComplement(seq)+'\n')
                else:
                    failedSeqs += 1
        fc.close()
        if failedSeqs > 0:
            print('WARNING: '+str(failedSeqs)+' couldn\'t be saved into the fasta file, because the sequences requests do not belong to the current genome!')

    def getSequences(self, genome, start, stop, width):
        """
        For each :class:`~stammp.obj.parclipsites.ParclipSites` between the
        indices *start* and *stop* genomic sequences +/- *width* nt from a
        :class:`~stammp.obj.genome.Genome` object are returned as a list of
        strings.

            Args:
                genome (Genome): :class:`~stammp.obj.genome.Genome` object
                start (int): start index
                stop (int): stop index
                width (int): number of nucleotides +/- around the crosslink position of a :class:`~stammp.obj.parclipsites.ParclipSites` entry.
        """
        if start < 0 or stop < start or start >= len(self.chrs):
            print('start and/or stop indices are out of range!')
            raise IndexError
        if stop >= len(self.chrs):
            stop = len(self.chrs)
            print('WARNING: stop index is higher than the number of ' +
                  'available PAR-CLIP sites. Stop was set to '+str(stop))
        seqs = []
        for i in range(start, stop):
            if self.strand[i] == '+':
                seq = genome.getSequence(self.chrs[i], (self.pos[i]-width),
                                         (self.pos[i]+width+1))
                if seq != -1:
                    seqs.append(seq)
            else:
                seq = genome.getSequence(self.chrs[i], (self.pos[i]-width-2),
                                         (self.pos[i]+width-1))
                if seq != -1:
                    seqs.append(functions.makeReverseComplement(seq))
        return seqs

    def getChromosomePositions(self):
        if len(self.chrs) == 0:
            return
        if self.sorted:
            print('PAR-CLIP sites have been sorted according to other values. Exiting!')
            raise Exception
        cur_chr = self.chrs[0]
        start = 0
        for i in range(self.size()):
            if cur_chr != self.chrs[i]:
                self.chrPositions[cur_chr] = [start, (i-1)]
                start = i
                cur_chr = self.chrs[i]
        self.chrPositions[cur_chr] = [start, (self.size()-1)]

    def exactSearch(self, chrname, position, strand, width=0):
        """
        Binary search which returns the index of the PAR-CLIP site entry, if a
        site is found at the given parameters and -1 otherwise.

        .. warning:: If you added or removed sites after loading you have to call :func:`~stammp.obj.parclipsites.getChromosomePositions`. Otherwise you won't get correct results.

        Args:
            chrname (str): chromosome identifier
            position (int): chromsome coordinate
            strand (str): strand identifier ['+','-']
            width (int): additional width. The functions returns the index, if the position is found +/- width around an entry.
        """
        count = 0
        if self.sorted:
            print('Search function can only be used if the sites are sorted according to the genomic position. PAR-CLIP sites have been sorted according to other values. Exiting!')
            raise Exception
        try:
            index_min = self.chrPositions[chrname][0]
            index_max = self.chrPositions[chrname][1]
            if position < self.pos[index_min]:
                return [index_min, False]
            if position > self.pos[index_max]:
                return [index_max, False]
            while index_min <= index_max:
                index_mid = index_min + math.trunc((index_max-index_min)/2)
                if position >= (self.pos[index_mid]-width) and position <= (self.pos[index_mid]+width):
                    if self.strand[index_mid] == strand:
                        return [index_mid, True]
                    else:
                        return [index_mid, False]
                if position < self.pos[index_mid]:
                    index_max = (index_mid-1)
                else:
                    index_min = (index_mid+1)
                count += 1
            return [index_mid, False]
        except:
            return [-1, False]

    def getValues(self, chrname, position, strand, sense, upstream, downstream):
        """
        Gets all occupancy values in the interval [position-upstream,
        position+downstream] for the given strand. If sense is true the
        PAR-CLIP values for the given strand are returned or for the opposite
        strand if false. If the strand is set to '-' the values are already
        flipped so that they are returned in 5\' to 3\' orientation.

        Args:
            chrname (str): chromosome identidifier
            position (int): postion within the chromosome
            strand (str): \'+\' or \'-\'
            sense (bool): If true values for the *strand* are returned or for the opposite strand if false
            upstream (int): number of nucleotides upstream of the position
            downstream (int): number of nucleotides downstream of the position

        Returns:
            list : list of length *upstream* + *downstream* +1 with the occupancy value at that genomic position
        """
        values = [0]*(upstream+downstream+1)
        strand2 = strand
        if strand == '+':
            subtract = position-upstream
            start = self.exactSearch(chrname, subtract, strand, width=0)
            stop = self.exactSearch(chrname, (position+downstream), strand, width=0)
            if not sense:
                strand2 = '-'
        else:
            subtract = position-downstream
            start = self.exactSearch(chrname, subtract, strand, width=0)
            stop = self.exactSearch(chrname, (position+upstream), strand, width=0)
            if not sense:
                strand2 = '+'
        if start[0] >= 0 and stop[0] >= 0 and (stop[0]+1) < self.size():
            for i in range(start[0], (stop[0]+1)):
                if self.strand[i] == strand2:
                    index = self.pos[i]-subtract
                    if index >= 0 and index < len(values):
                        values[index] = self.occ[i]
        else:
            return None
        if strand == '-':
            values = values[::-1]
        return values

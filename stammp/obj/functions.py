"""
:mod:`stammp.obj.functions`

Basic internal module that offers functions for loading, sorting, manipulating 
and saving of data files (mpileups, gff, fasta etc.) as well as PAR-CLIP data.
Most modules of STAMMP rely on these functions.
"""
import sys
import math
import copy

#################################################################
######## General file handling functions
#################################################################

def showProgress(count, total, message):
    """
    Displays a ASCII progress bar to stdout.
    
    Example:
        >>> showProgress(50, 100, 'This is a progress bar.')
        0% [=====.....] 100%	50.0%	This is a progress bar.
    
    Args:
        count (int): Current position of iteration
        total (int): Total number of iterataions / items etc.
        message (str): Message to display besides the progress bar
    
    """
    frac = math.trunc((count/total)*10)
    bar = '\t0% ['
    for i in range(1,11):
        if i <= frac:
            bar = bar+'='
        else:
            bar = bar+'.'
    frac = round((count/total)*100, 2)
    bar = bar+'] 100%\t'+str(frac)+'%\t'+message+'\r'
    sys.stdout.write(bar)

def getCounts(seqinfo, forward=True):
    """
    Returns the coverage, number of mutations for each nucleotide and total number of mutations for a given read base column of a line from a pileup file.
    
    Example:
        >>> s = 'chrI\t9605\tT\t11\tCC,..CCC,CC\t79H@I898H88'
        >>> s = s.split('\\t')[4]
        >>> s
        'CC,..CCC,CC'
        >>> getCounts(s)
        [9, {'T': 0, 'C': 7, 'A': 0, 'G': 0}, 7]
        >>> getCounts(s, forward=False)
        [2, {'T': 0, 'C': 0, 'A': 0, 'G': 0}, 0]
        
    
    Args:
        seqinfo (str): Read base column of a pileup line
        forward (bool): Defines if read and mutations counts are returned for the forward strand (default) or reverse strand
    Returns:
        list : 3-element list where [0] is the total coverage, [1] is a dictionary for the four bases {A,C,G,T} where each entry returns the number of observed mutations for that nucleotide, [2] total number of mutations
    """
    coverage  = 0
    totalmut  = 0
    mutcounts = {'A':0, 'C':0, 'G':0, 'T':0}
    if forward:
        for c in seqinfo:
            if c == '.':
                coverage +=1
            if c == 'A':
                mutcounts['A'] += 1
                coverage +=1
                totalmut += 1
            if c == 'C':
                mutcounts['C'] += 1
                coverage +=1
                totalmut += 1
            if c == 'G':
                mutcounts['G'] += 1
                coverage +=1
                totalmut += 1
            if c == 'T':
                mutcounts['T'] += 1
                coverage +=1
                totalmut += 1
    else:
        for c in seqinfo:
            if c == ',':
                coverage +=1
            if c == 'a':
                mutcounts['A'] += 1
                coverage +=1
                totalmut += 1
            if c == 'c':
                mutcounts['C'] += 1
                coverage +=1
                totalmut += 1
            if c == 'g':
                mutcounts['G'] += 1
                coverage +=1
                totalmut += 1
            if c == 't':
                mutcounts['T'] += 1
                coverage +=1
                totalmut += 1
    return([coverage, mutcounts, totalmut])

#includes anti sense stuff. selects every site within a annotation and 200upstream antisesne of tss and 200 upstream of pA
def selectSitesOverlappingWithGFF(sites, gff):
    chrs   = []
    pos    = []
    m      = []
    r      = []
    result = []
    strand = []
    occ    = []
    count_overlaps = 0
    chrstarts = {}
    curchr = gff['chr'][0]
    chrstarts[curchr] = 0
    for i in range(len(gff['chr'])):
        if gff['chr'][i] != curchr:
            curchr = gff['chr'][i]
            chrstarts[curchr] = i
    for i in range(len(sites['chr'])):
        if i%100 == 0:
            showProgress(i, len(sites['chr']), 'Filter progress: ')
        take = False
        try:
            j = chrstarts[sites['chr'][i]]
            run = True
        except:
            run = False
        while run:
#        for j in range(len(gff['chr'])):
            if gff['strand'][j] == sites['strand'][i] and sites['chr'][i] == gff['chr'][j] and sites['pos'][i] >= gff['start'][j] and sites['pos'][i] <= gff['stop'][j]:
                take = True
            #get antisense sites
            if gff['strand'][j] != sites['strand'][i] and sites['strand'][i] == '+' and sites['chr'][i] == gff['chr'][j] and sites['pos'][i] > (gff['start'][j]-200) and sites['pos'][i] < gff['start'][j]:
                take = True
            if gff['strand'][j] != sites['strand'][i] and sites['strand'][i] == '+' and sites['chr'][i] == gff['chr'][j] and sites['pos'][i] > (gff['stop'][j]-200) and sites['pos'][i] < gff['stop'][j]:
                take = True
            if gff['strand'][j] != sites['strand'][i] and sites['strand'][i] == '-' and sites['chr'][i] == gff['chr'][j] and sites['pos'][i] > gff['start'][j] and sites['pos'][i] < (gff['start'][j]+200):
                take = True
            if gff['strand'][j] != sites['strand'][i] and sites['strand'][i] == '-' and sites['chr'][i] == gff['chr'][j] and sites['pos'][i] > gff['stop'][j] and sites['pos'][i] < (gff['stop'][j]+200):
                take = True
            if take:
                count_overlaps += 1
                chrs.append(sites['chr'][i])
                pos.append(sites['pos'][i])
                m.append(sites['m'][i])
                r.append(sites['r'][i])
                result.append(sites['result'][i])
                strand.append(sites['strand'][i])
                occ.append(sites['occ'][i])
                break
            j += 1
            if j > len(gff['chr'])-2:
                break
            if gff['chr'][j] != sites['chr'][i]:
                break
    print(str(count_overlaps)+'/'+str(len(sites['chr'])))
    return({'chr':chrs, 'pos':pos, 'm':m, 'r':r, 'result':result, 'strand':strand, 'occ':occ})


#fehlt: nur sites die innerhalb eines GFFs sind, die um +/- left and right von start oder stop aus sind, 
#skript das alle overlapping sites entfernt und jeweils nur die stärkste site behält

##################################################################
######### Data handling for sequence data
##################################################################

##Reads a single fasta file and returns the sequence
#def readFasta(filename):
#    """
#    Returns the a single string, containing the sequence information of *filename*.
#    
#    Args:
#        filename (str): Input \*.fasta
#    Returns:
#        str : Content of *filename*
#    """
#    fc = open(filename, 'r')
#    line = fc.readline()
#    line = fc.readline()
#    seq = ''
#    while(line):
#        line = line.replace('\n', '')
#        seq = seq+line.upper()
#        line = fc.readline()
#    
#    fc.close()
#    return(seq)

##Reads all Fasta-Files from a directory
#def buildGenome(directory, verbose, chridentifier=None):
#    """
#    Reads all \*.fasta files of *directory* and returns a dictionary where each key is equal to a filename of *directory* and the value of a key:value pair contains the sequence returned by readFasta().
#    
#    Args:
#        directory (str): Directory which ONLY contains \*.fasta files.
#        verbose (bool): Verbose output
#    Returns:
#        dict : 
#    """
#    genome = {}
#    filenames = os.listdir(directory)
#    count = 0
#    for f in filenames:
#        count += 1
#        if chridentifier == None:
#            genome.update({f.split('.')[0] : readFasta(directory+f)})
#        else:
#            if f.split('.')[0] == chridentifier:
#                genome.update({f.split('.')[0] : readFasta(directory+f)})
#        if verbose and chridentifier == None: showProgress(count, len(filenames),'Reading Genome')
#    if verbose: print('')
#    return(genome)

def makeReverseComplement(s):
    """
    Returns the reverse complement of a nucleotide sequence
    
    Example:
        >>> makeReverseComplement('AACGTAGGCCT')
        'AGGCCTACGTT'
    
    Args:
        s (str): Nucleotide sequence ({A,C,G,T,-,N} are allowed)
    Returns:
        str : Reverse complement of s.
    """
    s_rev = s[::-1]
    rc = {'A':'T', 'G':'C', 'C':'G', 'T':'A', '-':'-','N':'N'}
    s_revcomp = ''
    for c in s_rev:
        s_revcomp = s_revcomp+rc[c]
    return(s_revcomp)

def findAllSubstrings(s, sub):
    """
    Returns a list of all start positions of the string *sub* in the string *s*.
    
    .. note:: In contrast to the regex used by python matches can overlap
    
    Example:
        >>> findAllSubstrings('TTTAATTTTATTT', 'TTT')
        [0, 5, 6, 10]
    
    Args:
        s (str): String you like to search in
        sub (str): Substring you like to find in *s*
    Returns:
        list : List of all start positions of *sub* in *s*. If *sub* does not occur in *s* an empty list is returned.
    """
    startpos = []
    
    for i in range(0,(len(s)-len(sub)+1)):
        if s[i:(i+len(sub))] == sub:
            startpos.append(i)
    return(startpos)

def makekmers(k,alphabet):
    """
    Returns a nested list of k-mer strings of an alphabet and a given number k.
    
    Example:
        >>> makekmers(2,['A','C','G','T'])
        [['A', 'C', 'G', 'T'], 
        ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']]
    
    Args:
        k (int): Defines the maximum length of a word
        alphabet (list): List of characters from which the k-mers are build
    Returns:
        list : Nested list of k-mers so that [0] returns the given alphabet, [1] all possible 2-mers, ..., [k-1] all words with length k
    """
    kmers = []
    kmers.append(alphabet)
    l = 2
    while l <= k:
        currentKmer = []
        for a in alphabet:
            for b in kmers[l-2]:
                currentKmer.append(a+b)
        kmers.append(currentKmer)
        l = l + 1
    return(kmers)

#die richtigen Sequenzen um annotationen aussuchen.

#################################################################
######## GFF data handling 
#################################################################

#gibt die Distanz zwischen zwei Annotationen auf einem Strang zurueck
def getDistanceToNextEntryASense(data, i):
    dist = 0
    run = True
    if data['strand'][i] == '+':
        j = i + 1
        while run:
            if j > (len(data['chr'])-1):
                run = False
                break
            if data['chr'][i] != data['chr'][j]:
                run = False
                break
            if data['strand'][i] != data['strand'][j] and data['chr'][i] == data['chr'][j]:
                dist = data['start'][j]-data['stop'][i]
                run = False
            j += 1
    else:
        j = i - 1
        while run:
            if j < 0:
                run = False
                break
            if data['chr'][i] != data['chr'][j]:
                run = False
                break
            if data['strand'][i] != data['strand'][j] and data['chr'][i] == data['chr'][j]:
                dist = data['start'][i]-data['stop'][j]
                run = False
            j -= 1
    return(dist)

def filterGFF_MAX(gff, maxdist=250, sense = True):
    if sense:
        d      = getDistanceToNextEntrySense(gff)
    else:
        d      = getDistanceToNextEntryASense(gff)
    chrs   = []
    start  = []
    stop   = []
    strand = []
    info   = []
    typeof = []
    count = 0
    for i in range(len(gff['chr'])):
        if d[i] > 0  and d[i] <= maxdist:
            count += 1
            chrs.append(gff['chr'][i])
            start.append(gff['start'][i])
            stop.append(gff['stop'][i])
            strand.append(gff['strand'][i])
            info.append(gff['info'][i])
            typeof.append(gff['type'][i])
    print(str(count)+' instances found')
    return({'chr':chrs,'start':start,'stop':stop,'strand':strand, 'info':info, 'type':typeof})

def filterGFF_MIN(gff, mindist=250, sense = True):
    if sense:
        d      = getDistanceToNextEntrySense(gff)
    else:
        d      = getDistanceToNextEntryASense(gff)
    chrs   = []
    start  = []
    stop   = []
    strand = []
    info   = []
    typeof = []
    count = 0
    for i in range(len(gff['chr'])):
        if d[i] >= mindist:
            count += 1
            chrs.append(gff['chr'][i])
            start.append(gff['start'][i])
            stop.append(gff['stop'][i])
            strand.append(gff['strand'][i])
            info.append(gff['info'][i])
            typeof.append(gff['type'][i])
    print(str(count)+' instances found')
    return({'chr':chrs,'start':start,'stop':stop,'strand':strand, 'info':info, 'type':typeof})

#sorts a GFF in a way that the the startpositions per chromosome are consecutive. Chrnames are sorted according to build in sort
def sortGFF(gff):
    chrnames = set(gff['chr'])
    for c in chrnames:
        tmp_start = 99999999
        #first find the first entry of a chromosome
        for i in range(len(gff['chr'])):
            if gff['chr'][i] == c:
                if gff['start'][i] < tmp_start:
                    tmp_start = gff['start'][i]
        

#The ugliest memory consumption method you ever saw!
def getInfo(pid, identifier='VmSize'):
    fc = open('/proc/'+str(pid)+'/status', 'r')
    line = fc.readline()
    while line:
        if line.split(':')[0] == identifier:
            return(int(line.split('\t')[1].split('kB')[0].strip()))
            break
        line = fc.readline()

def getQuantile(v,q):
    values = copy.deepcopy(v)
    values.sort()
    if int(len(values)*q) < len(values):
        return(values[int(len(values)*q)])
    else:
        return(values[len(values)-1])

def getQuantileIndex(length,q):
    if int(length*q) < length:
        return(int(length*q))
    else:
        return(length-1)




















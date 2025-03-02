class TELOMERE:
    def __init__(self, fileName, binSize):
        self.binSize = binSize

        fin = open(fileName)

        self.pCount_DICT = {}
        self.mCount_DICT = {}

        for line in fin:
            seqName, pos, strand = line.rstrip('\n').split('\t')

            pos = int(pos)

            binIDX = int(pos / self.binSize)

            if strand == '+':
                if not seqName in self.pCount_DICT: self.pCount_DICT[seqName] = {}
                if not binIDX in self.pCount_DICT[seqName]: self.pCount_DICT[seqName][binIDX] = 0
                self.pCount_DICT[seqName][binIDX] += 1
            elif strand == '-':
                if not seqName in self.mCount_DICT: self.mCount_DICT[seqName] = {}
                if not binIDX in self.mCount_DICT[seqName]: self.mCount_DICT[seqName][binIDX] = 0
                self.mCount_DICT[seqName][binIDX] += 1
            else:
                print('bug!')
    def get(self, seqName, strand, sPos, ePos):
        if strand == '+':
            count_DICT = self.pCount_DICT
        elif strand == '-':
            count_DICT = self.mCount_DICT
        else:
            print('bug')

        if not seqName in count_DICT: return 0

        bin_DICT = count_DICT[seqName]

        sbinIDX = int(sPos / self.binSize)
        ebinIDX = int(ePos / self.binSize)

        count = 0
        for binIDX in range(sbinIDX, ebinIDX + 1):
            if not binIDX in bin_DICT: continue
            count += bin_DICT[binIDX]

        return count


prefix = 'Sorghum_bicolor'

windowSize = 10000
slideSize = 1000


#####################################################################################
chromIDX_DICT = {}

fin = open(prefix + '.order')
for chromIDX, line in enumerate(fin):
    chromName = line.rstrip('\n').split('\t')[0]
    chromIDX_DICT[chromName] = chromIDX
fin.close()
#####################################################################################
chromSize_DICT = {}

fin = open(prefix + '.fa.fai')
for lienIDX, line in enumerate(fin):
    chromName,  chromSize = line.rstrip('\n').split('\t')[0:2]
    chromSize = int(chromSize)
    chromSize_DICT[chromName] = chromSize
fin.close()
#####################################################################################




telomere = TELOMERE(prefix + '.telomere', slideSize)

pCount = telomere.get('scaffold_17', '+', 0, 0 + windowSize)
mCount = telomere.get('scaffold_17', '-', 0, 0 + windowSize)

print(pCount, mCount)


telomere = TELOMERE(prefix + '.telomere', slideSize)
for chromName, chromIDX in chromIDX_DICT.items():
    chromSize = chromSize_DICT[chromName]
    #if chromName != 'scaffold_1': continue
    print(chromName)

    pMax = 0
    mMax = 0


    sPos = 0
    while True:
        pCount = telomere.get(chromName, '+', sPos, sPos + windowSize)
        mCount = telomere.get(chromName, '-', sPos, sPos + windowSize)

        pMax = max(pMax, pCount)
        mMax = max(mMax, mCount)

        sPos += slideSize
        if sPos > chromSize: break
    
    print(chromName, pMax, mMax)

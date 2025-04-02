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
#####################################################################################

from optparse import OptionParser
import sys
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
(opt, args) = parser.parse_args()
if opt.prefix == None:
    print("Error! usage: 04.teloCounter.py -p keumgang")
    sys.exit()

prefix = opt.prefix

chrSize_Ratio = 0.0000005
chromWidth = 30
chromGap = 10

#telomere
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

fout = open(prefix + '.Count', 'w')
for chromName, chromIDX in chromIDX_DICT.items():
    max_pCount = 0
    max_pPos = 0
    max_mCount = 0
    max_mPos = 0

    chromSize = chromSize_DICT[chromName]
    #if chromName != 'scaffold_1': continue
    print(chromName)

    sPos = 0
    while True:
        pCount = telomere.get(chromName, '+', sPos, sPos + windowSize)
        mCount = telomere.get(chromName, '-', sPos, sPos + windowSize)

        if max_pCount < pCount:
            max_pCount = pCount
            max_pPos = sPos
        if max_mCount < mCount:
            max_mCount = mCount
            max_mPos = sPos + windowSize

        sPos += slideSize

        if sPos > chromSize: break
    fout.write('\t'.join(map(str, [chromName, chromSize, max_pCount, max_mCount, max_pPos, max_mPos])) + '\n')
fout.close()

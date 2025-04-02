class element:
    def __init__(self,name, container):
        self.name = name
        self.attr_DICT = {}
        self.attrKey_LIST = []
        self.style_DICT = {}
        self.styleKey_LIST = []
        self.item_LIST = []

        if container != None:
            container.add(self)
        self.default()
    def default(self):
        self.style('box-sizing', 'border-box')
        self.style('float', 'left')
        self.style('font-family', 'Times New Roman')

    def add(self,item):
        self.item_LIST += [item]
    def __getitem__(self,key):
        return self.attr_DICT[key]
    def attr(self,key,value):
        if key in self.attr_DICT:
            self.attr_DICT[key] = value
        else:
            self.attr_DICT[key] = value
            self.attrKey_LIST += [key]
        return self
    def style(self,key,value):
        if key in self.style_DICT:
            self.style_DICT[key] = value
        else:
            self.style_DICT[key] = value
            self.styleKey_LIST += [key]
        return self
    def __repr__(self):
        tag = '<'+ self.name
        for attrKey in self.attrKey_LIST:
            tag += ' ' + attrKey + '=' + '"' + str(self.attr_DICT[attrKey]) + '"'
        style = ''
        for styleKey in self.styleKey_LIST:
            style += styleKey + ':' + str(self.style_DICT[styleKey]) + ';'
        if style != '':
            tag += ' ' + 'style' + '=' + '"'  + style + '"'
        tag += '>'

        for item in self.item_LIST:
            tag += str(item)
        tag += '</' + self.name + '>'
        return tag

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
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
(opt, args) = parser.parse_args()
if opt.prefix == None:
    print("Error! usage: 03.draw.py -p keumgang")
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

svg = element('svg', None)
svg.attr('viewBox', "0 0 1000 1000")


for chromName, chromIDX in chromIDX_DICT.items():
    chromSize = chromSize_DICT[chromName]

    chromosome = element('rect', svg)
    chromosome.attr('x', (chromGap + chromWidth)*chromIDX)
    chromosome.attr('y', '10')
    chromosome.attr('width', chromWidth)
    chromosome.attr('height', chromSize * chrSize_Ratio)
    chromosome.attr('rx', '0')
    chromosome.attr('fill', 'rgb(240,240,240)')
    chromosome.attr('stroke-width','0')
    chromosome.attr('stroke', 'gray')

    line = element('rect', svg)
    line.attr('x', (chromGap + chromWidth)*chromIDX + chromWidth/2)
    line.attr('y', '10')
    line.attr('width', 0.1)
    line.attr('height', chromSize * chrSize_Ratio)
    line.attr('rx', '0')
    line.attr('fill', 'black')
    line.attr('stroke-width','0')
    line.attr('stroke', 'gray')

fin = open(prefix + '.gap')
for line in fin:
    chromName, gapLength, sPos, ePos = line.rstrip('\n').split('\t')
    sPos = int(sPos)
    ePos = int(ePos)

    if not chromName in chromIDX_DICT: continue
    chromIDX = chromIDX_DICT[chromName]
    
    gap = element('rect', svg)

    gap.attr('x', (chromGap + chromWidth)*chromIDX)
    gap.attr('y', round(10 + sPos * chrSize_Ratio, 3))
    gap.attr('width', chromWidth)
    #gap.attr('height', round((ePos - sPos + 1) * chrSize_Ratio, 3))
    gap.attr('height', 0.1)
    gap.attr('fill', 'gray')
    gap.attr('opacity', '1')
    gap.attr('stroke-width','0')
    gap.attr('stroke', 'gray')



telomere = TELOMERE(prefix + '.telomere', slideSize)
for chromName, chromIDX in chromIDX_DICT.items():
    chromSize = chromSize_DICT[chromName]
    #if chromName != 'scaffold_1': continue
    print(chromName)

    sPos = 0
    while True:
        pCount = telomere.get(chromName, '+', sPos, sPos + windowSize)
        mCount = telomere.get(chromName, '-', sPos, sPos + windowSize)

        if pCount > 10:
            gap = element('rect', svg)
            gap.attr('x', (chromGap + chromWidth)*chromIDX + chromWidth/2)
            gap.attr('y', round(10 + sPos * chrSize_Ratio, 3))
            gap.attr('width', pCount/100)
            gap.attr('height', 0.5)
            gap.attr('fill', 'red')
            gap.attr('opacity', '1')

        if mCount > 10:
            gap = element('rect', svg)
            gap.attr('x', (chromGap + chromWidth)*chromIDX + chromWidth/2 - mCount/100)
            gap.attr('y', round(10 + sPos * chrSize_Ratio, 3))
            gap.attr('width', mCount/100)
            gap.attr('height', 0.5)
            gap.attr('fill', 'blue')
            gap.attr('opacity', '1')

        sPos += slideSize

        if sPos > chromSize: break


fout = open(prefix + '.html', 'w')
fout.write(str(svg))
fout.close()

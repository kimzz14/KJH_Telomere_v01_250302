class MFA:
    def __init__(self, fileName):
        self.seqName_LIST = []
        self.seq_DICT = {}

        fin = open(fileName)
        for line in fin:
            if line.startswith('>') == True:
                seqName = line.rstrip('\n').split('\t')[0].split(' ')[0][1:]
                self.seqName_LIST += [seqName]
                self.seq_DICT[seqName] = []
            else:
                sequence = line.rstrip('\n').upper()
                self.seq_DICT[seqName] += [sequence]
        fin.close()

        for seqName in self.seqName_LIST:
            self.seq_DICT[seqName] = ''.join(self.seq_DICT[seqName])
        
        print('[done] mfa file read - ' + fileName, flush=True)

seqLen = 0
gapLen = 0

prefix = 'keumgang'
fileName = prefix + '.fa'

fin = open(fileName)
for line in fin:
    if line.startswith('>') == True:
        pass
    else:
        sequence = line.rstrip('\n').upper()
        seqLen += len(sequence)
        gapLen += sequence.count('N')
fin.close()

print(seqLen)
print(gapLen)
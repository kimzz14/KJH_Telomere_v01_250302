from optparse import OptionParser
import sys
from joblib import Parallel, delayed
#option parser
parser = OptionParser(usage="""Run annotation.py \n Usage: %prog [options]""")
parser.add_option("-p","--prefix",action = 'store',type = 'string',dest = 'prefix',help = "")
parser.add_option("-t","--threadN",action = 'store',type = 'int',dest = 'threadN',help = "")
(opt, args) = parser.parse_args()
if opt.threadN == None or opt.prefix == None:
    print("Error! usage: 02.find_telomere.py -t 24 -p keumgang")
    sys.exit()

batchN = opt.threadN
prefix = opt.prefix
#prefix = 'keumgang'

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
    
def run_batch(seq_DICT, batchN):
    #motif = 'AAACCCT'
    #motif_rev = 'AGGGTTT'
    motif = 'AACCCT'
    motif_rev = 'AGGGTT'

    motifLen = len(motif)

    def run_single(seqName, sequence):
        result_LIST = []
        seqLen = len(sequence)
        for idx in range(seqLen - (motifLen - 1)):
            subSequence = sequence[idx:idx + motifLen]
            if subSequence == motif:
                result_LIST += [(seqName, idx + 1, '+')]
            elif subSequence == motif_rev:
                result_LIST += [(seqName, idx + 1, '-')]

        print('[done] find telomere - ' + seqName, flush=True)
        return (seqName, result_LIST)

    batchResult = Parallel(n_jobs=batchN)(delayed(run_single)(seqName, sequence) for seqName, sequence in seq_DICT.items())

    result_DICT = {}
    for seqName, result_LIST in batchResult:
        result_DICT[seqName] = result_LIST

    return result_DICT

ref = MFA(prefix + '.fa')
result_DICT = run_batch(ref.seq_DICT, batchN)

fout = open(prefix + '.telomere', 'w')
for seqName in ref.seqName_LIST:
    for result in result_DICT[seqName]:
        fout.write('\t'.join(map(str, result)) + '\n')
fout.close()

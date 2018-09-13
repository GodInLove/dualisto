import sys


class TranscriptInfo:
    def __init__(self):
        self.tran_num = -1
        self.tran_name = []
        self.tran_len = []


class KmerInfo:
    def __init__(self):
        self.contig_id = -1
        self.pos_in_contig = -1
        self.n_of_kmer_in_contig = -1
        self.sense_in_contig = -1


class ContigIncludeTrans:
    def __init__(self):
        self.tran_id = -1
        self.pos_in_tran = -1
        self.sense_in_tran = -1


class Contig:
    def __init__(self):
        self.id = -1
        self.n_of_kmer = -1
        self.ecs_id = -1
        self.seq = ""
        self.include_trans = []


class ContigFindTrans:
    def __init__(self):
        self.tran_id = -1
        self.start_in_tran = -1
        self.stop_in_tran = -1
        self.sense_in_tran = -1


def base_pair(base):
    if base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"
    elif base == "T":
        return "A"
    elif base == "N":
        return "N"


def str_pair(_str):
    i = 0
    while i < len(_str):
        _str[i] = base_pair(_str[i])
    return _str


def kmer_twin(dna_str):
    dna_str = reversed(dna_str)
    new_dna_str = ""
    for base in dna_str:
        print(base,end="")
        new_dna_str = new_dna_str + base_pair(base)
    print()
    return new_dna_str


def kmer_rep(dna_str):
    return min(dna_str, kmer_twin(dna_str))


def contig_bw(kmer):
    for x in 'ACGT':
        yield x + kmer[:-1]


def contig_fw(kmer):
    for x in 'ACGT':
        yield kmer[1:] + x


def program_stop(info):
    print(info)
    sys.exit(1)


def _int32(x):
    return int(0xFFFFFFFF & x)


class MT19937:
    def __init__(self, seed):
        self.mt = [0] * 624
        self.mt[0] = seed
        for i in range(1, 624):
            self.mt[i] = _int32(1812433253 * (self.mt[i - 1] ^ self.mt[i - 1] >> 30) + i)

    def extract_number(self):
        self.twist()
        y = self.mt[0]
        y = y ^ y >> 11
        y = y ^ y << 7 & 2636928640
        y = y ^ y << 15 & 4022730752
        y = y ^ y >> 18
        return _int32(y)

    def twist(self):
        for i in range(0, 624):
            y = _int32((self.mt[i] & 0x80000000) + (self.mt[(i + 1) % 624] & 0x7fffffff))
            self.mt[i] = y ^ self.mt[(i + 397) % 624] >> 1

            if y % 2 != 0:
                self.mt[i] = self.mt[i] ^ 0x9908b0df


def dna_base(num):
    dna = "ACGT"
    return dna[num & 0x03]


class DefArgs:
    def __init__(self):
        self.k = 31
        self.threads = 1
        self.index = "dualisto.idx"
        self.output = "dualisto_out"
        self.seed = 42
        self.fna_file = "dualisto_fna_file"
        self.fq_file_1 = "fq_file_1"
        self.fq_file_2 = "fq_file_2"
        self.single_mode = False
        self.fq_file_single = "fq_file_single"
        self.len_frag = 0.0
        self.sd = 0.0


transcript_info = TranscriptInfo()
ec_map = []
ecs_list = []
ec_inv_dict = {}
kmer_str_dict = {}
kmer_str_list = []
kmer_info_list = []
contig_list = []
skip = 1
counts = []

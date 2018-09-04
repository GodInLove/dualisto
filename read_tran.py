import gzip
from Bio import SeqIO

from _class import transcript_info, MT19937, dna_base, program_stop


def read_tran_files(file_name):
    seqs = []
    try:
        if ".gz" in file_name:
            fp = gzip.open(file_name, 'rt')
        else:
            fp = open(file_name, 'r')
    except IOError:
        print(file_name, "is not exist!")
        program_stop("read_tran.py")
    for seq_record in SeqIO.parse(fp, "fasta"):
        transcript_info.tran_name.append(seq_record.id)
        transcript_info.tran_len.append(len(seq_record.seq))
        seq = replace_base(seq_record.seq)
        seq = check_poly_a_tail(seq)
        seqs.append(seq)
    transcript_info.tran_num = len(transcript_info.tran_name)
    fp.close()
    return seqs


def replace_base(seq):
    seq_len = len(seq)
    i = 0
    seq = seq.upper()
    while i < seq_len:
        if seq[i] == 'U':
            seq[i] = 'T'
        elif seq[i] != 'A' and seq[i] != 'T' and seq[i] != 'C' and seq[i] != 'G':
            seq[i] = dna_base(MT19937(42).extract_number())
        i = i + 1
    return seq


def check_poly_a_tail(seq):
    if seq[-10:] == "AAAAAAAAAA":
        while seq[-1] == 'A':
            seq = seq[:-1]
    return seq


if __name__ == '__main__':
    test_seqs = read_tran_files("/home/lyd/PycharmProjects/dualisto/test_files/transcripts.fasta")
    print(test_seqs)

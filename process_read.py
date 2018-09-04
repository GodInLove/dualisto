import gzip
import datetime
from multiprocessing import Pool

import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import FastaIterator
from Bio.Seq import Seq

from _class import program_stop, counts, transcript_info
from em import em_run, write_em_tsv
from map import match, intersect_ecs_paired, map_pair, find_ec, compute_mean_flg_lens, get_each_tran_len

k = 0
new_ecs = []
f_lens = [0] * 1000
f_len_goal = 10000


def do_something_with_record(info):
    name, seq = info
    return name, len(seq)


def do_something_with_record_paired(info1, info2):
    name1, seq1 = info1
    name2, seq2 = info2
    v1 = match(seq1, len(seq1), k)
    v2 = match(seq2, len(seq2), k)
    u = intersect_ecs_paired(v1, v2)
    if len(v1) != 0 and len(v2) != 0:
        tl = map_pair(seq1, len(seq1), seq2, len(seq2), k)
    else:
        tl = -1
    return u, tl


def convert_to_fasta(in_handle):
    for rec_id, seq, _ in FastqGeneralIterator(in_handle):
        yield rec_id, str(Seq(seq))


def to_fasta(in_handle):
    for seq_record in FastaIterator(in_handle):
        yield seq_record.id, str(seq_record.seq)


def read_fq_fa_files(args):
    global k, f_len_goal
    k = args.k
    pool = Pool(processes=args.threads)
    begin = datetime.datetime.now()
    if args.single_mode:
        file_name = args.fq_file_single
        file_suffix = os.path.splitext(file_name)[1]
        try:
            if ".gz" == file_suffix:
                fp = gzip.open(file_name, 'rt')
                file_suffix_suffix = os.path.splitext(os.path.splitext(file_name)[0])[1]
            else:
                fp = open(file_name, 'r')
                file_suffix_suffix = file_suffix
        except IOError:
            print(file_name, "is not exist!")
            program_stop("process_read.py")
        if ".fasta" == file_suffix_suffix or ".fa" == file_suffix_suffix:
            res = pool.map(do_something_with_record, to_fasta(fp))
            pool.close()
            pool.join()
        elif ".fastq" == file_suffix_suffix or ".fq" == file_suffix_suffix:
            res = pool.map(do_something_with_record, convert_to_fasta(fp))
            pool.close()
            pool.join()
        else:
            print("cannot find .fa(.fasta) or .fq(.fastq).please check.")
            program_stop("process_read.py")
        print(len(res))
        print("please waite single mode.")
    else:
        # paired_mode
        file_name_1 = args.fq_file_1
        file_suffix_1 = os.path.splitext(file_name_1)[1]
        file_name_2 = args.fq_file_2
        file_suffix_2 = os.path.splitext(file_name_2)[1]
        try:
            if ".gz" in file_name_1:
                fp1 = gzip.open(file_name_1, 'rt')
                file_suffix_suffix_1 = os.path.splitext(os.path.splitext(file_name_1)[0])[1]
            else:
                fp1 = open(file_name_1, 'r')
                file_suffix_suffix_1 = file_suffix_1
            if ".gz" in file_name_2:
                fp2 = gzip.open(file_name_2, 'rt')
                file_suffix_suffix_2 = os.path.splitext(os.path.splitext(file_name_2)[0])[1]
            else:
                fp2 = open(file_name_2, 'r')
                file_suffix_suffix_2 = file_suffix_2
        except IOError:
            print(file_name_1, "or", file_name_2, "is not exist!")
            program_stop("process_read.py")
        if file_suffix_suffix_1 != file_suffix_suffix_2:
            print("suffix is not inconsistent.please check.")
            program_stop("process_read.py")
        else:
            if ".fasta" == file_suffix_suffix_1 or ".fa" == file_suffix_suffix_1:
                res = pool.starmap(do_something_with_record_paired, zip(to_fasta(fp1), to_fasta(fp2)))
                pool.close()
                pool.join()
            elif ".fastq" == file_suffix_suffix_1 or ".fq" == file_suffix_suffix_1:
                res = pool.starmap(do_something_with_record_paired, zip(convert_to_fasta(fp1), convert_to_fasta(fp2)))
                pool.close()
                pool.join()
            else:
                print("cannot find .fa(.fasta) or .fq(.fastq).please check.")
                program_stop("process_read.py")
            end = datetime.datetime.now()
            print("process_reads:", end - begin)
            begin = datetime.datetime.now()
            for each_u, each_tl in res:
                ec = find_ec(each_u)
                if ec == -1 or ec >= len(counts):
                    new_ecs.append(each_u)
                else:
                    counts[ec] = counts[ec] + 1
                if f_len_goal > 0 and 0 <= ec < transcript_info.tran_num:
                    if 0 < each_tl < len(f_lens):
                        f_lens[each_tl] = f_lens[each_tl] + 1
                        f_len_goal = f_len_goal - 1
            end = datetime.datetime.now()
            print("match_ecs:", end - begin)
            begin = datetime.datetime.now()
            mean_f_lens = compute_mean_flg_lens(f_lens)
            tran_lens_estimated = get_each_tran_len(mean_f_lens)
            alpha_list, eff_lens = em_run(counts, tran_lens_estimated)
            output_file = os.path.join(args.output, "dualisto_quant.tsv")
            write_em_tsv(output_file, alpha_list, eff_lens)
            end = datetime.datetime.now()
            print("em_run:", end - begin)

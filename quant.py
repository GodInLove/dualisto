import os
import datetime

from _class import program_stop, DefArgs
from index_write_load import load_idx
from process_read import read_fq_fa_files

def_args = DefArgs()


def start(args):
    if not os.path.exists(args.index_path):
        print("the index_file is not exist.please check.")
        program_stop("quant.py")
    else:
        begin = datetime.datetime.now()
        load_k = load_idx(args.index_path)
        end = datetime.datetime.now()
        print("load_idx:", end - begin)
        if args.kmer_len != load_k:
            print("the k value is inconsistent.")
            program_stop("quant.py")
    def_args.k = args.kmer_len
    def_args.index = args.index_path
    if args.single_mode:
        def_args.len_frag = args.fragment_len
        def_args.sd = args.sd_len
        def_args.fq_file_single = args.fa_or_fq_file
        if not os.path.exists(def_args.fq_file_single):
            print("fq_or_fa_files is required.please check.")
            program_stop("quant.py")
    else:
        def_args.fq_file_1 = args.fa_or_fq_files_1
        def_args.fq_file_2 = args.fa_or_fq_files_2
        if not (os.path.exists(def_args.fq_file_1) or os.path.exists(def_args.fq_file_2)):
            print("fq_or_fa_files is required.please check.")
            program_stop("quant.py")
    output_path_list = os.path.split(args.output_path)
    if len(output_path_list[0]) == 0:
        output_path = os.path.join(os.getcwd(), output_path_list[1])
    else:
        output_path = args.output_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    def_args.output = output_path
    def_args.threads = args.threads
    def_args.single_mode = args.single_mode
    read_fq_fa_files(def_args)


def args_handle(dualisto_quant):
    quant_args = dualisto_quant.add_argument_group("quant options")
    quant_args.add_argument("-k", action="store", dest="kmer_len", default=31, type=int,
                            help="kmer_length(default:31,must odd)")
    quant_args.add_argument("-L", action="store", dest="fa_or_fq_files_1", default="null",
                            help="input_fa/fq_files_1")
    quant_args.add_argument("-R", action="store", dest="fa_or_fq_files_2", default="null",
                            help="input_fa/fq_files_2")
    quant_args.add_argument("-i", action="store", dest="index_path", default="null",
                            help="index_path")
    quant_args.add_argument("-O", action="store", dest="output_path", default="dualisto_out",
                            help="output_path")
    quant_args.add_argument("-t", action="store", dest="threads", default=1, type=int,
                            help="threads number")
    quant_advanced_args = quant_args.add_argument_group("advanced options")
    quant_advanced_args.add_argument("-I", action="store", dest="fa_or_fq_file", default="null",
                                     help="input_fa/fq_file_single_mode")
    quant_advanced_args.add_argument("--single", action="store_true", dest="single_mode", default=False,
                                     help="single_mode(required -l,-s)")
    quant_advanced_args.add_argument("-l", action="store", dest="fragment_len", default=0.0, type=float,
                                     help="average fragment length")
    quant_advanced_args.add_argument("-s", action="store", dest="sd_len", default=0.0, type=float,
                                     help="standard deviation of fragment length")

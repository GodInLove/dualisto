import os

from _class import program_stop, DefArgs
from contig import build_dbg
from ecs import build_ecs
from index_write_load import write_idx
from read_tran import read_tran_files

def_args = DefArgs()


def start(args):
    if args.kmer_len % 2 == 0:
        print("k value must be odd.")
        program_stop("index.py")
    if args.input_fna_path == "null":
        print("the input_fna_path is required.")
        program_stop("index.py")
    else:
        if not os.path.isfile(args.input_fna_path):
            print("the input_fna_path is wrong.please check.")
            program_stop("index.py")
    index_out_dir = os.path.split(args.index_path)
    if len(index_out_dir[0]) == 0:
        index_path = os.path.join(os.getcwd(), index_out_dir[1])
    else:
        if not os.path.exists(index_out_dir[0]):
            print("the index_out_DIR is not exist.please check.")
            program_stop("index.py")
        else:
            index_path = args.index_path
    def_args.k = args.kmer_len
    def_args.index = index_path
    def_args.fna_file = args.input_fna_path
    seqs = read_tran_files(def_args.fna_file)
    build_dbg(seqs, def_args.k)
    build_ecs(seqs, def_args.k)
    write_idx(def_args.index, def_args.k)


def args_handle(dualisto_index):
    index_args = dualisto_index.add_argument_group("index options")
    index_args.add_argument("-k", action="store", dest="kmer_len", default=31, type=int,
                            help="kmer_length(default:31,must odd)")
    index_args.add_argument("-I", action="store", dest="input_fna_path", default="null",
                            help="tran_fna_file")
    index_args.add_argument("-i", action="store", dest="index_path", default="null",
                            help="index_path")

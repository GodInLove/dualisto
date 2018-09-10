from _class import program_stop, transcript_info, kmer_str_dict, kmer_info_list, ec_map, contig_list, ecs_list, \
    KmerInfo, Contig, ContigIncludeTrans, ec_inv_dict, counts, kmer_str_list
from version import version, kallisto_index_version


def write_idx(filename, k):
    try:
        fp = open(filename, 'w')
    except IOError:
        print(filename, "is not exist!")
        program_stop("index_write_load.py")
    # 1.write version
    fp.write(version + "\n")
    # 2.write k
    fp.write(str(k) + "\n")
    # 3.write num of trans
    fp.write(str(transcript_info.tran_num) + "\n")
    for i in range(transcript_info.tran_num):
        fp.write(str(transcript_info.tran_len[i]) + "\n")
    # 4.write kmer_str_dict and kmer_info_list
    for key in kmer_str_dict:
        kmer_str_list.append(key)
        kmer_info_list.append(kmer_str_dict[key])
    i = 0
    fp.write(str(len(kmer_str_list)) + "\n")
    while i < len(kmer_str_list):
        fp.write(kmer_str_list[i] + "\n")
        fp.write(str(kmer_info_list[i].contig_id) + "," + str(kmer_info_list[i].pos_in_contig) + "," + str(
            kmer_info_list[i].n_of_kmer_in_contig) + "," + str(kmer_info_list[i].sense_in_contig) + "\n")
        i = i + 1
    # 5.write num of ecs
    fp.write(str(len(ec_map)) + "\n")
    # 6.write each ecs
    i = 0
    while i < len(ec_map):
        # 6.1 write num of each ecs
        fp.write(str(len(ec_map[i])) + "\n")
        # 6.2 write each trans in each ecs
        j = 0
        while j < len(ec_map[i]):
            fp.write(str(ec_map[i][j]) + "\n")
            j = j + 1
        i = i + 1
    # 7.write trans_names
    for i in range(transcript_info.tran_num):
        fp.write(transcript_info.tran_name[i] + "\n")
    # 8.write contigs
    fp.write(str(len(contig_list)) + "\n")
    i = 0
    while i < len(contig_list):
        fp.write(str(contig_list[i].id) + "," + str(contig_list[i].n_of_kmer) + "," + str(contig_list[i].seq) + "\n")
        # 8.1 write contig_to_trans_info
        fp.write(str(len(contig_list[i].include_trans)) + "\n")
        j = 0
        while j < len(contig_list[i].include_trans):
            curr_include_trans = contig_list[i].include_trans[j]
            fp.write(str(curr_include_trans.tran_id) + "," + str(curr_include_trans.pos_in_tran) + "," + str(
                curr_include_trans.sense_in_tran) + "\n")
            j = j + 1
        i = i + 1
    # 9.write ecs info
    i = 0
    while i < len(ecs_list):
        fp.write(str(ecs_list[i]) + "\n")
        i = i + 1
    fp.flush()
    fp.close()


def load_idx(filename):
    try:
        fp = open(filename, 'r')
    except IOError:
        print(filename, "is not exist!")
        program_stop("index_write_load.py")
    # 1.read version
    read_version = fp.readline().strip()
    if read_version != kallisto_index_version:
        program_stop("index_write_load.py")
    # 2.read k
    k = int(fp.readline().strip())
    # 3.read num of trans
    transcript_info.tran_num = int(fp.readline().strip())
    transcript_info.tran_len = [0] * transcript_info.tran_num
    for i in range(transcript_info.tran_num):
        transcript_info.tran_len[i] = int(fp.readline().strip())
    # 4.read kmer_str_dict and kmer_info_list
    kmer_str_list_len = int(fp.readline().strip())
    for i in range(kmer_str_list_len):
        kmer_str_list.append("")
        kmer_info_list.append("")
    i = 0
    while i < kmer_str_list_len:
        kmer_str_list[i] = fp.readline().strip()
        kmer_info = KmerInfo()
        curr_kmer_info = fp.readline().strip().split(",")
        kmer_info.contig_id = int(curr_kmer_info[0])
        kmer_info.pos_in_contig = int(curr_kmer_info[1])
        kmer_info.n_of_kmer_in_contig = int(curr_kmer_info[2])
        kmer_info.sense_in_contig = int(curr_kmer_info[3])
        kmer_info_list[i] = kmer_info
        i = i + 1
    for i in range(kmer_str_list_len):
        kmer_str_dict[kmer_str_list[i]] = kmer_info_list[i]
    # 5.read num of ecs
    ec_map_len = int(fp.readline().strip())
    for i in range(ec_map_len):
        ec_map.append([])
        ec_inv_dict[i] = []
        counts.append(0)
    # 6.read each ecs
    i = 0
    while i < ec_map_len:
        # 6.1 read num of each ecs
        each_ecs_len = int(fp.readline().strip())
        for j in range(each_ecs_len):
            ec_map[i].append(0)
        # 6.2 read each trans in each ecs
        j = 0
        while j < each_ecs_len:
            ec_map[i][j] = int(fp.readline().strip())
            j = j + 1
        ec_inv_dict[i] = ec_map[i]
        i = i + 1
    # 7.read trans_names
    transcript_info.tran_name = [""] * transcript_info.tran_num
    for i in range(transcript_info.tran_num):
        transcript_info.tran_name[i] = fp.readline().strip()
    # 8.read contigs
    contig_list_len = int(fp.readline().strip())
    for i in range(contig_list_len):
        contig_list.append(0)
        ecs_list.append(0)
    i = 0
    while i < contig_list_len:
        load_contig = Contig()
        curr_contig_info = fp.readline().strip().split(",")
        load_contig.id = int(curr_contig_info[0])
        load_contig.n_of_kmer = int(curr_contig_info[1])
        load_contig.seq = curr_contig_info[2]
        # 8.1 read contig_to_trans_info
        curr_contig_include_trans_len = int(fp.readline().strip())
        load_contig.include_trans = [0] * curr_contig_include_trans_len
        j = 0
        while j < curr_contig_include_trans_len:
            load_include_trans = ContigIncludeTrans()
            curr_include_trans = fp.readline().strip().split(",")
            load_include_trans.tran_id = int(curr_include_trans[0])
            load_include_trans.pos_in_tran = int(curr_include_trans[1])
            load_include_trans.sense_in_tran = int(curr_include_trans[2])
            j = j + 1
        i = i + 1
    # 9.read ecs info
    i = 0
    while i < contig_list_len:
        ecs_list[i] = int(fp.readline().strip())
        i = i + 1
    fp.close()
    return k

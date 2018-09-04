from _class import ec_map, kmer_rep, kmer_str_list, program_stop, contig_list, kmer_info_list, \
    ContigFindTrans, kmer_twin, Contig, ecs_list, ContigIncludeTrans, KmerInfo

contig_find_trans_list = []
ec_inv_dict = {}


def build_ecs(seqs, k):
    contig_find_trans_dict = {}
    for i in range(len(contig_list)):
        contig_find_trans_dict[i] = []
    for i in range(len(seqs)):
        each_ec = [i]
        ec_map.append(each_ec)
        ec_inv_dict[i] = each_ec
        seq = seqs[i]
        n_of_kmers_of_seq = len(seq) - k + 1
        j = 0
        while j < n_of_kmers_of_seq:
            curr_kmer = seq[j: j + k]
            curr_kmer = str(curr_kmer)
            curr_rep = kmer_rep(curr_kmer)
            curr_twin = kmer_twin(curr_kmer)
            if curr_rep in kmer_str_list:
                curr_index = kmer_str_list.index(curr_rep)
                curr_contig_id = kmer_info_list[curr_index].contig_id
                curr_contig = contig_list[curr_contig_id]

                kmers_in_contig = []
                for m in range(len(curr_contig.seq) - k + 1):
                    curr_kmer_in_contig = curr_contig.seq[m: m + k]
                    kmers_in_contig.append(curr_kmer_in_contig)
                kmers_in_contig_len = len(kmers_in_contig)

                contig_find_trans = ContigFindTrans()
                contig_find_trans.tran_id = i
                if (curr_kmer == curr_rep) == kmer_info_list[curr_index].sense_in_contig:
                    contig_find_trans.sense_in_tran = 1
                    contig_find_trans.start_in_tran = kmers_in_contig.index(curr_kmer)
                    if kmers_in_contig_len - contig_find_trans.start_in_tran > n_of_kmers_of_seq - j:
                        # 说明这个contig不来自这个seq
                        contig_find_trans.stop_in_tran = contig_find_trans.start_in_tran + n_of_kmers_of_seq - j
                        # [start,stop)这里这个长度肯定是错误的，正好后面要用到
                        j = n_of_kmers_of_seq
                        # 到末尾，也就是结束
                    else:
                        contig_find_trans.stop_in_tran = kmers_in_contig_len
                        # contig属于这个seq，而且长度就是contig的长度
                        j = j + contig_find_trans.stop_in_tran - contig_find_trans.start_in_tran
                        # 下一个从紧接着当前contig的下一个contig开始
                else:
                    contig_find_trans.sense_in_tran = 0
                    contig_find_trans.stop_in_tran = kmers_in_contig.index(curr_twin) + 1
                    # 因为[start,stop) stop代表的是长度，取不到stop，实际取的是stop-1
                    if contig_find_trans.stop_in_tran > n_of_kmers_of_seq - j:
                        contig_find_trans.start_in_tran = contig_find_trans.stop_in_tran - n_of_kmers_of_seq + j
                        j = n_of_kmers_of_seq
                    else:
                        contig_find_trans.start_in_tran = 0
                        j = j + contig_find_trans.stop_in_tran - contig_find_trans.start_in_tran
                contig_find_trans_dict[curr_contig_id].append(contig_find_trans)
            else:
                print(curr_rep)
                print(j)
                program_stop("ecs.py-2")
    for i in range(len(contig_list)):
        contig_find_trans_list.append(contig_find_trans_dict[i])
    fix_contigs(k)
    for curr_contig_id in range(len(contig_find_trans_list)):
        curr_contig_list = []
        for curr_trans_info in contig_find_trans_list[curr_contig_id]:
            curr_contig_list.append(curr_trans_info.tran_id)
        curr_contig_list = sorted(list(set(curr_contig_list)))
        ec = -1
        for (key, value) in ec_inv_dict.items():
            if value == curr_contig_list:
                ec = key
        if ec == -1:
            ec = len(ec_inv_dict)
        ec_inv_dict[ec] = curr_contig_list
        ec_map.append(curr_contig_list)
        ecs_list[curr_contig_id] = ec
        contig_list[curr_contig_id].ecs_id = ec
    for i in range(len(seqs)):
        seq = seqs[i]
        n_of_kmers_of_seq = len(seq) - k + 1
        tmp_str = ""
        j = 0
        while j < n_of_kmers_of_seq:
            curr_kmer = seq[j: j + k]
            curr_kmer = str(curr_kmer)
            curr_rep = kmer_rep(curr_kmer)
            if curr_rep in kmer_str_list:
                curr_index = kmer_str_list.index(curr_rep)
                curr_kmer_info = kmer_info_list[curr_index]
                contig_include_trans = ContigIncludeTrans()
                contig_include_trans.tran_id = i
                contig_include_trans.pos_in_tran = j
                contig_include_trans.sense_in_tran = int(((curr_kmer == curr_rep) == curr_kmer_info.sense_in_contig))
                contig_list[curr_kmer_info.contig_id].include_trans.append(contig_include_trans)
                j = j + contig_list[curr_kmer_info.contig_id].n_of_kmer
                # if contig_include_trans.sense_in_tran:
                #     if contig_include_trans.pos_in_tran == 0:
                #         tmp_str = tmp_str + contig_list[curr_kmer_info.contig_id].seq
                #     else:
                #         tmp_str = tmp_str + contig_list[curr_kmer_info.contig_id].seq[k - 1]
                # else:
                #     new_tm_str = str_pair(contig_list[curr_kmer_info.contig_id].seq)
                #     if contig_include_trans.pos_in_tran == 0:
                #         tmp_str = tmp_str + new_tm_str
                #     else:
                #         tmp_str = tmp_str + new_tm_str[k - 1]
            else:
                program_stop("ecs.py-3")


def fix_contigs(k):
    perfect_contig_num = 0
    for contig_id in range(len(contig_find_trans_list)):
        all_ok = 1
        # 如前文所言，没有找错的
        n_of_kmer_in_curr_contig = contig_list[contig_id].n_of_kmer
        for each_trans_info in contig_find_trans_list[contig_id]:
            if each_trans_info.start_in_tran != 0 or each_trans_info.stop_in_tran != n_of_kmer_in_curr_contig:
                all_ok = 0
                # 说明这里有找错的
        if all_ok:
            perfect_contig_num = perfect_contig_num + 1
        else:
            tmp_list = []
            for each_trans_info in contig_find_trans_list[contig_id]:
                tmp_list.append(each_trans_info.start_in_tran)
                tmp_list.append(each_trans_info.stop_in_tran)
            break_points = sorted(list(set(tmp_list)))
            old_contig_seq = contig_list[contig_id].seq
            old_trans_info = contig_find_trans_list[contig_id]
            new_contig = Contig()
            for j in range(len(break_points) - 1):
                new_contig.seq = old_contig_seq[break_points[j]:break_points[j + 1] - break_points[j] + k - 1]
                new_contig.n_of_kmer = break_points[j + 1] - break_points[j]
                if j == 0:
                    new_contig.id = contig_id
                    contig_list[contig_id] = new_contig
                else:
                    new_contig.id = len(contig_list)
                    contig_list.append(new_contig)
                    ecs_list.append(-1)
                for m in range(len(new_contig.seq) - k + 1):
                    curr_kmer_in_contig = new_contig.seq[m: m + k]
                    if curr_kmer_in_contig in kmer_str_list:
                        curr_kmer_index = kmer_str_list.index(curr_kmer_in_contig)
                        kmer_info_list[curr_kmer_index].contig_id = new_contig.id
                        kmer_info_list[curr_kmer_index].pos_in_contig = m
                        kmer_info_list[curr_kmer_index].n_of_kmer_in_contig = new_contig.n_of_kmer
                    else:
                        program_stop("ecs.py-4")
                new_trans_info_list = []
                for each_trans_info in old_trans_info:
                    if not (each_trans_info.start_in_tran >= break_points[j + 1] or
                            each_trans_info.stop_in_tran <= break_points[j]):
                        new_trans_info = ContigFindTrans()
                        new_trans_info.sense_in_trans = each_trans_info.sense_in_trans
                        new_trans_info.tran_id = each_trans_info.tran_id
                        new_trans_info.start_in_tran = 0
                        new_trans_info.stop_in_tran = new_contig.n_of_kmer
                        new_trans_info_list.append(new_trans_info)
                if j == 0:
                    contig_find_trans_list[contig_id] = new_trans_info_list
                else:
                    contig_find_trans_list.append(new_trans_info_list)

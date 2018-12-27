from _class import kmer_rep, KmerInfo, kmer_str_dict, kmer_twin, contig_fw, contig_bw, contig_list, \
    Contig, ecs_list, program_stop

fw_step_kmer = ""


def build_dbg(seqs, k):
    global fw_step_kmer
    tmp_kmer_map = set()
    # 把序列打断,得到kmer的集合
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            curr_kmer = seq[i: i + k]
            curr_kmer = str(curr_kmer)
            curr_rep = kmer_rep(curr_kmer)
            tmp_kmer_map.add(curr_rep)

    # 给每一个kmer带上一个可以存储它信息的类
    for each_kmer in tmp_kmer_map:
        kmer_info = KmerInfo()
        kmer_str_dict[each_kmer] = kmer_info

    # 把kmer连起来,找contig
    for kmer_key in kmer_str_dict:
        curr_kmer = kmer_key
        curr_kmer_info = kmer_str_dict[kmer_key]
        if curr_kmer_info.contig_id == -1:
            self_loop = 0
            curr_twin = kmer_twin(curr_kmer)
            fw_list = [curr_kmer]
            last_kmer = curr_kmer
            fw_step_kmer = curr_kmer
            while fwstep(fw_step_kmer):
                if fw_step_kmer == curr_kmer:
                    # 该if判断是否成环
                    self_loop = 1
                    break
                    # pass
                    # example(3,"ACTGAC") (5,"ACCAACCA") (5,"TCTGTCTG") (5,"AACAAACA") (5,"CACACA") (5,"ACACAC")
                    # print("begin:",begin,"fw:", fw, "t_kmer", self.tmp_kmer, "seqs:", self.seqs)
                elif fw_step_kmer == curr_twin:
                    self_loop = (len(fw_list) > 1)
                    break
                elif fw_step_kmer == kmer_twin(last_kmer):
                    break
                fw_list.append(fw_step_kmer)
                last_kmer = fw_step_kmer
            fw_step_kmer = curr_twin
            bw_list = []
            first_kmer = curr_twin
            if not self_loop:
                while fwstep(fw_step_kmer):
                    if fw_step_kmer == curr_twin:
                        break
                    elif fw_step_kmer == curr_kmer:
                        break
                    elif fw_step_kmer == kmer_twin(first_kmer):
                        break
                    bw_list.append(fw_step_kmer)
                    first_kmer = fw_step_kmer
            curr_kmer_list = []
            for bw in reversed(bw_list):
                curr_kmer_list.append(kmer_twin(bw))
            for fw in fw_list:
                curr_kmer_list.append(fw)
            contig_list_len = len(curr_kmer_list)
            curr_contig_id = len(contig_list)
            contig = Contig()
            contig.id = curr_contig_id
            contig.n_of_kmer = contig_list_len
            # 存储kmer的信息
            for j in range(contig_list_len):
                each_kmer = curr_kmer_list[j]
                tmp_rep = kmer_rep(each_kmer)
                if tmp_rep in kmer_str_dict:
                    kmer_str_dict[tmp_rep].contig_id = curr_contig_id
                    kmer_str_dict[tmp_rep].pos_in_contig = j
                    kmer_str_dict[tmp_rep].n_of_kmer_in_contig = contig_list_len
                    if tmp_rep == each_kmer:
                        kmer_str_dict[tmp_rep].sense_in_contig = 1
                    else:
                        kmer_str_dict[tmp_rep].sense_in_contig = 0
                else:
                    program_stop("contig.py")
                if j == 0:
                    contig.seq = each_kmer
                elif j > 0:
                    contig.seq = contig.seq + each_kmer[-1]
            contig_list.append(contig)
            ecs_list.append(-1)


def fwstep(kmer):
    global fw_step_kmer
    fw_count = 0
    for fw_kmer in contig_fw(kmer):
        curr_rep = kmer_rep(fw_kmer)
        if curr_rep in kmer_str_dict:
            fw_count = fw_count + 1
            curr_fw = fw_kmer
        if fw_count > 1:
            return 0
    if fw_count != 1:
        return 0
    bw_count = sum(kmer_rep(tmp) in kmer_str_dict for tmp in contig_bw(curr_fw))
    if bw_count > 1:
        return 0
    if bw_count != 1:
        return 0
    else:
        if curr_fw != kmer:
            fw_step_kmer = curr_fw
            return 1
        else:
            return 0

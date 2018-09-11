import functools

from _class import kmer_str_dict, kmer_rep, skip, ecs_list, ec_map, ec_inv_dict, transcript_info

new_ecs = []
num_read = 0


def compare_in_v(m, n):
    if m[0].contig_id == n[0].contig_id:
        return m[1] - n[1]
    else:
        return m[0].contig_id - n[0].contig_id


def compare_in_test(m, n):
    test_list = [([1, 3, 4], 1), ([2, 4, 5], 3), ([2, 4, 5], 2), ([3, 3, 4], 1)]
    if m[0][0] == n[0][0]:
        return m[1] - n[1]
    else:
        return m[0][0] - n[0][0]


def intersect_ec_u(ec, u):
    v = ec_map[ec]
    res = []
    i = 0
    j = 0
    while i < len(v) and j < len(u):
        if v[i] < u[j]:
            i = i + 1
        elif u[j] < v[i]:
            j = j + 1
        else:
            res.append(v[i])
            i = i + 1
            j = j + 1
    return res


# def intersect_ec_u(ec, u):
#     res = set()
#     if ec < len(ec_map):
#         new_u = ec_map[ec]
#         for each_tran in u:
#             res.add(each_tran)
#         for each_tran in new_u:
#             res.add(each_tran)
#     return sorted(list(res))


def intersect_ecs(v):
    if len(v) == 0:
        return []
    else:
        sorted_v = sorted(v, key=functools.cmp_to_key(compare_in_v))
    ec = ecs_list[sorted_v[0][0].contig_id]
    last_ec = ec
    u = ec_map[ec]
    i = 1
    while i < len(sorted_v):
        if sorted_v[i][0].contig_id != sorted_v[i - 1][0].contig_id:
            ec = ecs_list[sorted_v[i][0].contig_id]
            if ec != last_ec:
                u = intersect_ec_u(ec, u)
                # 借此找到,一个read到底map到了哪几个序列上去.这一步是取交集.
                # 如果一个read只map到了一个contig,好的,没关系
                # map到了contig1,contig2,且他们属于一个ecs,好的,还是没关系
                # 如果map到的多个contig不属于同一个ecs,是对多个ecs中记录的转录本id取交集.
                last_ec = ec
            if len(u) == 0:
                return u
        i = i + 1
    return u


def intersect_u1_u2(u1, u2):
    res = []
    i = 0
    j = 0
    while i < len(u1) and j < len(u2):
        if u1[i] < u2[j]:
            i = i + 1
        elif u2[j] < u1[i]:
            j = j + 1
        else:
            res.append(u1[i])
            i = i + 1
            j = j + 1
    return res


# def intersect_u1_u2(u1, u2):
#     res = set()
#     for each_tran in u1:
#         res.add(each_tran)
#     for each_tran in u2:
#         res.add(each_tran)
#     return sorted(list(res))


def intersect_ecs_paired(v1, v2):
    u1 = intersect_ecs(v1)
    u2 = intersect_ecs(v2)
    # u1和u2分别是l_read和r_read分别准确找到的map到了哪几条转录本上
    if len(u1) == 0 and len(u2) == 0:
        return []
    if len(u1) == 0:
        if len(v1) == 0:
            return u2
        else:
            return []
    elif len(u2) == 0:
        if len(v2) == 0:
            return u1
        else:
            return []
    else:
        # 取交集,
        # 如果这一对reads, map到了一个转录本上, 美滋滋, 没问题
        # 如果是map到了多个转录本上, 但l和r map上去的一毛一样, 没问题
        # 如果map到了多个转录本, 且l和r之间有差异, 取交集
        return intersect_u1_u2(u1, u2)


def find_ec(each_u):
    if len(each_u) == 0:
        return -1
    elif len(each_u) == 1:
        return each_u[0]
    else:
        if each_u in ec_inv_dict.values():
            for (key, value) in ec_inv_dict.items():
                if value == each_u:
                    return key
        else:
            return -1


def compute_mean_flg_lens(f_lens):
    f_lens_len = len(f_lens)
    tmp_counts = [0] * f_lens_len
    tmp_mass = [0.0] * f_lens_len
    mean_f_lens = [0.0] * f_lens_len
    tmp_counts[0] = f_lens[0]
    i = 1
    while i < f_lens_len:
        tmp_mass[i] = f_lens[i] * i + tmp_mass[i - 1]
        tmp_counts[i] = f_lens[i] + tmp_counts[i - 1]
        if tmp_counts[i] > 0:
            mean_f_lens[i] = tmp_mass[i] / tmp_counts[i]
        i = i + 1
    return mean_f_lens


def get_each_tran_len(mean_f_lens):
    tran_num = transcript_info.tran_num
    tran_len = transcript_info.tran_len
    frag_len_means = [0.0] * tran_num
    final_mean_f_lens = mean_f_lens[-1]
    i = 0
    while i < tran_num:
        if tran_len[i] >= len(mean_f_lens):
            frag_len_means[i] = final_mean_f_lens
        else:
            frag_len_means[i] = mean_f_lens[tran_len[i]]
        i = i + 1
    return frag_len_means


# def count(u):
#     for each_u in u:
#         ec = find_ec(each_u)
#         if ec == -1 or ec >= len(counts):
#             new_ecs.append(each_u)
#         else:
#             counts[ec] = counts[ec] + 1


def map_pair(seq1, seq1_len, seq2, seq2_len, k):
    d1 = 1
    d2 = 1
    p1 = -1
    p2 = -1
    c1 = -1
    c2 = -1

    found1 = 0
    for i in range(seq1_len - k + 1):
        kit1 = seq1[i:i + k]
        kit1_rep = kmer_rep(kit1)
        if kit1_rep in kmer_str_dict:
            found1 = 1
            kit1_info = kmer_str_dict[kit1_rep]
            c1 = kit1_info.contig_id
            if (kit1_rep == kit1) == kit1_info.sense_in_contig:
                p1 = kit1_info.pos_in_contig - i
                d1 = 1
            else:
                p1 = kit1_info.pos_in_contig + k + i
                d1 = 0
            # print(seq1)
            # print(kit1_info.pos_in_contig)
            break
        else:
            pass
    if not found1:
        return -1
    found2 = 0
    for i in range(seq2_len - k + 1):
        kit2 = seq2[i:i + k]
        kit2_rep = kmer_rep(kit2)
        if kit2_rep in kmer_str_dict:
            found2 = 1
            kit2_info = kmer_str_dict[kit2_rep]
            c2 = kit2_info.contig_id
            if (kit2_rep == kit2) == kit2_info.sense_in_contig:
                p2 = kit2_info.pos_in_contig - i
                d2 = 1
            else:
                p2 = kit2_info.pos_in_contig + k + i
                d2 = 0
            # print(kit2_info.pos_in_contig)
            break
        else:
            pass
    if not found2:
        return -1
    if c1 != c2:
        return -1
    if d1 and d2:
        return -1
    if not d1 and not d2:
        return -1
    if p1 > p2:
        return p1 - p2
    else:
        return p2 - p1


def get_dist(is_fw, kmer_info):
    if is_fw == kmer_info.sense_in_contig:
        return kmer_info.n_of_kmer_in_contig - kmer_info.pos_in_contig - 1
    else:
        return kmer_info.pos_in_contig


def match(curr_read, curr_read_len, k):
    # global num_read
    # num_read = num_read + 1
    # if num_read %2 == 0:
    #     print(num_read)
    v = []
    i = 0
    again = 0
    while i < curr_read_len - k + 1:
        kit = curr_read[i:i + k]
        kit_rep = kmer_rep(kit)
        if kit_rep in kmer_str_dict:
            curr_kit_info = kmer_str_dict[kit_rep]
            pos = i
            v.append((curr_kit_info, pos))
            dist = get_dist(int((kit == kit_rep)), curr_kit_info)
            if dist > 2:
                pos2 = pos + dist
                if pos2 + k >= curr_read_len:
                    # 超界,检查该read的最后一个kmer
                    pos2 = curr_read_len - k
                kit2 = curr_read[pos2:pos2 + k]
                kit2_rep = kmer_rep(kit2)
                if kit2_rep in kmer_str_dict:
                    found2 = 0
                    found2pos = pos + dist
                    kit2_info = kmer_str_dict[kit2_rep]
                    if curr_kit_info.contig_id == kit2_info.contig_id:
                        found2 = 1
                        found2pos = pos + dist
                else:
                    found2 = 1
                    found2pos = pos
                if found2:
                    if found2pos >= curr_read_len - k:
                        v.append((curr_kit_info, curr_read_len - k))
                        return v
                    else:
                        v.append((curr_kit_info, found2pos))
                else:
                    found_middle = 0
                    if dist > 4:
                        middle_pos = int((pos + pos2) / 2)
                        found3pos = pos + dist
                        kit3 = curr_read[middle_pos:middle_pos + k]
                        kit3_rep = kmer_rep(kit3)
                        if kit3_rep in kmer_str_dict:
                            kit3_info = kmer_str_dict[kit3_rep]
                            middle_contig_id = kit3_info.contig_id
                            if middle_contig_id == curr_kit_info.contig_id:
                                found_middle = 1
                                found3pos = middle_pos
                            elif middle_contig_id == kit2_info.contig_id:
                                found_middle = 1
                                found3pos = pos + dist
                            else:
                                pass
                            if found_middle:
                                v.append((kit3_info, found3pos))
                                if pos2 >= curr_read_len - k:
                                    return v
                                else:
                                    i = pos2
                            else:
                                pass
                        else:
                            pass
                    else:
                        pass
                    if not found_middle:
                        i = i + 1
                        again = 1
                        if again:
                            j = 0
                            m = 0
                            while j < curr_read_len and m < curr_read_len:
                                if j == skip:
                                    j = 0
                                if j == 0:
                                    again_kit = curr_read[m:m + k]
                                    again_kit_rep = kmer_rep(again_kit)
                                    if again_kit_rep in kmer_str_dict:
                                        again_kit_info = kmer_str_dict[again_kit_rep]
                                        v.append((again_kit_info, m))
                                if m == pos2:
                                    again = 0
                                    return v
                                m = m + 1
                                j = j + 1
            else:
                pass
        else:
            pass
        i = i + 1
    return v

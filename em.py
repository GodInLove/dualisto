import math

from _class import transcript_info, ec_map, program_stop


def compute_eff_len(tran_lens_estimated):
    tran_len = transcript_info.tran_len
    tran_lens_estimated_len = len(tran_lens_estimated)
    eff_tran_lens = [0.0] * tran_lens_estimated_len
    i = 0
    while i < tran_lens_estimated_len:
        curr_tran_len = tran_len[i]
        curr_eff_len = curr_tran_len - tran_lens_estimated[i] + 1.0
        if curr_eff_len < 1.0:
            curr_eff_len = curr_tran_len
        eff_tran_lens[i] = curr_eff_len
        i = i + 1
    return eff_tran_lens


def compute_weight_map(a_counts, eff_tran_lens):
    ec_map_len = len(ec_map)
    weight_map = [[]] * len(ec_map)
    for ec in range(ec_map_len):
        v_list = ec_map[ec]
        for each_tran_id in v_list:
            weight_map[ec].append(a_counts[ec] / eff_tran_lens[each_tran_id])
    return weight_map


def init_alpha_list():
    alpha_list = []
    for i in range(transcript_info.tran_num):
        alpha_list.append(1.0 / transcript_info.tran_num)
    return alpha_list


def em_run(a_counts, tran_lens_estimated, max_round=10000, min_round=50, bias=0):
    eff_tran_lens = compute_eff_len(tran_lens_estimated)
    weight_map = compute_weight_map(a_counts, eff_tran_lens)
    alpha_list = init_alpha_list()
    next_alpha_list = [0.0] * len(alpha_list)
    alpha_limit = 1e-7
    alpha_change_limit = 1e-2
    alpha_change = 1e-2
    tolerance = 4.94066e-324
    final_round = 0
    i = 0
    while i < max_round and not final_round:
        for ec in range(transcript_info.tran_num):
            next_alpha_list[ec] = a_counts[ec]
        ec = transcript_info.tran_num
        while ec < len(ec_map):
            tmp_v = 0.0
            if a_counts[ec] > 0:
                wv = weight_map[ec]
                v_list = ec_map[ec]
                v_list_len = len(v_list)
                for t in range(v_list_len):
                    tmp_v = tmp_v + alpha_list[v_list[t]] * wv[t]
                if tmp_v > tolerance:
                    count_normal = a_counts[ec] / tmp_v
                    for t in range(v_list_len):
                        next_alpha_list[v_list[t]] = next_alpha_list[v_list[t]] + \
                                                     (wv[t] * alpha_list[v_list[t]]) * count_normal
                else:
                    ec = ec + 1
            else:
                ec = ec + 1
            ec = ec + 1
        stop_em = 0
        stop_num = 0
        for ec in range(transcript_info.tran_num):
            if next_alpha_list[ec] > alpha_change_limit and \
                    (abs(next_alpha_list[ec] - alpha_list[ec]) / next_alpha_list[ec]) > alpha_change:
                stop_num = stop_num + 1
            alpha_list[ec] = next_alpha_list[ec]
            next_alpha_list[ec] = 0.0
        if stop_num == 0 and i > min_round:
            stop_em = 1
        if stop_em:
            final_round = 1
            alpha_before_zero = [0.0] * len(alpha_list)
            for ec in range(transcript_info.tran_num):
                alpha_before_zero[ec] = alpha_list[ec]
                if alpha_list[ec] < alpha_limit / 10.0:
                    alpha_list[ec] = 0.0
        i = i + 1
    if i == max_round:
        alpha_before_zero = [0.0] * len(alpha_list)
        for ec in range(transcript_info.tran_num):
            alpha_before_zero[ec] = alpha_list[ec]
    return alpha_list, eff_tran_lens


def write_em_tsv(file_name, alpha_list, eff_lens):
    tran_name = transcript_info.tran_name
    tran_len = transcript_info.tran_len
    try:
        fp = open(file_name, 'w')
    except IOError:
        print("cannot open", file_name)
        program_stop("em.py")
    fp.write("tran_name\ttran_len\teff_len\tcounts" + "\n")
    for ec in range(len(alpha_list)):
        fp.write(tran_name[ec] + "\t" + str(tran_len[ec]) +
                 "\t" + str(round(eff_lens[ec], 4)) + "\t" + str(alpha_list[ec]) + "\n")
    fp.close()

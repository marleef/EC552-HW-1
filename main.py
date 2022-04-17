"""
EC552 HW 1: Genetic Circuit Design Program
Drew Gross and Marlee Feltham
"""

import sys
import json
import math


# circuit: NOT and NOR gate
# pLuxStar----- NOT P3_PhlF ________
#                      pTet ________|---NOR A1_AmtR--output
# not input: pluxstar
# not gate: P3_PhlF
# nor input a: not output
# nor input b: ptet
# nor gate: A1_AmtR

# ======================================================================================
# ====================================== FILE R/W ======================================
# ======================================================================================


def read_file(fname):
    # open and read .input and .UCF JSON files
    with open(fname, 'r') as file:
        content = file.read()
    data = json.loads(content)
    return data


def write_output(fname, data):
    # open and write to .output JSON file
    with open(fname, 'w') as file:
        json.dump(data, file)
    file.close()

# ======================================================================================
# ======================================= PARSER =======================================
# ======================================================================================


def parse_UCF(data, gate_not, gate_nor):
    ucf = {}
    name = []
    ymax = []
    ymin = []
    n = []
    k = []
    for i in range(len(data)):
        if data[i]["collection"] == "models":
            if data[i]['name'] == gate_not or data[i]['name'] == gate_nor:
                name.append(data[i]['name'])
                for j in range(len(data[i]['parameters'])):
                    if data[i]['parameters'][j]['name'] == 'ymax':
                        ymax.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'ymin':
                        ymin.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'n':
                        n.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'K':
                        k.append(data[i]['parameters'][j]['value'])

    ucf['name'] = name
    ucf['ymin'] = ymin
    ucf['ymax'] = ymax
    ucf['n'] = n
    ucf['K'] = k

    print(ucf)
    return ucf


def parse_input(data, not_prom, nor_prom):
    inputs = {}
    name = []
    ymax = []
    ymin = []
    n = []
    k = []
    for i in range(len(data)):
        if data[i]['collection'] == 'models':
            if data[i]['name'] == not_prom or data[i]['name'] == nor_prom:
                name.append(data[i]['name'])
                for j in range(len(data[i]['parameters'])):
                    if data[i]['parameters'][j]['name'] == 'ymax':
                        ymax.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'ymin':
                        ymin.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'alpha':
                        n.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'beta':
                        k.append(data[i]['parameters'][j]['value'])

    inputs['name'] = name
    inputs['ymin'] = ymin
    inputs['ymax'] = ymax
    inputs['n'] = n
    inputs['K'] = k

    print(inputs)
    return inputs

# ======================================================================================
# ===================================== OPERATIONS =====================================
# ======================================================================================


def find_idx(inputs, gate):
    for i in range(len(inputs['name'])):
        if inputs['name'][i] == gate:
            return i


def stretch(inputs, gate, x):
    new_inputs = inputs.copy()
    i = find_idx(inputs, gate)

    new_inputs['ymax'][i] = inputs['ymax'][i]*x
    new_inputs['ymin'][i] = inputs['ymin'][i]/x

    return new_inputs


def promoter(inputs, pick, gate, x):
    new_inputs = inputs.copy()
    i = find_idx(inputs, gate)

    if pick == 0:
        # weaker promoter
        new_inputs['ymax'][i] = inputs['ymax'][i]/x
        new_inputs['ymin'][i] = inputs['ymin'][i]/x

    elif pick == 1:
        # stronger promoter
        new_inputs['ymax'][i] = inputs['ymax'][i]*x
        new_inputs['ymin'][i] = inputs['ymin'][i]*x

    return inputs


def slope(inputs, pick, gate, x):
    new_inputs = inputs.copy()
    i = find_idx(inputs, gate)

    if pick == 0:
        # decrease slope
        new_inputs['n'][i] = inputs['n'][i]/x
    elif pick == 1:
        # increase slope
        new_inputs['n'][i] = inputs['n'][i]*x

    return new_inputs


def rbs(inputs, pick, gate, x):
    new_inputs = inputs.copy()
    i = find_idx(inputs, gate)

    if pick == 0:
        # weaker rbs
        new_inputs['k'][i] = inputs['k']*x
    elif pick == 1:
        # stronger rbs
        new_inputs['k'][i] = inputs['k']/x

    return new_inputs


# ======================================================================================
# =================================== SCORE CIRCUIT ====================================
# ======================================================================================
def nor_gate(ucf, inputs, not_output):
    x = [inputs['ymin'][0]+not_output[0],
         inputs['ymin'][0]+not_output[1],
         inputs['ymax'][0]+not_output[0],
         inputs['ymax'][0]+not_output[1]]

    ttable = [0]*4
    ttable[0] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[0]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[1] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[1]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[2] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[2]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[3] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[3]/ucf['K'][0]) ** ucf['n'][0]))

    on_min = ttable[1]
    off_max = ttable[3]
    score = math.log10(on_min/off_max)

    return ttable, score


def not_gate(ucf, inputs):
    ttable = [0]*2
    ttable[0] = float(ucf['ymin'][1] + (ucf['ymax'][1]-ucf['ymin'][1]) /
                      (1+(inputs['ymin'][1]/ucf['K'][1]) ** ucf['n'][1]))
    ttable[1] = float(ucf['ymin'][1] + (ucf['ymax'][1]-ucf['ymin'][1]) /
                      (1+(inputs['ymax'][1]/ucf['K'][1]) ** ucf['n'][1]))
    # print(ttable)
    return ttable

# score = math.log10(on_min/off_max)

# ======================================================================================
# ===================================== NEW VALUES =====================================
# ======================================================================================


# def compare(scores):
#     # returns the index of scores[] that returns the minimum value > 0
#     position = scores.index(min(i for i in scores if i > 0))

#     return position


# def y_decision(size, UCF, inputs):
#     # performs a combination of stretch and promoter operations
#     # returns the best (lowest) score and combination of operations that modify
#     # ymax and ymin

#     # (1) score the circuit without applying any operations
#     # (2) score the circuit after applying the stretch operation
#     # (3) compare the scores from (1) and (2)
#     # -> get best score and ymax_new and ymin_new value
#     # (4) use ymax_new and ymin_new values from (3) as inputs for promoter()
#     # (5) score the circuits
#     # (6) compare the scores and return the best score and ymax_new and ymin_new values

#     ymax_r0 = []
#     ymin_r0 = []
#     ymax_r0[0] = UCF['ymax']
#     ymin_r0[0] = UCF['ymin']

#     # score the circuit with the original ymax and ymin values
#     original_score = score_circuit(size, UCF, inputs)

#     # apply stretch operation and score the circuit
#     # # compare the original and stretch scores
#     # r1 indicates which score is the lowest

#     # change into dict
#     # ymax_r0[1], ymin_r0[1] = stretch(inputs, UCF)
#     stretch_score = score_circuit(size, UCF, inputs)

#     scores = [original_score, stretch_score]
#     r1 = compare(scores)

#     # ymax_r1 and ymin_r1 are the specified ymax and ymin values resulting from r0
#     ymax_r1 = ymax_r0[r1]
#     ymin_r1 = ymin_r0[r1]

#     ymax_r2 = []
#     ymin_r2 = []
#     scores = []

#     for i in range(1):
#         ymax_r2[i], ymin_r2[i] = promoter(x, ymax_r1, ymin_r1, i)
#         scores[i] = score_circuit(size, ymin_r2[i], ymax_r2[i], n, k, x)

#     # r2 = compare(scores[0], scores[1])
#     r2_pos = compare(scores)
#     ymax_new = ymax_r2[r2_pos]
#     ymin_new = ymin_r2[r2_pos]

#     return ymax_new, ymin_new, scores[r2_pos]


# def n_decision(size, ymin, ymax, n, k, x):
#     # performs slope operation
#     # returns the best (lowest) score and combination of operations that modify n
#     original_score = score_circuit(size, ymin, ymax, n, k, x)

#     n_vals = []
#     scores = []
#     n_vals[0] = n
#     scores[0] = original_score

#     for i in range(1):
#         n_vals[i+1] = slope(n, x, i)
#         scores[i+1] = score_circuit(size, ymin, ymax, n_vals[i+1], k, x)

#     pos = compare(scores)

#     return n_vals[pos], scores[pos]


# def k_decision(size, ymin, ymax, n, k, x):
#     # performs RBS operation
#     # returns the best (lowest) score and combination of operations that modify K
#     original_score = score_circuit(size, ymin, ymax, n, k, x)

#     k_vals = []
#     scores = []
#     k_vals[0] = k
#     scores[0] = original_score

#     for i in range(1):
#         k_vals[i+1] = rbs(k, x, i)
#         scores[i+1] = score_circuit(size, ymin, ymax, n, k_vals[i+1], x)

#     pos = compare(scores)

#     return k_vals[pos], scores[pos]


# def best_score(size, ymin, ymax, n, k, x):
#     # (1) retrieve outputs of y, n, and k_decision
#     # (2) add scores to list scores[]
#     # (3) generate scores from all combinations of each new parameter
#     # (4) find and return the best score
#     ymax_new, ymin_new, y_decision_score = y_decision(
#         size, ymin, ymax, n, k, x)
#     n_new, n_decision_score = n_decision(size, ymin, ymax, n, k, x)
#     k_new, k_decision_score = k_decision(size, ymin, ymax, n, k, x)

#     scores = [y_decision_score, n_decision_score, k_decision_score]

#     ymaxs = [ymax, ymax_new]
#     ymins = [ymin, ymin_new]
#     ns = [n, n_new]
#     ks = [k, k_new]

#     params = 3
#     num = 2**params  # 8 combinations

#     for i in range(num-1):
#         bin = format(i, "b")
#         scores[i+3] = score_circuit(size, ymins[bin[0]],
#                                     ymaxs[bin[0]], ns[bin[1]], ks[bin[2]])

#    # best_score = min(i for i in scores if i > 0)
#     pos = compare(scores)
#     best_score = scores[pos]

#     return best_score
def x_in():
    x = input("Define x value (0 < x <= 1.05)\n")
    if float(x) <= 0 or float(x) > 1.05:
        sys.exit("Invalid x value.\n")
    return x
# ======================================================================================
# ========================================MAIN==========================================
# ======================================================================================


def main():
    chassis_name = 'Eco1C1G1T1'
    in_ucf = f'{chassis_name}.UCF.json'
    input_sensor_file = f'{chassis_name}.input.json'
    output_device_file = f'{chassis_name}.output.json'
    in_param = read_file(input_sensor_file)
    UCF_param = read_file(in_ucf)
    # inputs = parse_input(in_param)
    gate_not = 'A1_AmtR_model'
    gate_nor = 'P3_PhlF_model'
    nor_prom = 'TetR_sensor_model'
    not_prom = 'LuxR_sensor_model'
    ucf = parse_UCF(UCF_param, gate_not, gate_nor)
    inputs = parse_input(in_param, not_prom, nor_prom)
    # not_output = not_gate(ucf, inputs)
    # nor_output = nor_gate(ucf, inputs, not_output)

    print("\n======== INPUT SIGNALS ======== \n")
    # ans = input(
    #     "Choose a NOT gate promoter:\n(a) LacI\n(b) TetR\n(c) AraC\n(d) LuxR\n(x) default\n")

    # if ans == 'a':
    #     not_prom = 'LacI_sensor_model'
    # elif ans == 'b':
    #     not_prom = 'TetR_sensor_model'
    # elif ans == 'c':
    #     not_prom = 'AraC_sensor_model'
    # elif ans == 'd':
    #     not_prom = 'LuxR_sensor_model'
    # elif ans == 'x':
    #     not_prom = 'LuxR_sensor_model'
    # else:
    #     sys.exit('Invalid  entry.\n')

    operation = input("Choose up to 4 operations from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
    operation = [i for i in operation]

    if len(operation) > 4:
        sys.exit("Invalid entry. Too many operations.\n")

    for i in range(len(operation)):
        if operation[i] == 'a':
            new_inputs = stretch(inputs, not_prom, float(x_in()))
        elif operation[i] == 'b':
            new_inputs = slope(inputs, 1, not_prom, float(x_in()))
        elif operation[i] == 'c':
            new_inputs = slope(inputs, 0, not_prom, float(x_in()))
        elif operation[i] == 'd':
            new_inputs = promoter(inputs, 1, not_prom, float(x_in()))
        elif operation[i] == 'e':
            new_inputs = promoter(inputs, 0, not_prom, float(x_in()))
        elif operation[i] == 'f':
            new_inputs = rbs(inputs, 1, not_prom, float(x_in()))
        elif operation[i] == 'g':
            new_inputs = rbs(inputs, 0, not_prom, float(x_in()))
        elif operation[i] == 'x':
            new_inputs = inputs.copy()
            break
    # inputs['ymin'][0] = inputs['ymin'][0]*2

    # print(inputs)
    # print(new_inputs)
    not_output = not_gate(ucf, new_inputs)
    nor_output, score = nor_gate(ucf, new_inputs, not_output)
    print(score)
    print('\n')

    # ans = input(
    #     "Choose a NOR gate promoter:\n(a) LacI\n(b) TetR\n(c) AraC\n(d) LuxR\n(x) default\n")

    # if ans == 'a':
    #     nor_prom = 'LacI_sensor_model'
    # elif ans == 'b':
    #     nor_prom = 'TetR_sensor_model'
    # elif ans == 'c':
    #     nor_prom = 'AraC_sensor_model'
    # elif ans == 'd':
    #     nor_prom = 'LuxR_sensor_model'
    # elif ans == 'x':
    #     nor_prom = 'TetR_sensor_model'
    # else:
    #     sys.exit('Invalid  entry\n')

    # print("\n============ GATES ============\n")
    # # NOR gate command line
    # ans = input(
    #     "Choose a NOR gate promoter:\n(a) A1_AmtR\n(b) B#_BM3R1\n(c) E1_BetI\n(d) F1_AmeR\n(e) H1_HlyIIR\n(f) I1_IcaRA\n(g) L1_LitR\n(h) N1_LmrA\n(i) P#_PhlF\n(j) Q2_QacR\n(k) R1_PsrA\n(l) S1_SrpR_model\n(x) default\n")
    # if ans == 'a':
    #     gate_nor = 'A1_AmtR_model'
    # elif ans == 'b':
    #     gate_num = input(
    #         '(1), (2), or (3)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_nor = 'B1_BM3R1_model'
    #     elif gate_num == '2':
    #         gate_nor = 'B2_BM3R1_model'
    #     elif gate_num == '3':
    #         gate_nor = 'B3_BM3R1_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'c':
    #     gate_nor = 'E1_BetI_model'
    # elif ans == 'd':
    #     gate_nor = 'F1_AmeR_model'
    # elif ans == 'e':
    #     gate_nor = 'H1_HlyIIR_model'
    # elif ans == 'f':
    #     gate_nor = 'I1_IcaRA_model'
    # elif ans == 'g':
    #     gate_nor = 'L1_LitR_model'
    # elif ans == 'h':
    #     gate_nor = 'N1_LmrA_model'
    # elif ans == 'i':
    #     gate_num = input('(1), (2), or (3)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_nor = 'P1_PhlF_model'
    #     elif gate_num == '2':
    #         gate_nor = 'P2_PhlF_model'
    #     elif gate_num == '3':
    #         gate_nor = 'P3_PhlF_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'j':
    #     gate_nor = 'Q2_QacR_model'
    # elif ans == 'k':
    #     gate_nor = 'R1_PsrA_model'
    # elif ans == 'l':
    #     gate_num = input('(1), (2), (3), or (4)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_nor = 'S1_SrpR_model'
    #     elif gate_num == '2':
    #         gate_nor = 'S2_SrpR_model'
    #     elif gate_num == '3':
    #         gate_nor = 'S3_SrpR_model'
    #     elif gate_num == '4':
    #         gate_nor = 'S4_SrpR_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'x':
    #     gate_nor = 'P3_PhlF_model'

    # # NOT gate command line
    # ans = input(
    #     "Choose a NOT gate promoter:\n(a) A1_AmtR\n(b) B#_BM3R1\n(c) E1_BetI\n(d) F1_AmeR\n(e) H1_HlyIIR\n(f) I1_IcaRA\n(g) L1_LitR\n(h) N1_LmrA\n(i) P#_PhlF\n(j) Q2_QacR\n(k) R1_PsrA\n(l) S1_SrpR_model\n(x) default\n")
    # if ans == 'a':
    #     gate_not = 'A1_AmtR_model'
    # elif ans == 'b':
    #     gate_num = input(
    #         '(1), (2), or (3)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_not = 'B1_BM3R1_model'
    #     elif gate_num == '2':
    #         gate_not = 'B2_BM3R1_model'
    #     elif gate_num == '3':
    #         gate_not = 'B3_BM3R1_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'c':
    #     gate_not = 'E1_BetI_model'
    # elif ans == 'd':
    #     gate_not = 'F1_AmeR_model'
    # elif ans == 'e':
    #     gate_not = 'H1_HlyIIR_model'
    # elif ans == 'f':
    #     gate_not = 'I1_IcaRA_model'
    # elif ans == 'g':
    #     gate_not = 'L1_LitR_model'
    # elif ans == 'h':
    #     gate_not = 'N1_LmrA_model'
    # elif ans == 'i':
    #     gate_num = input('(1), (2), or (3)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_not = 'P1_PhlF_model'
    #     elif gate_num == '2':
    #         gate_not = 'P2_PhlF_model'
    #     elif gate_num == '3':
    #         gate_not = 'P3_PhlF_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'j':
    #     gate_not = 'Q2_QacR_model'
    # elif ans == 'k':
    #     gate_not = 'R1_PsrA_model'
    # elif ans == 'l':
    #     gate_num = input('(1), (2), (3), or (4)? Enter the corresponding #\n')
    #     if gate_num == '1':
    #         gate_not = 'S1_SrpR_model'
    #     elif gate_num == '2':
    #         gate_not = 'S2_SrpR_model'
    #     elif gate_num == '3':
    #         gate_not = 'S3_SrpR_model'
    #     elif gate_num == '4':
    #         gate_not = 'S4_SrpR_model'
    #     else:
    #         sys.exit("Invalid entry\n")
    # elif ans == 'x':
    #     gate_not = 'A1_AmtR_model'


if __name__ == "__main__":
    main()

"""
EC552 HW 1: Genetic Circuit Design Program
Drew Gross and Marlee Feltham
"""

import sys
import json
import math
import copy


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

    # print('ucf parameters: ', ucf)
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

    print('input parameters: ', inputs)
    return inputs

# ======================================================================================
# ===================================== OPERATIONS =====================================
# ======================================================================================


def find_idx(ucf, gate):
    for i in range(len(ucf['name'])):
        if ucf['name'][i] == gate:
            return i


def stretch(ucf, gate, x):
    new_ucf = copy.deepcopy(ucf)
    i = find_idx(ucf, gate)

    new_ucf['ymax'][i] = ucf['ymax'][i]*x
    new_ucf['ymin'][i] = ucf['ymin'][i]/x

    # print('stretch: ', new_ucf)
    return new_ucf


def promoter(ucf, pick, gate, x):
    new_ucf = copy.deepcopy(ucf)
    i = find_idx(ucf, gate)

    if pick == 0:
        # weaker promoter
        new_ucf['ymax'][i] = ucf['ymax'][i]/x
        new_ucf['ymin'][i] = ucf['ymin'][i]/x

    elif pick == 1:
        # stronger promoter
        new_ucf['ymax'][i] = ucf['ymax'][i]*x
        new_ucf['ymin'][i] = ucf['ymin'][i]*x

    # print('promoter: ', new_ucf)

    return new_ucf


def slope(ucf, pick, gate, x):
    new_ucf = copy.deepcopy(ucf)
    i = find_idx(ucf, gate)

    if pick == 0:
        # decrease slope
        new_ucf['n'][i] = ucf['n'][i]/x
    elif pick == 1:
        # increase slope
        new_ucf['n'][i] = ucf['n'][i]*x

    # print('slope: ', new_ucf)

    return new_ucf


def rbs(ucf, pick, gate, x):
    new_ucf = copy.deepcopy(ucf)
    i = find_idx(ucf, gate)

    if pick == 0:
        # weaker rbs
        new_ucf['k'][i] = ucf['k']*x
    elif pick == 1:
        # stronger rbs
        new_ucf['k'][i] = ucf['k']/x

    # print('rbs: ', new_ucf)

    return new_ucf


# ======================================================================================
# =================================== SCORE CIRCUIT ====================================
# ======================================================================================
def nor_gate(ucf, inputs, not_output):
    x = [ucf['ymin'][0]+not_output[0],
         ucf['ymin'][0]+not_output[1],
         ucf['ymax'][0]+not_output[0],
         ucf['ymax'][0]+not_output[1]]

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


def score_circuit(ucf, inputs):
    not_output = not_gate(ucf, inputs)
    [ttable, score] = nor_gate(ucf, inputs, not_output)
    return score


# ======================================================================================
# ===================================== NEW VALUES =====================================
# ======================================================================================
def merge(ucf, param):
    cpy = copy.deepcopy(ucf[0])
    cpy2 = copy.deepcopy(ucf[1])

    for i in range(len(param)):
        cpy[param[i]][1] = cpy2[param[i]][1]
        # print('cpy param:', cpy[param[i]])
    return cpy


def y_decision(ucf, inputs, gate_nor, gate_not):
    # performs a combination of stretch and promoter operations
    # returns the best (lowest) score and combination of operations that modify
    # ymax and ymin

    gates = [gate_nor, gate_not]
    original_score = score_circuit(ucf, inputs)

    x = .85
    scores = [original_score]
    ucfs = []

    for i in range(len(gates)):
        ucfs.append(stretch(ucf, gates[i], x))
        scores.append(score_circuit(ucfs[i], inputs))

    ucfs.append(merge(ucfs, ['ymax', 'ymin']))
    scores.append(score_circuit(ucfs[i+1], inputs))
    ucfs.append(ucf)
    idx = scores.index(min(scores))  # best ucfs index
    # print('ucfs: ', ucfs)
    prom_ucf = []
    scores = [original_score]

    for i in range(len(gates)):
        prom_ucf.append(promoter(ucfs[idx], i, gates[i], x))
        scores.append(score_circuit(prom_ucf[i], inputs))

    prom_ucf.append(merge(prom_ucf, ['ymax', 'ymin']))
    scores.append(score_circuit(prom_ucf[i+1], inputs))

    ucfs = [ucf]
    for i in range(len(prom_ucf)):
        ucfs.append(prom_ucf[i])

    idx = scores.index(min(scores))  # best inputs index

    return ucfs[idx], scores[idx]


def n_decision(ucf, inputs, gate_nor, gate_not):
    # performs slope operation
    # returns the best (lowest) score and combination of operations that modify n
    gates = [gate_nor, gate_not]
    original_score = score_circuit(ucf, inputs)
    x = .85
    scores = [original_score]
    slope_ucf = []

    for i in range(len(gates)):
        slope_ucf.append(slope(ucf, i, gates[i], x))
        scores.append(score_circuit(slope_ucf[i], inputs))

    slope_ucf.append(merge(slope_ucf, 'n'))
    scores.append(score_circuit(slope_ucf[i+1], inputs))

    ins = [ucf]
    for i in range(len(slope_ucf)):
        ins.append(slope_ucf[i])
    # print('n ins: ', ins)

    idx = scores.index(min(scores))

    return ins[idx], scores[idx]


def k_decision(ucf, inputs, gate_nor, gate_not):
    # performs RBS operation
    # returns the best (lowest) score and combination of operations that modify K
    gates = [gate_nor, gate_not]
    original_score = score_circuit(ucf, inputs)
    x = .85
    scores = [original_score]
    rbs_ucf = []

    for i in range(len(gates)):
        rbs_ucf.append(slope(ucf, i, gates[i], x))
        scores.append(score_circuit(rbs_ucf[i], inputs))

    rbs_ucf.append(merge(rbs_ucf, 'K'))
    scores.append(score_circuit(rbs_ucf[i+1], inputs))

    ins = [ucf]
    for i in range(len(rbs_ucf)):
        ins.append(rbs_ucf[i])
    # print('n ins: ', ins)

    idx = scores.index(min(scores))

    return ins[idx], scores[idx]


def best_score(ucf, inputs, gate_nor, gate_not):
    y_ucf, y_score = y_decision(ucf, inputs, gate_nor, gate_not)
    n_ucf, n_score = n_decision(ucf, inputs, gate_nor, gate_not)
    k_ucf, k_score = k_decision(ucf, inputs, gate_nor, gate_not)

    scores = [y_score, n_score, k_score]
    # print('scores: ', scores)
    ucfs = [ucf, y_ucf, n_ucf, k_ucf]

    cpy = copy.deepcopy(y_ucf)
    cpy['n'] = n_ucf['n']
    ucfs.append(cpy)

    cpy = copy.deepcopy(y_ucf)
    cpy['K'] = k_ucf['K']
    ucfs.append(cpy)

    cpy['n'] = n_ucf['n']
    ucfs.append(cpy)

    cpy = copy.deepcopy(k_ucf)
    cpy['n'] = n_ucf['n']
    ucfs.append(cpy)
    # [inputs, y, n, k, y+n, y+k, y+k+n, n+k]

    for i in range(len(ucfs)-4):
        scores.append(score_circuit(ucfs[i+4], inputs))

    idx = scores.index(min(scores))

    return scores[idx]


def x_in():
    x = input("Define x value (0 < x <= 1.05): \n")
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

    gate_not = 'A1_AmtR_model'
    gate_nor = 'P3_PhlF_model'
    nor_prom = 'TetR_sensor_model'
    not_prom = 'LuxR_sensor_model'
    ucf = parse_UCF(UCF_param, gate_not, gate_nor)
    inputs = parse_input(in_param, not_prom, nor_prom)

    print('====== INPUT SIGNALS ================================ \n')
    operation = input("Choose up to 4 operations to perform on the NOT gate input from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
    operation = [i for i in operation]

    if len(operation) > 4:
        sys.exit("Invalid entry. Too many operations.\n")

    is_op = 1
    for i in range(len(operation)):
        if operation[i] == 'a':
            new_ucf = stretch(ucf, gate_nor, float(x_in()))
        elif operation[i] == 'b':
            new_ucf = slope(ucf, 1, gate_nor, float(x_in()))
        elif operation[i] == 'c':
            new_ucf = slope(ucf, 0, gate_nor, float(x_in()))
        elif operation[i] == 'd':
            new_ucf = promoter(ucf, 1, gate_nor, float(x_in()))
        elif operation[i] == 'e':
            new_ucf = promoter(ucf, 0, gate_nor, float(x_in()))
        elif operation[i] == 'f':
            new_ucf = rbs(ucf, 1, gate_nor, float(x_in()))
        elif operation[i] == 'g':
            new_ucf = rbs(ucf, 0, gate_nor, float(x_in()))
        elif operation[i] == 'x':
            is_op = 0
            break

    print('\n')
    print('====== Assignment ===================================')
    print('pTetR        0011')
    print('pLuxStar     0101')

    print('P3_PhlF      pLuxStar')
    print('A1_AmtR      pLuxStar    pTetR')

    print('Output       A1_AmtR')
    print('\n')
    if is_op == 0:
        print('====== Default Score (No Operations) ================')
        print(score_circuit(ucf, inputs))
        print('\n')
        print('====== Optimized Score ==============================')
        print(best_score(ucf, inputs, gate_nor, gate_not))
    elif is_op == 1:
        print('====== Score with Operations ========================')
        print(score_circuit(new_ucf, inputs))
    print('\n')


if __name__ == "__main__":
    main()

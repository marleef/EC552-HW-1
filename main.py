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
    new_inputs = copy.deepcopy(inputs)
    i = find_idx(inputs, gate)

    new_inputs['ymax'][i] = inputs['ymax'][i]*x
    new_inputs['ymin'][i] = inputs['ymin'][i]/x

    print('stretch: ', new_inputs)
    return new_inputs


def promoter(inputs, pick, gate, x):
    new_inputs = copy.deepcopy(inputs)
    i = find_idx(inputs, gate)

    if pick == 0:
        # weaker promoter
        new_inputs['ymax'][i] = inputs['ymax'][i]/x
        new_inputs['ymin'][i] = inputs['ymin'][i]/x

    elif pick == 1:
        # stronger promoter
        new_inputs['ymax'][i] = inputs['ymax'][i]*x
        new_inputs['ymin'][i] = inputs['ymin'][i]*x

    print('promoter: ', new_inputs)

    return new_inputs


def slope(inputs, pick, gate, x):
    new_inputs = copy.deepcopy(inputs)
    i = find_idx(inputs, gate)

    if pick == 0:
        # decrease slope
        new_inputs['n'][i] = inputs['n'][i]/x
    elif pick == 1:
        # increase slope
        new_inputs['n'][i] = inputs['n'][i]*x

    print('slope: ', new_inputs)

    return new_inputs


def rbs(inputs, pick, gate, x):
    new_inputs = copy.deepcopy(inputs)
    i = find_idx(inputs, gate)

    if pick == 0:
        # weaker rbs
        new_inputs['k'][i] = inputs['k']*x
    elif pick == 1:
        # stronger rbs
        new_inputs['k'][i] = inputs['k']/x

    print('rbs: ', new_inputs)

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
    print(ttable)
    return ttable


def score_circuit(ucf, inputs):
    not_output = not_gate(ucf, inputs)
    [ttable, score] = nor_gate(ucf, inputs, not_output)
    return score


# ======================================================================================
# ===================================== NEW VALUES =====================================
# ======================================================================================
def merge(inputs, param):
    cpy = copy.deepcopy(inputs[0])
    cpy2 = copy.deepcopy(inputs[1])

    for i in range(len(param)):
        cpy[param[i]][1] = cpy2[param[i]][1]
        print('cpy param:', cpy[param[i]])
    return cpy


def y_decision(ucf, inputs, nor_prom, not_prom):
    # performs a combination of stretch and promoter operations
    # returns the best (lowest) score and combination of operations that modify
    # ymax and ymin

    proms = [nor_prom, not_prom]
    original_score = score_circuit(ucf, inputs)

    x = .85
    scores = [original_score]
    ins = []

    for i in range(len(proms)):
        ins.append(stretch(inputs, proms[i], x))
        scores.append(score_circuit(ucf, ins[i]))

    ins.append(merge(ins, ['ymax', 'ymin']))
    scores.append(score_circuit(ucf, ins[i+1]))
    ins.append(inputs)
    idx = scores.index(min(scores))  # best inputs index
    print('ins: ', ins)
    prom_inputs = []
    scores = [original_score]

    for i in range(len(proms)):
        prom_inputs.append(promoter(ins[idx], i, proms[i], x))
        scores.append(score_circuit(ucf, prom_inputs[i]))

    prom_inputs.append(merge(prom_inputs, ['ymax', 'ymin']))
    scores.append(score_circuit(ucf, prom_inputs[i+1]))

    ins = [inputs]
    for i in range(len(prom_inputs)):
        ins.append(prom_inputs[i])

    idx = scores.index(min(scores))  # best inputs index

    return ins[idx], scores[idx]


def n_decision(ucf, inputs, nor_prom, not_prom):
    # performs slope operation
    # returns the best (lowest) score and combination of operations that modify n
    proms = [nor_prom, not_prom]
    original_score = score_circuit(ucf, inputs)
    x = .85
    scores = [original_score]
    slope_inputs = []

    for i in range(len(proms)):
        slope_inputs.append(slope(inputs, i, proms[i], x))
        scores.append(score_circuit(ucf, slope_inputs[i]))

    slope_inputs.append(merge(slope_inputs, 'n'))
    scores.append(score_circuit(ucf, slope_inputs[i+1]))

    # ins = [inputs, slope_inputs]
    ins = [inputs]
    for i in range(len(slope_inputs)):
        ins.append(slope_inputs[i])
    print('n ins: ', ins)

    idx = scores.index(min(scores))

    return ins[idx], scores[idx]


def k_decision(ucf, inputs, nor_prom, not_prom):
    # performs RBS operation
    # returns the best (lowest) score and combination of operations that modify K
    proms = [nor_prom, not_prom]
    original_score = score_circuit(ucf, inputs)
    x = .85
    scores = [original_score]
    rbs_inputs = []

    for i in range(len(proms)):
        rbs_inputs.append(slope(inputs, i, proms[i], x))
        scores.append(score_circuit(ucf, rbs_inputs[i]))

    rbs_inputs.append(merge(rbs_inputs, 'K'))
    scores.append(score_circuit(ucf, rbs_inputs[i+1]))

    # ins = [inputs, rbs_inputs]
    ins = [inputs]
    for i in range(len(rbs_inputs)):
        ins.append(rbs_inputs[i])
    idx = scores.index(min(scores))

    return ins[idx], scores[idx]


def best_score(ucf, inputs, nor_prom, not_prom):
    y_inputs, y_score = y_decision(ucf, inputs, nor_prom, not_prom)
    n_inputs, n_score = n_decision(ucf, inputs, nor_prom, not_prom)
    k_inputs, k_score = k_decision(ucf, inputs, nor_prom, not_prom)

    scores = [y_score, n_score, k_score]
    print('scores: ', scores)
    ins = [inputs, y_inputs, n_inputs, k_inputs]

    cpy = copy.deepcopy(y_inputs)
    cpy['n'] = n_inputs['n']
    ins.append(cpy)

    cpy = copy.deepcopy(y_inputs)
    cpy['K'] = k_inputs['K']
    ins.append(cpy)

    cpy['n'] = n_inputs['n']
    ins.append(cpy)

    cpy = copy.deepcopy(n_inputs)
    cpy['K'] = n_inputs['n']
    ins.append(cpy)
    # [inputs, y, n, k, y+n, y+k, y+k+n, n+k]

    for i in range(len(ins)-4):
        scores.append(score_circuit(ucf, ins[i+4]))

    idx = scores.index(min(scores))

    return scores[idx]


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

    gate_not = 'A1_AmtR_model'
    gate_nor = 'P3_PhlF_model'
    nor_prom = 'TetR_sensor_model'
    not_prom = 'LuxR_sensor_model'
    ucf = parse_UCF(UCF_param, gate_not, gate_nor)
    inputs = parse_input(in_param, not_prom, nor_prom)

    print("\n======== INPUT SIGNALS ======== \n")
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

    print(best_score(ucf, inputs, nor_prom, not_prom))
    print('\n')


if __name__ == "__main__":
    main()

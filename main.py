"""
EC552 HW 1: Genetic Circuit Design Program
Drew Gross and Marlee Feltham
"""

import os
import sys
import math
import json


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


def stretch(inputs, ucf, gate, x):
    # apply stretch operation
    for i in len(inputs['name']):
        if ucf['name'][i] == gate:
            break

    inputs['ymax'][i] = inputs['ymax'][i]*x
    inputs['ymin'][i] = inputs['ymin'][i]/x

    return inputs


def promoter(inputs, ucf, pick, signal, x):
    # apply weaker or stronger promoter operation based on pick value
    # pick == 0 -> weaker promoter
    # pick == 1 -> stronger promoter
    for i in len(inputs['name']):
        if inputs['name'][i] == signal:
            break

    if pick == 0:
        inputs['ymax'][i] = inputs['ymax']/x
        inputs['ymin'][i] = inputs['ymin']/x
    elif pick == 1:
        inputs['ymax'][i] = inputs['ymax'][i]*x
        inputs['ymin'][i] = inputs['ymin'][i]*x

    return inputs


def slope(inputs, pick, signal, x):
    # apply increase or decrease slope operation based on pick value
    # pick == 0 -> decrease slope
    # pick == 1 -> increase slope
    for i in len(inputs['name']):
        if inputs['name'][i] == signal:
            break

        if pick == 0:
            inputs['n'][i] = inputs['n']/x
        elif pick == 1:
            inputs['n'][i] = inputs['n']*x

    return inputs


def rbs(inputs, pick, signal, x):
    # apply weaker or stronger RBS operation based on pick value
    # pick == 0 -> weaker rbs
    # pick == 1 -> stronger rbs
    for i in len(inputs['name']):
        if inputs['name'][i] == signal:
            break

    if pick == 0:
        inputs['k'][i] = inputs['k']*x
    elif pick == 1:
        inputs['k'][i] = inputs['k']/x

    return inputs


# ======================================================================================
# =================================== SCORE CIRCUIT ====================================
# ======================================================================================
def nor_gate(ucf, inputs, not_output):
    x = [inputs['ymin'][0]+not_output[0],
         inputs['ymin'][0]+not_output[1],
         inputs['ymax'][0]+not_output[1],
         inputs['ymax'][0]+not_output[0]]
    print(x)

    ttable = [0]*4

    ttable[0] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[0]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[1] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[1]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[2] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[2]/ucf['K'][0]) ** ucf['n'][0]))
    ttable[3] = float(ucf['ymin'][0] + (ucf['ymax'][0]-ucf['ymin']
                      [0]) / (1+(x[3]/ucf['K'][0]) ** ucf['n'][0]))
    return ttable


def not_gate(ucf, inputs):
    ttable = [0]*2

    ttable[0] = float(ucf['ymin'][1] + (ucf['ymax'][1]-ucf['ymin'][1]) /
                      (1+(inputs['ymin'][0]/ucf['K'][1]) ** ucf['n'][1]))
    ttable[1] = float(ucf['ymin'][1] + (ucf['ymax'][1]-ucf['ymin'][1]) /
                      (1+(inputs['ymax'][0]/ucf['K'][1]) ** ucf['n'][1]))
    print(ttable)
    return ttable


# score = math.log10(on_min/off_max)


# ======================================================================================
# ===================================== NEW VALUES =====================================
# ======================================================================================

def compare(scores):
    # returns the index of scores[] that returns the minimum value > 0
    position = scores.index(min(i for i in scores if i > 0))

    return position


def y_decision(size, UCF, inputs):
    # performs a combination of stretch and promoter operations
    # returns the best (lowest) score and combination of operations that modify
    # ymax and ymin

    # (1) score the circuit without applying any operations
    # (2) score the circuit after applying the stretch operation
    # (3) compare the scores from (1) and (2)
    # -> get best score and ymax_new and ymin_new value
    # (4) use ymax_new and ymin_new values from (3) as inputs for promoter()
    # (5) score the circuits
    # (6) compare the scores and return the best score and ymax_new and ymin_new values

    ymax_r0 = []
    ymin_r0 = []
    ymax_r0[0] = UCF['ymax']
    ymin_r0[0] = UCF['ymin']

    # score the circuit with the original ymax and ymin values
    original_score = score_circuit(size, UCF, inputs)

    # apply stretch operation and score the circuit
    # # compare the original and stretch scores
    # r1 indicates which score is the lowest

    # change into dict
    # ymax_r0[1], ymin_r0[1] = stretch(inputs, UCF)
    stretch_score = score_circuit(size, UCF, inputs)

    scores = [original_score, stretch_score]
    r1 = compare(scores)

    # ymax_r1 and ymin_r1 are the specified ymax and ymin values resulting from r0
    ymax_r1 = ymax_r0[r1]
    ymin_r1 = ymin_r0[r1]

    ymax_r2 = []
    ymin_r2 = []
    scores = []

    for i in range(1):
        ymax_r2[i], ymin_r2[i] = promoter(x, ymax_r1, ymin_r1, i)
        scores[i] = score_circuit(size, ymin_r2[i], ymax_r2[i], n, k, x)

    # r2 = compare(scores[0], scores[1])
    r2_pos = compare(scores)
    ymax_new = ymax_r2[r2_pos]
    ymin_new = ymin_r2[r2_pos]

    return ymax_new, ymin_new, scores[r2_pos]


def n_decision(size, ymin, ymax, n, k, x):
    # performs slope operation
    # returns the best (lowest) score and combination of operations that modify n
    original_score = score_circuit(size, ymin, ymax, n, k, x)

    n_vals = []
    scores = []
    n_vals[0] = n
    scores[0] = original_score

    for i in range(1):
        n_vals[i+1] = slope(n, x, i)
        scores[i+1] = score_circuit(size, ymin, ymax, n_vals[i+1], k, x)

    pos = compare(scores)

    return n_vals[pos], scores[pos]


def k_decision(size, ymin, ymax, n, k, x):
    # performs RBS operation
    # returns the best (lowest) score and combination of operations that modify K
    original_score = score_circuit(size, ymin, ymax, n, k, x)

    k_vals = []
    scores = []
    k_vals[0] = k
    scores[0] = original_score

    for i in range(1):
        k_vals[i+1] = rbs(k, x, i)
        scores[i+1] = score_circuit(size, ymin, ymax, n, k_vals[i+1], x)

    pos = compare(scores)

    return k_vals[pos], scores[pos]


def best_score(size, ymin, ymax, n, k, x):
    # (1) retrieve outputs of y, n, and k_decision
    # (2) add scores to list scores[]
    # (3) generate scores from all combinations of each new parameter
    # (4) find and return the best score
    ymax_new, ymin_new, y_decision_score = y_decision(
        size, ymin, ymax, n, k, x)
    n_new, n_decision_score = n_decision(size, ymin, ymax, n, k, x)
    k_new, k_decision_score = k_decision(size, ymin, ymax, n, k, x)

    scores = [y_decision_score, n_decision_score, k_decision_score]

    ymaxs = [ymax, ymax_new]
    ymins = [ymin, ymin_new]
    ns = [n, n_new]
    ks = [k, k_new]

    params = 3
    num = 2**params  # 8 combinations

    for i in range(num-1):
        bin = format(i, "b")
        scores[i+3] = score_circuit(size, ymins[bin[0]],
                                    ymaxs[bin[0]], ns[bin[1]], ks[bin[2]])

   # best_score = min(i for i in scores if i > 0)
    pos = compare(scores)
    best_score = scores[pos]

    return best_score

# ======================================================================================
# ========================================MAIN==========================================
# ======================================================================================


def main():
    # Set our directory variables.
    in_dir = os.path.join(os.getcwd(), 'input')
    out_dir = os.path.join(os.getcwd(), 'output')

    # Set our input files.
    chassis_name = 'Eco1C1G1T1'
    in_ucf = f'{chassis_name}.UCF.json'
    # options = 'options.csv'
    input_sensor_file = f'{chassis_name}.input.json'
    output_device_file = f'{chassis_name}.output.json'

    # Open then parse .json files
    in_param = read_file(input_sensor_file)
    UCF_param = read_file(in_ucf)

    # Get user input of operations.
    operation = input("Choose up to 4 operations from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
    operation = operation.split()
    while len(operation) > 4 or len(operation) == 0:
        print("Incorrect entry of operations. Try again.\n")
        operation = input("Choose up to 4 operations from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
        operation = operation.split()

    # xhi = in_param[1]
    # xlow = in_param[2]
    # x = [xhi, xlow]
    # ymax = UCF_param[1]
    # ymin = UCF_param[2]
    # k = UCF_param[3]
    # n = UCF_param[4]
    if x > 1.05:
        sys.exit("Invalid x value\n")
    # Call functions based on truth table.
    for i in range(len(operation)):
        match operation[i]:
            case 'a':
                stretch(x, ymax, ymin)
            case 'b':
                slope(n, x, 1)
            case 'c':
                slope(n, x, 0)
            case 'd':
                promoter(x, ymax, ymin, 1)
            case 'e':
                promoter(x, ymax, ymin, 0)
            case 'f':
                rbs(k, x, 1)
            case 'g':
                rbs(k, x, 0)
            case 'x':
                break


if __name__ == "__main__":
    main()


# TO DO:
# fix parse
# match input to fxn call
# our scoring
# api scoring
# decide what we want to print (scores, sequence)

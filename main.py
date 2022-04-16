"""
EC552 HW 1: Genetic Circuit Design Program
Drew Gross and Marlee Feltham
"""

import os

import sys
import math
# from celloapi2 import CelloQuery, CelloResult
import json


# circuit: NOT and NOR gate
# pLuxStar----- NOT ________
#              pTet ________|---NOR--output
# not input: pluxstar
# not gate: P3_PhlF
# nor input a: not output
# nor input b: ptet
# nor gate: A1_AmtR

# ======================================================================================
# ====================================== FILE R/W ======================================
# ======================================================================================


def read_file(fname, chassis_name):
    # open and read .input and .UCF JSON files
    with open('input/' + fname, 'r') as file:
        content = file.read()
    data = json.loads(content)
    if fname == f'{chassis_name}.input.json':
        inputs = parse_input(data)
        file.close()
        return inputs
        # data_list = parse_UCF(data)
    elif fname == f'{chassis_name}.UCF.json':
        UCF = {}
        UCF = parse_UCF(data)
        file.close()
        return UCF
        # data_list = parse_input(data)


def write_output(fname, data):
    # open and write to .output JSON file
    with open(fname, 'w') as file:
        json.dump(data, file)
    file.close()

# ======================================================================================
# ======================================= PARSER =======================================
# ======================================================================================


def parse_UCF(data):
    UCF = {}
    for i in range(len(data)):
        if data[i]["collection"] == 'models':
            if data[i]['name'] == "P3_PhlF_model" or data[i]['name'] == "A1_AmtR_model":
                UCF['name'] = data[i]['name']
                if (data[i]["collection"] == 'parameters'):
                    for j in range(len(data[i]['name'])):
                        if data[i]['name'] == 'ymax':
                            UCF['ymax'] = data[i]['name']['value']
                        elif data[i]['name'] == 'ymin':
                            UCF['ymax'] = data[i]['name']['value']
                        elif data[i]['name'] == 'K':
                            UCF['k'] = data[i]['name']['value']
                        elif data[i]['name'] == 'n':
                            UCF['n'] = data[i]['name']['value']
    return UCF


def parse_input(data):
    inputs = {}
    name = []
    ymax = []
    ymin = []
    for i in range(len(data)):
        if data[i]['collection'] == 'models':
            if data[i]['name'] == 'LuxR_sensor_model' or data[i]['name'] == 'TetR_sensor_model':
                name.append(data[i]['name'])
                for j in range(len(data[i]['parameters'])):
                    if data[i]['parameters'][j]['name'] == 'ymax':
                        ymax.append(data[i]['parameters'][j]['value'])
                    elif data[i]['parameters'][j]['name'] == 'ymin':
                        ymin.append(data[i]['parameters'][j]['value'])
    inputs['name'] = name
    inputs['ymin'] = ymin
    inputs['ymax'] = ymax
    return inputs

# ======================================================================================
# ===================================== OPERATIONS =====================================
# ======================================================================================


def stretch(inputs, UCF):
    # apply stretch operation
    newUCF = UCF.copy()

    newUCF['ymax'] = UCF['ymax']*inputs['x']
    newUCF['ymin'] = UCF['ymin']/inputs['x']

    return newUCF


def promoter(inputs, UCF, pick):
    # apply weaker or stronger promoter operation based on pick value
    # pick == 0 -> weaker promoter
    # pick == 1 -> stronger promoter
    newUCF = UCF.copy()

    if pick == 0:
        newUCF['ymax'] = UCF['ymax']/inputs['x']
        newUCF['ymin'] = UCF['ymin']/inputs['x']
    elif pick == 1:
        newUCF['ymax'] = UCF['ymax']*inputs['x']
        newUCF['ymin'] = UCF['ymin']*inputs['x']

    return newUCF


def slope(UCF, inputs, pick):
    # apply increase or decrease slope operation based on pick value
    # pick == 0 -> decrease slope
    # pick == 1 -> increase slope
    newUCF = UCF.copy()

    if inputs['x'] <= 1.05:
        if pick == 0:
            newUCF['n'] = UCF['n']/inputs['x']
        elif pick == 1:
            newUCF['n'] = UCF['n']*inputs['x']
    else:
        sys.exit("Invalid x value\n")

    return newUCF


def rbs(UCF, inputs, pick):
    # apply weaker or stronger RBS operation based on pick value
    # pick == 0 -> weaker rbs
    # pick == 1 -> stronger rbs
    newUCF = UCF.copy()

    if pick == 0:
        newUCF['k'] = UCF['k']*inputs['x']
    elif pick == 1:
        newUCF['k'] = UCF['k']/inputs['x']

    return newUCF


# ======================================================================================
# =================================== SCORE CIRCUIT ====================================
# ======================================================================================

def response_function(UCF, inputs):
    # generate response function given parameters
    response = UCF['ymin'] + (UCF['ymax']-UCF['ymin']) / \
        (1+(inputs['x']/UCF['k']) ^ UCF['n'])

    return response


def score_circuit(size, UCF, inputs):
    # construct truth table for gate
    ttable = [0]*size
    for i in range(0, size-1):
        ttable[i] = response_function(
            UCF['ymin'], UCF['ymax'], UCF['n'], UCF['k'], inputs['x'])

    on_min = ttable[0]

    # score the gate based on the truth table
    if size == 2:
        off_max = max(ttable[1])
    else:
        off_max = max(ttable[1, 3])

    score = math.log10(on_min/off_max)

    return score


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
    # v_file = 'logic.v'
    options = 'options.csv'
    input_sensor_file = f'{chassis_name}.input.json'
    output_device_file = f'{chassis_name}.output.json'
    # q = CelloQuery(
    #     input_directory=in_dir,
    #     output_directory=out_dir,
    #     verilog_file=None,
    #     compiler_options=options,
    #     input_ucf=in_ucf,
    #     input_sensors=input_sensor_file,
    #     output_device=output_device_file,
    # )

    # Open then parse .json files
    in_param = read_file(input_sensor_file, chassis_name)
    UCF_param = read_file(in_ucf, chassis_name)

    # Get user input of operations.
    operation = input("Choose up to 4 operations from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
    operation = operation.split()
    while len(operation) > 4 or len(operation) == 0:
        print("Incorrect entry of operations. Try again.\n")
        operation = input("Choose up to 4 operations from the following list:\n(a) Stretch\n(b) Increase slope\n(c) Decrease slope\n(d) Stronger promoter\n(e) Weaker promoter\n(f) Stronger RBS\n(g) Weaker RBS\n(x) done\n")
        operation = operation.split()

    xhi = in_param[1]
    xlow = in_param[2]
    x = [xhi, xlow]
    ymax = UCF_param[1]
    ymin = UCF_param[2]
    k = UCF_param[3]
    n = UCF_param[4]

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
    # if logic_gate = 'NOT':
    #     size = 2
    # elif logic_gate = 'NOR':
    #     size = 4
    our_score = score_circuit(size, ymin, ymax, n, k, x)
    # Submit our query to Cello. This might take a second.
    # q.get_results()
    # # Fetch our Results.
    # res = CelloResult(results_dir=out_dir)
    # print(res.circuit_score)


if __name__ == "__main__":
    main()


# TO DO:
# fix parse
# match input to fxn call
# our scoring
# api scoring
# decide what we want to print (scores, sequence)

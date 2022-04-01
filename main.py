import os
import sys
import math
from celloapi2 import CelloQuery, CelloResult
import json

# ======================================================================================
# ====================================== FILE R/W ======================================
# ======================================================================================

def get_file(fname):
    # open and read .input and .UCF JSON files
    with open(fname, 'r') as file:
        content = file.read()

    data = json.loads(content)
    # do something here
    # send data to another function
    file.close()


def write_output(fname, data):
    # open and write to .output JSON file
    with open(fname, 'w') as file:
        json.dump(data, file)
    file.close()

# ======================================================================================
# ======================================= PARSER =======================================
# ======================================================================================


def parse_UCF(data):
    # parse .UCF JSON and store parameters in corresponding lists
    name = []
    ymax = []
    ymin = []
    K = []
    n = []
    alpha = []
    beta = []

    for i in range(len(data)):
        if data[i]["collections"] == 'models':
            name.append(data[i]['name'])
            if (data[i]["collections"] == 'parameters'):
                for j in range(len(data[i][name])):
                    if data[i][name][j]['name'] == 'ymax':
                        ymax.append(data[i][name][j]['value'])
                    elif data[i][name][j]['name'] == 'ymin':
                        ymin.append(data[i][name][j]['value'])
                    elif data[i][name][j]['name'] == 'K':
                        K.append(data[i][name][j]['value'])
                    elif data[i][name][j]['name'] == 'n':
                        n.append(data[i][name][j]['value'])
                    elif data[i][name][j]['name'] == 'alpha':
                        alpha.append(data[i][name][j]['value'])
                    elif data[i][name][j]['name'] == 'beta':
                        beta.append(data[i][name][j]['value'])

    return name, ymax, ymin, K, n, alpha, beta


def parse_input(data):
    # parse .input JSON and store parameters in corresponding lists
    name = []
    ymax = []
    ymin = []
    alpha = []
    beta = []

    for i in range(len(data)):
        if data[i]['collections'] == 'models':
            name.append(data[i]['name'])
            for j in range(len(data[i][name])):
                if data[i][name][j]['name'] == 'ymax':
                    ymax.append(data[i][name][j]['value'])
                elif data[i][name][j]['name'] == 'ymin':
                    ymin.append(data[i][name][j]['value'])
                elif data[i][name][j]['name'] == 'alpha':
                    alpha.append(data[i][name][j]['value'])
                elif data[i][name][j]['name'] == 'beta':
                    beta.append(data[i][name][j]['value'])

    return name, ymax, ymin, alpha, beta


# ======================================================================================
# ===================================== OPERATIONS =====================================
# ======================================================================================

def stretch(x, ymax, ymin):
    # apply stretch operation
    ymax_new = ymax*x
    ymin_new = ymin/x
    return ymax_new, ymin_new


def promoter(x, ymax, ymin, pick):
    # apply weaker or stronger promoter operation based on pick value
    # pick == 0 -> weaker promoter
    # pick == 1 -> stronger promoter
    if pick == 0:
        ymax_new = ymax/x
        ymin_new = ymin/x
    elif pick == 1:
        ymax_new = ymax*x
        ymin_new = ymin*x
    return ymax_new, ymin_new


def slope(n, x, pick):
    # apply increase or decrease slope operation based on pick value
    # pick == 0 -> decrease slope
    # pick == 1 -> increase slope
    if x <= 1.05:
        if pick == 0:
            n_new = n/x
        elif pick == 1:
            n_new = n*x
    else:
        sys.exit("Invalid x value\n")
    return n_new


def rbs(k, x, pick):
    # apply weaker or stronger RBS operation based on pick value
    # pick == 0 -> weaker rbs
    # pick == 1 -> stronger rbs
    if pick == 0:
        k_new = k*x
    elif pick == 1:
        k_new = k/x
    return k_new


# ======================================================================================
# =================================== SCORE CIRCUIT ===================================
# ======================================================================================

def response_function(ymin, ymax, n, k, x):
    # generate response function given parameters
    response = ymin + (ymax-ymin)/(1+(x/k) ^ n)
    return response


size = 2  # change later


def score_circuit(size, ymin, ymax, n, k, x):
    # construct truth table for gate
    ttable = [0]*size
    for i in range(0, size-1):
        ttable[i] = response_function(ymin, ymax, n, k, x[i])

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


def y_decision(size, ymin, ymax, n, k, x):
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
    ymax_r0[0] = ymax
    ymin_r0[0] = ymin
    # score the circuit with the original ymax and ymin values
    original_score = score_circuit(size, ymin, ymax, n, k, x)

    # apply stretch operation and score the circuit
    # # compare the original and stretch scores
    # r1 indicates which score is the lowest
    ymax_r0[1], ymin_r0[1] = stretch(x, ymax, ymin)
    stretch_score = score_circuit(size, ymin_r0[1], ymax_r0[1], n, k, x)

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
# ======================================================================================


def main():
    # Set our directory variables.
    in_dir = os.path.join(os.getcwd(), 'input')
    out_dir = os.path.join(os.getcwd(), 'output')

    # Set our input files.
    chassis_name = 'Eco1C1G1T1'
    in_ucf = f'{chassis_name}.UCF.json'
    v_file = 'and.v'
    options = 'options.csv'
    input_sensor_file = f'{chassis_name}.input.json'
    output_device_file = f'{chassis_name}.output.json'
    q = CelloQuery(
        input_directory=in_dir,
        output_directory=out_dir,
        verilog_file=v_file,
        compiler_options=options,
        input_ucf=in_ucf,
        input_sensors=input_sensor_file,
        output_device=output_device_file,
    )

    # Submit our query to Cello. This might take a second.
    q.get_results()
    # Fetch our Results.
    res = CelloResult(results_dir=out_dir)
    print(res.circuit_score)


if __name__ == "__main__":
    main()

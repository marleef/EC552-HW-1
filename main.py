import os
import sys
import math
from celloapi2 import CelloQuery, CelloResult
import json


def get_file(fname):
    with open(fname, 'r') as file:
        content = file.read()

    data = json.loads(content)
    # do something here
    # send data to another function
    file.close()

# ======================================================================================
# ======================================= PARSER =======================================
# ======================================================================================


def parse_UCF(data):
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
    ymax_new = ymax*x
    ymin_new = ymin/x
    return ymax_new, ymin_new


def promoter(x, ymax, ymin, pick):
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
    response = ymin + (ymax-ymin)/(1+(x/k) ^ n)
    return response


size = 2  # change later


def score_circuit(size, ymin, ymax, n, k, x):
    ttable = [0]*size
    for i in range(0, size-1):
        ttable[i] = response_function(ymin, ymax, n, k, x[i])

    on_min = ttable[0]
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
    position = scores.index(min(scores))
    return position


def y_decision(size, ymin, ymax, n, k, x):
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

    return ymax_new, ymin_new


def n_decision(size, ymin, ymax, n, k, x):
    original_score = score_circuit(size, ymin, ymax, n, k, x)  # global?

    n_vals = []
    scores = []
    n_vals[0] = n
    scores[0] = original_score

    for i in range(1):
        n_vals[i+1] = slope(n, x, i)
        scores[i+1] = score_circuit(size, ymin, ymax, n_vals[i+1], k, x)

    pos = compare(scores)

    return n_vals[pos]


def k_decision(size, ymin, ymax, n, k, x):
    original_score = score_circuit(size, ymin, ymax, n, k, x)  # global?

    k_vals = []
    scores = []
    k_vals[0] = k
    scores[0] = original_score

    for i in range(1):
        k_vals[i+1] = rbs(k, x, i)
        scores[i+1] = score_circuit(size, ymin, ymax, k_vals[i+1], k, x)

    pos = compare(scores)

    return k_vals[pos]

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

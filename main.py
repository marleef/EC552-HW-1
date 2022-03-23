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
# ======================================================================================
# individual/single list?


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
# ======================================================================================


def stretch(x, ymax, ymin):
    ymax_new = ymax*x
    ymin_new = ymin/x
    return ymax_new, ymin_new


def increase_slope(n, x):
    if x <= 1.05:
        n_new = n*x
    else:
        sys.exit("Invalid x value\n")
    return n_new


def decrease_slope(n, x):
    if x <= 1.05:
        n_new = n/x
    else:
        sys.exit("Invalid x value\n")
    return n_new


def stronger_promoter(x, ymax, ymin):
    ymax_new = ymax*x
    ymin_new = ymin*x
    return ymax_new, ymin_new


def weaker_promoter(x, ymax, ymin):
    ymax_new = ymax/x
    ymin_new = ymin/x
    return ymax_new, ymin_new


def strong_rbs(k, x):
    k_new = k/x
    return k_new


def weak_rbs(k, x):
    k_new = k*x
    return k_new

# ======================================================================================
# ======================================================================================


def response_function(ymin, ymax, n, k, x):
    response = ymin + (ymax-ymin)/(1+(x/k) ^ n)
    return response


size = 2


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

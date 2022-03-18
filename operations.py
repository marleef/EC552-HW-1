# gate operations

def stretch(x, ymax, ymin):
    ymax_new = ymax*x
    ymin_new = ymin/x

    return ymax_new, ymin_new

def increase_slope(n, x):
    if x <= 1.05:
        n_new = n*x

def decrease_slope(n, x):
    if x <= 1.05:
        n_new = n/x

def stronger_promoter(x, ymax, ymin):
    ymax_nw = ymax*x
    ymin_new = ymin*x

def weaker_promoter(x, ymax, ymin):
    ymax_nw = ymax/x
    ymin_new = ymin/x 

def strong_rbs(k, x):
    k_new = k/x

def weak_rbs(k, x):
    k_nw = k*x


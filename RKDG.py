"""
Solve eqs using P^1 RKDG(RK2)
"""

import matplotlib.pyplot as plt
import sympy as sy
import numpy as np 
from scipy import integrate 

x = sy.Symbol('x')

def Euler_forward(step_len, matrices, u_list, deri_func):
    u_list_new = []
    deri_vec = deri_func(matrices, u_list)
    for num in range(len(u_list)):
        u_list_new.append(u_list[num] + step_len * deri_vec[num])
    return u_list_new

def SSP_RK_2nd(step_len, matrices, u_list, deri_func):
    u_temp_list = Euler_forward(step_len, matrices, u_list, deri_func)
    u_temp_list = Euler_forward(step_len, matrices, u_temp_list, deri_func)
    u_list_new = [0.5 * u_temp_list[num] + 0.5 * u_list[num] for num in range(len(u_list))]
    return u_list_new

def cut(left_bdry, mesh_len, num_cuts):
    len_interval = mesh_len / num_cuts
    cuts = [(left_bdry + num * len_interval, left_bdry + (num + 1) * len_interval) for num in range(num_cuts)]
    return cuts 

def compute_basis(left_bdry, mesh_len, num_cuts):
    basis = []
    len_interval = mesh_len / num_cuts
    for num in range(num_cuts):
        mid_point = left_bdry + (num + 0.5)*len_interval
        basis.append([1.0 + 0 * x, x - mid_point])
    return basis

def project(cuts, basis, f):
    return [np.mat([sy.integrate(phi * f, (x, cuts[num][0], cuts[num][1])).evalf() / sy.integrate(phi**2, (x, cuts[num][0], cuts[num][1]))
        for phi in basis[num]]).T for num in range(len(cuts))] 

def product(f, g, interval):
    return float(sy.integrate(f * g, (x, interval[0], interval[1])))


def combine_expr(u_list, basis):
    u_expr_list = []
    for num in range(len(u_list)):
        base_vec = np.mat(basis[num]).T
        u_expr = (u_list[num].T * base_vec).tolist()[0][0]
        u_expr_list.append(u_expr)
    return u_expr_list

def compute_error(u_target, u_expr_list, cuts):
    error_square = 0.0
    for num in range(len(cuts)):
        interval = cuts[num]
        error_square_num = integrate.quad(sy.lambdify(x, (u_expr_list[num] - u_target)**2), interval[0], interval[1])
        error_square += error_square_num[0]
    return np.sqrt(float(error_square))

def RKDG(u_0, u_1, num_cuts, matrices_func, deri_func, time_total, time_step_len):
    cuts = cut(num_cuts)
    basis = compute_basis(num_cuts)
    matrices = matrices_func(cuts, basis)
    u_list = project(cuts, basis, u_0)


    time = 0.0
    u_expr_list = combine_expr(u_list, basis) 
    u_t = sy.sin(2.0*sy.pi*(x-time))
    error = compute_error(u_t, u_expr_list, cuts)    
    #print "time = " +  str(time) + "\n"
    #print "L2 err of u is " + str(error) + "\n"
    
    while time + time_step_len  < time_total:
        u_list = SSP_RK_2nd(time_step_len, matrices, u_list, deri_func) 
        u_expr_list = combine_expr(u_list, basis) 
        time += time_step_len 
        u_t = sy.sin(2.0*sy.pi*(x-time))
        error = compute_error(u_t, u_expr_list, cuts)    
        #print "time = " +  str(time) + "\n"
        #print "L2 err of u is " + str(error) + "\n"
    
    time_left = time_total - time
    if time_left > 0:
        u_list = SSP_RK_2nd(time_left, matrices, u_list, deri_func)
    
    u_expr_list = combine_expr(u_list, basis) 
    error = compute_error(u_1, u_expr_list, cuts)    
     
    return (num_cuts, u_expr_list, error)

     

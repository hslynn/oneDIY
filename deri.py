import getopt
import sys
from RKDG import *
def compute_matrices(cuts, basis):
    product_base_base_inv = []
    product_deri_base_base = []
    product_deri_deri = []
    boundary_right_plus = []
    boundary_right_minus = []
    boundary_left_plus = []
    boundary_left_minus = []
    
    for num in range(len(cuts)):
        interval = cuts[num]
        interval_left = cuts[num - 1]
        interval_right = cuts[(num+1)%len(cuts)]
        base = basis[num]
        base_left = basis[num - 1]
        base_right = basis[(num+1)%len(cuts)]

        product_base_base_inv.append(np.linalg.inv(np.mat([[product(phi_row, phi_col, interval) for phi_col in base] for phi_row in base])))
    
        deri_base = [sy.diff(phi, x) for phi in base]
        product_deri_base_base.append(np.mat([[product(phi_x, phi, interval) for phi in base] for phi_x in deri_base]))
        product_deri_deri.append(np.mat([[product(phi_x_row, phi_x_col, interval) for phi_x_col in deri_base] for phi_x_row in deri_base]))
        
        boundary_right_plus.append(np.mat([[float(phi_row.subs(x, interval[1]) * phi_col.subs(x, interval_right[0])) for phi_col in base_right] for phi_row in base]))
        boundary_right_minus.append(np.mat([[float(phi_row.subs(x, interval[1]) * phi_col.subs(x, interval[1])) for phi_col in base] for phi_row in base]))
        
        boundary_left_plus.append(np.mat([[float(phi_row.subs(x, interval[0]) * phi_col.subs(x, interval[0])) for phi_col in base] for phi_row in base]))
        boundary_left_minus.append(np.mat([[float(phi_row.subs(x, interval[0]) * phi_col.subs(x, interval_left[1])) for phi_col in base_left] for phi_row in base]))

    return (product_base_base_inv, product_deri_base_base, product_deri_deri, boundary_right_plus, boundary_right_minus, boundary_left_plus, boundary_left_minus)

def get_deri_plus(matrices, u_list):
    deri_list = []
    num_cuts = len(u_list)
    for num in range(num_cuts):
        u = u_list[num]
        if num == num_cuts - 1:
            deri_vec = matrices[0][num] * (-matrices[1][num]*u + matrices[4][num]*u - matrices[5][num]*u)
        else:
            u_right = u_list[num+1]
            deri_vec = matrices[0][num] * (-matrices[1][num]*u + matrices[3][num]*u_right - matrices[5][num]*u)
        deri_list.append(deri_vec)
    return deri_list

def get_deri_minus(matrices, u_list):
    deri_list = []
    num_cuts = len(u_list)
    for num in range(num_cuts):
        u = u_list[num]
        if num == 0:
            deri_vec = matrices[0][num] * (-matrices[1][num]*u + matrices[4][num]*u - matrices[5][num]*u)
        else:
            u_left = u_list[num-1]
            deri_vec = matrices[0][num] * (-matrices[1][num]*u + matrices[4][num]*u - matrices[6][num]*u_left)
        deri_list.append(deri_vec)
    return deri_list

def main():
    mesh_num = 100
    DG_degree = 1
    mesh_len = 3.0
    left_bdry = 0.5
    opts, dumps = getopt.getopt(sys.argv[1:], "-m:-d:-l:")
    for opt, arg in opts:
        if opt == "-m":
            mesh_num = int(arg)
        if opt == "-d":
            DG_degree = int(arg)
        if opt == "-i":
            left_bdry = float(arg)
        if opt == "-l":
            mesh_len = float(arg)


    num_cuts = mesh_num
    u = 2/x 
    p = -2/x**2 

    cuts = cut(left_bdry, mesh_len, num_cuts)
    basis = compute_basis(left_bdry, mesh_len, num_cuts)
    matrices = compute_matrices(cuts, basis)

    u_list = project(cuts, basis, u)
    u_expr_list = combine_expr(u_list, basis)  

    p_list = project(cuts, basis, p)
    p_expr_list = combine_expr(p_list, basis)

    p_plus_list = get_deri_plus(matrices, u_list)
    p_plus_expr_list = combine_expr(p_plus_list, basis)

    p_minus_list = get_deri_minus(matrices, u_list)
    p_minus_expr_list = combine_expr(p_minus_list, basis)

    dif_plus_expr_list = [p_plus_expr_list[idx] - p_expr_list[idx] for idx in range(len(u_expr_list))]
    dif_minus_expr_list = [p_minus_expr_list[idx] - p_expr_list[idx] for idx in range(len(u_expr_list))]

    #plot error
    for num in range(len(cuts)):
        interval = cuts[num]
        t = np.arange(interval[0], interval[1], 0.002) 
        dif_plus_expr = dif_plus_expr_list[num]
        dif_minus_expr = dif_plus_expr_list[num]
        p_expr = p_expr_list[num]
        p_plus_expr = p_plus_expr_list[num]
        p_minus_expr = p_minus_expr_list[num]
    
        dif_deri_expr = p_plus_expr - p_minus_expr

        plt.subplot(3, 2, 1)
        plt_error_plus = plt.plot(t, map(sy.lambdify(x, dif_plus_expr), t), 'b')

        plt.subplot(3, 2, 2)
        plt_error_minus = plt.plot(t, map(sy.lambdify(x, dif_minus_expr), t), 'b')

        plt.subplot(3, 2, 3)
        plt_p_plus = plt.plot(t, map(sy.lambdify(x, p_plus_expr), t), 'r')
        plt_p = plt.plot(t, map(sy.lambdify(x, p_expr), t), 'y')


        plt.subplot(3, 2, 4)
        plt_p_minus = plt.plot(t, map(sy.lambdify(x, p_minus_expr), t), 'r')
        plt_p = plt.plot(t, map(sy.lambdify(x, p_expr), t), 'y')

        plt.subplot(3, 1, 3)
        plt_dif_deri = plt.plot(t, map(sy.lambdify(x, dif_deri_expr), t), 'r')

    plt.show()

main()

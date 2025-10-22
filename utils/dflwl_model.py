# this file includes functions to implement the modified dFLWL model to utilize the BMA-PDFs of forecasted inflow within the forecast horizon (i.e., I2) 

from sympy.solvers import solve
from sympy import Symbol
import sympy


#%%%%%% define the Partial Differential Equation (PDE) of the water conservation loss, 
#                  which can be considered as the Marginal Value (MV) of carryover storage
def f_L1_W1(W1, W_N, w, m):
    f1_W1 = -(-(1-w)*m/W_N*((W_N-W1)/W_N)**(m-1))
    return(f1_W1)

#%%%%%% define functions to obtain the equal Marginal Value (MV) solution
# W1_c: critical W1 value [note: the regression coeff. is respect to W1 instead of I2]
W1_c = Symbol('W1_c')

# a supporting function (used in solve_equal_MV function)
def f_L1(W_N, w, m):
    return((1-w)*m/W_N*((W_N-W1_c)/W_N)**(m-1))

# a supporting function (used in solve_equal_MV function)
def reg_pdf(coeffcients, w):
    return(w*(coeffcients[0]*W1_c + coeffcients[1]*W1_c**2 + coeffcients[2]*W1_c**3 + coeffcients[3]))

# solve the optimal solutions of the equal marginal values
def solve_equal_MV(reg_coef, I2_expected, delta_min, I2_995, date_R2_max, W_N, w, m):
    W1_eq = solve(reg_pdf(reg_coef, w) - f_L1(W_N, w, m), W1_c)
    #print("solve_equal_MV-W1")
    if sympy.im(W1_eq[0])!= 0:
        W1_eq[0] = complex_filter(W1_eq, I2_expected, delta_min, I2_995, date_R2_max)
    delta_op = date_R2_max - I2_expected - W1_eq[0]
    return W1_eq[0], delta_op



#%%%%%% define a supporting function to ignore the imaginary part of 'solve' function solutions
def complex_filter(W1_i, I2_expected, delta_min, I2_995, date_R2_max):
    final_W1 = -9999
    W1_up_bound = date_R2_max - (I2_expected+delta_min) +5 # +5 to account for approximation errors
    W1_lo_bound = date_R2_max - I2_995

    for item in W1_i:
        if sympy.im(item)<1e-6:
            if (sympy.re(item)<W1_up_bound) and (sympy.re(item)>W1_lo_bound):
                final_W1 = sympy.re(item)
    if final_W1 == -9999:
        if (sympy.re(W1_i[-1])<W1_lo_bound):
            final_W1 = W1_lo_bound
        if (sympy.re(W1_i[0])>W1_lo_bound):
            final_W1 = W1_up_bound
    return(final_W1)



#%%%%%% compute daily-updated critical inflow levels (I2_a, I2_0, I2_e, I2_h) 
# used to determine which optimal case of KKT conditions to apply
# 1. critical inflow where delta = delta_min
def calculate_I2_a(date_R2_max, delta_min):
    I2_a = date_R2_max - delta_min
    return(I2_a)

# 2. critical inflow where f2(delta_0) = f1(0) 
# implying f2(delta_min)> f1(0)> f1(W_max)
def calculate_I2_0(reg_coef, f1_W1_0, I2_expected, delta_min, I2_995, date_R2_max, w):
    W1_0 = solve(reg_pdf(reg_coef, w) - f1_W1_0, W1_c)
    if len(W1_0) == 0:
        I2_0 = date_R2_max-(I2_995-I2_expected)
        return(I2_0)
    
    if sympy.im(W1_0[0])!= 0: # complex filter
        W1_0[0] = complex_filter(W1_0, I2_expected, delta_min, I2_995, date_R2_max)
    I2_0 = date_R2_max-(date_R2_max-W1_0[0]-I2_expected)-0 #W1=0
    return(I2_0)

# 3. critical inflow where f2(delta_e) = f1(W_max)
# implying f2(delta_min)> f1(W_max)
def calculate_I2_e(reg_coef, f1_W1_max, I2_expected, delta_min, I2_995, f2_delta_995, W_max, date_R2_max, w):
    # determine whether the f2(I_99.5th) >= f1_W1_max
    # if so, let I2_e = I_99.5th, and no need to consider larger inflow
    if f2_delta_995 >= f1_W1_max:
        I2_e = date_R2_max-(I2_995-I2_expected)-W_max
    else:
        W1_e = solve(reg_pdf(reg_coef, w) - f1_W1_max, W1_c)
        
        if sympy.im(W1_e[0])!= 0: # complex filter
            W1_e[0] = complex_filter(W1_e, I2_expected, delta_min, I2_995, date_R2_max)
        if W1_e[0] >= date_R2_max: # no feasible solution for f2(delta_e) = f1(W_max)
            I2_e = I2_expected-1
        else:
            I2_e = date_R2_max-(date_R2_max-W1_e[0]-I2_expected)-W_max
    return(I2_e)

# 4. critical inflow where f2(delta_min) = f1(W1_h) 
# implying f1(W_max) <=f2(delta_min) <=f1(0)
def calculate_I2_h(f2_delta_min, delta_min, date_R2_max, W_N, w, m):
    W1_h = solve(f_L1(W_N, w, m) - f2_delta_min, W1_c)
    I2_h = date_R2_max-delta_min-W1_h[0]
    return(I2_h)


#%%%%%% apply KKT conditions to derive daily optimal hedging rules
def daily_optimal_KKT(I2_expected, delta_min, f2_delta_min, I2_995, f2_delta_995, reg_coef, R2_max, W_max, W_N, w, m, f1_W1_0, f1_W1_max):
    #############################
    # prepare: compute I2_a, I2_0, I2_e, I2_h at the begining of each time step
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # case 1: max MV of water conservation is smaller than max MV of flood control
    if f2_delta_min >= f1_W1_0:
        I2_0 = calculate_I2_0(reg_coef, f1_W1_0, I2_expected, delta_min, I2_995, R2_max, w)
        if I2_0 > 0:
            I2_e = calculate_I2_e(reg_coef, f1_W1_max, I2_expected, delta_min, I2_995, f2_delta_995, W_max, R2_max, w)
        
        if I2_expected > I2_0:
            case_num = '1.1'
            W1_optimal = 0
            delta_optimal = R2_max - I2_expected
        elif I2_expected < I2_e:
            case_num = '1.3' 
            W1_optimal = W_max
            delta_optimal = R2_max - I2_expected - W_max
        else:
            case_num = '1.2'
            if f_L1_W1((R2_max-I2_995), W_N, w, m)<f2_delta_995:
                W1_optimal = R2_max-I2_995
                delta_optimal = I2_995-I2_expected
            else:    
                W1_optimal, delta_optimal = solve_equal_MV(reg_coef,I2_expected, delta_min, I2_995, R2_max, W_N, w, m)
                
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # case 2:
    elif f1_W1_0 >= f2_delta_min and f1_W1_max <= f2_delta_min:
        I2_a = calculate_I2_a(R2_max, delta_min)
        if I2_expected > (I2_a-5): # -5 to account for numerical approximation error
            # under this case, I2_h and I2_e cannot be calculated
            case_num = '2.1'
            W1_optimal = 0
            delta_optimal = R2_max - I2_expected
            return (W1_optimal, delta_optimal, case_num)
        
        I2_h = calculate_I2_h(f2_delta_min, delta_min, R2_max, W_N, w, m)
        I2_e = calculate_I2_e(reg_coef, f1_W1_max, I2_expected, delta_min, I2_995, f2_delta_995, W_max, R2_max, w)
        
        if (I2_expected > I2_h) and (I2_expected < I2_a):
            case_num = '2.2'
            W1_optimal = R2_max - I2_expected - delta_min
            delta_optimal = delta_min
        elif (I2_expected > I2_e) and (I2_expected <= I2_h):
            case_num = '2.3'
            W1_optimal, delta_optimal = solve_equal_MV(reg_coef,I2_expected, delta_min, I2_995, R2_max, W_N, w, m)
        else:
            case_num = '2.4'
            W1_optimal = W_max
            delta_optimal = R2_max - I2_expected - W_max
            

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    # case 3: MV of water conservation is always larger than MV of flood control 
    elif f1_W1_max >= f2_delta_min:
        I2_a = calculate_I2_a(R2_max, delta_min)
        
        if I2_expected > I2_a:
            W1_optimal = 0
            delta_optimal = R2_max - I2_expected
            case_num = '3.1'
        elif I2_expected < I2_a - W_max:
            W1_optimal = W_max
            delta_optimal = R2_max - I2_expected - W_max
            case_num = '3.3'
        else:
            W1_optimal = R2_max - I2_expected - delta_min
            delta_optimal = delta_min
            case_num = '3.2'
            

    ############################# 
    return (W1_optimal, delta_optimal, case_num)
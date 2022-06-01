import numpy as np
import pandas as pd
import sympy as sp
from sympy import Interval
from sympy.calculus.util import minimum, maximum

# x = sp.Symbol('x', float = True)
#
# interv = Interval(0.1, 0.2)
# f = sp.atan(0.8*x + 0.2)
# g = sp.exp(2*x + 1)
# res_min_f = minimum(f, x, interv)
# res_max_f = maximum(f, x, interv)
# res_min_g = minimum(g, x, interv)
# res_max_g = maximum(g, x, interv)
#
# print(res_min_f, res_max_f)
# print(res_min_g, res_max_g)
#
# x = sp.Symbol('x', float = True)
# u_, v_ = sp.symbols('u_ v_')
# f_zv = (1+u_)**0.5*v_
#
# c_1 = sp.diff(f_zv, u_)
# c_2 = sp.diff(f_zv, v_)
#
# c_1 = c_1.subs(u_, res_min_f)
# c_1 = c_1.subs(v_, res_max_g)
# c_2 = c_2.subs(u_, res_max_f)
# c_2 = c_2.subs(v_, res_min_g)
#
# rez = ((1+f)**0.5)*g
#
# # Находим эпсилоны
# e_1 = 10**(-6)/(3*c_1)
# e_2 = 10**(-6)/(3*c_2)
# e_3 = 10**(-6)/3
#
# porog = 10**(-8)
# # ---------------------------------------------------------------------
# interv = Interval(2*0.1+1, 2*0.2+1)
# k_exp = 1
# chlen_exp = x**k_exp/sp.factorial(k_exp)
# myexp = 1+chlen_exp
# while maximum(chlen_exp, x, interv) - e_2 > porog: # Считаем номер члена в раложении exp для необходимой точности
#     k_exp = k_exp + 1
#     chlen_exp=chlen_exp*x/k_exp
#     myexp=myexp+chlen_exp
#
#
# interv = Interval(0.8*0.1+0.2, 0.8*0.2+0.2)
# k_arc = 0
# myarc = x**(2*k_arc+1)/(2*k_arc+1)
# chlen_arc=x**(2*k_arc+1)/(2*k_arc+1)
# while maximum(chlen_arc, x, interv) - e_1 > porog: # Считаем номер члена в раложении arctan для необходимой точности
#     k_arc = k_arc+1
#     chlen_arc=chlen_arc*x**2*(2*(k_arc-1)+1)/(2*k_arc+1)
#     myarc=myarc+(-1)**k_arc*chlen_arc
#
# # ---------------------------------------------------------------------
def my_Geron(x, eps):
  w_1 = 1
  w_2 = (w_1 + x/w_1)/2
  while np.abs(w_1 - w_2) > eps:
    w_1 = w_2
    w_2 = (w_1 + x/w_1)/2
  return w_2
  i=10
# # ---------------------------------------------------------------------
# eps=e_3
# for i in range(10, 20, 1):
#    while (abs(my_Geron(1+myarc.subs(x, 0.8 * i / 100 + 0.2), eps) * myexp.subs(x, 2 * i / 100 + 1) - z(i / 100)) - e_3 > porog):
#    eps=eps/10
# # ---------------------------------------------------------------------
def myfinarc(x):
    return x ** 13 / 13 - x ** 11 / 11 + x ** 9 / 9 - x ** 7 / 7 + x ** 5 / 5 - x ** 3 / 3 + x
def myfinexp(x):
    return x**12/479001600 + x**11/39916800 + x**10/3628800 + x**9/362880 + x**8/40320 + x**7/5040 + x**6/720 + x**5/120 + x**4/24 + x**3/6 + x**2/2 + x + 1
def z(x):
    z = np.power((1 + np.arctan(0.8*x + 0.2)), 1/2) * np.exp(2*x + 1)
    return z
# ---------------------------------------------------------------------

vect_approx = np.zeros(10)
for i in range(10, 20, 1):
    vect_approx[i-10] = my_Geron(1+myfinarc(0.8 * i / 100 + 0.2), 10**(-6)/3) * myfinexp(2 * i / 100 + 1)

vect_exact = np.zeros(10)
for i in range(10, 20, 1):
    vect_exact[i-10] = z(i / 100)

vect_dif = np.zeros(10)
for i in range(10, 20, 1):
    vect_dif[i-10] = abs(vect_approx[i-10] - vect_exact[i-10])

vect=np.zeros(10)
for i in range(10, 20, 1):
    vect[i-10] = i/100

result = pd.DataFrame({'x' : vect, 'f_exact' : vect_exact, 'f_approx' : vect_approx, 'error': vect_dif})
print(result)
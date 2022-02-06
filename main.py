import sympy as sym
import tkinter as tk
import time

x = sym.Symbol('x')

#evaluate x root using bisection method
def Bisection_eval(exp, xl, xu, es, imax):
    start = time.time()
    ea = None
    result = ''
    if exp.evalf(10,subs={x:xl})*exp.evalf(10,subs={x:xu}) > 0:
        print("no bisection")
        return
    for i in range(imax):
        xr_new = (xl + xu) / 2  #compute new root
        if i != 0:
            ea = abs((xr_new-xr_old)/xr_new) #compute relative error percentage
        xr_old = xr_new
        if exp.evalf(10,subs={x:xl})*exp.evalf(10,subs={x:xr_new}) > 0:
            print(str(i) + 'iteration:')
            st_i = str(i) + 'iteration:'
            print(' xlower: ' + str(xl))
            st_xl = ' xlower: ' + str(xl)
            print(' xupper: ' + str(xu))
            st_xu = ' xupper: ' + str(xu)
            xl = xr_new
            print(' xroot: ' + str(xr_new))
            st_xr = ' xroot: ' + str(xr_new)
            print(' eabs: ' + str(ea))
            st_ea = ' eabs: ' + str(ea)
            print("-------------------------")
            st_sp = "-------------------------" + "\n"
            result = result + st_i + st_xl + st_xu + st_xr + st_ea + st_sp
        elif exp.evalf(10,subs={x:xl})*exp.evalf(10,subs={x:xr_new}) < 0:
            print(str(i) + 'iteration:')
            st_i = str(i) + 'iteration:'
            print(' xlower: ' + str(xl))
            st_xl = ' xlower: ' + str(xl)
            print(' xupper: ' + str(xu))
            st_xu = ' xupper: ' + str(xu)
            xu = xr_new
            print(' xroot: ' + str(xr_new))
            st_xr = ' xroot: ' + str(xr_new)
            print(' eabs: ' + str(ea))
            st_ea = ' eabs: ' + str(ea)
            print("-------------------------")
            st_sp = "-------------------------" + "\n"
            result = result + st_i + st_xl + st_xu + st_xr + st_ea + st_sp
        else:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
        if ea is not None and ea < es:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
    end = time.time()
    result = result + "\n" + str(end - start)
    return result


#compute root of function using false position method


def False_Position(exp,xl,xu,es,imax):
    start = time.time()
    ea = None
    result = ''
    if exp.evalf(10,subs={x:xl}) > 0 and exp.evalf(10,subs={x:xu}) < 0:
        print("function has same sign at endpoints")
        return
    for i in range(imax):
        y = xl * exp.evalf(10,subs={x:xu}) - xu * exp.evalf(10,subs={x:xl})
        z = exp.evalf(10,subs={x:xu}) - exp.evalf(10,subs={x:xl})
        xr_new = y/z   #compute new root x = xl*f(xu)-xu*f(xl)/f(xu)-f(xl)
        if i != 0:    #if i is zero dont compute relative error because we have only 1 root
            ea = abs((abs(xr_new)-xr_old)/xr_new)
        xr_old = xr_new
        print(' f(XR): ' + str(exp.evalf(10,subs={x:xr_new})))
        if exp.evalf(10,subs={x:xr_new}) > 0:
            print(str(i) + 'iteration:')
            st_i = str(i) + 'iteration:'
            print(' xlower: ' + str(xl))
            st_xl = ' xlower: ' + str(xl)
            print(' xupper: ' + str(xu))
            st_xu = ' xupper: ' + str(xu)
            xu = xr_new
            print(' xroot: ' + str(xr_new))
            st_xr = ' xroot: ' + str(xr_new)
            print(' eabs: ' + str(ea))
            st_ea = ' eabs: ' + str(ea)
            print("-------------------------")
            st_sp = "-------------------------" + "\n"
            result = result + st_i + st_xl + st_xu + st_xr + st_ea + st_sp
        elif exp.evalf(10,subs={x:xr_new}) < 0:
            print(str(i) + 'iteration:')
            st_i = str(i) + 'iteration:'
            print(' xlower: ' + str(xl))
            st_xl = ' xlower: ' + str(xl)
            print(' xupper: ' + str(xu))
            st_xu = ' xupper: ' + str(xu)
            xl = xr_new
            print(' xroot: ' + str(xr_new))
            st_xr = ' xroot: ' + str(xr_new)
            print(' eabs: ' + str(ea))
            st_ea = ' eabs: ' + str(ea)
            print("-------------------------")
            st_sp = "-------------------------" + "\n"
            result = result + st_i + st_xl + st_xu + st_xr + st_ea + st_sp
        else:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
        if ea is not None and ea < es:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
    end = time.time()
    result = result + "\n" + str(end - start)
    return result


#compute root of function using fixed point method
def Fixed_Point(xr, expG, es, imax):
    result = ''
    start = time.time()
    ea = None
    for i in range(imax):
        xr_old = xr
        xr = expG.evalf(10, subs={x:xr}) #compute new root
        if xr != 0:
            ea = abs((xr - xr_old) / xr)
        print(str(i) + ' iteration:')
        st_i = str(i) + ' iteration:'
        print(' xroot: ' + str(xr))
        st_xr = ' xroot: ' + str(xr)
        print(' eabs: ' + str(ea))
        st_ea = ' eabs: ' + str(ea)
        print("-------------------------")
        st_sp = "-------------------------" + "\n"
        result = result + st_i + st_xr + st_ea + st_sp
        if xr == 0 or ea < es:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
    end = time.time()
    result = result + "\n" + str(end - start)
    return result


#compute root of function using newton raphson method
def Newton_Raphson(expr, xr, es, imax):
    start = time.time()
    result = ''
    diff_expr = sym.diff(expr)   #get f'(x)
    ea = None
    for i in range(imax):
        xr_old = xr
        xr = xr - (expr.evalf(10, subs={x:xr})) / (diff_expr.evalf(10, subs={x:xr})) #compute x(i+1) = x(i) - f(x)/f'(x)
        ea = abs((xr_old - xr)/xr)
        st_i = 'iteration' + str(i) + "\n"
        print(str(i) + ' iteration:')
        print(' xroot: ' + str(xr))
        st_xr ="xroot: " + str(xr)
        print(' eabs: ' + str(ea))
        st_ea = ' eabs: ' + str(ea) + "\n"
        print("-------------------------")
        st_sp = "-------------------------"
        result = result + st_i + st_xr + st_ea + st_sp
        if xr == 0 or ea < es:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
    end = time.time()
    result = result + "\n" + str(end - start)
    return result


#compute root of function using secant method
def Secant_Eval(expr, x0, x1, es, imax):
    ea = None
    start = time.time()
    result = ''
    xr_i1 = x1 # equivalent to x(i-1)
    xr_i0 = x0 # equivalent to x(i)
    for i in range(imax):
        #compute x(i+1) = x(i) - f(x(i-1)) * x(i-1) - x(i) / f(x(i-1) - f(x(i))
        xr = xr_i0 - (expr.evalf(10, subs={x:xr_i0}) * (xr_i1 - xr_i0)) / (expr.evalf(10, subs={x:xr_i1}) - expr.evalf(10, subs={x:xr_i0}))
        ea = abs((xr-xr_i0)/xr)
        xr_i1 = xr_i0 # make x(i-1) of next iteration = x(i) of current iteration
        xr_i0 = xr # make x(i) of next iteration = xr of current iteration
        print(str(i) + ' iteration:')
        st_i = str(i) + ' iteration:'
        print(' xroot: ' + str(xr))
        st_xr = ' xroot: ' + str(xr)
        print(' eabs: ' + str(ea))
        st_ea = ' eabs: ' + str(ea)
        print("-------------------------")
        st_sp = "-------------------------" + "\n"
        result = result + st_i + st_xr + st_ea + st_sp
        if xr == 0 or ea < es:
            end = time.time()
            result = result + "\n" + str(end - start)
            return result
    end = time.time()
    result = result + "\n" + str(end - start)
    return result

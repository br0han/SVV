import numpy as np

def trap_int (f, l_limit, u_limit): ### function to integrate, upper limit, lower limit

    result = ((u_limit - l_limit)/2)*(f(u_limit)+f(l_limit)) #trapezoidal rule, exact for linear stuff

    return (result)


def Simpson_int (f, l_limit, u_limit): ### function to integrate, upper limit, lower limit

    result = ((u_limit - l_limit)/6)*(f(u_limit)+ 4*f((u_limit + l_limit)/2) +f(l_limit)) #Simpson rule

    return result

def Simpson38_int (f, l_limit, u_limit): ### function to integrate, upper limit, lower limit
    h = (u_limit - l_limit)/3

    result = (3/8)*h*(f(l_limit) + 3*f(l_limit + h)+ 3*f(l_limit + 2*h) +f(u_limit)) #Simpson rule

    return (result)



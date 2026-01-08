using Plots 
using BenchmarkTools

"""
We consider the BVP 
    u''(x) - u(x) = f(x), 0 <= x <= 1 
    u(0) = c, u(1) = d 
"""


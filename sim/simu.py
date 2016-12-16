from scipy.optimize import brentq
from math import log
from math import sqrt

def NSbest_l(n):
    def f(l, k):
        return l + (log(1.3)+log(l) - log(log(l)))/2 - k     
    return brentq(f, 2, n-1, args=(n // 2))


def mybest_l(n):
    def f(l, k):
        return l + (log(n - l))/2 - k
    return brentq(f, 2, n-1, args=(n//2))

def newbest_l(n):
    def f(l, k):
        return l + (log(4*(n-l)) - log(3))/2 - k
    return brentq(f, 2, n-1, args=(n//2))

print("l : NS vs new")
print(NSbest_l(64), newbest_l(64), mybest_l(64))
print(NSbest_l(128), newbest_l(128), mybest_l(32))
print(NSbest_l(256), newbest_l(256), mybest_l(16))
print(NSbest_l(512), newbest_l(512), mybest_l(512))

def NScst(n):
    l = NSbest_l(n)
    return sqrt(1.3 * l / log(l))

def mycst(n):
    l = mybest_l(n)
    return sqrt(n-l)

def newcst(n):
    l = newbest_l(n)
    return sqrt(4/3 *(n-l))

def Jouxcst(n):
    return sqrt(n//2)

print("sqrt(q) : NS vs Joux vs new")
print(NScst(64), newcst(64), Jouxcst(64), mycst(64))
print(NScst(128), newcst(128), Jouxcst(128), mycst(128))
print(NScst(256), newcst(256), Jouxcst(256), mycst(256))
print(NScst(512), newcst(512), Jouxcst(512), mycst(512))

    





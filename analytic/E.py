#!/usr/bin/env python
import numpy as np
from math import pi

N = 3

def hat(p):
    return 2. * np.sin(p / 2)

def ring(p):
    return np.sin(p)

def check(p):
    return 2. * np.sin(p / 4)

def shat(p, x):
    return np.sin(p * x / 2)

def chat(p, x):
    return np.cos(p * (x + 0.5) / 2)

def tsqEs(L, T, x, c):
    n = [0,0,0,0]
    result = 0
    for n[0] in range(-T, T):
      for n[1] in range(-L/2, L/2):
        for n[2] in range(-L/2, L/2):
          for n[3] in range(-L/2, L/2):
            p = np.array([ni * 2 * pi / LL for ni, LL in zip(n, (T, L, L, L))])
            psqd = np.sum((hat(i)**2 for i in p[1:])) + check(p[0])**2
            tmp = np.exp(- L**2 * c**2 / 4 * psqd)
            prsq = np.sum((ring(i)**2 for i in p[1:]))
            num = 0
            for i in (1,2,3):
                num += prsq * np.cos(p[i] / 2)**2 \
                        - (ring(p[i]) * np.cos(p[i] / 2)**2)**2
            num *= shat(p[0], x)**2
            if psqd:
                result += tmp * num / psqd
    rho = float(T) / L
    return (N**2 - 1) * c**4 * result / 128 / rho

if __name__ == "__main__":
    print tsqEs(4, 4, 2, 0.5)

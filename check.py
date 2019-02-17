#!/usr/bin/env python
import sys

small_primes = [2, 3, 5, 7, 11, 13, 17]

def smooth_part(x):
    res = 1
    factorization = {}
    for p in small_primes:
        while x % p == 0:
            factorization[p] = factorization.get(p, 0) + 1
            res *= p
            x /= p
    return (res, factorization)

if __name__ == '__main__':
    n = int(sys.argv[1])
    print len('{0:b}'.format(n)), smooth_part(n - 1)

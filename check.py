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

def process_one(n):
    print len('{0:b}'.format(n)), smooth_part(n - 1)

def check_line(s):
    cols = s[7:].split(' , ')
    process_one(int(cols[0]))
    process_one(int(cols[1]))


if __name__ == '__main__':
    if len(sys.argv[1]) < 20:
        for line in open(sys.argv[1]):
            check_line(line)
            print
    else:
        n = int(sys.argv[1])
        print len('{0:b}'.format(n)), smooth_part(n - 1)

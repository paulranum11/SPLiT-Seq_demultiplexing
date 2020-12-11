#!/usr/bin/env python


# Hamming script test

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


x = 'ACTGACTG'
green = ['GGCGCATG', 'ACTGAAAT', 'ATGCCCGT', 'ACTGAGTG']

for barcode in green:
	res = hamming(x, barcode)
	print(res)

res2 = [s for s in green if hamming(s, x) <= 1]

print(res2[0])
print(len(res2))
"""
    Contain commonly used helper functions.
"""


import json
import decimal
import sympy


def generatePrime(n):
    """
    Generate a prime number of n bits
    """
    
    p = 2 ** n - 1
    while not sympy.isprime(p):
        p = p - 1

    return p


def norm(v):
    """
    Return the 2-norm of vector
    """
    res = decimal.Decimal(0)
    for vi in v:
        res += decimal.Decimal(vi) * decimal.Decimal(vi)
    return res.sqrt()


def mod(value, q):
    """
    The result is in (-q//2, q//2]
    """
    value = value % q

    if value > q // 2:
        value -= q

    return value


def saveVariable(filename, v):
    with open(filename, 'a') as file:
        file.write(str(v))
        file.write('\n')


def loadVariable(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    last_line = int(lines[-1])
    return last_line


def saveList(lst, filename):
    try:
        with open(filename, 'w') as file:
            json.dump(lst, file)
    except IOError:
        print("An error occurred while saving the list.")


def loadList(filename):
    try:
        with open(filename, 'r') as file:
            lst = json.load(file)
        return lst
    except IOError:
        print("An error occurred while loading the list.")
        return []


def btoi(b):
    return int.from_bytes(b, "big")


def itob(i, len):
    return int.to_bytes(int(i), length=len, byteorder="big")


def computeHNPDataLSB(line, lsb, leakage, q, baselen):
    const1 = pow(2 ** leakage, -1, q)
    h, sig = line.strip().split()
    h = int(h, 16)
    r = int(sig[:2 * baselen], 16)
    s = int(sig[2 * baselen:], 16)
    s_inverse = pow(s, -1, q)
    t = const1 * s_inverse * r % q
    a = const1 * (lsb - s_inverse * h) % q
    return t, a


def computeHNPDataMSB(line, msb, leakage, q, baselen):
    h, sig = line.strip().split()
    h = int(h, 16)
    r = int(sig[:2 * baselen], 16)
    s = int(sig[2 * baselen:], 16)
    s_inverse = pow(s, -1, q)
    t = s_inverse * r % q
    a = (msb * pow(2, (q - 1).bit_length()-leakage) - s_inverse * h) % q
    return t, a
import argparse
import json

q = 115792089210356248762697446949407573529996955224135760342422259061068512044369

def computeHNPDataMSB(line, msb, leakage, q, baselen):
    h, sig = line.strip().split()
    h = int(h, 16)
    r = int(sig[:2 * baselen], 16)
    s = int(sig[2 * baselen:], 16)
    s_inverse = pow(s, -1, q)
    t = s_inverse * r % q
    a = (msb * pow(2, (q - 1).bit_length()-leakage) - s_inverse * h) % q
    return t, a


def loadList(filename):
    try:
        with open(filename, 'r') as file:
            lst = json.load(file)
        return lst
    except IOError:
        print("An error occurred while loading the list.")
        return []


def mod(value, q):
    """
    The result is in (-q//2, q//2]
    """
    value = value % q

    if value > q // 2:
        value -= q

    return value


def boundConstraint(t0_inv, ti, const):
        """
        The bound constraint.
        """
        return abs(mod(t0_inv * ti, q)) <= const


# def outputFormat():
    
    

parser = argparse.ArgumentParser()

parser.add_argument('-f1', type=str, help='The file of lines')
parser.add_argument('-f2', type=str, help='The file of msb')
parser.add_argument('-f3', type=str, help='The file of sk')
parser.add_argument('-s', type=int, help='The leakage')


args = parser.parse_args()
lines = loadList(args.f1)
msb = loadList(args.f2)
sk = loadList(args.f3)

x = 4

const1 = q // 2 ** (args.s + x + 4)
const2 = q // 2 ** (args.s + x - 2)

k_replaced = []
new_vec = []


num1 = num2 = num3 = num4 = 0

t0_inv = 0
for i in range(len(msb)):
    ti, ai = computeHNPDataMSB(lines[i], msb[i], args.s, q, 32)
    if i == 0:
        t0_inv = pow(ti, -1, q)
        continue
    if boundConstraint(t0_inv, ti, const1):
        num1 += 1
    elif boundConstraint(t0_inv, ti, 2 * const1):
        num2 += 1
    elif boundConstraint(t0_inv, ti, const2):
        num3 += 1
    else:
        num4 += 1

print(num1, num2, num3, num4)

    
    

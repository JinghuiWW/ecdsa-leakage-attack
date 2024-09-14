"""
Implement the ECDSA class, which samples ECDSA instances with nonce leakage, specifically where the leakage is restricted to the least significant bits of the nonce.
The ECDSASolver class accepts an ECDSA instance and performs key recovery.
"""
 
import ecdsa
import random
from util import *
from HNPSolver import HNP, HNPSolver
from multiprocessing import Pool


class ECDSA:
    def __init__(self, bit_size, curve=None):
        """
        Initialize the curve according to the bit size of modulus
        """
        curves = {str.lower(str(c)): c for c in ecdsa.curves.curves}
        if curve is None:
            if bit_size == 112:
                curve = "secp112r1"
            elif bit_size == 128:
                curve = "secp128r1"
            elif bit_size == 160:
                curve = "secp160r1"
            elif bit_size == 192:
                curve = "nist192p"
            elif bit_size == 224:
                curve = "nist224p"
            elif bit_size == 256:
                curve = "nist256p"
            elif bit_size == 320:
                curve = "brainpoolp320r1"
            elif bit_size == 384:
                curve = "nist384p"
            elif bit_size == 512:
                curve = "brainpoolp512r1"
            elif bit_size == 521:
                curve = "nist521p"
            else:
                raise NotImplementedError("bit size={bit_size} is not implemented".format(bit_size=bit_size))

        if str.lower(curve) in curves.keys():
            self.curve = curves[curve]
            self.baselen = self.curve.baselen
            self.bit_size = self.curve.order.bit_length()
            self.G = self.curve.generator
            self.q = self.G.order()
        else:
            raise NotImplementedError("curve={curve} is not implemented".format(curve=curve))

    def sign(self, h, sk, leakage):
        """
        Assume that the least significant bits of nonce is leaked
        """
        d = btoi(sk.to_string())
        hi = btoi(h)
        k = random.randint(0, self.q - 1)

        r = (self.G * k).x() % self.q
        s = pow(k, -1, self.q) * ((hi + d * r) % self.q) % self.q
        sig = itob(r, self.baselen) + itob(s, self.baselen)

        return sig, k % (2 ** leakage)

    def sampling(self, sk, leakage):
        """
        In a side-channel attack against ECDSA, the adversary may obtain
        some of the least significant bits of the signature nonce k.
        """
        h = random.randint(0, 2 ** self.bit_size)
        hb = itob(h, self.baselen)
        sig, k1 = self.sign(hb, sk, leakage)
        line = "%s %s" % (bytes.hex(hb), bytes.hex(sig))
        return line, k1

    def boundConstraint(self, t0_inverse, ti, const):
        """
        The bound constraint.
        """
        return abs(mod(t0_inverse * ti, self.q)) <= const

    def parallelSampling(self, args):
        sk, leakage, x, t0_inverse, const1, const2, error_rate, goal = args

        while True:
            line, k1 = self.sampling(sk, leakage)

            t, a = computeHNPDataLSB(line, k1, leakage, self.q, self.baselen)

            if x != 0:
                if goal == "construct":
                    const = const1
                elif goal == "reduce":
                    const = const1 * 2
                elif goal == "narrow":
                    const = const2
                    if (t0_inverse * t % self.q) <= const:
                        break
                    else:
                        continue
                elif goal == "check":
                    const = self.q // 2
                else:
                    const = self.q // 2
                    print("Error type! Sample randomly")

                if self.boundConstraint(t0_inverse, t, const):
                    break
            else:
                break
        
        if random.random() <= error_rate:
            k1 = (random.randint(0, 2 ** leakage - 1)) % (2 ** leakage)

        return [[line, k1]]

    def generate(self, num_samples, leakage, x=0, error_rate=0):
        def parallelGen(num, args, goal):
            if num == 0:
                return [], []

            args = args + (goal,)

            lines = []
            lsb = []
            pool = Pool(num)
            results = []
            for _ in range(num):
                results.append(pool.apply_async(self.parallelSampling, args=([args])))
            pool.close()
            pool.join()
            for res in results:
                for pair in res.get():
                    lines.append(pair[0])
                    lsb.append(pair[1])

            return lines, lsb

        const1 = self.q // 2 ** (leakage + x + 4)
        const2 = self.q // 2 ** (leakage + x - 2)

        sk = ecdsa.SigningKey.generate(curve=self.curve)

        lines = []
        lsb = []

        line, k1 = self.sampling(sk, leakage)
        lines.append(line)
        lsb.append(k1)

        t0, a0 = computeHNPDataLSB(line, k1, leakage, self.q, self.baselen)
        t0_inverse = pow(t0, -1, self.q)

        args = sk, leakage, x, t0_inverse, const1, const2, error_rate

        lines_lattice, lsb_lattice = parallelGen(num_samples[0] - 1, args, goal="construct") 
        lines_reduce, lsb_reduce = parallelGen(num_samples[1], args, goal="reduce")
        lines_narrow, lsb_narrow = parallelGen(num_samples[2], args, goal="narrow")
        lines_check, lsb_check = parallelGen(num_samples[3], args, goal="check")

        lines += lines_lattice + lines_reduce + lines_narrow + lines_check
        lsb += lsb_lattice + lsb_reduce + lsb_narrow + lsb_check

        return lines, lsb, btoi(sk.to_string())


class ECDSASolver:
    """
    Recover the signing key of ECDSA from leaked signatures
    """

    def __init__(self, ECDSA_instance, lines, leakedbits, leakage, num_samples, x, error_rate, sieve_algorithm = "gpu", threads = 4, gpus = 1, model=1):
        self.ECDSA_instance = ECDSA_instance
        self.leakage = leakage
        self.num_samples = num_samples
        self.x = x 
        self.error_rate = error_rate
        self.sieve_algorithm = sieve_algorithm
        self.threads = threads
        self.gpus = gpus
        self.a = []
        self.t = []

        for i in range(len(lines)):
            if model == 1:
                ti, ai = computeHNPDataLSB(lines[i], leakedbits[i], leakage, self.ECDSA_instance.q, self.ECDSA_instance.baselen)
            else:
                ti, ai = computeHNPDataMSB(lines[i], leakedbits[i], leakage, self.ECDSA_instance.q, self.ECDSA_instance.baselen)

            self.t.append(ti)
            self.a.append(ai)

        self.k1 = leakedbits[0]
        h, sig = lines[0].strip().split()
        self.r0 = int(sig[:2 * self.ECDSA_instance.baselen], 16)

    def recoverKey(self):
        """
        Treat the ECDSA instance as a HNP instance, and recover the hidden number.
        Finally recover the signing key
        """

        # Treat the ECDSA instance as a HNP instance
        HNP_instance = HNP(self.ECDSA_instance.q, self.leakage, self.t, self.a)

        # Solve HNP. Recover the hidden number if the algorithm works
        Solver = HNPSolver(HNP_instance, self.num_samples, self.r0, self.k1, self.ECDSA_instance.G, x=self.x, error_rate=self.error_rate, quickCheckTech=True, threads=self.threads, sieve_algorithm=self.sieve_algorithm, gpus=self.gpus)

        # Repeat strategy
        round = 1
        if self.x == 0:
            for _ in range(round):
                if Solver.alpha == 0:
                    Solver.solve("eliminateAlpha")

        else:
            for _ in range(round):
                if Solver.alpha == 0:
                    Solver.solve("increaseVolume")

        return Solver.alpha

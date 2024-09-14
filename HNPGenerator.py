"""
    Implement the HNPGenerator class, which generates a HNP instance.
"""


from random import randint
from util import *
from multiprocessing import Pool


class HNPGenerator:
    """
    Generate a HNP instance.
    """

    def __init__(self, q, s, num_lattice, x=0, num_reduce=0, num_narrow=0, num_check=0):
        self.q = q  # The modulus q
        self.n = (q - 1).bit_length()  # The bit-size
        self.s = s  # The number of leaked bits
        self.l = self.n - self.s
        self.x = x

        # Store some constants to avoid repetitive calculations
        self.const1 = q // 2 ** (s + 1)
        self.const2 = q // 2 ** (s + x + 4)
        self.const3 = q // 2 ** (s + x - 2)
        self.const = [self.const2, self.const2 * 2, self.const3]

        self.num_lattice = num_lattice
        self.num_reduce = num_reduce
        self.num_narrow = num_narrow
        self.num_check = num_check

        self.t = []
        self.a = []
        self.alpha = None
        self.t0_inverse = None

    def sampling(self):
        """
        Sampling randomly
        """
        t = randint(1, self.q - 1)
        a = ((t * self.alpha % self.q) >> self.l) << self.l
        return t, a

    def parallelSampling(self, goal):
        """
        Parallel sampling with some constraints
        """
        t_tmp, a_tmp = self.sampling()
        if self.x != 0:
            while not self.boundConstraint(t_tmp, self.const[goal - 1]):
                t_tmp, a_tmp = self.sampling()

        return [[t_tmp, a_tmp]]

    def boundConstraint(self, ti, const):
        """
        The bound constraint
        """
        return abs(mod(self.t0_inverse * ti, self.q)) < const

    def generate(self):
        """
        Generate (t, a)
        """

        def parallelGen(num, goal):
            if num == 0:
                return [], []
                
            t = []
            a = []
            pool = Pool(num)
            results = []
            for _ in range(num):
                results.append(pool.apply_async(self.parallelSampling, args=([goal])))
            pool.close()
            pool.join()
            for res in results:
                for pair in res.get():
                    t.append(pair[0])
                    a.append(pair[1])
            return t, a

        self.alpha = randint(0, self.q - 1)

        t_tmp, a_tmp = self.sampling()
        self.t.append(t_tmp)
        self.a.append(a_tmp)

        self.t0_inverse = pow(self.t[0], -1, self.q)

        # For the lattice construction
        t_lattice, a_lattice = parallelGen(self.num_lattice - 1, 1)

        # For the prescreening technique
        t_reduce, a_reduce = parallelGen(self.num_reduce, 2)

        # For the interval reduction algorithm
        t_narrow, a_narrow = parallelGen(self.num_narrow, 3)

        # For the linear predicate
        t_check = []
        a_check = []
        for _ in range(self.num_check):
            t_tmp, a_tmp = self.sampling()
            t_check.append(t_tmp)
            a_check.append(a_tmp)

        self.t += t_lattice + t_reduce + t_narrow + t_check
        self.a += a_lattice + a_reduce + a_narrow + a_check

    def checkHiddenNumber(self, recovered_alpha):
        """
        Check whether the recovered hidden number is right
        """
        return self.alpha == recovered_alpha
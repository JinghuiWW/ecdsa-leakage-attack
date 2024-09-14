"""
    The primary module for solving the Hidden Number Problem. 
    Implement the HNPSolver class, which includes a new lattice construction, linear predicate, interval reduction algorithm, and pre-screening technique.
"""


import math
from decimal import Decimal
from time import time
from multiprocessing import Pool
from util import *
from fpylll import IntegerMatrix, BKZ
from fpylll.util import gaussian_heuristic


class HNP:
    def __init__(self, q, s, t, a):
        self.q = q  # The modulus q
        self.n = (q - 1).bit_length()  # The bit-size
        self.s = s  # The number of leakage
        self.l = self.n - self.s  # l is a parameter to HNP
        self.t = t
        self.a = a


class HNPSolver:
    def handleData(self, num_samples):
        len1, len2, len3, _ = num_samples

        for i in range(sum(num_samples)):
            t_mod = self.HNP_instance.t[i] % self.HNP_instance.q
            a_mod = self.HNP_instance.a[i] % self.HNP_instance.q

            # Distribute samples into respective lists
            if i < len1:
                self.t_lattice.append(t_mod)
                self.a_lattice.append(a_mod)
            elif len1 <= i < len1 + len2:
                self.t_reduce.append(t_mod)
                self.a_reduce.append(a_mod)
            elif len1 + len2 <= i < len1 + len2 + len3:
                self.t_narrow.append(t_mod)
                self.a_narrow.append(a_mod)
            else:
                self.t_check.append(t_mod)
                self.a_check.append(a_mod)

    def __init__(self, HNP_instance, num_samples, r=None, k1=None, G=None, x=0, error_rate=0, quickCheckTech=True, threads=4, sieve_algorithm="bdgl", gpus=0):
        """
        Initialize the instance for hidden number recovery.

        Args:
            HNP_instance: An instance containing parameters for the hidden number problem.
            num_samples (tuple): A tuple representing the number of samples for different stages 
                                (num_lattice, num_reduce, num_narrow, num_check).
            r (optional): An optional parameter representing a specific known value (default: None).
            k1 (optional): Another optional parameter (default: None).
            G (optional): The generator point (default: None).
            x (int, optional): An optional parameter affecting lattice construction (default: 0).
            error_rate (float, optional): The error rate for calculations (default: 0).
            quickCheckTech (bool, optional): Enables quick check technique for optimization (default: True).
            threads (int, optional): Number of threads to use for sieving (default: 4).
            sieve_algorithm (str, optional): The sieving algorithm to use (default: "gpu").
            gpus (int, optional): The number of GPUs to use for sieving (default: 1).
        """

        self.HNP_instance = HNP_instance

        # Unpack number of samples for different stages
        self.num_lattice, self.num_reduce, self.num_narrow, self.num_check = num_samples
        # For constructing the lattice
        self.t_lattice, self.a_lattice = [], []
        # For pre-sceening technique
        self.t_reduce, self.a_reduce = [], []
        # For interval reduction algorithm
        self.t_narrow, self.a_narrow = [], []
        # For linear predicate
        self.t_check, self.a_check = [], []

        self.handleData(num_samples)

        # Precompute constants to avoid repetitive calculations
        self.const1 = self.HNP_instance.q // 2 ** (self.HNP_instance.s + 1)
        self.const2 = self.const1 + self.HNP_instance.q // 2 ** (self.HNP_instance.s + 4)

        # Calculate the inverse of the first t_lattice element modulo q
        self.t0_inverse = pow(self.t_lattice[0], -1, self.HNP_instance.q)
        
        # Replace arrays to store adjusted sample values
        self.t_lattice_replaced = [0] * (self.num_lattice - 1)
        self.a_lattice_replaced = [0] * (self.num_lattice - 1)
        self.t_reduce_replaced = [0] * self.num_reduce
        self.a_reduce_replaced = [0] * self.num_reduce
        self.t_narrow_replaced = [0] * self.num_narrow
        self.a_narrow_replaced = [0] * self.num_narrow
        if quickCheckTech:
            for j in range(len(self.t_reduce_replaced)):
                self.t_reduce_replaced[j] = int(self.t0_inverse * self.t_reduce[j] % self.HNP_instance.q)
                self.a_reduce_replaced[j] = int(((self.a_reduce[j] + self.const1) - self.t_reduce_replaced[j] * (
                            self.a_lattice[0] + self.const1)) % self.HNP_instance.q)
            for j in range(len(self.t_narrow_replaced)):
                self.t_narrow_replaced[j] = int(self.t0_inverse * self.t_narrow[j] % self.HNP_instance.q)
                self.a_narrow_replaced[j] = int(self.a_narrow[j] - self.a_lattice[0] * self.t_narrow_replaced[j] % self.HNP_instance.q)

        # Set other parameters
        self.error_rate = error_rate
        self.max_error = [2 * math.ceil(error_rate * num) for num in num_samples]
        self.x = x
        self.threads = threads
        self.sieve_algorithm = sieve_algorithm
        self.gpus = gpus
        self.quickCheckTech = quickCheckTech
        self.embedding_number = int(self.const1 / math.sqrt(3))
        self.lattice = None
        self.gh = None
        self.last_column = None
        self.alpha = 0
        saveList(self.alpha, "hidden_number.txt")
        self.r = r
        self.k1 = k1
        self.G = G

    def constructLatticeEliminateAlpha(self):
        """
        Construct lattice based on [AH21]. The lattice dimension is "num_lattice + 1"
        """
        
        self.lattice = IntegerMatrix(self.num_lattice + 1, self.num_lattice + 1)

        for i in range(self.num_lattice - 1):
            self.lattice[i, i] = self.HNP_instance.q

        for j in range(1, self.num_lattice):
            self.t_lattice_replaced[j - 1] = int(self.t0_inverse * self.t_lattice[j] % self.HNP_instance.q)
            self.lattice[self.num_lattice - 1, j - 1] = self.t_lattice_replaced[j - 1]
            self.a_lattice_replaced[j - 1] = int(((self.a_lattice[j] + self.const1) - self.t_lattice_replaced[j - 1] * (self.a_lattice[0] + self.const1)) % self.HNP_instance.q)
            self.lattice[self.num_lattice, j - 1] = self.a_lattice_replaced[j - 1]

        self.lattice[self.num_lattice - 1, self.num_lattice - 1] = 1
        self.lattice[self.num_lattice, self.num_lattice] = self.embedding_number

    def constructLatticeIncreaseVolume(self):
        """
        Construct lattice with an increased volume. The lattice dimension is "num_lattice + 1"
        """
        self.lattice = IntegerMatrix(self.num_lattice + 1, self.num_lattice + 1)

        for i in range(self.num_lattice - 1):
            self.lattice[i, i] = self.HNP_instance.q

        for j in range(1, self.num_lattice):
            self.t_lattice_replaced[j - 1] = mod(self.t0_inverse * self.t_lattice[j], self.HNP_instance.q)
            self.lattice[self.num_lattice - 1, j - 1] = self.t_lattice_replaced[j - 1] << self.x

            self.a_lattice_replaced[j - 1] = int(((self.a_lattice[j] + self.const1) - self.t_lattice_replaced[j - 1] * (
                        self.a_lattice[0] + self.const1)) % self.HNP_instance.q)
            self.lattice[self.num_lattice, j - 1] = self.a_lattice_replaced[j - 1]

        self.lattice[self.num_lattice - 1, self.num_lattice - 1] = 1 << self.x
        self.lattice[self.num_lattice, self.num_lattice] = self.embedding_number

    def checkLattice(self):
        """
        Check if the lattice contains the hidden number by using an improved linear predicate.
        Iterates through the lattice rows to find a valid solution.
        """
        
        for i in range(self.lattice.nrows):
            target = self.lattice[i, self.lattice.nrows - 2]
            tau = self.lattice[i, self.lattice.nrows - 1]
            if self.improvedLinearPredicate(target, tau):
                self.alpha = int(loadVariable("hidden_number.txt"))
                return True

        return False

    def solve(self, lattice_construction_method):
        """
        Recover the hidden number using lattice-based algorithms
        """

        # Construct the lattice based on the provided method
        if lattice_construction_method == "eliminateAlpha":
            self.constructLatticeEliminateAlpha()
        
        if lattice_construction_method == "increaseVolume":
            self.constructLatticeIncreaseVolume()

        # BKZ pre-processing
        BKZ.reduction(self.lattice, BKZ.Param(20))
        
        # Check if the lattice contains the hidden number
        if self.checkLattice():
            print("Find the target vector in the lattice generated by BKZ")
            return

        # Initialize sieving with G6K and configure the siever parameters    
        if self.gpus > 0:
            from g6k.siever import Siever, SieverParams
            g6k = Siever(self.lattice, SieverParams(default_sieve=self.sieve_algorithm, threads=self.threads, gpus=self.gpus))
        else:
            from g6k import Siever, SieverParams
            g6k = Siever(self.lattice, SieverParams(default_sieve=self.sieve_algorithm, threads=self.threads))

        # Calculate Gaussian heuristic for the lattice
        self.gh = Decimal(gaussian_heuristic([g6k.M.get_r(i, i) for i in range(self.lattice.nrows)]))
        self.gh = self.gh.sqrt()

        # Calculate the expected length of the target vector
        expected_length = Decimal(self.num_lattice + 1) * Decimal(self.const1 ** 2) / Decimal(3) * Decimal(1 + self.error_rate * (pow(2, 2 * self.HNP_instance.s) - 1))
        
        # Calculate the ratio between the expected length and the Gaussian heuristic
        ratio = expected_length.sqrt() / Decimal(self.gh)

        if ratio <= 1.1547:
            print("The ratio between the expected length of target vector and gaussian heuristic: {:.4f} <= sqrt(4/3) = 1.1547".format(ratio))
        else:
            print("The ratio between the expected length of target vector and gaussian heuristic: {:.4f} > sqrt(4/3) = 1.1547".format(ratio))
        
        g6k.resize_db(0)
        g6k.initialize_local(0, self.lattice.nrows // 2, self.lattice.nrows)

        while g6k.l > 0:
            t1 = time()
            g6k.extend_left(1)
            try:
                g6k()
            except Exception as e:
                print(f"Error during sieving: {e}")
            t2 = time()

            # Output the sieving process information
            if g6k.l <= 10:
                print("g6k.l", g6k.l, "g6k.db_size", g6k.db_size(), "Sieving Dim", self.lattice.nrows - g6k.l, "Spent Time {:.1f} s".format(t2 - t1), flush=True)

            # Once sieving is complete, check the database for the hidden number
            if g6k.l == 0:
                # Store the last two columns of lattice
                self.last_column = self.lattice.submatrix(range(self.lattice.nrows), [self.lattice.ncols - 2, self.lattice.ncols - 1])
                for i in range(self.lattice.nrows):
                    self.last_column[i, 1] = self.last_column[i, 1] // self.embedding_number
                
        self.checkDatabaseBatch(g6k.itervalues(), int(g6k.db_size()))

    def checkDatabase(self, db):
        """
        Check all the vectors in the database to see if they contain information about the hidden number
        """
        
        search_lst = [([v]) for v in db]

        if self.threads == 1:
            flag = False
            t1 = time()
            
            for args in search_lst:
                if self.x == 0:
                    if self.checkEliminateAlpha(args):
                        flag = True
                        break
                else:
                    if self.checkIncreaseVolume(args):
                        flag = True
                        break
            
            t2 = time()
            print("Time of searching the database with one thread:{:.2f} s".format(t2 - t1))
            
            return flag
        
        # Check the database in parallel
        else:
            t1 = time()
            
            pool = Pool(processes=self.threads)

            if self.x == 0:
                pool.map(self.checkEliminateAlpha, search_lst)
            else:
                pool.map(self.checkIncreaseVolume, search_lst)

            pool.close()
            pool.join()
            
            t2 = time()
            print("Time of searching the database in parallel: {:.2f} s".format(t2 - t1))
            
            return (int(loadList("hidden_number.txt")) > 0)

    def checkDatabaseBatch(self, whole_db, db_size):
        """
        Check the vectors in the database in batches to prevent program crashes due to memory issues
        """
        t1 = time()
        batch_size = int(5e5)

        # Check the vectors in the database in batches
        if db_size <= batch_size:
            if self.checkDatabase(list(whole_db)):
                self.alpha = int(loadList("hidden_number.txt"))
        else:
            check_round = math.ceil(db_size / batch_size)
            print("batch_size:", batch_size, "Round of checkDatabase:", check_round)
            db = []
            i = 0
            for item in whole_db:
                db.append(item)
                if len(db) == batch_size:
                    vlen_first = norm(self.lattice.multiply_left(db[0])) / self.gh
                    vlen_last = norm(self.lattice.multiply_left(db[-1])) / self.gh
                    print("Searching for vectors with lengths between {:.4f}".format(vlen_first), "and {:.4f}".format(vlen_last))

                    if self.checkDatabase(db):
                        self.alpha = int(loadList("hidden_number.txt"))
                        db = []
                        break

                    db = []
                    i += 1
                    print("Finished, process of checking the database ", i, "/", check_round)

            if len(db) != 0:
                vlen_first = norm(self.lattice.multiply_left(db[0])) / self.gh
                vlen_last = norm(self.lattice.multiply_left(db[-1])) / self.gh
                print("Searching for vectors with lengths between {:.4f}".format(vlen_first), "and {:.4f}".format(vlen_last))
                self.checkDatabase(db)
                print("The database checking has finished")

        t2 = time()
        print("checkDatabase Time:{:.2f} min".format((t2 - t1) / 60))

    def checkEliminateAlpha(self, v):
        """
        Check whether the vector contains information about the hidden number using an improved linear predicate.
        """
        
        # Only require two vector-vector multiplications
        v = self.last_column.multiply_left(v[0])
        
        if abs(v[1]) != 1:
            return False

        # if self.nonlinearPredicate(v[0], v[1] * self.embedding_number):
        if self.improvedLinearPredicate(v[0], v[1] * self.embedding_number):
            return True
        
        return False

    def checkIncreaseVolume(self, v): 
        """
        Check function for our new construction of lattice
        """
        
        # Only require two vector-vector multiplications
        v = self.last_column.multiply_left(v[0])
        
        if abs(v[1]) != 1:
            return False

        if self.linearPredicateIncreaseVolume(v[0], v[1] * self.embedding_number):
            return True

        return False

    def improvedLinearPredicate(self, target, tau):
        """
        Improved linear predicate.
        """
        
        # Check for early exit conditions
        if target == 0 or abs(target) > self.const1:
            return False

        # Adjust k0 based on the value of tau
        if tau == self.embedding_number:
            k0 = self.const1 - target
        elif tau == - self.embedding_number:
            k0 = target + self.const1
        else:
            return False

        # Compute the potential hidden number
        alpha = self.t0_inverse * int(self.a_lattice[0] + k0) % self.HNP_instance.q

        # Use additional samples to ensure the uniqueness of the solution
        error_count = 0
        for i in range(len(self.t_check)):
            if int(alpha * self.t_check[i] - self.a_check[i]) % self.HNP_instance.q > 2 * self.const1:
                error_count += 1
                # Early exit if errors exceed the allowed threshold
                if error_count > self.max_error[3]:
                    return False

        print("The hidden number is found:", alpha)
        saveList(alpha, "hidden_number.txt")

        return True # self.nonlinearPredicate(target, tau)

    def linearPredicateIncreaseVolume(self, target, tau):
        """
        Improved linear predicate for our new lattice construction
        """

        if target == 0 or abs(tau) != self.embedding_number:
            return False

        if tau == self.embedding_number:
            tau, target = -tau, -target

        # Use quick check technique if enabled, otherwise perform enumeration
        if self.quickCheckTech:
            return self.quickCheck(target, tau) 
        else:
            return self.enumerateCheck(target, tau)

    def enumerateCheck(self, target, tau):
        """
        Enumerate integers in the range [-2^{x-1}, 2^{x-1}) and check the predicate for the range [target - 2^{x-1}, target + 2^{x-1}).
        """
        for num in range(target - pow(2, self.x - 1), target + pow(2, self.x - 1)):
            if self.improvedLinearPredicate(num, tau):
                return True

        return False

    def presceening(self, target):
        error_count = 0
        for i in range(self.num_reduce):
            tmp = (self.t_reduce_replaced[i] * target - self.a_reduce_replaced[i] + self.const2) % self.HNP_instance.q - self.const2
            if abs(tmp) > self.const2:
                error_count += 1
                if error_count > self.max_error[1]:
                    return False

        return True

    def quickCheck(self, target, tau):
        """
            Quick check: first reduce the number of possible target by t_reduce_replaced and a_reduce_replaced
            then narrow the possible range of each target by t_narrow_replaced and a_narrow_replaced
        """

        def intersection(intervals_a, intervals_b):
            if intervals_a == [] or intervals_b == []:
                return []

            res = []
            i, j = 0, 0
            while i < len(intervals_a) and j < len(intervals_b):
                a_start, a_end = intervals_a[i]
                b_start, b_end = intervals_b[j]
                if a_end < b_start:
                    i += 1
                elif b_end < a_start:
                    j += 1
                elif a_end >= b_end:
                    res.append([max(a_start, b_start), b_end])
                    j += 1
                else:
                    res.append([max(a_start, b_start), a_end])
                    i += 1
            return res

        def remove(R):
            if R == []:
                return []
            R_new = []
            i = 0
            x, y = R[0][0], R[0][1]
            while i < len(R) - 1:
                i += 1
                if y == R[i][0] or y == R[i][0] - 1:
                    y = R[i][1]
                else:
                    R_new.append([x, y])
                    x, y = R[i][0], R[i][1]
            R_new.append([x, y])
            return R_new

        def divideInterval(t, r):
            interval = []
            for ri in r:
                tmp1 = ri[1] // t
                tmp2 = (ri[0] - 1) // t + 1
                if tmp1 >= tmp2:
                    interval.append([tmp2, tmp1])
            if interval != []:
                return remove(interval)

            return []

        def intervals(a, t, q, l1, l2, alpha1, alpha2):
            n1, n2 = (t * alpha1 - a - l2 - 1) // q + 1, (t * alpha2 - a - l1) // q

            if n1 > n2:
                return []

            if n1 == n2:
                return divideInterval(t, [[max(a + n1 * q + l1, t * alpha1), min(a + n2 * q + l2, t * alpha2)]])

            r = []
            r.append([int(max(a + n1 * q + l1, t * alpha1)), int(a + n1 * q + l2)])
            for i in range(n1 + 1, n2):
                r.append([int(i * q + a + l1), int(i * q + a + l2)])
            r.append([int(a + n2 * q + l1), int(min(t * alpha2, a + n2 * q + l2))])

            return divideInterval(t, r)

        def narrowRange(low, high, arguments):
            t, a, q, const1 = arguments
            r = [[low, high]]
            for i in range(len(t)):
                if t[i] > q // 2:
                    continue
                r_new = intervals(a[i], t[i], q, 0, 2 * const1 - 1, low, high)
                r = intersection(r, r_new)
            r = remove(r)
            return r

        # Pre-sceening technique
        if not self.presceening(target):
            return False

        # Subsampling technique
        repeat = 1
        if self.error_rate != 0:
            repeat = 10
        for _ in range(repeat):
            indices = random.sample(list(range(self.num_narrow)), self.x)
            perfect_t_narrow_replace = [self.t_narrow_replaced[i] for i in indices]
            perfect_a_narrow_replaced = [self.a_narrow_replaced[i] for i in indices]
            argus = perfect_t_narrow_replace, perfect_a_narrow_replaced, self.HNP_instance.q, self.const1
            new_range = narrowRange(target + self.const1 - pow(2, self.x + 1), target + self.const1 + pow(2, self.x + 1), argus)
            
            for r in new_range:
                for k0_replaced in range(r[0] - self.const1, r[1] + 1 - self.const1):
                    if self.improvedLinearPredicate(k0_replaced, tau):
                        return True

        return False

    def nonlinearPredicate(self, target, tau):
        """
        Nonlinear predicate used in [AH21]. This function checks whether `r` is equal to ([k]G)_x.
        """
        
        # Check for early exit conditions
        if target == 0 or abs(target) > self.const1:
            return False

        # Adjust k0 based on the value of tau
        if tau == self.embedding_number:
            k0 = self.const1 - target
        elif tau == -self.embedding_number:
            k0 = target + self.const1
        else:
            return False

        # Calculate k
        k = k0 * pow(2, self.HNP_instance.s) + self.k1

        # Check if the result matches the given condition
        if self.r == (self.G * k).x() % self.HNP_instance.q:
            alpha = self.t0_inverse * int(self.a_lattice[0] + k0) % self.HNP_instance.q
            print("The hidden number is found:", alpha)
            saveList(alpha, "hidden_number.txt")
            return True

        return False
"""
    Solve a HNP instance.
"""


from HNPGenerator import HNPGenerator
from HNPSolver import HNP, HNPSolver
from util import *
import argparse
from time import time, process_time

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Solving Hidden Number Problem.")
    parser.add_argument('-n', type=int, help='Bit size of the modulus')
    parser.add_argument('-s', type=int, help='Leakage (LSB)')
    parser.add_argument('-m', type=int, help='Number of equations for constructing the lattice')
    parser.add_argument('-x', type=int, help='Parameter for increasing volume in the lattice')
    parser.add_argument('-g', type=int, help='Number of GPUs to use')
    parser.add_argument('-t', type=int, help='Number of CPU threads to use')

    args = parser.parse_args()

    # Generate a prime number
    q = generatePrime(args.n)
    
    num_samples = args.m, 2 * args.x, args.x, 2 * args.n

    # Create HNP instance by HNPGenerator and recover the hidden number by HNPSolver
    HNPGen = HNPGenerator(q, args.s, args.m, x=args.x, num_reduce=num_samples[1], num_narrow=num_samples[2], num_check=num_samples[3])
    t1 = time()
    HNPGen.generate()
    t2 = time()
    print("Sampling Time:{:.2f} min".format((t2 - t1) / 60))

    HNP_instance = HNP(q, args.s, HNPGen.t, HNPGen.a)
    
    if args.g > 0:
        sieve_alg = "gpu"
    else:
        sieve_alg = "bdgl"
        
    Solver = HNPSolver(HNP_instance, num_samples, x=args.x, quickCheckTech=True, threads=args.t, sieve_algorithm=sieve_alg, gpus=args.g)

    walltime_1 = time()
    cputime_1 = process_time()

    if args.x == 0:
        Solver.solve("eliminateAlpha")
    else:
        Solver.solve("increaseVolume")

    cputime_2 = process_time()
    walltime_2 = time()
    saveVariable("walltime.txt", walltime_2 - walltime_1)
    saveVariable("cputime.txt", cputime_2 - cputime_1)

    print("HNPSolver Spent Time:{:.2f} min".format((walltime_2 - walltime_1) / 60))

    # Check the result
    if HNPGen.checkHiddenNumber(Solver.alpha):
        print("Success")
    else:
        print("Failed to recover the hidden number")
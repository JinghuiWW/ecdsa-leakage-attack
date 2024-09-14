"""
    Attack targeting ECDSA nonce leakage from the most significant bits.
"""

from ECDSASolver import ECDSA, ECDSASolver
from util import *
import argparse
from time import time, process_time

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ECDSA Nonce Leakage Lattice-based Attack.")
    parser.add_argument('-curve', type=str, help='Elliptic curve', default=None)
    parser.add_argument('-n', type=int, help='Bit size of the modulus', required=True)
    parser.add_argument('-s', type=int, help='Leakage (MSB)', required=True)
    parser.add_argument('-m1', type=int, help='Number of equations for constructing the lattice', required=True)
    parser.add_argument('-m2', type=int, help='Number of equations for Pre-screening', default=0)
    parser.add_argument('-m3', type=int, help='Number of equations for interval reduction algorithm', default=0)
    parser.add_argument('-m4', type=int, help='Number of equations for linear predicate', required=True)
    parser.add_argument('-x', type=int, help='Parameter for increasing volume in the lattice', default=0)
    parser.add_argument('-g', type=int, help='Number of GPUs to use', default=0)
    parser.add_argument('-t', type=int, help='Number of CPU threads to use', default=4)
    parser.add_argument('-e', type=float, help='Error rate', default=0.0)
    parser.add_argument('-f1', type=str, help='File containing data and public key', required=True)
    parser.add_argument('-f2', type=str, help='File containing MSB leakage data', required=True)
    parser.add_argument('-f3', type=str, help='File containing the secret key', default=None)

    args = parser.parse_args()

    # Create an ECDSA instance
    ECDSA_instance = ECDSA(args.n, args.curve)

    lines = loadList(args.f1)
    msb = loadList(args.f2)
    
    # Solve ECDSA instance using ECDSASolver
    num_samples = args.m1, args.m2, args.m3, args.m4
    
    lines = lines[:sum(num_samples)]
    msb = msb[:sum(num_samples)]

    if args.g > 0:
        sieve_alg = "gpu"
    else:
        sieve_alg = "bdgl"
        
    Solver = ECDSASolver(ECDSA_instance, lines, msb, args.s, num_samples, args.x, args.e, sieve_algorithm=sieve_alg, threads=args.t, gpus=args.g, model=0)

    walltime_1 = time()
    cputime_1 = process_time()

    recovered_sk = Solver.recoverKey()

    cputime_2 = process_time()
    walltime_2 = time()
    print("ECDSASolver Spent Time:{:.2f} min".format((walltime_2 - walltime_1) / 60))

    saveVariable("walltime.txt", walltime_2 - walltime_1)
    saveVariable("cputime.txt", cputime_2 - cputime_1)
    
    # Check the result
    sk_file = args.f3
    if sk_file is not None:
        sk = loadList(args.f3)
        if sk == recovered_sk:
            print("Success")
        else:
            print("Failed to recover the secret key")
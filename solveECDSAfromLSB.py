"""
    Attack targeting ECDSA nonce leakage from the least significant bits.
"""

from ECDSASolver import ECDSA, ECDSASolver
from util import *
import argparse
from time import time, process_time

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="ECDSA Nonce Leakage Lattice-based Attack.")
    parser.add_argument('-curve', type=str, help='Elliptic curve', default=None)
    parser.add_argument('-n', type=int, help='Bit size of the modulus', required=True)
    parser.add_argument('-s', type=int, help='Leakage (LSB)', required=True)
    parser.add_argument('-m', type=int, help='Number of equations for constructing the lattice', required=True)
    parser.add_argument('-x', type=int, help='Parameter for increasing volume in the lattice', default=0)
    parser.add_argument('-g', type=int, help='Number of GPUs to use', default=0)
    parser.add_argument('-t', type=int, help='Number of CPU threads to use', default=4)
    parser.add_argument('-e', type=float, help='Error rate', default=0.0)
    parser.add_argument('-f', type=int, help='Data source: 0 to read from files, 1 to generate randomly and save', default=None)
    parser.add_argument('-f1', type=str, help='File containing data and public key', default=None)
    parser.add_argument('-f2', type=str, help='File containing LSB leakage data', default=None)
    parser.add_argument('-f3', type=str, help='File containing the secret key', default=None)


    args = parser.parse_args()

    # Create an ECDSA instance
    ECDSA_instance = ECDSA(args.n, args.curve)
    
    num_samples = args.m, 2 * args.x, args.x, 2 * args.n
    if args.f1 is None or args.f == 1:
        t1 = time()
        lines, lsb, sk = ECDSA_instance.generate(num_samples, args.s, x=args.x, error_rate=args.e)
        t2 = time()
        print("ECDSA samples generation:{:.2f} min".format((t2 - t1) / 60))

        if args.f == 1:
            saveList(lines, args.f1)
            saveList(lsb, args.f2)
            saveList(sk, args.f3)
    else:
        lines = loadList(args.f1)
        lsb = loadList(args.f2)
        sk = loadList(args.f3)

    # Solve ECDSA instance using ECDSASolver    
    if args.g > 0:
        sieve_alg = "gpu"
    else:
        sieve_alg = "bdgl"
        
    Solver = ECDSASolver(ECDSA_instance, lines, lsb, args.s, num_samples, args.x, args.e, sieve_algorithm=sieve_alg, threads=args.t, gpus=args.g, model=1)

    walltime_1 = time()
    cputime_1 = process_time()

    recovered_sk = Solver.recoverKey()

    cputime_2 = process_time()
    walltime_2 = time()
    print("ECDSASolver Spent Time:{:.2f} min".format((walltime_2 - walltime_1) / 60))

    saveVariable("walltime.txt", walltime_2 - walltime_1)
    saveVariable("cputime.txt", cputime_2 - cputime_1)

    # Check the result
    if sk == recovered_sk:
        print("Success")
    else:
        print("Failed to recover the secret key")
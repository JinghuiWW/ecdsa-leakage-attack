# ECDSA Nonce Leakage Lattice-based Attack

This repository contains Python source code for attacking ECDSA with nonce leakage using lattice-based algorithms. The specific algorithm is detailed in the paper:

** Yiming Gao, Jinghui Wang, Honggang Hu and Binang He, Attacking ECDSA with Nonce Leakage by Lattice Sieving: Bridging the Gap with Fourier Analysis-based Attacks. ** 

In this paper, we aim to give a solution to an open question: Can lattice-based attacks be enhanced by utilizing more samples? Using this repository, we can break 160-bit ECDSA with 1-bit leakage using approximately $2^{25}$ ECDSA samples. In addition, our new algorithms
for solving the HNP are extended to address the case of erroneous input, increasing the robustness of lattice-based attacks.

## Key Recovery of ECDSA with Nonce Leakage
You can perform the attack on an instance using multiple CPU cores. For example:
``` shell
python solveECDSAfromLSB.py -n 256 -s 4 -m 65 -t 4
```
This command will solve an ECDSA (256, 4) instance using 4 CPU threads.
- -n specifies the bit-size of modulus.
- -s indicates the leakage.
- -m is the number of samples for constructing lattice.
- -t defines the number of CPU threads to use.

To accelerate the attack using GPUs, specify the number of GPUs with the -g option. For detailed information on all parameters, run:
``` shell
python solveECDSAfromLSB.py -h
```

Another example, for solving an ECDSA(128, 1) instance, run:
 ``` shell
python solveECDSAfromLSB.py -n 128 -s 1 -m 117 -x 15 -t 24 -g 2 -f 0 -f1 "Instances/128_1/lines.txt" -f2 "Instances/128_1/lsb.txt" -f3 "Instances/128_1/sk.txt"
```

You can also perform attacks on public datasets such as [minerva](https://github.com/crocs-muni/minerva/tree/master/data). Note that the dataset indicates the number of leading-zero bits (i.e., most significant bits) of the nonce. For example:
``` shell
python solveECDSAfromMSB.py -n 256 -s 3 -m1 90 -m4 512 -t 16 -f1 "minerva-data/athena/256_3/lines.txt" -f2 "minerva-data/athena/256_3/msb.txt"
```
This command will solve an ECDSA instance from the Athena dataset, where:

- -n specifies the bit-size of modulus.
- -s indicates the leakage.
- -m1 is the number of samples for constructing lattice.
- -m4 is the number of samples for the linear predicate.
- -t defines the number of CPU threads to use.
- -f1 provides the path to the file containing the lines (including data and public key (r, s)).
- -f2 the path to the file containing the MSBs.

For convience, we provide two bash scripts, autoRunECDSALSB.sh and autoRunECDSAMSB.sh. Depending on the target of the attack, you could modify the parameters in the corresponding script file.

To execute the attack for the LSB situation, use the following command:

``` shell
bash autoRunECDSALSB.sh
```

To execute the attack for the MSB situation, use the following command:

``` shell
bash autoRunECDSAMSB.sh
```

## Environment
To conduct the attacks, ensure the following environment is properly set up on your machine:
- [FPLLL](https://github.com/fplll/fplll) and [FPyLLL](https://github.com/fplll/fpylll) for data structures and BKZ algorithm.
- [G6K](https://github.com/fplll/g6k) for lattice sieving.
- [G6K-GPU-Tensor](https://github.com/WvanWoerden/G6K-GPU-Tensor) for lattice sieving with GPUs accelerating.


## New Records of Lattice-based Attacks against ECDSA
We achieved several new records of lattice-based attacks against ECDSA. For reproducibility, the successfully broken instances have been stored in the Instances folder. These attacks were conducted using an Intel Xeon Platinum 8480+ CPU and four GeForce RTX 4090 GPUs. Specific details about time and memory consumption can be seen in the following table: 

**4-bit leakage**
| **Curve**        | **Leakage** | **d**   | **x**     | **Expected Sample Size** | **Wall time** | **Mem GiB** |
|------------------|-------------|---------|-----------|-------------|---------------|-------------|
| brainpoolp512r1   | 4           | 130     | 0         | $2^{10}$    | 96min         | 254         |

**1-bit leakage**

| **Curve**        | **Leakage** | **d**   | **x**     | **Expected Sample Size** | **Wall time** | **Mem GiB** |
|------------------|-------------|---------|-----------|-------------|---------------|-------------|
| secp128r1        | 1           | 131     | 0         | $2^8$       | 72min         | 294         |
| secp128r1        | 1           | 118     | 15  | $2^{26}$    | 8min          | 53          |
| secp160r1        | 1           | 144     | 14  | $2^{25}$    | 824min        | 1939        |
| secp160r1        | 1           | 138     | 25  | $2^{36}$    | 279min        | 850         |

**less than 1-bit leakage**

| **Curve**        | **Error rate** | **d**   | **x**       | **Expected Sample Size** | **Wall time** | **Mem GiB** |
|------------------|----------------|---------|-------------|-------------|---------------|-------------|
|secp128r1 |  0.1            | 140     | 20    | $2^{31}$    | 370min        | 1090        |
| secp160r1        | 0.02           | 144     | 14    | $2^{25}$    | 1009min       | 1960        |

Note that the definition of $x$ in this implementation differs slightly from that in the paper. Specifically, $2^{x}$ in the implementation corresponds to $x$ in the paper.

## Acknowledgements
This work was supported by National Natural Science Foundation of China (Grant No. 62472397) and Innovation Program for Quantum Science and Technology (Grant No. 2021ZD0302902).

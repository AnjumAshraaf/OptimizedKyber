#  Optimized Kyber Key Encapsulation Mechanism (KEM)

This repository provides an optimized C implementation of the [Kyber](https://pq-crystals.org/kyber/) post-quantum key encapsulation algorithm, focusing on improving the performance of its polynomial multiplication operations.

The optimization fuses three multiplication routines into a single efficient function, significantly reducing computational overhead and improving overall execution speed compared to the reference implementation.

---

##  Project Objectives

- Optimize polynomial multiplication in the Kyber KEM.
- Maintain algorithmic correctness and cryptographic equivalence.
- Achieve performance gains over the reference implementation.
- Ensure compatibility with Linux systems using GCC.

---

---

##  Build Instructions

### Prerequisites

- Linux OS (tested on Ubuntu 18.04+)
- GCC (version 7 or higher)
- `make` utility

### Build Steps

```bash
git clone https://github.com/AnjumAshraaf/OptimizedKyber.git
cd OptimizedKyber
make speed





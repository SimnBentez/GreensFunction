# GreensFunction

GreensFunction is a high-performance C++ computational engine with Python bindings
designed to solve beam problems using Green's functions as an alternative to classical stiffness-based methods.

This project aims to provide a faster and more flexible framework for beam analysis,
allowing arbitrary load configurations and easy extension to different beam theories
such as Euler-Bernoulli, Timoshenko, and Winkler foundation models.

---

## Motivation

Classical beam solvers rely heavily on stiffness matrix assembly and global system inversion.
While effective, these methods can become computationally expensive and less flexible
when dealing with arbitrary distributed loads or alternative constitutive models.

GreensFunction explores a convolution-based approach using fundamental solutions,
reducing matrix assembly overhead and enabling more direct physical modeling.

The primary objective is to significantly reduce computation time
(from minutes to seconds) while maintaining numerical accuracy.

---

## Key Features (Planned)

- Green's function kernels for:
  - Euler-Bernoulli beams
  - Timoshenko beams
  - Winkler foundation models
- High-performance convolution engine in C++
- Python interface via pybind11
- Support for arbitrary load distributions
- Benchmark comparison against classical stiffness methods

---

## Mathematical Background

The method is based on the convolution between the load distribution and
the corresponding Green's function of the governing differential equation:

$$
v_f(x) = \int G(x, \xi) \cdot q(\xi) \text{d}\xi
$$

Where:
- G(x, ξ) is the fundamental solution
- q(ξ) is the applied load
- w(x) is the beam deflection

---

## Project Status

Early development stage.

Current focus:
- Core C++ kernel implementation
- Performance benchmarking
- Python bindings integration

---

## Roadmap

- [ ] Implement 1D Euler-Bernoulli Green kernel
- [ ] Add numerical integration engine
- [ ] Benchmark vs stiffness matrix method
- [ ] Add Timoshenko formulation
- [ ] Publish first PyPI release

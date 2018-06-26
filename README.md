# NSGA-III
An implementation of NSGA-III algorithm in C++ according to the article:

* Kalyanmoy Deb, & Himanshu Jain, An Evolutionary Many-Objective Optimization Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part I: Solving Problems With Box Constraints. IEEE Transactions on Evolutionary Computation, Vol. 18, No. 4, pp. 577–601, Aug 2014. doi:10.1109/TEVC.2013.2281535.

This code was tested using DTLZ and WFG test problems and the obtained results were quite similar to those reported by the authors.

Contributions and bug fixes are welcome.

## Demonstration
The demo "src/main_nsga3.cpp" consists of a scalability test of the algorithm NSGA-III by varying the number of objectives from three to ten considering the DTLZ2 problem.

Compile the specific target file:
```bash
make Makefile
```

# sharpSMT

## What is it?
sharpSMT is a tool designed for measuring the solution space of an SMT(LA) formula,
It is licensed under the [GNU General Public License](COPYING).
Intuitively, sharpSMT subdivides the solution spaces of SMT(LA) formulas into a series of polytopes,
then it calls tools on polytopes and finaly sums the volume or counts up.
In detail, it calls various polytope subroutines as the following.

| Polytope Subroutine | Functions | Command-line Option |
|  ------------- | ------------- | ------------- |
| PolyVest | Volume Approximation | -p |
| Vinci | Volume Computation | -v |
| LattE | Exact Integer Solution Counting | -l |
| Barvinok | Exact Integer Solution Counting | -b |
| ApproxLatCount | Approximate Integer Solution Counting | -a |
| Vol2Lat | Approximate Integer Solution Counting | -v |

This release includes the source codes of PolyVest, Vinci, ApproxLatCount (ALC) and Vol2Lat (V2L).
The tool LattE and Barvinok are not included.


## Directories
| Name           | Description   |
|  ------------- | ------------- |
| bin/ | Polytope Subroutines |
| src/ | Sorce codes of sharpSMT |
| test/ | Formulas for Testing |
| ApproxLatCount.zip | Source codes of ALC |
| benchmarks.zip | The complete set of benchmarks |
| build.sh | Shell for building sharpSMT |
| COPYING | GNU GPL lincese |
| makefile | makefile |
| vinci-1.0.5.zip | Source codes for a modified version of Vinci |
| z3-master.zip | Source codes of Z3 |

## Build status
This release of sharpSMT has been successfully built on the following operating systems:
* Ubuntu 22.04 on 64-bit with g++ 11.4.0

## Building sharpSMT
* Step 1: Make sure that g++ is installed on your machine (you can type "g++ -v" to check this).
* Step 2: The functionality of sharpSMT is dependent on some other libraries: [boost](http://www.boost.org/), [glpk](http://www.gnu.org/software/glpk/), and [Armadillo](http://arma.sourceforge.net/).
* Step 3: Execute:
```bash
sh build.sh
```
* Step 4: Build and install [LattE](https://www.math.ucdavis.edu/~latte/). Then move the executable files (*count* and *scdd*) into directory *bin/*. It is optional if you do not use LattE (-L) to compute the integer solution counts.

* Step 5: Build and install [barvinok](https://barvinok.sourceforge.io/). Then move the executable files (*barvinok_count*) into directory *bin/*. It is optional if you do not use barvinok to compute the integer solution counts.

Quick test, simply execute:
```bash
./sharpSMT -v test/f1_lra.smt2
./sharpSMT -a -w=2 test/coloring.smt
```

sharpSMT should pop up the help menu by "-h" command.
```bash
./sharpSMT -h
```


### Quick guide for building on Ubuntu

Execute:

```bash
sudo apt install g++ libglpk-dev libboost-dev libarmadillo-dev
sh build.sh
```


### Questions/Feedback/Comments ###
Please contact:

  1. Cunjing Ge ([gecj@ios.ac.cn](mailto:gecunjing@nju.edu.cn))


Enjoy!




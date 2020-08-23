# History 

1. Around 1997, Mike Stillman and Harry Tsai created some code to compute SAGBI bases. Two versions of the SAGBI code were written, one just using the front-end of Macaulay2 and the other built into the C/C++ engine of Macaulay2. 
2. Around 2006, Mike Sillman created a package called SubalgebraBases to provide documentation and the front-end package for computing SAGBI bases. Sometime before 2018, a change to the way monomial orders are represented broke at least one of the two implementations of SAGBI bases. 
3. In the Summer of 2019 (at the David Cox retirement conference), Michael Burr, Tim Duff, and Mike Stillman fixed the error in the engine due to the change in monomial orders for computing SAGBI bases. 
4. In the Summer of 2020 (at the Cleveland M2 workshop), Michael Burr, Tim Duff, and Elise Walker worked on modernizing and cleaning the code. Some future projects were also planned. 
5. Summer of 2020 (Warwick M2 workshop) - Work was done on refactoring the code and implementing a few new algorithms.


# Where's the Code?

1. [Version on the master branch](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases.m2)
2. However, the code was still a bit opaque. A stripped-down, working rewrite is available on this repository. Try [installing](./install.m2) it!

# What are (Khovanskii bases / SAGBI bases / Canonical Subalgebra Bases)?

1. See [here](./Subalgebra-Basics/StuCh11.pdf) for a first answer to this question.

## 

1. It is **not** reccommended to work on the documenation because the package is still changing. 
2. The package's "official" tests (which are automatically run when the package installs) are broken and should be ignored for now. Instead, the package's working tests can be run by executing the file benchmark.m2. This is done for several reasons:
    - If the official tests fail when you are trying to debug, it would normally prevent the package from installing and therefore prevent debugging.
    - The package's official tests have to run in a very short amount of time because they contribute to the total time it takes to build Macaulay2. Putting the tests in the file benchmark.m2 allows us to have a more comprehensive set of tests that we can use during development. 
    - Not everyone is an Emacs wizard
    

## Design notes

1. The Subring type:
    - The function "subring" is the canonical constructor of the Subring type.
    - An instance of Subring should be associated with a particular set of generators. Any operation that modifies the generators should return a new Subring instance.
    - The function "debugPrintAllMaps" is the easiest way to understand the data contained inside of an instance of Subring.
    - Every subring instance has a cache. Only use the cache when it can be done in a way that does not cause side effects. The state of the cache should never effect the result of a function. Also, all of the data that is considered the result of a function should be returned instead of stored in the cache.
    
    
 

## Recommended tasks for new contributors

Using the package is a good way to learn about the code as well as the underlying theory. Here are some suggestions for new contributors:

1. Write code that tests the package's newer features:
    - Partial Sagbi bases (Especially examples of using them when a finite Sagbi basis doesn't exist)
    - Modules over subrings (See the file Khovanskii-Group/subring_modules.m2) 
    - Groebner bases of ideals inside of subrings (Usage of the function extrinsicBuchberger)
    - Monomial syzygies (The function toricSyz)
    
2. Find specific computations that involve features that aren't implemented yet:
    - non-monomial syzygies, resolutions
    - Khovanskii bases
    - Newton–Okounkov bodies
    

## References with examples worth looking at

1. [Introductory reference](./Subalgebra-Basics/StuCh11.pdf).

2. [Khovanskii bases / valuations from weight matrices](https://arxiv.org/pdf/1610.00298.pdf)

3. (related to 2, paper associated to the original package) [SAGBI bases for quotient rings](https://www.sciencedirect.com/science/article/pii/S0022404999000158)

4. Other examples worth pursuing [here](https://homepages.ecs.vuw.ac.nz/foswiki/pub/Users/Donelan/WebHome/multiscrews.pdf), and [here](https://arxiv.org/abs/1809.01026), and [here](https://arxiv.org/pdf/1702.05480.pdf), and [here](https://arxiv.org/abs/0803.0892), and [here](https://arxiv.org/abs/1612.03838), and ...

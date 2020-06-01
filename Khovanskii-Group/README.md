# History 

1. Around 1997, Mike Stillman and Harry Tsai created some code to compute SAGBI bases. Two versions of the SAGBI code were written, one just using the front-end of Macaulay2 and the other built into the C/C++ engine of Macaulay2. 
2. Around 2006, Mike Sillman created a package called SubalgebraBases to provide documentation and the front-end package for computing SAGBI bases. Sometime before 2018, a change to the way monomial orders are represented broke at least one of the two implementations of SAGBI bases. 
3. In the Summer of 2019 (at the David Cox retirement conference), Michael Burr, Tim Duff, and Mike Stillman fixed the error in the engine due to the change in monomial orders for computing SAGBI bases. 
4. In the Summer of 2020 (at the Cleveland M2 workshop), Michael Burr, Tim Duff, and Elise Walker worked on modernizing and cleaning the code. Some future projects were also planned. 

# Where's the Code?

1. [Version on the master branch](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases.m2). All tests pass since the release of Macaulay2 1.15.
2. However, the code was still a bit opaque. A stripped-down, working rewrite is available on this repository. Try installing it!

# What are (Khovanskii bases / SAGBI bases / Canonical Subalgebra Bases)?

1. See [here](./Subalgebra-Basics/StuCh11.pdf) for a first answer to this question.

# Potential Projects 

A work in progress. If you have other ideas, include them in your response to the group survey.

## New to M2

1. Write and export additional functions (eg. an exported interface to subduction, or NOBodies and NOBody volume) to compute interesting aspects of SAGBI/Khovanskii bodies. This should tie into other packages where appropriate (Tropical, NormalToricVarieties, FourTiTwo.)

2. Add / clean up the documentation, add (additional) [tests](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases/sagbi-tests.m2), and prepare the package for future publication. 

3. (Related to 2) Develop documented examples that give an introduction at the level of [this introductory chapter](./Subalgebra-Basics/StuCh11.pdf).

## Some experience with M2

1. At the Cleveland M2 workshop, we rewrote and commented the old code. There are two functions related to auto-reduction that remain somewhat mysterious. Figure out what these functions are trying to do and make sure that they are sane. Then, rewrite them to be simpler. 

2. Our working template derives from the [engine-based code](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases/sagbieng.m2). Most of the [top level code](https://github.com/Macaulay2/M2/blob/master/M2/Macaulay2/packages/SubalgebraBases/sagbitop.m2) is duplicated, but there are some differences. Try to re-incorporate any novel elements of the top-level code.

3. (Related to 1 and 2). Performance-tuning. Develop different strategies for intensive computations---for instance, when computing binomial syzygies, the package "FourTiTwo" gives a potential optimization.

## Expert at M2

1. Sort out what the basic data types should be. On this repository, the current datatype is called "Subring". The intended usage is for k-subalgebras of polynomial rings over a field k. Compared to the old package (which had no dedicated datatype), we now stash various results of the computation inside of subrings. The interface is still primitive. Currently, Subring inherits from HashTable. Down the road, it may make sense to inherit from Ring --- however, there are some subtleties. Ideally, this type will support *intrinsic computations* (those requiring a subalgebra basis) and *extrinsic* computations (those that use Groebner bases for some presentation of the subring.)

2. (Related to 1, for extension to [Khovanskii bases](https://arxiv.org/pdf/1610.00298.pdf)) Create a good valuation type. There are many design choices that should be considered here. Should there be a single valuation type or should all common valuation constructions have their own types (inheriting from an abstract base class)? There are advantages to either approach that should be discussed. Some code has been started in this direction, but much more work (and a design) remains. This should also tie in with the TropicalGeometry package.

3. Implement additional algorithms (eg. intrinsic Groebner, toric syzygies as described [here](./Subalgebra-Basics/StuCh11.pdf).)

## References with examples worth looking at

1. [Introductory reference](./Subalgebra-Basics/StuCh11.pdf).

2. [Khovanskii bases / valuations from weight matrices](https://arxiv.org/pdf/1610.00298.pdf)

3. (related to 2, paper associated to the original package) [SAGBI bases for quotient rings](https://www.sciencedirect.com/science/article/pii/S0022404999000158)

4. Other examples worth pursuing [here](https://homepages.ecs.vuw.ac.nz/foswiki/pub/Users/Donelan/WebHome/multiscrews.pdf), and [here](https://arxiv.org/abs/1809.01026), and [here](https://arxiv.org/pdf/1702.05480.pdf), and [here](https://arxiv.org/abs/0803.0892), and [here](https://arxiv.org/abs/1612.03838), and ...

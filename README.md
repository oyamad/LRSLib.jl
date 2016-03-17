# LRSLib

LRSLib.jl is a wrapper for [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html). This module can either be used in a "lower level" using the API of lrs or using the higher level interface of [Polyhedra.jl](https://github.com/blegat/Polyhedra.jl).

As written in the [user guide of lrs](http://cgm.cs.mcgill.ca/~avis/C/lrslib/USERGUIDE.html#Introduction):
> A polyhedron can be described by a list of inequalities (H-representation) or as by a list of its vertices and extreme rays (V-representation). lrs is a C program that converts a H-representation of a polyhedron to its V-representation, and vice versa.  These problems are known respectively at the vertex enumeration and convex hull problems.

The current version of `lrs` only works with files as input and output. I have have added the support for creating H-representation and V-representation directly from a matrix on a [fork of lrs](https://github.com/blegat/lrslib).
I still need to rewrite `lrs_main` to support the conversions from H-representation to V-representation, and vice versa.
If you want to do those conversions in Julia, use the cdd wrapper [CDDLib.jl](https://github.com/blegat/Polyhedra.jl) which is fully working.
If your code use the higher level interface of [Polyhedra.jl](https://github.com/blegat/Polyhedra.jl), you will only have to change 3 characters to use `LRSLib` instead of `CDDLib` once `LRSLib` is fully working.

[![Build Status](https://travis-ci.org/blegat/LRSLib.jl.svg?branch=master)](https://travis-ci.org/blegat/LRSLib.jl)

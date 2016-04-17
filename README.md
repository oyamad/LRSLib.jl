# LRSLib

LRSLib.jl is a wrapper for [lrs](http://cgm.cs.mcgill.ca/~avis/C/lrs.html). This module can either be used in a "lower level" using the API of lrs or using the higher level interface of [Polyhedra.jl](https://github.com/blegat/Polyhedra.jl).

As written in the [user guide of lrs](http://cgm.cs.mcgill.ca/~avis/C/lrslib/USERGUIDE.html#Introduction):
> A polyhedron can be described by a list of inequalities (H-representation) or as by a list of its vertices and extreme rays (V-representation). lrs is a C program that converts a H-representation of a polyhedron to its V-representation, and vice versa.  These problems are known respectively at the vertex enumeration and convex hull problems.

I have have [fork of lrs](https://github.com/blegat/lrslib) to add a few functions to help doing the wrapper.
These changes are not upstream yet so this version is used instead of the upstream version.

[![Build Status](https://travis-ci.org/blegat/LRSLib.jl.svg?branch=master)](https://travis-ci.org/blegat/LRSLib.jl)

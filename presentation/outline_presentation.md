# Outline PhD presentation (max 30 min)



 1. Short outline of what I did, and did not. Mention the unsuccessful attempt to produce hyper-Kähler varieties. The smoothing components of `C(dP6)`. The construction of new Calabi-Yau varieties (which I will explain here). (~5 min)
 3. Explain SR-schemes: as unions of `P^d`'s. Mention join (that join of spheres are spheres). Show the hexagon (and mention its ideal). (~5 min) (DROP THIS? NOT EXTREMELY RELEVANT)
 3. Introduce Calabi-Yau manifolds. (~ 3 min)
 4. Mention mirror symmetry (strange duality between different CY manifolds). Explain orbifolding heuristic and Roan's formula (<- very constructive). (~ 5 min)
 5. Motivate my construction: comes from the two embeddings of `dP6 -> P^6` and `dP6 -> P^7`. (~ 3 min)
 6. Explain construction (is smoothing of `E6*E6`). Explain calculation of Euler characteristic, and also heuristic for Hodge numbers. (~ 7 min)
 7. Explain orbifold construction. (~ 5 min)
 8. Finish (1 min)

Sum tid: 35 min (for mye!)

---

## Short outline

The thesis consists of three main themes. The first chapter describes an unsuccessful attempt to find new families of hyper-Kähler varieties. In the second chapter I study a special toric singularity: the cone over a del Pezzo surface of degree 6. It is a 3-dimensional toric variety, having two smoothing components. I describe the smoothing components, and I calculate the singular homology groups of the two smoothings.

In this presentation I will focus on my results from the last chapter, namely the construction of new Calabi-Yau manifolds, and candidates for their mirrors.

## Introduction to Calabi-Yau manifolds

Here is the definition of Calabi-Yau manifolds (for us, they all have dimension 3):

 - Smooth and projective (sometimes allow mild singularities).
 - `H^0(X,OO_X)=H^3(X,OO_X)=k` and `h^i(X,OO_X)=0` for `i=1,2`.
 - Trivial canonical bundle. `omega_X=OO_X`.

Thousands of them have been found as anticanonical sections of toric varieties, and there are many other (mostly quite ad hoc) constructions. The canonical first example of a Calabi-Yau manifold is the quintic threefold in `P^4`. It is the hypersurface defined by a general degree 5 polynomial in 5 variables.

It has Hodge numbers `(1,101)`. One can find the Euler characteristic by using Chern classes. The one counts the number of moduli parameters: the number of degree 5 polynomials is 125. But we can subtract by automorphisms of `P^4`, of which there `5*5-1=24`. Hence `h^11=125-24=101`. The same kind of reasoning was used to predict Hodge numbers for my constructions.

## Mirror symmetry and orbifolding

Connected to Calabi-Yau threefolds is the phenomenon of mirror symmetry. In its simplest form, it states that for each Calabi-Yau variety `X`, there is a "mirror partner" `X'`, whose Hodge diamond is the Hodge diamond of `X` mirrored diagonally.

This simple form is "almost true" (there are counter-examples, coming from cases with `h^11=0`). There are more precise formulations (nowadays mirror symmetry is mostly considered as an equivalence of two derived categories).

Constructing mirrors from given Calabi-Yau manifolds is interesting because it might give us some insight into mirror symmetry. Let us describe one heuristic which often produces a mirror, called _orbifolding_.

The idea is this: suppose you have a family of Calabi-Yau varieties, where the general member `X` is smooth. Suppose further that there is a subfamily invariant under some finite group, where the general member `X_0` have isolated singularities. On this family there might be some action of a finite subgroup `H` of the big torus of `P^n`. One then considers `X_0/H`, and a crepant resolution `X'` of it. In many cases of interest, `X` and `X'` are mirror manifolds.

This is for example the case for the quintic: here one considers `f=x_0x_1x_2x_3x_4 + ∑x_i^5`, which is invariant under `S_5`. There are isolated singularities, but using toric geometry, one can resolve them. The resolved manifold is a mirror to the quintic.

Given such a crepant resolution, one wants to verify that the mirror candidate at least have the correct Euler characteristic (minus that of the original). For this, we have Roan's formula, which lets us compute the Euler characteristic in terms of the number of fixed points for elements from the torus group acting. (show the formula).

## Motivate my constructions

The singularity `C(dP6)` have two smoothing components, coming from different ways of writing its set of equations. This is because `dP6` lies naturally in both `P^5` and in `P^6`. In the first presentation, the equations are given by the minors of a `3x3` matrix (show it). In the second presentation, the equations are given by the "minors" of a `2x2x2` cube.

## My constructions

The idea is to join (literally) the equations of the smoothings of `C(dP6)` in the three possible ways (giving varieties in `P^15`, `P^16`, and `P^17`).


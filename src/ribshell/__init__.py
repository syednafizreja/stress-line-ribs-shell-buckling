"""Stress-line ribs vs optimised thickness (Project 2).

Modules:
    constants  model data fixed by the brief
    form       form-found surface z(x, y)
    mesh       Gmsh plan mesh lifted onto the surface, with embedded rib curves
    fe         CalculiX decks (static G, perturbation *BUCKLE with a_h) and result readers
    trace      tensor-line ribs (uses the shared tracer in cecore)
    opt        variant E, CMA-ES over thickness zones

Units are SI throughout: m, N, Pa, kg.
"""
__version__ = "0.1.0"

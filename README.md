# Gravitational N-Body Simulation with Barnes-Hut Algorithm

This repository contains an efficient implementation of the gravitational N-body problem using the Barnes-Hut algorithm. The N-body problem deals with simulating the motion of celestial bodies like planets and stars interacting through gravitational forces.

## Problem Description

The simulation employs the Plummer sphere model, utilizing Newton's law for gravitational forces. Instead of a simple time integration, the Velocity Verlet method is employed for improved numerical accuracy. The Barnes-Hut algorithm is utilized to reduce computational complexity from O(n^2) to O(n log n) by approximating forces exerted on distant groups of objects.

## Solution Method

The code implements a quadtree structure for efficient force calculations. Each node in the quadtree represents a region of space containing particles. Forces are recursively calculated between particles or approximated using group centers of mass. The implementation optimizes memory access for improved performance and parallelizes force calculations using OpenMP.

## Experiments

Various experiments were conducted to measure runtime performance and optimize parameters such as Θmax and timestep (∆t). Results indicate significant runtime improvements with optimized parameter values.

## Conclusions

Future optimizations could include modifying the tree structure instead of rebuilding it entirely at each step and reducing memory indirection for faster access to data. These enhancements could further improve the efficiency of the simulation.

For detailed information on implementation, experiments, and conclusions, refer to the provided documentation.

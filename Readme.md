# Bee Foraging Path Optimization Using Genetic Algorithm

## Overview
This project applies a Genetic Algorithm to solve a combinatorial optimization problem similar to the **traveling salesman problem**. 

The objective is to optimize the foraging path of bees, ensuring they collect nectar from all flowers in a field in the least time possible, starting and ending at their hive located at coordinates (500, 500).

## Key Concepts
- **Genes**: Represent points (flower and hive coordinates) in a bee’s flying path. Hive coordinates are fixed at the start and end of the chromosome.
- **Chromosomes**: Sequences of genes depicting a complete foraging path.
- **Population**: A group of various foraging paths taken by the bees.
- **Fitness Function**: Evaluates path efficiency based on the total distance. Lower distance implies higher efficiency.
- **Euclidean Distance**: Used to calculate distances between points on the path. The distance `(d)` between two points `((x_1, y_1)\)` and `((x_2, y_2)\)` in 2D Cartesian coordinates is given by the formula:

$$
d = \sqrt{(x_2 - x_1)^2 + (y_2 - y_1)^2}
$$

## Genetic Algorithm Components
- **Pairs Selection**: Top 50 fittest bees are chosen for crossover. Roulette wheel selection is used, where a bee's chance of selection is proportional to its fitness.
- **Crossover ([Partially Mapped Crossover](https://github.com/ruta-tamosiunaite/partially-mapped-crossover) - PMX)**: Two parent chromosomes are combined to create offspring. Segments between two chosen crossover points are swapped and duplicates are resolved using a mapping approach.
- **Mutation (Place Change)**: Random changes in a path to explore new possibilities. Mutation rate determines the magnitude of change. Higher rates move a gene more significantly to the left or right in the chromosome.
- **Updating Population**: Population is updated by replacing less fit ancestors with fitter offspring. The combined population and offspring are sorted, and the top 100 individuals are selected for the next generation.

## Parameters
- `min_relative_separation`: Controls inbreeding by setting the minimum generational gap required for pairing. Higher values avoid close lineage connections.

## Files
- **Main.py**: Controls the algorithm flow and contains parameters.
- **Beehive.py**: Hosts the `Bee` and `Beehive` classes. The Beehive class manages the bee population and oversees their evolution through generations.

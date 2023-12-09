from beehive import Bee, Beehive, BeesArchive
import sys
import pandas as pd
import matplotlib.pyplot as plt


# H I V E  &  F L O W E R S
excel_file_path = 'Champ de pissenlits et de sauge des pres.xlsx'
df = pd.read_excel(excel_file_path)
flowers = list(zip(df['x'], df['y']))
hive_location = (500, 500)

# P A R A M E T E R S
population_size = 100
num_generations = 10
select_bees = 100
children_per_bee = 2
mutation = 1 # Mutation rate. 1 means a gene changes the place by 1 step.
mutation_frequency = 3 # Mutation frequency. Means how many genes will change its places.

# A L G O R I T H M
average_fitness_per_generation = [] # Saving average fitnesses of all generations
fittest_bee_per_generation = [] # Saving fittest bees fitnesses of all generations
bees_archive = BeesArchive() # Saving all bees to track the lineage


beehive = Beehive(flowers, hive_location, population_size)
beehive.initialize_population(method='random', randomness_ratio=0.99, swap_count=25, swap_method='neighbour_swap') # METHOD: 'random' or 'nearest_flower' or 'hybrid'. For hybrid: randomness_ratio (what percentage of population is random), swap_count (in nearest flower paths, how many genes to swap), swap_method (how to swap genes). SWAP method: 'neighbour_swap' or 'random_swap'.

for bee in beehive.population:
    bees_archive.add_bee(bee)
fitness_results = beehive.analyze_fitness()
average_fitness_per_generation.append(fitness_results["average_fitness"])
fittest_bee_per_generation.append(fitness_results["fittest_bee_fitness"])

# Analyzing differences in bees paths in initial population
if False:
    initial_population = beehive.population
    for bee in initial_population:
        x_coords, y_coords = zip(*bee.chromosome)  # Unpack the list of tuples into two lists
        plt.plot(x_coords, y_coords, marker='o')
    plt.title("Paths of All Initial Population Bees")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.show()

for generation in range(num_generations):
    #beehive.run_generation(select_bees=select_bees, number_of_children=children_per_bee, crossover_method='partially_mapped', mutation_method='place_change', mutation_rate=mutation, mutation_frequency=mutation_frequency)
    print(f'\n{generation + 1}')
    
    # Select pairs for crossover
    selected_pairs = beehive.select_pairs_for_crossover(number_of_bees=100) # Selection method = roulette wheel
    
    # Perform crossover on the selected pairs
    offspring = beehive.perform_crossover(selected_pairs, number_of_children=2, current_generation=generation+1)

    # Mutation
    if generation > 1:
        if average_fitness_per_generation[generation-1] == average_fitness_per_generation[generation-2]:
            mutation_frequency = mutation_frequency + 1
            print(f'New mutation frequency: {mutation_frequency}')
        elif fittest_bee_per_generation[generation-1] == fittest_bee_per_generation[generation-2]:
            mutation_frequency = mutation_frequency + 1
            mutation = mutation + 1
            print(f'New mutation rate: {mutation}')
            print(f'New mutation frequency: {mutation_frequency}')
    mutated_offspring = beehive.perform_mutation(offspring, mutation_method='place_change', mutation_rate=mutation, mutation_frequency=mutation_frequency)

    # Update Population
    for bee in offspring:
        bees_archive.add_bee(bee)
    beehive.update_population(mutated_offspring)

    # Calculate and Analyze Fitness
    fitness_results = beehive.analyze_fitness()
    average_fitness_per_generation.append(fitness_results['average_fitness'])
    fittest_bee_per_generation.append(fitness_results['fittest_bee_fitness'])
    if generation == num_generations - 1:
        fittest_bee_path = fitness_results['fittest_bee_chromosome']
        fittest_bee_id = fitness_results['fittest_bee_id']

print(f'Fittest bee fitness: {fittest_bee_per_generation[-1]}')
if fittest_bee_per_generation[-1] < 8500:
    print(fittest_bee_path)
    x_coords, y_coords = zip(*fittest_bee_path)  # Unpack the list of tuples into two lists
    plt.plot(x_coords, y_coords, marker='o')
    plt.title("Fittest bee path")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.show()

    # Plotting the timeline of fitness averages
    plt.figure(figsize=(12, 6))
    plt.plot(average_fitness_per_generation, label='Average Fitness', marker='o')
    plt.plot(fittest_bee_per_generation, label='Fittest Bee Fitness', marker='x')
    plt.title("Fitness Over Generations")
    plt.xlabel("Generation")
    plt.ylabel("Fitness")
    plt.legend()
    plt.grid(True)
    plt.show()

print()
for bee in bees_archive.all_bees:
    # print(bee) this prints all individuals
    pass
print()

# Lineage of the fittest bee
fittest_bee = bees_archive.get_bee_by_id(fittest_bee_id)
print(fittest_bee.unique_id)
print(fittest_bee.parent1_id)
print(fittest_bee.parent2_id)
if fittest_bee:
    parent1 = bees_archive.get_bee_by_id(fittest_bee.parent1_id)
    parent2 = bees_archive.get_bee_by_id(fittest_bee.parent2_id)
    print(parent1)
    print(parent2)
    # Now you have the fittest bee and its parents

def build_lineage(bee_id, tree, bees_archive, depth=0):
    bee = bees_archive.get_bee_by_id(bee_id)
    if bee is None:
        return

    # Add the current bee to the tree
    tree[bee_id] = {'depth': depth, 'parent1_id': bee.parent1_id, 'parent2_id': bee.parent2_id}

    # Recursively add parents
    if bee.parent1_id is not None:
        build_lineage(bee.parent1_id, tree, bees_archive, depth + 1)
    if bee.parent2_id is not None:
        build_lineage(bee.parent2_id, tree, bees_archive, depth + 1)

# Initialize an empty dictionary to hold the tree structure
family_tree = {}

# Build the lineage starting from the fittest bee
build_lineage(fittest_bee_id, family_tree, bees_archive)


import networkx as nx
import matplotlib.pyplot as plt

# Create a directed graph
G = nx.DiGraph()

# Add nodes and edges based on the family tree
for bee_id, info in family_tree.items():
    G.add_node(bee_id)
    if info['parent1_id'] is not None:
        G.add_edge(info['parent1_id'], bee_id)
    if info['parent2_id'] is not None:
        G.add_edge(info['parent2_id'], bee_id)

# Generate positions for each node using a hierarchical layout
pos = nx.spring_layout(G, k=0.5)  # k is a spacing parameter

# Draw the graph
nx.draw(G, pos, with_labels=True, arrows=True)
plt.title("Family Tree of the Fittest Bee")
plt.show()

for bee_id, info in family_tree.items():
    print(f"Bee ID: {bee_id}, Parent 1 ID: {info['parent1_id']}, Parent 2 ID: {info['parent2_id']}")

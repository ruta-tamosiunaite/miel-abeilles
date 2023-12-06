from beehive import Beehive, BeesArchive
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
num_generations = 100
mutation = 1 # Mutation rate. 1 means a gene changes the place by 1 step.
mutation_frequency = 3 # Mutation frequency. Means how many genes will change its places.

# A L G O R I T H M
average_fitness_per_generation = [] # Saving average fitnesses of all generations
fittest_bee_per_generation = [] # Saving fittest bees fitnesses of all generations
bees_archive = BeesArchive() # Saving all bees to track the lineage

beehive = Beehive(flowers, hive_location, population_size)
beehive.initialize_population(method='random', randomness_ratio=0.99, swap_count=25, swap_method='neighbour_swap') # METHOD: 'random' or 'nearest_flower' or 'hybrid'. For hybrid: randomness_ratio (what percentage of population is random), swap_count (in nearest flower paths, how many genes to swap), swap_method (how to swap genes). SWAP method: 'neighbour_swap' or 'random_swap'.

bees_archive.update_archive(beehive.population)
fitness_results = beehive.analyze_fitness()
average_fitness_per_generation.append(fitness_results["average_fitness"])
fittest_bee_per_generation.append(fitness_results["fittest_bee_fitness"])

initial_population = beehive.population
# Analyzing differences in bees paths in initial population
if False:
    for bee in initial_population:
        x_coords, y_coords = zip(*bee.chromosome)  # Unpack the list of tuples into two lists
        plt.plot(x_coords, y_coords, marker='o')
    plt.title("Paths of All Initial Population Bees")
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.show()

for generation in range(num_generations):
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
    bees_archive.update_archive(offspring)
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

# Lineage of the fittest bee



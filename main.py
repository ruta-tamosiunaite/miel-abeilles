from beehive import Beehive
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
num_generations = 50
mutation = 1 # Mutation rate. 1 means 1 gene changes the place by 1

# A L G O R I T H M
beehive = Beehive(flowers, hive_location, population_size)
beehive.initialize_population(method='hybrid', randomness_ratio=0, swap_count=5, swap_method='neighbour_swap') # METHOD: 'random' or 'nearest_flower' or 'hybrid'. SWAP method: 'neighbour_swap' or 'random_swap'.

initial_population = beehive.population
# Analyzing differences in bees paths in initial population
for bee in initial_population:
    x_coords, y_coords = zip(*bee.chromosome)  # Unpack the list of tuples into two lists
    plt.plot(x_coords, y_coords, marker='o')

plt.title("Paths of All Initial Population Bees")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

average_fitness_per_generation = [] # Saving average fitnesses of all generations
fittest_bee_per_generation = [] # Saving fittest bees fitnesses of all generations
fitness_results = beehive.analyze_fitness()
average_fitness_per_generation.append(fitness_results["average_fitness"])
fittest_bee_per_generation.append(fitness_results["fittest_bee_fitness"])

for generation in range(num_generations):
    print(f'\n{generation + 1}')
    
    # Select pairs for crossover
    selected_pairs = beehive.select_pairs_for_crossover(number_of_bees=50) # Selection method = roulette wheel
    
    # Perform crossover on the selected pairs
    offspring = beehive.perform_crossover(selected_pairs, crossover_method="partially_mapped", number_of_children=1, current_generation=generation+1)

    # Mutation
    if generation > 1:
        #print(f'current generation fitness = {average_fitness_per_generation[generation]}')
        #print(f'previous generation fitness = {average_fitness_per_generation[generation-1]}')
        #print(f'before previous generation fitness = {average_fitness_per_generation[generation-2]}')
        if average_fitness_per_generation[generation-1] == average_fitness_per_generation[generation-2]:
            mutation = mutation + 1
            print(f'New mutation rate: {mutation}')
    mutated_offspring = beehive.perform_mutation(offspring, mutation_method='place_change', mutation_rate=mutation)

    # Update Population
    beehive.update_population(mutated_offspring)

    # Calculate and Analyze Fitness
    fitness_results = beehive.analyze_fitness()
    average_fitness_per_generation.append(fitness_results['average_fitness'])
    fittest_bee_per_generation.append(fitness_results['fittest_bee_fitness'])
    if generation == num_generations - 1:
        fittest_bee_path = fitness_results['fittest_bee_chromosome']

for generation, fitness in enumerate(average_fitness_per_generation):
    #print(generation)
    #print(fitness)
    pass

print(f'Fittest bee fitness: {fittest_bee_per_generation[-1]}')
x_coords, y_coords = zip(*fittest_bee_path)  # Unpack the list of tuples into two lists
plt.plot(x_coords, y_coords, marker='o')
plt.title("Fittest bee path")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.show()

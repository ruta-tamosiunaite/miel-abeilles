from beehive import Beehive
import sys
import pandas as pd
import matplotlib.pyplot as plt

# P A R A M E T E R S
min_fittest_bee_distance_in_initial_population = None
min_initial_population_average_fitness = None
select_TOP_bees_for_crosover = 50
pairing_method = None
pairing_rules = [True, None, None] # [avoid_inbreeding, bee_reproduction_times_limit, pair_reproduction_times_limit]
crossover_method = None
crossover_rate = 0.5 # 1 child from pair. 1 would be 2 childs from pair
mutation_rate_is_evolutive = False
some_condition_for_mutation = None
new_mutation_rate = None
mutation_rate = None
update_method = None

# Load flower coordinates and hive location
excel_file_path = 'Champ de pissenlits et de sauge des pres.xlsx'
df = pd.read_excel(excel_file_path)
flowers = list(zip(df['x'], df['y']))
hive_location = (500, 500)

# Set parameters for the genetic algorithm
population_size = 100
num_generations = 55

# A L G O R I T H M
beehive = Beehive(flowers, hive_location, population_size)
beehive.initialize_population(min_fittest_bee_distance_in_initial_population=20000, min_initial_population_average_fitness=25500)
sys.exit()
average_fitness_per_generation = []
fitness_results = beehive.analyze_fitness()
average_fitness_per_generation.append(fitness_results["average_fitness"])

for generation in range(num_generations):
    print(f'\n{generation + 1}')
    
    # Select pairs for crossover
    selected_pairs = beehive.select_pairs_for_crossover(selection_method="roulette_wheel", min_relative_separation=5)
    
    # Perform crossover on the selected pairs
    offspring = beehive.perform_crossover(selected_pairs, crossover_method="partially_mapped", number_of_children=1, current_generation=generation+1)

    # Mutation
    if mutation_rate_is_evolutive == True:
        if some_condition_for_mutation:
            mutation_rate = new_mutation_rate
    mutated_offspring = beehive.perform_mutation(offspring, mutation_method='place_change', mutation_rate=1)

    # Update Population
    beehive.update_population(mutated_offspring)

    # Calculate and Analyze Fitness
    fitness_results = beehive.analyze_fitness()
    average_fitness_per_generation.append(fitness_results["average_fitness"])
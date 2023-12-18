from beehive import Beehive, BeesArchive
import matplotlib.pyplot as plt
import networkx as nx

def collect_lineage_data(bee_id, generation=0, max_generations=2, lineage_data=None):
    if generation > max_generations or bee_id is None:
        return

    if lineage_data is None:
        lineage_data = []

    bee = bees_archive.get_bee_by_id(bee_id)
    bee_data = {
        'id': bee.unique_id,
        'parent1_id': getattr(bee, 'parent1_id', None),
        'parent2_id': getattr(bee, 'parent2_id', None)
    }
    lineage_data.append(bee_data)

    collect_lineage_data(bee.parent1_id, generation + 1, max_generations, lineage_data)
    collect_lineage_data(bee.parent2_id, generation + 1, max_generations, lineage_data)

    return lineage_data

# H I V E  &  F L O W E R S
flowers = [(796, 310), (774, 130), (116, 69), (908, 534), (708, 99), (444, 428), (220, 307), (501, 287), (345, 560), (628, 311), (901, 639), (436, 619), (938, 646), (45, 549), (837, 787), (328, 489), (278, 434), (704, 995), (101, 482), (921, 964), (493, 970), (494, 898), (929, 389), (730, 742), (528, 794), (371, 429), (98, 711), (724, 631), (573, 903), (964, 726), (213, 639), (549, 329), (684, 273), (273, 105), (897, 324), (508, 31), (758, 405), (862, 361), (898, 898), (2, 897), (951, 209), (189, 739), (602, 68), (437, 601), (330, 410), (3, 517), (643, 404), (875, 407), (761, 772), (276, 666)]
hive_location = (500, 500)

# P A R A M E T E R S
population_size = 100
num_generations = 500
select_bees = 100 # How many bees to select for crossover
mutation_rate = 1 # Mutation rate. 1 means a random gene changes the place by 1 step.
mutation_frequency = 3 # Mutation frequency. Means how many genes will change its places.
random_seed = None

# A L G O R I T H M
average_fitness_per_generation = [] # Saving average fitnesses of all generations
fittest_bee_per_generation = [] # Saving fittest bees fitnesses of all generations
bees_archive = BeesArchive() # Saving all bees to track the lineage

beehive = Beehive(flowers, hive_location, population_size, bees_archive, seed=random_seed)

fitness_results = beehive.analyze_fitness()
average_fitness_per_generation.append(fitness_results["average_fitness"])
fittest_bee_per_generation.append(fitness_results["fittest_bee_fitness"])

for generation in range(num_generations):
    print(f'\n{generation + 1}')
    new_population = beehive.run_generation(bees_archive, select_bees=select_bees, mutation_rate=mutation_rate, mutation_frequency=mutation_frequency, current_generation=generation+1)
    
    if generation > 1:
        if average_fitness_per_generation[generation-1] == average_fitness_per_generation[generation-2]:
            mutation_frequency = mutation_frequency + 1
            print(f'New mutation frequency: {mutation_frequency}')
        elif fittest_bee_per_generation[generation-1] == fittest_bee_per_generation[generation-2]:
            mutation_frequency = mutation_frequency + 1 # If 2 last fittest bees' fitnesses are equal, increase mutation rate and mutation frequency by one
            mutation_rate = mutation_rate + 1
            print(f'New mutation rate: {mutation_rate}')
            print(f'New mutation frequency: {mutation_frequency}')

        # Check if the first and last bees have identical chromosomes
        if beehive.population[0].chromosome == beehive.population[-1].chromosome:
            print("First and last bees have identical chromosomes. Stopping further generations.")
            break

    average_fitness_per_generation.append(new_population['average_fitness'])
    fittest_bee_per_generation.append(new_population['fittest_bee_fitness'])

fittest_bee_id = new_population['fittest_bee_id']
fittest_bee_path = new_population['fittest_bee_chromosome']
print(f'Fittest bee ID: {fittest_bee_id}')
print(f'Fittest bee fitness: {new_population['fittest_bee_fitness']}')
print(fittest_bee_path)

# Displaying the fittest bee path
x_coords, y_coords = zip(*fittest_bee_path)
plt.plot(x_coords, y_coords, marker='o')
plt.title(f"Fittest bee path - Fitness {int(new_population['fittest_bee_fitness'])} - Seed {beehive.seed}")
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

# Displaying the lineage of the fittest bee
bee_lineage_data = collect_lineage_data(fittest_bee_id)
G = nx.DiGraph()
for bee in bee_lineage_data:
    G.add_node(bee['id'])
    if bee['parent1_id']:
        G.add_edge(bee['parent1_id'], bee['id'])
    if bee['parent2_id']:
        G.add_edge(bee['parent2_id'], bee['id'])
pos = nx.spring_layout(G)  # positions for all nodes
nx.draw(G, pos, with_labels=True, node_color='skyblue', node_size=1500, edge_color='black', linewidths=1, font_size=10)
plt.show()
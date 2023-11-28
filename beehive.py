import sys
import math
import random


class Bee:
    def __init__(self, chromosome):
        self.unique_id = None
        self.chromosome = chromosome
        self.fitness = self.calculate_fitness(chromosome)
        self.chromosome_before_mutation = None
        self.fitness_before_mutation = None
        self.parent1 = None
        self.parent2 = None
        self.normalized_fitness = None # variable depending on what population bee is at the moment of calculation
        
    def calculate_fitness(self, chromosome):
        def calculate_distance(point1, point2):
            return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

        total_distance = 0
        for i in range(len(chromosome) - 1):
            total_distance += calculate_distance(chromosome[i], chromosome[i+1])
        return total_distance
    
class Beehive:
    def __init__(self, flowers, hive_location, population_size):
        self.flowers = flowers
        self.hive_location = hive_location
        self.population_size = population_size
        self.population = []  # List of Bee objects
        self.archive = []  # Archive for lineage tracking

    def update_archive(self, bee):
        # Add an entry to the archive for the given bee
        self.archive[bee.unique_id] = {
            'parent1': bee.parent1,
            'parent2': bee.parent2
        }

    def check_lineage(self, bee1_id, bee2_id, min_relative_separation):
        """Check if two bees share common ancestors."""
        ancestors_bee1 = self.get_ancestors(bee1_id, min_relative_separation)
        ancestors_bee2 = self.get_ancestors(bee2_id, min_relative_separation)
        return not ancestors_bee1.isdisjoint(ancestors_bee2)

    def get_ancestors(self, bee_id, min_relative_separation):
        """Get ancestors of a given bee up to a certain number of generations back."""
        #print(f'getting ancestors for bee {bee_id}')
        ancestors = set()
        current_id = bee_id
        generation_count = 0

        while current_id in self.archive and generation_count < min_relative_separation:
            parent1 = self.archive[current_id]['parent1']
            parent2 = self.archive[current_id]['parent2']
            if parent1:
                ancestors.add(parent1)
                current_id = parent1
            elif parent2:
                ancestors.add(parent2)
                current_id = parent2
            else:
                break
            generation_count += 1
        return ancestors

    def initialize_population(self, min_fittest_bee_distance_in_initial_population=20000, min_initial_population_average_fitness=18300):        
        print('Initializing population')
        while True:
            self.population.clear()
            for _ in range(self.population_size):
                chromosome = self.create_random_chromosome()
                new_bee = Bee(chromosome=chromosome)
                self.population.append(new_bee)

            results = self.analyze_fitness()
            average_fitness = results['average_fitness']
            fittest_bee_fitness = results['fittest_bee_fitness']
            if fittest_bee_fitness <= min_fittest_bee_distance_in_initial_population and average_fitness <= min_initial_population_average_fitness:
                break

        # Sort the bees based on fitness and give unique IDs
        self.population.sort(key=lambda bee: bee.fitness)
        for index, bee in enumerate(self.population, start=1):
            bee.unique_id = f"0-{index}"
            print(f"Bee ID: {bee.unique_id}, Fitness: {bee.fitness}")
    
    def create_random_chromosome(self):
        path = self.flowers.copy()
        random.shuffle(path)
        return [self.hive_location] + path + [self.hive_location]
    
    def select_pairs_for_crossover(self, selection_method="roulette_wheel", min_relative_separation=0):
        selected_bees = sorted(self.population, key=lambda bee: bee.fitness)[:50]
        '''print('Bees selected for crossover:')
        for bee in selected_bees:
            print(f"Bee ID: {bee.unique_id}, Fitness: {bee.fitness}")
        print(f'Selecting pairs (min_relative_separation={min_relative_separation}):')'''
        pairs = []
        already_chosen = set()  # To keep track of bees that have been paired
        while len(pairs) < 25 and selected_bees:
            pair = []
            while len(pair) < 2:
                if not pair:
                    chosen_bee = self.roulette_wheel_selection(selected_bees, already_chosen)
                    if chosen_bee:
                        pair.append(chosen_bee)
                        already_chosen.add(chosen_bee)
                else:
                    # Find possible pairs for the first chosen bee
                    possible_pairs = [bee for bee in selected_bees if bee not in already_chosen and not self.check_lineage(pair[0].unique_id, bee.unique_id, min_relative_separation)]
                    
                    # If no possible pairs are found, adjust min_relative_separation
                    min_relative_separation_change = 0
                    while not possible_pairs:
                        min_relative_separation_change += 1
                        # Code to adjust min_relative_separation or handle no possible pairs
                        print(f'No possible pairs for bee {chosen_bee.unique_id}. Lowering min_relative_separation to {min_relative_separation-min_relative_separation_change}')
                        possible_pairs = [bee for bee in selected_bees if bee not in already_chosen and not self.check_lineage(pair[0].unique_id, bee.unique_id, min_relative_separation - min_relative_separation_change)]
                        
                    # Select second bee from possible pairs
                    second_bee = self.roulette_wheel_selection(possible_pairs, already_chosen)
                    if second_bee:
                        pair.append(second_bee)
                        already_chosen.add(second_bee)
            if pair:
                pairs.append(tuple(pair))
        '''print('\nPairs selected:')
        for no, pair in enumerate(pairs, start=1):
            print(no)
            print(f'{pair[0].unique_id}: {pair[0].fitness}')
            print(f'{pair[1].unique_id}: {pair[1].fitness}')
            print()'''
        return pairs

    def roulette_wheel_selection(self, bees, already_chosen):
        # C H E C K  H E R E  F O R  T H E  E R R O R - need to use normalized fitness. Recalculate normalized fitness for possible pairs before each wheel turn. 
        S = sum(bee.fitness for bee in bees if bee not in already_chosen)
        r = random.uniform(0, S)
        P = 0
        for bee in bees:
            if bee not in already_chosen:
                P += bee.fitness
                if P >= r:
                    return bee
        return None


    def perform_crossover(self, pairs, crossover_method="partially_mapped", number_of_children=1, current_generation=0):
        #print(f'Performing crossover of pairs. Number of children per bee: {number_of_children}')
        offspring = []

        for parent1, parent2 in pairs:
            if crossover_method == "partially_mapped":
                children = self.partially_mapped_crossover(parent1.chromosome, parent2.chromosome, number_of_children)

                if number_of_children == 1:
                    offspring.append(self.create_bee(children, parent1.unique_id, parent2.unique_id))
                else:
                    for child_chromosome in children:
                        offspring.append(self.create_bee(child_chromosome, parent1.unique_id, parent2.unique_id))
            else:
                raise NotImplementedError(f"Crossover method '{crossover_method}' is not implemented.")
        # Sort the offspring by fitness and assign unique IDs
        offspring.sort(key=lambda bee: bee.fitness)
        
        for rank, bee in enumerate(offspring, start=1):
            bee.unique_id = f"{current_generation}-{rank}"
            '''print(bee.unique_id)
            print(bee.fitness)
            print()'''
        return offspring

    def partially_mapped_crossover(self, parent1, parent2, number_of_children=1):
        size = len(parent1)
        # Step 1: Select crossover range at random
        start, end = sorted(random.sample(range(1, size - 1), 2))  # Avoid the first and last gene (the hive)

        # Step 2: Create offspring by exchanging the selected range
        child1 = parent1[:start] + parent2[start:end] + parent1[end:]
        child2 = parent2[:start] + parent1[start:end] + parent2[end:]

        # Step 3: Determine the mapping relationship to legalize offspring
        mapping1 = {parent2[i]: parent1[i] for i in range(start, end)}
        mapping2 = {parent1[i]: parent2[i] for i in range(start, end)}

        # Step 4: Legalize children with the mapping relationship
        for i in list(range(start)) + list(range(end, size)):
            if child1[i] in mapping1:
                while child1[i] in mapping1:
                    child1[i] = mapping1[child1[i]]
            if child2[i] in mapping2:
                while child2[i] in mapping2:
                    child2[i] = mapping2[child2[i]]
        if number_of_children == 1:
            if random.random() < 0.5:
                return child1
            else:
                return child2
        else:
            return child1, child2
    
    def create_bee(self, chromosome, parent1_id, parent2_id):
        new_bee = Bee(chromosome=chromosome)
        new_bee.parent1 = parent1_id
        new_bee.parent2 = parent2_id
        return new_bee
    
    def perform_mutation(self, bees, mutation_method='place_change', mutation_rate=1):
        for bee in bees:
            if mutation_method == 'place_change':
                self.mutate_bee_place_change(bee, mutation_rate)
            else:
                raise NotImplementedError(f"Mutation method '{mutation_method}' is not implemented.")
        # Sort the offspring by fitness and assign unique IDs
        bees.sort(key=lambda bee: bee.fitness)
        '''print(f'Performing mutation. Mutation method: {mutation_method}. Mutation rate: {mutation_rate}')
        for child in bees:
            print(child.unique_id)
            print(child.fitness)
            print()'''
        return bees
    
    def mutate_bee_place_change(self, bee, mutation_rate):
        # Store the original chromosome and fitness
        bee.chromosome_before_mutation = bee.chromosome[:]
        bee.fitness_before_mutation = bee.fitness

        # Select a gene to mutate (excluding the hive locations)
        gene_index = random.randint(1, len(bee.chromosome) - 2)

        # Determine the move direction and magnitude
        move_direction = random.choice([-1, 1])
        move_steps = min(mutation_rate, gene_index if move_direction == -1 else len(bee.chromosome) - 1 - gene_index)

        # Perform the mutation by moving the gene
        new_position = gene_index + move_direction * move_steps
        gene_to_move = bee.chromosome[gene_index]
        bee.chromosome.pop(gene_index)
        bee.chromosome.insert(new_position, gene_to_move)

        # Update fitness after mutation and revert if not better
        new_fitness = bee.calculate_fitness(bee.chromosome)
        if new_fitness > bee.fitness_before_mutation:
            # Revert to pre-mutation chromosome if fitness has worsened
            bee.chromosome = bee.chromosome_before_mutation
            bee.fitness = bee.fitness_before_mutation
        else:
            bee.fitness = new_fitness

    def update_population(self, offspring):
        """Update the population with fitter offspring."""
        # Combine current population and offspring
        combined_population = self.population + offspring

        # Sort combined population by fitness (lower is better)
        combined_population.sort(key=lambda bee: bee.fitness)

        # Keep only the fittest bees up to the population size
        self.population = combined_population[:self.population_size]

        '''print('Updating population with fitter offspring.')
        for bee in self.population:
            print(f"Bee ID: {bee.unique_id}, Fitness: {bee.fitness}")
        print(f'Population size: {len(self.population)}')'''

    def analyze_fitness(self):
        """Calculate and analyze the fitness of the population."""
        #print(f'\nAnalyzing fitness of the population')

        # Calculate the average fitness of the population
        total_fitness = sum(bee.fitness for bee in self.population)
        average_fitness = total_fitness / len(self.population)

        # Find the fittest bee (lowest fitness)
        fittest_bee = min(self.population, key=lambda bee: bee.fitness)

        # Print results
        print(f"Average Fitness of Population: {average_fitness:.2f}")
        print(f"Fittest Bee ID: {fittest_bee.unique_id}, Fitness: {fittest_bee.fitness}")

        # Return a dictionary with the results
        return {
            "average_fitness": average_fitness,
            "fittest_bee_id": fittest_bee.unique_id,
            "fittest_bee_chromosome:": fittest_bee.chromosome,
            "fittest_bee_fitness": fittest_bee.fitness
        }
    
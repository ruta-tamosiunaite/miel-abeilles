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

    def ensure_hive_location(self, chromosome):
        """Ensure that the chromosome starts and ends with the hive location."""
        if chromosome[0] != self.hive_location:
            print('0 is NOT THE HIVE')
        if chromosome[-1] != self.hive_location:
            print('LAST is NOT THE HIVE')

    def initialize_population(self, method='hybrid', randomness_ratio=0.3, swap_count=2, swap_method='random_swap'):
        """
        Initialize the population of bees.
        
        :param method: Method to initialize ('random', 'nearest_flower', or 'hybrid').
        :param randomness_ratio: Proportion of the population to initialize with added randomness.
        :param swap_count: Number of swaps to introduce randomness in the path.
        """
        print('Initializing population')
        num_random_bees = int(self.population_size * randomness_ratio)

        if method in ['nearest_flower', 'hybrid']:
            sorted_flowers = sorted(self.flowers, key=lambda point: calculate_distance(point, self.hive_location))

        for i in range(self.population_size):
            if method == 'random' or (method == 'hybrid' and i < num_random_bees):
                chromosome = self.create_random_chromosome()
            elif method in ['nearest_flower', 'hybrid']:
                start_point = sorted_flowers[i % len(self.flowers)]
                chromosome = self.create_path_using_nearest_flower(start_point)

                if method == 'hybrid':
                    chromosome = self.introduce_randomness_to_path(chromosome, swap_count, swap_method)

            new_bee = Bee(chromosome=chromosome)
            self.population.append(new_bee)

        # Sort the bees based on fitness and give unique IDs
        self.population.sort(key=lambda bee: bee.fitness)
        for index, bee in enumerate(self.population, start=1):
            bee.unique_id = f"0-{index}"
            print(f'Bee ID: {bee.unique_id}, Fitness: {bee.fitness}')
            #self.ensure_hive_location(bee.chromosome)

    def introduce_randomness_to_path(self, path, swap_count, method='neighbour_swap'):
        for _ in range(swap_count):
            if method == 'random_swap':
                idx1, idx2 = random.sample(range(1, len(path) - 2), 2)  # Avoid swapping the hive location
                path[idx1], path[idx2] = path[idx2], path[idx1]
            elif method == 'neighbour_swap':
                gene_index = random.randint(1, len(path) - 2)
                # Determine move direction (backward or forward)
                if gene_index == 1:  # If the gene is next to the starting hive, it can only move forward
                    move_direction = 1
                elif gene_index == len(path) - 2:  # If the gene is next to the ending hive, it can only move backward
                    move_direction = -1
                else:  # Otherwise, randomly choose to move backward or forward
                    move_direction = random.choice([-1, 1])
                # Perform the move by swapping the selected gene with its neighbor
                neighbor_index = gene_index + move_direction
                path[gene_index], path[neighbor_index] = path[neighbor_index], path[gene_index]
        return path
        
    def create_random_chromosome(self):
        path = self.flowers.copy()
        random.shuffle(path)
        return [self.hive_location] + path + [self.hive_location]
    
    def create_path_using_nearest_flower(self, start_point):
        current_point = start_point
        visited = {start_point}
        path = [self.hive_location, current_point]

        while len(visited) < len(self.flowers):
            closest_flower = self.find_closest_point(current_point, self.flowers, visited)
            if closest_flower:
                path.append(closest_flower)
                visited.add(closest_flower)
                current_point = closest_flower

        path.append(self.hive_location)
        return path
        
    def find_closest_point(self, current_point, points, visited):
        min_distance = float('inf')
        closest_point = None
        for point in points:
            if point not in visited:
                distance = calculate_distance(current_point, point)
                if distance < min_distance:
                    min_distance = distance
                    closest_point = point
        return closest_point
    
    def select_pairs_for_crossover(self, selection_method="roulette_wheel", number_of_bees=50):
        selected_bees = sorted(self.population, key=lambda bee: bee.fitness)[:50]
        '''print('Bees selected for crosover:')
        for bee in selected_bees:
            if bee is None or bee.chromosome is None:
                print("Error: Found a bee or its chromosome as None!")
            else:
                print(f'{bee.unique_id}: {bee.fitness}; len: {len(bee.chromosome)}')'''
        pairs = []
        while len(pairs) < 25 and selected_bees:
            pair = []
            for _ in range(2):
                chosen_bee = self.roulette_wheel_selection(selected_bees)
                if chosen_bee:
                    pair.append(chosen_bee)
                    selected_bees.remove(chosen_bee)  # Remove chosen bee to avoid selecting it again
            if pair:
                pairs.append(tuple(pair))

        return pairs

    def roulette_wheel_selection(self, bees):
        # Invert fitness values to make lower fitness more likely to be selected
        inverted_fitness = [1.0 / bee.fitness for bee in bees]
        total_inverted_fitness = sum(inverted_fitness)

        r = random.uniform(0, total_inverted_fitness)
        P = 0
        for bee, inv_fit in zip(bees, inverted_fitness):
            P += inv_fit
            if P >= r:
                return bee
        return None
        
    def perform_crossover(self, pairs, crossover_method="partially_mapped", number_of_children=1, current_generation=0):
        offspring = []
        #print('parents:')
        
        for parent1, parent2 in pairs:
            if crossover_method == "partially_mapped":
                children = self.partially_mapped_crossover(parent1, parent2, number_of_children)
                if number_of_children == 1:
                    new_bee1 = self.create_bee(children, parent1.unique_id, parent2.unique_id)
                else:
                    new_bee1 = self.create_bee(children[0], parent1.unique_id, parent2.unique_id)
                    new_bee2 = self.create_bee(children[1], parent1.unique_id, parent2.unique_id)
                    offspring.append(new_bee2)
                    #self.ensure_hive_location(new_bee2.chromosome)
                    #print(new_bee2.chromosome)
                offspring.append(new_bee1)
                #print(new_bee1.chromosome)
                #self.ensure_hive_location(new_bee1.chromosome)
                '''if parent1.chromosome == parent2.chromosome:
                    print('identical parents')
                    if parent1.chromosome == children:
                        print('identical child')'''
            else:
                raise NotImplementedError(f"Crossover method '{crossover_method}' is not implemented.")
        # Sort the offspring by fitness and assign unique IDs
        offspring.sort(key=lambda bee: bee.fitness)
        
        for rank, bee in enumerate(offspring, start=1):
            bee.unique_id = f"{current_generation}-{rank}"
            '''print(bee.unique_id)
            print(bee.chromosome)
            print()'''
        #print('/parents')
        return offspring

    def partially_mapped_crossover(self, parent1, parent2, number_of_children=1):
        parent1_chromosome = parent1.chromosome
        parent2_chromosome = parent2.chromosome 
        '''self.ensure_hive_location(parent1_chromosome)
        self.ensure_hive_location(parent2_chromosome)'''
        '''print('parents chromosomes:')
        print(parent1_chromosome)
        print(parent2_chromosome)'''
        size = len(parent1_chromosome)
        # Step 1: Select crossover range at random
        start, end = sorted(random.sample(range(1, size - 2), 2))  # Avoid the first and last gene (the hive)
        
        # Step 2: Create offspring by exchanging the selected range
        child1 = parent1_chromosome[:start] + parent2_chromosome[start:end] + parent1_chromosome[end:]
        child2 = parent2_chromosome[:start] + parent1_chromosome[start:end] + parent2_chromosome[end:]

        # Step 3: Determine the mapping relationship to legalize offspring
        mapping1 = {parent2_chromosome[i]: parent1_chromosome[i] for i in range(start, end)}
        mapping2 = {parent1_chromosome[i]: parent2_chromosome[i] for i in range(start, end)}

        # Step 4: Legalize children with the mapping relationship
        for i in list(range(start)) + list(range(end, size)):
            if child1[i] in mapping1:
                while child1[i] in mapping1:
                    child1[i] = mapping1[child1[i]]
            if child2[i] in mapping2:
                while child2[i] in mapping2:
                    child2[i] = mapping2[child2[i]]

        '''if child1[0] != self.hive_location:
            print('0 is NOT THE HIVE')
        if child1[-1] != self.hive_location:
            print('Child 1 LAST is NOT THE HIVE')'''
        # Step 5: Return children based on the number of children requested
        if number_of_children == 1:
            return child1 if random.random() < 0.5 else child2
        else:
            '''if child2[0] != self.hive_location:
                print('0 is NOT THE HIVE')
            if child2[-1] != self.hive_location:
                print('Child 2 LAST is NOT THE HIVE')'''
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
            print(child.chromosome)
            print()'''
        return bees
    
    def mutate_bee_place_change(self, bee, mutation_rate):
        # Store the original chromosome and fitness
        self.ensure_hive_location(bee.chromosome)
        chromosome_before_mutation = bee.chromosome[:]
        fitness_before_mutation = bee.fitness

        # Create a copy for the mutated chromosome
        mutated_chromosome = chromosome_before_mutation[:]

        # Select a gene to mutate (excluding the hive locations)
        gene_index = random.randint(1, len(mutated_chromosome) - 2)

        # Determine the move direction and magnitude
        move_direction = random.choice([-1, 1])
        move_steps = min(mutation_rate, gene_index if move_direction == -1 else len(mutated_chromosome) - 2 - gene_index)

        # Calculate the new position, ensuring it stays within valid range
        new_position = (gene_index + move_direction * move_steps - 1) % (len(mutated_chromosome) - 2) + 1

        # Perform the mutation by moving the gene
        gene_to_move = mutated_chromosome[gene_index]
        mutated_chromosome.pop(gene_index)
        mutated_chromosome.insert(new_position, gene_to_move)

        # Update the bee's chromosome with the mutated one
        bee.chromosome = mutated_chromosome
        self.ensure_hive_location(mutated_chromosome)

        # Update fitness and revert if not better
        new_fitness = bee.calculate_fitness(bee.chromosome)
        if new_fitness < fitness_before_mutation:
            # If the new fitness is better, update fitness
            bee.fitness = new_fitness
            #print(f'{bee.unique_id} mutated is better')
        else:
            # Revert to the original chromosome and fitness
            #print(f'{bee.unique_id} not mutated is better')
            bee.chromosome = chromosome_before_mutation
            bee.fitness = fitness_before_mutation


    def update_population(self, offspring):
        """Update the population with fitter offspring."""
        combined_population = self.population + offspring
        combined_population.sort(key=lambda bee: bee.fitness)
        self.population = combined_population[:self.population_size]

    def analyze_fitness(self):
        """Calculate and analyze the fitness of the population."""
        total_fitness = sum(bee.fitness for bee in self.population)
        average_fitness = total_fitness / len(self.population)
        fittest_bee = min(self.population, key=lambda bee: bee.fitness)
        print(f"Average Fitness of Population: {average_fitness:.2f}")
        print(f"Fittest Bee ID: {fittest_bee.unique_id}, Fitness: {fittest_bee.fitness}")
        return {
            "average_fitness": average_fitness,
            "fittest_bee_id": fittest_bee.unique_id,
            "fittest_bee_chromosome": fittest_bee.chromosome,
            "fittest_bee_fitness": fittest_bee.fitness
        }
    
def calculate_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
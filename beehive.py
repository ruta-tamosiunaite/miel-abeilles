import sys
import math
import random


class Bee:
    def __init__(self, chromosome):
        self.unique_id = None
        self.chromosome = chromosome
        self.fitness = self.calculate_fitness()
        self.chromosome_before_mutation = None
        self.fitness_before_mutation = None
        self.parent1 = None
        self.parent2 = None
        self.normalized_fitness = None

    def calculate_fitness(self, chromosome=None):
        if chromosome is None:
            chromosome = self.chromosome

        total_distance = 0
        for i in range(len(chromosome) - 1):
            total_distance += calculate_distance(chromosome[i], chromosome[i + 1])
        return total_distance

    def mutate(self, mutation_rate, mutation_frequency, hive_location):
        best_fitness = self.fitness
        best_chromosome = self.chromosome[:]
        
        for _ in range(mutation_frequency):
            gene_index = random.randint(1, len(self.chromosome) - 2)
            move_direction = random.choice([-1, 1])
            move_steps = min(mutation_rate, gene_index if move_direction == -1 else len(self.chromosome) - 2 - gene_index)
            new_position = (gene_index + move_direction * move_steps - 1) % (len(self.chromosome) - 2) + 1

            # Perform mutation
            self.chromosome[gene_index], self.chromosome[new_position] = self.chromosome[new_position], self.chromosome[gene_index]
            self.ensure_hive_location(hive_location)
            self.fitness = self.calculate_fitness()

            # Check if this mutation improved fitness
            if self.fitness < best_fitness:
                best_fitness = self.fitness
                best_chromosome = self.chromosome[:]
            else:
                # Revert mutation if it didn't improve fitness
                self.chromosome = best_chromosome[:]
                self.fitness = best_fitness

    def ensure_hive_location(self, hive_location):
        if self.chromosome[0] != hive_location or self.chromosome[-1] != hive_location:
            print('Hive location is not correct !')

    @staticmethod
    def crossover(parent1, parent2, crossover_method="partially_mapped", number_of_children=1, current_generation=0):
        offspring = []
        if crossover_method == "partially_mapped":
            children = Bee.partially_mapped_crossover(parent1, parent2, number_of_children)
            if number_of_children == 1:
                offspring.append(Bee(children))
            else:
                offspring.extend([Bee(children[0]), Bee(children[1])])
        else:
            raise NotImplementedError(f"Crossover method '{crossover_method}' is not implemented.")
        return offspring

    @staticmethod
    def partially_mapped_crossover(parent1, parent2, number_of_children):
        parent1_chromosome = parent1.chromosome
        parent2_chromosome = parent2.chromosome 
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

        # Step 5: Return children based on the number of children requested
        if number_of_children == 1:
            return child1 if random.random() < 0.5 else child2
        else:
            return child1, child2
    
class Beehive:
    def __init__(self, flowers, hive_location, population_size):
        self.flowers = flowers
        self.hive_location = hive_location
        self.population_size = population_size
        self.population = []  # List of Bee objects

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

        # NEED TO DO THE SAME FOR ALL BEES...
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
        for parent1, parent2 in pairs:
            new_offspring = Bee.crossover(parent1, parent2, crossover_method, number_of_children, current_generation)
            offspring.extend(new_offspring)

        # Assign unique IDs to offspring and sort them
        for rank, bee in enumerate(offspring, start=1):
            bee.unique_id = f"{current_generation}-{rank}"
            bee.ensure_hive_location(self.hive_location)
        return offspring
    
    def perform_mutation(self, bees, mutation_method='place_change', mutation_rate=1, mutation_frequency=2):
        for bee in bees:
            if mutation_method == 'place_change':
                bee.mutate(mutation_rate, mutation_frequency, self.hive_location)
            else:
                raise NotImplementedError(f"Mutation method '{mutation_method}' is not implemented.")
            bee.ensure_hive_location(self.hive_location)
        bees.sort(key=lambda bee: bee.fitness)
        return bees

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
    
class BeesArchive:
    def __init__(self):
        self.all_bees = []  # List of all bee objects

    def update_archive(self, bees):
        for bee in bees:
            self.all_bees.append(bee)
    
def calculate_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
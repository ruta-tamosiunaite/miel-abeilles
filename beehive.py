import sys
import math
import random
import pandas as pd


class Bee:
    def __init__(self, chromosome, parent1_id='0-0', parent2_id='0-0'):
        self.unique_id = None
        self.chromosome = chromosome
        self.fitness = self.calculate_fitness()
        self.parent1_id = parent1_id
        self.parent2_id = parent2_id

    def calculate_fitness(self, chromosome=None):
        if chromosome is None:
            chromosome = self.chromosome

        total_distance = 0
        for i in range(len(chromosome) - 1):
            total_distance += calculate_distance(chromosome[i], chromosome[i + 1])
        return total_distance

    def mutate(self, mutation_rate, mutation_frequency):
        best_fitness = self.fitness
        best_chromosome = self.chromosome[:]
        
        for _ in range(mutation_frequency):
            gene_index = random.randint(1, len(self.chromosome) - 2)
            move_direction = random.choice([-1, 1])
            move_steps = min(mutation_rate, gene_index if move_direction == -1 else len(self.chromosome) - 2 - gene_index)
            new_position = (gene_index + move_direction * move_steps - 1) % (len(self.chromosome) - 2) + 1
    
            # Perform mutation
            self.chromosome[gene_index], self.chromosome[new_position] = self.chromosome[new_position], self.chromosome[gene_index]
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
            print('Error: Hive location is not correct !')
            sys.exit()
    
    @staticmethod
    def partially_mapped_crossover(parent1, parent2):
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

        # Step 5: Return the fitter child
        temp_bee1 = Bee(child1)
        temp_bee2 = Bee(child2)
        child1_fitness = temp_bee1.calculate_fitness()
        child2_fitness = temp_bee2.calculate_fitness()
        return child1 if child1_fitness < child2_fitness else child2 # Select the fitter child if only one child is needed

class Beehive:
    def __init__(self, flowers, hive_location, population_size, bees_archive, seed=None):
        self.flowers = flowers
        self.hive_location = hive_location
        self.population_size = population_size
        self.population = []  # List of Bee objects
        self.distances_table = self.calculate_distances_table()

        # Generate initial population
        if seed is None:
            seed = random.randint(0, 1000000)
        self.seed = seed
        random.seed(seed)  # Set the seed for reproducibility
        print(f"Random seed used for this run: {seed}")

        for i in range(self.population_size):
            chromosome = self.flowers[:]  # Create a copy of flowers
            random.shuffle(chromosome)  # Shuffle the copy
            new_bee = Bee(chromosome=[self.hive_location] + chromosome + [self.hive_location])
            self.population.append(new_bee)
            bees_archive.add_bee(new_bee)
    
    def calculate_distances_table(self):
        # Create a DataFrame to store distances
        locations = [(f'{self.hive_location}', self.hive_location)] + [(f'{flower}', flower) for i, flower in enumerate(self.flowers)]
        distances = {}

        for loc1_name, loc1_coords in locations:
            row_distances = {}
            for loc2_name, loc2_coords in locations:
                if loc1_name == loc2_name:
                    row_distances[loc2_name] = '-'
                else:
                    distance = math.sqrt((loc1_coords[0] - loc2_coords[0])**2 + (loc1_coords[1] - loc2_coords[1])**2)
                    row_distances[loc2_name] = round(distance, 2)
            distances[loc1_name] = row_distances

        distances_table = pd.DataFrame.from_dict(distances, orient='index')
        print(distances_table)
        return distances_table
    
    def run_generation(self, bees_archive, select_bees=50, mutation_rate=1, mutation_frequency=1, current_generation=0):
        # S E L E C T  P A I R S  F O R  C R O S S O V E R
        selected_bees = sorted(self.population, key=lambda bee: bee.fitness)[:select_bees]
        pairs = []
        while len(pairs) < (select_bees / 2) and selected_bees:
            pair = []
            for _ in range(2):
                chosen_bee = self.roulette_wheel_selection(selected_bees)
                if chosen_bee:
                    pair.append(chosen_bee)
                    selected_bees.remove(chosen_bee)  # Remove chosen bee to avoid selecting it again
            if pair:
                pairs.append(tuple(pair))

        # P E R F O R M  C R O S S O V E R
        offspring = []
        for parent1, parent2 in pairs:
            new_offspring = Bee.partially_mapped_crossover(parent1, parent2)
            offspring.append(Bee(new_offspring, parent1_id=parent1.unique_id, parent2_id=parent2.unique_id))

        # P E R F O R M  M U T A T I O N
        average_fitness = sum(bee.fitness for bee in self.population) / len(self.population)
        for bee in offspring:
            bee.mutate(mutation_rate, mutation_frequency)
            if bee.fitness > average_fitness:
                pass
                #bee.mutate(mutation_rate, self.distances_table)
            # If the new bee fitness is bigger than the average fitness of the population, then launch mutation.
            # mutation rate 1 will mean that a random gene will be selected. 
            # The two closest points to the gene will be taken from the already created table of all distances. 
            # These fllowers will be moved to be from the left and right of the randomly selected gene (trying both and choosing the option where bee fitness is better - smaller). 
            # If the hive happens to be the closest point, then the not the hive but the randomly selected gene moves to the beginning or the end of the chromosome, and the second closest point to the gene go before of after the randomly selected gene, but always keeping the hive the first and the last element of the chromosome (trying both - putting gene as a first or one before last element of the chrosomosome and choosing the option where bee fitness is better) .
        
        # Sort offspring and assign unique IDs to them
        offspring.sort(key=lambda bee: bee.fitness)
        for rank, bee in enumerate(offspring, start=1):
            bee.unique_id = f"{current_generation}-{rank}"
            bee.ensure_hive_location(self.hive_location)

        # U P D A T E  P O P U L A T I O N
        combined_population = self.population + offspring
        combined_population.sort(key=lambda bee: bee.fitness)
        self.population = combined_population[:self.population_size]
        for bee in combined_population: # Add new bees to the archive
            if bee in offspring:
                bees_archive.add_bee(bee)

        # A N A L Y Z E  F I T N E S S
        results = self.analyze_fitness()
        return results
    
    def roulette_wheel_selection(self, bees):
        # Invert fitness values to make lower fitness more likely to be selected
        inverted_fitness = [1000 / bee.fitness for bee in bees]
        total_inverted_fitness = sum(inverted_fitness)

        r = random.uniform(0, total_inverted_fitness)
        P = 0
        for bee, inv_fit in zip(bees, inverted_fitness):
            P += inv_fit
            if P >= r:
                return bee
        return None

    def analyze_fitness(self):
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
        self.all_bees = []

    def add_bee(self, bee):
        self.all_bees.append(bee)

    def get_bee_by_id(self, unique_id):
        for bee in self.all_bees:
            if bee.unique_id == unique_id:
                return bee
        return None
    
def calculate_distance(point1, point2):
    return math.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)

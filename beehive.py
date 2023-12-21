import sys
import math
import random


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

    def mutate(self, distances_table, mutation_rate):
        beehive = self.chromosome[0]
        available_genes = self.chromosome[:]
        copy_of_the_chromosome = self.chromosome[:]
        chromosome_length = len(self.chromosome)
        for _ in range(min(mutation_rate, len(available_genes))):
            point_selected = random.choice(available_genes)
            available_genes.remove(point_selected)
            random_gene = self.chromosome.index(point_selected)

            closest_points = distances_table[point_selected] # closest points in desc order for the random point selected
            closest_point_1 = closest_points[0]
            closest_point_2 = closest_points[1]

            neighbour_from_left = self.chromosome[random_gene-1] # Find neighbours
            try: neighbour_from_right = self.chromosome[random_gene+1]
            except: neighbour_from_right = beehive # In case when beehive selected

            # Neighbour to the left  -> closest_point_2
            # Neighbour to the right -> closest_point_1
            # left_or_right = random.randint(0, 1)
            if True:
                if closest_point_2 == beehive:
                    self.chromosome.remove(point_selected)
                    self.chromosome.insert(1, point_selected)
                    neighbour_from_left = beehive
                    neighbour_from_right = self.chromosome[2]

                if neighbour_from_left != closest_point_2:
                    self.chromosome.remove(closest_point_2)
                    index_of_closest_point_2 = self.chromosome.index(point_selected)
                    if index_of_closest_point_2 == 0:
                        self.chromosome.insert(chromosome_length-2, closest_point_2)
                    else:
                        self.chromosome.insert(index_of_closest_point_2, closest_point_2)
            if True:
                if closest_point_1 == beehive:
                    self.chromosome.remove(point_selected)
                    self.chromosome.insert(chromosome_length-2, point_selected)
                    neighbour_from_left = self.chromosome[chromosome_length-3]
                    neighbour_from_right = beehive

                if neighbour_from_right != closest_point_1:
                    self.chromosome.remove(closest_point_1)
                    index_of_closest_point_1 = self.chromosome.index(point_selected) + 1
                    self.chromosome.insert(index_of_closest_point_1, closest_point_1)
            previous_fitness = self.calculate_fitness(copy_of_the_chromosome)
            fitness_after_mutation = self.calculate_fitness(self.chromosome)
            if previous_fitness < fitness_after_mutation:
                self.chromosome = copy_of_the_chromosome
                self.fitness = previous_fitness
        return self.chromosome

    def ensure_hive_location(self, hive_location):
        if self.chromosome[0] != hive_location or self.chromosome[-1] != hive_location:
            print('Error: Hive location is not correct !')
            print(self.chromosome)
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
        self.distances_table = self.calculate_closest_points()

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
    
    def calculate_closest_points(self):
        points = self.flowers + [self.hive_location]
        closest_points = {} # Dictionary to hold the closest points for each location

        for i, point in enumerate(points):
            distances = []
            for j, other_point in enumerate(points):
                if i != j:
                    distance = calculate_distance(point, other_point)
                    distances.append((distance, other_point))

            distances.sort() # Sort the distances and get the four closest points
            closest_points[point] = [location for _, location in distances[:4]]

        return closest_points
    
    def run_generation(self, bees_archive, select_bees=50, mutation_rate=1, current_generation=0):
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
            bee.ensure_hive_location(self.hive_location)
            if bee.fitness > average_fitness:
                bee.mutate(self.distances_table, mutation_rate)
        
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

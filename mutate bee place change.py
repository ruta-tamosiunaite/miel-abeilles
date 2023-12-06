import random
import math

def mutate_bee_place_change(bee, mutation_rate):
    # Store the original chromosome and fitness
    bee_chromosome_before_mutation = bee[:]

    # Select a gene to mutate (excluding the hive locations)
    gene_index = random.randint(1, len(bee_chromosome_before_mutation) - 2)
    print(f'gene_index: {gene_index}')

    move_direction = random.choice([-1, 1])
    print(f'move_direction: {move_direction}')

    move_steps = min(mutation_rate, gene_index if move_direction == -1 else len(bee_chromosome_before_mutation) - 1 - gene_index)
    print(f'move_steps: {move_steps}')

    #new_position = gene_index + move_direction * move_steps
    new_position = (gene_index + move_direction * move_steps - 1) % len(bee_chromosome_before_mutation) + 1
    if new_position < 1:
        new_position = 1
    elif new_position >= len(bee) - 1:
        new_position = len(bee) - 2
    print(f'new_position: {new_position}')

    gene_to_move = bee_chromosome_before_mutation[gene_index]
    print(f'gene_to_move: {gene_to_move}')

    bee.pop(gene_index)
    bee.insert(new_position, gene_to_move)
    print(f'bee: {bee}')

    return bee

#bee = [(500, 500), (730, 742), (898, 898), (921, 964), (573, 903), (837, 787), (549, 329), (501, 287), (951, 209), (98, 711), (276, 666), (444, 428), (437, 601), (704, 995), (494, 898), (761, 772), (938, 646), (684, 273), (796, 310), (774, 130), (964, 726), (602, 68), (220, 307), (628, 311), (897, 324), (724, 631), (278, 434), (908, 534), (371, 429), (708, 99), (508, 31), (643, 404), (929, 389), (189, 739), (345, 560), (101, 482), (213, 639), (273, 105), (45, 549), (116, 69), (758, 405), (3, 517), (330, 410), (2, 897), (328, 489), (875, 407), (862, 361), (493, 970), (528, 794), (436, 619), (901, 639), (500, 500)]
bee = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0]

print(f'bee: {bee}')

mutation = 5

print(f'mutation: {mutation}')
mutate_bee_place_change(bee, mutation)
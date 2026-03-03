import numpy as np

def initialize_population(pop_size, chromosome_length):
    return np.random.randint(0, 2, size=(pop_size, chromosome_length))

def fitness(population):
    return population.sum(axis=1)

def roulette_wheel_selection(population, fitness_scores):
    total_fitness = fitness_scores.sum()
    if total_fitness == 0:
        idx = np.random.randint(len(population), size=2)
        return population[idx[0]], population[idx[1]]

    selection_probs = fitness_scores / total_fitness
    indices = np.random.choice(len(population), size=2, p=selection_probs)
    return population[indices[0]], population[indices[1]]

def single_point_crossover(parent1, parent2, crossover_rate):
    if np.random.rand() < crossover_rate:
        point = np.random.randint(1, len(parent1))
        child1 = np.concatenate((parent1[:point], parent2[point:]))
        child2 = np.concatenate((parent2[:point], parent1[point:]))
        return child1, child2
    return parent1.copy(), parent2.copy()

def mutate(chromosome, mutation_rate):
    mutation_mask = np.random.rand(len(chromosome)) < mutation_rate
    chromosome[mutation_mask] ^= 1
    return chromosome

def run_ga(crossover_rate, runs=20, pop_size=100, chromosome_length=100, mutation_rate=0.001):
    generations_to_solution = []

    for run in range(runs):
        population = initialize_population(pop_size, chromosome_length)
        generation = 0

        while True:
            fit = fitness(population)

            if np.any(fit == chromosome_length):
                print(f"Run {run+1}: Solution found in {generation} generations")
                generations_to_solution.append(generation)
                break

            new_population = []

            while len(new_population) < pop_size:
                parent1, parent2 = roulette_wheel_selection(population, fit)
                child1, child2 = single_point_crossover(parent1, parent2, crossover_rate)
                child1 = mutate(child1, mutation_rate)
                child2 = mutate(child2, mutation_rate)
                new_population.extend([child1, child2])

            population = np.array(new_population[:pop_size])
            generation += 1

    return generations_to_solution


# ---------------- RUN AND PRINT RESULTS ----------------

results = run_ga(crossover_rate=0.7, runs=20)

print("\nFinal Results:")
print("Generations per run:", results)
print("Average generations:", np.mean(results))
print("Best run:", np.min(results))
print("Worst run:", np.max(results))
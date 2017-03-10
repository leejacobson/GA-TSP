package tsp;

import java.util.Arrays;

public class GeneticAlgorithm {
    private int populationSize;
    private double mutationRate;
    private double crossoverRate;
    private int elitismCount;
    private int tournamentSize;
    
    /**
     * Initalize genetic algorithm
     * 
     * @param populationSize The size of the GA population
     * @param mutationRate GA mutation rate
     * @param crossoverRate GA crossover rate
     * @param elitismCount Number of elite individuals in population
     * @param tournamentSize The number of individual per tournament
     */
    public GeneticAlgorithm(int populationSize, double mutationRate, 
            double crossoverRate, int elitismCount, int tournamentSize){
        this.populationSize = populationSize;
        this.mutationRate = mutationRate;
        this.crossoverRate = crossoverRate;
        this.elitismCount = elitismCount;
        this.tournamentSize = tournamentSize;
    }
    
    /**
     * Initialize population
     * 
     * @param chromosomeLength The length of the individuals chromosome
     * @return population The initial population generated
     */
    public Population initPopulation(int chromosomeLength){
        // Initialize population
        Population population = new Population(this.populationSize, chromosomeLength);
        return population;
    }
    
    /**
     * Calculate individuals fitness value
     * 
     * @param individual the individual to evaluate
     * @param cities the cities being referenced
     * @return double The fitness value for individual
     */
    public double calcFitness(Individual individual, City cities[]){        
        // Get fitness
        Route route = new Route(individual, cities);
        double fitness = 1 / route.getDistance();
                
        // Store fitness
        individual.setFitness(fitness);
        
        return fitness;
    }
    
    /**
     * Evaluate population
     * 
     * @param population the population to evaluate
     * @param cities the cities being referenced
     */
    public void evalPopulation(Population population, City cities[]){
        double totalPopulationFitness = 0;
        
        // Loop over population evaluating individuals and suming population fitness
        for (Individual individual : population.getIndividuals()) {
            totalPopulationFitness += this.calcFitness(individual, cities);
        }
        
        double avgFitness = totalPopulationFitness / population.size();
        population.setAvgFitness(avgFitness);
    }
    
    /**
     * Check if population has met termination condition
     * 
     * @param generationsCount Number of generations passed
     * @param maxGenerations Number of generations to terminate after
     * @return boolean True if termination condition met, otherwise, false
     */
    public boolean isTerminationConditionMet(int generationsCount, int maxGenerations){     
        return (generationsCount > maxGenerations);
    }
    
    /** 
     * Selects parent for crossover using tournament selection
     * 
     * @param population
     * @return The individual selected as a parent
     */
    public Individual tournamentSelection(Population population){
        // Create tournament
        Population tournament = new Population(this.tournamentSize);

        // Add random individuals to the tournament
        population.shuffle();
        for (int i=0; i < this.tournamentSize; i++) {
            Individual tournamentIndividual = population.getIndividual(i);
            tournament.setIndividual(tournamentIndividual, i);
        }
        
        // Return the best
        return tournament.getFittest(0);
    }
    
    /**
     * Select parent for crossover
     * 
     * @param population The population to select parent from
     * @return The individual selected as a parent
     */
    public Individual selectParent(Population population) {
        // Get individuals
        Individual individuals[] = population.getIndividuals();
       
        // Spin roulette wheel
        double populationFitness = population.getPopulationFitness();
        double rouletteWheelPosition = Math.random() * populationFitness;

        // Find parent
        double spinWheel = 0;
        for (Individual individual : individuals) {
            spinWheel += individual.getFitness();
            if (spinWheel >= rouletteWheelPosition) {
                return individual;
            }
        }

        return individuals[population.size()-1];
    }
    
    /**
     * Crossover population using single point crossover
     * 
     * @param population Population to crossover
     * @return Population The new population
     */
    public Population singlePointCrossover(Population population){
        // Create new population
        Population newPopulation = new Population(population.size());
        
        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);
            
            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex > this.elitismCount) {
                // Initialize offspring
                int offspringChromosome[] = new int[parent1.getChromosomeLength()];
                
                 // Find second parent
                Individual parent2 = this.tournamentSelection(population);
                
                // Get random swap point
                int swapPoint = (int) (Math.random() * (parent1.getChromosomeLength()+1));

                // Loop over genome
                for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    // Use half of parent1's genes and half of parent2's genes 
                    if (geneIndex < swapPoint) {
                        offspringChromosome[geneIndex] = parent1.getGene(geneIndex);
                    } else {
                        offspringChromosome[geneIndex] = parent2.getGene(geneIndex);
                    }
                }

                // Add offspring to new population
                newPopulation.setIndividual(new Individual(offspringChromosome), populationIndex);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(parent1, populationIndex);
            }
        }
        
        return newPopulation;
    }
    
    /**
     * Apply crossover to population
     * 
     * @param population The population to apply crossover to
     * @return The new population
     */
    public Population crossoverPopulation(Population population){
        // Create new population
        Population newPopulation = new Population(population.size());
        
        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual parent1 = population.getFittest(populationIndex);
            
            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex > this.elitismCount) {
                // Initialize offspring
                int offspringChromosome[] = new int[parent1.getChromosomeLength()];
                
                 // Find second parent
                Individual parent2 = selectParent(population);

                // Loop over genome
                for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
                    // Use half of parent1's genes and half of parent2's genes 
                    if (0.5 > Math.random()) {
                        offspringChromosome[geneIndex] = parent1.getGene(geneIndex);
                    } else {
                        offspringChromosome[geneIndex] = parent2.getGene(geneIndex);
                    }
                }

                // Add offspring to new population
                newPopulation.setIndividual(new Individual(offspringChromosome), populationIndex);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(parent1, populationIndex);
            }
        }
        
        return newPopulation;
    }
    
    /**
     * Apply mutation to population
     * 
     * @param population The population to apply mutation to
     * @return The mutated population
     */
    public Population mutatePopulation(Population population){
        // Initialize new population
        Population newPopulation = new Population(this.populationSize);
        
        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            Individual individual = population.getFittest(populationIndex);

            // Skip mutation if this is an elite individual
            if (populationIndex > this.elitismCount) {   
                // Loop over individual's genes
                for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {            
                    // Does this gene need mutation?
                    if (this.mutationRate > Math.random()) {
                        // Get new gene position
                        int newGenePos = (int) (Math.random() * individual.getChromosomeLength());
                        // Get genes to swap
                        int gene1 = individual.getGene(newGenePos);
                        int gene2 = individual.getGene(geneIndex);
                        // Swap genes
                        individual.setGene(gene1, geneIndex);
                        individual.setGene(gene2, newGenePos);
                    }
                }
            }
            
            // Add individual to population
            newPopulation.setIndividual(individual, populationIndex);
        }
        
        // Return mutated population
        return newPopulation;
    }
    
    /**
     * Order crossover mutation to population
     * 
     * @param population The population to apply crossover to
     * @return The new population
     */
    public Population orderedCrossover(Population population){
        // Create new population
        Population newPopulation = new Population(population.size());
        
        // Loop over current population by fitness
        for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
            // Get parent1
            Individual parent1 = population.getFittest(populationIndex);
            
            // Apply crossover to this individual?
            if (this.crossoverRate > Math.random() && populationIndex > this.elitismCount) {
                // Find parent2 with tournament selection
                Individual parent2 = this.tournamentSelection(population);

                // Create blank offspring chromosome
                int offspringChromosome[] = new int[parent1.getChromosomeLength()];
                Arrays.fill(offspringChromosome, -1);
                Individual offspring = new Individual(offspringChromosome);

                // Get subset of parent chromosomes
                int substrPos1 = (int) (Math.random() * parent1.getChromosomeLength());
                int substrPos2 = (int) (Math.random() * parent1.getChromosomeLength());

                // make the smaller the start and the larger the end
                final int startSubstr = Math.min(substrPos1, substrPos2);
                final int endSubstr = Math.max(substrPos1, substrPos2);

                // Loop and add the sub tour from parent1 to our child
                for (int i = startSubstr; i < endSubstr; i++) {
                    offspring.setGene(parent1.getGene(i), i);
                }

                // Loop through parent2's city tour
                for (int i = 0; i < parent2.getChromosomeLength(); i++) {
                    int parent2Gene = i + endSubstr;
                    if (parent2Gene >= parent2.getChromosomeLength()) {
                        parent2Gene -= parent2.getChromosomeLength();
                    }

                    // If offspring doesn't have the city add it
                    if (offspring.containsGene(parent2.getGene(parent2Gene)) == false) {
                        // Loop to find a spare position in the child's tour
                        for (int ii = 0; ii < offspring.getChromosomeLength(); ii++) {
                            // Spare position found, add city
                            if (offspring.getGene(ii) == -1) {
                                offspring.setGene(parent2.getGene(parent2Gene), ii);
                                break;
                            }
                        }
                    }
                }

                // Add child
                newPopulation.setIndividual(offspring, populationIndex);
            } else {
                // Add individual to new population without applying crossover
                newPopulation.setIndividual(parent1, populationIndex);
            }
        }
        
        return newPopulation;
    }
    
}

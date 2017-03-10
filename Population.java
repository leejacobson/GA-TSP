package tsp;

import java.util.Arrays;
import java.util.Comparator;

public class Population {
    private Individual population[];
    private double avgFitness = -1;
    
    /**
     * Initializes blank population of individuals
     * 
     * @param populationSize The size of the population
     */
    public Population(int populationSize){
        // Initial population
        this.population = new Individual[populationSize];
    }
    
    /**
     * Initializes population of individuals
     * 
     * @param populationSize The size of the population
     * @param chromosomeLength The length of the individuals chromosome
     */
    public Population(int populationSize, int chromosomeLength){
        // Initial population
        this.population = new Individual[populationSize];
        
        // Loop over population size
        for (int individualCount = 0; individualCount < populationSize; individualCount++) {
            // Create individual
            Individual individual = new Individual(chromosomeLength);
            // Add individual to population
            this.population[individualCount] = individual;
        }
    }
    
    /**
     * Get individuals from the population
     * 
     * @return individuals Individuals in population
     */
    public Individual[] getIndividuals(){
        return this.population;
    }
    
    /**
     * Find fittest individual in the population
     * 
     * @param offset
     * @return individual Fittest individual at offset
     */
    public Individual getFittest(int offset){ 
        // Order population by fitness
        Arrays.sort(this.population, new Comparator<Individual>() {
            @Override
            public int compare(Individual o1, Individual o2) {
                if (o1.getFitness() > o2.getFitness()) {
                    return -1;
                }
                else if (o1.getFitness() < o2.getFitness()) {
                    return 1;
                }
                return 0;
            }
        });

        // Return the fittest individual
        return this.population[offset];
    }
    
    /**
     * Randomly order individuals in the population
     */
    public void shuffle(){
        // Run through each individual
        for (int individualIndex=0; this.size() > individualIndex; individualIndex++) {
          // Get new, random, index for individual
          int index = (int) (Math.random() * this.size());
          // Swap individuals
          Individual tempIndividual = this.population[index];
          this.population[index] = this.population[individualIndex];
          this.population[individualIndex] = tempIndividual;
        }
    }
    
    /**
     * Sets the average population's fitness
     * 
     * @param fitness The population's average fitness
     */
    public void setAvgFitness(double fitness){
        this.avgFitness = fitness;
    }
    
    /**
     * Get the population's average fitness
     * 
     * @return avgFitness The population's average fitness
     */
    public double getPopulationFitness(){
        return this.avgFitness;
    }
    
    /**
     * Get population's size
     * 
     * @return size The population's size
     */
    public int size(){
        return this.population.length;
    }
    
    /**
     * Set individual at offset
     * 
     * @param individual
     * @param offset
     * @return individual
     */
    public Individual setIndividual(Individual individual, int offset){
        return population[offset] = individual;
    }
    
    /**
     * Get individual at offset
     * 
     * @param offset
     * @return individual
     */
    public Individual getIndividual(int offset){
        return population[offset];
    }
}
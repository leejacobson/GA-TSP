package tsp;

public class Individual {
    private int chromosome[];
    private double fitness = -1;
      
    /**
     * Initializes individual with specific chromosome
     * 
     * @param chromosome The chromosome to give individual
     */
    public Individual(int chromosome[]){
        // Create individual chromosome
        this.chromosome = chromosome;
    }
    
    /**
     * Initializes random individual
     * 
     * @param chromosomeLength The length of the individuals chromosome
     */
    public Individual(int chromosomeLength){
        // Create random individual
        int chromosome[] = new int[chromosomeLength];
        for (int gene = 0; gene < chromosomeLength; gene++) {
            chromosome[gene] = gene;
        }
        
        for (int geneIndex=0; geneIndex < chromosomeLength; geneIndex++) {
          // Shuffle genes
          int newGeneIndex = (int) (Math.random() * chromosomeLength);
          // Swap individuals
          int tempIndividual = chromosome[newGeneIndex];
          chromosome[newGeneIndex] = chromosome[geneIndex];
          chromosome[geneIndex] = tempIndividual;
        }
        
        this.chromosome = chromosome;
    }
    
    /**
     * Gets individual's chromosome
     * 
     * @return The individual's chromosome
     */
    public int[] getChromosome(){
        return this.chromosome;
    }
    
    /**
     * Gets individual's chromosome length
     * 
     * @return The individual's chromosome length
     */
    public int getChromosomeLength(){
        return this.chromosome.length;
    }
    
    /**
     * Check if gene is contained in this individual
     * 
     * @param gene Gene to look for
     * @return boolean
     */
    public boolean containsGene(int gene){
        // Loop over genes
        for (int i=0; i < this.getChromosomeLength(); i++) {
            // Check for gene
            if (this.getGene(i) == gene) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Set gene at offset
     * 
     * @param gene
     * @param offset
     * @return gene
     */
    public void setGene(int gene, int offset){
        this.chromosome[offset] = gene;
    }
    
    /**
     * Get gene at offset
     * 
     * @param offset
     * @return gene
     */
    public int getGene(int offset){
        return this.chromosome[offset];
    }
    
    /**
     * Store individual's fitness
     * 
     * @param fitness The individuals fitness
     */
    public void setFitness(double fitness){
        this.fitness = fitness;
    }
    
    /**
     * Gets individual's fitness
     * 
     * @return The individual's fitness
     */
    public double getFitness(){
        return this.fitness;
    }
}

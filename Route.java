package tsp;

public class Route {
    private City route[];
    private double distance = 0;
    
    /**
     * Initialize Route
     * 
     * @param individual A GA individual
     * @param cities The cities referenced
     */
    public Route(Individual individual, City cities[]){
        // Get individual's chromosome
        int chromosome[] = individual.getChromosome();
        // Create route
        this.route = new City[cities.length];
        for(int geneIndex=0; geneIndex < chromosome.length; geneIndex++){
            this.route[geneIndex] = cities[chromosome[geneIndex]];
        }
    }
    
    /**
     * Get route distance
     * 
     * @return distance The route's distance
     */
    public double getDistance(){
        if (this.distance > 0) {
            return this.distance;
        }
        
        // Loop over cities in route and calculate route distance
        double totalDistance = 0;
        for(int cityIndex=0; cityIndex+1<this.route.length; cityIndex++){
            totalDistance += this.route[cityIndex].distanceFrom(this.route[cityIndex+1]);
        }
        
        totalDistance += this.route[this.route.length-1].distanceFrom(this.route[0]);
        this.distance = totalDistance;
        
        return totalDistance;
    }
    
    /**
     * Print route
     */
    public String toString(){
        String route = "";
        for (City route1 : this.route) {
            route += "(" + route1.getX() + ", " + route1.getY() + ")";
        }
        return route;
    }
}

//
//  NEATPlusAlgorithm.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__NEATPlusAlgorithm__
#define __SystemGenerator__NEATPlusAlgorithm__

#include <vector>
#include <fstream>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>

template <class IndividualType>
struct SystemInfo {
    IndividualType individual;
    double fitness;
    
    SystemInfo(IndividualType i):individual(i) { fitness = -INFINITY;}
};

template <class IndividualType>
class NEATSpecies {
public:
// each generation the members are cleared out, repopulated, and the representative is updated
    IndividualType *representative;
    std::vector<SystemInfo<IndividualType> >members;
    double totalSharedFitness;
    int speciesNumber; // for data collection
};

template <class IndividualType, class InnovationType>
class NEATPlusAlgorithm {
protected:
    std::vector<InnovationType> newConnections;
    std::vector<SystemInfo<IndividualType> > population;
    std::vector<NEATSpecies<IndividualType> >speciesList;

    int populationSize;
    int stagnantGenerations;
    int maxStagnation;
    int maxGenerations;
    
    long generations;
    IndividualType *_bestIndividual;
    double allTimeBestFitness;
    double lastBestFitness; // before sharing
 
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one

    void mutateSystem(IndividualType& original); // does not make a copy

    virtual IndividualType combineSystems(IndividualType &sys1, IndividualType &sys2); // assumes the fitter individual is first
    virtual double genomeDistance( IndividualType& sys1,  IndividualType& sys2);
    virtual std::pair<SystemInfo<IndividualType>, SystemInfo<IndividualType> > selectParents(double fitnessSum);
    virtual void assignInnovationNumberToAttachment(IndividualType& individual, InnovationType i);

    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
        
    void logPopulationStatistics();
    NEATSpecies<IndividualType> chooseBreedingSpecies(double totalFitness);
public:
    NEATPlusAlgorithm(int populationSize, int maxGenerations, int maxStagnation);
    ~NEATPlusAlgorithm();
    
    // notice that this cannot be overriden.
    bool tick();
#pragma mark - Function pointers - MUST BE SET OR ELSE
    // please set these
    double (*evaluationFunc)(IndividualType& individual); // provide a fitness for each individual
    bool (*stopFunc)(double bestFitness); // set a stopping fitness
    IndividualType (*createInitialIndividual)(void); // the smallest/starting individual
#pragma mark - 
    
    IndividualType *bestIndividual();
    long getNumberOfIterations();
    
#pragma mark - Tuning parameters
    // probability of recombination taking from weaker parent
    double p_c = 0.3;
    
    // mutation probabilities
    double p_m_node_ins = 0.1; // how likely we are to add a new node
    double p_m_conn_ins = 0.2;   // how likely we are to add a new connection
    
    double p_m_node = 0.3; // how likely we are to mutate each existing node
    double p_m_conn = 0.3; // how likely are we to mutate each existing connection
    
    double d_threshold = 3.0; // threshold distance - a difference bigger than this will exclude an individual from a species
    
    double w_excess = 0.25; // how heavily the number of excess genes are weighed
    double w_disjoint = 0.25; // how heavily the number of disjoint genes are weighed
    double w_matching = 0.5; // how heavily the total difference between matching genes is weighed
    
    double w_matching_node = 0.5; // how heavily we weigh the difference between matching nodes
#pragma mark -
    
private:
    std::ofstream currentLogFile;
    IndividualType *origin; // to find distance of species from start
    int nextInnovationNumber;
    int nextSpeciesNumber;
    void speciate(); // divide everything into species
};


#endif /* defined(__SystemGenerator__NEATPlusAlgorithm__) */

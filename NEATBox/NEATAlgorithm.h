//
//  NEATAlgorithm.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__NEATAlgorithm__
#define __SystemGenerator__NEATAlgorithm__

#include <vector>
#include <fstream>
#include "NEATTypes.h"

template <class IndividualType>
struct SystemInfo {
    IndividualType individual;
    double fitness;
    
    SystemInfo(IndividualType i)
    :individual(i) { fitness = -INFINITY;}
};

template <class IndividualType>
class NEATSpecies {
// each generation the members are cleared out, repopulated, and the representative is updated
    IndividualType representative;
    std::vector<SystemInfo<IndividualType> *>members;
    double totalSharedFitness;
    int speciesNumber; // for data collection
};

template <class IndividualType, class InnovationType>
class NEATAlgorithm {
protected:
    std::vector<InnovationType> newConnections;
    std::vector<SystemInfo<IndividualType> > population;
    std::vector<NEATSpecies<IndividualType> >speciesList;

    // recombination probability
    double p_c;
    
    // mutation probabilities
    double p_m_node; // how likely we are to add a new node
    double p_m_conn;   // how likely we are to add a new connection
    double p_m_attach; // how likely are we to mutate each existing attachment
    
    // threshold distance - a difference bigger than this will exclude an individual from a species
    
    double d_threshold = 1.0;

    double w_excess = 0.5;
    double w_disjoint = 0.25;
    double w_matching = 0.25;

    int populationSize;
    int stagnantGenerations;
    int maxStagnation;
    int maxGenerations;
    
    long generations;
    IndividualType _bestIndividual;
    double allTimeBestFitness;
    double lastBestFitness; // before sharing
 
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one

    void mutateSystem(IndividualType original); // does not make a copy

#pragma mark - Override these functions in NEAT variants only
    virtual IndividualType combineSystems(IndividualType sys1, IndividualType sys2); // assumes the fitter individual is first
    virtual double genomeDistance(const IndividualType& sys1, const IndividualType& sys2);
    virtual void selectParents(SystemInfo<IndividualType> *parent1, SystemInfo<IndividualType> *parent2, double fitnessSum);
    virtual void assignInnovationNumberToAttachment(IndividualType individual, InnovationType i);

#pragma mark - Override these functions any time you need
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
        
    void logPopulationStatistics();
public:
    NEATAlgorithm(int populationSize, int maxGenerations, int maxStagnation, float p_c=0.5, float p_m_attach=0.1, float p_m_node=0.1, float p_m_conn=0.1);
    ~NEATAlgorithm();
    
    // notice that this cannot be overriden.
    bool tick();
    
    // please set these
    double (*evaluationFunc)(IndividualType& individual); // provide a fitness for each individual
    bool (*stopFunc)(double bestFitness); // set a stopping fitness
    
    const IndividualType& bestIndividual();
    long getNumberOfIterations();
    
    
private:
    std::ofstream currentLogFile;
    IndividualType origin; // to find distance of species from start
    int nextInnovationNumber;
    int nextSpeciesNumber;
    void speciate(); // divide everything into species
};

#endif /* defined(__SystemGenerator__NEATAlgorithm__) */

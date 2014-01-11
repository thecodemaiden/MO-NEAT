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
class NEATAlgorithm {
protected:
    std::vector<InnovationType> newConnections;
    std::vector<SystemInfo<IndividualType> > population;
    std::vector<NEATSpecies<IndividualType> >speciesList;

    // probability of recombination taking from weaker parent
    double p_c = 0.3;
    
    // mutation probabilities
    double p_m_node_ins = 0.1; // how likely we are to add a new node
    double p_m_conn_ins = 0.2;   // how likely we are to add a new connection
    
    double p_m_node = 0.3; // how likely we are to mutate each existing node
    double p_m_conn = 0.3; // how likely are we to mutate each existing connection
    
    // threshold distance - a difference bigger than this will exclude an individual from a species
    
    double d_threshold = 1.0;

    double w_excess = 0.25;
    double w_disjoint = 0.25;
    double w_matching = 0.5;

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

#pragma mark - Override these functions in NEAT variants only
    virtual IndividualType combineSystems(IndividualType &sys1, IndividualType &sys2); // assumes the fitter individual is first
    virtual double genomeDistance( IndividualType& sys1,  IndividualType& sys2);
    virtual std::pair<SystemInfo<IndividualType>, SystemInfo<IndividualType> > selectParents(double fitnessSum);
    virtual void assignInnovationNumberToAttachment(IndividualType& individual, InnovationType i);

#pragma mark -
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
        
    void logPopulationStatistics();
public:
    NEATAlgorithm(int populationSize, int maxGenerations, int maxStagnation);
    ~NEATAlgorithm();
    
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
    
    
private:
    std::ofstream currentLogFile;
    IndividualType *origin; // to find distance of species from start
    int nextInnovationNumber;
    int nextSpeciesNumber;
    void speciate(); // divide everything into species
};


#endif /* defined(__SystemGenerator__NEATAlgorithm__) */

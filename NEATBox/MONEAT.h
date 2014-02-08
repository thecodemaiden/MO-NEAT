//
//  MONEAT.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__MONEAT__
#define __SystemGenerator__MONEAT__

#include <vector>
#include <fstream>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>

#include "MNEdge.h"
#include "MNIndividual.h"

struct InnovationInfo;

struct SystemInfo {
    MNIndividual *individual;
    std::vector<double> fitnesses;
    double rankFitness;
    int dominationCount;
    
    SystemInfo(MNIndividual *i):individual(i) {rankFitness = -INFINITY; dominationCount=0;}
    SystemInfo(const SystemInfo &other)
     :individual(other.individual->clone()),
      fitnesses(other.fitnesses),
      rankFitness(other.rankFitness),
      dominationCount(0)
    {}
    
    ~SystemInfo() {
        delete individual;
    }
};

class NEATSpecies {
public:
// each generation the members are cleared out, repopulated, and the representative is updated
    MNIndividual *representative;
    std::vector<SystemInfo *>members;
    double totalSharedFitness;
    int speciesNumber; // for data collection
};


class MONEAT {
protected:
    std::vector<InnovationInfo *> newConnections;
    std::vector<SystemInfo *> population;
    std::vector<NEATSpecies>speciesList;

    int populationSize;
    int stagnantGenerations;
    int maxStagnation;
    int maxGenerations;
    
    long generations;
    double allTimeBestFitness;
    double lastBestFitness; // before sharing
 
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one

    void mutateSystem(MNIndividual  *original); // does not make a copy
    
    // finds the Pareto fronts of individuals
    void rankSystems();

    virtual MNIndividual *combineSystems(SystemInfo *sys1, SystemInfo *sys2);
    virtual double genomeDistance( MNIndividual *sys1,  MNIndividual *sys2);

    virtual void assignInnovationNumberToAttachment(InnovationInfo *i);
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
        
    void logPopulationStatistics();
    NEATSpecies &chooseBreedingSpecies(double totalFitness); // fitness proportionate selection
    
    std::vector<InnovationInfo *> orderedConnectionGenome(std::vector<MNEdge *> genome);
    
public:
    MONEAT(int populationSize, int maxGenerations, int maxStagnation);
    ~MONEAT();
    
    // notice that this cannot be overriden.
    bool tick();
#pragma mark - Function pointers - MUST BE SET OR ELSE
    // MUST SET
    
    typedef  double (*EvaluationFunction)(MNIndividual  *individual);
    std::vector<EvaluationFunction> evaluationFunctions;
    
    MNIndividual *(*createInitialIndividual)(void); // the smallest/starting individual
#pragma mark - 
    
    std::vector<SystemInfo *> bestIndividuals; // the optimal front
    long getNumberOfIterations();
    
#pragma mark - Tuning parameters
    // probability of recombination taking from weaker parent
    double p_c = 0.01;
    
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
    MNIndividual *origin; // to find distance of species from start
    int nextInnovationNumber;
    int nextSpeciesNumber;
    void speciate(); // divide everything into species
    void updateSharedFitnesses();
};


#endif /* defined(__SystemGenerator__MONEAT__) */

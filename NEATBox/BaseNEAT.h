//
//  BaseNEAT.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/10/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef NEATBox_BaseNEAT_h
#define NEATBox_BaseNEAT_h

#include "NEATCommon.h"
#include "MNIndividual.h"
#include <vector>


class BaseNEAT {
protected:
    // define your population yourself...
    
    std::vector<InnovationInfo *> newConnections;
    
    long populationSize;
    long maxGenerations;
    
    long generations;
    
    MNIndividual *origin; // to find distance of species from start
    int nextInnovationNumber;

    std::vector<SystemInfo *> bestIndividuals; // the optimal front
    
    virtual void spawnNextGeneration() = 0; // recombine species to get enough children then mutate each one
    
    virtual void mutateSystem(MNIndividual  *original); // does not make a copy
    
    // finds the Pareto fronts of individuals
    virtual void rankSystems() = 0;
    
    virtual MNIndividual *combineSystems(SystemInfo *sys1, SystemInfo *sys2);
    virtual double genomeDistance( MNIndividual *sys1,  MNIndividual *sys2);
    
    virtual void assignInnovationNumberToGene(InnovationInfo *i);
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation() = 0;
    
    void logPopulationStatistics();
    
    std::vector<InnovationInfo *> orderedConnectionGenome(std::vector<MNEdge *> genome);
    
   
    
public:
    BaseNEAT(long populationSize, long maxGenerations);
    virtual ~BaseNEAT(){};
    
    virtual bool tick() = 0;
    virtual std::vector<SystemInfo *> optimalSolutions();
#pragma mark - Function pointers - MUST BE SET OR ELSE
    // MUST SET
    
    typedef  double (*EvaluationFunction)(MNIndividual  *individual);
    std::vector<EvaluationFunction> evaluationFunctions;
    
    MNIndividual *(*createInitialIndividual)(void); // the smallest/starting individual
#pragma mark -
    
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

};

#endif

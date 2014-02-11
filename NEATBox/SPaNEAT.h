//
//  SPaNEAT.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/9/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__SPaNEAT__
#define __NEATBox__SPaNEAT__

#include <iostream>


//
//  SPaNEAT.h
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#ifndef __SystemGenerator__SPaNEAT__
#define __SystemGenerator__SPaNEAT__

#include <vector>
#include <fstream>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>

#include "MNEdge.h"
#include "MNIndividual.h"

#include "NEATCommon.h"

struct SPSystemInfo : public SystemInfo
{
    std::vector<double> distances;
    
    SPSystemInfo(MNIndividual *i):SystemInfo(i) {}
    
    SPSystemInfo(const SPSystemInfo &other)
    :SystemInfo(other),
    distances(other.distances)
    {}
};

class SPaNEAT {
protected:
    std::vector<InnovationInfo *> newConnections;
    std::vector<SPSystemInfo *> population;
    std::vector<SPSystemInfo *> archive;
    
    long populationSize;
//    int stagnantGenerations;
//    int maxStagnation;
    long maxGenerations;
    
    long generations;
    
    long archiveSize;
    
    void spawnNextGeneration(); // recombine species to get enough children then mutate each one
    
    void mutateSystem(MNIndividual  *original); // does not make a copy
    
    // finds the Pareto fronts of individuals
    void rankSystems();
    
    virtual MNIndividual *combineSystems(SystemInfo *sys1, SystemInfo *sys2);
    virtual double genomeDistance( MNIndividual *sys1,  MNIndividual *sys2);
    
    virtual void assignInnovationNumberToGene(InnovationInfo *i);
    
    // overriden functions cannot be called in constructors, so this is called on the first tick();
    virtual void prepareInitialPopulation();
    
    void logPopulationStatistics();
    
    std::vector<InnovationInfo *> orderedConnectionGenome(std::vector<MNEdge *> genome);
    
public:
    SPaNEAT(long populationSize, long archiveSize, long maxGenerations=500);
    ~SPaNEAT();
    
    // notice that this cannot be overriden.
    bool tick();
    std::vector<SystemInfo *> optimalSolutions();

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
#pragma mark -
    
private:
    MNIndividual *origin; // to find distance of species from start
    int nextInnovationNumber;
};


#endif /* defined(__SystemGenerator__SPaNEAT__) */



#endif /* defined(__NEATBox__SPaNEAT__) */

//
//  SPaNEAT.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/9/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "SPaNEAT.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>
#include <assert.h>
#include <map>
#include <iostream>

#include "NBUtils.h"

#define USE_NODE_DIFF 1
#define SEPARATE_ARCHIVE 1

SPaNEAT::SPaNEAT(long populationSize, long archiveSize, long maxGenerations)
:populationSize(populationSize), archiveSize(archiveSize), maxGenerations(maxGenerations)

{
    nextInnovationNumber = 1;
    generations = 0;
    origin = NULL;
}

SPaNEAT::~SPaNEAT()
{
    //empty the archive
#if SEPARATE_ARCHIVE
    for (std::vector<SPSystemInfo *>::iterator it = archive.begin(); it != archive.end(); it++) {
        delete *it;
    }
    
#endif
    //empty the population
    for (std::vector<SPSystemInfo *>::iterator it = population.begin(); it != population.end(); it++) {
        delete *it;
    }
 
    delete origin;
}

// return true if sys1 comes before sys2
// smaller rank fitness is better
static bool compareIndividuals(SystemInfo *sys1, SystemInfo *sys2)
{
    return sys1->rankFitness < sys2->rankFitness;
}

static bool compareInnovationNumbers(const InnovationInfo *a1, const InnovationInfo *a2)
{
    return a1->innovationNumber < a2->innovationNumber;
}


std::vector<InnovationInfo*> SPaNEAT::orderedConnectionGenome(std::vector<MNEdge *> genome)
{
    std::vector<InnovationInfo *> newGenome;
    
    std::vector<MNEdge *>::iterator it = genome.begin();
    while (it != genome.end()) {
        InnovationInfo *innov = new InnovationInfo(*it);
        // assign the number to innov
        assignInnovationNumberToGene(innov);
        newGenome.push_back(innov);
        it++;
    }
    
    std::sort(newGenome.begin(), newGenome.end(), compareInnovationNumbers);
    
    return newGenome;
}


MNIndividual *SPaNEAT::combineSystems(SystemInfo *sys1, SystemInfo *sys2)
{
    
    if (sys1 == sys2) {
        // breeding with myself?
        return (sys1->individual->clone());
    }
    
    MNIndividual *newChild = createInitialIndividual();
    std::vector<InnovationInfo *> genome1 = orderedConnectionGenome(sys1->individual->connectionGenome());
    std::vector<InnovationInfo *> genome2 = orderedConnectionGenome(sys2->individual->connectionGenome());
    
    // holy crap, c++ has some good standard methods, like *set difference*
    std::vector<InnovationInfo *> matchingGenes1;
    std::vector<InnovationInfo *> matchingGenes2;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers);
    
    
    
    // first we handle the matching genes
    for (int i=0; i<matchingGenes1.size(); i++){
        InnovationInfo *gene = matchingGenes1[i];
        MNIndividual *parent = sys1->individual;
        // randomly choose 1 to add to the child
        double selector = (double )rand()/RAND_MAX;
        
        if (selector <= p_c) {
            parent = sys2->individual;
            gene = matchingGenes2[i];
        }
        newChild->addGeneFromParentSystem(parent, gene->data);
    }
    
    
    std::vector<InnovationInfo *> disjointandExcess1;
    std::vector<InnovationInfo *> disjointandExcess2;
    
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointandExcess1), compareInnovationNumbers);
    std::set_difference(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(), std::back_inserter(disjointandExcess2), compareInnovationNumbers);
    
    
    // then the disjoint genes
    if (sys1->rankFitness <= sys2->rankFitness) {
        for (int i=0; i<disjointandExcess1.size(); i++) {
            newChild->addGeneFromParentSystem(sys1->individual, disjointandExcess1[i]->data);
        }
    }
    
    // if rankFitness is equal, include all genes
    if (sys2->rankFitness <= sys1->rankFitness) {
        for (int i=0; i<disjointandExcess2.size(); i++) {
            newChild->addGeneFromParentSystem(sys2->individual, disjointandExcess2[i]->data);
        }
    }
    
    // clean out the allocated InnovationInfos
    for (std::vector<InnovationInfo *>::iterator dit = genome1.begin(); dit != genome1.end();) {
        delete *dit;
        dit = genome1.erase(dit);
    }
    
    for (std::vector<InnovationInfo *>::iterator dit = genome2.begin(); dit != genome2.end();) {
        delete *dit;
        dit = genome2.erase(dit);
    }
    
    
    return newChild;
}

void SPaNEAT::spawnNextGeneration()
{
    std::vector<SPSystemInfo *>newGeneration;
    
    // create population children through recombination and mutation of archive members
    while (newGeneration.size() < populationSize-archiveSize) {
        const int tournamentSize = 3;
        // select parent one
        SPSystemInfo *p1 = archive[uniformlyDistributed(archiveSize)];
        for (int i=0; i<tournamentSize-1; i++) {
             SPSystemInfo *contender =  archive[uniformlyDistributed(archiveSize)];
            if (contender->rankFitness < p1->rankFitness)
                p1 = contender;
        }
        
        SPSystemInfo *p2 = archive[uniformlyDistributed(archiveSize)];
        for (int i=0; i<tournamentSize-1; i++) {
            SPSystemInfo *contender =  archive[uniformlyDistributed(archiveSize)];
            if (contender->rankFitness < p2->rankFitness)
                p2 = contender;
        }
        
        MNIndividual *i = combineSystems(p1, p2);
        mutateSystem(i);

        SPSystemInfo *child = new SPSystemInfo(i);
        newGeneration.push_back(child);
    }
    
    // add the archive too (for now)
#if SEPARATE_ARCHIVE
    population = newGeneration;
#else
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());
#endif
}



// unlike base NEAT, also take node differences into account
double  SPaNEAT::genomeDistance( MNIndividual *sys1,  MNIndividual *sys2)
{
    std::vector<InnovationInfo *> genome1 = orderedConnectionGenome(sys1->connectionGenome());
    std::vector<InnovationInfo *> genome2 = orderedConnectionGenome(sys2->connectionGenome());
    
    // find the longer one
    size_t longerSize = (genome1.size() > genome2.size() ? genome1.size()  : genome2.size());
    
    size_t nDisjoint = 0;
    size_t nExcess = 0;
    size_t nMatching = 0;
    size_t matchingDiff = 0;
    
    // yay stdlib!
    std::vector<InnovationInfo *> matchingGenes1;
    std::vector<InnovationInfo *> matchingGenes2;
    std::vector<InnovationInfo *> disjointFromSys1;
    std::vector<InnovationInfo *> disjointFromSys2;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers);
    
    // difference takes things in the first range that are not in the second - we will have to check for excess ourself
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointFromSys1), compareInnovationNumbers);
    std::set_difference(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(), std::back_inserter(disjointFromSys2), compareInnovationNumbers);
    
    
    nMatching = matchingGenes2.size();
    // first find the distance between matching genes
    for (int i=0; i<nMatching; i++) {
        MNEdge *a1 = matchingGenes1[i]->data;
        MNEdge *a2 = matchingGenes2[i]->data;
        
        matchingDiff += sys1->connectionDifference(a1, a2);
    }
    
    // now determine the excess vs disjoint
    // if one of the disjoint sets is empty, then it's all excess
    if (disjointFromSys2.size() == 0 || disjointFromSys1.size() == 0) {
        nExcess = disjointFromSys1.size() + disjointFromSys2.size();
    } else {
        // else chalk everything up to disjoint
        nDisjoint = disjointFromSys1.size() + disjointFromSys2.size();
        
        // then find the number of excess genes
        int last1 = disjointFromSys1.back()->innovationNumber;
        int last2 = disjointFromSys2.back()->innovationNumber;
        
        std::vector<InnovationInfo *> shorter;
        std::vector<InnovationInfo *> longer;
        
        if (last1 != last2) {
            if (last1 > last2) {
                // genome 1 has excess
                shorter = disjointFromSys2;
                longer = disjointFromSys1;
            } else if (last1 < last2) {
                // genome 2 has excess, switch last1 and last2
                shorter=disjointFromSys1;
                longer = disjointFromSys2;
                
                int temp = last1;
                last1 = last2;
                last2 = temp;
            }
            // go back from end of longer until we run into a value <= the end of shorter
            size_t i;
            for (i=longer.size(); i; i--) {
                if (longer[i-1]->innovationNumber < shorter.back()->innovationNumber)
                    break;
            }
            nExcess = i;
            nDisjoint -= nExcess;
        }
    }
    
    
    double diff_node = sys1->nodeDifference(sys2);
    long moreNodes = std::min(sys1->numberOfNodes(), sys2->numberOfNodes());
    
   // double  d = w_disjoint*nDisjoint/longerSize + w_excess*nExcess/longerSize + w_matching*matchingDiff/nMatching;
    
        double  d = w_disjoint*nDisjoint + w_excess*nExcess + w_matching*matchingDiff;
#if USE_NODE_DIFF
    //d+= + w_matching_node*diff_node/moreNodes;
    d += w_matching_node*diff_node;
#endif
    
    // clean out the allocated InnovationInfos
    for (std::vector<InnovationInfo *>::iterator dit = genome1.begin(); dit != genome1.end();) {
        delete *dit;
        dit = genome1.erase(dit);
    }
    
    for (std::vector<InnovationInfo *>::iterator dit = genome2.begin(); dit != genome2.end();) {
        delete *dit;
        dit = genome2.erase(dit);
    }
    
    
    return fabs(d);
}


void SPaNEAT::prepareInitialPopulation()
{
    // we need to index the genome of the 'base' system
    origin = createInitialIndividual();
    std::vector<InnovationInfo *> baseGenome = orderedConnectionGenome(origin->connectionGenome());
    
    // clean out the allocated InnovationInfos
    for (std::vector<InnovationInfo *>::iterator dit = baseGenome.begin(); dit != baseGenome.end();) {
        delete *dit;
        dit = baseGenome.erase(dit);
    }
    
    
    // mutate the initial system to get an initial population
    while (population.size() < populationSize) {
        MNIndividual *newSystem = createInitialIndividual();
        mutateSystem(newSystem);
        SPSystemInfo *info = new SPSystemInfo(newSystem);
        population.push_back(info);
    }
}

// does s1 dominate s2?
// i.e. all s1 fitness values are greater than s2?
static int domination(SPSystemInfo *s1, SPSystemInfo *s2)
{
    bool firstDominates = true;
    bool secondDominates = true;
    
    std::vector<double>::iterator fitnessIter1 = s1->fitnesses.begin();
    std::vector<double>::iterator fitnessIter2 = s2->fitnesses.begin();
    
    while (fitnessIter1 != s1->fitnesses.end()) {
        if (*fitnessIter1 > *fitnessIter2)
            firstDominates = false;
        
        if (*fitnessIter2 > *fitnessIter1)
            secondDominates = false;
        fitnessIter2++;
        fitnessIter1++;
    }
    
    if (firstDominates == secondDominates)
        return 0;
    
    if (firstDominates)
        return 1;
    
    if (secondDominates)
        return -1;
    
    return 0;
}

// using SPEA2
void SPaNEAT::rankSystems()
{
    std::vector<SPSystemInfo *>::iterator popIter;
    
    std::map<SPSystemInfo *, std::vector<SPSystemInfo *> > dominationMap;
    
    int k = sqrt(populationSize+archiveSize);
#if SEPARATE_ARCHIVE
    std::vector<SPSystemInfo *> toRank = std::vector<SPSystemInfo *>(population);
    toRank.insert(toRank.end(), archive.begin(), archive.end());
#else
    std::vector<SPSystemInfo *> toRank = population;
#endif
    // we're assigning a strength number based on # of dominated individuals
    // domination map and count take  opposite meanings from NSGA-II
    
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        std::vector<SPSystemInfo *> dominators;
        SPSystemInfo *sys = *popIter;
        
        sys->dominationCount = 0;
        sys->distances.clear();
        
        typename std::vector<SPSystemInfo *>::iterator innerIter;
        for (innerIter = toRank.begin(); innerIter != toRank.end(); innerIter++) {
            if (*innerIter == *popIter)
                continue;
            SPSystemInfo *otherSys = *innerIter;
            
            int dominationScore = domination( otherSys, sys);
            
            if (dominationScore == 1) {
                // *innerIter dominates
                dominators.push_back(otherSys);
            }
            if (dominationScore == -1) {
                // *popIter dominates
                sys->dominationCount += 1;
            }
            sys->distances.push_back(genomeDistance(sys->individual, otherSys->individual));
        }
        dominationMap[sys] = dominators;
        
        // the density calculation can come here, since we have all the distances any way
        std::sort(sys->distances.begin(), sys->distances.end());
        double kDist = sys->distances[k];
        sys->rankFitness = 1.0/(kDist + 2);
    }
    
    
    // now assign raw fitness based on the strengths of dominators, plus population neighborhood density
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        SPSystemInfo *sys = *popIter;
        std::vector<SPSystemInfo *> dominators = dominationMap[sys];
        std::vector<SPSystemInfo *>::iterator it;
        for (it = dominators.begin(); it != dominators.end(); it++) {
            sys->rankFitness += (*it)->dominationCount;
        }
    }
}

bool SPaNEAT::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
        first_run = true;
    }
    
    std::vector<SPSystemInfo *>::iterator popIter;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        SPSystemInfo *sys = *popIter;
        sys->fitnesses.clear();
        
        std::vector<EvaluationFunction>::iterator funcIter;
        for (funcIter =  evaluationFunctions.begin(); funcIter != evaluationFunctions.end(); funcIter++) {
            double val = (*funcIter)(sys->individual);
            sys->fitnesses.push_back(val);
        }
    }
    
    rankSystems();
    
    bool last_run = (generations >= maxGenerations) ;
    
    std::vector<SPSystemInfo *> newArchive;
    
#if SEPARATE_ARCHIVE
    std::vector<SPSystemInfo *> toRank = std::vector<SPSystemInfo *>(population);
    toRank.insert(toRank.end(), archive.begin(), archive.end());
#else
    std::vector<SPSystemInfo *> toRank = population;
#endif
    
    std::sort(toRank.begin(), toRank.end(), compareIndividuals);
    
    std::vector<SPSystemInfo *>::iterator lastArchived = population.begin();
    for (popIter = toRank.begin(); popIter != toRank.end(); popIter++) {
        if ((*popIter)->rankFitness < 1.0) {
            newArchive.push_back(new SPSystemInfo(**popIter));
            lastArchived = popIter;
        }
    }
    
    std::cout << newArchive.size() << " non dominated solutions (" << generations << ").\n";

    
    long archiveExcess = newArchive.size() - archiveSize;
    
    if (archiveExcess < 0) {
        // add in the next best -archiveExcess solutions
        for (popIter = lastArchived+1; popIter != lastArchived - archiveExcess +1; popIter++) {
            newArchive.push_back(new SPSystemInfo(**popIter));
        }
    }
    
    if (archiveExcess > 0) {
        // remove by crowding distance...
        std::sort(newArchive.begin(), newArchive.end(), [](const SPSystemInfo * elem1, const SPSystemInfo *elem2){return elem1->distances < elem2->distances;});
        
        popIter = newArchive.begin();
        for (long i=0; i<archiveExcess; i++) {
            delete *popIter;
            popIter = newArchive.erase(popIter);
        }
    }

    generations++;

    // delete everything in the old population
    // XXX: this includes the old archive (for now)
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        delete  *popIter;
    }
    population.clear();
    
#if SEPARATE_ARCHIVE
    for (popIter = archive.begin(); popIter != archive.end(); popIter++) {
        delete  *popIter;
    }
#endif
    
    archive = newArchive;
    
    if (!last_run) {
        spawnNextGeneration();
    }
    
    return  last_run;
}


long SPaNEAT::getNumberOfIterations()
{
    return generations;
}


void SPaNEAT::mutateSystem(MNIndividual *original)
{
    
    std::vector<MNEdge *> allAttachments = original->connectionGenome();
    typename std::vector<MNEdge *>::iterator it = allAttachments.begin();
    while (it != allAttachments.end()) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_conn)
            original->mutateConnectionWeight();
        it++;
    }
    
    long nNodes = original->numberOfNodes();
    for (long i=0; i<nNodes; i++) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_node)
            original->mutateNode(i);
    }
    
    // add attachment or add node - or both
    double selector1 = (double )rand()/RAND_MAX;
    if (selector1 < p_m_conn_ins) {
        // insert a new connection
        // update the innovation number
        MNEdge *created = original->createConnection();
        InnovationInfo *newInfo = new InnovationInfo(created, false);
        assignInnovationNumberToGene(newInfo);
        delete newInfo;
    }
    
    double selector2 = (double )rand()/RAND_MAX;
    if (selector2 < p_m_node_ins) {
        // add a node somewhere if possible
        std::vector<MNEdge *> new_edges = original->createNode();
        for (int i=0; i<new_edges.size(); i++) {
            InnovationInfo *newInfo = new InnovationInfo(new_edges[i]);
            assignInnovationNumberToGene(newInfo);
            delete newInfo;
        }
    }
    
}



void SPaNEAT::assignInnovationNumberToGene(InnovationInfo *i){
    // have we already created this 'innovation' in this generation?
    pointer_values_equal<InnovationInfo> pred = {i};
    std::vector<InnovationInfo *>::iterator it = std::find_if(newConnections.begin(), newConnections.end(), pred);
    
    InnovationInfo *found = NULL;
    
    if (it== newConnections.end()) {
        // not already existing
        found =  new InnovationInfo(i->data); // clone the data
        found->innovationNumber = nextInnovationNumber++;
        newConnections.push_back(found);
    } else {
        found = *it;
    }
    i->innovationNumber = found->innovationNumber;
}

void SPaNEAT::logPopulationStatistics()
{
    
}

std::vector<SystemInfo *> SPaNEAT::optimalSolutions()
{
    // some suboptimal solutions might be in here: make sure we only get the nondominated ones
        std::vector<SystemInfo *>converted;
    
    std::vector<SPSystemInfo *>::iterator it = std::find_if(archive.begin(), archive.end(), [](const SystemInfo * elem){return elem->rankFitness < 1.0;});
    for (; it != archive.end(); it = std::find_if(++it, archive.end(), [](const SystemInfo * elem){return elem->rankFitness < 1.0;})) {
        converted.push_back(*it);
    }

    return converted;
}
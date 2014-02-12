//
//  BaseNEAT.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/10/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "BaseNEAT.h"
#include "NBUtils.h"

#define USE_NODE_DIFF 1

BaseNEAT::BaseNEAT(long populationSize, long maxGenerations)
:populationSize(populationSize), maxGenerations(maxGenerations),
nextInnovationNumber(1), generations(0), origin(NULL)
{
}

std::vector<SystemInfo *> BaseNEAT::optimalSolutions()
{
    return bestIndividuals;
}



void BaseNEAT::mutateSystem(MNIndividual  *original)
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
            InnovationInfo *newInfo = new InnovationInfo(new_edges[i], false);
            assignInnovationNumberToGene(newInfo);
            delete newInfo;
        }
    }
}


static bool compareInnovationNumbers(const InnovationInfo *a1, const InnovationInfo *a2)
{
    return a1->innovationNumber < a2->innovationNumber;
}

MNIndividual * BaseNEAT::combineSystems(SystemInfo *sys1, SystemInfo *sys2)
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
        double selector = uniformProbability();
        
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

double BaseNEAT::genomeDistance(MNIndividual *sys1,  MNIndividual *sys2)
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
    
    double  d = w_disjoint*nDisjoint/longerSize + w_excess*nExcess/longerSize + w_matching*matchingDiff/nMatching;
    
    //    double  d = w_disjoint*nDisjoint + w_excess*nExcess + w_matching*matchingDiff;
#if USE_NODE_DIFF
    d+= + w_matching_node*diff_node/moreNodes;
    // d += w_matching_node*diff_node
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
    
    
    return d;
}

void BaseNEAT::assignInnovationNumberToGene(InnovationInfo *i) {
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

void BaseNEAT::logPopulationStatistics() {}

std::vector<InnovationInfo *> BaseNEAT::orderedConnectionGenome(std::vector<MNEdge *> genome)
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

long BaseNEAT::getNumberOfIterations()
{
    return generations;
}

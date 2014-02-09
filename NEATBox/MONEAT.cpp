
//
//  MONEAT.cpp
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#include "MONEAT.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>
#include <assert.h>
#include <map>
#include <iostream>

#include "NBUtils.h"

#define RESTRICT_SPECIES 1
#define ROULETTE_SELECT 0
#define INTERSPECIES 0

#define TIGHT_CLUSTERS 0

#define USE_NODE_DIFF 1

// SO DIRTY
struct InnovationInfo {
//private:
    InnovationInfo(const InnovationInfo& other)
    :data(other.data), innovationNumber(other.innovationNumber), wasCloned(false)
    {
        
    }
//public:
    MNEdge *data;
    bool wasCloned;
    int innovationNumber;
    
    InnovationInfo(MNEdge *e)
    :data(e->clone()),innovationNumber(-1), wasCloned(true){}
    
    InnovationInfo *clone()
    {
        InnovationInfo *newOne = new InnovationInfo(data->clone());
        newOne->innovationNumber = innovationNumber;
        newOne->wasCloned = true;
        return newOne;
    }
    
    
    ~InnovationInfo(){
        if (wasCloned)
            delete data;
    }
    
    bool operator==(const InnovationInfo& other)const {
        return *(data) == *(other.data);
    }
    
    bool operator==(const MNEdge& other)const {
        return *data == other;
    }
};

/// OH GOD I HATE C++ FUNCTORSSSS
// http://stackoverflow.com/questions/258871/how-to-use-find-algorithm-with-a-vector-of-pointers-to-objects-in-c
template <typename T>
struct pointer_values_equal
{
    const T* to_find;
    
    bool operator()(const T* other) const
    {
        return *to_find == *other;
    }
};

MONEAT::MONEAT(int populationSize, int maxGenerations, int maxStagnation)
:populationSize(populationSize), maxStagnation(maxStagnation), maxGenerations(maxGenerations)

{
    nextInnovationNumber = 1;
    nextSpeciesNumber = 1;
    allTimeBestFitness = -INFINITY;
    generations = 0;
    stagnantGenerations = 0;
    origin = NULL;
}

MONEAT::~MONEAT()
{
    //empty the species lists
    std::vector<NEATSpecies>::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo *> members = speciesIterator->members;
         std::vector<SystemInfo *>::iterator it = members.begin();
        while (it != members.end()) {
            it = members.erase(it);
        }
        speciesIterator = speciesList.erase(speciesIterator);
    }
    delete origin;
}

// return true if sys1 comes before sys2
static bool compareIndividuals(SystemInfo *sys1, SystemInfo *sys2)
{
    return sys1->rankFitness > sys2->rankFitness;
}

static bool compareInnovationNumbers(const InnovationInfo *a1, const InnovationInfo *a2)
{
    return a1->innovationNumber < a2->innovationNumber;
}


std::vector<InnovationInfo*> MONEAT::orderedConnectionGenome(std::vector<MNEdge *> genome)
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


MNIndividual *MONEAT::combineSystems(SystemInfo *sys1, SystemInfo *sys2)
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
    if (sys1->rankFitness >= sys2->rankFitness) {
        for (int i=0; i<disjointandExcess1.size(); i++) {
            newChild->addGeneFromParentSystem(sys1->individual, disjointandExcess1[i]->data);
        }
    }
    
    // if rankFitness is equal, include all genes
    if (sys2->rankFitness >= sys1->rankFitness) {
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

void MONEAT::spawnNextGeneration()
{
    double  rankFitnessSum = 0.0;
    typename std::vector<SystemInfo *>::iterator it = population.begin();
    while (it != population.end()) {
        rankFitnessSum += (*it)->rankFitness;
        it++;
    }

    // create population children through recombination
    std::vector<SystemInfo *>newGeneration;
    
     std::vector<NEATSpecies>::iterator speciesIt = speciesList.begin();
    while(speciesIt != speciesList.end()) {
        // selectionnnnn
        NEATSpecies breedingSpecies = *speciesIt;//chooseBreedingSpecies(rankFitnessSum);
        std::vector<SystemInfo *>newMembers;
        
        std::vector<SystemInfo *>individuals = breedingSpecies.members;
        
        while (newMembers.size() < individuals.size()) {
            SystemInfo *p1;
            SystemInfo *p2;
            
            MNIndividual *child;
            
            int parentPosition = arc4random_uniform((int)individuals.size());
            p1 = (individuals[parentPosition]);
#if INTERSPECIES
            double interspecies_mate = (double)rand()/RAND_MAX;
            if (stagnantGenerations > maxStagnation/3 || interspecies_mate < 0.1) {
                breedingSpecies = chooseBreedingSpecies(rankFitnessSum);
                individuals = breedingSpecies.members;
                
                int parentPosition = arc4random_uniform((int)individuals.size());
                p2 = (individuals[parentPosition]);
            } else {
#endif
                if (individuals.size() > 1) {
                    
                    int parentPosition2;
                    do {
                        parentPosition2 = arc4random_uniform((int)individuals.size());
                    } while (parentPosition2 == parentPosition);
                    p2 = (individuals[parentPosition2]);
                } else {
                    p2 = p1;
                }
#if INTERSPECIES
            }
#endif
            
            child = combineSystems(p1, p2);
            // mutate the new individual (maybe)
            mutateSystem(child);
            
            SystemInfo *i = new SystemInfo(child);
            newGeneration.push_back(i);
            newMembers.push_back(i);
        }
        speciesIt->members.insert(speciesIt->members.end(), newMembers.begin(), newMembers.end());
        speciesIt++;
    }
    // add the original population too
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());

//    population = newGeneration;
}

NEATSpecies &MONEAT::chooseBreedingSpecies(double totalFitness)
{
    
    // choose a species according to species fitness
    // then chose parents at random from within that species
    
    // Assume sorted population
    double selector = ((double )rand()/RAND_MAX);
    
    typename std::vector<NEATSpecies >::iterator it = speciesList.begin();
    double  cumulativeFitness = 0.0;
    while (it != speciesList.end()) {
        cumulativeFitness += it->totalSharedFitness/totalFitness;
        if (cumulativeFitness >= selector) {
            break;
        }
        it++;
    }
    
    if (it == speciesList.end())
        it--;
    
    return *it;
}


void MONEAT::speciate()
{
     std::vector<SystemInfo *>::iterator populationIter = population.begin();
    while (populationIter != population.end()) {

        // assign to species
        bool added= false;
#if TIGHT_CLUSTERS
        double smallestDiff = INFINITY;
        std::vector<NEATSpecies>::iterator chosenSpecies = speciesList.end();
#endif
        
         std::vector<NEATSpecies >::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            double  distance = fabs(genomeDistance((*populationIter)->individual, speciesIterator->representative));
            if (distance < d_threshold) {
#if TIGHT_CLUSTERS
                if (distance < smallestDiff) {
                    smallestDiff = distance;
                    chosenSpecies = speciesIterator;
                }
#else
                added = true;
                speciesIterator->members.push_back((*populationIter));
                break;
#endif
            }
            speciesIterator++;
        }
#if TIGHT_CLUSTERS
        if (chosenSpecies != speciesList.end()) {
            added = true;
            chosenSpecies->members.push_back((*populationIter));
        }
#endif
        
        if (!added) {
            // new species!!
            NEATSpecies s = NEATSpecies();
            s.representative = (*populationIter)->individual->clone();
            s.members.push_back(*populationIter);
            speciesList.push_back(s);
            s.speciesNumber = nextSpeciesNumber++;
        }
        
        populationIter++;
    }
}

void MONEAT::updateSharedFitnesses()
{
    long nSpecies = speciesList.size();
    fprintf(stderr, "%ld species\n", nSpecies);
    
    // for each species, copy a random current member to be the representative for the next generation, and adjust the fitnesses for sharing
    // kill off a species with no members
    typename std::vector<NEATSpecies >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        // pick a new rep
        std::vector<SystemInfo *> members = speciesIterator->members;
        if (members.size() == 0) {
            NEATSpecies extinctSpecies = *speciesIterator;
            speciesIterator = speciesList.erase(speciesIterator);
        } else {
            long index = uniformlyDistributed(members.size());
            MNIndividual *rep = (members[index]->individual);
            delete speciesIterator->representative;
            speciesIterator->representative = rep->clone();
            
            // fprintf(stderr, "\tSpecies %d: %ld members\n", ++i, members.size());
            
            // update the fitnesses
            double  totalFitness = 0.0;
            typename std::vector<SystemInfo *>::iterator memberIterator = members.begin();
            while (memberIterator != members.end()) {
                assert((*memberIterator)->rankFitness != -INFINITY);
                (*memberIterator)->rankFitness /= members.size();
                totalFitness += (*memberIterator)->rankFitness;
                memberIterator++;
            }
            
            speciesIterator->totalSharedFitness = totalFitness;
            speciesIterator++;
        }
    }
}

// unlike base NEAT, also take node differences into account
double  MONEAT::genomeDistance( MNIndividual *sys1,  MNIndividual *sys2)
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


void MONEAT::prepareInitialPopulation()
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
        SystemInfo *info = new SystemInfo(newSystem);
        population.push_back(info);
    }
}

// does s1 dominate s2?
// i.e. all s1 fitness values are greater than s2?
static int domination(SystemInfo *s1, SystemInfo *s2)
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


void MONEAT::rankSystems()
{
    typename std::vector<SystemInfo *>::iterator popIter;
    
    std::map<SystemInfo *, std::vector<SystemInfo *> > dominationMap;
    
    int counted = 0;

    std::vector<SystemInfo *> bestFront;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        std::vector<SystemInfo *> dominated;
        SystemInfo *sys = *popIter;
        sys->dominationCount = 0;
        typename std::vector<SystemInfo *>::iterator innerIter;
        for (innerIter = population.begin(); innerIter != population.end(); innerIter++) {
            if (*innerIter == *popIter)
                continue;
            
            int dominationScore = domination( *innerIter, *popIter);
            
            if (dominationScore == -1) {
                // *popIter dominates
                dominated.push_back(*innerIter);
            }
            if (dominationScore == 1) {
                // *innerIter dominates
                sys->dominationCount += 1;
            }
        }
        if (sys->dominationCount == 0) {
            counted++;
            sys->rankFitness = 1.0;
            bestFront.push_back(sys);
        }
        dominationMap[sys] = dominated;
        
    }
    
    int currentRank = 1;
    std::vector<SystemInfo *> currentFront = bestFront;
    while (!currentFront.empty()) {
        std::vector<SystemInfo *> nextFront;
        typename std::vector<SystemInfo *>::iterator frontIterator;
        
        for (frontIterator = currentFront.begin(); frontIterator != currentFront.end(); frontIterator++) {
            std::vector<SystemInfo *> dominated = dominationMap[*frontIterator];
            typename std::vector<SystemInfo *>::iterator it;
            for (it = dominated.begin(); it != dominated.end(); it++) {
                SystemInfo *s = *it;
                s->dominationCount -= 1;
               // assert(s->dominationCount >=0);
                if (s->dominationCount == 0) {
                    s->rankFitness = 1.0/(currentRank+1);
                    nextFront.push_back(s);
                    counted++;
                }
            }
        }

        currentRank += 1;
        currentFront = nextFront;
    }
    
    assert(counted == population.size());

    std::cout << currentRank << " ranks examined.\n";
    
    //return bestFront;
}

// descending, not ascending, order
static bool compareSpeciesFitness(const NEATSpecies& s1, const NEATSpecies& s2)
{
    return s1.totalSharedFitness > s2.totalSharedFitness;
}

bool MONEAT::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
     
        first_run = true;
    }
    
    typename std::vector<SystemInfo *>::iterator popIter;
    for (popIter = population.begin(); popIter != population.end(); popIter++) {
        SystemInfo *sys = *popIter;
        sys->fitnesses.clear();
        
        typename std::vector<EvaluationFunction>::iterator funcIter;
        for (funcIter =  evaluationFunctions.begin(); funcIter != evaluationFunctions.end(); funcIter++) {
            double val = (*funcIter)(sys->individual);
            sys->fitnesses.push_back(val);
        }
    }
    
    rankSystems();
    
    bool last_run =  (generations >= maxGenerations) ;

    
#if !RESTRICT_SPECIES
   {
       // clear species membership lists
       typename std::vector<NEATSpecies >::iterator speciesIterator = speciesList.begin();
       while (speciesIterator != speciesList.end()) {
           speciesIterator->members.clear();
           speciesIterator++;
       }
   }
#endif
    
#if RESTRICT_SPECIES
    if (first_run) {
#endif
        speciate();
#if RESTRICT_SPECIES
    } else {
    
          bool regroup = false;
          if (speciesList.size() < 10 || stagnantGenerations > maxStagnation/5) {
              regroup = true;
              d_threshold /= 1.05;
          }
      
          if (speciesList.size() > 20) {
              regroup = true;
              d_threshold *= 1.05;
          }
      
          if (regroup) {
              // empty the species member lists and remake them
              speciesList.clear();
              speciate();
          }
          
    }
#endif
    updateSharedFitnesses();

    
    double  bestFitness = -INFINITY; // we want zero fitness
    // figure out how many of each species to save
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += speciesList[i].totalSharedFitness;
    }
    
    std::vector<SystemInfo *> individualsToSave;
    // sort the species to ensure we keep the best
    std::sort(speciesList.begin(), speciesList.end(), compareSpeciesFitness);
    
    std::vector<NEATSpecies >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo *> members =speciesIterator->members;
        size_t numMembers = members.size();

        double  proportionToSave = (speciesIterator->totalSharedFitness)/sharedFitnessSum;
        long numToSave = std::min((size_t)(proportionToSave*populationSize), numMembers);
        
        assert(numToSave >= 0);
        
        std::sort(members.begin(), members.end(), compareIndividuals);
        
        // don't let population grow unchecked
        if (individualsToSave.size() > populationSize)
            numToSave = 0;
        
        
        if (numToSave > 0) {
            std::vector<SystemInfo *> newMembers;
            std::vector<SystemInfo *>::iterator it = members.begin();

#if ROULETTE_SELECT
            //stochastic universal selection

            double choiceSep = (speciesIterator->totalSharedFitness/numToSave);
            double pointer = ((double)rand()/RAND_MAX)*choiceSep;
            std::vector<long> savedPositions;
            double cumFitness = (*it)->rankFitness;
            while (newMembers.size() < numToSave) {
                if (cumFitness >= pointer) {
                    pointer += choiceSep;
                    
                    SystemInfo *individual = *it;
                    
                    SystemInfo *saved = new SystemInfo(**it);
                    double  rankFitness = individual->rankFitness * numMembers;
                    saved->rankFitness = rankFitness;
                    newMembers.push_back(saved);
                    individualsToSave.push_back(saved);
                    if (rankFitness > bestFitness) {
                        bestFitness = rankFitness;
                        bestFound = copied;
                    }
                    
                } else {
                    it++;
                    assert(it != members.end());
                    cumFitness += (*it)->rankFitness;
                }
            }
#else
            while (newMembers.size() < numToSave) {
                
                SystemInfo *individual = *it;
                
                SystemInfo *saved = new SystemInfo(**it);
                double  rankFitness = individual->rankFitness * numMembers;
                saved->rankFitness = rankFitness;
                newMembers.push_back(saved);
                individualsToSave.push_back(saved);
                if (rankFitness > bestFitness) {
                    bestFitness = rankFitness;
                }
                it++;
            }
#endif
            
            std::vector<SystemInfo *>::iterator deleteIt =  members.begin();
            
            while (deleteIt !=  members.end()) {
                delete *deleteIt;
                deleteIt =  members.erase(deleteIt);
            }
            
            speciesIterator->members = newMembers;
            speciesIterator++;

        } else {
            speciesIterator = speciesList.erase(speciesIterator);
 
        }
        
    }
 
    
    population = individualsToSave;
    
    if (bestFitness > allTimeBestFitness)
        allTimeBestFitness = bestFitness;
    
    // stagnation if fitnesses are within 1% of each other
    
    if (fabs( 1 - (lastBestFitness/bestFitness)) < .01)
        stagnantGenerations++;
    else
        stagnantGenerations = 0;
    
    generations++;
    
    lastBestFitness = bestFitness;
   
    bestIndividuals.clear();
    // find the top ranking individuals
     std::vector<SystemInfo *>::iterator it = std::find_if(individualsToSave.begin(), individualsToSave.end(), [](const SystemInfo * elem){return elem->rankFitness == 1.0;});
    for (; it != individualsToSave.end(); it = std::find_if(++it, individualsToSave.end(), [](const SystemInfo * elem){return elem->rankFitness == 1.0;})) {
        bestIndividuals.push_back(*it);
    }
    
    
    if (!last_run) {
        spawnNextGeneration();
    }
    
    return  last_run;
}


long MONEAT::getNumberOfIterations()
{
    return generations;
}


void MONEAT::mutateSystem(MNIndividual *original)
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
        InnovationInfo *newInfo = new InnovationInfo(created);
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



void MONEAT::assignInnovationNumberToGene(InnovationInfo *i){
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

void MONEAT::logPopulationStatistics()
{
   
}


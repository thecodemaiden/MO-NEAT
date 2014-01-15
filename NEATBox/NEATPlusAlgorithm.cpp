
//
//  NEATPlusAlgorithm.cpp
//  SystemGenerator
//
//  Created by Adeola Bannis on 11/4/13.
//
//

#include "NEATPlusAlgorithm.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <wordexp.h>



template <class IndividualType, class InnovationType>
NEATPlusAlgorithm <IndividualType, InnovationType>::NEATPlusAlgorithm(int populationSize, int maxGenerations, int maxStagnation)
:populationSize(populationSize), maxStagnation(maxStagnation), maxGenerations(maxGenerations)

{
    _bestIndividual = NULL;
    nextInnovationNumber = 1;
    nextSpeciesNumber = 1;
    allTimeBestFitness = -INFINITY;
    generations = 0;
    stagnantGenerations = 0;
    origin = NULL;
}

template <class IndividualType, class InnovationType>
NEATPlusAlgorithm <IndividualType, InnovationType>::~NEATPlusAlgorithm()
{
    //empty the species lists
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        std::vector<SystemInfo<IndividualType> *> members = speciesIterator->members;
        typename std::vector<SystemInfo<IndividualType> *>::iterator it = members.begin();
        while (it != members.end()) {
            it = members.erase(it);
        }
        speciesIterator = speciesList.erase(speciesIterator);
    }
    currentLogFile.close();
    delete origin;
}

// return true if sys1 comes before sys2
template <class IndividualType>
static bool compareIndividuals(SystemInfo<IndividualType> *sys1, SystemInfo<IndividualType> *sys2)
{
    return sys1->fitness > sys2->fitness;
}

template <class InnovationType>
static bool compareInnovationNumbers(const InnovationType a1, const InnovationType a2)
{
    return a1.innovationNumber < a2.innovationNumber;
}

template <class IndividualType, class InnovationType>
IndividualType NEATPlusAlgorithm <IndividualType, InnovationType>::combineSystems(IndividualType &sys1, IndividualType &sys2)
{
    
//    if (&sys1 == &sys2) {
//        // breeding with myself?
//        return IndividualType(sys1);
//    }
    
    IndividualType newChild = createInitialIndividual();
    std::vector<InnovationType> genome1 = sys1.connectionGenome();
    std::vector<InnovationType> genome2 = sys2.connectionGenome();
    
    // holy crap, c++ has some good standard methods, like *set difference*
    std::vector<InnovationType> matchingGenes1;
    std::vector<InnovationType> matchingGenes2;
    std::vector<InnovationType> disjointandExcess;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers<InnovationType>);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers<InnovationType>);
    
    // difference takes things in the first range that are not in the second - perfect if we assume parent 1 is more fit
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointandExcess), compareInnovationNumbers<InnovationType>);
    
    
    // first we handle the matching genes
    for (int i=0; i<matchingGenes1.size(); i++){
        InnovationType gene = matchingGenes1[i];
        IndividualType parent = sys1;
        // randomly choose 1 to add to the child
        double selector = (double )rand()/RAND_MAX;
        
        if (selector <= p_c) {
            parent = sys2;
            gene = matchingGenes2[i];
        }
        newChild.addGeneFromParentSystem(parent, gene);
    }
    
    // then the disjoint genes
    for (int i=0; i<disjointandExcess.size(); i++) {
        newChild.addGeneFromParentSystem(sys1, disjointandExcess[i]);
    }
    
    return newChild;
}

//fitness proportionate selection
template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::spawnNextGeneration()
{
    double  fitnessSum = 0.0;
    typename std::vector<SystemInfo<IndividualType> >::iterator it = population.begin();
    while (it != population.end()) {
        fitnessSum += it->fitness;
        it++;
    }
    
    // create population children through recombination
    std::vector<SystemInfo<IndividualType> >newGeneration;
    
    do {

        // selectionnnnn
        NEATSpecies<IndividualType> &breedingSpecies  = chooseBreedingSpecies(fitnessSum);
        
        std::vector<SystemInfo<IndividualType> *>individuals = breedingSpecies.members;
        
        SystemInfo<IndividualType> *p1;
        SystemInfo<IndividualType> *p2;
        
        int parentPosition = arc4random_uniform((int)individuals.size());
        int parentPosition2 = parentPosition;

        p1 = (individuals[parentPosition]);
        
        if (individuals.size() > 1) {
            do {
                parentPosition2 = arc4random_uniform((int)individuals.size());
            } while (parentPosition2 == parentPosition);
        }
        p2 = (individuals[parentPosition2]);
        
        if (p2->fitness < 0.01 || p2->fitness > 1000)
            fprintf(stderr, "WTF");
        
        IndividualType child = combineSystems(p1->individual, p2->individual);
        // mutate the new individual (maybe)
        mutateSystem(child);
        
        SystemInfo<IndividualType> i = SystemInfo<IndividualType>(child);
        newGeneration.push_back(i);
        
        SystemInfo<IndividualType> *temp = &newGeneration.back();
        
        breedingSpecies.members.push_back(temp);
    } while (newGeneration.size() < populationSize);
    
    // add the original population too
    population.insert(population.end(), newGeneration.begin(), newGeneration.end());
    // population = newGeneration;
}

////fitness proportionate selection
//template <class IndividualType, class InnovationType>
//void NEATPlusAlgorithm <IndividualType, InnovationType>::spawnNextGeneration()
//{
//    double  fitnessSum = 0.0;
//    typename std::vector<SystemInfo<IndividualType> >::iterator it = population.begin();
//    while (it != population.end()) {
//        fitnessSum += it->fitness;
//        it++;
//    }
//    
//    // create population children through recombination
//    std::vector<SystemInfo<IndividualType> >newGeneration;
//    
//    do {
//        // selectionnnnn
//        NEATSpecies<IndividualType> breedingSpecies = chooseBreedingSpecies(fitnessSum);
//        
//        std::vector<SystemInfo<IndividualType> >individuals = breedingSpecies.members;
//        
//        SystemInfo<IndividualType> p1 = (individuals[0]);
//        SystemInfo<IndividualType> p2 = (individuals[0]);
//        
//        int parentPosition = arc4random_uniform((int)individuals.size());
//        p1 = SystemInfo<IndividualType>(individuals[parentPosition]);
//        
//        double interspecies_mate = (double)rand()/RAND_MAX;
//        if (stagnantGenerations > maxStagnation/2 || interspecies_mate < 0.05) {
//            breedingSpecies = chooseBreedingSpecies(fitnessSum);
//            individuals = breedingSpecies.members;
//            
//            int parentPosition = arc4random_uniform((int)individuals.size());
//            p2 = SystemInfo<IndividualType>(individuals[parentPosition]);
//        } else
//            if (individuals.size() > 1) {
//                
//                int parentPosition2;
//                do {
//                    parentPosition2 = arc4random_uniform((int)individuals.size());
//                } while (parentPosition2 == parentPosition);
//                p2 = SystemInfo<IndividualType>(individuals[parentPosition2]);
//            }
//        
//        IndividualType child = combineSystems(p1.individual, p2.individual);
//        // mutate the new individual (maybe)
//        mutateSystem(child);
//        
//        SystemInfo<IndividualType> i = SystemInfo<IndividualType>(child);
//        newGeneration.push_back(i);
//    } while (newGeneration.size() < populationSize);
//    
//    // add the original population too
//   population.insert(population.end(), newGeneration.begin(), newGeneration.end());
//   // population = newGeneration;
//}

template <class IndividualType, class InnovationType>
NEATSpecies<IndividualType> &NEATPlusAlgorithm<IndividualType, InnovationType>::chooseBreedingSpecies(double totalFitness)
{
    
    // choose a species according to species fitness
    // then chose parents at random from within that species
    
    // Assume sorted population
    double selector = ((double )rand()/RAND_MAX);
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator it = speciesList.begin();
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

//template <class IndividualType, class InnovationType>
//std::pair<SystemInfo<IndividualType>, SystemInfo<IndividualType> > NEATPlusAlgorithm <IndividualType, InnovationType>::selectParents(double  fitnessSum)
//{
//    NEATSpecies<IndividualType> breedingSpecies = chooseBreedingSpecies(fitnessSum);
//  
//        std::vector<SystemInfo<IndividualType> >individuals = breedingSpecies.members;
//        
//        SystemInfo<IndividualType> parent1 = (individuals[0]);
//        SystemInfo<IndividualType> parent2 = (individuals[0]);
//    
//         int parentPosition = arc4random_uniform(individuals.size());
//            parent1 = SystemInfo<IndividualType>(individuals[parentPosition]);
//    
//    double interspecies_mate = (double)rand()/RAND_MAX;
//    if (stagnantGenerations > maxStagnation/2 || interspecies_mate < 0.05) {
//        breedingSpecies = chooseBreedingSpecies(fitnessSum);
//        individuals = breedingSpecies.members;
//
//        int parentPosition = arc4random_uniform(individuals.size());
//        parent2 = SystemInfo<IndividualType>(individuals[parentPosition]);
//    } else
//        if (individuals.size() > 1) {
//   
//            int parentPosition2;
//            do {
//                parentPosition2 = arc4random_uniform(individuals.size());
//            } while (parentPosition2 == parentPosition);
//            parent2 = SystemInfo<IndividualType>(individuals[parentPosition2]);
//        }
//
//        return std::pair<SystemInfo<IndividualType> , SystemInfo<IndividualType> >(parent1, parent2);
//    
//    
//}


template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::speciate()
{
    typename std::vector<SystemInfo<IndividualType> >::iterator populationIter = population.begin();
    
    while (populationIter != population.end()) {
        // assign to species
        bool added= false;
        
        typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
        while (speciesIterator != speciesList.end()) {
            double  distance = genomeDistance(populationIter->individual, *speciesIterator->representative);
            if (fabs(distance) < d_threshold) {
                added = true;
                speciesIterator->members.push_back(&(*populationIter));
                break;
            }
            speciesIterator++;
        }
        
        if (!added) {
            // new species!!
            NEATSpecies<IndividualType> s = NEATSpecies<IndividualType>();
            s.representative = new IndividualType(populationIter->individual);
            s.members.push_back(&(*populationIter));
            speciesList.push_back(s);
            s.speciesNumber = nextSpeciesNumber++;
        }
        
        populationIter++;
    }
    
   
}

template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::updateSharedFitnesses()
{
    long nSpecies = speciesList.size();
    fprintf(stderr, "%ld species\n", nSpecies);
    
    // for each species, copy a random current member to be the representative for the next generation, and adjust the fitnesses for sharing
    // kill off a species with no members
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    while (speciesIterator != speciesList.end()) {
        // pick a new rep
        std::vector<SystemInfo<IndividualType> *> members = speciesIterator->members;
        if (members.size() == 0) {
            NEATSpecies<IndividualType> extinctSpecies = *speciesIterator;
            speciesIterator = speciesList.erase(speciesIterator);
        } else {
            int index = arc4random_uniform(members.size());
            IndividualType rep = (members[index]->individual);
            delete speciesIterator->representative;
            speciesIterator->representative = new IndividualType(rep);
            
            // fprintf(stderr, "\tSpecies %d: %ld members\n", ++i, members.size());
            
            // update the fitnesses
            double  totalFitness = 0.0;
            typename std::vector<SystemInfo<IndividualType> *>::iterator memberIterator = members.begin();
            while (memberIterator != members.end()) {
                (*memberIterator)->fitness /= members.size();
                totalFitness += (*memberIterator)->fitness;
                memberIterator++;
            }
            if (totalFitness < 0) {
                fprintf(stderr, "WOAH");
            }
            
            speciesIterator->totalSharedFitness = totalFitness;
            speciesIterator++;
        }
    }
}

// unlike base NEAT, also take node differences into account
template <class IndividualType, class InnovationType>
double  NEATPlusAlgorithm <IndividualType, InnovationType>::genomeDistance( IndividualType& sys1,  IndividualType& sys2)
{
    std::vector<InnovationType> genome1 = sys1.connectionGenome();
    std::vector<InnovationType> genome2 = sys2.connectionGenome();
    
    // find the longer one
    size_t longerSize = (genome1.size() > genome2.size() ? genome1.size()  : genome2.size());
    
    size_t nDisjoint = 0;
    size_t nExcess = 0;
    size_t nMatching = 0;
    size_t matchingDiff = 0;
    
    // yay stdlib!
    std::vector<InnovationType> matchingGenes1;
    std::vector<InnovationType> matchingGenes2;
    std::vector<InnovationType> disjointFromSys1;
    std::vector<InnovationType> disjointFromSys2;
    
    // set intersection takes from the first range given - I do it twice so I have the matches from parent 1 and from parent 2.
    std::set_intersection(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(matchingGenes1), compareInnovationNumbers<InnovationType>);
    std::set_intersection(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(),std::back_inserter(matchingGenes2), compareInnovationNumbers<InnovationType>);
    
    // difference takes things in the first range that are not in the second - we will have to check for excess ourself
    std::set_difference(genome1.begin(), genome1.end(), genome2.begin(), genome2.end(), std::back_inserter(disjointFromSys1), compareInnovationNumbers<InnovationType>);
    std::set_difference(genome2.begin(), genome2.end(), genome1.begin(), genome1.end(), std::back_inserter(disjointFromSys2), compareInnovationNumbers<InnovationType>);
    
    
    nMatching = matchingGenes2.size();
    // first find the distance between matching genes
    for (int i=0; i<nMatching; i++) {
        InnovationType a1 = matchingGenes1[i];
        InnovationType a2 = matchingGenes2[i];
        
        matchingDiff += IndividualType::connectionDifference(a1, a2);
    }
    
    // now determine the excess vs disjoint
    // if one of the disjoint sets is empty, then it's all excess
    if (disjointFromSys2.size() == 0 || disjointFromSys1.size() == 0) {
        nExcess = disjointFromSys1.size() + disjointFromSys2.size();
    } else {
        // else chalk everything up to disjoint
        nDisjoint = disjointFromSys1.size() + disjointFromSys2.size();
        
        // then find the number of excess genes
        int last1 = disjointFromSys1.back().innovationNumber;
        int last2 = disjointFromSys2.back().innovationNumber;
        
        std::vector<InnovationType> shorter;
        std::vector<InnovationType> longer;
        
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
                if (longer[i].innovationNumber < shorter.back().innovationNumber)
                    break;
            }
            nExcess = i;
            nDisjoint -= nExcess;
        }
    }
    
    
    double diff_node = sys1.nodeDifference(sys2);
    long moreNodes = std::min(sys1.numberOfNodes(), sys2.numberOfNodes());
    
    double  d = w_disjoint*nDisjoint/longerSize + w_excess*nExcess/longerSize + w_matching*matchingDiff/nMatching;

//    double  d = w_disjoint*nDisjoint + w_excess*nExcess + w_matching*matchingDiff;
    
  //  d+= + w_matching_node*diff_node/moreNodes;
    
    //assert(!isUnreasonable(d));
    return d;
}


template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::prepareInitialPopulation()
{
    // mutate the initial system to get an initial population
    while (population.size() < populationSize) {
        IndividualType newSystem = createInitialIndividual();
        mutateSystem(newSystem);
        SystemInfo<IndividualType> info = SystemInfo<IndividualType>(newSystem);
        population.push_back(info);
    }
}


// descending, not ascending, order
template <class IndividualType>
static bool compareSpeciesFitness(const NEATSpecies<IndividualType>& s1, const NEATSpecies<IndividualType>& s2)
{
    return s1.totalSharedFitness > s2.totalSharedFitness;
}

template <class IndividualType, class InnovationType>
bool NEATPlusAlgorithm <IndividualType, InnovationType>::tick()
{
    bool first_run = false;
    if (population.size() == 0) {
        prepareInitialPopulation();
        
        origin = new IndividualType(population.front().individual);
        first_run = true;
    }
    
    double  bestFitness = -INFINITY; // we want zero fitness
    
    for (size_t popIter = 0; popIter <population.size(); popIter++) {
        population[popIter].fitness = evaluationFunc(population[popIter].individual);
    }
    
    if (first_run)
        speciate(); // later on we stay in our parents' species! OOOOH
    
    updateSharedFitnesses();
    
    // figure out how many of each species to save
    double  sharedFitnessSum = 0.0;
    for (int i=0; i<speciesList.size(); i++) {
        sharedFitnessSum += speciesList[i].totalSharedFitness;
    }
    
    std::vector<SystemInfo<IndividualType> > individualsToSave;
    // sort the species to ensure we keep the best
    std::sort(speciesList.begin(), speciesList.end(), compareSpeciesFitness<IndividualType>);
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIterator = speciesList.begin();
    
    while (speciesIterator != speciesList.end()) {
        size_t numMembers = speciesIterator->members.size();
        
        double  proportionToSave = (speciesIterator->totalSharedFitness)/sharedFitnessSum;
        int numToSave = std::min((int)(proportionToSave*populationSize), (int)numMembers);
      //  numToSave = std::max(numToSave, 1); // save 1 of each species
        
        if (numToSave < 0) {
            fprintf(stderr, "WOAH");
        }

        std::sort(speciesIterator->members.begin(), speciesIterator->members.end(), compareIndividuals<IndividualType>);
        
        // don't let population grow unchecked
        if (individualsToSave.size() > populationSize)
            numToSave = 0;
        
        // fprintf(stderr, "Species %d: %d/%ld\n", i++, numToSave,  (*speciesIterator)->members.size());
        
        
        for (int j=0; j<numToSave; j++) {
            SystemInfo<IndividualType> *individual = speciesIterator->members[j];
            double  rawFitness = individual->fitness * numMembers;
            if (rawFitness > bestFitness) {
                bestFitness = rawFitness;
                _bestIndividual = new IndividualType(individual->individual);
            }
            individualsToSave.push_back(*individual);
        }
        
        typename std::vector<SystemInfo<IndividualType> *>::iterator deleteIt =  speciesIterator->members.begin()+numToSave;
        while (deleteIt !=  speciesIterator->members.end()) {
            deleteIt =  speciesIterator->members.erase(deleteIt);
        }
        
        if (numToSave == 0)
            speciesIterator = speciesList.erase(speciesIterator);
        else
            speciesIterator++;
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
    
    bool stop =  (generations >= maxGenerations) || (stopFunc && stopFunc(bestFitness));
    logPopulationStatistics();
    
    // empty the innovation list before spawning more
  //  newConnections.clear();
    fprintf(stderr, "BEST FITNESS: %f\n", bestFitness);
    if (stop) {
        fprintf(stderr, "ALL TIME BEST FITNESS: %f\n", allTimeBestFitness);
    } else {
        spawnNextGeneration();
    }
    
    bool regroup = false;
    if (speciesList.size() < 10) {
        regroup = true;
        d_threshold /= 1.05;
    }
    
    if (speciesList.size() > 20) {
        regroup = true;
        d_threshold *= 1.05;
    }
    
    if (regroup) {
        // empty the species lists
        speciesList.clear();
        speciate();
    }
    

    
    
    return  stop;
}

template <class IndividualType, class InnovationType>
IndividualType  *NEATPlusAlgorithm <IndividualType, InnovationType>::bestIndividual()
{
    return _bestIndividual;
}

template <class IndividualType, class InnovationType>
long NEATPlusAlgorithm <IndividualType, InnovationType>::getNumberOfIterations()
{
    return generations;
}


template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::mutateSystem(IndividualType& original)
{
    
    std::vector<InnovationType> allAttachments = original.connectionGenome();
    typename std::vector<InnovationType>::iterator it = allAttachments.begin();
    while (it != allAttachments.end()) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_conn)
            original.mutateConnectionWeight();
        it++;
    }
    
    long nNodes = original.numberOfNodes();
    for (long i=0; i<nNodes; i++) {
        double selector = (double )rand()/RAND_MAX;
        if (selector < p_m_node)
            original.mutateNode(i);
    }
        
    // add attachment or add node - or both
    double selector1 = (double )rand()/RAND_MAX;
    if (selector1 < p_m_conn_ins) {
        // insert a new connection
        // update the innovation number
        InnovationType created = original.createConnection();
        assignInnovationNumberToAttachment(original, created);
    }
    
    double selector2 = (double )rand()/RAND_MAX;
    if (selector2 < p_m_node_ins) {
        // add a node somewhere if possible
        std::vector<InnovationType> new_edges = original.insertNode();
        for (int i=0; i<new_edges.size(); i++) {
            assignInnovationNumberToAttachment(original,  new_edges[i]);
        }
    }
    
}

template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::assignInnovationNumberToAttachment(IndividualType& individual, InnovationType i){
    // have we already created this 'innovation' in this generation?
    typename std::vector<InnovationType>::iterator it = std::find(newConnections.begin(), newConnections.end(), i);
    if (it!= newConnections.end()) {
        i.innovationNumber = it->innovationNumber;
    } else {
        i.innovationNumber = nextInnovationNumber++;
        newConnections.push_back(i);
    }
    individual.updateInnovationNumber(i);
}

template <class IndividualType, class InnovationType>
void NEATPlusAlgorithm <IndividualType, InnovationType>::logPopulationStatistics()
{
    if (!currentLogFile.is_open()) {
        time_t now = time(NULL);
        std::stringstream s;
        wordexp_t directory;
        wordexp("~/temp/neatNN/", &directory, 0);
        s << directory.we_wordv[0];
        s << "neat-log" << now << ".log";
        std::string filename = s.str();
        currentLogFile.open(filename.c_str());
    }
    
    // 1 species per line : # in species, best fitness, worst fitness, mean fitness, std. dev, distance to starting genome
    // then an empty line
    
    typename std::vector<NEATSpecies<IndividualType> >::iterator speciesIter = speciesList.begin();
    
    while (speciesIter != speciesList.end()) {
        std::vector<SystemInfo<IndividualType> *>members = speciesIter->members;
        size_t n = members.size();
        currentLogFile << speciesIter->speciesNumber << " ";
        currentLogFile << n << " ";
        SystemInfo<IndividualType> *best = members.front();
        currentLogFile << best->fitness*n << " ";
        SystemInfo<IndividualType> *worst = members.back();
        currentLogFile << worst->fitness*n << " ";
        double  meanFitness =speciesIter->totalSharedFitness;
        currentLogFile << meanFitness << " ";
        // standard dev - from http://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
        
        double  sq_sum = 0;
        typename std::vector<SystemInfo<IndividualType> *>::iterator it = members.begin();
        while (it != members.end()) {
            sq_sum += ((*it)->fitness - meanFitness)*((*it)->fitness - meanFitness);
            it++;
        }
        
        double  stdev = sqrt(sq_sum /n);
        currentLogFile << stdev << " ";
        
        double  speciesDist = genomeDistance(*speciesIter->representative, *origin);
        currentLogFile << speciesDist << "\n";
        speciesIter++;
        
    }
    currentLogFile << "\n";
    currentLogFile.flush();
}


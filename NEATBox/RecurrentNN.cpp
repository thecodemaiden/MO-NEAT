//
//  BaseNetwork.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "RecurrentNN.h"
#include <algorithm>
#include <sstream>
#include <cmath>
#include <assert.h>

// locate incoming and outgoing edges for a node number
struct has_input {
    has_input(long n) : n(n) { }
    bool operator()(DelayEdge e) const { return e.nodeTo == n && !e.disabled; }
private:
    long n;
};

struct has_output {
    has_output(long n) : n(n) { }
    bool operator()(DelayEdge e) const { return e.nodeFrom == n && !e.disabled;}
private:
    long n;
};

// need a way to detect cycles so that delays can be enforced in cycles

RecurrentNN::RecurrentNN(int nInputs, int nOutputs)
:MNIndividual(), CycledNN()
{
    for (int i=0; i<nInputs; i++) {
        nodes.push_back(Node());
        inputNodes.insert(i);
    }
    
    // connect every input to every output
    for (int j=0; j<nOutputs; j++) {
        Node n = Node();
        int nextNode = (int)nodes.size();
        outputNodes.insert(nextNode);
        
        for (int i=0; i<nInputs; i++) {
            DelayEdge e = DelayEdge(i, nextNode);
            nodes[i].outdegree++;
            n.indegree++;
            edges.push_back(e);
        }
        nodes.push_back(n);
    }
    
}

void RecurrentNN::addGeneFromParentSystem(MNIndividual *p, MNEdge *g)
{
    
    RecurrentNN *parent = dynamic_cast<RecurrentNN *>(p);
    DelayEdge *gene = dynamic_cast<DelayEdge *>(g);
    
    if (!parent || !gene)
        return; // can't do nothin with this
    
    if (gene->nodeFrom < 0 || gene->nodeTo < 0)
        return; // invalid
    
    std::vector<DelayEdge>::iterator found = std::find(edges.begin(), edges.end(), *gene);
    
    bool didDisable = false;
    bool wasDisabled = false;
    
    
    if (found != edges.end()) {
        found->weight = gene->weight;
        wasDisabled = found->disabled;
        if (found->disabled == gene->disabled)
            return; // nothing more to do, don't modify in/outdegrees
        found->disabled = gene->disabled;
        didDisable = gene->disabled;
    } else {
        wasDisabled = true;
        edges.push_back(*gene);
        
        while (gene->nodeFrom >= nodes.size() || gene->nodeTo >= nodes.size()) {
            // add nodes up to the required number, copying from the parent
            // hope that we don't have floating nodes...
            Node newNode = Node(parent->nodes[nodes.size()]);
            newNode.indegree = newNode.outdegree = 0;
            nodes.push_back(newNode);
        }
    }
    if (!gene->disabled) {
        nodes.at(gene->nodeFrom).outdegree++;
        nodes.at(gene->nodeTo).indegree++;
    }
    
    if (didDisable) {
        nodes.at(gene->nodeFrom).outdegree--;
        nodes.at(gene->nodeTo).indegree--;
    }
}

std::vector<MNEdge *> RecurrentNN::connectionGenome()
{
    std::vector<MNEdge *> genome;
    for (std::vector<DelayEdge>::iterator it = edges.begin(); it != edges.end(); it++) {
        genome.push_back(&*it);
    }
    return genome;
}

void RecurrentNN::mutateConnectionWeights(double p_m)
{
    // pick a random edge (disabled is fine?) and mutate its weight
    for (size_t pos = 0; pos < edges.size(); pos++) {
        if (uniformProbability() < p_m) {
            double var = normallyDistributed();
            DelayEdge &toMutate = edges.at(pos);
            
            toMutate.weight = var;
        }
    }
}

void RecurrentNN::mutateNodes(double p_m)
{
    double selector = (double)rand()/RAND_MAX;
    
    for (size_t i = 0; i<nodes.size(); i++) {
        Node &n  = nodes[i];
        if (uniformProbability() < p_m) {
            if (selector < 0.33)
                n.bias = normallyDistributed(n.bias, 1);
            // else if (selector < 0.67)
            //   n.param1 = normallyDistributed(n.param1, 1);
            else
                n.type = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
        }
    }
}


std::vector<MNEdge *> RecurrentNN::createNode()
{
    // if all edges are disabled, that's BAD and we are uselesssssss
    long activeEdges = std::count_if(edges.begin(), edges.end(), [](DelayEdge e) {return !e.disabled;});
    if (activeEdges == 0)
        return std::vector<MNEdge *>();
    
    uint pos;
    do {
        pos = arc4random_uniform((uint)edges.size());
        
    } while (edges.at(pos).disabled);
    
    return insertNodeOnEdge(pos);
}

std::vector<MNEdge *> RecurrentNN::insertNodeOnEdge(long ePos)
{
    
    uint nextNode = (uint)nodes.size();
    Node newNode = Node();
    
    bool wasDisabled = edges[ePos].disabled;
    
    if (!wasDisabled) {
        newNode.indegree = newNode.outdegree = 1;
    }
    
    nodes.push_back(newNode);
    
    DelayEdge e1 = DelayEdge(edges[ePos].nodeFrom, nextNode);
    DelayEdge e2 = DelayEdge(nextNode, edges[ePos].nodeTo);
    
    e1.disabled = edges[ePos].disabled;
    e2.disabled = edges[ePos].disabled;
    edges[ePos].disabled = true;
    
    
    edges.push_back(e1);
    DelayEdge *first = &edges.back();
    
    edges.push_back(e2);
    DelayEdge *second = &edges.back();
    
    std::vector<MNEdge *> toReturn;
    
    toReturn.push_back(first);
    toReturn.push_back(second);
    return toReturn;
}

std::vector<std::vector<double> > RecurrentNN::simulateSequence(const std::vector<std::vector<double> > &inputValues, int delay)
{
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<std::vector<double> > allOuputs;
    
    std::vector<std::vector<double> > memory;
    
    for (long i=0; i<inputValues.size(); i++) {
        std::vector<double> currentOutputs = std::vector<double>(nNodes);
        std::vector<double> currentInputs = inputValues[i];
        for (long steps = -1; steps<delay; steps++) {
            std::vector<double> lastOutputs = currentOutputs;
            currentOutputs = nodeOutputsForInputs(currentInputs, currentOutputs, memory);
            if (lastOutputs == currentOutputs) {
                break;
            }
        }
        
        std::vector<double> finalOutputs;
        for (std::set<long>::iterator it = outputNodes.begin(); it != outputNodes.end(); it++) {
            long node = *it;
            finalOutputs.push_back(currentOutputs[node]);
        }
        allOuputs.push_back(finalOutputs);
    }
    
    return allOuputs;
}


DelayEdge *RecurrentNN::createConnection()
{
    // pick two existing, unconnected nodes and connect them - can connect self to self!
    long maxConnections = nodes.size()* nodes.size();
    
    long activeEdges = std::count_if(edges.begin(), edges.end(), [](DelayEdge e) {return !e.disabled;});
    
    DelayEdge newEdge = DelayEdge(-1,-1);
    
    std::vector<DelayEdge>::iterator edgePos = edges.end();
    
    DelayEdge *toReturn = NULL;
    
    if (activeEdges < maxConnections) {
        // chose a source node at random
        long maxDegree = nodes.size();
        Node *source = NULL;
        uint pos;
        do {
            pos = arc4random_uniform((uint)nodes.size());
            source = &nodes.at(pos);
        } while (!source || source->outdegree >= maxDegree);
        newEdge.nodeFrom = pos;
        nodes.at(pos).outdegree++;
        
        // choose one it's not connected to yet -
        // reduce the search space
        std::vector<DelayEdge> existingEdges;
        
        std::vector<DelayEdge>::iterator it = edges.begin();
        while (it != edges.end()) {
            if (it->nodeFrom == newEdge.nodeFrom)
                existingEdges.push_back(*it);
            it++;
        }
        
        Node *sink = NULL;
        bool found = false;
        
        do {
            pos = arc4random_uniform((uint)nodes.size());
            
            sink = &nodes.at(pos);
            newEdge.nodeTo = pos;
            edgePos = std::find(existingEdges.begin(), existingEdges.end(), newEdge);
            
            found = (edgePos != existingEdges.end() && !edgePos->disabled); // already enabled edge?
            
            
        } while (!sink || sink->indegree >= maxDegree || found);
        nodes.at(pos).indegree++;
        
        if (edgePos != existingEdges.end() && edgePos->disabled) {
            edgePos = std::find(edges.begin(), edges.end(), newEdge);
            edgePos->disabled = false;
        } else {
            edges.push_back(newEdge);
            edgePos = edges.end()-1;
        }
        toReturn = &(*edgePos);
        
        // if this created a cycle with no delay, put a delay onto this edge
        std::vector<std::vector<long> > cycles = cycleSort();
        
        // determines if the edge is in a cycle with no delays
        std::vector<long> cycle;
        for (int i=0; i<cycles.size(); i++) {
            std::vector<long>::const_iterator nodeFromPos;
            std::vector<long>::const_iterator nodeToPos;
            
            nodeFromPos = std::find(cycles[i].begin(), cycles[i].end(), edgePos->nodeTo);
            nodeToPos = std::find(cycles[i].begin(), cycles[i].end(), edgePos->nodeFrom);
            
            if (nodeFromPos != cycles[i].end() && nodeToPos != cycles[i].end()) {
                cycle = cycles[i];
                break;
            }
        }
        
        
        bool tightCycle = false;
        if (cycle.size() > 1) {
            std::sort(cycle.begin(), cycle.end());
            tightCycle = true;
            // check that at least one edge in the cycle has delay >= 1
            std::vector<DelayEdge>::iterator it = edges.begin();
            do {
                it = std::find_if(it, edges.end(), [cycle](const DelayEdge &e){
                    bool fromPresent = std::binary_search(cycle.begin(), cycle.end(), e.nodeFrom);
                    bool toPresent = std::binary_search(cycle.begin(), cycle.end(), e.nodeTo);
                    if (fromPresent && toPresent) {
                        // this is an edge cycle
                        return true;
                    }
                    return false;
                });
                if (it != edges.end()) {
                    if (it->delay > 0) {
                        tightCycle = false;
                        // we only need one delayed edge
                        break;
                    }
                    ++it;
                } else {
                    break;
                }
            } while (1);
        }
        
        if (tightCycle) {
            toReturn->delay = 1;
        }
    }
    
    
    
    return toReturn;
}

double RecurrentNN::connectionDifference(MNEdge *e1, MNEdge *e2)
{
    
    DelayEdge *first = dynamic_cast<DelayEdge *>(e1);
    DelayEdge *second = dynamic_cast<DelayEdge *>(e2);
    
    if (!first || !second)
        return INFINITY;
    
    return fabs(first->weight - second->weight);
}

double RecurrentNN::nodeDifference(MNIndividual *ind)
{
    RecurrentNN *other = dynamic_cast<RecurrentNN *>(ind);
    if (!other)
        return INFINITY;
    
    
    
    long length = std::max(other->nodes.size(), nodes.size());
    double d = 0;
    for (long i=0; i<length; i++) {
        
        double d11 = 0;
        double d21 = 0;
        
        double d12 = 0;
        double d22 = 0;
        
        if (i < nodes.size()) {
            d11 = applyActivationFunc(nodes[i], 0) - applyActivationFunc(nodes[i], 1);
            d21 = applyActivationFunc(nodes[i], -1) - applyActivationFunc(nodes[i], 0);
        }
        if (i < other->nodes.size()) {
            d12 = applyActivationFunc(other->nodes[i], 0) - applyActivationFunc(other->nodes[i], 1);
            d22 = applyActivationFunc(other->nodes[i], -1) - applyActivationFunc(other->nodes[i], 0);
        }
        
        d += fabs((d11 - d21)/2 - (d12 - d22)/2);
    }
    return d;
}

long RecurrentNN::numberOfNodes()
{
    return nodes.size();
}

long RecurrentNN::numberOfEdges()
{
    return edges.size();
}


double RecurrentNN::visitNode(long i, std::set<long> &visitedNodes, std::vector<double> &lastOutputs)
{
    Node n = nodes[i];
    visitedNodes.insert(i);
    std::vector<DelayEdge> inputEdges = inputsToNode(i);
    
    double inputSum = n.bias;
    
    for (long j=0; j<inputEdges.size(); j++) {
        long nodeN = inputEdges[j].nodeFrom;
        if (visitedNodes.find(nodeN) == visitedNodes.end()) {
            lastOutputs[nodeN] = visitNode(nodeN, visitedNodes, lastOutputs);
        }
        double rawVal = lastOutputs[nodeN];
        
        inputSum += inputEdges[j].weight*rawVal;
        
    }
    
    return inputSum;
}

std::vector<double> RecurrentNN::nodeOutputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs, std::vector<std::vector<double> > &memory)
{
    std::vector<double> newOutputs;
    
    long inputNodeN = 0; // used to map values from the inputs vector to input nodes
    std::set<long>visitedNodes;
    for (long i=0; i<nodes.size(); i++) {
        // collect the input values
        double inputSum = visitNode(i, visitedNodes, lastOutputs);
        
        bool isInput = (inputNodes.find(i) != inputNodes.end());
        
        if (isInput) {
            inputSum += (inputs[inputNodeN++]);
        }
        
        double output = applyActivationFunc(nodes[i], inputSum);
        assert(!isUnreasonable(output));
        
        newOutputs.push_back(output);
    }
    return newOutputs;
}



std::vector<DelayEdge> RecurrentNN::inputsToNode(long n)
{
    
    std::vector<DelayEdge> found;
    
    has_input pred = has_input(n);
    
    std::vector<DelayEdge>::iterator it = std::find_if(edges.begin(), edges.end(), pred);
    for (; it != edges.end(); it = std::find_if(++it, edges.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

std::vector<DelayEdge> RecurrentNN::outputsFromNode(long n)
{
    
    std::vector<DelayEdge> found;
    
    has_output pred = has_output(n);
    
    std::vector<DelayEdge>::iterator it = std::find_if(edges.begin(), edges.end(), pred);
    for (; it != edges.end(); it = std::find_if(++it, edges.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

std::string RecurrentNN::display()
{
    // sort the edges by source node and display
    
    std::ostringstream ss;
    std::sort(edges.begin(), edges.end());
    
    std::vector<DelayEdge>::iterator it = edges.begin();
    for (; it != edges.end(); it++) {
        ss << it->nodeFrom << " -> " << it->nodeTo << ": " << it->weight;
        if (it->disabled)
            ss << " (disabled)";
        ss << "\n";
    }
    
    // display the nodes too
    for (int i=0; i<nodes.size(); i++) {
        Node n = nodes[i];
        ss << i << ": " << activationFuncName(n.type);
        ss << " P: " << n.param1;
        ss << " B: " << n.bias << "\n";
    }
    
    return ss.str();
}

std::string RecurrentNN::dotFormat(std::string graphName)
{
    std::ostringstream ss;
    
    ss << "digraph " << graphName << "{\n" ;
    ss << "size=\"7.75,10.75\";\n";
    ss << "\trankdir=LR;\n";
    ss << "\tcenter=true;\n";
    ss << "\torientation=landscape;\n";
    
    ss << "{rank=source;\n";
    int i = 0;
    for (std::set<long>::iterator it = inputNodes.begin(); it != inputNodes.end(); it++, i++) {
        ss << "\tin_" << i <<";\n";
    }
    ss << "}\n\n";
    
    i=0;
    for (std::set<long>::iterator it = inputNodes.begin(); it != inputNodes.end(); it++, i++) {
        ss << "\tin_" << i << " -> " << *it <<";\n";
    }
    
    ss << "{rank=sink;\n";
    i = 0;
    for (std::set<long>::iterator it = outputNodes.begin(); it != outputNodes.end(); it++, i++) {
        ss << "\t" << "out_" << i << ";\n";
    }
    ss <<"}\n\n";
    
    i=0;
    for (std::set<long>::iterator it = outputNodes.begin(); it != outputNodes.end(); it++, i++) {
        ss << "\t" << *it << " -> " << "out_" << i << ";\n";
    }
    
    for (long j = 0; j<nodes.size(); j++) {
        Node n = nodes[j];
        ss << "\t" << j << "[label = \"" << activationFuncName(nodes[j].type) << "\\n" << n.param1 << "\"];\n";
        if (n.bias != 0) {
            ss << "\tb_" << j << " -> " << j << " [label = \"" << n.bias <<"\n" <<"\"];\n";
        }
        ss << "\n";
    }
    
    for (long j=0; j<edges.size(); j++) {
        DelayEdge e = edges[j];
        ss << "\t" << e.nodeFrom << " -> " << e.nodeTo << " [label =\"" << e.weight;
        if (e.disabled) {
            ss << ", style=dotted";
        } else {
            ss <<  "\n" << e.delay;
        }
        ss << "\"];\n";
    }
    
    ss <<"}\n";
    
    
    return ss.str();
}

// helper struct and function for cycleSort
struct RecurrentNN::CycleNode {
    long index;
    long lowlink;
    long node;
    CycleNode(long n, long i=-1, long l=INFINITY)
     :node(n), index(i),lowlink(l){}
};

void RecurrentNN::strongConnect(long v, long index, std::vector<long> &stack,  std::vector<RecurrentNN::CycleNode> &nodes,  std::vector<std::vector<long> > &components)
{
    std::vector<long> connectedComponent;
    
    CycleNode &cn = nodes[v];
    cn.index = index;
    cn.lowlink = index;
    index = index+1;
    stack.push_back(v);
    
    std::vector<DelayEdge> succesors = outputsFromNode(v);
    for (long i=0; i<succesors.size(); i++) {
        long toNode = succesors[i].nodeTo;
        CycleNode &nextN = nodes[toNode];
        if (nextN.index < 0) {
            strongConnect(toNode, index, stack, nodes, components);
            cn.lowlink = std::min(nextN.lowlink, cn.lowlink);
        } else if (std::find(stack.begin(), stack.end(), toNode) != stack.end()) {
            // Successor w is in stack S and hence in the current SCC
            cn.lowlink = std::min(nextN.index, cn.lowlink);
        }
    }
    
    if (cn.lowlink == cn.index) {
        long nodeFrom = v;
        long nodeTo = -1;
        do {
            nodeTo = stack.back();
            stack.pop_back();
            connectedComponent.push_back(nodeTo);
        } while (nodeFrom != nodeTo);
        components.push_back(connectedComponent);
    }
}

//http://en.wikipedia.org/wiki/Tarjan’s_strongly_connected_components_algorithm
std::vector<std::vector<long> > RecurrentNN::cycleSort()
{
    std::vector<std::vector<long> > sorted;
    std::vector<long> stack;
    
    std::vector<RecurrentNN::CycleNode> cycleNodes;

    long index = 0;
    // make a CycleNode for each actual node
    for (long i=0; i<nodes.size(); i++) {
        CycleNode cn = CycleNode(i);
        cycleNodes.push_back(cn);
    }
    
    for (long i=0; i<cycleNodes.size(); i++) {
        if (cycleNodes[i].index < 0) {
            strongConnect(i, index, stack, cycleNodes, sorted);
        }
    }
    
    
    return sorted;
}

//void RecurrentNN::cleanup()
//{
//    // in, out
//    std::vector<std::pair<int, int> > nodeDegrees(nodes.size());
//    for (long i=0; i< edges.size(); i++) {
//        if (!edges[i].disabled) {
//            nodeDegrees[edges[i].nodeFrom].second++;
//            nodeDegrees[edges[i].nodeTo].first++;
//        }
//    }
//    
//    for (long i=0; i<nodes.size(); i++) {
//        assert(nodes[i].indegree == nodeDegrees[i].first);
//        assert(nodes[i].outdegree == nodeDegrees[i].second);
//    }
//    
//}

//void RecurrentNN::cleanup()
//{
//    bool done;
//    do {
//       done = true;
//        std::set<long> disabledNodes;
//        for (int i=0; i< nodes.size(); i++) {
//            Node &n = nodes[i];
//            // if our indegree is 0, disable our inputs
//            if (inputNodes.count(i) == 0 && n.indegree == 0 && n.outdegree > 0) {
//                disabledNodes.insert(i);
//            }
//        }
//        if (!done) {
//            for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); it++) {
//                if (disabledNodes.count(it->nodeFrom) > 0) {
//                    if (!it->disabled) {
//                       done=false;
//                        it->disabled = true;
//                       nodes[it->nodeTo].indegree--;
//                        nodes[it->nodeFrom].outdegree--;
//                    }
//                }
//            }
//        }
//    }  while (!done);
//}
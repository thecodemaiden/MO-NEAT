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
#include <deque>
#include <set>

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
:MNIndividual(), CycledNN(),
maxDelay(0)
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
        if (found->nodeFrom == found->nodeTo)
            assert(found->delay > 0);
    } else {
        wasDisabled = true;
        edges.push_back(*gene);
        found = edges.end()-1;
        
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
 //   fixCycles(*found);
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
    
    // put the delay on the first edge
    e1.delay = edges[ePos].delay;
    
    edges.push_back(e1);
    DelayEdge *first = &edges.back();
    
    edges.push_back(e2);
    DelayEdge *second = &edges.back();
    
    std::vector<MNEdge *> toReturn;
    
    toReturn.push_back(first);
    toReturn.push_back(second);
    return toReturn;
}

std::vector<std::vector<double> > RecurrentNN::simulateSequence(const std::vector<std::vector<double> > &inputValues)
{
    fixCycles(edges.front());
    
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<std::vector<double> > allOuputs;
    
    std::vector<std::vector<double> > memory(nNodes);
    
    for (long i=0; i<inputValues.size(); i++) {
        std::vector<double> currentOutputs = std::vector<double>(nNodes);
        std::vector<double> currentInputs = inputValues[i];
        for (long steps = -1; steps<maxDelay; steps++) {
            std::vector<double> lastOutputs = currentOutputs;
            currentOutputs = nodeOutputsForInputs(currentInputs, memory);
            if (lastOutputs == currentOutputs) {
                // settled
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



std::vector<double> RecurrentNN::nodeOutputsForInputs(std::vector<double> inputs, std::vector<std::vector<double> > &memory)
{
    std::vector<double> newOutputs(nodes.size());
    
    std::deque<long> toVisitNodes;
    for (int i=0; i<nodes.size(); i++)
        toVisitNodes.push_back(i);
    
    std::set<long> resolvedNodes;
    
    int runs = 0;
    
    while (toVisitNodes.size() > 0) {
        long n = toVisitNodes.front();
        
        std::vector<DelayEdge> inputEdges = inputsToNode(n);
        double inputSum = nodes[n].bias;
        bool tryAgain = false;
        for (long j=0; j<inputEdges.size(); j++) {
            DelayEdge &e = inputEdges[j];
            long from = e.nodeFrom;
            long delay = e.delay;
            std::vector<double> history = memory[from];
            
            bool isSelfInspection = false;
            if (from == n) {
                // we're looking at ourselves
                assert(delay > 0);
                isSelfInspection = true;
            }
            
            bool longEnoughHistory = (history.size() > delay);
            
            if (longEnoughHistory) {
                if (std::count(resolvedNodes.begin(), resolvedNodes.end(), from) > 0 || delay > 0) {
                    std::vector<double>::iterator valuePos = history.end() - (delay +1);
                    assert(!isUnreasonable(*valuePos));
                    inputSum += *valuePos;
                } else {
                    tryAgain = true;
                    break;
                }
            }
            
//            if (std::count(resolvedNodes.begin(), resolvedNodes.end(), from) > 0 && longEnoughHistory) {
//                // we've already resolved the output of this node
//                // the values can be taken from memory
//                std::vector<double>::iterator valuePos = history.end() - (delay+1);
//                inputSum += *valuePos;
//            } else {
//                if (longEnoughHistory) {
//                    // we should have an answer, but we don't - reschedule for later
//                    tryAgain = true;
//                    break;
//                }
//            }
        }
        
        if (!tryAgain) {
            // all input nodes resolved
            resolvedNodes.insert(n);
            
            
            bool isInput = (inputNodes.find(n) != inputNodes.end());
            
            if (isInput) {
                inputSum += (inputs[n]); // the first [inputNodes.size()] nodes are the inputs
            }
            
            double output = applyActivationFunc(nodes[n], inputSum);
            assert(!isUnreasonable(output));
            newOutputs[n] = output;
            std::vector<double> &thisMemory = memory[n];
            thisMemory.push_back(output);
        } else {
            // save this for later
            toVisitNodes.push_back(n);
        }
        toVisitNodes.pop_front();
        runs++;
    }
    
    return newOutputs;
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
        
    }
    
  //  fixCycles(*edgePos);
    
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


std::vector<DelayEdge> RecurrentNN::inputsToNode(long n, std::vector<DelayEdge> edgeSet)
{
    std::vector<DelayEdge> found;
    
    has_output pred = has_output(n);
    
    std::vector<DelayEdge>::iterator it = std::find_if(edgeSet.begin(), edgeSet.end(), pred);
    for (; it != edgeSet.end(); it = std::find_if(++it, edgeSet.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

std::vector<DelayEdge> RecurrentNN::outputsFromNode(long n, std::vector<DelayEdge> edgeSet)
{
    std::vector<DelayEdge> found;
    
    has_output pred = has_output(n);
    
    std::vector<DelayEdge>::iterator it = std::find_if(edgeSet.begin(), edgeSet.end(), pred);
    for (; it != edgeSet.end(); it = std::find_if(++it, edgeSet.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}


std::vector<DelayEdge> RecurrentNN::inputsToNode(long n)
{
    return inputsToNode(n, edges);
   
}

std::vector<DelayEdge> RecurrentNN::outputsFromNode(long n)
{
    return outputsFromNode(n, edges);
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
            ss << "\tb_" << j << " -> " << j << " [label = \"" << n.bias <<"\\n" <<"\"];\n";
        }
        ss << "\n";
    }
    
    for (long j=0; j<edges.size(); j++) {
        DelayEdge e = edges[j];
        ss << "\t" << e.nodeFrom << " -> " << e.nodeTo << " [label =\"" << e.weight;
        if (e.disabled) {
            ss << "\", style=dotted";
        } else {
            ss <<  "\\n" << e.delay << "\"";
        }
        ss << "];\n";
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

// may make this part of RecurrentNN
struct Path {
    std::vector<long> nodes;

public:
    Path(){}
    Path(long startNode) { nodes.push_back(startNode); }
    
    long totalDelay = 0;
    
    void extend(long next) {
        nodes.push_back(next);
    }
    
    bool containsNode(long n) const {
        return (std::find(nodes.begin(), nodes.end(), n) != nodes.end());
    }
    
    bool endsOnNode(long n)  const {
        return !nodes.empty() && nodes.back() == n;
    }
    
    bool containsEdge(long from, long to) {
        std::vector<long>::iterator fromPos = std::find(nodes.begin(), nodes.end(), from);
        std::vector<long>::iterator toPos = std::find(nodes.begin(), nodes.end(), to);
        
        if (fromPos != nodes.end() && toPos != nodes.end()) {
            return (toPos - fromPos == 1);
        }
            
        return false;
    }
    
    
};



void RecurrentNN::fixCycles(DelayEdge &e)
{
   // collect all paths - if one turns out to be a cycle, put a delay on the last edge we encountered
    std::vector<Path> paths;
    std::deque<DelayEdge> stack;

    std::vector<DelayEdge> amendedEdges = std::vector<DelayEdge>(edges);
    
    
    for (std::set<long>::iterator inputIt = inputNodes.begin(); inputIt != inputNodes.end(); inputIt++) {
        paths.push_back(Path(*inputIt));
        
        DelayEdge e = DelayEdge(-1, *inputIt);
        stack.push_back(e);
        amendedEdges.push_back(e);
    }
    
    while (stack.size() > 0) {
        DelayEdge current = stack.front();
        stack.pop_front();

        long from = current.nodeTo;
        std::vector<DelayEdge> outgoing = outputsFromNode(from);
        
        // find the paths that ended on this node, so we can extend them
        std::vector<Path>::iterator it = paths.begin();
        std::vector<Path> newPaths;
        std::vector<std::pair<long,long> > delayedEdges;
        
            do {
               // bool created_cycle = false;
                it = std::find_if(it, paths.end(), [from](const Path &p){return p.endsOnNode(from);});
                
                if (it != paths.end()) {
                    for (std::vector<DelayEdge>::iterator outIt = outgoing.begin(); outIt != outgoing.end(); outIt++) {
                        // spawn new paths that fork out here, then get rid of the original
                        // EXPONENTIAL SPACE - WOOOOOOO
                        if (!it->containsNode(outIt->nodeTo)) {
                            // not making a cycle
                            Path newPath = Path(*it);
                            newPath.extend(outIt->nodeTo);
                            newPath.totalDelay = it->totalDelay + outIt->delay;
                            newPaths.push_back(newPath);
                            stack.push_back(*outIt);
                        } else {
                            // if the delay is 0, there is a PROBLEM
                            // we need to put a delay in this path
                            if (it->totalDelay == 0) {
                                std::vector<DelayEdge>::iterator realEdge = std::find(edges.begin(), edges.begin(), *outIt);
                                realEdge->delay = 1;
                                delayedEdges.push_back(std::pair<long, long>(realEdge->nodeFrom, realEdge->nodeTo));
                                 it->totalDelay = 1;
                            }
                            
                            // either way, no need to extend this path... right?
                            // definitely no need to put this edge on the stack
                            
                        }
                    }
                    it = paths.erase(it); // get rid of this path - it either has descendants or it ended.
                }
            } while (it!=paths.end());
       
        // update paths containing delayed edges
        if (!delayedEdges.empty()) {
            for (long i=0; i<paths.size(); i++) {
                for (int j=0; j<delayedEdges.size(); j++) {
                    if (paths[i].containsEdge(delayedEdges[j].first, delayedEdges[j].second)) {
                        paths[i].totalDelay +=1;
                    }
                }
            }
        }
        
        paths.insert(paths.end(), newPaths.begin(), newPaths.end());
    }
    
    
}

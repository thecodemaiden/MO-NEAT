//
//  BaseNetwork.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "BasicNN.h"
#include <algorithm>
#include <sstream>
#include <cmath>
#include <assert.h>

// locate incoming and outgoing edges for a node number
struct has_input {
    has_input(long n) : n(n) { }
    bool operator()(Edge e) const { return e.nodeTo == n && !e.disabled; }
private:
    long n;
};

struct has_output {
    has_output(long n) : n(n) { }
    bool operator()(Edge e) const { return e.nodeFrom == n && !e.disabled;}
private:
    long n;
};

BasicNN::BasicNN(int nInputs, int nOutputs)
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
            Edge e = Edge(i, nextNode);
            nodes[i].outdegree++;
            n.indegree++;
            edges.push_back(e);
        }
        nodes.push_back(n);
    }
    
}

void BasicNN::addGeneFromParentSystem(MNIndividual *p, MNEdge *g)
{
    
    BasicNN *parent = dynamic_cast<BasicNN *>(p);
    Edge *gene = dynamic_cast<Edge *>(g);
    
    if (!parent || !gene)
        return; // can't do nothin with this
    
    if (gene->nodeFrom < 0 || gene->nodeTo < 0)
        return; // invalid
    
    std::vector<Edge>::iterator found = std::find(edges.begin(), edges.end(), *gene);
    
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
    cleanup();
}

std::vector<MNEdge *> BasicNN::connectionGenome()
{
    std::vector<MNEdge *> genome;
    for (std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); it++) {
        genome.push_back(&*it);
    }
    return genome;
}

void BasicNN::mutateConnectionWeights(double p_m)
{
    // pick a random edge (disabled is fine?) and mutate its weight
    for (size_t pos = 0; pos < edges.size(); pos++) {
        if (uniformProbability() < p_m) {
            double var = normallyDistributed();
            Edge &toMutate = edges.at(pos);
            
            toMutate.weight = var;
        }
    }
}

void BasicNN::mutateNodes(double p_m)
{
    double selector = (double)rand()/RAND_MAX;
    
    for (size_t i = 0; i<nodes.size(); i++) {
        Node &n  = nodes[i];
        if (uniformProbability() < p_m) {
            if (selector < 0.33)
                n.bias = normallyDistributed(n.bias, 1);
            else if (selector < 0.67)
                n.param1 = normallyDistributed(n.param1, 1);
            else
                n.type = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
        }
    }
}


std::vector<MNEdge *> BasicNN::createNode()
{
    // if all edges are disabled, that's BAD and we are uselesssssss
    long activeEdges = std::count_if(edges.begin(), edges.end(), [](Edge e) {return !e.disabled;});
        if (activeEdges == 0)
            return std::vector<MNEdge *>();
    
        uint pos;
        do {
            pos = arc4random_uniform((uint)edges.size());
            
        } while (edges.at(pos).disabled);
        Edge &toSplit = edges.at(pos);

        return insertNodeOnEdge(pos);
}

std::vector<MNEdge *> BasicNN::insertNodeOnEdge(long ePos)
{
    
    uint nextNode = (uint)nodes.size();
    Node newNode = Node();

    bool wasDisabled = edges[ePos].disabled;
    
    if (!wasDisabled) {
        newNode.indegree = newNode.outdegree = 1;
    }
    
    nodes.push_back(newNode);
    
    Edge e1 = Edge(edges[ePos].nodeFrom, nextNode);
    Edge e2 = Edge(nextNode, edges[ePos].nodeTo);
    
    e1.disabled = edges[ePos].disabled;
    e2.disabled = edges[ePos].disabled;
    edges[ePos].disabled = true;

    
    edges.push_back(e1);
    edges.push_back(e2);
    
    cleanup();

    
    Edge *first = new Edge(e1);
    Edge *second = new Edge(e2);
    
    std::vector<MNEdge *> toReturn;

    toReturn.push_back(first);
    toReturn.push_back(second);
    return toReturn;
}


Edge *BasicNN::createConnection()
{

    // pick two existing, unconnected nodes and connect them - can connect self to self!
    long maxConnections = nodes.size()*( nodes.size() - 1);

    long activeEdges = std::count_if(edges.begin(), edges.end(), [](Edge e) {return !e.disabled;});
    
    Edge newEdge = Edge(-1,-1);
    
    std::vector<Edge>::iterator edgePos = edges.end();
    
    if (activeEdges < maxConnections) {
        // chose a source node at random
        long maxDegree = nodes.size() - 1;
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
        std::vector<Edge> existingEdges;
        
        std::vector<Edge>::iterator it = edges.begin();
        while (it != edges.end()) {
            if (it->nodeFrom == newEdge.nodeFrom)
                existingEdges.push_back(*it);
            it++;
        }
        
        Node *sink = NULL;
        bool found = false;
        
        do {
            pos = arc4random_uniform((uint)nodes.size());
         
            if (pos == newEdge.nodeFrom)
                continue; // don't connect to ourselves!
            
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
    }
    cleanup();
    return new Edge(*edgePos);

}

double BasicNN::connectionDifference(MNEdge *e1, MNEdge *e2)
{
    
    Edge *first = dynamic_cast<Edge *>(e1);
    Edge *second = dynamic_cast<Edge *>(e2);
    
    if (!first || !second)
        return INFINITY;
    
    return fabs(first->weight - second->weight);
}

double BasicNN::nodeDifference(MNIndividual *ind)
{
    BasicNN *other = dynamic_cast<BasicNN *>(ind);
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

long BasicNN::numberOfNodes()
{
    return nodes.size();
}

long BasicNN::numberOfEdges()
{
    return edges.size();
}



std::vector<double> BasicNN::simulateSequence(const std::vector<double> &inputValues)
{
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<std::vector<double> > allOuputs;
    
        std::vector<double> currentOutputs = std::vector<double>(nNodes);
        std::vector<double> currentInputs = inputValues;
    
        for (long steps = -1; steps<1; steps++) {
            std::vector<double> lastOutputs = currentOutputs;
            currentOutputs = nodeOutputsForInputs(currentInputs, currentOutputs);
            if (lastOutputs == currentOutputs) {
                break;
            }
        }
        
        std::vector<double> finalOutputs;
        for (std::set<long>::iterator it = outputNodes.begin(); it != outputNodes.end(); it++) {
            long node = *it;
            finalOutputs.push_back(currentOutputs[node]);
        }
    
    
    return finalOutputs;
}

double BasicNN::visitNode(long i, std::set<long> &visitedNodes, std::vector<double> &lastOutputs)
{
    Node n = nodes[i];
    visitedNodes.insert(i);
    std::vector<Edge> inputEdges = inputsToNode(i);
    
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

 std::vector<double> BasicNN::nodeOutputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs)
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



std::vector<Edge> BasicNN::inputsToNode(long n)
{
    
    std::vector<Edge> found;
    
    has_input pred = has_input(n);
    
    std::vector<Edge>::iterator it = std::find_if(edges.begin(), edges.end(), pred);
    for (; it != edges.end(); it = std::find_if(++it, edges.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

std::vector<Edge> BasicNN::outputsFromNode(long n)
{
    
    std::vector<Edge> found;
    
    has_output pred = has_output(n);
    
    std::vector<Edge>::iterator it = std::find_if(edges.begin(), edges.end(), pred);
    for (; it != edges.end(); it = std::find_if(++it, edges.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

std::string BasicNN::display()
{
    // sort the edges by source node and display
    
    std::ostringstream ss;
    std::sort(edges.begin(), edges.end());
    
    std::vector<Edge>::iterator it = edges.begin();
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

std::string BasicNN::dotFormat(std::string graphName)
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
        Edge e = edges[j];
        ss << "\t" << e.nodeFrom << " -> " << e.nodeTo << " [label =\"" << e.weight << "\"";
        if (e.disabled) {
            ss << ", style=dotted";
        }
        ss << "];\n";
    }
    
    ss <<"}\n";
    
    
    return ss.str();
}

void BasicNN::cleanup()
{
    // in, out
    std::vector<std::pair<int, int> > nodeDegrees(nodes.size());
    for (long i=0; i< edges.size(); i++) {
        if (!edges[i].disabled) {
            nodeDegrees[edges[i].nodeFrom].second++;
            nodeDegrees[edges[i].nodeTo].first++;
        }
    }
    
    for (long i=0; i<nodes.size(); i++) {
        assert(nodes[i].indegree == nodeDegrees[i].first);
        assert(nodes[i].outdegree == nodeDegrees[i].second);
    }
    
}
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
    bool operator()(Edge e) const { return e.nodeFrom == n;} // disabled connections count for insertAsNode
private:
    long n;
};

std::string activationFuncName(ActivationFunc f)
{
  //  std::ostringstream name;
    std::string name;
    switch (f) {
        case SIN_FUNC:
            name = "SIN";
            break;
        case GAUSSIAN_FUNC:
            name = "GAUSSIAN";
            break;
        case TANH_FUNC:
            name = "TANH";
            break;
        case STEP_FUNC:
            name = "THRESHOLD";
            break;
        default:
            break;
    }
   // return name.str();
    return name;
}

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

void BasicNN::addGeneFromParentSystem(BasicNN parent, Edge gene)
{
    if (gene.nodeFrom < 0 || gene.nodeTo < 0)
        return; // invalid
    
    std::vector<Edge>::iterator found = std::find(edges.begin(), edges.end(), gene);
    
    if (found != edges.end()) {
        found->weight = gene.weight;
        if (found->disabled == gene.disabled)
            return; // nothing more to do, don't modify in/outdegrees
        found->disabled = gene.disabled;
    } else {
        while (gene.nodeFrom >= nodes.size() || gene.nodeTo >= nodes.size()) {
            // add nodes up to the required number, copying from the parent
            // hope that we don't have floating nodes...
            Node newNode = Node(parent.nodes[nodes.size()]);
            newNode.indegree = newNode.outdegree = 0;
            nodes.push_back(newNode);
        }
        edges.push_back(gene);
    }
    if (!gene.disabled) {
        nodes.at(gene.nodeFrom).outdegree++;
        nodes.at(gene.nodeTo).indegree++;
    }
    
    for(std::set<long>::iterator it = outputNodes.begin(); it!=outputNodes.end(); it++) {
        assert(nodes[*it].indegree > 0);
    }
    
}

static bool compareInnovationNumbers(const Edge &e1, const Edge &e2)
{
    return e1.innovationNumber < e2.innovationNumber;
}

std::vector<Edge> BasicNN::connectionGenome()
{
    std::vector<Edge> genome = std::vector<Edge>(edges);
    std::sort(genome.begin(), genome.end(), compareInnovationNumbers);
    return genome;
}

void BasicNN::mutateConnectionWeight()
{
    // pick a random edge (disabled is fine?) and mutate its weight
    double var = normallyDistributed();
    
    uint pos = arc4random_uniform((uint)edges.size());
    Edge &toMutate = edges.at(pos);
    
    toMutate.weight = var;
}

void BasicNN::mutateNode(long n)
{
    double var = normallyDistributed();
    nodes[n].bias = var;
    
    ActivationFunc t = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
    nodes[n].type = t;
}


std::vector<Edge> BasicNN::insertNode()
{
    bool insertNormal = true;
    if (insertNormal) {
        uint pos;
        pos = arc4random_uniform((uint)edges.size());
        
        Edge &toSplit = edges.at(pos);
        
        return insertNodeOnEdge(toSplit);
    } else {
        uint pos;
        pos = arc4random_uniform((uint)nodes.size());
        
        return insertNodeAsNode(pos);
    }
}

std::vector<Edge> BasicNN::insertNodeOnEdge(Edge &e)
{
    uint nextNode = (uint)nodes.size();
    Node newNode = Node();
    newNode.indegree = newNode.outdegree = 1;
    nodes.push_back(newNode);
    
    Edge e1 = Edge(e.nodeFrom, nextNode);
    Edge e2 = Edge(nextNode, e.nodeTo);
        e1.disabled = e.disabled;
        e2.disabled = e.disabled;
    e.disabled = true;
    
    edges.push_back(e1);
    edges.push_back(e2);
    
    std::vector<Edge> toReturn;
    toReturn.push_back(e1);
    toReturn.push_back(e2);
    
    return toReturn;
}

std::vector<Edge>  BasicNN::insertNodeAsNode(int n)
{
    Node newNode = Node();
    Node oldNode = nodes[n];
    
    int newNodeNum = (int)nodes.size();
    
    Edge bridge = Edge(n, newNodeNum);
    
    std::vector<Edge> toReturn;

    // gather the existing outputs from the old node
    std::vector<Edge> outputEdges = outputsFromNode(n);
    
    // add the bridge
    toReturn.push_back(bridge);
    edges.push_back(bridge);
    
    // now move the existing outputs to the new node
    for (long i=0; i<outputEdges.size(); i++) {
        std::vector<Edge>::iterator it = std::find(edges.begin(), edges.end(), outputEdges[i]);
        it->nodeFrom = newNodeNum;
        toReturn.push_back(*it);
    }
    
    newNode.outdegree = oldNode.outdegree;
    oldNode.outdegree = 1;
    
    nodes.push_back(newNode);
    
    return toReturn;
}


Edge BasicNN::createConnection()
{
    // pick two existing, unconnected nodes and connect them - can connect self to self!
    int maxConnections = (int)(nodes.size()*nodes.size());
   // int maxConnections = (int)((nodes.size() - 1)*nodes.size());

    Edge newEdge = Edge(-1,-1);
    
    if (edges.size() < maxConnections) {
        // chose a source node at random
   //     long maxDegree = nodes.size() -1;
        long maxDegree = nodes.size();
        Node source;
        uint pos;
        do {
            pos = arc4random_uniform((uint)nodes.size());
            source = nodes.at(pos);
        } while (source.outdegree >= maxDegree);
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
        
        Node sink;
        bool found = false;
        std::vector<Edge>::iterator existing = existingEdges.end();
        
        do {
            pos = arc4random_uniform((uint)nodes.size());
         
          //  if (pos == newEdge.nodeFrom)
          //      continue;
            
            sink = nodes.at(pos);
            newEdge.nodeTo = pos;
            existing = std::find(existingEdges.begin(), existingEdges.end(), newEdge);
            
            found = (existing != existingEdges.end() && !existing->disabled); // already enabled edge?
            
            
        } while (sink.indegree >= maxDegree || found);
        nodes.at(pos).indegree++;

        if (existing != existingEdges.end() && existing->disabled)
            existing->disabled = false;
        else
            edges.push_back(newEdge);
    }
    
    return newEdge;
}

void BasicNN::updateInnovationNumber(const Edge &info)
{
    //find the edge and update it
    std::vector<Edge>::iterator it = std::find(edges.begin(), edges.end(), info);
    
    if (it != edges.end()) {
        it->innovationNumber = info.innovationNumber;
    }
    
}

double BasicNN::connectionDifference(const Edge &c1, const Edge &c2)
{
    double weightDiff = fabs(c1.weight - c2.weight);
    
    return weightDiff;
}

double BasicNN::nodeDifference(BasicNN other)
{
    long length = std::min(other.nodes.size(), nodes.size());
    double d = 0;
    for (long i=0; i<length; i++) {
        if (other.nodes[i].type != nodes[i].type)
            d++;
        d+= (nodes[i].bias - other.nodes[i].bias);
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


SimReturn BasicNN::simulateTillEquilibrium(std::vector<double> inputValues, int maxSteps)
{
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<double> currentOutputs = std::vector<double>(nNodes);
    
    bool steady = false;
    long steps;
    for (steps = -1; steps<maxSteps; steps++) {
        std::vector<double> lastOutputs = currentOutputs;
        currentOutputs = nodeOuputsForInputs(inputValues, currentOutputs);
        if (lastOutputs == currentOutputs) {
            steady = true;
            break;
        }
    }
    
    std::vector<double> finalOutputs;
    for (std::set<long>::iterator it = outputNodes.begin(); it != outputNodes.end(); it++) {
        long node = *it;
        finalOutputs.push_back(currentOutputs[node]);
    }
    
    return SimReturn(finalOutputs, steady, steps);
}


static double applyActivationFunc(Node n, double inputSum)
{
    double output;
    
    switch (n.type) {
        case STEP_FUNC:
            if (inputSum > n.threshold)
                output = n.activatedVal;
            else
                output = n.deactivatedVal;
            break;
        case TANH_FUNC:
            output = tanh(inputSum);
            break;
        case SIN_FUNC:
            output = sin(inputSum);
            break;
        case GAUSSIAN_FUNC:
            output = 2*exp(-inputSum*inputSum) - 1;
            break;
        default:
            output = inputSum;
            break;
    }
    
    return output;
}


 std::vector<double> BasicNN::nodeOuputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs)
{
    std::vector<double> newOutputs;
    
    long inputNodeN = 0; // used to map values from the inputs vector to input nodes
    
    for (long i=0; i<nodes.size(); i++) {
        // collect the input values
        Node n = nodes[i];
        std::vector<Edge> inputEdges = inputsToNode(i);
        double inputSum = n.bias;
        
        for (long i=0; i<inputEdges.size(); i++) {
            int nodeN = inputEdges[i].nodeFrom;
            double rawVal = lastOutputs[nodeN];
            inputSum += inputEdges[i].weight*rawVal;
        }
        
        bool isInput = (inputNodes.find(i) != inputNodes.end());
        
        if (isInput) {
            inputSum += (inputs[inputNodeN++]);
        }
        double output = applyActivationFunc(n, inputSum);
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
        if (n.type == STEP_FUNC)
            ss << " T: " << n.threshold;
        ss << " B: " << n.bias << "\n";
    }
    
    return ss.str();
}

std::string BasicNN::dotFormat(std::string graphName)
{
    std::ostringstream ss;

    ss << "digraph " << graphName << "{\n" ;
    ss << "rankdir=LR;\n";
    ss << "center=true;\n";
    ss << "orientation=landscape;\n";
    
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
        Node n = nodes[i];
        ss << j << " [label = \"" << activationFuncName(nodes[j].type) << "\"];\n";
        if (n.bias != 0) {
            ss << "b_" << j << " -> " << j << " [label = \"" << n.bias <<"\"];\n";
        }
        ss << "\n";
    }
    
    for (long j=0; j<edges.size(); j++) {
        Edge e = edges[j];
        ss << e.nodeFrom << " -> " << e.nodeTo << " [label =\"" << e.weight << "\"";
        if (e.disabled) {
            ss << ", style=dotted";
        }
        ss << "];\n";
    }
    
    ss <<"}\n";
    
    
    return ss.str();
}
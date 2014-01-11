//
//  BaseNetwork.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 1/7/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "ExampleNetwork.h"
#include <algorithm>
#include <sstream>
#include <cmath>


ExampleNetwork::ExampleNetwork(int nInputs, int nOutputs)
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

void ExampleNetwork::addGeneFromParentSystem(ExampleNetwork parent, Edge gene)
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
            // add nodes up to the required number
            // hope that we don't have floating nodes...
            nodes.push_back(Node());
        }
        edges.push_back(gene);
    }
    if (!gene.disabled) {
        nodes.at(gene.nodeFrom).outdegree++;
        nodes.at(gene.nodeTo).indegree++;
    }
}

static bool compareInnovationNumbers(const Edge &e1, const Edge &e2)
{
    return e1.innovationNumber < e2.innovationNumber;
}

std::vector<Edge> ExampleNetwork::connectionGenome()
{
    std::vector<Edge> genome = std::vector<Edge>(edges);
    std::sort(genome.begin(), genome.end(), compareInnovationNumbers);
    return genome;
}

void ExampleNetwork::mutateConnectionWeight()
{
    // pick a random edge (disabled is fine?) and mutate its weight
    double var = ((double)rand()/RAND_MAX - 0.5)*2;
    
    uint pos = arc4random_uniform((uint)edges.size());
    Edge &toMutate = edges.at(pos);
    
    toMutate.weight = var;
}

void ExampleNetwork::mutateNode(long n)
{
    double var = ((double)rand()/RAND_MAX - 0.5)*2;
    nodes[n].bias = var;
    
    // should I just change the func?
    ActivationFunc t = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
    nodes[n].type = t;
}


std::pair<Edge, Edge> ExampleNetwork::insertNode()
{
    // pick an existing, enabled edge at random: split it and insert the node there
    uint pos;
   // do {
        pos = arc4random_uniform((uint)edges.size());
  //  } while (edges.at(pos).disabled);
    
    Edge &toSplit = edges.at(pos);
    
    return insertNodeOnEdge(toSplit);
}

std::pair<Edge, Edge> ExampleNetwork::insertNodeOnEdge(Edge &e)
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
    
    return std::pair<Edge, Edge>(e1, e2);
}

Edge ExampleNetwork::createConnection()
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

void ExampleNetwork::updateInnovationNumber(const Edge &info)
{
    //find the edge and update it
    std::vector<Edge>::iterator it = std::find(edges.begin(), edges.end(), info);
    
    if (it != edges.end()) {
        it->innovationNumber = info.innovationNumber;
    }
    
}

double ExampleNetwork::connectionDifference(const Edge &c1, const Edge &c2)
{
    return 4.0*fabs(c1.weight - c2.weight);
}


std::string ExampleNetwork::display()
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
    
    
    return ss.str();
}

long ExampleNetwork::numberOfNodes()
{
    return nodes.size();
}

long ExampleNetwork::numberOfEdges()
{
    return edges.size();
}

//long ExampleNetwork::numberOfInputNodes()
//{
//    // all nodes with indegree == 0
//    return std::count_if(nodes.begin(), nodes.end(), isInput);
//}
//
//static bool isOutput(Node n) { return n.outdegree == 0;}
//long ExampleNetwork::numberOfOutputNodes()
//{
//    return std::count_if(nodes.begin(), nodes.end(), isOutput);
//}


std::pair<std::vector<double>, bool> ExampleNetwork::simulateTillEquilibrium(std::vector<double> inputValues, int maxSteps)
{
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<double> currentOutputs = std::vector<double>(nNodes);
    
    bool steady = false;
    for (int steps = 0; steps<maxSteps; steps++) {
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
    
    return std::pair<std::vector<double>, bool>(finalOutputs, steady);
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
        case ABS_SIG_FUNC:
            output = 1.0/(1 + fabs(inputSum));
            break;
        case GAUSSIAN_FUNC:
            output = exp(-inputSum*inputSum);
            break;
        default:
            output = inputSum;
            break;
    }
    
    return output;
}


 std::vector<double> ExampleNetwork::nodeOuputsForInputs(std::vector<double> inputs, std::vector<double> lastOutputs)
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
        
        newOutputs.push_back(output);
    }
    return newOutputs;
}

struct has_input {
    has_input(long n) : n(n) { }
    bool operator()(Edge e) const { return e.nodeTo == n && !e.disabled; }
private:
    long n;
};

std::vector<Edge> ExampleNetwork::inputsToNode(long n)
{
    
    std::vector<Edge> found;
    
    has_input pred = has_input(n);
    
    std::vector<Edge>::iterator it = std::find_if(edges.begin(), edges.end(), pred);
    for (; it != edges.end(); it = std::find_if(++it, edges.end(), pred)) {
        found.push_back(*it);
    }
    return found;
}

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
    
//    for(std::set<long>::iterator it = outputNodes.begin(); it!=outputNodes.end(); it++) {
//        assert(nodes[*it].indegree > 0);
//    }
    
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
    double selector = (double)rand()/RAND_MAX;
    
    if (selector < 0.33)
        nodes[n].bias = normallyDistributed();
    else if (selector < 0.67)
        nodes[n].param1 = normallyDistributed(0,2);
    else
        nodes[n].type = (ActivationFunc)arc4random_uniform(FUNC_SENTINEL);
}


std::vector<Edge> BasicNN::insertNode()
{
        uint pos;
        do {
            pos = arc4random_uniform((uint)edges.size());
            
        } while (edges.at(pos).disabled);
        Edge &toSplit = edges.at(pos);

        return insertNodeOnEdge(toSplit);
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

//std::vector<Edge>  BasicNN::insertNodeAsNode(int n)
//{
//    Node newNode = Node();
//    Node oldNode = nodes[n];
//    
//    int newNodeNum = (int)nodes.size();
//    
//    Edge bridge = Edge(n, newNodeNum);
//    
//    std::vector<Edge> toReturn;
//
//    // gather the existing outputs from the old node
//    std::vector<Edge> outputEdges = outputsFromNode(n);
//    
//    // add the bridge
//    toReturn.push_back(bridge);
//    edges.push_back(bridge);
//    
//    // now move the existing outputs to the new node
//    for (long i=0; i<outputEdges.size(); i++) {
//        std::vector<Edge>::iterator it = std::find(edges.begin(), edges.end(), outputEdges[i]);
//        it->nodeFrom = newNodeNum;
//        toReturn.push_back(*it);
//    }
//    
//    newNode.outdegree = oldNode.outdegree;
//    oldNode.outdegree = 1;
//    
//    nodes.push_back(newNode);
//    
//    return toReturn;
//}


Edge BasicNN::createConnection()
{
    // pick two existing, unconnected nodes and connect them - can connect self to self!
    long maxConnections = (nodes.size()*nodes.size());

    long activeEdges = std::count_if(edges.begin(), edges.end(), [](Edge e) {return !e.disabled;});
    
    Edge newEdge = Edge(-1,-1);
    
    if (activeEdges < maxConnections) {
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



std::vector<std::vector<double> > BasicNN::simulateSequence(std::vector<std::vector<double> > inputValues, int delay)
{
    // ensure that we have the right number of input values
    long nNodes = nodes.size();
    
    std::vector<std::vector<double> > allOuputs;
    
    for (long i=0; i<inputValues.size(); i++) {
        std::vector<double> currentOutputs = std::vector<double>(nNodes);
        std::vector<double> currentInputs = inputValues[i];
        for (long steps = -1; steps<delay; steps++) {
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
        allOuputs.push_back(finalOutputs);
    }
    
    return allOuputs;
}


static double applyActivationFunc(Node n, double inputSum)
{
    double output;
    
    switch (n.type) {
        case STEP_FUNC:
            if (inputSum > n.param1)
                output = n.activatedVal;
            else
                output = n.deactivatedVal;
            break;
        case TANH_FUNC:
            output = tanh(inputSum/exp(n.param1));
            break;
        case SIN_FUNC:
            output = sin(inputSum*exp(n.param1));
            break;
        case GAUSSIAN_FUNC:
        {
            double x = inputSum - n.param1;
            output = 2*exp(-x*x) - 1;
            break;
        }
        default:
            output = inputSum;
            break;
    }
    
    return output;
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
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

ExampleNetwork::ExampleNetwork()
{
    nodes.push_back(Node());
    nodes.push_back(Node());
    
    Edge e = Edge(0, 1);
    nodes[0].outdegree += 1;
    nodes[1].indegree += 1;
    edges.push_back(e);
}

void ExampleNetwork::addGeneFromParentSystem(ExampleNetwork parent, Edge gene)
{
    
}

std::vector<Edge> ExampleNetwork::attachmentGenome()
{
    std::vector<Edge> genome;
    
    return genome;
}

void ExampleNetwork::mutateAttachmentWeight()
{
    
}

std::pair<Edge, Edge> ExampleNetwork::insertNode()
{
    // pick an existing edge at random: split it and insert the node there
    uint pos = arc4random_uniform((uint)edges.size());
    Edge toSplit = edges.at(pos);
    
    uint nextNode = (uint)nodes.size();
    nodes.push_back(Node());
    Edge e1 = Edge(toSplit.nodeFrom, nextNode);
    Edge e2 = Edge(nextNode, toSplit.nodeTo);
    
    toSplit.disabled = true;
    
    edges.push_back(e1);
    edges.push_back(e2);
    
    return std::pair<Edge, Edge>(e1, e2);
}

Edge ExampleNetwork::createConnection()
{
    // pick two existing, unconnected nodes and connect them - can connect self to self!
    int maxConnections = (int)(nodes.size()*nodes.size());
    
    Edge newEdge = Edge(-1,-1);
    
    if (edges.size() < maxConnections) {
        // chose a source node at random
        int maxDegree = maxConnections/2;
        Node source;
        uint pos;
        do {
            pos = arc4random_uniform((uint)nodes.size());
            source = nodes.at(pos);
        } while (source.outdegree < maxDegree);
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
        std::vector<Edge>::iterator found;
        do {
            pos = arc4random_uniform((uint)nodes.size());
            sink = nodes.at(pos);
            newEdge.nodeTo = pos;
            found = std::find(existingEdges.begin(), existingEdges.end(), newEdge);
        } while (sink.indegree < maxDegree && found != edges.end());
        nodes.at(pos).indegree++;
    }
    
    edges.push_back(newEdge);
    
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

double ExampleNetwork::attachmentDifference(const Edge &c1, const Edge &c2)
{
    return 1.0;
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

//
//  BayesNet.cpp
//  NEATBox
//
//  Created by Adeola Bannis on 2/23/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#include "BayesNet.h"


BayesNet::BayesNet(std::vector<std::string>nodeNames, std::map<std::vector<int>, double> jdt)
:lastEdgeCount(0)
{
    int i=0;
    for (std::vector<std::string>::iterator it = nodeNames.begin(); it != nodeNames.end(); it++) {
        BayesNode bn = BayesNode(*it);
        if (it != nodeNames.begin()) {
            BayesEdge e = BayesEdge(i-1,i);
            edges.push_back(e);
        }
        nodes.push_back(bn);
        i++;
    }
    // create an arbitrary structure - a chain from one to the other
}

std::vector<MNEdge *> BayesNet::connectionGenome()
{
    std::vector<MNEdge *> toReturn;
    
    return toReturn;
}

void BayesNet::mutateConnectionWeights(double p_m)
{
    
}

void BayesNet::mutateNodes(double p_m)
{
    
}

MNEdge * BayesNet::createConnection()
{
    return NULL;
}

std::vector<MNEdge *> BayesNet::createNode()
{
    std::vector<MNEdge *> toReturn;
    
    return toReturn;
}

long BayesNet::numberOfNodes()
{
    return nodes.size();
}

long BayesNet::numberOfEdges()
{
    return lastEdgeCount;
}

double BayesNet::nodeDifference(MNIndividual *other)
{
    return 0.0;
}

void BayesNet::addGeneFromParentSystem(MNIndividual *parent, MNEdge *gene)
{
    
}

double BayesNet::connectionDifference(MNEdge *e1, MNEdge *e2)
{
    return 0.0;
}
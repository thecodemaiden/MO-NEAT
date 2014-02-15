//
//  LayeredNN.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/13/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef __NEATBox__LayeredNN__
#define __NEATBox__LayeredNN__

#include <iostream>
#include "BasicNN.h"

class LayeredNN :public BasicNN {
    LayeredNN(int nInputs, int nOutputs):BasicNN(nInputs, nOutputs){}
    
    virtual LayeredNN *clone() const {return new LayeredNN(*this);}
protected:
    virtual bool isRecurrent() {return true;}
    virtual std::vector<MNEdge *> insertNodeOnEdge(Edge &e);
    virtual  Edge *createConnection();
};

#endif /* defined(__NEATBox__LayeredNN__) */

//
//  MNEdge.h
//  NEATBox
//
//  Created by Adeola Bannis on 2/3/14.
//  Copyright (c) 2014 Adeola Bannis. All rights reserved.
//

#ifndef NEATBox_MNEdge_h
#define NEATBox_MNEdge_h

class MNEdge {
    
public:
    virtual MNEdge * clone() const = 0;
    virtual ~MNEdge(){};
    
    virtual bool operator==(const MNEdge& other)const = 0;

};

#endif

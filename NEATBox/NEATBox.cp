/*
 *  NEATBox.cp
 *  NEATBox
 *
 *  Created by Adeola Bannis on 1/7/14.
 *  Copyright (c) 2014 Adeola Bannis. All rights reserved.
 *
 */

#include <iostream>
#include "NEATBox.h"
#include "NEATBoxPriv.h"

void NEATBox::HelloWorld(const char * s)
{
	 NEATBoxPriv *theObj = new NEATBoxPriv;
	 theObj->HelloWorldPriv(s);
	 delete theObj;
};

void NEATBoxPriv::HelloWorldPriv(const char * s) 
{
	std::cout << s << std::endl;
};


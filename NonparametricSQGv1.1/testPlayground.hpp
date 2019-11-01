//
//  testPlayground.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef testPlayground_hpp
#define testPlayground_hpp

#include <stdio.h>
#include "NonparametricSQG_dataType.hpp"
#include "nonparametricSQG.hpp"
#include "ioNonparametricDB.hpp"
#include "ioNonparametricModel.hpp"
#include "ioStochastic.hpp"


// Target: ioStochastic
// test functions
void testStochastic();
void testStochasticOneStagePiecewiseLP();
void testStochasticUnorderedMap(const oneStageParameters& parameters);
#endif /* testPlayground_hpp */

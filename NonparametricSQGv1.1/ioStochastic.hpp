//
//  ioStochastic.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef ioStochastic_hpp
#define ioStochastic_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include "NonparametricSQG_dataType.hpp"

// input funtions for random vectors
secondStageRHS readStochastic(const std::string& stochasticPath);
oneStagePiecewiseLinearObjective readStochasticOneStagePiecewiseLP(const std::string& stochasticPath);

// functions for merging nonparametericDB and randomVector
secondStageRHSDB mergeDB_randomVector(const std::vector<std::vector<dataPoint>>& be_DB, const std::vector<std::vector<dataPoint>>& bi_DB, const secondStageRHS& bRHS);

// functions for mergeing nonparametricDB and random objetive function
oneStagePiecewiseLinearObjectiveDB mergeDB_randomObjective(const std::vector<std::vector<dataPoint>>& DB, oneStagePiecewiseLinearObjective parameters_objective, long x_size);

// output functions
void printStochastic(const secondStageRHS& bRHS);
void printStochasticOneStagePiecewiseLP(oneStagePiecewiseLinearObjective parameters_objective, oneStageParameters parameters);


#endif /* ioStochastic_hpp */

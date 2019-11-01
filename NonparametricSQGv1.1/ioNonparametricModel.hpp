//
//  ioNonparametricModel.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef ioNonparametricModel_hpp
#define ioNonparametricModel_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include "NonparametricSQG_dataType.hpp"

// input functions for two stage linear programming
twoStageParameters readTwoStageParameters(const std::string& parameterPath);
// input functions for one stage piecewise linear programming
oneStageParameters readOneStageParameters(const std::string& parameterPath);
// print parameters
void printTwoStageParameters(const twoStageParameters& parameters);
void printOneStageParameters(oneStageParameters parameters);
// test
void inputModelTest();
void inputOneStageModelTest();
#endif /* ioNonparametricModel_hpp */

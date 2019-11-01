//
//  testPlayground.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "testPlayground.hpp"

// Target: ioStochastic
// test functions
void testStochastic() {
    std::string filePath = "/Users/sonny/Documents/predictorDB/stochastic01.txt";
    secondStageRHS bRHS = readStochastic(filePath);
    printStochastic(bRHS);
} // end testStochastic
void testStochasticOneStagePiecewiseLP() {
    std::string filePath1 = "/Users/sonny/Documents/DebuggingPlayground/oneStagePiecewiseLP/model.txt";
    std::string filePath2 = "/Users/sonny/Documents/DebuggingPlayground/oneStagePiecewiseLP/sto.txt";
    oneStageParameters parameters = readOneStageParameters(filePath1);
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(filePath2);
    printStochasticOneStagePiecewiseLP(parameters_objective, parameters);
    // test on constant reference
    std::cout << "Test on constant reference" << std::endl;
    testStochasticUnorderedMap(parameters);
}
void testStochasticUnorderedMap(const oneStageParameters& parameters) {
    // size of A
    long A_rowsize = parameters.A.size();
    //
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        for (auto a: parameters.A[row_index]) {
            std::cout << "row: " << row_index << " col: " << a.first << " value:" << a.second << std::endl;
        }
    }
}

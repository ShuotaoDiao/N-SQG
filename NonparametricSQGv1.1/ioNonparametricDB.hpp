//
//  ioNonparametricDB.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef ioNonparametricDB_hpp
#define ioNonparametricDB_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "NonparametricSQG_dataType.hpp"

std::vector<std::vector<dataPoint>> readNonparametricDB(std::string readPath); // read database from a text file
void printNonparametricDB(const std::vector<std::vector<dataPoint>>& dataPointDB);
// print dataPoint
void printDataPoint(const dataPoint& dataPoint01);
void inputDBTest(); // test on input functions of nonparametric DB

void dataPointTest01(); // test on robustness of dataPoint
#endif /* ioNonparametricDB_hpp */

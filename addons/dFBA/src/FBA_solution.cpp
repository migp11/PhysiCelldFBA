/*
 * Reaction.cpp
 *
 *  Created on: jun. 2022
 *      Author: mponce
 */

#include "FBA_solution.h"

#include <iostream>
#include <vector>
#include <map>


FBA_solution::FBA_solution()
{

}

FBA_solution::FBA_solution(float objective_value, std::string status, std::map<std::string,float> fluxes)
{
    this->objective_value = objective_value;
    this->status = status;
    this->fluxes = fluxes;
}

FBA_solution::~FBA_solution()
{

}

/*
 * Metabolite.cpp
 *
 *  Created on: 13 jun. 2019
 *      Author: mponce
 */


#include "PhysiFBA_metabolite.h"

#include <iostream>
#include <vector>
#include <map>



PhysiFBA_metabolite::PhysiFBA_metabolite(std::string id) {
   this->id = id;
}

PhysiFBA_metabolite::~PhysiFBA_metabolite() {

}

const std::string& PhysiFBA_metabolite::getId() const {
  return this->id;
}

void PhysiFBA_metabolite::setName(std::string value) {
  this->name = value;
}

const std::string& PhysiFBA_metabolite::getName() const {
  return this->name;
}



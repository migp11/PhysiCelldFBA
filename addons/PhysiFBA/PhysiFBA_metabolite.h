/*
 * Metabolite.h
 *
 *  Created on: 13 jun. 2019
 *      Author: mponce
 */

#ifndef SRC_METABOLITE_H_
#define SRC_METABOLITE_H_

#include <iostream>
#include <vector>
#include <map>

class PhysiFBA_metabolite
{
	private:
		std::string id;
		std::string name;

	public:
		PhysiFBA_metabolite(std::string id);
		~PhysiFBA_metabolite();

		const std::string& getId() const;

		void setName(std::string value);
		const std::string& getName() const;

};



#endif /* SRC_METABOLITE_H_ */

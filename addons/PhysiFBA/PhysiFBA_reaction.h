/*
 * Reaction.h
 *
 *  Created on: 13 jun. 2019
 *      Author: mponce
 */

#ifndef SRC_REACTION_H_
#define SRC_REACTION_H_

#include <iostream>
#include <vector>
#include <map>

#include "PhysiFBA_metabolite.h"

class PhysiFBA_reaction
{
	private:
		std::string id;
		std::string name;

		double lowerBound;
		double upperBound;
		double objectiveCoefficient;
		double fluxValue;

		std::map<const PhysiFBA_metabolite*, double> metabolites;
		std::map<std::string, const PhysiFBA_metabolite*> idMetaboliteMap;

	public:
		PhysiFBA_reaction(std::string id);
		~PhysiFBA_reaction();

		const std::string& getId() const;

		void setName(std::string name);
		const std::string& getName() const;

		int getNumberOfMetabolites();

		void setLowerBound(double lowerBound);
		double getLowerBound();

		void setUpperBound(double upperBound);
		double getUpperBound();

		void setObjectiveCoefficient(double ojectiveCoefficient);
		double getObjectiveCoefficient();

		void setFluxValue(double flux_value);
		double getFluxValue();

		bool reversible();
		bool hasMetabolite(std::string mId);
		void addMetabolite(const PhysiFBA_metabolite* met, double stoich);

		std::vector<std::string> getReactants();
		std::vector<std::string> getProducts();
		double getStoichCoefficient(std::string mId);

		std::string getReactionString();

		const std::map<const PhysiFBA_metabolite*, double>& getMetabolites() const;
};



#endif /* SRC_REACTION_H_ */

/* This file is part of MODEL, the Mechanics Of Defect Evolution Library.
 *
 * Copyright (C) 2015 by Giacomo Po <gpo@ucla.edu>.
 *
 * model is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_aLoopIndividualSpecification_H_
#define model_aLoopIndividualSpecification_H_

#include <string>
#include <vector>
#include <Eigen/Dense>

#include <MicrostructureSpecificationBase.h>

namespace model
{

struct aLoopIndividualSpecification : public MicrostructureSpecificationBase
{
    std::vector<int> planeIDs;
    std::vector<double> loopRadii;
    Eigen::Matrix<double,Eigen::Dynamic,3> loopCenters;
    std::vector<int> loopSides;
    std::vector<int> isVacancyLoop;
    
    aLoopIndividualSpecification();
    aLoopIndividualSpecification(const std::string& fileName);
};

}
#endif

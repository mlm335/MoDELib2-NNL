/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoise_cpp
#define model_GlidePlaneNoise_cpp

#include <filesystem>
#include <GlidePlaneNoise.h>

namespace model
{

    GlidePlaneNoise::GlidePlaneNoise(const PolycrystallineMaterialBase& mat)
    {
        
        const auto noiseFiles(TextFileParser(mat.materialFile).readStringVector("glidePlaneNoise"));
        for(const auto& pair : noiseFiles)
        {
            const std::string noiseFileName(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser::removeSpaces(pair.first));
            TextFileParser parser(noiseFileName);
            const std::string type(TextFileParser::removeSpaces(parser.readString("type",true)));
            const std::string  tag(TextFileParser::removeSpaces(parser.readString("tag",true)));
            const int seed(parser.readScalar<int>("seed",true));
            const Eigen::Matrix<int,1,3> gridSize(parser.readMatrix<int,1,3>("gridSize",true));
            const Eigen::Matrix<double,1,3> gridSpacing(parser.readMatrix<double,1,3>("gridSpacing_SI",true)/mat.b_SI);

            if(type=="AnalyticalSolidSolutionNoise")
            {
                const double a(parser.readScalar<double>("spreadLstress_SI",true)/mat.b_SI);      // spreading length for stresses [AA]
                const double a_Cai(parser.readScalar<double>("a_cai_SI",true)/mat.b_SI);
                const double MSSS(parser.readScalar<double>("MSSS_SI",true)/std::pow(mat.mu_SI,2));
                const auto success(solidSolutionNoise().emplace(tag,new AnalyticalSolidSolutionNoise(tag,seed,gridSize,gridSpacing,a,a_Cai,MSSS)));
                if(!success.second)
                {
                    throw std::runtime_error("Could not insert noise "+tag);
                }
            }
            if(type=="MDSolidSolutionNoise")
            {
                //                solidSolutionNoise().emplace_back(new AnalyticalSolidSolutionNoise(mat,
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<int,1,2>("gridSize",true),
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<double,1,2>("gridSpacing_SI",true)/mat.b_SI));
            }
            if(type=="MDStackingFaultNoise")
            {
                //                stackingFaultNoise().emplace_back(new AnalyticalSolidSolutionNoise(mat,
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<int,1,2>("gridSize",true),
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<double,1,2>("gridSpacing_SI",true)/mat.b_SI));
            }
            if(type=="MDShortRangeOrderNoise")
            {
                //                stackingFaultNoise().emplace_back(new AnalyticalSolidSolutionNoise(mat,
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<int,1,2>("gridSize",true),
                //                                                                                  TextFileParser(mat.materialFile).readMatrix<double,1,2>("gridSpacing_SI",true)/mat.b_SI));
            }
        }
        
        for(auto& pair : solidSolutionNoise())
        {
            pair.second->computeRealNoise();
            pair.second->computeRealNoiseStatistics(mat);
        }
        
        for(auto& pair : stackingFaultNoise())
        {
            pair.second->computeRealNoise();
            pair.second->computeRealNoiseStatistics(mat);
        }
        
        //        if(solidSolution)
        //        {
        //            Eigen::VectorXd rowsAvr0(Eigen::VectorXd::Zero(solidSolution->gridSize(0)));
        //            Eigen::VectorXd colsAvr0(Eigen::VectorXd::Zero(solidSolution->gridSize(1)));
        //            Eigen::VectorXd rowsAvr1(Eigen::VectorXd::Zero(solidSolution->gridSize(0)));
        //            Eigen::VectorXd colsAvr1(Eigen::VectorXd::Zero(solidSolution->gridSize(1)));
        //
        //            for(size_t k=0; k<solidSolution->size(); ++k)
        //            {
        //                const GridSizeType rowCol(rowAndColIndices(k));
        //                rowsAvr0(rowCol(0)) += this->solidSolution->operator[](k)(0);
        //                colsAvr0(rowCol(1)) += this->solidSolution->operator[](k)(0);
        //                rowsAvr1(rowCol(0)) += this->solidSolution->operator[](k)(1);
        //                colsAvr1(rowCol(1)) += this->solidSolution->operator[](k)(1);
        //            }
        //
        //            const auto rowsNorm0(rowsAvr0.norm());
        //            const auto colsNorm0(colsAvr0.norm());
        //            const auto rowsNorm1(rowsAvr1.norm());
        //            const auto colsNorm1(colsAvr1.norm());
        //
        //            std::cout<<"rowsNorm0= "<<rowsNorm0<< " ,rowsAvr0.size= " <<rowsAvr0.size() << std::endl;
        //            std::cout<<"colsNorm0= "<<colsNorm0<< " ,colsAvr0.size= " <<colsAvr0.size() << std::endl;
        //            std::cout<<"rowsNorm1= "<<rowsNorm1<< " ,rowsAvr1.size= " <<rowsAvr1.size() << std::endl;
        //            std::cout<<"colsNorm1= "<<colsNorm1<< " ,colsAvr1.size= " <<colsAvr1.size() << std::endl;
        //        }
    }

    //    typename GlidePlaneNoise::GridSizeType GlidePlaneNoise::rowAndColIndices(const int& storageIdx) const
    //    {
    //        return GridSizeType(storageIdx/this->gridSize(1), storageIdx%this->gridSize(1));
    //    }

    //    int GlidePlaneNoise::storageIndex(const int& i,const int& j) const
    //    {/*!\param[in] localPos the  position vector on the grid
    //      * \returns The grid index periodically wrapped within the gridSize bounds
    //      */
    //        return gridSize(1)*i+j;
    //    }

    std::tuple<double,double,double> GlidePlaneNoise::gridInterp(const Eigen::Matrix<double,2,1>& localPos) const
    {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
        
    //    std::cout<<"gridInterp 1"<<std::endl;
        
        double effsolNoiseXZ(0.0);
        double effsolNoiseYZ(0.0);
        for(const auto& noise : solidSolutionNoise())
        {
            const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(localPos));
            for(size_t p=0;p<idxAndWeights.first.size();++p)
            {
                const int storageID(noise.second->storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
                effsolNoiseXZ+=noise.second->operator[](storageID)(0)*idxAndWeights.second[p];
                effsolNoiseYZ+=noise.second->operator[](storageID)(1)*idxAndWeights.second[p];
            }
        }
        
        double effsfNoise(0.0);
        for(const auto& noise : stackingFaultNoise())
        {
            const auto idxAndWeights(noise.second->posToPeriodicCornerIdxAndWeights(localPos));
            for(size_t p=0;p<idxAndWeights.first.size();++p)
            {
                const int storageID(noise.second->storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
                effsfNoise+=noise.second->operator[](storageID)*idxAndWeights.second[p];
            }
        }
        
        //        for(size_t p=0;p<idxAndWeights.first.size();++p)
        //        {
        //            const int storageID(storageIndex(idxAndWeights.first[p](0),idxAndWeights.first[p](1)));
        //            if(solidSolution)
        //            {
        //                effsolNoiseXZ+=solidSolution->operator[](storageID)(0)*idxAndWeights.second[p];
        //                effsolNoiseYZ+=solidSolution->operator[](storageID)(1)*idxAndWeights.second[p];
        //            }
        //            if(stackingFault)
        //            {
        //                effsfNoise+=stackingFault->operator[](storageID)*idxAndWeights.second[p];
        //            }
        //        }
    //    std::cout<<"gridInterp 10"<<std::endl;
        return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise);
    }

    std::tuple<double,double,double> GlidePlaneNoise::gridVal(const Eigen::Array<int,2,1>& idx) const
    {   // Added by Hyunsoo (hyunsol@g.clemson.edu)
        
        double effsolNoiseXZ(0.0);
        double effsolNoiseYZ(0.0);
        for(const auto& noise : solidSolutionNoise())
        {
            const Eigen::Array<int,2,1> pidx(noise.second->idxToPeriodicIdx(idx));
            const int storageID(noise.second->storageIndex(pidx(0),pidx(1)));
            effsolNoiseXZ+=noise.second->operator[](storageID)(0);
            effsolNoiseYZ+=noise.second->operator[](storageID)(1);
            
        }
        
        double effsfNoise(0.0);
        for(const auto& noise : stackingFaultNoise())
        {
            const Eigen::Array<int,2,1> pidx(noise.second->idxToPeriodicIdx(idx));
            const int storageID(noise.second->storageIndex(pidx(0),pidx(1)));
            effsfNoise+=noise.second->operator[](storageID);
        }
        
        return std::make_tuple(effsolNoiseXZ,effsolNoiseYZ,effsfNoise);
    }

    const typename GlidePlaneNoise::SolidSolutionNoiseContainer& GlidePlaneNoise::solidSolutionNoise() const
    {
        return *this;
    }

    typename GlidePlaneNoise::SolidSolutionNoiseContainer& GlidePlaneNoise::solidSolutionNoise()
    {
        return *this;
    }

    const typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::stackingFaultNoise() const
    {
        return *this;
    }

    typename GlidePlaneNoise::StackingFaultNoiseContainer& GlidePlaneNoise::stackingFaultNoise()
    {
        return *this;
    }

}
#endif


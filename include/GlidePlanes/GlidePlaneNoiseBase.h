/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoiseBase_H
#define model_GlidePlaneNoiseBase_H

#include <cmath>
#include <random>
#include <Eigen/Dense>

#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
#include <fftw3.h>
#endif

#include <NoiseTraits.h>
#include <UniformPeriodicGrid.h>

namespace model
{


    template <int N>
    struct GlidePlaneNoiseBase : public UniformPeriodicGrid<2>
    /*                        */,public NoiseTraits<N>::NoiseContainerType
    {
        typedef typename NoiseTraitsBase::REAL_SCALAR REAL_SCALAR;
        typedef typename NoiseTraitsBase::COMPLEX COMPLEX;
        typedef typename NoiseTraitsBase::GridSizeType GridSizeType;
        typedef typename NoiseTraitsBase::GridSpacingType GridSpacingType;
        typedef typename NoiseTraits<N>::NoiseType NoiseType;
        typedef typename NoiseTraits<N>::NoiseContainerType NoiseContainerType;

        const std::string tag;
        const int seed;
        const GridSizeType gridSize;
        const GridSpacingType gridSpacing;
//        const size_t NK;
        // lattice basis1
        // lattice basis2
        
        GlidePlaneNoiseBase(const std::string& tag_in,
                            const int& seed_in,
                            const NoiseTraitsBase::GridSizeType& gridSize_in,
                            const NoiseTraitsBase::GridSpacingType& gridSpacing_SI_in);
        virtual ~GlidePlaneNoiseBase(){};

        void computeRealNoise();
        void computeRealNoiseStatistics() const;
        
        GridSizeType rowAndColIndices(const int& storageIndex) const;
        
        int storageIndex(const int& i,const int& j) const;
        
        const NoiseContainerType& noiseVector() const;
        NoiseContainerType& noiseVector();
        
        virtual std::array<COMPLEX,N> kCorrelations(const Eigen::Matrix<double,3,1>& k,const Eigen::Matrix<int,3,1>& kID) const = 0;


    };


}
#endif


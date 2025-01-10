/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_GlidePlaneNoiseBase_cpp
#define model_GlidePlaneNoiseBase_cpp

#include <chrono>
#include <GlidePlaneNoiseBase.h>
#include <TerminalColors.h>

namespace model
{

    template <int N>
    GlidePlaneNoiseBase<N>::GlidePlaneNoiseBase(const std::string& tag_in,
                                                const int& seed_in,
                                                const NoiseTraitsBase::GridSizeType& gridSize_in,
                                                const NoiseTraitsBase::GridSpacingType& gridSpacing_in):
    /*init*/ UniformPeriodicGrid<2>(gridSize_in,gridSpacing_in)
    /*init*/,tag(tag_in)
    /*init*/,seed(seed_in)
    /*init*/,gridSize(gridSize_in)
    /*init*/,gridSpacing(gridSpacing_in)
//    /*init*/,NK(gridSize(0)*gridSize(1)*(gridSize(2)/2+1))
    {
//        std::cout<<"gridSize="<<gridSize<<std::endl;
//        std::cout<<"gridSize(2)="<<gridSize(2)<<std::endl;
        
//        THIS IS WRONG BECAUSE gridSize(2) must be 64!!
        
        std::cout<<greenBoldColor<<"Creating GlidePlaneNoiseBase "+tag<<defaultColor<<std::endl;
        
    }

    template <int N>
    const typename GlidePlaneNoiseBase<N>::NoiseContainerType& GlidePlaneNoiseBase<N>::noiseVector() const
    {
        return *this;
    }

    template <int N>
    typename GlidePlaneNoiseBase<N>::NoiseContainerType& GlidePlaneNoiseBase<N>::noiseVector()
    {
        return *this;
    }

    template <int N>
    int GlidePlaneNoiseBase<N>::storageIndex(const int& i,const int& j) const
    {/*!\param[in] localPos the  position vector on the grid
      * \returns The grid index periodically wrapped within the gridSize bounds
      */
        return gridSize(1)*i+j;
    }

    template <int N>
    void GlidePlaneNoiseBase<N>::computeRealNoise()
    {
        std::cout<<tag<<" computing noise "<<std::flush;
        const auto t0= std::chrono::system_clock::now();

        const int NX(gridSize(0));
        const int NY(gridSize(1));
        const int NZ(64);
        const int NK(NX*NY*(NZ/2+1));
        const int NR(NX*NY*NZ);
        const double LX(NX*gridSpacing(0));
        const double LY(NY*gridSpacing(1));
        const double LZ(NZ*0.5*(gridSpacing(0)+gridSpacing(1)));


#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_
        

        
    // Apply Gaussian noise to kCorrelations
    std::default_random_engine generator(seed);
    std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);
    
//        std::cout<<"here 1"<<std::endl;
        
    std::vector<COMPLEX*> kNoisyCorrelations;
    std::vector<REAL_SCALAR*> rNoisyCorrelations;
//    std::vector<fftw_plan> fftPlans;
    for(int n=0;n<N;++n)
    {
        kNoisyCorrelations.push_back((COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK));
        rNoisyCorrelations.push_back((REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR));
    }
        
//        std::cout<<"here 2"<<std::endl;
    
//    for(int n=0;n<N;++n)
//    {
//        fftPlans.push_back(fftw_plan_dft_c2r_3d(gridSize(0), gridSize(1), gridSize(2), reinterpret_cast<fftw_complex*>(kNoisyCorrelations[n]), rNoisyCorrelations[n], FFTW_ESTIMATE));
//    }
        
//        std::cout<<"here 3"<<std::endl;
    
//        std::cout<<", randomizing "<<std::flush;
    for(int i=0; i<NX; i++)
    {
        for(int j=0; j<NY; j++)
        {
            for(int k=0; k<(NZ/2+1); k++)
            {
                const int ind = NY*(NZ/2+1)*i + j*(NZ/2+1) + k;
                
                REAL_SCALAR kx = 2.*M_PI/LX*REAL_SCALAR(i);
                if(i>NX/2)
                {
                    kx = 2.*M_PI/LX*REAL_SCALAR(i-NX);
                }
                
                REAL_SCALAR ky = 2*M_PI/LY*REAL_SCALAR(j);
                if(j>NY/2)
                {
                    ky = 2.*M_PI/LY*REAL_SCALAR(j-NY);
                }
                
                REAL_SCALAR kz = 2.*M_PI/LZ*REAL_SCALAR(k);
                
                const Eigen::Matrix<double,3,1> kv((Eigen::Matrix<double,3,1>()<<kx,ky,kz).finished());
                const Eigen::Matrix<int,3,1> kvID((Eigen::Matrix<int,3,1>()<<i,j,k).finished());
                
                // random numbers
                REAL_SCALAR Nk_yz = distribution(generator);
                REAL_SCALAR Mk_yz = distribution(generator);
                REAL_SCALAR Nk_xz, Mk_xz;
                if(kx*ky>=0)
                {
                    Nk_xz = Nk_yz;
                    Mk_xz = Mk_yz;
                }
                else
                {
                    Nk_xz = -Nk_yz;
                    Mk_xz = -Mk_yz;
                }
                
                const double kCorrFactor((k==0 || k==NZ/2)? 1.0 : 2.0); // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
                const auto kCorr(kCorrelations(kv,kvID));
                for(int n=0;n<N;++n)
                {
//                    std::cout<<kCorr[n]/kCorrFactor<<std::endl;
//                    kNoisyCorrelations[n][ind]=sqrt(kCorr[n]/kCorrFactor);
                    kNoisyCorrelations[n][ind]=sqrt(kCorr[n]/kCorrFactor)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    std::cout<<kNoisyCorrelations[n][ind]<<std::endl;
                }
//                if(k==0 || k==NZ/2) // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
//                {
//                    for(int n=0;n<N;++n)
//                    {
//                        kNoisyCorrelations[n][ind]=sqrt(kCorr[n])*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    }
//                }
////                else if(k==NZ/2)
////                {
////                    for(int n=0;n<N;++n)
////                    {
////                        kNoisyCorrelations[n][ind]=sqrt(kCorr[n])*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
////                    }
////                }
//                else
//                {
//                    for(int n=0;n<N;++n)
//                    {
//                        kNoisyCorrelations[n][ind]=sqrt(kCorr[n]/2.0)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    }
//                }
            }
        }
    }
    
//        std::cout<<"here 4"<<std::endl;
        
    // Compute Real noise
        std::cout<<", ifft "<<std::flush;
    for(int n=0;n<N;++n)
    {
        kNoisyCorrelations[n][0]=0;
//        std::cout<<rNoisyCorrelations[n][120000]<<std::endl;;
//        std::cout<<kNoisyCorrelations[n][120000]<<std::endl;;
//        std::cout<<"here 41"<<std::endl;
        fftw_plan nPlan = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(kNoisyCorrelations[n]), rNoisyCorrelations[n], FFTW_ESTIMATE);
        //        fftPlans.push_back(fftw_plan_dft_c2r_3d(gridSize(0), gridSize(1), gridSize(2), reinterpret_cast<fftw_complex*>(kNoisyCorrelations[n]), rNoisyCorrelations[n], FFTW_ESTIMATE));

        
        fftw_execute(nPlan);
    }
        
//        std::cout<<"here 5"<<std::endl;
    
    // Store Real noise
    this->reserve(NX*NY);
    const int k=0;
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            const int ind = NY*NZ*i + j*NZ + k;
            
            Eigen::Matrix<REAL_SCALAR,1,N> temp;
            for(int n=0;n<N;++n)
            {
                temp(n)=rNoisyCorrelations[n][ind];
            }
            this->push_back(NoiseTraits<N>::fromMatrix(temp));
        }
    }
//        std::cout<<"here 6"<<std::endl;
#else
        
//        std::cout<<"here 7"<<std::endl;
    this->resize(NX*NY,NoiseTraits<N>::Zero());
//        std::cout<<"here 8"<<std::endl;
#endif
        std::cout<<greenColor<<" ["<<(std::chrono::duration<double>(std::chrono::system_clock::now()-t0)).count()<<" sec]"<<defaultColor<<std::endl;
    }

    template <int N>
    void GlidePlaneNoiseBase<N>::computeRealNoiseStatistics() const
    {
        // Compute Statistics
        NoiseType ave(NoiseTraits<N>::Zero());
        for(const auto& valArr: noiseVector())
        {
            ave+=valArr;
        }
        ave/=noiseVector().size();
        
        NoiseType var(NoiseTraits<N>::Zero());
        for(const auto& valArr: noiseVector())
        {
            var+=NoiseTraits<N>::squared(valArr-ave);
//            var+= ((valArr-ave).array()*(valArr-ave).array()).matrix();
        }
        var/=noiseVector().size();
        
        std::cout<<"gridSize= "<<gridSize.transpose()<<std::endl;
        std::cout<<"gridSpacing= "<<gridSpacing.transpose()<<std::endl;
        std::cout<<"noiseAverage="<<ave<<std::endl;
        std::cout<<"noiseVariance="<<var<<std::endl;
    }

    template struct GlidePlaneNoiseBase<1>;
    template struct GlidePlaneNoiseBase<2>;

}
#endif


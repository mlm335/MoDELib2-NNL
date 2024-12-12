/* This file is part of MoDELib, the Mechanics Of Defects Evolution Library.
 *
 *
 * MoDELib is distributed without any warranty under the
 * GNU General Public License (GPL) v2 <http://www.gnu.org/licenses/>.
 */

#ifndef model_AnalyticalSolidSolutionNoise_cpp
#define model_AnalyticalSolidSolutionNoise_cpp

#include <AnalyticalSolidSolutionNoise.h>

namespace model
{

AnalyticalSolidSolutionNoiseExpression::AnalyticalSolidSolutionNoiseExpression(const REAL_SCALAR& a_in,const REAL_SCALAR& a_cai_in,const GridSpacingType& gridSpacing,const double& MSSS):
/* init */ a(a_in)
/* init */,a_cai(a_cai_in)
/* init */,stressPrefactor(MSSS*120.0*M_PI*sqrt(M_PI)*a*a*a/(gridSpacing.prod()))
{
//    std::cout<<"stressPrefactor="<<stressPrefactor<<std::endl;
//    std::cout<<"a="<<a<<std::endl;
//    std::cout<<"a_cai="<<a_cai<<std::endl;    
}

// Cai doubly-convoluted spreading function in Fourier space
typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::Wk_Cai(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz, REAL_SCALAR a)
{
    REAL_SCALAR k = sqrt(kx*kx + ky*ky + kz*kz);
    if(k>0)
    {
        return a*k*sqrt(0.5*boost::math::cyl_bessel_k(2,a*k));
//            return a*k*sqrt(0.5*std::cyl_bessel_k(2,a*k));
    }
    else
    {
        return 1.;
    }
}

// Cai spreading function
typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::W_Cai(REAL_SCALAR r2, REAL_SCALAR a)
{
    return 15.*a*a*a*a/(8.*M_PI*pow(r2+a*a,7./2.));
}

typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::W_t_Cai(REAL_SCALAR r2, REAL_SCALAR a)
{
    return 0.3425*W_Cai(r2,0.9038*a) + 0.6575*W_Cai(r2,0.5451*a);
}

// normalized auto-correlation function in Fourier space for sigma_xy
typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::S_xy_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
{
    REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
    return stressPrefactor*(kx*kx*ky*ky)/(k2*k2)*exp(-a*a*k2);
}

// normalized auto-correlation function in Fourier space for sigma_xz
typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::S_xz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
{
    REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
    return stressPrefactor*(kx*kx*kz*kz)/(k2*k2)*exp(-a*a*k2);
}

// normalized auto-correlation function in Fourier space for sigma_yz
typename AnalyticalSolidSolutionNoiseExpression::REAL_SCALAR AnalyticalSolidSolutionNoiseExpression::S_yz_k(REAL_SCALAR kx, REAL_SCALAR ky, REAL_SCALAR kz) const
{
    REAL_SCALAR k2 = kx*kx + ky*ky + kz*kz;
    return stressPrefactor*(ky*ky*kz*kz)/(k2*k2)*exp(-a*a*k2);
}



//#ifdef _MODEL_GLIDE_PLANE_NOISE_GENERATOR_

AnalyticalSolidSolutionNoise::AnalyticalSolidSolutionNoise(const PolycrystallineMaterialBase& mat,
                                                           const std::string& tag,
                                                           const int& seed,
                                                           const GridSizeType& gridSize,
                                                           const GridSpacingType& gridSpacing,
                                                           const double& a_in,
                                                           const double& a_Cai_in,
                                                           const double& MSSS) :
/*init*/ AnalyticalSolidSolutionNoiseExpression(a_in,a_Cai_in,gridSpacing,MSSS)
/*init*/,GlidePlaneNoiseBase<2>("AnalyticalSolidSolutionNoise"+tag,seed,gridSize,gridSpacing)
///*init*/,NX(_gridSize(0))     // dimension along x
///*init*/,NY(_gridSize(1))     // dimension along y
///*init*/,NZ(64)      // dimension along z
///*init*/,DX(_gridSpacing_A(0))     // grid spacing [AA]
///*init*/,DY(_gridSpacing_A(1))     // grid spacing [AA]
///*init*/,DZ(_gridSpacing_A(1))     // grid spacing [AA]
///*init*/,a(a_in)      // spreading length for stresses [AA]
///*init*/,a_cai(a_Cai_in)
///*init*/,a_cai(DislocationFieldBase<3>::a*mat.b_SI*1e10)  // spreading length for non-singular dislocaion theory [AA]
///*init*/,seed(TextFileParser(mat.materialFile).readScalar<double>("seed",true))  // random seed
///*init*/,LX(NX*DX)
///*init*/,LY(NY*DY)
///*init*/,LZ(NZ*DZ)
///*init*/,DV(DX*DY*DZ)
///*init*/,NR(NX*NY*NZ)
///*init*/,NK(NX*NY*(NZ/2+1))
///*init*/,Norm(1./REAL_SCALAR(NR))
{
    
    this->computeRealNoise();
    this->computeRealNoiseStatistics();

    
//    std::cout<<"Computing SolidSolutionNoise..."<<std::endl;
//
//    std::cout << "a_cai = " << a_cai << std::endl;
//    std::cout << "a = " << a << std::endl;
//        const double MSSS_SI(TextFileParser(mat.materialFile).readScalar<double>("MSSS_SI",true));
//        const double MSS(std::sqrt(MSSS_SI)/mat.mu_SI);
//
//    
//    // fftw_plan plan_R_yz_r2c, plan_R_xz_r2c;  // fft plan
//    fftw_plan plan_R_yz_c2r, plan_R_xz_c2r;  // fft plan
//
//    // allocate
//    REAL_SCALAR *Rr_yz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//    COMPLEX *Rk_yz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
//    REAL_SCALAR *Rr_xz = (REAL_SCALAR*) fftw_malloc(sizeof(REAL_SCALAR)*NR);
//    COMPLEX *Rk_xz = (COMPLEX*) fftw_malloc(sizeof(COMPLEX)*NK);
//
//
//    // prepare plans
//    //    plan_W_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Wr, reinterpret_cast<fftw_complex*>(Wk), FFTW_ESTIMATE);
//    //    plan_W_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Wk), Wr, FFTW_ESTIMATE);
//    // plan_R_yz_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Rr_yz, reinterpret_cast<fftw_complex*>(Rk_yz), FFTW_ESTIMATE);
//    // plan_R_xz_r2c = fftw_plan_dft_r2c_3d(NX, NY, NZ, Rr_xz, reinterpret_cast<fftw_complex*>(Rk_xz), FFTW_ESTIMATE);
//    plan_R_yz_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Rk_yz), Rr_yz, FFTW_ESTIMATE);
//    plan_R_xz_c2r = fftw_plan_dft_c2r_3d(NX, NY, NZ, reinterpret_cast<fftw_complex*>(Rk_xz), Rr_xz, FFTW_ESTIMATE);
//
//
//    std::default_random_engine generator(seed);
//    std::normal_distribution<REAL_SCALAR> distribution(0.0,1.0);
//
//
//    /////////////////////////////////////////////////////
//    // generate yz and xz correlated component //
//    /////////////////////////////////////////////////////
//    for(int i=0; i<NX; i++)
//    {
//        for(int j=0; j<NY; j++)
//        {
//            for(int k=0; k<(NZ/2+1); k++)
//            {
//                const int ind = NY*(NZ/2+1)*i + j*(NZ/2+1) + k;
//
//                REAL_SCALAR kx = 2.*M_PI/LX*REAL_SCALAR(i);
//                if(i>NX/2)
//                {
//                    kx = 2.*M_PI/LX*REAL_SCALAR(i-NX);
//                }
//
//                REAL_SCALAR ky = 2*M_PI/LY*REAL_SCALAR(j);
//                if(j>NY/2)
//                {
//                    ky = 2.*M_PI/LY*REAL_SCALAR(j-NY);
//                }
//
//                REAL_SCALAR kz = 2.*M_PI/LZ*REAL_SCALAR(k);
//
//                // random numbers
//                REAL_SCALAR Nk_yz = distribution(generator);
//                REAL_SCALAR Mk_yz = distribution(generator);
//                REAL_SCALAR Nk_xz, Mk_xz;
//                if(kx*ky>=0)
//                {
//                    Nk_xz = Nk_yz;
//                    Mk_xz = Mk_yz;
//                }
//                else
//                {
//                    Nk_xz = -Nk_yz;
//                    Mk_xz = -Mk_yz;
//                }
//
//                if(k==0) // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
//                {
//                    Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                }
//                else if(k==NZ/2)
//                {
//                    Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                }
//                else
//                {
//                    Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz)/2.)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//                    Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz)/2.)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//                }
//
//                if(a_cai>0)
//                {
//                    Rk_yz[ind] = Rk_yz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//                    Rk_xz[ind] = Rk_xz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//                }
//                
//            }
//        }
//    }
//    Rk_yz[0] = 0;
//    Rk_xz[0] = 0;
//
//
//    // FFT back to REAL_SCALAR space
//    fftw_execute(plan_R_yz_c2r);
//    fftw_execute(plan_R_xz_c2r);
//
//
//    this->reserve(NX*NY);
//    for(int i=0;i<NX;i++)
//    {
//        for(int j=0;j<NY;j++)
//        {
//            const int k=0;
//            const int ind = NY*NZ*i + j*NZ + k;
//            this->push_back(MSS*(NoiseType()<<Rr_xz[ind],Rr_yz[ind]).finished());
//            // std::cout<<"Rr_xz[ind]="<<Rr_xz[ind]<<", Rr_yz[ind]="<<Rr_yz[ind]<<std::endl;
//        }
//    }
//    
//    // ouput vtk files
//    const std::string fileName_xz(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser(mat.materialFile).readString("solidSolutionNoiseFile_xz",true));
//    const std::string fileName_yz(std::filesystem::path(mat.materialFile).parent_path().string()+"/"+TextFileParser(mat.materialFile).readString("solidSolutionNoiseFile_yz",true));
//    std::cout<<"Writing noise file "<<fileName_xz<<std::endl;
//    Write_field_slice(Rr_xz, fileName_xz.c_str());
//    std::cout<<"Writing noise file "<<fileName_yz<<std::endl;
//    Write_field_slice(Rr_yz, fileName_yz.c_str());
}



std::array<AnalyticalSolidSolutionNoise::COMPLEX,2> AnalyticalSolidSolutionNoise::kCorrelations(const Eigen::Matrix<double,3,1>& kv,const Eigen::Matrix<int,3,1>&) const
{
    
//    if(k==0) // /!\ special case for k=0 and k==NZ/2 because of folding of C2R Fourier transform
//    {
//        Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//        Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//    }
//    else if(k==NZ/2)
//    {
//        Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz))*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//        Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz))*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//    }
//    else
//    {
//        Rk_yz[ind] = sqrt(S_yz_k(kx,ky,kz)/2.)*(Nk_yz+Mk_yz*COMPLEX(0.0,1.0));
//        Rk_xz[ind] = sqrt(S_xz_k(kx,ky,kz)/2.)*(Nk_xz+Mk_xz*COMPLEX(0.0,1.0));
//    }
//
//    if(a_cai>0)
//    {
//        Rk_yz[ind] = Rk_yz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//        Rk_xz[ind] = Rk_xz[ind]*Wk_Cai(kx, ky, kz, a_cai);
//    }
    
    std::array<AnalyticalSolidSolutionNoise::COMPLEX,2> temp{this->S_xz_k(kv(0),kv(1),kv(2)),this->S_yz_k(kv(0),kv(1),kv(2))};
//    std::cout<<"this->S_xz_k(kv(0),kv(1),kv(2))="<<this->S_xz_k(kv(0),kv(1),kv(2))<<std::endl;
//    std::cout<<"this->S_yz_k(kv(0),kv(1),kv(2))="<<this->S_yz_k(kv(0),kv(1),kv(2))<<std::endl;
//    std::cout<<"kv="<<kv.transpose()<<", temp="<<temp[0]<<","<<temp[1]<<std::endl;
    if(a_cai>0.0)
    {
//        const double wkc(Wk_Cai(kv(0),kv(1),kv(2), a_cai));
        const double wkc2(std::pow(Wk_Cai(kv(0),kv(1),kv(2), a_cai),2)); // using the square because this is before the square root

        temp[0]*=wkc2;
        temp[1]*=wkc2;
    }
    return temp;
}


//    void AnalyticalSolidSolutionNoise::Write_field_slice(REAL_SCALAR *F, const char *fname)
//    {
//        FILE *OutFile=fopen(fname,"w");
//        
//        fprintf(OutFile,"# vtk DataFile Version 2.0\n");
//        fprintf(OutFile,"iter %d\n",0);
//        fprintf(OutFile,"BINARY\n");
//        fprintf(OutFile,"DATASET STRUCTURED_POINTS\n");
//        fprintf(OutFile,"ORIGIN \t %f %f %f\n",0.,0.,0.);
//        fprintf(OutFile,"SPACING \t %f %f %f\n", DX, DY, DZ);
//        fprintf(OutFile,"DIMENSIONS \t %d %d %d\n", NX, NY, 1);
//        fprintf(OutFile,"POINT_DATA \t %d\n",NX*NY);
//        fprintf(OutFile,"SCALARS \t volume_scalars double 1\n");
//        fprintf(OutFile,"LOOKUP_TABLE \t default\n");
//        
//        for(int i=0;i<NX;i++)
//        {
//            for(int j=0;j<NY;j++)
//            {
//                const int k=0;
//                const int ind = NY*NZ*i + j*NZ + k;
//                const double temp=NoiseTraitsBase::ReverseDouble(double(F[ind]));
//                fwrite(&temp, sizeof(double), 1, OutFile);
//            }
//        }
//        
//        fclose(OutFile);
//    }


}
#endif


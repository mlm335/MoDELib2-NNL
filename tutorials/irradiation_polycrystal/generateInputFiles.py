import sys
sys.path.append("../../python/")
from modlibUtils import *
import random
from datetime import datetime

# Create folder structure
folders=['evl','F','inputFiles']
for x in folders:
    if not os.path.exists(x):
        os.makedirs(x)

# Make a local copy of DD parameters file and modify that copy if necessary
DDfile='DD.txt'
DDfileTemplate='../../Library/DislocationDynamics/'+DDfile
print("\033[1;32mCreating  DDfile\033[0m")
shutil.copy2(DDfileTemplate,'inputFiles/'+DDfile)
setInputVariable('inputFiles/'+DDfile,'useFEM','1')
setInputVariable('inputFiles/'+DDfile,'useDislocations','1')
setInputVariable('inputFiles/'+DDfile,'useInclusions','0')
setInputVariable('inputFiles/'+DDfile,'useElasticDeformation','1')
setInputVariable('inputFiles/'+DDfile,'useClusterDynamics','1')
setInputVariable('inputFiles/'+DDfile,'Nsteps','6')  # number of simulation steps
setInputVariable('inputFiles/'+DDfile,'timeSteppingMethod','fixed') # adaptive or fixed
setInputVariable('inputFiles/'+DDfile,'dtMax','26.47935e18')
setInputVariable('inputFiles/'+DDfile,'dxMax','5') # max nodal displacement for when timeSteppingMethod=adaptive
setInputVariable('inputFiles/'+DDfile,'use_velocityFilter','0') # don't filter velocity if noise is enabled
setInputVariable('inputFiles/'+DDfile,'use_stochasticForce','0') # Langevin thermal noise enabled
setInputVariable('inputFiles/'+DDfile,'glideSolverType','none')  # type of glide solver, or none
setInputVariable('inputFiles/'+DDfile,'climbSolverType','Galerkin')  # type of clim solver, or none
setInputVariable('inputFiles/'+DDfile,'quadPerLength','0.0001')  # quadrature points per unit length (in Burgers vector)
setInputVariable('inputFiles/'+DDfile,'alphaLineTension','1.0') # dimensionless scale factor in for line tension forces
setInputVariable('inputFiles/'+DDfile,'remeshFrequency','1')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmin','450')  # min segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'Lmax','500')  # max segment length (in Burgers vector units)
setInputVariable('inputFiles/'+DDfile,'outputFrequency','1')  # output frequency
setInputVariable('inputFiles/'+DDfile,'outputQuadraturePoints','1')  # output quadrature data

# Make a local copy of material file, and modify that copy if necessary
materialFile='Zr4.txt';
materialFileTemplate='../../Library/Materials/'+materialFile;
print("\033[1;32mCreating  materialFile\033[0m")
shutil.copy2(materialFileTemplate,'inputFiles/'+materialFile)
setInputVariable('inputFiles/'+materialFile,'enabledSlipSystems','fullBasal fullPrismatic')

# Make a local copy of ElasticDeformation file, and modify that copy if necessary
elasticDeformatinoFile='ElasticDeformation.txt';
elasticDeformatinoFileTemplate='../../Library/ElasticDeformation/'+elasticDeformatinoFile;
print("\033[1;32mCreating  elasticDeformatinoFile\033[0m")
shutil.copy2(elasticDeformatinoFileTemplate,'inputFiles/'+elasticDeformatinoFile)
setInputVector('inputFiles/'+elasticDeformatinoFile,'ExternalStress0',np.array([0.0,0.0,0.0,0.0,0.0,0.0]),'applied stress')

# Create polycrystal.txt using local material file
meshFile='poly30_100K.msh';
meshFileTemplate='../../Library/Meshes/'+meshFile;
print("\033[1;32mCreating  polycrystalFile\033[0m")
shutil.copy2(meshFileTemplate,'inputFiles/'+meshFile)
pf=PolyCrystalFile(materialFile);
pf.absoluteTemperature=553;
pf.meshFile=meshFile

pf.singleCrystal=False
pf.f_param=np.array([0.63,0.32])
pf.numberGrains=30

#pf.grain1globalX1=np.array([1,0,0])     # global x1 axis. Overwritten if alignToSlipSystem0=true
#pf.grain1globalX3=np.array([0,0,1])    # global x3 axis. Overwritten if alignToSlipSystem0=true
#pf.boxEdges=np.array([[1,0,0],[0,1,0],[0,0,1]]) # i-throw is the direction of i-th box edge

pf.boxScaling=np.array([3093,3093,3093]) # must be a vector of integers
pf.X0=np.array([0,0,0]) # Centering unitCube mesh. Mesh nodes X are mapped to x=F*(X-X0)
pf.periodicFaceIDs=np.array([])
pf.write('inputFiles')


# Generate a Texture for Polycrystals
#f_param = [0.63,0.32] #PWR radial, transverse, axial f-params
##f_param = [0.32,0.61] #CANDU radial, transverse, axial f-params
#def random_output(numbers):
#    random.seed(datetime.now().timestamp())
#    j = random.random()
#    if j < f_param[0]:
#            return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) #crystal c align with global x
#            #return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) #All crystal c align with global x
#    elif j > f_param[0] and j < f_param[0]+f_param[1]:
#            return np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]) #crystal c align with global y
#            #return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) #All crystal c align with global x
#    else:
#            return np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]) #crystal c align with global z
#            #return np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]]) #All crystal c align with global x
#N = 30
#C2G = []
#for i in range(N):
#    R = random_output(i)
#    C2G.append(R)
#with open('poly'+str(N)+'.txt', 'w') as fID:
#    for i in range(N):
#        fID.write(f'C2G{i + 1} =\n')
#        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][0]]) + '\n')
#        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][1]]) + '\n')
#        fID.write(' '.join([f'{val:.15f}' for val in C2G[i][2]]) + ';\n\n')



# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile1='frankLoopsDensity.txt';
microstructureFileTemplate1='../../Library/Microstructures/'+microstructureFile1;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate1,'inputFiles/'+microstructureFile1) # target filename is /dst/dir/file.ext
setInputVariable('inputFiles/'+microstructureFile1,'targetDensity','1e13')
setInputVector('inputFiles/'+microstructureFile1,'planeIDs',np.array([0]),'')
setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionMean','50e-9') # [m] mean of loop radii
setInputVariable('inputFiles/'+microstructureFile1,'radiusDistributionStd','0e-8') # [m] std of loop radii
setInputVariable('inputFiles/'+microstructureFile1,'numberOfSides','20') # [-] number of sides in each loop polygon
setInputVariable('inputFiles/'+microstructureFile1,'burgersFactor','0.5') # [-] scaling of burgers vector
setInputVariable('inputFiles/'+microstructureFile1,'areVacancyLoops','1') # 1=vacancy-type loops, 0=interstiatial-type loops

# make a local copy of microstructure file, and modify that copy if necessary
microstructureFile2='aLoopsDensity.txt';
microstructureFileTemplate2='../../Library/Microstructures/'+microstructureFile2;
print("\033[1;32mCreating  microstructureFile\033[0m")
shutil.copy2(microstructureFileTemplate2,'inputFiles/'+microstructureFile2) # target filename is /dst/dir/file.ext
setInputVector('inputFiles/'+microstructureFile2,'slipSystemIDs',np.array([6,8,10,6,8,10]),'')
setInputVector('inputFiles/'+microstructureFile2,'targetDensity',np.array([1e20,1e20,1e20,1e20,1e20,1e20]),'')
setInputVector('inputFiles/'+microstructureFile2,'loopRadiusMean',np.array([10e-9,10e-9,10e-9,10e-9,10e-9,10e-9]),'')
setInputVector('inputFiles/'+microstructureFile2,'loopRadiusStd',np.array([0,0,0,0,0,0]),'')
setInputVector('inputFiles/'+microstructureFile2,'numberOfSides',np.array([36,36,36,36,36,36]),'')
setInputVector('inputFiles/'+microstructureFile2,'areVacancyLoops',np.array([1,1,1,0,0,0]),'')
setInputVector('inputFiles/'+microstructureFile2,'ellipticityFactor',np.array([1.5,1.5,1.5,1.5,1.5,1.5]),'')

print("\033[1;32mCreating  initialMicrostructureFile\033[0m")
with open('inputFiles/initialMicrostructure.txt', "w") as initialMicrostructureFile:
    initialMicrostructureFile.write('microstructureFile='+microstructureFile1+';\n')
    initialMicrostructureFile.write('microstructureFile='+microstructureFile2+';\n')

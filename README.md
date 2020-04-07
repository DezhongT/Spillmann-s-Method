# Spillmann-s-Method
DER codes for simulating overhand knots with frictional contact based on spillmann's method

Instructions for Ubuntu:
(1) To run this code you need Eigen, OpenGL, Pardiso and Lapack. Lapack is usually preinstalled on your computer. 
Eigen can be found at http://eigen.tuxfamily.org/index.php?title=Main_Page

(2) Compile Command:
e.g.  g++ -I /usr/local/include/eigen3/ main.cpp world.cpp elasticRod.cpp elasticStretchingForce.cpp elasticBendingForce.cpp elasticTwistingForce.cpp externalGravityForce.cpp inertialForce.cpp dampingForce.cpp timeStepper.cpp setInput.cpp collision.cpp -llapack -lGL -lglut -lGLU -Ofast -o simDER  

(3) To start the simulation, run the command "./simDER option_ploy.txt" or "./simDER option_PRL.txt"  (without the quotes). More on option.txt later.  

(4) Before running the simulation, the initial configuration file should be putted into the folder. The initial configuration file should be renamed as 'simDER.txt'. Initial configurations for overhands with different unknotting number can be found in the \intial_configurations.  

Simulating Physical process:
Holding two ends of a overhand knots to keep static for a while, then pulling two ends away slowly. During this process, the configurations and relatinonship between forces and shortening will be recoreded.  

 Details on the options (we use SI units): 
    "RodLength" is the length of a rod for tying the overhand knot.  
    "rodRadius" is the cross-sectional radius of the rod.  
    "youngM" is the young's modulus.  
    "Poisson" is the Poisson ratio.  
    "deltaTime" is the time step size.  
    "totalTime" is the time at which the simulation ends.  
    "tol" and "stol" are small numbers used in solving the linear system. Fraction of a percent, e.g. 1.0e-3, is often a good choice.  
    "maxIter" is the maximum number of iterations allowed before the solver quits.  
    "density" is the mass per unit volume.  
    "gVector" is the vector specifying acceleration due to gravity.  
    "viscosity" is the viscosity of the fluid medium.  
    "numVertices" is the number of nodes on the rod.  
    "render" (0 or 1) indicates whether OpenGL visualization should be displayed.  
    "saveData" (0 or 1) indicates whether the location of the head should be saved in "datafiles/" folder (this folder will be created by the program).  
    "pulltime" is the duration for pulling two ends.  
    "friction" is the coefficient of friction.  
    "statictime" is the duration to keep the rod static at the begining of the simulation.   



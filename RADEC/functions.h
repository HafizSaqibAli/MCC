/* 
** Created 13/02/2014
** Modified 15/07/2016
**
** List of all included functions to aid Main.cpp readability
**
*/

//AMBER Input functions 
#include "functions/readTopology.h"
#include "functions/readTrajectory.h"
#include "functions/readForces.h"
#include "functions/getSigmas.h"

//Housekeeping functions
#include "functions/outputVector.h"
#include "functions/createFragments.h"
#include "functions/fragmentTypes.h"
#include "functions/createMolecules.h"
#include "functions/moleculeTypes.h"
#include "functions/createMoleculeShells.h"
#include "functions/outputMolecules.h"
#include "functions/sortMolecules.h"
#include "functions/sortFragments.h"

//Force-Torque Covariance Functions
#include "functions/localiseForces.h"
#include "functions/FTCMOperations.h"
#include "functions/FTCMOutputs.h"

//g(r) functions
//#include "functions/nrCalculation.h"
//#include "functions/nrShellCalculation.h"
//#include "functions/grCalculation.h"
//#include "functions/smoothCurve.h"
//#include "functions/outputgrCurve.h"

//Shell functions
#include "functions/calculateRADFragmentShells.h"
#include "functions/calculateRADMoleculeShells.h"
//#include "functions/cutoffShells.h"
//#include "functions/updateFragmentShellTypes.h"
//#include "functions/updateMoleculeShellTypes.h"
#include "functions/updateFragmentShellTypeForces.h"
#include "functions/updateMoleculeShellTypeForces.h"
#include "functions/outputFragmentShellForces.h"
#include "functions/outputMoleculeShellForces.h"
//#include "functions/outputShells.h"

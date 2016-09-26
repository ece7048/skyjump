/* -------------------------------------------------------------------------- *
 *                  OpenSim:  JumpingOptimizationDouble.cpp                   *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Carmichael Ong, Matt Titchenal                                  *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/* 
 *  This was created for a class project at Stanford for ME 485. 
 *  
 *  Look in JumpingOptimization.cpp for more information. This version of the code
 *  is similar to JumpingOptimization.cpp, but it loops through each of the muscles
 *  and doubles MaxIsometricForce and doubles MaxContractionVelocity before sending
 *  the model to the optimizer. To do either one of these tasks, the code for the other
 *  loop was commented out and the .exe was rebuilt.
 *  
 *  
 */

//==============================================================================
//==============================================================================
#include <OpenSim/OpenSim.h>
#include <OpenSim/Common/PiecewiseConstantFunction.h>
#include <ctime>  // clock(), clock_t, CLOCKS_PER_SEC

using namespace OpenSim;
using namespace SimTK;
using namespace std;


// Global variables
int stepCount = 0;
double initialTime = 0.0;
double finalTime = 1.0;
double bestSoFar = Infinity;
double forceThreshold = 1.0;
double settleTime = 0.3;

ofstream optLog("optLog.txt", ofstream::out);
ofstream bestSoFarLog("bestSoFarLog.txt", ofstream::out);
ofstream verboseLog("verboseLog.txt", ofstream::out);
ofstream debugLog("debugLog.txt", ofstream::out);
std::clock_t optimizerStartTime;

//=============================================================================
// EVENT HANDLER: TerminateSimulation
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulation : public SimTK::TriggeredEventHandler {
public:
    
    // CONSTRUCTOR
    TerminateSimulation(const Model& m, double threshold) : TriggeredEventHandler(Stage::Dynamics), 
        _model(m), _forceThreshold(threshold)
    {
        getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    // WITNESS FUNCTION
    SimTK::Real getValue(const SimTK::State& s) const
	{

		// keep returning a positive value while the contact settles out
		double simTime = s.getTime();
		if (simTime < settleTime) {
			return 100.0;
		}

		// total vertical force applied to both feet
		Array<double> force_foot_r = _model.getForceSet().get("foot_r").getRecordValues(s);
		Array<double> force_foot_l = _model.getForceSet().get("foot_l").getRecordValues(s);
		double netGRFVertical = -force_foot_r[1] - force_foot_l[1];
		return netGRFVertical - _forceThreshold;
	}

    // EVENT HANDLER FUNCTION
    void handleEvent (SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
	{
		//cout << "terminated because vertical force passed threshold" << endl;
		terminate = true;
	}


private:
    const Model& _model;
    double _forceThreshold;
};


//=============================================================================
// OPTIMIZER SYSTEM: JumpingOptimizationSystem
// Defines a constructor and objective function for the optimization
//=============================================================================
class JumpingOptimizationSystem : public OptimizerSystem {
   public:

	   /* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */
	   JumpingOptimizationSystem(int numParameters, Model& aModel): 
             OptimizerSystem(numParameters), osimModel(aModel) {
			 
			// Print a heading for the udpates.
			bestSoFarLog << "iter. count" << "\t" << "ObjectiveValue" << "\t"
			<< "Parameters" << "\t" << "TimeSoFar" << endl;

			verboseLog << "iter. count" << "\t" << "ObjectiveValue" << "\t"
            << "Parameters" << "\t" << "TimeSoFar" << endl;

			optimizerStartTime = std::clock();	
		}
			 	
	/* Objective Function. */
	int objectiveFunc(  const Vector &newControllerParameters, bool new_coefficients, Real& f ) const {

		// Grab actuator and controller sets
		const Set<Actuator> &actuatorSet = osimModel.getActuators();
		int numActuators = actuatorSet.getSize();
		const ControllerSet &controllerSet = osimModel.getControllerSet();
		int numControllers = controllerSet.getSize();
		// Update the control values
		for (int i = 0; i < numControllers/2; i++) {

			// define a piecewise constant function for both sides
			PiecewiseConstantFunction bangBangControl;
			bangBangControl.addPoint(initialTime, 0);
			bangBangControl.addPoint(newControllerParameters[2*i],1);
			bangBangControl.addPoint(newControllerParameters[2*i]+newControllerParameters[2*i+1],0);

			// Update the right side
			PrescribedController* controller_r = dynamic_cast<PrescribedController*>( &controllerSet.get(2*i) );
			controller_r->prescribeControlForActuator(actuatorSet.get(i).getName(), bangBangControl.clone());

			// Update the left side
			PrescribedController* controller_l = dynamic_cast<PrescribedController*>( &controllerSet.get(2*i+1) );
			controller_l->prescribeControlForActuator(actuatorSet.get(i + numActuators/2).getName(), bangBangControl.clone());
			
		}

		// Initialize the system. initSystem() cannot be used here because adding the event handler
		// must be done between buildSystem() and initializeState().
		osimModel.buildSystem();
		TerminateSimulation *terminate = new TerminateSimulation(osimModel, forceThreshold);
		osimModel.updMultibodySystem().addEventHandler(terminate);
		State &osimState = osimModel.initializeState();

		// Set the initial muscle activations 
		const Set<Muscle> &muscleSet = osimModel.getMuscles();
     	for(int i=0; i< muscleSet.getSize(); i++ ){
			muscleSet[i].setActivation(osimState, 0.05);
		}
	
		// Make sure the muscles states are in equilibrium
		osimModel.equilibrateMuscles(osimState);

		// Create the integrator for the simulation.
		RungeKuttaMersonIntegrator integrator(osimModel.getMultibodySystem());
		integrator.setAccuracy(1.0e-6); //1.0e-6 seems to converge faster (both objFunc calls and time) than 1.0e-5

		// Create a manager to run the simulation. Can change manager options to save run time and memory or print more information
		Manager manager(osimModel, integrator);
		//manager.setWriteToStorage(false);
		manager.setPerformAnalyses(false);

		// Integrate from initial time to final time and integrate
		manager.setInitialTime(initialTime);
		manager.setFinalTime(finalTime);
		manager.integrate(osimState);

		/* Calculate the scalar quantity we want to minimize or maximize. 
		*  In this case, we’re maximizing the height of the COM of the jumper
		*  so to maximize, calculate (position + velocity^2)/(2g) when simulation ends.
		*  Then take the negative of that quantity (since optimizers minimize).
		*/
		osimModel.getMultibodySystem().realize(osimState, Stage::Velocity);
		Vec3 COM_position = osimModel.getMultibodySystem().getMatterSubsystem().calcSystemMassCenterLocationInGround(osimState);
		Vec3 COM_velocity = osimModel.getMultibodySystem().getMatterSubsystem().calcSystemMassCenterVelocityInGround(osimState);
		double g = -osimModel.getGravity()[1];
		double maxHeight = COM_position[1] + pow(COM_velocity[1], 2.0)/(2.0*g);
		f = -maxHeight;

		stepCount++;
		
		// Store and print the results of a "random sample"
		if( stepCount == 23 ){
			//manager.getStateStorage().print("Jumper_randomSample_states.sto");
			osimModel.print("Jumper_randomSample.osim");
			//osimModel.printControlStorage("Arm26_randomSample_controls.sto");
		}
		// Store and print the  results of the first step.
		else if( stepCount == 1){ 
			manager.getStateStorage().print("Jumper_noActivation_states.sto");
			osimModel.print("Jumper_noActivation.osim");
			osimModel.printControlStorage("Jumper_noActivation_controls.sto");
		}
		// Use an if statement to only store and print the results of an 
		//  optimization step if it is better than a previous result.
		else if( f < bestSoFar){
			manager.getStateStorage().print("Jumper_bestSoFar_states.sto");
			osimModel.print("Jumper_bestSoFar.osim");
			//osimModel.printControlStorage("Arm26_bestSoFar_controls.sto");
			bestSoFar = f;
			bestSoFarLog << stepCount << "\t" << f << "\t"
				<< newControllerParameters << "\t" << (std::clock()-optimizerStartTime)/CLOCKS_PER_SEC << endl;
			
			verboseLog << stepCount << "\t" << f << "\t"
				<< newControllerParameters << "\t" << (std::clock()-optimizerStartTime)/CLOCKS_PER_SEC << endl;
		}
		// Print every 100 objective function calls.
		else if( stepCount % 100 == 0) {
			verboseLog << stepCount << "\t" << f << "\t"
				<< newControllerParameters << "\t" << (std::clock()-optimizerStartTime)/CLOCKS_PER_SEC << endl;
		}

      return(0);

   }	

private:
	Model& osimModel;
 };



//______________________________________________________________________________
/**
 * Define an optimization problem that finds a set of muscle controls to maximize 
 * the vertical jump height of a model.
 */
int main()
{
	try {
		std::clock_t startTime = std::clock();	

		// Ensures all components are printed out to the model file.
		Object::setSerializeAllDefaults(true);
	
		// Create a new OpenSim model from file
		Model osimModel("jumper10dof24musc.osim");

		/* Get the muscle set, then double MaxIsometricForce and/or MaxContractionVelocity for each muscle. */
		const Set<Muscle> &muscleSet = osimModel.getMuscles();

		//// Double maximum isometric force of all muscles: COMMENT OUT TWO LINES BELOW IF YOU DO NOT WANT TO DO THIS
		//for (int i=0; i < muscleSet.getSize(); i++) {
		//	muscleSet[i].setMaxIsometricForce(muscleSet[i].getMaxIsometricForce()*2);
		//}

		// Double maximum contraction velocity of all muscles: COMMENT OUT TWO LINES BELOW IF YOU DO NOT WANT TO DO THIS
		for (int i=0; i < muscleSet.getSize(); i++) {
			muscleSet[i].setMaxContractionVelocity(muscleSet[i].getMaxContractionVelocity()*2);
		}

		// The number of parameters is the number of actuators
		const Set<Actuator> &actuatorSet = osimModel.getActuators();
		int numActuators = actuatorSet.getSize();
		optLog << "numActuators = " << numActuators << endl;
		int numParameters = numActuators;

		/* Define initial values for controllerParameters. Each controller has two parameters, an
		   initial time, and a duration for bang-bang control. There are as many controls as
		   there are actuators because we assume symmetry (but each controller has two parameters).
		   controllerParameters = [ti_0 duration_0 ti_1 duration_1 ... ]' */
		
		Vector controllerParameters(numParameters, 0.05);

		// initialize initial times
		controllerParameters[0] = 0.107863; //hamstrings
		controllerParameters[2] = 0.00069462; //bifmesh
		controllerParameters[4] = 0.296889; //glut_max
		controllerParameters[6] = 0.0533121; //iliopsoas
		controllerParameters[8] = 0.166427; //rect_fem
		controllerParameters[10] = 0.277592; //vasti
		controllerParameters[12] = 0.344144; //gastroc
		controllerParameters[14] = 0.315376; //soleus
		controllerParameters[16] = 0.0857591; //tib_ant
		controllerParameters[18] = 0.149644; //ercspn
		controllerParameters[20] = 0.468741; //intobl
		controllerParameters[22] = 0.455184; //extobl

		// initialize durations
		controllerParameters[1] = 0.429072; //hamstrings
		controllerParameters[3] = 0.120626; //bifemsh
		controllerParameters[5] = 0.39246; //glut_max
		controllerParameters[7] = 0.0303192; //iliopsoas
		controllerParameters[9] = 0.370385; //rect_fem
		controllerParameters[11] = 0.3; //vasti
		controllerParameters[13] = 0.343613; //gastroc
		controllerParameters[15] = 0.350808; //soleus
		controllerParameters[17] = 0.0277077; //tib_ant
		controllerParameters[19] = 0.243332; //ercspn
		controllerParameters[21] = 0.160016; //intobl
		controllerParameters[23] = 0.15192; //extobl

		// Add prescribed controllers to each muscle. Need to only loop numActuators/2 times since we are enforcing symmetry.
		// It is assumed that all actuators are listed for one side of the model first, then the other side, in the same order.
		for(int i=0; i < numActuators/2; i++) {

			// make a piecewise constant function for both sides
			PiecewiseConstantFunction bangBangControl;
			bangBangControl.addPoint(initialTime,0);
			bangBangControl.addPoint(controllerParameters[2*i],1);
			bangBangControl.addPoint(controllerParameters[2*i]+controllerParameters[2*i+1],0);

			// add controller to right side
			PrescribedController *actuatorController_r = new PrescribedController();
			Set<Actuator> actuator_r;
			actuator_r.insert(0,osimModel.getActuators().get(i));
			actuatorController_r->setName(actuatorSet.get(i).getName());
			actuatorController_r->setActuators(*actuator_r.clone());
			actuatorController_r->prescribeControlForActuator(osimModel.getActuators().get(i).getName(), bangBangControl.clone());
			osimModel.addController(actuatorController_r);

			// add controller to left side
			PrescribedController *actuatorController_l = new PrescribedController();
			Set<Actuator> actuator_l;
			actuator_l.insert(0,osimModel.getActuators().get(i + numActuators/2));
			actuatorController_l->setName(actuatorSet.get(i + numActuators/2).getName());
			actuatorController_l->setActuators(*actuator_l.clone());
			actuatorController_l->prescribeControlForActuator(osimModel.getActuators().get(i + numActuators/2).getName(), bangBangControl.clone());
			osimModel.addController(actuatorController_l);

		}

		// Create the OptimizationSystem. Initialize the objective function value "f".
		JumpingOptimizationSystem sys(numParameters, osimModel);
		Real f = NaN;
		
		// Set lower and upper bounds.
		Vector lower_bounds(numParameters, initialTime);
		Vector upper_bounds(numParameters, finalTime);

		// Limit the duration of the "bang" to be at least a certain value
		for (int i = 1; i<numParameters; i+=2) {
			lower_bounds[i] = 0.0001;
		}

		sys.setParameterLimits( lower_bounds, upper_bounds );
		
		// Create an optimizer. Pass in our OptimizerSystem
		// and the name of the optimization algorithm.
		optLog << "using LBFGSB" << endl; Optimizer opt(sys, SimTK::LBFGSB); //LBFGSB was found to be beter for this problem
		//optLog << "using IPOPT" << endl; Optimizer opt(sys, InteriorPoint);

		// Specify settings for the optimizer
		opt.setConvergenceTolerance(0.01);
		opt.useNumericalGradient(true);
		opt.setMaxIterations(1000);
		opt.setLimitedMemoryHistory(500);
		//opt.setDiagnosticsLevel(4); // First level that gives outer loop information for IPOPT
			
		// Optimize it!
		optLog << "Optimizing!" << endl;
		optLog << "Initial controllerParameters: " << controllerParameters << endl;
		optLog << "lower_bounds: " << lower_bounds << endl;
		optLog << "upper_bounds: " << upper_bounds << endl;
		f = opt.optimize(controllerParameters);
			
		optLog << "Elapsed time = " << (std::clock()-startTime)/CLOCKS_PER_SEC << "s" << endl;


        optLog << "OpenSim example completed successfully. Look in bestSoFarLog.txt for best solution.\n";

		optLog.close();
		bestSoFarLog.close();
		verboseLog.close();
		debugLog.close();
	}
    catch (const std::exception& ex)
    {
        std::cout << ex.what() << std::endl;
        return 1;
    }
	
	// End of main() routine.
	return 0;
}

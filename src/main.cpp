#include <OpenSim/OpenSim.h>
#include <OpenSim/Common/PiecewiseConstantFunction.h>
#include <OpenSim/Analyses/ForceReporter.h>
#include <OpenSim/Analyses/BodyKinematics.h>
#include <ctime>

#include "Settings.h"
#include "INIReader.h"

using namespace OpenSim;
using namespace SimTK;
using namespace std;

// Global variables
INIReader ini;
int stepCount = 0;
double initialTime;
double finalTime;
double bestSoFar = Infinity;
double forceThreshold = 1.0;
double settleTime = 0.2;
double penaltyWeight;
string resultDir;
ForceReporter* forceAnalysis;
BodyKinematics* bodyKinematics;

ofstream optLog;
ofstream bestSoFarLog;
ofstream verboseLog;
ofstream debugLog;
std::clock_t optimizerStartTime;

//=============================================================================
// EVENT HANDLER: TerminateSimulation
// Triggered when the GRF normal falls below some threshold
//=============================================================================
class TerminateSimulation : public SimTK::TriggeredEventHandler
{
public:

    TerminateSimulation(const Model& m, double threshold) : TriggeredEventHandler(Stage::Dynamics),
        model(m), forceThreshold(threshold)
    {
        getTriggerInfo().setTriggerOnFallingSignTransition(true);
    }

    SimTK::Real getValue(const SimTK::State& s) const
    {
        // keep returning a positive value while the contact settles out
        double simTime = s.getTime();
        if (simTime < settleTime)
        {
            return 100.0;
        }

        // total vertical force applied to both feet
        Array<double> force_foot_r = model.getForceSet().get("foot_r").getRecordValues(s);
        Array<double> force_foot_l = model.getForceSet().get("foot_l").getRecordValues(s);
        double netGRFVertical = -force_foot_r[1] - force_foot_l[1];
        return netGRFVertical - forceThreshold;
    }

    void handleEvent(SimTK::State& s, SimTK::Real accuracy, bool& terminate) const
    {
        //cout << "terminated because vertical force passed threshold" << endl; //debug line
        terminate = true;
    }

private:
    const Model& model;
    double forceThreshold;
};

//=============================================================================
// OPTIMIZER SYSTEM: JumpingOptimizationSystem
// Defines a constructor and objective function for the optimization
//=============================================================================
class JumpingOptimizationSystem : public OptimizerSystem
{
public:

    /* Constructor class. Parameters passed are accessed in the objectiveFunc() class. */
    JumpingOptimizationSystem(int numParameters, Model& aModel) :
        OptimizerSystem(numParameters), osimModel(aModel)
    {
        // Print a heading for the udpates.
        bestSoFarLog << "iter. count" << "\t" << "ObjectiveValue" << "\t"
            << "Parameters" << "\t" << "TimeSoFar" << endl;

        verboseLog << "iter. count" << "\t" << "ObjectiveValue" << "\t"
            << "Parameters" << "\t" << "TimeSoFar" << endl;

        optimizerStartTime = std::clock();
    }

    /* Objective Function. */
    int objectiveFunc(const Vector &newControllerParameters, bool new_coefficients, Real& f) const
    {
        // Grab actuator and controller sets
        const Set<Actuator> &actuatorSet = osimModel.getActuators();
        int numActuators = actuatorSet.getSize();
        const ControllerSet &controllerSet = osimModel.getControllerSet();
        int numControllers = controllerSet.getSize();

        // Update the control values
        for (int i = 0; i < numControllers / 2; i++)
        {
            // define a piecewise constant function for both sides
            PiecewiseConstantFunction bangBangControl;
            bangBangControl.addPoint(initialTime, 0);
            bangBangControl.addPoint(newControllerParameters[2 * i], 1);
            bangBangControl.addPoint(newControllerParameters[2 * i] + newControllerParameters[2 * i + 1], 0);

            // Update the right side
            PrescribedController* controller_r = dynamic_cast<PrescribedController*>(&controllerSet.get(2 * i));
            controller_r->prescribeControlForActuator(actuatorSet.get(i).getName(), bangBangControl.clone());

            // Update the left side
            PrescribedController* controller_l = dynamic_cast<PrescribedController*>(&controllerSet.get(2 * i + 1));
            controller_l->prescribeControlForActuator(actuatorSet.get(i + numActuators / 2).getName(), bangBangControl.clone());
        }

        // Initialize the system. initSystem() cannot be used here because adding the event handler
        // must be done between buildSystem() and initializeState().
        osimModel.buildSystem();
        TerminateSimulation *terminate = new TerminateSimulation(osimModel, forceThreshold);
        osimModel.updMultibodySystem().addEventHandler(terminate);

        // reset analysis
        forceAnalysis->updForceStorage().reset(0);
        bodyKinematics->getPositionStorage()->reset(0);
        bodyKinematics->getAccelerationStorage()->reset(0);
        bodyKinematics->getVelocityStorage()->reset(0);

        State &osimState = osimModel.initializeState();

        // Set the initial muscle activations
        const Set<Muscle> &muscleSet = osimModel.getMuscles();
        for (int i = 0; i < muscleSet.getSize(); i++)
        {
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
        //manager.setPerformAnalyses(false);

        // Integrate from initial time to final time and integrate
        manager.setInitialTime(initialTime);
        manager.setFinalTime(finalTime);
        manager.integrate(osimState);

        //// compute cost function
        //Storage* kin = bodyKinematics->getPositionStorage();
        //Storage& force = forceAnalysis->updForceStorage();

        //Array<double> time, comY, lumbarLimit, kneeRLimit, kneeLLimit;
        //kin->getTimeColumn(time);
        //kin->getDataColumn("center_of_mass_Y", comY);
        //force.getDataColumn("LumbarExtensionLimit", lumbarLimit);
        //force.getDataColumn("KneeLimit_r", kneeRLimit);
        //force.getDataColumn("KneeLimit_l", kneeLLimit);

        //double maxHeight = 0;
        //int maxIndex = 0;
        //for (int i = 0; i < time.size(); i++)
        //{
        //    if (comY[i] > maxHeight)
        //    {
        //        maxHeight = comY[i];
        //        maxIndex = i;
        //    }
        //}

        //// trapezoid integration 0.5 sum (ti - ti-1) * (fi - fi-1)
        //double trapezSum = 0;
        //for (int i = 1; i <= maxIndex; i++)
        //{
        //    double fi = pow(lumbarLimit[i], 2) + pow(kneeRLimit[i], 2) + pow(kneeLLimit[i], 2);
        //    double fminus = pow(lumbarLimit[i - 1], 2) + pow(kneeRLimit[i - 1], 2) + pow(kneeLLimit[i - 1], 2);
        //    trapezSum += (time[i] - time[i - 1]) * (fi + fminus);
        //}
        //// multiply by - because optimization tries to minimize
        //f = -(maxHeight - penaltyWeight * 0.5 * trapezSum);

        //cout << "Max height: " << maxHeight << endl;
        //cout << "Ligaments: " << 0.5 * trapezSum << endl;
        //cout << "Cost: " << -f << endl;

        /* Calculate the scalar quantity we want to minimize or maximize.
        *  In this case, we’re maximizing the height of the COM of the jumper
        *  so to maximize, calculate (position + velocity^2)/(2g) when simulation ends.
        *  Then take the negative of that quantity (since optimizers minimize).
        */
        osimModel.getMultibodySystem().realize(osimState, Stage::Velocity);
        Vec3 COM_position = osimModel.getMultibodySystem().getMatterSubsystem().calcSystemMassCenterLocationInGround(osimState);
        Vec3 COM_velocity = osimModel.getMultibodySystem().getMatterSubsystem().calcSystemMassCenterVelocityInGround(osimState);
        double g = -osimModel.getGravity()[1];
        double maxHeight = COM_position[1] + pow(COM_velocity[1], 2.0) / (2.0*g);
        f = -maxHeight;

        stepCount++;

        // Store and print the results of a "random sample"
        if (stepCount == 23)
        {
            //manager.getStateStorage().print("Jumper_randomSample_states.sto");
            osimModel.print(resultDir + "/Jumper_randomSample.osim");
            //osimModel.printControlStorage("Arm26_randomSample_controls.sto");
        }
        // Store and print the results of the first step.
        else if (stepCount == 1)
        {
            manager.getStateStorage().print(resultDir + "/Jumper_noActivation_states.sto");
            osimModel.print(resultDir + "/Jumper_noActivation.osim");
            osimModel.printControlStorage(resultDir + "/Jumper_noActivation_controls.sto");
            forceAnalysis->printResults("Jumper_no", resultDir);
            bodyKinematics->printResults("Jumper_no", resultDir);
        }
        // Use an if statement to only store and print the results of an
        //  optimization step if it is better than a previous result.
        else if (f < bestSoFar)
        {
            manager.getStateStorage().print(resultDir + "/Jumper_bestSoFar_states.sto");
            osimModel.print(resultDir + "/Jumper_bestSoFar.osim");
            osimModel.printControlStorage(resultDir + "/Jumper_bestSoFar_controls.sto");
            forceAnalysis->printResults("Jumper_best", resultDir);
            bodyKinematics->printResults("Jumper_best", resultDir);
            bestSoFar = f;
            bestSoFarLog << stepCount << "\t" << f << "\t"
                << newControllerParameters << "\t" << (std::clock() - optimizerStartTime) / CLOCKS_PER_SEC << endl;

            verboseLog << stepCount << "\t" << f << "\t"
                << newControllerParameters << "\t" << (std::clock() - optimizerStartTime) / CLOCKS_PER_SEC << endl;
        }
        // Print every 100 objective function calls.
        else if (stepCount % 100 == 0)
        {
            verboseLog << stepCount << "\t" << f << "\t"
                << newControllerParameters << "\t" << (std::clock() - optimizerStartTime) / CLOCKS_PER_SEC << endl;
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
    try
    {
        ini = INIReader(INI_FILE);
        string modelPath = BASE_DIR + ini.Get("PATH", "MODEL", "");
        resultDir = string(BASE_DIR + ini.Get("PATH", "RESULT_DIR", ""));

        initialTime = ini.GetReal("SIMULATION", "START_TIME", 0);
        finalTime = ini.GetReal("SIMULATION", "END_TIME", 1);
        penaltyWeight = ini.GetReal("OPTIMIZATION", "PENALTY_WEIGHT", 0.001);

        optLog = ofstream(resultDir + "/optLog.txt", ofstream::out);
        bestSoFarLog = ofstream(resultDir + "/bestSoFarLog.txt", ofstream::out);
        verboseLog = ofstream(resultDir + "/verboseLog.txt", ofstream::out);
        debugLog = ofstream(resultDir + "/debugLog.txt", ofstream::out);

        std::clock_t startTime = std::clock();

        // Ensures all components are printed out to the model file.
        Object::setSerializeAllDefaults(true);

        // Create a new OpenSim model from file
        Model osimModel(modelPath);

        forceAnalysis = new ForceReporter(&osimModel);
        osimModel.addAnalysis(forceAnalysis);

        bodyKinematics = new BodyKinematics(&osimModel);
        osimModel.addAnalysis(bodyKinematics);

        // The number of parameters is the number of actuators
        const Set<Actuator> &actuatorSet = osimModel.getActuators();
        int numActuators = actuatorSet.getSize();
        optLog << "numActuators = " << numActuators << endl;
        int numParameters = numActuators;

        /* Define initial values for controllerParameters. Each controller has two parameters, an
           initial time, and a duration for bang-bang control. There are as many controls as
           there are actuators because we assume symmetry (but each controller has two parameters).
           controllerParameters = [ti_0 duration_0 ti_1 duration_1 ... ]' */

           // initialize controller parameters
        Vector controllerParameters(numParameters, 0.05);

        // initialize initial time for each excitation
        controllerParameters[0] = 0.107829; //hamstrings
        controllerParameters[2] = 0.000694389; //bifmesh
        controllerParameters[4] = 0.296862; //glut_max
        controllerParameters[6] = 0.0533367; //iliopsoas
        controllerParameters[8] = 0.166411; //rect_fem
        controllerParameters[10] = 0.277515; //vasti
        controllerParameters[12] = 0.34415; //gastroc
        controllerParameters[14] = 0.315352; //soleus
        controllerParameters[16] = 0.0857547; //tib_ant
        controllerParameters[18] = 0.149619; //ercspn
        controllerParameters[20] = 0.468751; //intobl
        controllerParameters[22] = 0.455196; //extobl

        // initialize durations
        controllerParameters[1] = 0.429037; //hamstrings
        controllerParameters[3] = 0.12063; //bifemsh
        controllerParameters[5] = 0.392451; //glut_max
        controllerParameters[7] = 0.0303019; //iliopsoas
        controllerParameters[9] = 0.370407; //rect_fem
        controllerParameters[11] = 0.3; //vasti
        controllerParameters[13] = 0.343607; //gastroc
        controllerParameters[15] = 0.350789; //soleus
        controllerParameters[17] = 0.0276857; //tib_ant
        controllerParameters[19] = 0.243365; //ercspn
        controllerParameters[21] = 0.160011; //intobl
        controllerParameters[23] = 0.151908; //extobl

        // Add prescribed controllers to each muscle. Need to only loop numActuators/2 times since we are enforcing symmetry.
        // It is assumed that all actuators are listed for one side of the model first, then the other side, in the same order.
        for (int i = 0; i < numActuators / 2; i++)
        {
            // make a piecewise constant function for both sides
            PiecewiseConstantFunction bangBangControl;
            bangBangControl.addPoint(initialTime, 0);
            bangBangControl.addPoint(controllerParameters[2 * i], 1);
            bangBangControl.addPoint(controllerParameters[2 * i] + controllerParameters[2 * i + 1], 0);

            // add controller to right side
            PrescribedController *actuatorController_r = new PrescribedController();
            Set<Actuator> actuator_r;
            actuator_r.insert(0, osimModel.getActuators().get(i));
            actuatorController_r->setName(actuatorSet.get(i).getName());
            actuatorController_r->setActuators(*actuator_r.clone());
            actuatorController_r->prescribeControlForActuator(osimModel.getActuators().get(i).getName(), bangBangControl.clone());
            osimModel.addController(actuatorController_r);

            // add controller to left side
            PrescribedController *actuatorController_l = new PrescribedController();
            Set<Actuator> actuator_l;
            actuator_l.insert(0, osimModel.getActuators().get(i + numActuators / 2));
            actuatorController_l->setName(actuatorSet.get(i + numActuators / 2).getName());
            actuatorController_l->setActuators(*actuator_l.clone());
            actuatorController_l->prescribeControlForActuator(osimModel.getActuators().get(i + numActuators / 2).getName(), bangBangControl.clone());
            osimModel.addController(actuatorController_l);
        }

        // Create the OptimizationSystem. Initialize the objective function value "f".
        JumpingOptimizationSystem sys(numParameters, osimModel);
        Real f = NaN;

        // Set lower and upper bounds.
        Vector lower_bounds(numParameters, initialTime);
        Vector upper_bounds(numParameters, finalTime);

        // Limit the duration of the "bang" to be at least a certain value to decrease search space.
        for (int i = 1; i < numParameters; i += 2)
        {
            lower_bounds[i] = 0.0001;
        }

        sys.setParameterLimits(lower_bounds, upper_bounds);

        // Create an optimizer. Pass in our OptimizerSystem
        // and the name of the optimization algorithm.
        optLog << "using LBFGSB" << endl; Optimizer opt(sys, SimTK::LBFGSB); //LBFGSB was found to be beter for this problem
        //optLog << "using IPOPT" << endl; Optimizer opt(sys, InteriorPoint);

        // Specify settings for the optimizer
        opt.setConvergenceTolerance(0.01);
        opt.useNumericalGradient(true);
        opt.setMaxIterations(1000);
        opt.setLimitedMemoryHistory(500);
        //opt.setDiagnosticsLevel(4); // Lowest level that gives outer loop information for IPOPT

        // Optimize it!
        optLog << "Optimizing!" << endl;
        optLog << "Initial controllerParameters: " << controllerParameters << endl;
        optLog << "lower_bounds: " << lower_bounds << endl;
        optLog << "upper_bounds: " << upper_bounds << endl;
        f = opt.optimize(controllerParameters);

        optLog << "Elapsed time = " << (std::clock() - startTime) / CLOCKS_PER_SEC << "s" << endl;

        optLog << "OpenSim example completed successfully. Look in bestSoFarLog.txt for best solution.\n";

        optLog.close();
        bestSoFarLog.close();
        verboseLog.close();
        debugLog.close();
    }
    catch (const std::exception& ex)
    {
        std::cout << ex.what() << std::endl;
        //return 1;
    }

    // End of main() routine.
    system("pause");
    return 0;
}
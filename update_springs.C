#include "update_springs.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/IBSpringForceSpec.h>
#include <ibamr/IBBeamForceSpec.h>
#include <ibtk/LData.h>
#include <stdio.h>
#include <stdlib.h>

void
update_springs(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt, 
    ParameterFile & pf)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    //static const double pi = 4*atan(1);
    static const double L1 = 1; // length of computational domain (meters)
    static const int N1 = 512; // number of cartesian grid meshwidths at the finest level of the AMR grid
    static const int Nend = 10;
    static const double Let = pf.Let;
    static const double diameter = pf.tdiameter;
    static const double R2 = pf.tR2;
    static const double R1 = R2+diameter;
    static const double freq = pf.freq; // frequency of heart beats
    static const double period = 1/freq;  // Period of single heart beat
    static const double s_ramp = 0.02;    // time it takes to pinch or unpinch tube
    static const double speed = 2.0;  // Speed of contraction wave
    static const double pinch_time = Let/speed;  // time it takes for the pinch to travel

    static const double kappa1 = 3000.0;
    static const double kappa2 = 100*kappa1;
    static const double rest2 = 0.0014;
    static const double rest1 = rest2/10;
    //static const double bend1 = 5.0e8;
    //static const double bend2 = 5.0e6;
    static const int numPts = 860;
    static const double cutoffhigh = -R2-0.005*diameter;
    static const double cutofflow = -R1+0.005*diameter;

    // Normalize the current time of the simulation to the period of the heart beat
    double num = current_time/period; // Normalizes current time to the period of heart beat
    int intpart = (int)num;   // Extracts the integer part 
    double loop_time = num-intpart; // Subtracts the integer part (we are only interested in the decimal part, the unfinished portion of the beat, not the number of beats that have been completed)


    // Normalize the different periods to the period of a heart beat.
    double rampup = s_ramp/period; 
    double pinchend = (pinch_time+s_ramp)/period;
    double rampdown = (pinch_time+2*s_ramp)/period;

    // Find out the Lagrangian index ranges.
    const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

    // Update the spring lengths in their associated spring specs.
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
      {
        LNode* node_idx = *it;
       	IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
	IBSpringForceSpec* spring_spec = node_idx->getNodeDataItem<IBSpringForceSpec>();
	//	IBBeamForceSpec* beam_spec = node_idx->getNodeDataItem<IBBeamForceSpec>();
        if (spring_spec == NULL) continue;  // skip to next node

        // Here we update the resting length of the spring
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        resting_length    	is the resting length of the current Lagrangian point that is the "master index'
	//		  						Since there may be more than one spring associated with each point, it's a vector.
        //        resting_length[0]  	is the resting length of the first spring
        //        resting_length[1]  	would be the resting length of the second spring associated with that Lagrangian point.
        //
        // In this example, the resting length is increased by 0.01*dt each time step.

        const int lag_idx = node_idx->getLagrangianIndex();
	Point& X_target = force_spec->getTargetPointPosition();
	double& X_spring = force_spec->getStiffness();
	//std::vector<double>& resting_length = spring_spec->getRestingLengths();
	//The above is the old way of doing this.
        double resting_length = spring_spec->getParameters()[0][1];
	//	double& beam_stiff = beam_spec->getBendingRigidities()[0][0];
	//Note that you can also getStiffnesses
	double spring_stiffness = spring_spec->getParameters()[0][0];
	//resting_length+=0.01*dt;
	if (loop_time>=0.05*rampup && loop_time<=period){
	  if (lag_idx>=(lag_idxs.first+Nend) && lag_idx<(numPts/2-Nend) && X_target[1]<cutoffhigh) {
	    //resting_length = rest2;
	    X_spring = kappa2*10;
	    spring_stiffness = kappa2*10;
	    //  beam_stiff = bend2;
	  } else if (lag_idx<=(lag_idxs.second-Nend) && lag_idx>=(numPts/2+Nend) && X_target[1]>cutofflow){
	    //resting_length = rest2;
	    X_spring = kappa2*10;
	    spring_stiffness = kappa2*10;
	    //beam_stiff = bend2;
	  } else if (lag_idx>=(lag_idxs.first+Nend) && lag_idx<=(lag_idxs.second-Nend)) {
	    //resting_length = rest1;
	    X_spring = 100*kappa1;
	    spring_stiffness = 100*kappa1;
	    //beam_stiff = bend1;
	  } else {}
	} else {}
      }
    return;
}// update_springs

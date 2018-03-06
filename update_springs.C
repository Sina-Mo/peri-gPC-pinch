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

	// Initialize variables (some are pulled from the parameter file, those with "pf." in front of them)
    static const int numPts = pf.numPts;			// Number of target points on the heart tube
    static const int Nend = pf.Nendspring;			// Number of points to ignore on each of the four ends of the heart tube
    static const double diameter = pf.tdiameter;	// Diameter of the heart tube (m)
	static const double R2 = pf.tR2;				// Starting position of the top of the tube
    static const double R1 = R2+diameter;			// Starting position of the bottom of the tube
	static const double s_ramp = 0.02;   			// time it takes to pinch or unpinch tube
    static const double kappa1 = pf.kappa1;  		// Value that works for single-pinch: 3000.0
    static const double kappa2 = 750*kappa1; 		// Value that works for single-pinch: 100*kappa1
    //static const double rest2 = 0.0014;			// Value for resting length 1
    //static const double rest1 = rest2/10;			// Value for resting length 2
    static const double cutoffhigh = -R2-0.005*diameter;  // Cutoff value on y-axis for top of tube
    static const double cutofflow = -R1+0.005*diameter;	  // Cutoff value on y-axis for bottom of tube

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
        
    	const int lag_idx = node_idx->getLagrangianIndex(); 		// Get info on target points
		Point& X_target = force_spec->getTargetPointPosition();		// Get target point positions
		double& X_spring = force_spec->getStiffness();				// Get stiffness of target point tethering springs
       	//double resting_length = spring_spec->getParameters()[0][1]; // Get resting lengths of springs between vertex points
		//double& beam_stiff = beam_spec->getBendingRigidities()[0][0];	// Get stiffness of beams, not working right now
		double spring_stiffness = spring_spec->getParameters()[0][0]; // Get stiffness of springs between vertex points

		if (current_time>=0.05*s_ramp){  // Make sure time is greater than some portion of the ramp up time
	  
		  	if (lag_idx>=(lag_idxs.first+Nend) && lag_idx<(numPts/2-Nend) && X_target[1]<cutoffhigh) { // Top target points,
	    		// If the point y value on the top tube is less than the cutoff, make the springs stiffer:
		    	//resting_length = rest2;		// Change resting length of springs
	    		X_spring = kappa2*50;			// Change spring constant of target point tethering springs
		    	spring_stiffness = kappa2*10;	// Change spring constant of springs between vertex points
		    
		  	} else if (lag_idx<=(lag_idxs.second-Nend) && lag_idx>=(numPts/2+Nend) && X_target[1]>cutofflow){ // Bottom target points
	    		// If the point y value on the bottom tube is more than the cutoff, make the springs stiffer:
			    //resting_length = rest2;		// Change resting length of springs
	    		X_spring = kappa2*50;			// Change spring constant of target point tethering springs
		    	spring_stiffness = kappa2*10;	// Change spring constant of springs between vertex points
	    
			} else if (lag_idx>=(lag_idxs.first+Nend) && lag_idx<=(lag_idxs.second-Nend)) { // If not top or bottom but still in the tube    
				// If the point is on the tube but not in the stiff cutoff range, make the springs very flexible:
			    //resting_length = rest1;		// Change resting length of springs
				X_spring = 10*kappa1;			// Change spring constant of target point tethering springs
		    	spring_stiffness = 100*kappa1;	// Change spring constant of springs between vertex points
	    
			} else {} 	 // If not in tube, do nothing.
		} else {}		 // If before percent of ramp up time, do nothing.
	} // End main loop over points
    return;
}// update_springs

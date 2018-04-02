#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibtk/LData.h>
#include <stdio.h>
#include <stdlib.h>

void update_target_point_positions(
				   tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
				   LDataManager* const l_data_manager,
				   const double current_time,
				   const double dt,
				   ParameterFile & pf)
{
  // Find finest grid level number in simulation
  const int finest_ln = hierarchy->getFinestLevelNumber();

	// Initialize variables (some are pulled from the parameter file, those with "pf." in front of them)
	static const double freq = pf.freq;			// frequency of heart beats, decoupled from contraction speed
	static const double s_ramp = 0.02;			// time it takes to initially pinch tube
	static const double Let = pf.Let;			// Length of heart tube
	static const double speed = pf.speed;		// Speed of contraction wave
	static const int numPts = pf.numPts;		// number of points in geometry
	static const int Nend = pf.Nendtarget;		// Number of target points on the end to ignore
	static const double centery = 0.0;			// Center of domain y
	static const double diameter = pf.tdiameter;// Diameter of the tube
	static const double R2 = pf.tR2;			// distance from middle of domain to inner wall
	static const double R1 = R2+diameter;		// distance from middle of domain to outer wall
	static const double pamp = pf.pamp;			// percent occlusion of tube
	static const double sigma = pf.sigma;		// pointiness of pinch
	static const double mu = -0.5*Let;			// how far from x-center pinch starts
	static const double offset = speed/freq;	//	Offset distance between pinches
	static const double stop_pt = mu+3*offset;	//	Stopping distance
	double crampup = 0.0;						// Initial parameter controlling how far to pinch	
	double mu1 = mu; 							// Initial location of pinch #1
	double mu2 = mu+offset;						// Initial location of pinch #2
	double mu3 = mu+2*offset;					// Initial location of pinch #3
	double num;									// Needed for looping of pinches
	int intpart;								// Needed for looping of pinches
	double addto;								// Needed for looping of pinches
  
  	// Define the way the pinches should move by altering the values mu1, mu2, and mu3.
	if(current_time <= s_ramp)  // If the time is less than the time it takes to start the pinch
	 {  // Start the pinch by ramping it up

		 crampup = current_time/s_ramp;  // Ramp up the pinch based on current time

	 } else {  	   // Make the pinch travel down the tube 	
	 
	   	crampup = 1.0;  // Set the ramp up to 100% 
	   	mu1 += (current_time-s_ramp)*speed;	// Move pinch #1 forward by one step 
	   	mu2 += (current_time-s_ramp)*speed;	// Move pinch #2 forward by one step 
	   	mu3 += (current_time-s_ramp)*speed;  // Move pinch #3 forward by one step 
	   
   	// These loops control what happens when each pinch reaches the stopping point on the tube, 
   	// as soon as the position exceeds the stop point, it forces the pinch to start again from mu. 
   	// It then calculates the position as a fraction of where between mu and the stop point the pinch is, 
   	// independent of how many times it has cycled through the heart tube. 
	   if (mu1 >= stop_pt) { 	// Pinch #1
	     num = (mu1-mu)/(stop_pt-mu);
	   		intpart = (int)num;
	  	 	addto = num-intpart;
	   		mu1 = mu + (addto*3*offset); } 
	   else { }
	   if (mu2 >= stop_pt) { 	// Pinch #2
	     num = (mu2-mu)/(stop_pt-mu);
	   		intpart = (int)num;
	   		addto = num-intpart;
	   		mu2 = mu + (addto*3*offset); } 
	   else { }
	   if (mu3 >= stop_pt) { 	// Pinch #3
	     num = (mu3-mu)/(stop_pt-mu);
	   		intpart = (int)num;
		   	addto = num-intpart;
		   	mu3 = mu + (addto*3*offset); } 
	   else { }
	} // Ends loop over time
	
  // Find out the Lagrangian index ranges.
    const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

  // Get LMesh associated with finest level of patch.   
   Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
   vector<LNode*> nodes;
   nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
   nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

  //
  // Update target point positions

   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);

  // Loop over target point nodes
  for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        LNode* node_idx = *it;
        IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
	if (force_spec == NULL) continue;  // skip to next node

    const int lag_idx = node_idx->getLagrangianIndex();
    Point& X_target = force_spec->getTargetPointPosition();

	// Tell the target points how to move: 
	if (lag_idx >=(lag_idxs.first+Nend) && lag_idx<(numPts/2-Nend)) // If it is on the top (inner) part of the heart tube
    	{
            // Create the pinches at the new positions mu1, mu2, and mu3.                   
	  	 	X_target[1] = centery-R2 - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu1)/sigma,2.0)) - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu2)/sigma,2.0)) - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu3)/sigma,2.0));
	       
	    } else if (lag_idx >=(numPts/2+Nend) && lag_idx <(lag_idxs.second-Nend))  // If it is on the bottom (outer) part of the heart tube 
	    { 
			// Create the pinches at the new positions mu1, mu2, and mu3.
	    	X_target[1] = centery-R1 + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu1)/sigma,2.0)) + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu2)/sigma,2.0)) + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu3)/sigma,2.0));
	     
	    } else { } // Otherwise do nothing 

       
	} // Ends loop over target points

	return;

} //update_target_point_positions

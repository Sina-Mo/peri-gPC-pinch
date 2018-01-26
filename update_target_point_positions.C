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

  // Initialize a bunch of parameters
  static const double freq = pf.freq; // frequency of heart beats
  static const double Let = pf.Let;  // Length of heart tube
  static const double speed = 2.0;  // Speed of contraction wave
  static const int numPts = 860;  // number of points in geometry
  static const double s_ramp = 0.02;    // time it takes to pinch or unpinch tube
  static const double pinch_time = Let/speed;  // time it takes for the pinch to travel
  static const double centery = 0.0;  // Center of domain y
  static const double diameter = pf.tdiameter;  // Diameter of the tube
  static const double R2 = pf.tR2; // distance from middle of domain to inner wall
  static const double R1 = R2+diameter; // distance from middle of domain to outer wall
  static const double pamp = pf.pamp;  // percent occlusion of tube
  static const double sigma = 0.007;  // pointiness of pinch
  static const double mu = 0.18;   // how far from x-center pinch starts
  double mu1 = 0.18; 
  static const double period = 1/freq;  // Period of single heart beat
  // static const double pi = 4*atan(1);  // The number Pi  

  // Normalize the current time of the simulation to the period of the heart beat
  double num = current_time/period; // Normalizes current time to the period of heart beat
  int intpart = (int)num;   // Extracts the integer part 
  double loop_time = num-intpart; // Subtracts the integer part (we are only interested in the decimal part, the unfinished portion of the beat, not the number of beats that have been completed)

  // Normalize the different periods to the period of a heart beat.
  double rampup = s_ramp/period; 
  double pinchend = (pinch_time+s_ramp)/period;
  double rampdown = (pinch_time+2*s_ramp)/period;

  // Find out the Lagrangian index ranges.
  //
    const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

  //
  // Get LMesh associated with finest level of patch. 
  
   Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
   vector<LNode*> nodes;
   nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
   nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

  //
  // Update target point positions
  //
  
   tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);

  // Loop over target point nodes
  for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        LNode* node_idx = *it;
        IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec == NULL) continue;  // skip to next node

      // Update the target point positions here
      // NOETS: lag_idx is the index of the Lagrangian point (lag_idx = 0, 1, ... N-1 where N is the number of Lagrangian points. 
      //  X_target     is the target position fo the target point
      //  X_target[0]  is the x-component of the target position
      //  X_target[1]  is the y-component of the target position
      //  X_target[2]  is the z-component of the target position (for 3D simulations)

      const int lag_idx = node_idx->getLagrangianIndex();

      // Depending on the version of IBAMR, you'll need to select one fo the ways of accessing target point positions
      // FOR KD MODULE:
       Point& X_target = force_spec->getTargetPointPosition();

      //FOR NEMOS / KD (NOT MODULE):
      //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();

      //OLD:
      //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();

      // Tell the target points how to move: 
       if(loop_time<=rampup)  // If the normalized time is less than the time it takes to start the pinch
	 {  // Start the pinch

	   if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2)) // If it is on the top (inner) part of the heart tube
             {
               // Move points into a pinch, ramp this up based on the loop time                      
	       X_target[1] = centery-R2-(loop_time/rampup)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));

	     } else if (lag_idx >=(numPts/2) && lag_idx <numPts)  // If it is on the bottom (outer) part of the heart tube 
	     { 

	     X_target[1] = centery-R1+(loop_time/rampup)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));

	     } else { } // Otherwise do nothing 

	 } else if(loop_time>rampup && loop_time<=pinchend) // If the normalized time is greater than starting the pinch, less than the end of the pinch traveling along the tube,
	 { 
	   // Make the pinch travel down the tube 
	   double mu1 = mu-2*mu*((loop_time-rampup)/(pinchend-rampup));  // Define the offset of the pinch center based on the normalized time
	   
	   if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2)) // If it is on the top part of the tube
	     { 

	       X_target[1] = centery-R2-(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0)); 

	     } else if (lag_idx >=(numPts/2) && lag_idx < numPts) {  // If it is on the bottom part of the tube

	     X_target[1] = centery-R1+(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0));

	   } else { }  // Otherwise do nothing

	 } else if (loop_time>pinchend && loop_time<=rampdown) //If the normalized time is greater than the end of travel but less than the end of the unpinch 
	   { 
	     // Unpinch the tube
	     if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2))  // If it is on the top part of the tube
	       {
		
		 X_target[1] = centery-R2-(1-(loop_time-pinchend)/rampup)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu)/sigma,2.0)); 
	       
	       } else if (lag_idx >=(numPts/2) && lag_idx <numPts) // If it is on the bottom part of the tube 
	       { 
	       
	         X_target[1] = centery-R1+(1-(loop_time-pinchend)/rampup)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu)/sigma,2.0)); 
	     
	       } else { }  // Otherwise do nothing
	   
	   } else {  // Remaining time steps are a break with no target points being updated
      
       } //Ends time loop

       
  } // Ends loop over target points

  return;

} //update_target_point_positions

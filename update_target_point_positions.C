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
  static const double period = 1/freq;  // Period of single heart beat   
  static const double s_ramp = 0.02;    // time it takes to pinch or unpinch tube
  
  static const double Let = pf.Let;  // Length of heart tube
  static const double speed = 0.3;  // Speed of contraction wave
  static const int numPts = 860;  // number of points in geometry
  static const int Nend = 54;
  static const double acc_ramp = 0.01;  // time to accelerate the pinch movement
  static const double centery = 0.0;  // Center of domain y
  static const double diameter = pf.tdiameter;  // Diameter of the tube
  static const double R2 = pf.tR2; // distance from middle of domain to inner wall
  static const double R1 = R2+diameter; // distance from middle of domain to outer wall
  static const double pamp = pf.pamp;  // percent occlusion of tube
  static const double sigma = 0.0085;  // pointiness of pinch
  static const double mu = -0.5*Let;   // how far from x-center pinch starts
  static const double offset = speed/freq; //Offset bewteen peaks
  double dx = speed*dt;	//Distance to move each time step
  static const double stop_pt = mu+3*offset; //Stopping distance
  double crampup = 0.0;
  double mu1 = mu; 
  double mu2 = mu+offset;
  double mu3 = mu+2*offset;
  double num;
  int intpart;
  double addto;
  
  if(current_time <= s_ramp)  // If the normalized time is less than the time it takes to start the pinch
	 {  // Start the pinch
		 crampup = current_time/s_ramp;
	 
	// } else if (current_time > s_ramp && current_time <= (s_ramp+acc_ramp)) {
	 
	//	 camp = 1.0;
		
	//	 mu1 += (current_time/(s_ramp+acc_ramp))*dx;
	//	 mu2 += (current_time/(s_ramp+acc_ramp))*dx;
	//	 mu3 += (current_time/(s_ramp+acc_ramp))*dx;

       //if (mu1 >= stop_pt) { mu1 = mu; } else { }
	   //if (mu2 >= stop_pt) { mu2 = mu; } else { }
	   //if (mu3 >= stop_pt) { mu3 = mu; } else { }
		 
	 } else {
	 
	   crampup = 1.0;
	   mu1 += (current_time-s_ramp)*speed;
	   mu2 += (current_time-s_ramp)*speed;
	   mu3 += (current_time-s_ramp)*speed;
	   
	   // Make the pinch travel down the tube 	   
	   if (mu1 >= stop_pt) { 
	   	num = mu1/stop_pt;
	   	intpart = (int)num;
	   	addto = num-intpart;
	   	mu1 = mu + (addto*stop_pt); } 
	   	else { }
	   if (mu2 >= stop_pt) { 
	   	num = mu2/stop_pt;
	   	intpart = (int)num;
	   	addto = num-intpart;
	   	mu2 = mu + (addto*stop_pt); } 
	   	else { }
	   if (mu3 >= stop_pt) { 
	   	num = mu3/stop_pt;
	   	intpart = (int)num;
	   	addto = num-intpart;
	   	mu3 = mu + (addto*stop_pt); } 
	   	else { }
	} // Ends time loop
	
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
       //double X_stiff = force_spec->getStiffness();
      //FOR NEMOS / KD (NOT MODULE):
      //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();

      //OLD:
      //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();

      // Tell the target points how to move: 

	if (lag_idx >=(lag_idxs.first+Nend) && lag_idx<(numPts/2-Nend)) // If it is on the top (inner) part of the heart tube
    	{
               // Move points into a pinch, ramp this up based on the loop time                      
	  	 	X_target[1] = centery-R2 - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu1)/sigma,2.0)) - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu2)/sigma,2.0)) - crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu3)/sigma,2.0));
	       
	    } else if (lag_idx >=(numPts/2+Nend) && lag_idx <(lag_idxs.second-Nend))  // If it is on the bottom (outer) part of the heart tube 
	    { 

	    	X_target[1] = centery-R1 + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu1)/sigma,2.0)) + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu2)/sigma,2.0)) + crampup*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu3)/sigma,2.0));
	     
	    } else { } // Otherwise do nothing 

       
	  } // Ends loop over target points

	return;

} //update_target_point_positions

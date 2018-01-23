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
  const int finest_ln = hierarchy->getFinestLevelNumber();

  static const int numPts = 860;  // number of points in geometry
  static double s_ramp = 0.02;    // time it takes to pinch or unpinch tube
  static double pinch_time = 0.20;  // time it takes for the pinch to travel
  static double centery = 0.0;  // Center of domain y
  static double diameter = pf.tdiameter;  // Diameter of the tube
  static double R2 = pf.tR2; // distance from middle of domain to inner wall
  static double R1 = R2+diameter; // distance from middle of domain to outer wall
  static double pamp = pf.pamp;  // percent occlusion of tube
  static double sigma = 0.007;  // pointiness of pinch
  static double mu = 0.18;   // how far from x-center pinch starts
  double mu1 = 0.18; 

  // static const double pi = 4*atan(1);  // The number Pi

  //
  // Find out the Lagrangian index ranges.
  //
    const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

  //
  // Get LMesh associated with finest level of patch. 
  
   Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
   vector<LNode*> nodes;
   nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
   nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());

  // Update target point positions
  //
  tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);

  // Loop over nodes
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
       if(current_time<=s_ramp)
	 {
	   // X_target[0] = xP1[lag_idx] + (xP2[lag_idx]-xP1[lag_idx])*(current_time/s_ramp);
	   // X_target[1] = yP1[lag_idx] + (yP2[lag_idx]-yP1[lag_idx])*(current_time/s_ramp);  //This works.

	   if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2))
             {
	       //X_target[0] = xP2[lag_idx] + (xP3[lag_idx]-xP2[lag_idx])*(current_time/s_ramp);                                                              
	       X_target[1] = centery-R2-(current_time/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));
	       //X_target[1]=X_target[1]+0.00001;
	     } else if (lag_idx >=(numPts/2) && lag_idx <numPts) {
	     X_target[1] = centery-R1+(current_time/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));
	     //X_target[1] = X_target[1]-0.00001;
	   } else { }

	 } else if(current_time>s_ramp && current_time<=(pinch_time+s_ramp))
	 {
	   double mu1 = mu-2*mu*((current_time-s_ramp)/(pinch_time));
	   
	   if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2))
	     { 
	   //X_target[0] = xP2[lag_idx] + (xP3[lag_idx]-xP2[lag_idx])*(current_time/s_ramp);
	       X_target[1] = centery-R2-(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0));
	     } else if (lag_idx >=(numPts/2) && lag_idx < numPts) {
	     X_target[1] = centery-R1+(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0));
	   } else { }
	 } else if (current_time>(pinch_time+s_ramp) && current_time<=(2*s_ramp+pinch_time)) 
	   {
	   
	     if (lag_idx >=lag_idxs.first && lag_idx<(numPts/2))
	       {
		 X_target[1] = centery-R2-(1-(current_time-(pinch_time+s_ramp))/(s_ramp))*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu)/sigma,2.0)); 
	       } else if (lag_idx >=(numPts/2) && lag_idx <numPts) {
	       X_target[1] = centery-R1+(1-(current_time-(pinch_time+s_ramp))/(s_ramp))*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]-mu)/sigma,2.0)); 
	     } else { }
	   
	   } else {
	   }

       
  } // Ends loop over target points

  return;

} //update_target_point_positions

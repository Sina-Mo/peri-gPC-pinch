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
  static double s_ramp = 0.05;    // time it takes to pinch or unpinch tube
  static double pinch_time = 0.50;  // time it takes for the pinch to travel
  static double centery = 0.0;  // Center of domain y
  static double diameter = pf.tdiameter;  // Diameter of the tube
  static double R2 = pf.tR2; // distance from middle of domain to inner wall
  static double R1 = R2+diameter; // distance from middle of domain to outer wall
  static double pamp = pf.pamp;  // percent occlusion of tube
  static double sigma = 0.007;  // pointiness of pinch
  static double mu = 0.18;   // how far from x-center pinch starts
  double mu1 = 0.18; 

  // static const double pi = 4*atan(1);  // The number Pi

  //Storing the (x,y) positions in all three states
  static double xP1[numPts];  //Stores x-values for state 1
  static double yP1[numPts];  //Stores y-values for state 1
  static double xP2[numPts];  //Stores x-values for state 2
  static double yP2[numPts];  //Stores y-values for state 2
  static double xP3[numPts];  //Stores x-values for state 3
  static double yP3[numPts];  //Stores y-values for state 3

  FILE *fP1x, *fP1y, *fP2x, *fP2y, *fP3x, *fP3y;

  fP1x = fopen("xP1.txt", "r");
  fP1y = fopen("yP1.txt", "r");
  
  fP2x = fopen("xP2.txt", "r");
  fP2y = fopen("yP2.txt", "r");
  
  fP3x = fopen("xP3.txt", "r");
  fP3y = fopen("yP3.txt", "r");

  for(int i=0; i<numPts; i++){

    fscanf(fP1x,"%lf\n",&xP1[i]);
    fscanf(fP1y,"%lf\n",&yP1[i]);

    fscanf(fP2x,"%lf\n",&xP2[i]);
    fscanf(fP2y,"%lf\n",&yP2[i]);

    fscanf(fP3x,"%lf\n",&xP3[i]);
    fscanf(fP3y,"%lf\n",&yP3[i]);

  }

  fclose(fP1x);
  fclose(fP1y);
  fclose(fP2x);
  fclose(fP2y);
  fclose(fP3x);
  fclose(fP3y);
  
  
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

	   if (lag_idx <(numPts/2))
             {
	       //X_target[0] = xP2[lag_idx] + (xP3[lag_idx]-xP2[lag_idx])*(current_time/s_ramp);                                                              
	       X_target[1] = centery-R2-(current_time/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));
	       //X_target[1]=X_target[1]+0.00001;
	     } else {
	     X_target[1] = centery-R1+(current_time/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0));
	     //X_target[1] = X_target[1]-0.00001;
	   }

	 } else if(current_time>s_ramp && current_time<=(pinch_time+s_ramp))
	 {
	   double mu1 = mu-2*mu*((current_time-s_ramp)/(pinch_time+s_ramp));
	   
	   if (lag_idx <(numPts/2))
	     { 
	   //X_target[0] = xP2[lag_idx] + (xP3[lag_idx]-xP2[lag_idx])*(current_time/s_ramp);
	       X_target[1] = centery-R2-(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0));
	     } else {
	     X_target[1] = centery-R1+(diameter*pamp/2)*exp(-0.5*pow((X_target[0]+mu1)/sigma,2.0));
	   }
	 } else if (current_time>(pinch_time+s_ramp) && current_time<=(2*s_ramp+pinch_time)) 
	   {
	   
	     if (lag_idx <(numPts/2))
	       {
		 X_target[1] = centery-R2-(1-(current_time-(pinch_time+s_ramp))/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0)); 
	       } else {
	         X_target[1] = centery-R1+(1-(current_time-(pinch_time+s_ramp))/s_ramp)*(diameter*pamp/2.0)*exp(-0.5*pow((X_target[0]+mu)/sigma,2.0)); 
	     }
	   
	   } else {
	   }

       
  } // Ends loop over target points

  return;

} //update_target_point_positions

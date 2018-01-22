#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <stdio.h>
#include <stdlib.h>

void update_target_point_positions(
    tbox::Pointer<hier::PatchHierarchy<NDIM> > hierarchy,
    LDataManager* const l_data_manager,
    const double current_time,
    const double dt)
{
    const int finest_ln = hierarchy->getFinestLevelNumber();

    static const int numPts = 1287; //number of points in geometry

 	//
    // Storing all the (X,Y) positions for both PHASES
    //
	static double xP1[numPts]; //Stores X-values for Phase 1
    static double yP1[numPts]; //Stores Y-values for Phase 1
    static double xP2[numPts]; //Stores X-values for Phase 2
    static double yP2[numPts]; //Stores Y-values for Phase 2
	
	FILE *fP1x, *fP1y, *fP2x, *fP2y;
    
	fP1x = fopen("xP1.txt", "r");
    fP1y = fopen("yP1.txt", "r");

    fP2x = fopen("xP2.txt", "r");
    fP2y = fopen("yP2.txt", "r");

	for(int i=0; i<numPts; i++){
		
	  fscanf(fP1x,"%lf\n",&xP1[i]);
	  fscanf(fP1y,"%lf\n",&yP1[i]);
		
	  fscanf(fP2x,"%lf\n",&xP2[i]);
	  fscanf(fP2y,"%lf\n",&yP2[i]);
		
	  //cout << i << "\n";
	}
    fclose(fP1x);
    fclose(fP1y);
    fclose(fP2x);
    fclose(fP2y);













    //
    // Find out the Lagrangian index ranges.
    //
    //const std::pair<int,int>& lag_idxs = l_data_manager->getLagrangianStructureIndexRange(0, finest_ln);

    //
    // Get the LMesh (which we assume to be associated with the finest level of
    // the patch hierarchy).  Note that we currently need to update both "local"
    // and "ghost" node data.
    
    Pointer<LMesh> mesh = l_data_manager->getLMesh(finest_ln);
    vector<LNode*> nodes;
    nodes.insert(nodes.end(), mesh->getLocalNodes().begin(), mesh->getLocalNodes().end());
    nodes.insert(nodes.end(), mesh->getGhostNodes().begin(), mesh->getGhostNodes().end());


    //
    // Update the target point positions in their associated target point force specs.
    //
    tbox::Pointer<hier::PatchLevel<NDIM> > level = hierarchy->getPatchLevel(finest_ln);
	
 
    // Loop over ALL nodes
    //
    for (vector<LNode*>::iterator it = nodes.begin(); it != nodes.end(); ++it)
    {
        LNode* node_idx = *it;
        IBTargetPointForceSpec* force_spec = node_idx->getNodeDataItem<IBTargetPointForceSpec>();
        if (force_spec == NULL) continue;  // skip to next node


        // Here we update the position of the target point.
        //
        // NOTES: lag_idx      is the "index" of the Lagrangian point (lag_idx = 0, 1, ..., N-1, where N is the number of Lagrangian points)
        //        X_target     is the target position of the target point
        //        X_target[0]  is the x component of the target position
        //        X_target[1]  is the y component of the target position
        //        X_target[2]  is the z component of the target position (for a 3D simulation)
        //



        // Get lag_idx associated with lag_pt (i.e, index in array)
        //
        const int lag_idx = node_idx->getLagrangianIndex();
		

        //
        //Depending on the version of IBAMR, you need to select one of the ways of accessing target point positions
        //
        //FOR KD MODULE:
        Point& X_target = force_spec->getTargetPointPosition();
        //
        //FOR NEMOs / KD (NOT MODULE)
        //TinyVector<double,NDIM>& X_target = force_spec->getTargetPointPosition();
        //
        //OLD:
        //IBTK::Vector<double,NDIM>& X_target = force_spec->getTargetPointPosition();



        //
        // BLOCK WHERE YOU TELL YOUR TARGET POINTS HOW TO MOVE!!!!
        //
        // Here we wish to interpolate between State 1 -> State 2
        //
        
    
        //Simple linear interpolation between states! (Can change this, i.e. tanh(t), poly, etc.)
        X_target[0] = xP1[lag_idx] + current_time*(xP2[lag_idx]-xP1[lag_idx]);
		X_target[1] = yP1[lag_idx] + current_time*(yP2[lag_idx]-yP1[lag_idx]);




		

	// Move x-Pts over by dt to the right each time-step	
        //
	// Write that HERE!		
		

		//cout << lag_idx << "\n"
    
    } //Ends for-loop over all lag_pts
    
    return;

}// update_target_point_positions

#include "update_springs.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibamr/IBSpringForceSpec.h>
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
    static const double Let = pf.Let;
    static const double diameter = pf.tdiameter;
    static const double R2 = pf.tR2;
    static const double R1 = R2+diameter;
    static const double kappa1 = 2.0;
    static const double kappa2 = 40.0;
    static const int numPts = 860;
    static const double cutoffhigh = R2-0.1*diameter;
    static const double cutofflow = R1+0.1*diameter;

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
	//std::vector<double>& resting_length = spring_spec->getRestingLengths();
	//The above is the old way of doing this.
	//double resting_length = spring_spec->getParameters()[0][1];

	//Note that you can also getStiffnesses
	double spring_stiffness = spring_spec->getParameters()[0][0];
	//resting_length+=0.01*dt;
	if (lag_idx<(numPts/2) && X_target[1]<cutoffhigh) {
	  spring_stiffness = kappa2;
	} else if (lag_idx>=(numPts/2) && X_target[1]>cutofflow){
	  spring_stiffness = kappa2;
	} else {
	  spring_stiffness = kappa1;
	}
	
      }
    return;
}// update_springs

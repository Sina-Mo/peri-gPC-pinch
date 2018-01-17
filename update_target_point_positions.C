#include "update_target_point_positions.h"
#include <ibamr/IBTargetPointForceSpec.h>
#include <ibtk/LData.h>
#include <stdio.h>
#include <stdlib.h>

void update_target_point_positions(
				   tbox::Pointer<hier::PatchHierarch<NDIM> > hierarchy, 
				   LDataManager* const l_data_manager, 
				   const double current_time, 
				   const double dt, 
				   ParameterFile & pf)
{
  const int finest_ln = hierarchy->getFinestLevelNumber();
 
  static const double pi = 4*atan(1);

} //update_target_point_positions

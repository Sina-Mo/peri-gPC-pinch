#ifndef included_update_target_point_positions
#define included_update_target_point_positions

#include <ibamr/app_namespaces.h>
#include <ibtk/LDataManager.h>
#include "parameterfile.h"

/*
 * Update the positions of the target point specifications.
 */
void
update_target_point_positions(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* l_data_manager,
    const double current_time,
    const double dt,
    ParameterFile &);

#endif //#ifndef included_update_target_point_positions


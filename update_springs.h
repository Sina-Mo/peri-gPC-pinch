#ifndef included_update_springs
#define included_update_springs

#include <ibamr/app_namespaces.h>
#include <ibtk/LDataManager.h>

/*
 * Update the positions of the target point specifications.
 */
void
update_springs(
    Pointer<PatchHierarchy<NDIM> > hierarchy,
    LDataManager* l_data_manager,
    const double current_time,
    const double dt,
    ParameterFile &);

#endif //#ifndef included_update_target_point_positions

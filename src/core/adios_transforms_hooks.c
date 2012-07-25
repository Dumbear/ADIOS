
#include <stdint.h>
#include <assert.h>

#include "adios_logger.h"
#include "adios_internals.h"
#include "adios_transforms_common.h"
#include "adios_transforms_hooks.h"
#include "adios_transforms_util.h"
#include "public/adios_selection.h"

/*
 * The transform method registry, containing all necessary function pointers
 */
adios_transform_method TRANSFORM_METHODS[num_adios_transform_types];

int adios_transforms_initialized = 0;

void adios_transform_read_init() {
    if (adios_transforms_initialized)
        return;

    ASSIGN_READ_FNS(none, adios_transform_none);
    ASSIGN_READ_FNS(identity, adios_transform_identity);
    ASSIGN_READ_FNS(alacrity, adios_transform_alacrity);

    adios_transforms_initialized = 1;
}

uint64_t adios_transform_calc_vars_transformed_size(enum ADIOS_TRANSFORM_TYPE transform_type, uint64_t orig_size, int num_vars) {
    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_METHODS[transform_type].transform_calc_vars_transformed_size(orig_size, num_vars);
}
enum ADIOS_ERRCODES adios_transform_retrieve_subvolume(
        enum ADIOS_TRANSFORM_TYPE transform_type,
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *subv_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {

    assert(transform_type >= adios_transform_none && transform_type < num_adios_transform_types);
    return TRANSFORM_METHODS[transform_type].transform_retrieve_subvolume(var, global_out, time_index, subv_spec, read_state, read_delegate, swap_endianness);
}

// Error function used by unimplemented transport methods
void unimplemented_transform_function(const char *tmethod, const char *func) {
    log_error("Transport method %s, function %s is not available " \
              "in this build of ADIOS\n", tmethod, func);
    assert(0); // To generate an error
}

// Add error message stubs for the "none" method (it should never be called)
DEFINE_FNS_UNIMPL(none);

// Implementation of the "identity" transform, which does nothing, but
// exercises the transform framework for testing.
// (adios_transform_identity_apply is in adios_transforms_hooks_write.c)

uint64_t adios_transform_identity_calc_vars_transformed_size(uint64_t orig_size, int num_vars) {
    return orig_size;
}

// TODO: implement, need more powerful helper functions
enum ADIOS_ERRCODES adios_transform_identity_retrieve_subvolume(
        struct adios_index_var_struct_v1 *var, void *global_out, int time_index,
        const adios_subvolume_copy_spec *copy_spec,
        void *read_state, adios_transform_var_read_delegate read_delegate,
        enum ADIOS_FLAG swap_endianness) {

    assert(var->characteristics[time_index].transform.transform_type == adios_transform_identity);

    // Read the variable in its entirety
    void *buf = read_delegate(var, time_index, 0, adios_transform_var_get_transformed_size(var, time_index), read_state);

    copy_subvolume_with_spec(global_out, buf, copy_spec, adios_transform_get_var_original_type(var), swap_endianness);

    return err_no_error;
}
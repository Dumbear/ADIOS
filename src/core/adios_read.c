/* 
 * ADIOS is freely available under the terms of the BSD license described
 * in the COPYING file in the top level directory of this source distribution.
 *
 * Copyright (c) 2008 - 2009.  UT-BATTELLE, LLC. All rights reserved.
 */

#include <stdlib.h>
#include <string.h>
#include "public/adios_read_v2.h"
#include "public/adios_error.h"
#include "core/adios_logger.h"
#include "core/common_read.h"
#define BYTE_ALIGN 8

typedef struct _mask_read_request {
    char *var_name;
    ADIOS_SELECTION *sel;
    int from_step;
    int n_steps;
    void *data;
    struct _mask_read_request *next;
} mask_read_request;

static int g_n_dims;
static int g_n_blocks;
static int g_mask_id;

static int **g_mask_bits = NULL;

static int *g_mask_bits_offset;
static int *g_mask_bits_length;
static int *g_mask_offset;
static int *g_mask_length;

static int **g_block_offset;
static int **g_block_length;

static mask_read_request *g_requests = NULL;

static void add_mask_read_request(mask_read_request *req) {
    req->next = g_requests;
    g_requests = req;
}

const char *adios_errmsg ()
{
    return adios_get_last_errmsg();
}

int adios_read_init_method (enum ADIOS_READ_METHOD method, MPI_Comm comm, const char * parameters)
{
    return common_read_init_method (method,comm,parameters);
} 

int adios_read_finalize_method(enum ADIOS_READ_METHOD method)
{
    int retval = common_read_finalize_method(method);
    log_debug ("adios_read_finalize_method completed\n");
    return retval;
}

ADIOS_FILE * adios_read_open (const char * fname,
                              enum ADIOS_READ_METHOD method,
                              MPI_Comm comm,
                              enum ADIOS_LOCKMODE lock_mode,
                              float timeout_sec)
{
    return common_read_open (fname, method, comm, lock_mode, timeout_sec);
}

ADIOS_FILE * adios_read_open_file (const char * fname,
                                   enum ADIOS_READ_METHOD method,
                                   MPI_Comm comm)
{
    return common_read_open_file (fname, method, comm);
}

int adios_read_close (ADIOS_FILE *fp) 
{
    return common_read_close (fp);
}

int adios_advance_step (ADIOS_FILE *fp, int last, float timeout_sec)
{
    return common_read_advance_step (fp, last, timeout_sec);
}

void adios_release_step (ADIOS_FILE *fp)
{
    common_read_release_step (fp);
}

ADIOS_VARINFO * adios_inq_var (ADIOS_FILE  *fp, const char * varname) 
{
    return common_read_inq_var (fp, varname);
}

ADIOS_VARINFO * adios_inq_var_byid (ADIOS_FILE  *fp, int varid)
{
    return common_read_inq_var_byid (fp, varid);
}

void adios_free_varinfo (ADIOS_VARINFO *vp)
{
    common_read_free_varinfo (vp);
}

int adios_inq_var_stat (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo,
                                    int per_step_stat, int per_writer_stat)
{
    return common_read_inq_var_stat (fp, varinfo, per_step_stat, per_writer_stat);
}

int adios_inq_var_blockinfo (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    return common_read_inq_var_blockinfo (fp, varinfo);
}

ADIOS_MESH * adios_inq_mesh_byid (ADIOS_FILE *fp, int meshid)
{
    return common_read_inq_mesh_byid (fp, meshid);
}

int adios_complete_meshinfo (ADIOS_FILE *datafile, ADIOS_FILE *meshfile, ADIOS_MESH *meshinfo)
{
    return common_read_complete_meshinfo (datafile, meshfile, meshinfo);
}

void adios_free_meshinfo (ADIOS_MESH *meshinfo)
{
    common_read_free_meshinfo (meshinfo);
}

int adios_inq_var_meshinfo (ADIOS_FILE *fp, ADIOS_VARINFO * varinfo)
{
    return common_read_inq_var_meshinfo (fp, varinfo);
}

static int bits_decompress(uint64_t mask_bits_length, int *mask_bits, uint64_t *mask_length, char **mask) {
    *mask_length = mask_bits_length * sizeof(int) * 8;
    *mask = malloc(*mask_length * sizeof(char) );
    if (*mask == NULL) {
        return 1;
    }

    uint64_t i, j, index = 0;
    for (i = 0; i < mask_bits_length; ++i) {
        for (j = 0; j < sizeof(int) * 8; ++j) {
            (*mask)[index] = (mask_bits[i] >> j) & 1;
            ++index;
        }
    }
    return err_no_error;
}

int adios_read_set_mask(ADIOS_FILE *fp, int mask_id, int n_dims, const char **offset_var_names, const char **length_var_names) {
    int i, j;
    char s[4][64];

    g_n_dims = n_dims;
    g_mask_id = mask_id;
    sprintf(s[0], "mask_vars/mask_%d/mask_bits", g_mask_id);
    ADIOS_VARINFO *mask_info = adios_inq_var(fp, s[0]);
    g_n_blocks = mask_info->nblocks[0];
    adios_free_varinfo(mask_info);

    sprintf(s[0], "mask_vars/mask_%d/mask_bits_offset", g_mask_id);
    sprintf(s[1], "mask_vars/mask_%d/mask_bits_length", g_mask_id);
    sprintf(s[2], "mask_vars/mask_%d/mask_offset", g_mask_id);
    sprintf(s[3], "mask_vars/mask_%d/final_length", g_mask_id);

    g_mask_bits_offset = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_bits_length = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_offset = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_length = (int *)malloc(sizeof(int) * g_n_blocks);

    g_block_offset = (int **)malloc(sizeof(int *) * g_n_blocks);
    g_block_length = (int **)malloc(sizeof(int *) * g_n_blocks);
    for (i = 0; i < g_n_blocks; ++i) {
        g_block_offset[i] = (int *)malloc(sizeof(int) * g_n_dims);
        g_block_length[i] = (int *)malloc(sizeof(int) * g_n_dims);
    }

    for (i = 0; i < g_n_blocks; ++i) {
        ADIOS_SELECTION *sel = adios_selection_writeblock(i);
        adios_schedule_read(fp, sel, s[0], 0, 1, &g_mask_bits_offset[i]);
        adios_schedule_read(fp, sel, s[1], 0, 1, &g_mask_bits_length[i]);
        adios_schedule_read(fp, sel, s[2], 0, 1, &g_mask_offset[i]);
        adios_schedule_read(fp, sel, s[3], 0, 1, &g_mask_length[i]);
        for (j = 0; j < g_n_dims; ++j) {
            adios_schedule_read(fp, sel, offset_var_names[j], 0, 1, &g_block_offset[i][j]);
            adios_schedule_read(fp, sel, length_var_names[j], 0, 1, &g_block_length[i][j]);
        }
        adios_selection_delete(sel);
    }
    adios_perform_reads(fp, 1);

    g_mask_bits = (int **)malloc(sizeof(int *) * g_n_blocks);
    for (i = 0; i < g_n_blocks; ++i) {
        g_mask_bits[i] = NULL;
    }

    for (i = 0; i < g_n_blocks; ++i) {
        printf("-------------------- Block %d --------------------\n", i);
        printf("Mask bits range: [%d, %d]\n", g_mask_bits_offset[i], g_mask_bits_offset[i] + g_mask_bits_length[i] - 1);
        printf("Mask range: [%d, %d]\n", g_mask_offset[i], g_mask_offset[i] + g_mask_length[i] - 1);
        printf("Original domain range: ");
        for (j = 0; j < g_n_dims; ++j) {
            if (j == 0) printf("[");
            printf("%d", g_block_offset[i][j]);
            printf(j + 1 == g_n_dims ? "]" : ", ");
        }
        printf(" - ");
        for (j = 0; j < g_n_dims; ++j) {
            if (j == 0) printf("[");
            printf("%d", g_block_offset[i][j] + g_block_length[i][j] - 1);
            printf(j + 1 == g_n_dims ? "]" : ", ");
        }
        printf("\n");
    }

    return err_no_error;
}

int adios_read_unset_mask(const ADIOS_FILE *fp) {
    int i, j;

    for (i = 0; i < g_n_blocks; ++i) {
        if (g_mask_bits[i] != NULL) {
            free(g_mask_bits[i]);
        }
    }
    free(g_mask_bits);
    g_mask_bits = NULL;

    free(g_mask_bits_offset);
    free(g_mask_bits_length);
    free(g_mask_offset);
    free(g_mask_length);
    g_mask_bits_offset = NULL;
    g_mask_bits_length = NULL;
    g_mask_offset = NULL;
    g_mask_length = NULL;

    for (i = 0; i < g_n_blocks; ++i) {
        free(g_block_offset[i]);
        free(g_block_length[i]);
    }
    free(g_block_offset);
    free(g_block_length);
    g_block_offset = NULL;
    g_block_length = NULL;

    /* TODO: free g_requests*/
    g_requests = NULL;

    return err_no_error;
}

int adios_boundingbox_touched(const ADIOS_SELECTION *sel, int *offset, int *length) {
    int i;
    for (i = 0; i < sel->u.bb.ndim; ++i) {
        if (sel->u.bb.start[i] >= offset[i] + length[i] || sel->u.bb.start[i] + sel->u.bb.count[i] <= offset[i]) {
            return 0;
        }
    }
    return 1;
}

int adios_schedule_read (const ADIOS_FILE      * fp,
                             const ADIOS_SELECTION * sel,
                             const char            * varname,
                             int                     from_steps,
                             int                     nsteps,
                             void                  * data)
{
    if (strncmp(varname, "mask_vars/", 10) == 0 || g_mask_bits == NULL) {
        return common_read_schedule_read (fp, sel, varname, from_steps, nsteps, NULL, data);
    }

    int i;

    char s[64];
    sprintf(s, "mask_vars/mask_%d/mask_bits", g_mask_id);
    for (i = 0; i < g_n_blocks; ++i) {
        if(g_mask_bits[i] == NULL && adios_boundingbox_touched(sel, g_block_offset[i], g_block_length[i])) {
            g_mask_bits[i] = (int *)malloc(sizeof(int) * g_mask_bits_length[i]);
            uint64_t offset = g_mask_bits_offset[i], length = g_mask_bits_length[i];
            ADIOS_SELECTION *mask_bits_sel = adios_selection_boundingbox(1, &offset, &length);
            adios_schedule_read(fp, mask_bits_sel, s, from_steps, nsteps, g_mask_bits[i]);

            mask_read_request *req = (mask_read_request *)malloc(sizeof(mask_read_request));
            req->var_name = strdup(varname);
            req->sel = copy_selection(sel);
            req->from_step = from_steps;
            req->n_steps = nsteps;
            req->data = data;
            add_mask_read_request(req);
        }
    }

    return err_no_error;
}

int adios_schedule_read_byid (const ADIOS_FILE      * fp,
                                  const ADIOS_SELECTION * sel,
                                  int                     varid,
                                  int                     from_steps,
                                  int                     nsteps,
                                  void                  * data)
{
    return common_read_schedule_read_byid (fp, sel, varid, from_steps, nsteps, NULL, data);
}

int adios_schedule_read_param (const ADIOS_FILE * fp,
                               const ADIOS_SELECTION * sel,
                               const char            * varname,
                               int                     from_steps,
                               int                     nsteps,
                               const char            * param,
                               void                  * data) {
    return common_read_schedule_read (fp, sel, varname, from_steps, nsteps, param, data);
}

int adios_schedule_read_byid_param (const ADIOS_FILE * fp,
                                    const ADIOS_SELECTION * sel,
                                    int                     varid,
                                    int                     from_steps,
                                    int                     nsteps,
                                    const char            * param,
                                    void                  * data) {
    return common_read_schedule_read_byid (fp, sel, varid, from_steps, nsteps, param, data);
}

static int adios_read_mask_var(mask_read_request *req) {
    printf("    Should read var: %s\n", req->var_name);
    return err_no_error;
}

int adios_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    if (g_mask_bits == NULL) {
        return common_read_perform_reads (fp, blocking);
    }

    common_read_perform_reads (fp, blocking);

    mask_read_request *req;
    for (req = g_requests; req != NULL; req = req->next) {
        adios_read_mask_var(req);
    }
}

int adios_check_reads (const ADIOS_FILE * fp, ADIOS_VARCHUNK ** chunk)
{
    return common_read_check_reads (fp, chunk);
}

void adios_free_chunk (ADIOS_VARCHUNK *chunk)
{
    common_read_free_chunk (chunk);
}

int adios_get_attr (ADIOS_FILE  * fp, const char * attrname, enum ADIOS_DATATYPES * type,
                    int * size, void ** data)
{
    return common_read_get_attr (fp, attrname, type, size, data);
}

int adios_get_attr_byid (ADIOS_FILE  * fp, int attrid, 
                    enum ADIOS_DATATYPES * type, int * size, void ** data)
{
    return common_read_get_attr_byid (fp, attrid, type, size, data);
}

const char * adios_type_to_string (enum ADIOS_DATATYPES type)
{
    return common_read_type_to_string (type);
}

int adios_type_size(enum ADIOS_DATATYPES type, void *data)
{
    return common_read_type_size(type, data);
}

int adios_get_grouplist (ADIOS_FILE  *fp, char ***group_namelist)
{
    return common_read_get_grouplist (fp, group_namelist);
}

int adios_group_view (ADIOS_FILE  *fp, int groupid)
{
    return common_read_group_view (fp, groupid);
}

void adios_print_fileinfo (ADIOS_FILE *fp) 
{
    common_read_print_fileinfo(fp);
}


ADIOS_SELECTION * adios_selection_boundingbox (int ndim, const uint64_t *start, const uint64_t *count)
{
    return common_read_selection_boundingbox (ndim, start, count);
}

ADIOS_SELECTION * adios_selection_points (int ndim, uint64_t npoints, const uint64_t *points)
{
    return common_read_selection_points (ndim, npoints, points);
}

ADIOS_SELECTION * adios_selection_writeblock (int index)
{
    return common_read_selection_writeblock (index);
}

ADIOS_SELECTION * adios_selection_auto (char *hints)
{
    return common_read_selection_auto (hints);
}

void adios_selection_delete (ADIOS_SELECTION *sel)
{
    common_read_selection_delete (sel);
}



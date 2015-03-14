#define DUMB_DEBUG
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
    char *data;
    struct _mask_read_request *next;
} mask_read_request;

static const int G_MAX_MASK_NUM = 8;

static int g_f_using_dataspaces = 0;

static int g_n_dims;
static int g_n_blocks;
static int g_mask_id = -1;

static int g_f_mask_done[G_MAX_MASK_NUM] = {};
static int **g_mask_bits[G_MAX_MASK_NUM];
static int **g_mask_index[G_MAX_MASK_NUM];

static int *g_mask_bits_offset[G_MAX_MASK_NUM];
static int *g_mask_bits_length[G_MAX_MASK_NUM];
static int *g_mask_offset[G_MAX_MASK_NUM];
static int *g_mask_length[G_MAX_MASK_NUM];

static int **g_block_offset[G_MAX_MASK_NUM];
static int **g_block_length[G_MAX_MASK_NUM];

static mask_read_request *g_requests = NULL;

static void add_mask_read_request(mask_read_request *req) {
    req->next = g_requests;
    g_requests = req;
}

const char *adios_errmsg ()
{
    return adios_get_last_errmsg();
}

int adios_set_using_dataspaces()
{
    int i;
    g_f_using_dataspaces = 1;
    g_mask_id = -1;
    memset(g_f_mask_done, 0, sizeof(g_f_mask_done));
    for ( i = 0; i < G_MAX_MASK_NUM; i++ )
        g_mask_bits[i] = NULL;
    return err_no_error;
}

int adios_read_init_method (enum ADIOS_READ_METHOD method, MPI_Comm comm, const char * parameters)
{
    if (method == ADIOS_READ_METHOD_DATASPACES) {
        g_f_using_dataspaces = 1;
    }
    g_mask_id = -1;
    memset(g_f_mask_done, 0, sizeof(g_f_mask_done));
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

// static int bits_decompress(uint64_t mask_bits_length, int *mask_bits, uint64_t *mask_length, char **mask) {
//     *mask_length = mask_bits_length * sizeof(int) * 8;
//     *mask = malloc(*mask_length * sizeof(char));
//     if (*mask == NULL) {
//         return 1;
//     }

//     uint64_t i, j, index = 0;
//     for (i = 0; i < mask_bits_length; ++i) {
//         for (j = 0; j < sizeof(int) * 8; ++j) {
//             (*mask)[index] = (mask_bits[i] >> j) & 1;
//             ++index;
//         }
//     }
//     return err_no_error;
// }

int adios_read_set_mask(ADIOS_FILE *fp, int mask_id, int n_dims, const char **offset_var_names, const char **length_var_names) {
    int i, j;
    char s[4][64];

    g_n_dims = n_dims;
    g_mask_id = mask_id;
    if (g_f_using_dataspaces) {
        sprintf(s[0], "const/proc_size");
        ADIOS_VARINFO *proc_size_info = adios_inq_var(fp, s[0]);
        g_n_blocks = *(int *)proc_size_info->value;
        adios_free_varinfo(proc_size_info);
    } else {
        sprintf(s[0], "mask_vars/mask_%d/mask_bits", g_mask_id);
        ADIOS_VARINFO *mask_info = adios_inq_var(fp, s[0]);
        g_n_blocks = mask_info->nblocks[0];
        adios_free_varinfo(mask_info);
    }

    if (g_f_mask_done[g_mask_id]) {
        return err_no_error;
    }

    if (g_f_using_dataspaces) {
        sprintf(s[0], "const/mask_%d/mask_bits_offset", g_mask_id);
        sprintf(s[1], "const/mask_%d/mask_bits_length", g_mask_id);
        sprintf(s[2], "const/mask_%d/mask_offset", g_mask_id);
        sprintf(s[3], "const/mask_%d/final_length", g_mask_id);
    } else {
        sprintf(s[0], "mask_vars/mask_%d/mask_bits_offset", g_mask_id);
        sprintf(s[1], "mask_vars/mask_%d/mask_bits_length", g_mask_id);
        sprintf(s[2], "mask_vars/mask_%d/mask_offset", g_mask_id);
        sprintf(s[3], "mask_vars/mask_%d/final_length", g_mask_id);
    }

    g_mask_bits_offset[g_mask_id] = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_bits_length[g_mask_id] = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_offset[g_mask_id] = (int *)malloc(sizeof(int) * g_n_blocks);
    g_mask_length[g_mask_id] = (int *)malloc(sizeof(int) * g_n_blocks);

    g_block_offset[g_mask_id] = (int **)malloc(sizeof(int *) * g_n_blocks);
    g_block_length[g_mask_id] = (int **)malloc(sizeof(int *) * g_n_blocks);
    for (i = 0; i < g_n_blocks; ++i) {
        g_block_offset[g_mask_id][i] = (int *)malloc(sizeof(int) * g_n_dims);
        g_block_length[g_mask_id][i] = (int *)malloc(sizeof(int) * g_n_dims);
    }
    if (g_f_using_dataspaces) {
        uint64_t t_start = 0, t_count = g_n_blocks;
        ADIOS_SELECTION *sel = adios_selection_boundingbox(1, &t_start, &t_count);
        adios_schedule_read(fp, sel, s[0], 0, 1, g_mask_bits_offset[g_mask_id]);
        adios_schedule_read(fp, sel, s[1], 0, 1, g_mask_bits_length[g_mask_id]);
        adios_schedule_read(fp, sel, s[2], 0, 1, g_mask_offset[g_mask_id]);
        adios_schedule_read(fp, sel, s[3], 0, 1, g_mask_length[g_mask_id]);
        int *offset_buffer = (int *)malloc(sizeof(int) * g_n_blocks * g_n_dims);
        int *length_buffer = (int *)malloc(sizeof(int) * g_n_blocks * g_n_dims);
        for (j = 0; j < g_n_dims; ++j) {
            adios_schedule_read(fp, sel, offset_var_names[j], 0, 1, &offset_buffer[j * g_n_blocks]);
            adios_schedule_read(fp, sel, length_var_names[j], 0, 1, &length_buffer[j * g_n_blocks]);
        }
        adios_selection_delete(sel);
        adios_perform_reads(fp, 1);
        for (i = 0; i < g_n_blocks; ++i) {
            for (j = 0; j < g_n_dims; ++j) {
                g_block_offset[g_mask_id][i][j] = offset_buffer[j * g_n_blocks + i];
                g_block_length[g_mask_id][i][j] = length_buffer[j * g_n_blocks + i];
            }
        }
        free(offset_buffer);
        free(length_buffer);
    } else {
        for (i = 0; i < g_n_blocks; ++i) {
            ADIOS_SELECTION *sel;
            sel = adios_selection_writeblock(i);
            adios_schedule_read(fp, sel, s[0], 0, 1, &g_mask_bits_offset[g_mask_id][i]);
            adios_schedule_read(fp, sel, s[1], 0, 1, &g_mask_bits_length[g_mask_id][i]);
            adios_schedule_read(fp, sel, s[2], 0, 1, &g_mask_offset[g_mask_id][i]);
            adios_schedule_read(fp, sel, s[3], 0, 1, &g_mask_length[g_mask_id][i]);
            for (j = 0; j < g_n_dims; ++j) {
                adios_schedule_read(fp, sel, offset_var_names[j], 0, 1, &g_block_offset[g_mask_id][i][j]);
                adios_schedule_read(fp, sel, length_var_names[j], 0, 1, &g_block_length[g_mask_id][i][j]);
            }
            adios_selection_delete(sel);
        }
        adios_perform_reads(fp, 1);
    }
    g_mask_bits[g_mask_id] = (int **)malloc(sizeof(int *) * g_n_blocks);
    g_mask_index[g_mask_id] = (int **)malloc(sizeof(int *) * g_n_blocks);
    for (i = 0; i < g_n_blocks; ++i) {
        g_mask_bits[g_mask_id][i] = NULL;
    }

    // for (i = 0; i < g_n_blocks; ++i) {
    //     printf("-------------------- Block %d --------------------\n", i);
    //     printf("Mask bits range: [%d, %d]\n", g_mask_bits_offset[g_mask_id][i], g_mask_bits_offset[g_mask_id][i] + g_mask_bits_length[g_mask_id][i] - 1);
    //     printf("Mask range: [%d, %d]\n", g_mask_offset[g_mask_id][i], g_mask_offset[g_mask_id][i] + g_mask_length[g_mask_id][i] - 1);
    //     printf("Original domain range: ");
    //     for (j = 0; j < g_n_dims; ++j) {
    //         if (j == 0) printf("[");
    //         printf("%d", g_block_offset[g_mask_id][i][j]);
    //         printf(j + 1 == g_n_dims ? "]" : ", ");
    //     }
    //     printf(" - ");
    //     for (j = 0; j < g_n_dims; ++j) {
    //         if (j == 0) printf("[");
    //         printf("%d", g_block_offset[g_mask_id][i][j] + g_block_length[g_mask_id][i][j] - 1);
    //         printf(j + 1 == g_n_dims ? "]" : ", ");
    //     }
    //     printf("\n");
    // }

    return err_no_error;
}

int adios_read_unset_mask(const ADIOS_FILE *fp) {
    g_mask_id = -1;

    mask_read_request *req;
    for (req = g_requests; req != NULL; req = g_requests) {
        g_requests = req->next;
        free(req->var_name);
        adios_selection_delete(req->sel);
        free(req);
    }

    return err_no_error;
}

int adios_read_remove_mask() {
    int i, j, k;

    for (k = 0; k < G_MAX_MASK_NUM; ++k) {
        if (g_mask_bits[k] == NULL) {
            continue;
        }

        for (i = 0; i < g_n_blocks; ++i) {
            if (g_mask_bits[k][i] != NULL) {
                free(g_mask_bits[k][i]);
                free(g_mask_index[k][i]);
            }
        }
        free(g_mask_bits[k]);
        free(g_mask_index[k]);
        g_mask_bits[k] = NULL;

        free(g_mask_bits_offset[k]);
        free(g_mask_bits_length[k]);
        free(g_mask_offset[k]);
        free(g_mask_length[k]);
        g_mask_bits_offset[k] = NULL;
        g_mask_bits_length[k] = NULL;
        g_mask_offset[k] = NULL;
        g_mask_length[k] = NULL;

        for (i = 0; i < g_n_blocks; ++i) {
            if (g_block_offset[k][i] != NULL ) {
                free(g_block_offset[k][i]);
            }
            if (g_block_length[k][i] != NULL ) {
                free(g_block_length[k][i]);
            }
        }
        free(g_block_offset[k]);
        free(g_block_length[k]);
        g_block_offset[k] = NULL;
        g_block_length[k] = NULL;
    }
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
    if (strncmp(varname, "mask_vars/", 10) == 0 || g_mask_id == -1 || g_mask_bits[g_mask_id] == NULL) {
        return common_read_schedule_read (fp, sel, varname, from_steps, nsteps, NULL, data);
    }

    int i;

    char s[64];
    sprintf(s, "mask_vars/mask_%d/mask_bits", g_mask_id);
    for (i = 0; i < g_n_blocks; ++i) {
        if (g_mask_bits[g_mask_id][i] == NULL && adios_boundingbox_touched(sel, g_block_offset[g_mask_id][i], g_block_length[g_mask_id][i])) {
            // printf("    Variable %s touched block %d\n", varname, i);
            g_mask_bits[g_mask_id][i] = (int *)malloc(sizeof(int) * g_mask_bits_length[g_mask_id][i]);
            g_mask_index[g_mask_id][i] = (int *)malloc(sizeof(int) * g_mask_bits_length[g_mask_id][i] * sizeof(int) * 8);
            uint64_t offset = g_mask_bits_offset[g_mask_id][i], length = g_mask_bits_length[g_mask_id][i];
            ADIOS_SELECTION *mask_bits_sel = adios_selection_boundingbox(1, &offset, &length);
            adios_schedule_read(fp, mask_bits_sel, s, from_steps, nsteps, g_mask_bits[g_mask_id][i]);
            adios_selection_delete(mask_bits_sel);
        }
    }

    mask_read_request *req = (mask_read_request *)malloc(sizeof(mask_read_request));
    req->var_name = strdup(varname);
    req->sel = copy_selection(sel);
    req->from_step = from_steps;
    req->n_steps = nsteps;
    req->data = data;
    add_mask_read_request(req);

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

static int adios_read_mask_var(const ADIOS_FILE *fp, mask_read_request *req) {
    printf("    Should read var: %s\n", req->var_name);

    char s[64];
    if( req->var_name[0] != '/' ) sprintf(s, "mask_vars/%s", req->var_name);
    else sprintf(s, "mask_vars%s", req->var_name);

    ADIOS_VARINFO *var_info = adios_inq_var((ADIOS_FILE *)fp, s);
    uint64_t element_size = adios_get_type_size(var_info->type, var_info->value);
    adios_free_varinfo(var_info);

    int i, j, k;
    uint64_t index;

    // char **mask = (char **)malloc(sizeof(char *) * g_n_blocks);
    // for (i = 0; i < g_n_blocks; ++i) {
    //     if (g_mask_bits[i] == NULL) {
    //         mask[i] = NULL;
    //     } else {
    //         uint64_t mask_length = 0;
    //         bits_decompress(g_mask_bits_length[i], g_mask_bits[i], &mask_length, &mask[i]);
    //     }
    // }

    double t1, t2, t3;
#ifdef DUMB_DEBUG
    t1 = MPI_Wtime();
#endif
    char **mask_vars = (char **)malloc(sizeof(char *) * g_n_blocks);
    for (i = 0; i < g_n_blocks; ++i) {
        mask_vars[i] = NULL;
        if (adios_boundingbox_touched(req->sel, g_block_offset[g_mask_id][i], g_block_length[g_mask_id][i])) {
            mask_vars[i] = (char *)malloc(sizeof(char) * element_size * g_mask_length[g_mask_id][i]); /* TODO: take care of steps */
            uint64_t offset = g_mask_offset[g_mask_id][i], length = g_mask_length[g_mask_id][i];
            ADIOS_SELECTION *mask_sel = adios_selection_boundingbox(1, &offset, &length);
            adios_schedule_read(fp, mask_sel, s, req->from_step, req->n_steps, mask_vars[i]);
            adios_selection_delete(mask_sel);
        }
    }
#ifdef DUMB_DEBUG
    t2 = MPI_Wtime();
#endif
    common_read_perform_reads(fp, 1);
#ifdef DUMB_DEBUG
    t3 = MPI_Wtime();
#endif
#ifdef DUMB_DEBUG
    printf("[%s]: schedule var %s all data bits: %.8f\n", __func__, req->var_name, t2 - t1);
    printf("[%s]: read var %s all data bits: %.8f\n", __func__, req->var_name, t3 - t2);
#endif

    for (i = 0; i < g_n_blocks; ++i) {
        if (mask_vars[i] == NULL) {
            continue;
        }
        uint64_t var_length = 1;
        for (j = 0; j < g_n_dims; ++j) {
            var_length *= g_block_length[g_mask_id][i][j];
        }
        // char *var_data = (char *)malloc(sizeof(char) * element_size * var_length);

        int from[10], to[10];
        int pos[10];
        uint64_t dim_offset_block[11];
        uint64_t dim_offset_req[11];
        dim_offset_block[g_n_dims] = 1;
        dim_offset_req[g_n_dims] = 1;
        uint64_t index_block = 0, index_req = 0;
        for (k = g_n_dims - 1; k >= 0; --k) {
            from[k] = g_block_offset[g_mask_id][i][k];
            int t1 = req->sel->u.bb.start[k];
            if (t1 > from[k]) {
                from[k] = t1;
            }
            pos[k] = from[k];
            to[k] = g_block_offset[g_mask_id][i][k] + g_block_length[g_mask_id][i][k];
            int t2 = req->sel->u.bb.start[k] + req->sel->u.bb.count[k];
            if (t2 < to[k]) {
                to[k] = t2;
            }
            dim_offset_block[k] = dim_offset_block[k + 1] * g_block_length[g_mask_id][i][k];
            dim_offset_req[k] = dim_offset_req[k + 1] * req->sel->u.bb.count[k];
            index_block += (pos[k] - g_block_offset[g_mask_id][i][k]) * dim_offset_block[k + 1];
            index_req += (pos[k] - req->sel->u.bb.start[k]) * dim_offset_req[k + 1];
        }
#ifdef DUMB_DEBUG
        t1 = MPI_Wtime();
#endif
        while (true) {
            /*int index_block = 0;
            for (k = 0; k < g_n_dims; ++k) {
                index_block = index_block * g_block_length[g_mask_id][i][k] + pos[k] - g_block_offset[g_mask_id][i][k];
            }
            int index_req = 0;
            for (k = 0; k < g_n_dims; ++k) {
                index_req = index_req * req->sel->u.bb.count[k] + pos[k] - req->sel->u.bb.start[k];
            }*/
            int p = index_block / (sizeof(int) * 8), q = index_block % (sizeof(int) * 8);
            if (((g_mask_bits[g_mask_id][i][p] >> q) & 1) == 1) {
                memcpy(req->data + element_size * index_req, (char *)mask_vars[i] + element_size * g_mask_index[g_mask_id][i][index_block], element_size);
            } else {
                memset(req->data + element_size * index_req, 0, element_size);
            }
            for (k = 0; k < g_n_dims; ++k) {
                if (++pos[k] == to[k]) {
                    index_block -= (pos[k] - from[k] - 1) * dim_offset_block[k + 1];
                    index_req -= (pos[k] - from[k] - 1) * dim_offset_req[k + 1];
                    pos[k] = from[k];
                } else {
                    index_block += dim_offset_block[k + 1];
                    index_req += dim_offset_req[k + 1];
                    break;
                }
            }
            if (k == g_n_dims) {
                break;
            }
        }
#ifdef DUMB_DEBUG
        t2 = MPI_Wtime();
#endif
#ifdef DUMB_DEBUG
        printf("[%s]: reorganize block %d: %.8f\n", __func__, i, t2 - t1);
#endif

        // index = 0;
        // for (j = 0; j < var_length; ++j) {
        //     uint64_t p = j / (sizeof(int) * 8), q = j % (sizeof(int) * 8);
        //     if (((g_mask_bits[g_mask_id][i][p] >> q) & 1) == 1) {
        //         memcpy(var_data + element_size * j, (char *)mask_vars[i] + index, element_size);
        //         index += element_size;
        //     } else {
        //         memset(var_data + element_size * j, 0, element_size);
        //     }
        // }

        // uint64_t start[10], end[10];
        // for (j = 0; j < g_n_dims; ++j) {
        //     start[j] = req->sel->u.bb.start[j];
        //     if (g_block_offset[i][j] > start[j]) {
        //         start[j] = g_block_offset[i][j];
        //     }
        //     end[j] = req->sel->u.bb.start[j] + req->sel->u.bb.count[j] - 1;
        //     if (g_block_offset[i][j] + g_block_length[i][j] - 1 < end[j]) {
        //         end[j] = g_block_offset[i][j] + g_block_length[i][j] - 1;
        //     }
        // }
        // for (j = 0; j < var_length; ++j) {
        //     int offset[10];
        //     int pos = j;
        //     int f_in_domain = 1;
        //     for (k = g_n_dims - 1; k >= 0; --k) {
        //         offset[k] = g_block_offset[g_mask_id][i][k] + pos % g_block_length[g_mask_id][i][k] - req->sel->u.bb.start[k];
        //         if (offset[k] < 0 || offset[k] >= req->sel->u.bb.count[k]) {
        //             f_in_domain = 0;
        //             break;
        //         }
        //         pos /= g_block_length[g_mask_id][i][k];
        //     }
        //     if (f_in_domain) {
        //         pos = 0;
        //         for (k = 0; k < g_n_dims; ++k) {
        //             pos = pos * req->sel->u.bb.count[k] + offset[k];
        //         }
        //         memcpy(req->data + element_size * pos, var_data + element_size * j, element_size);
        //     }
        // }

        // free(var_data);
    }

    for (i = 0; i < g_n_blocks; ++i) {
        // if (mask[i] != NULL) {
        //     free(mask[i]);
        // }
        if (mask_vars[i] != NULL) {
            free(mask_vars[i]);
        }
    }
    // free(mask);
    free(mask_vars);

    return err_no_error;
}

int adios_perform_reads (const ADIOS_FILE *fp, int blocking)
{
    if (g_mask_id == -1 || g_mask_bits[g_mask_id] == NULL) {
        return common_read_perform_reads (fp, blocking);
    }

    double t1, t2, t3;
    int i, j, k;
#ifdef DUMB_DEBUG
    t1 = MPI_Wtime();
#endif
    common_read_perform_reads (fp, blocking);
#ifdef DUMB_DEBUG
    t2 = MPI_Wtime();
#endif
    for (i = 0; i < g_n_blocks; ++i) {
        if (g_mask_bits[g_mask_id][i] == NULL) {
            continue;
        }
        int index = 0;
        g_mask_index[g_mask_id][i][index] = 0;
        for (j = 0; j < g_mask_bits_length[g_mask_id][i]; ++j) {
            for (k = 0; k < sizeof(int) * 8; ++k) {
                g_mask_index[g_mask_id][i][index + 1] = g_mask_index[g_mask_id][i][index] + ((g_mask_bits[g_mask_id][i][j] >> k) & 1);
                ++index;
            }
        }
    }
#ifdef DUMB_DEBUG
    t3 = MPI_Wtime();
#endif
    g_f_mask_done[g_mask_id] = 1;

    mask_read_request *req;
    for (req = g_requests; req != NULL; req = req->next) {
        adios_read_mask_var(fp, req);
    }
#ifdef DUMB_DEBUG
    printf("[%s]: read all mask bits: %.8f\n", __func__, t2 - t1);
    printf("[%s]: make index for all mask bits: %.8f\n", __func__, t3 - t2);
#endif
    return err_no_error;
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



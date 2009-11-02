#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "adios.h"

/*************************************************************/
/*          Example of writing arrays in ADIOS               */
/*************************************************************/
int main (int argc, char ** argv) 
{
    char        filename [256];
    int         rank, size, i, j;
    int         NX = 10, NY = 100; 
    double      t[NX][NY];
    int         p[NX];
    MPI_Comm    comm = MPI_COMM_WORLD;

    int         adios_err;
    uint64_t    adios_groupsize, adios_totalsize;
    int64_t     adios_handle;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (comm, &rank);

    for (i = 0; i < NX; i++)
        for (j = 0; j< NY; j++)
            t[i][j] = rank * 100 + 10 * i + j;

    for (i = 0; i < NX; i++)
        p[i] = rank * 100 + 10 * i;

    strcpy (filename, "restart.bp");
    adios_init ("config_array.xml");
    adios_open (&adios_handle, "my_group", filename, "w");

    adios_groupsize = sizeof (int)              \
                    + sizeof (int)              \
                    + NX * NY * sizeof (double) \
                    + NX * sizeof (int);
 
    adios_group_size (adios_handle, adios_groupsize, &adios_totalsize, &comm);

    adios_write (adios_handle, "NX", &NX);
    adios_write (adios_handle, "NY", &NY);
    adios_write (adios_handle, "var_double_array", t);
    adios_write (adios_handle, "var_int_array", p);

    adios_close (adios_handle);

    MPI_Barrier (comm);

    adios_finalize (rank);

    MPI_Finalize ();
    return 0;
}

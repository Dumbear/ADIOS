#!/usr/bin/env python
from adios_mpi import *
from scipy.io import netcdf
import numpy as np
import sys
import os
import operator

def usage():
    print os.path.basename(sys.argv[0]), "netcdf_file"

if len(sys.argv) < 2:
    usage()
    sys.exit(0)

##fname = "MERRA100.prod.assim.tavg3_3d_mst_Cp.19791010.SUB.nc"

fname = sys.argv[1]
fout = '.'.join(fname.split('.')[:-1]) + ".bp"

tname = "time"
if len(sys.argv) > 2:
    tname = sys.argv[2]

## Open NetCDF file
f = netcdf.netcdf_file(fname, 'r')

## Check dimension
assert (all(map(lambda x: x is not None,
                [ val for k, val in f.dimensions.items()
                  if k != tname])))

## Two types of variables : time-dependent or time-independent
dimvar = {n:v for n,v in f.variables.items() if n in f.dimensions.keys()}
var = {n:v for n,v in f.variables.items() if n not in f.dimensions.keys()}
tdepvar = {n:v for n,v in var.items() if tname in v.dimensions}
tindvar = {n:v for n,v in var.items() if tname not in v.dimensions}

## Time dimension
assert (len(set([v.dimensions.index(tname) for v in tdepvar.values()]))==1)
tdx = tdepvar.values()[0].dimensions.index(tname)

assert (all([v.data.shape[tdx] for v in tdepvar.values()]))
tdim = tdepvar.values()[0].shape[tdx]

## Init ADIOS without xml
init_noxml()
allocate_buffer(BUFFER_ALLOC_WHEN.NOW, 10)
gid = declare_group ("group", tname, FLAG.YES)
select_method (gid, "MPI", "", "")

d1size = 0
for name, val in f.dimensions.items():
    if name == tname:
        continue
    print "Dimension : %s (%d)" % (name, val)
    define_var (gid, name, "", DATATYPE.integer, "", "", "")
    d1size += 4

"""
d2size = 0
for name, var in dimvar.items():
    if name == tname:
        continue
    if name in f.dimensions.keys():
        name = "v_" + name
    print "Variable : %s (%s)" % (name, ','.join(var.dimensions))
    define_var (gid, name, "", DATATYPE.double,
                ','.join(var.dimensions),
                "",
                "")
    d2size += reduce(operator.mul, var.shape) * 8

v1size = 0
for name, var in tindvar.items():
    print "Variable : %s (%s)" % (name, ','.join(var.dimensions))
    define_var (gid, name, "", DATATYPE.double,
                ','.join(var.dimensions),
                "",
                "")
    v1size += reduce(operator.mul, var.shape) * 8
"""
    
v2size = 0
for name, var in tdepvar.items():
    print "Variable : %s (%s)" % (name, ','.join(var.dimensions))
    define_var (gid, name, "", DATATYPE.double,
                ','.join(var.dimensions),
                ','.join([dname for dname in var.dimensions
                          if dname != tname]),
                "0,0,0")
    v2size += reduce(operator.mul, var.shape) / tdim * 8

print "Count (dim, var) : ", (d1size, v2size)

## Clean old file
if os.access(fout, os.F_OK):
    os.remove(fout)

for it in range(tdim):
    print 
    print "Time step : %d" % (it)
    
    fd = open("group", fout, "a")
    groupsize = d1size + v2size
    set_group_size(fd, groupsize)

    for name, val in f.dimensions.items():
        if name == tname:
            continue
        print "Dimension writing : %s (%d)" % (name, val)
        write_int(fd, name, val)
        
    for name, var in tdepvar.items():
        try:
            arr = np.array(var.data.take([it], axis=tdx),
                           dtype=np.float64)
            print "Variable writing : %s %s" % (name, arr.shape)
            write(fd, name, arr)
        except ValueError:
            print "Skip:", name

    close(fd)
    
f.close()
finalize()

print
print "Done. Saved:", fout

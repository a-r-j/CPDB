# pdb_parser.pyx

cimport cython
from libc.stdio cimport FILE, fopen, fclose, fgets
from libc.stdlib cimport malloc, free
from libc.stdlib cimport atof, atoi

import numpy as np
cimport numpy as np

ctypedef np.float32_t FLOAT32_t
ctypedef np.int32_t INT32_t

@cython.boundscheck(False)
@cython.wraparound(False)
def parse_pdb(filename, int max_atoms=1000000):
    cdef FILE *fp = NULL
    cdef char *line = NULL
    cdef int n_atoms = 0
    cdef int i = 0

    # Preallocate arrays for parsed data
    cdef np.ndarray[FLOAT32_t, ndim=1] x = np.empty(max_atoms, dtype=np.float32)
    cdef np.ndarray[FLOAT32_t, ndim=1] y = np.empty(max_atoms, dtype=np.float32)
    cdef np.ndarray[FLOAT32_t, ndim=1] z = np.empty(max_atoms, dtype=np.float32)
    cdef np.ndarray[FLOAT32_t, ndim=1] occupancy = np.empty(max_atoms, dtype=np.float32)
    cdef np.ndarray[FLOAT32_t, ndim=1] tempFactor = np.empty(max_atoms, dtype=np.float32)

    cdef np.ndarray[INT32_t, ndim=1] atom_serial = np.empty(max_atoms, dtype=np.int32)
    cdef np.ndarray[INT32_t, ndim=1] resNum = np.empty(max_atoms, dtype=np.int32)


    cdef list record_type = [''] * max_atoms
    cdef list atomName = [''] * max_atoms
    cdef list altLoc = [''] * max_atoms
    cdef list residueName = [''] * max_atoms
    cdef list chainID = [''] * max_atoms
    cdef list iCode = [''] * max_atoms
    cdef list elementSymbol = [''] * max_atoms
    cdef list charge = [''] * max_atoms

    line = <char *> malloc(82 * cython.sizeof(char))  # 82 chars to include newline and null terminator
    fp = fopen(filename.encode('utf-8'), "r")
    if fp == NULL:
        raise FileNotFoundError("Could not open file")

    while fgets(line, 82, fp) != NULL:
        if line[:4] == b'ATOM' or line[:6] == b'HETATM':
            record_type[n_atoms] = line[:6].decode('utf-8').strip()
            atom_serial[n_atoms] = atoi(line[6:11])
            atomName[n_atoms] = line[12:16].decode('utf-8').strip()
            altLoc[n_atoms] = line[16:17].decode('utf-8').strip()
            residueName[n_atoms] = line[17:20].decode('utf-8').strip()
            chainID[n_atoms] = line[21:22].decode('utf-8').strip()
            resNum[n_atoms] = atoi(line[22:26])
            x[n_atoms] = atof(line[30:38])
            y[n_atoms] = atof(line[38:46])
            z[n_atoms] = atof(line[46:54])
            occupancy[n_atoms] = atof(line[54:60])
            tempFactor[n_atoms] = atof(line[60:66])
            elementSymbol[n_atoms] = line[76:78].decode('utf-8').strip()
            charge[n_atoms] = line[78:80].decode('utf-8').strip()
            n_atoms += 1

            if n_atoms >= max_atoms:
                break

    fclose(fp)
    free(line)

    return {
        'record': np.array(record_type[:n_atoms], dtype='U6'),
        'atom_serial': atom_serial[:n_atoms],
        'atom_name': np.array(atomName[:n_atoms], dtype='U4'),
        'altLoc': np.array(altLoc[:n_atoms], dtype='U1'),
        'residueName': np.array(residueName[:n_atoms], dtype='U3'),
        'chainID': np.array(chainID[:n_atoms], dtype='U1'),
        'resNum': resNum[:n_atoms],
        'iCode': np.array(iCode[:n_atoms], dtype='U1'),
        'x': x[:n_atoms],
        'y': y[:n_atoms],
        'z': z[:n_atoms],
        'occupancy': occupancy[:n_atoms],
        'tempFactor': tempFactor[:n_atoms],
        'elementSymbol': np.array(elementSymbol[:n_atoms], dtype='U2'),
    }

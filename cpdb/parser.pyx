cimport cython
from libc.stdio cimport FILE, fopen, fclose, fgets
from libc.stdlib cimport malloc, free
from libc.stdlib cimport atof, atoi
from libc.string cimport strncpy

cdef struct AtomData:
    char[7] record_type
    int atom_serial
    char[4] atom_name
    char[2] alt_loc
    char[4] residue_name
    char[2] chain_id
    int res_num
    char[2] iCode
    float x
    float y
    float z
    float occupancy
    float temp_factor
    char[3] element_symbol
    char[2] charge

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def parse_pdb(filename, int max_atoms=1000000):
    cdef FILE *fp = NULL
    cdef char *line = NULL
    cdef int n_atoms = 0
    cdef int i = 0

    cdef AtomData *atom_data = <AtomData *> malloc(max_atoms * cython.sizeof(AtomData))

    line = <char *> malloc(82 * cython.sizeof(char))  # 82 chars to include newline and null terminator
    fp = fopen(filename.encode('utf-8'), "r, ccs=UTF-8")
    if fp == NULL:
        raise FileNotFoundError("Could not open file")

    while fgets(line, 82, fp) != NULL: # 82
        if line[:4] == b'ATOM' or line[:6] == b'HETATM':
            strncpy(atom_data[n_atoms].record_type, line, 6)
            atom_data[n_atoms].record_type[6] = '\0'
            atom_data[n_atoms].atom_serial = atoi(line[6:11])
            strncpy(atom_data[n_atoms].atom_name, line[12:16], 4)
            strncpy(atom_data[n_atoms].alt_loc, line[16:17], 1)
            strncpy(atom_data[n_atoms].residue_name, line[17:20], 3)
            strncpy(atom_data[n_atoms].chain_id, line[21:22], 1)
            atom_data[n_atoms].res_num = atoi(line[22:26])
            strncpy(atom_data[n_atoms].iCode, line[26:27], 1)
            atom_data[n_atoms].x = atof(line[30:38])
            atom_data[n_atoms].y = atof(line[38:46])
            atom_data[n_atoms].z = atof(line[46:54])
            atom_data[n_atoms].occupancy = atof(line[54:60])
            atom_data[n_atoms].temp_factor = atof(line[60:66])
            strncpy(atom_data[n_atoms].element_symbol, line[76:78], 2)
            strncpy(atom_data[n_atoms].charge, line[78:79], 2)
            n_atoms += 1

            if n_atoms >= max_atoms:
                break
    fclose(fp)
    free(line)

    # Convert the C struct data to Python objects and return
    return {
        'record_name': [atom_data[i].record_type.decode('ascii', errors='replace').strip() for i in range(n_atoms)],
        'atom_number': [atom_data[i].atom_serial for i in range(n_atoms)],
        'atom_name': [atom_data[i].atom_name.decode('utf-8').strip() for i in range(n_atoms)],
        'alt_loc': [atom_data[i].alt_loc.decode('utf-8').strip() for i in range(n_atoms)],
        'residue_name': [atom_data[i].residue_name.decode('utf-8', errors='replace').strip() for i in range(n_atoms)],
        'chain_id': [atom_data[i].chain_id.decode('utf-8', errors='replace').strip() for i in range(n_atoms)],
        'residue_number': [atom_data[i].res_num for i in range(n_atoms)],
        "insertion": [atom_data[i].iCode.decode('utf-8', errors='replace').strip() for i in range(n_atoms)],
        'x_coord': [atom_data[i].x for i in range(n_atoms)],
        'y_coord': [atom_data[i].y for i in range(n_atoms)],
        'z_coord': [atom_data[i].z for i in range(n_atoms)],
        'occupancy': [atom_data[i].occupancy for i in range(n_atoms)],
        'b_factor': [atom_data[i].temp_factor for i in range(n_atoms)],
        'element_symbol': [atom_data[i].element_symbol.decode('utf-8', errors='replace').strip() for i in range(n_atoms)],
        'charge': [atom_data[i].charge.decode('utf-8').strip() for i in range(n_atoms)],
    }

    # Free the allocated memory for the atom_data array
    free(atom_data)

cimport cython
cimport numpy as cnp
import numpy as np
from libc.stdio cimport FILE, fopen, fclose, fgets
from libc.stdlib cimport malloc, free
from libc.stdlib cimport atof, atoi
from libc.string cimport strncpy, strlen
from libc.stdint cimport uint8_t


cdef struct AtomData:
    char[7] record_type
    int atom_serial
    char[5] atom_name
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
    int model_idx


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def parse_pdb(filename, int max_atoms=1000000) -> object:
    cdef FILE *fp = NULL
    cdef char *line = NULL
    cdef int n_atoms = 0
    cdef int i = 0
    cdef int n_model =1

    cdef AtomData *atom_data = <AtomData *> malloc(max_atoms * cython.sizeof(AtomData))

    line = <char *> malloc(82 * cython.sizeof(char))  # 82 chars to include newline and null terminator
    fp = fopen(filename.encode('utf-8'), "r, ccs=UTF-8")
    if fp == NULL:
        raise FileNotFoundError("Could not open file")

    while fgets(line, 82, fp) != NULL: # 82
        if line[:4] == b'ATOM' or line[:6] == b'HETATM':
            strncpy(atom_data[n_atoms].record_type, line, 6)
            atom_data[n_atoms].record_type[6] = b'\0' # Add Null terminator

            atom_data[n_atoms].atom_serial = atoi(line[6:11])

            strncpy(atom_data[n_atoms].atom_name, line[12:16], 4)
            atom_data[n_atoms].atom_name[4] = b'\0'

            strncpy(atom_data[n_atoms].alt_loc, line[16:17], 1)
            atom_data[n_atoms].alt_loc[1] = b'\0'

            strncpy(atom_data[n_atoms].residue_name, line[17:20], 3)
            atom_data[n_atoms].residue_name[3] = b'\0'

            strncpy(atom_data[n_atoms].chain_id, line[21:22], 1)
            atom_data[n_atoms].chain_id[1] = b'\0'
            
            atom_data[n_atoms].res_num = atoi(line[22:26])
            
            strncpy(atom_data[n_atoms].iCode, line[26:27], 1)
            atom_data[n_atoms].iCode[1] = b'\0'

            atom_data[n_atoms].x = atof(line[30:38])
            atom_data[n_atoms].y = atof(line[38:46])
            atom_data[n_atoms].z = atof(line[46:54])

            atom_data[n_atoms].occupancy = atof(line[54:60])
            atom_data[n_atoms].temp_factor = atof(line[60:66])
            
            strncpy(atom_data[n_atoms].element_symbol, line[76:78], 2)
            atom_data[n_atoms].element_symbol[2] = b'\0'

            strncpy(atom_data[n_atoms].charge, line[78:79], 2)
            atom_data[n_atoms].model_idx = n_model
            n_atoms += 1

            if n_atoms >= max_atoms:
                break
        elif line[:6] == b'ENDMDL':
            n_model += 1

    fclose(fp)
    free(line)

    # Create memory view arrays
    cdef cnp.ndarray[float, ndim=1] x_coords = np.empty(n_atoms, dtype=np.float32)
    cdef cnp.ndarray[float, ndim=1] y_coords = np.empty(n_atoms, dtype=np.float32)
    cdef cnp.ndarray[float, ndim=1] z_coords = np.empty(n_atoms, dtype=np.float32)
    cdef cnp.ndarray[float, ndim=1] occupancies = np.empty(n_atoms, dtype=np.float32)
    cdef cnp.ndarray[float, ndim=1] b_factors = np.empty(n_atoms, dtype=np.float32)
    cdef cnp.ndarray[int, ndim=1] atom_number = np.empty(n_atoms, dtype=np.int32)
    cdef cnp.ndarray[int, ndim=1] residue_number = np.empty(n_atoms, dtype=np.int32)
    cdef cnp.ndarray[int, ndim=1] model_idx = np.empty(n_atoms, dtype=np.int32)
    
    # Create NumPy arrays for string columns
    cdef cnp.ndarray[object, ndim=1] record_name = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] atom_name = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] alt_loc = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] residue_name = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] chain_id = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] insertion = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] element_symbol = np.empty(n_atoms, dtype=object)
    cdef cnp.ndarray[object, ndim=1] charge = np.empty(n_atoms, dtype=object)

    # Fill memory view arrays
    for i in range(n_atoms):
        x_coords[i] = atom_data[i].x
        y_coords[i] = atom_data[i].y
        z_coords[i] = atom_data[i].z
        occupancies[i] = atom_data[i].occupancy
        b_factors[i] = atom_data[i].temp_factor
        model_idx[i] = atom_data[i].model_idx
        atom_number[i] = atom_data[i].atom_serial
        residue_number[i] = atom_data[i].res_num
        atom_name[i] = ATOM_NUMBERING.get()

        # Fill NumPy arrays for string columns
        record_name[i] = atom_data[i].record_type[:strlen(atom_data[i].record_type)].decode('ascii', errors='replace').strip()
        atom_name[i] = atom_data[i].atom_name[:strlen(atom_data[i].atom_name)].decode('utf-8').strip()
        alt_loc[i] = atom_data[i].alt_loc[:strlen(atom_data[i].alt_loc)].decode('utf-8').strip()
        residue_name[i] = atom_data[i].residue_name[:strlen(atom_data[i].residue_name)].decode('utf-8', errors='replace').strip()
        chain_id[i] = atom_data[i].chain_id[:strlen(atom_data[i].chain_id)].decode('utf-8', errors='replace').strip()
        insertion[i] = atom_data[i].iCode[:strlen(atom_data[i].iCode)].decode('utf-8', errors='replace').strip()
        element_symbol[i] = atom_data[i].element_symbol[:strlen(atom_data[i].element_symbol)].decode('utf-8', errors='replace').strip()
        charge[i] = atom_data[i].charge[:strlen(atom_data[i].charge)].decode('utf-8').strip()


    # Convert the C struct data to Python objects and return
    result = {
        'record_name': record_name,
        'atom_number': atom_number,
        'atom_name': atom_name,
        'alt_loc': alt_loc,
        'residue_name': residue_name,
        'chain_id': chain_id,
        'residue_number': residue_number,
        "insertion": insertion,
        'x_coord': x_coords,
        'y_coord': y_coords,
        'z_coord': z_coords,
        'occupancy': occupancies,
        'b_factor': b_factors,
        'element_symbol': element_symbol,
        'charge': charge,
        'model_idx': model_idx
    }

    # Free the allocated memory for the atom_data array
    free(atom_data)

    return result
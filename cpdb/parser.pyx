cimport cython
cimport numpy as cnp

import numpy as np

from libc.stdint cimport uint64_t
from libc.stdlib cimport atof, atoi, free, malloc, realloc
from libc.string cimport memchr, memcmp, memcpy, memset


cdef struct AtomData:
    bint hetatm
    int atom_serial
    float x
    float y
    float z
    float occupancy
    float temp_factor
    int res_num
    int model_idx
    # Raw fixed-width text fields (see TEXT_FIELDS), packed into ints so they can be dict keys
    uint64_t text[7]


# (column name, start, width) of each fixed-width text field, in AtomData.text order
TEXT_FIELDS = (
    ('atom_name', 12, 4),
    ('alt_loc', 16, 1),
    ('residue_name', 17, 3),
    ('chain_id', 21, 1),
    ('insertion', 26, 1),
    ('element_symbol', 76, 2),
    ('charge', 78, 2),
)

cdef int _TEXT_START[7]
cdef int _TEXT_WIDTH[7]
for _i, (_, _start, _width) in enumerate(TEXT_FIELDS):
    _TEXT_START[_i] = _start
    _TEXT_WIDTH[_i] = _width

COLUMNS = ('record_name', 'atom_number', 'atom_name', 'alt_loc', 'residue_name', 'chain_id',
           'residue_number', 'insertion', 'x_coord', 'y_coord', 'z_coord', 'occupancy',
           'b_factor', 'element_symbol', 'charge', 'model_idx')


cdef inline uint64_t _field_key(const char *line, Py_ssize_t n, int start, int width) noexcept nogil:
    # Pack columns [start, start + width) into an int; columns past the end of a short line count as blank
    cdef uint64_t key = 0
    cdef int i
    for i in range(width):
        key = (key << 8) | (<unsigned char> line[start + i] if start + i < n else 32)
    return key


cdef double[19] _POW10 = [1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12,
                          1e13, 1e14, 1e15, 1e16, 1e17, 1e18]


@cython.cdivision(True)
cdef inline double _field_float(const char *line, Py_ssize_t n, int start, int width) noexcept nogil:
    # Fast path for PDB's fixed-point "  -12.345"; anything else (exponents, junk) falls back to atof.
    # mantissa / 10**k divides two exactly representable doubles, so it rounds the same as atof.
    cdef char buf[16]
    cdef int w = min(width, max(0, <int> n - start))
    cdef int i = 0
    cdef int n_digits = 0
    cdef int n_frac = 0
    cdef bint neg = False
    cdef bint dot = False
    cdef long long mantissa = 0
    cdef char c
    while i < w and line[start + i] == c' ':
        i += 1
    if i < w and (line[start + i] == c'-' or line[start + i] == c'+'):
        neg = line[start + i] == c'-'
        i += 1
    while i < w:
        c = line[start + i]
        if c'0' <= c <= c'9':
            mantissa = mantissa * 10 + (c - c'0')
            n_digits += 1
            n_frac += dot
        elif c == c'.' and not dot:
            dot = True
        else:
            break
        i += 1
    while i < w and line[start + i] == c' ':
        i += 1
    if i == w and 0 < n_digits <= 18:
        return -(mantissa / _POW10[n_frac]) if neg else mantissa / _POW10[n_frac]
    memcpy(buf, line + start, w)
    buf[w] = 0
    return atof(buf)


cdef inline int _field_int(const char *line, Py_ssize_t n, int start, int width) noexcept nogil:
    cdef char buf[16]
    cdef int w = min(width, max(0, <int> n - start))
    memcpy(buf, line + start, w)
    buf[w] = 0
    return atoi(buf)


cdef inline Py_ssize_t _slot(uint64_t key, Py_ssize_t mask) noexcept nogil:
    return <Py_ssize_t> ((key * <uint64_t> 0x9E3779B97F4A7C15) >> 32) & mask


@cython.boundscheck(False)
@cython.wraparound(False)
cdef object _decode_column(AtomData *atom_data, Py_ssize_t n_atoms, int field, int width):
    # Text columns have few distinct values: map each raw key to a code with a C open-addressing
    # table (slot_code -1 = empty), then decode + strip only the distinct keys
    codes = np.empty(n_atoms, dtype=np.intp)
    cdef Py_ssize_t[::1] cv = codes
    cdef Py_ssize_t capacity = 256
    cdef Py_ssize_t mask, s, code
    cdef Py_ssize_t i = 0
    cdef Py_ssize_t n_unique = 0
    cdef uint64_t key
    cdef uint64_t *unique_keys = NULL
    cdef uint64_t *grown_keys
    cdef Py_ssize_t *slot_code = NULL
    cdef int j
    cdef char raw[8]
    try:
        while True:
            # (Re)build an empty table at the current capacity and re-insert the keys seen so far
            free(slot_code)
            slot_code = <Py_ssize_t *> malloc(capacity * sizeof(Py_ssize_t))
            grown_keys = <uint64_t *> realloc(unique_keys, (capacity // 2) * sizeof(uint64_t))
            if slot_code == NULL or grown_keys == NULL:
                raise MemoryError()
            unique_keys = grown_keys
            mask = capacity - 1
            memset(slot_code, 0xFF, capacity * sizeof(Py_ssize_t))
            for code in range(n_unique):
                s = _slot(unique_keys[code], mask)
                while slot_code[s] != -1:
                    s = (s + 1) & mask
                slot_code[s] = code

            with nogil:
                while i < n_atoms:
                    key = atom_data[i].text[field]
                    s = _slot(key, mask)
                    while slot_code[s] != -1 and unique_keys[slot_code[s]] != key:
                        s = (s + 1) & mask
                    if slot_code[s] == -1:
                        if n_unique * 2 >= capacity:
                            break  # keep load <= 1/2: grow, then resume at atom i
                        slot_code[s] = n_unique
                        unique_keys[n_unique] = key
                        n_unique += 1
                    cv[i] = slot_code[s]
                    i += 1
            if i == n_atoms:
                break
            capacity *= 2

        values = []
        for code in range(n_unique):
            for j in range(width):
                raw[j] = <char> (unique_keys[code] >> (8 * (width - 1 - j)))
            values.append(raw[:width].decode('utf-8', errors='replace').strip())
    finally:
        free(slot_code)
        free(unique_keys)
    out = np.empty(len(values), dtype=object)
    out[:] = values
    return out[codes]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef dict _parse(const char *buf, Py_ssize_t size, Py_ssize_t max_atoms):
    cdef const char *line = buf
    cdef const char *end = buf + size
    cdef const char *nl
    cdef Py_ssize_t n
    cdef Py_ssize_t n_atoms = 0
    cdef Py_ssize_t capacity = 0
    cdef int n_model = 1
    cdef Py_ssize_t i
    cdef AtomData *atom_data = NULL
    cdef AtomData *a
    cdef AtomData *grown
    cdef int f
    cdef float[::1] xv, yv, zv, ov, bv
    cdef int[::1] anv, rnv, mv
    cdef unsigned char[::1] hv

    try:
        while line < end and (max_atoms < 0 or n_atoms < max_atoms):
            nl = <const char *> memchr(line, b'\n', end - line)
            n = (nl if nl != NULL else end) - line
            if n > 0 and line[n - 1] == b'\r':
                n -= 1

            if (n >= 4 and memcmp(line, b"ATOM", 4) == 0) or (n >= 6 and memcmp(line, b"HETATM", 6) == 0):
                if n_atoms == capacity:
                    capacity = max(1024, capacity * 2)
                    grown = <AtomData *> realloc(atom_data, capacity * sizeof(AtomData))
                    if grown == NULL:
                        raise MemoryError()
                    atom_data = grown
                a = &atom_data[n_atoms]
                a.hetatm = line[0] == b'H'
                a.atom_serial = _field_int(line, n, 6, 5)
                a.res_num = _field_int(line, n, 22, 4)
                a.x = _field_float(line, n, 30, 8)
                a.y = _field_float(line, n, 38, 8)
                a.z = _field_float(line, n, 46, 8)
                a.occupancy = _field_float(line, n, 54, 6)
                a.temp_factor = _field_float(line, n, 60, 6)
                for f in range(7):
                    a.text[f] = _field_key(line, n, _TEXT_START[f], _TEXT_WIDTH[f])
                a.model_idx = n_model
                n_atoms += 1
            elif n >= 6 and memcmp(line, b"ENDMDL", 6) == 0:
                n_model += 1

            line = nl + 1 if nl != NULL else end

        x_coords = np.empty(n_atoms, dtype=np.float32)
        y_coords = np.empty(n_atoms, dtype=np.float32)
        z_coords = np.empty(n_atoms, dtype=np.float32)
        occupancies = np.empty(n_atoms, dtype=np.float32)
        b_factors = np.empty(n_atoms, dtype=np.float32)
        atom_number = np.empty(n_atoms, dtype=np.int32)
        residue_number = np.empty(n_atoms, dtype=np.int32)
        model_idx = np.empty(n_atoms, dtype=np.int32)

        xv, yv, zv, ov, bv = x_coords, y_coords, z_coords, occupancies, b_factors
        anv, rnv, mv = atom_number, residue_number, model_idx
        hetatm = np.empty(n_atoms, dtype=np.uint8)
        hv = hetatm

        for i in range(n_atoms):
            xv[i] = atom_data[i].x
            yv[i] = atom_data[i].y
            zv[i] = atom_data[i].z
            ov[i] = atom_data[i].occupancy
            bv[i] = atom_data[i].temp_factor
            anv[i] = atom_data[i].atom_serial
            rnv[i] = atom_data[i].res_num
            mv[i] = atom_data[i].model_idx
            hv[i] = atom_data[i].hetatm

        result = {
            'record_name': np.array(['ATOM', 'HETATM'], dtype=object)[hetatm],
            'atom_number': atom_number,
            'residue_number': residue_number,
            'x_coord': x_coords,
            'y_coord': y_coords,
            'z_coord': z_coords,
            'occupancy': occupancies,
            'b_factor': b_factors,
            'model_idx': model_idx,
        }
        for f, (name, _, width) in enumerate(TEXT_FIELDS):
            result[name] = _decode_column(atom_data, n_atoms, f, width)
        return {k: result[k] for k in COLUMNS}
    finally:
        free(atom_data)


# max_atoms < 0 means no limit
def parse_pdb_file(filename, Py_ssize_t max_atoms=-1) -> object:
    with open(filename, 'rb') as f:
        data = f.read()
    return parse_pdb_bytes(data, max_atoms)


def parse_pdb_string(str input_str, Py_ssize_t max_atoms=-1) -> object:
    return parse_pdb_bytes(input_str.encode('utf-8'), max_atoms)


def parse_pdb_bytes(bytes data, Py_ssize_t max_atoms=-1) -> object:
    return _parse(data, len(data), max_atoms)

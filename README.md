# CPDB
Cython implementation of PDB -> DataFrame parsing


## Usage


### To Dictionary
```python
# To dictionary
from cpdb import parse

data = parse("path_to_pdb.pdb")
```

```
{'record_name': array(['ATOM', 'ATOM', 'ATOM', ..., 'HETATM', 'HETATM', 'HETATM'],
      dtype=object), 'atom_number': array([   1,    2,    3, ..., 1773, 1774, 1775], dtype=int32), 'atom_name': array(['N', 'CA', 'C', ..., 'O', 'O', 'O'], dtype=object), 'alt_loc': array(['', '', '', ..., '', '', ''], dtype=object), 'residue_name': array(['GLY', 'GLY', 'GLY', ..., 'HOH', 'HOH', 'HOH'], dtype=object), 'chain_id': array(['A', 'A', 'A', ..., 'A', 'A', 'A'], dtype=object), 'residue_number': array([  30,   30,   30, ..., 2276, 2277, 2278], dtype=int32), 'insertion': array(['', '', '', ..., '', '', ''], dtype=object), 'x_coord': array([31.203, 32.02 , 33.358, ..., 44.665, 41.786, 38.498], dtype=float32), 'y_coord': array([26.31 , 27.046, 26.387, ..., 13.172, 10.059, 12.491], dtype=float32), 'z_coord': array([ 6.06 ,  5.069,  4.79 , ..., 18.445, 22.316, 15.004], dtype=float32), 'occupancy': array([0.5, 0.5, 0.5, ..., 1. , 1. , 1. ], dtype=float32), 'b_factor': array([26.27, 29.29, 30.21, ..., 24.67, 34.64, 41.14], dtype=float32), 'element_symbol': array(['N', 'C', 'C', ..., 'O', 'O', 'O'], dtype=object), 'charge': array(['', '', '', ..., '', '', ''], dtype=object), 'model_idx': array([1, 1, 1, ..., 1, 1, 1], dtype=int32)}
```

### To Pandas DataFrame

```python
from cpdb import parse

data = parse("path_to_pdb.pdb", df=True)
```

```
     record_name  atom_number atom_name alt_loc residue_name chain_id  residue_number insertion    x_coord    y_coord    z_coord  occupancy   b_factor element_symbol charge  model_idx
0           ATOM            1         N                  GLY        A              30            31.202999  26.309999   6.060000       0.50  26.270000              N                 1
1           ATOM            2        CA                  GLY        A              30            32.020000  27.046000   5.069000       0.50  29.290001              C                 1
2           ATOM            3         C                  GLY        A              30            33.358002  26.386999   4.790000       0.50  30.209999              C                 1
3           ATOM            4         O                  GLY        A              30            33.810001  25.535999   5.552000       0.50  29.299999              O                 1
4           ATOM            5         N                  GLY        A              31            33.987000  26.789000   3.684000       0.50  31.889999              N                 1
...          ...          ...       ...     ...          ...      ...             ...       ...        ...        ...        ...        ...        ...            ...    ...        ...
1769      HETATM         1771         O                  HOH        A            2274            42.688999  61.925999  29.589001       1.00  39.950001              O                 1
1770      HETATM         1772         O                  HOH        A            2275            32.055000  62.648998  30.961000       0.66  15.680000              O                 1
1771      HETATM         1773         O                  HOH        A            2276            44.665001  13.172000  18.445000       1.00  24.670000              O                 1
1772      HETATM         1774         O                  HOH        A            2277            41.785999  10.059000  22.316000       1.00  34.639999              O                 1
1773      HETATM         1775         O                  HOH        A            2278            38.498001  12.491000  15.004000       1.00  41.139999              O                 1
```


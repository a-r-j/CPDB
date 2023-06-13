
import os
import pathlib

import numpy as np
import pandas as pd
import pyximport
from biopandas.pdb import PandasPdb
from pandas.testing import assert_frame_equal

pyximport.install(setup_args={
                              "include_dirs": np.get_include()},
                  reload_support=True)
from cpdb import parse_pdb

#from cpdb.parser import parse_pdb as parse

COLUMNS_TO_CHECK = ['record_name', 'atom_number', 'atom_name', 'alt_loc',
                    'residue_name', 'chain_id', 'residue_number', 'insertion',
                    'x_coord', 'y_coord', 'z_coord', 'occupancy', 'b_factor',
                    'element_symbol']

DATA_DIR = pathlib.Path(__file__).parent / "test_data"


def test():
    TEST_FILES = os.listdir(DATA_DIR)
    print(TEST_FILES)

    for f in TEST_FILES:
        print(f)
        ppdb = PandasPdb().read_pdb(str(DATA_DIR / f))
        cp = parse_pdb(str(DATA_DIR / f))
        cp = pd.DataFrame.from_dict(cp)
        cp = cp[COLUMNS_TO_CHECK]
        print(len(cp))

        atom_df = ppdb.df['ATOM'][COLUMNS_TO_CHECK].reset_index(drop=True)
        print(len(atom_df))
        print(cp.record_name.value_counts())
        atom_c_df = cp.loc[cp.record_name == "ATOM"]#.reset_index(drop=True)
        atom_c_df = atom_c_df.reset_index(drop=True)
        print(len(atom_c_df))
        print(atom_c_df["chain_id"])
        assert_frame_equal(atom_df, atom_c_df, check_dtype=False)

        hetatm_df = ppdb.df['HETATM'][COLUMNS_TO_CHECK].reset_index(drop=True)
        assert_frame_equal(hetatm_df, cp.loc[cp.record_name == 'HETATM'].reset_index(drop=True), check_dtype=False)



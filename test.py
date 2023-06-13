import os
import time

import numpy
import pandas as pd
import pyximport

pyximport.install(setup_args={
                              "include_dirs": numpy.get_include()},
                  reload_support=True)
import numpy as np
from cpdb import parse_pdb

DATA_DIR = "./tests/test_data"
TEST_FILES = os.listdir(DATA_DIR)

for f in TEST_FILES:
    print(f)
    start = time.time()
    df = parse_pdb(f"{DATA_DIR}/{f}")
    #df = pd.DataFrame.from_dict(df)
    print(df)
    end = time.time()
    #print(df)
    print(end - start)
#print(df)
#print(pd.DataFrame.from_dict(df))
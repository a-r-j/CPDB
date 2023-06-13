import os
import time

import numpy
import pyximport

pyximport.install(setup_args={
                              "include_dirs": numpy.get_include()},
                  reload_support=True)
from cpdb import parse_pdb

DATA_DIR = "./tests/test_data"
TEST_FILES = os.listdir(DATA_DIR)

for f in TEST_FILES:
    print(f)
    start = time.time()
    df = parse_pdb(f"{DATA_DIR}/{f}")
    print(df)
    end = time.time()
    print(end - start)
from typing import Dict, Union

import numpy as np
import pandas as pd

from .parser import parse_pdb


def parse(fname: str, df: bool = True) -> Union[Dict[str, np.ndarray], pd.DataFrame]:
    d = parse_pdb(fname)
    return pd.DataFrame(d) if df else d

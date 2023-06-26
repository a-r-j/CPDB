import gzip
import os
import pathlib
import sys
from typing import Dict, Optional, Union
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

import numpy as np
import pandas as pd

from .parser import parse_pdb_file, parse_pdb_string

AF2_VERSION: int = 4


def parse(
    fname: Optional[os.PathLike] = None,
    pdb_str: Optional[str] = None,
    pdb_code: Optional[str] = None,
    uniprot_id: Optional[str] = None,
    df: bool = True,
) -> Union[Dict[str, np.ndarray], pd.DataFrame]:
    if fname is not None:
        if isinstance(fname, pathlib.Path):
            fname = str(fname)
        if fname.endswith(("pdb.gz", ".ent.gz")):
            with gzip.open(fname, "rb") as f:
                pdb_str = f.read()
            pdb_str = (
                pdb_str.decode("utf-8")
                if sys.version_info[0] >= 3
                else pdb_str.encode("ascii")
            )
        else:
            d = parse_pdb_file(fname)

    if pdb_str is not None:
        if isinstance(pdb_str, list):
            pdb_str = "".join(pdb_str)
        d = parse_pdb_string(pdb_str)

    if pdb_code is not None:
        pdb_str = _fetch_pdb(pdb_code)
        d = parse_pdb_string(pdb_str)

    if uniprot_id is not None:
        pdb_str = _fetch_af2(uniprot_id)
        d = parse_pdb_string(pdb_str)

    return pd.DataFrame(d) if df else d


def _fetch_pdb(pdb_code: str) -> str:
    """Load PDB file from rcsb.org."""
    txt = None
    url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    try:
        response = urlopen(url)
        txt = response.read()
        txt = txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
    except HTTPError as e:
        print(f"HTTP Error {e.code}")
    except URLError as e:
        print(f"URL Error {e.args}")
    return txt


def _fetch_af2(uniprot_id: str) -> str:
    """Load PDB file from https://alphafold.ebi.ac.uk/."""
    txt = None
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v{AF2_VERSION}.pdb"
    try:
        response = urlopen(url)
        txt = response.read()
        txt = txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
    except HTTPError as e:
        print(f"HTTP Error {e.code}")
    except URLError as e:
        print(f"URL Error {e.args}")
    return txt

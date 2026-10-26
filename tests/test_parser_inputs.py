import gzip
import pathlib

import numpy as np

from cpdb import parse
from cpdb.parser import parse_pdb_bytes, parse_pdb_file, parse_pdb_string

DATA_DIR = pathlib.Path(__file__).parent / "test_data"


def test_entry_points_agree(tmp_path):
    f = DATA_DIR / "1crn.pdb"
    gz = tmp_path / "1crn.pdb.gz"
    gz.write_bytes(gzip.compress(f.read_bytes()))
    expected = parse_pdb_file(str(f))
    for d in (parse_pdb_string(f.read_text()), parse_pdb_bytes(f.read_bytes()), parse(str(gz), df=False)):
        assert all(np.array_equal(expected[k], d[k]) for k in expected)


def test_short_and_crlf_lines():
    d = parse_pdb_string("ATOM      1  N   MET A   1      1.000   2.000   3.000\r\nENDMDL\r\nHETATM")
    assert d["record_name"].tolist() == ["ATOM", "HETATM"]
    assert d["z_coord"].tolist() == [3.0, 0.0]
    assert d["element_symbol"].tolist() == ["", ""]
    assert d["model_idx"].tolist() == [1, 2]


# Columns 1-30, then x (8), y (8), z (8), occupancy (6), b-factor (6), 10 blank, element (2), charge (2)
ATOM_LINE = "HETATM    1 ZN    ZN A 101    " + "{x}" + "  -2.000   3.000  1.00 20.00" + " " * 10 + "ZN2+\n"


def test_two_char_charge():
    assert parse_pdb_string(ATOM_LINE.format(x="   1.000"))["charge"].tolist() == ["2+"]


def test_no_default_atom_limit():
    d = parse_pdb_string(ATOM_LINE.format(x="   1.000") * 1_000_001)
    assert len(d["x_coord"]) == 1_000_001
    assert len(parse_pdb_string(ATOM_LINE.format(x="   1.000") * 10, max_atoms=3)["x_coord"]) == 3


def test_float_fields_match_python_float():
    for x in ("-999.999", "   0.001", "  -0.000", "  12.5  ", "1.5e2   ", "        ", "   -    "):
        got = parse_pdb_string(ATOM_LINE.format(x=x))["x_coord"][0]
        expected = np.float32(float(x) if x.strip() not in ("", "-") else 0.0)
        assert got == expected and np.signbit(got) == np.signbit(expected), x


def test_many_distinct_text_values():
    # Enough distinct atom names to force the text-column hash table to grow several times
    names = [f"{i:04X}" for i in range(5000)]
    pdb = "".join(f"ATOM  {i:5d} {n} ALA A   1       1.000   2.000   3.000  1.00 20.00           C\n" for i, n in enumerate(names))
    d = parse_pdb_string(pdb * 2)
    assert d["atom_name"].tolist() == names * 2
    assert set(d["residue_name"]) == {"ALA"}

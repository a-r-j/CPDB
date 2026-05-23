import gzip
from pathlib import Path
from urllib.error import HTTPError, URLError

import cpdb
from cpdb import _fetch_af2, _fetch_pdb, parse


DATA_DIR = Path(__file__).parent / "test_data"
SAMPLE_PDB = DATA_DIR / "1crn.pdb"


class DummyResponse:
    def __init__(self, payload: bytes):
        self._payload = payload

    def read(self) -> bytes:
        return self._payload


def test_parse_accepts_pathlib_file_path():
    result = parse(fname=SAMPLE_PDB, df=False)
    assert isinstance(result, dict)
    assert len(result["atom_number"]) > 0


def test_parse_accepts_gzipped_pdb_file(tmp_path):
    gz_path = tmp_path / "copy.pdb.gz"
    gz_path.write_bytes(gzip.compress(SAMPLE_PDB.read_bytes()))

    result = parse(fname=str(gz_path), df=False)
    assert isinstance(result, dict)
    assert len(result["atom_number"]) > 0


def test_parse_accepts_list_of_lines():
    pdb_lines = SAMPLE_PDB.read_text().splitlines(keepends=True)
    result = parse(pdb_str=pdb_lines, df=False)
    assert isinstance(result, dict)
    assert len(result["atom_number"]) > 0


def test_parse_uses_pdb_code_source(monkeypatch):
    fake_dict = {"atom_number": [1], "record_name": ["ATOM"]}

    def fake_fetch(code):
        assert code == "1ABC"
        return "fake pdb contents"

    monkeypatch.setattr(cpdb, "_fetch_pdb", fake_fetch)
    monkeypatch.setattr(cpdb, "parse_pdb_string", lambda _: fake_dict)

    assert parse(pdb_code="1ABC", df=False) == fake_dict


def test_parse_uses_uniprot_source(monkeypatch):
    fake_dict = {"atom_number": [1], "record_name": ["ATOM"]}

    def fake_fetch(uniprot_id, af2_version):
        assert uniprot_id == "Q8W3K0"
        assert af2_version == 4
        return "fake af2 contents"

    monkeypatch.setattr(cpdb, "_fetch_af2", fake_fetch)
    monkeypatch.setattr(cpdb, "parse_pdb_string", lambda _: fake_dict)

    assert parse(uniprot_id="Q8W3K0", af2_version=4, df=False) == fake_dict


def test_fetch_pdb_success_and_lowercases_code(monkeypatch):
    called = {}

    def fake_urlopen(url):
        called["url"] = url
        return DummyResponse(b"pdb")

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_pdb("1ABC") == "pdb"
    assert called["url"] == "https://files.rcsb.org/download/1abc.pdb"


def test_fetch_pdb_handles_http_error(monkeypatch, capsys):
    def fake_urlopen(_):
        raise HTTPError("https://example.org", 404, "Not Found", None, None)

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_pdb("1ABC") is None
    assert "HTTP Error 404" in capsys.readouterr().out


def test_fetch_pdb_handles_url_error(monkeypatch, capsys):
    def fake_urlopen(_):
        raise URLError("network down")

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_pdb("1ABC") is None
    assert "URL Error" in capsys.readouterr().out


def test_fetch_af2_success_uses_uppercase_uniprot(monkeypatch):
    called = {}

    def fake_urlopen(url):
        called["url"] = url
        return DummyResponse(b"af2")

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_af2("q8w3k0", af2_version=3) == "af2"
    assert (
        called["url"]
        == "https://alphafold.ebi.ac.uk/files/AF-Q8W3K0-F1-model_v3.pdb"
    )


def test_fetch_af2_handles_http_error(monkeypatch, capsys):
    def fake_urlopen(_):
        raise HTTPError("https://example.org", 500, "Server Error", None, None)

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_af2("Q8W3K0") is None
    assert "HTTP Error 500" in capsys.readouterr().out


def test_fetch_af2_handles_url_error(monkeypatch, capsys):
    def fake_urlopen(_):
        raise URLError("timeout")

    monkeypatch.setattr(cpdb, "urlopen", fake_urlopen)
    assert _fetch_af2("Q8W3K0") is None
    assert "URL Error" in capsys.readouterr().out

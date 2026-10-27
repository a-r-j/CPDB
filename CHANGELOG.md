# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Source builds failing with clang (macOS) and gcc >= 14 when the build-time NumPy was older than the one
  the shipped `parser.c` was generated against (`call to undeclared function 'PyDataType_TYPEOBJ'`)

### Changed

- The extension no longer uses the NumPy C-API (the `cimport numpy` was unused), so it builds without
  NumPy headers and is not tied to the NumPy version present at build time
- The sdist no longer ships a pre-generated `parser.c`; it is regenerated from `parser.pyx` at build time
- Build requires `Cython>=3.0`; NumPy and `wheel` dropped from build requirements
- Cython is no longer a runtime dependency
- Releases publish wheels (Linux x86_64/aarch64, macOS x86_64/arm64, Windows) built with cibuildwheel
- CI tests on macOS and Windows, and against NumPy 1.x

## [0.3.0] - 2026-10-26

### Added

- [uv](https://docs.astral.sh/uv/) for dependency management and local development (`pyproject.toml`, `uv.lock`)
- GitHub Actions CI: tests on Python 3.10–3.13 and package build verification
- Automated PyPI publishing (sdist) on version tags (`v*`)
- Dependabot for GitHub Actions and uv dependencies
- This changelog
- `parse_pdb_bytes()` for parsing raw PDB bytes

### Changed

- Project metadata and dependencies consolidated in `pyproject.toml`
- README badges and development install instructions updated for uv
- Parser rewritten as a single C-level scanner: ~4x faster on large files (980k atoms: 515 ms → 125 ms)
- `max_atoms` now defaults to no limit (previously 1,000,000; larger files were silently truncated)
- Text fields are decoded with `errors="replace"` instead of raising on invalid UTF-8

### Removed

- Support for Python 3.8 and 3.9 (end of life); `requires-python` is now `>=3.10`

### Fixed

- `charge` read only column 79, truncating values such as `2+` to `2`; now reads columns 79-80
- Out-of-bounds reads on lines shorter than 80 columns; missing columns are now blank
- `\r\n` line endings
- Memory leak when the input file does not exist

## [0.2.1] - 2025-12-23

### Added

- `af2_version` argument to `parse()` to select the AlphaFold DB model version ([#6])

### Changed

- Default AlphaFold DB model version updated from v4 to v6 ([#6])

### Fixed

- AlphaFold structure downloads failing after EBI deprecated model v4 ([#6])

### Packaging

- Include `LICENSE` in published package metadata

## [0.2.0] - 2023-06-26

### Added

- Support for gzipped PDB files (`.pdb.gz`, `.ent.gz`) ([#3])
- Parse structures from RCSB by PDB code (`pdb_code`) ([#2], [#4])
- Parse structures from AlphaFold DB by UniProt ID (`uniprot_id`) ([#2], [#4])
- Parse from in-memory PDB strings or line lists (`pdb_str`) ([#2], [#4])
- Optional dict-of-arrays output via `df=False` ([#4])

## [0.1.0] - 2023-06-13

### Added

- Initial PyPI release of `cpdb-protein`
- Packaging via `setup.py` with Cython extension build
- Tests comparing parser output to [BioPandas](https://github.com/BioPandas/biopandas)

## [0.0.2] - 2023-06-13

### Fixed

- Installation and packaging issues ([#1])

## [0.0.1] - 2023-04-03

### Added

- Cython-accelerated PDB parser (`cpdb.parser`)
- High-level `parse()` API returning a pandas `DataFrame`
- Benchmark notebook and test data under `tests/test_data/`

[Unreleased]: https://github.com/a-r-j/CPDB/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/a-r-j/CPDB/releases/tag/v0.3.0
[0.2.1]: https://github.com/a-r-j/CPDB/releases/tag/v0.2.1
[0.2.0]: https://pypi.org/project/cpdb-protein/0.2.0/
[0.1.0]: https://pypi.org/project/cpdb-protein/0.1.0/
[0.0.2]: https://pypi.org/project/cpdb-protein/0.0.2/

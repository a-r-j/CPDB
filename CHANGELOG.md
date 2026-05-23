# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- [uv](https://docs.astral.sh/uv/) for dependency management and local development (`pyproject.toml`, `uv.lock`)
- GitHub Actions CI: tests on Python 3.9–3.13 and package build verification
- Automated PyPI publishing on version tags (`v*`)
- Dependabot for GitHub Actions and uv dependencies
- This changelog

### Changed

- Project metadata and dependencies consolidated in `pyproject.toml`
- README badges and development install instructions updated for uv

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

[Unreleased]: https://github.com/a-r-j/CPDB/compare/v0.2.1...HEAD
[0.2.1]: https://github.com/a-r-j/CPDB/releases/tag/v0.2.1
[0.2.0]: https://pypi.org/project/cpdb-protein/0.2.0/
[0.1.0]: https://pypi.org/project/cpdb-protein/0.1.0/
[0.0.2]: https://pypi.org/project/cpdb-protein/0.0.2/

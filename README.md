# Limits of defect-insertion bounds for slow-phonon pair order in the Holstein model

Davis Rens  
SuperMatics Inc.  
[ORCID 0000-0001-5839-5197](https://orcid.org/0000-0001-5839-5197)

Paper DOI: [10.5281/zenodo.24385569](https://doi.org/10.5281/zenodo.24385569). Supplement DOI: [10.5281/zenodo.24385573](https://doi.org/10.5281/zenodo.24385573).  
Source repository: https://github.com/SuperMatics-Inc/holstein-defect-insertion-limits

## Scope and evidence

The preprint proves an exponential estimate for centered defect-insertion observables in the spinful Holstein model under a stated low-temperature parameter regime. The estimate is uniform over finite periodic square lattices, chemical potential, insertion times, and deterministic insertion sets, with explicit restrictions on the coupling, hopping, temperature, and support sizes. The paper also identifies limits of converting this estimate into an uncentered envelope bound or a connected polymer expansion. These method-specific results do not prove or disprove a linear-in-phonon-frequency bound for pair order in the full Gibbs state. No material prediction, measurement, or claim about realized pair order is made.

The supplement replays finite rational calculations and guarded mutants. These checks do not formalize the analytic arguments. This preprint has not been peer reviewed.

## Contents and replay

- `paper/` contains the publication-metadata-filled manuscript source and its unchanged bibliography.
- `replay/` contains the dispatcher and the reviewed finite replay scripts for the three defect-insertion arguments, together with the source pins.
- `ADMITTED.sha256` records the raw-byte hashes of the complete admitted source whitelist, including the original PDF retained in release custody but not included in this repository.
- `MANIFEST.sha256` pins each file in this repository except itself.

With Python 3.12+, install the pinned replay dependencies and run both modes:

```sh
python -m pip install numpy==2.5.5 scipy==1.18.1 mpmath==1.4.1 python-flint==0.9.0
(cd replay && python -B replay.py)
(cd replay && python -B replay.py --all-mutants)
sha256sum -c MANIFEST.sha256
```

The manuscript can be built from `paper/` with Tectonic:

```sh
tectonic --keep-logs main.tex
```

The GitHub Actions workflow checks the source manifest and replay pins, runs the base and mutant modes, verifies the manuscript prose and URLs, and builds the PDF with Tectonic.

## Licenses and provenance

The manuscript text in `paper/` is distributed under CC-BY-4.0 through the paper record. The replay code and repository documentation use Apache-2.0. The Zenodo supplement is software licensed under Apache-2.0.

AI agents generated the proofs, code, and manuscript drafts under D. Rens' direction. Separate AI agents performed independent reviews.

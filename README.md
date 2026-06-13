# Partner-Midword Joint Weight Enumerators

[![Verification](https://github.com/matheorist/vscode/actions/workflows/verify.yml/badge.svg)](https://github.com/matheorist/vscode/actions/workflows/verify.yml)

This repository contains the verification scripts accompanying the manuscript

> Yong Zhang and Guangzhou Chen, "Partner-Midword Joint Weight Enumerators of
> Nearly Perfect 1-Covering Codes and Diamond Completely Regular Codes."

## Contents

- `scripts/verify_n8.py`: exhaustive length-8 checks for the structural
  identities, the canonical Type A and Type B constructions, and the associated
  length-9 diamond code.
- `scripts/verify_n16.py`: exact length-16 verification of the product
  construction, the strict-refinement witness, and the classification of the
  57 tested Vasil'ev constructions into 12 Type C joint-spectrum classes.

## Requirements

- Python 3.12 or later
- NumPy 2.4.3 for `verify_n16.py`

Install the numerical dependency with:

```console
python -m pip install -r requirements.txt
```

## Reproduction

Run the scripts from the repository root:

```console
python scripts/verify_n8.py
python scripts/verify_n16.py
```

The first script uses only the Python standard library. The second script also
uses NumPy. Both scripts terminate with `DONE.` after all assertions pass.

## Verified Results

The computations check:

- the partner involution and the endpoint-weight gap;
- the rigidity and marginal identities;
- the closed formulas for the canonical Type A and Type B constructions;
- the diamond-code reconstruction identity;
- the product trichotomy and the intersection formula for the Type I spectrum;
- a pair of length-16 NP1CCs with the same weight distribution and different
  joint spectra;
- exactly 12 Type C joint-spectrum classes among the 57 tested Vasil'ev
  constructions.

## Citation

Please cite the accompanying manuscript when using these scripts. GitHub can
also generate a software citation from `CITATION.cff`. A DOI or journal
reference can be added when the manuscript or a software release is archived.

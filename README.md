# Schoenberg validation code

Computational certificates for *Kushel–Tyaglov's problem and Sharp Schoenberg type inequalities* by Teng Zhang.

The `certificates/` directory contains the Python checkers and the allocation table used by the cubic-tail verifier. The filenames match those in the paper. The directory also includes supplementary checks and alternative filenames from the original code bundle.

## Run the cubic-tail checker

From the repository root, install the dependencies and run the short check of four allocation rows.

```sh
python -m pip install -r requirements.txt
python certificates/verify_odd_ge_three.py --smoke
```

To reconstruct all 104 polynomial identities, run the same checker without `--smoke`.

```sh
python certificates/verify_odd_ge_three.py
```

# PPARF — Plant Poly(A) Site Prediction

PPARF scans plant (e.g. *Arabidopsis*) genomic sequences and predicts **poly(A) / cleavage sites (CS)**.
Every position in a sequence is treated as a candidate site. The program scores the cis-elements
around that position and combines them into a single weighted score:

| Feature | Region (relative to candidate) | Method |
|---|---|---|
| NUE (near upstream element, AAUAAA-like) | −30 … −6 | best 6-mer PWM match |
| FUE (far upstream element) | −140 … −41 | best 10-mer PWM match |
| NUE / FUE / downstream composition | −50…−1, −150…−51, 0…+59 | pentamer probabilities |
| Cleavage site profile | −10 … +9 | 20-position PWM |
| NUE–CS distance | from NUE hit | binned distance prior |
| NUE hexamer weight | NUE motif | lookup table |

A candidate is reported if the combined score is **≥ 68.043**. At the end the program prints how many
sequences had a prediction inside the expected site window (positions 140–160) and how many had any prediction.

## Requirements

- Java JDK 8 or newer (`javac -version` to check)

## Setup

The source lives in `pparf code.txt` (package `ppa_program`). Put it where Java expects it and compile:

```bash
mkdir -p src/ppa_program
cp "pparf code.txt" src/ppa_program/Main.java
javac -d out src/ppa_program/Main.java
```

### Input files

File paths are **hardcoded** near the top of `main()` and at the `FileWriter`/`FileReader` lines.
Before compiling, edit them to match your machine:

| Variable / location | Default path | Contents |
|---|---|---|
| `nuePath` | `C:\ppa_matrices\nue_ji_final.txt` | NUE PWM (4 × 6) |
| `fueLongPath` | `C:\ppa_matrices\fue_long.res` | FUE PWM (4 × 10) |
| `cePath`, `fueShortPath`, `ndePath` | `C:\ppa_matrices\*.res` | loaded but not used in the score; they must still exist |
| `pentaFilePath` | `C:\ppa_matrices\pentamer_prob_cds.txt` | pentamer probabilities |
| `nueWeightFilePath` | `C:\ppa_matrices\weighted_hexamers_NUE_new.txt` | NUE hexamer counts/weights |
| input FASTA | `D:\Loke_8k_train.seq` | sequences to scan |
| output | `D:\pparf_output.txt` | predictions |

## File formats (with sample values)

**PWM files** (`.res`, `nue_ji_final.txt`): 3 header lines (skipped), then 4 rows (A, C, G, T)
with one space-separated column per motif position:

```
NUE matrix
rows: A C G T
--
1.20 1.85 0.10 1.90 1.95 1.70
-0.50 -1.20 -0.80 -1.10 -1.30 -0.90
-0.70 -1.10 -0.90 -1.20 -1.40 -0.60
0.30 -0.40 1.50 -0.30 -0.60 0.20
```

**Pentamer file**: 5-mer, 3 spaces, then the downstream, NUE and FUE probabilities:

```
AAAAA   0.00213 0.00412 0.00187
AAAAT   0.00198 0.00385 0.00201
```

**NUE hexamer weight file**: 6-mer, 5 spaces, count, then weight:

```
AATAAA     1520 0.8731
AATGAA     410 0.2214
```

**Input sequences** (FASTA, DNA alphabet; each sequence at least ~212 nt, with the true site near position 150):

```
>AT1G01010_polyA
ATGCTTAAATTTGTATATAATAAATTTGTTGAAATATAAAAA...
```

## Run

```bash
java -cp out ppa_program.Main
```

Sample console output:

```
 Analyzed : 8000 NUE = ... FUElong = ...
count : 6123 totalPos = 7410
Processs Completed :
```

Sample `pparf_output.txt`:

```
>AT1G01010_polyA
Pos: 153.0     Score: 71.284
Pos: 154.0     Score: 69.902
```

## Known limitations

- Paths are hardcoded Windows paths; edit them before compiling.
- Every position above the threshold is listed. Nearby hits are not merged into one peak.

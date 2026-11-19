# PPARF — Plant Poly(A) Site Prediction

PPARF scans plant (e.g. *Arabidopsis*) genomic sequences and predicts **poly(A) / cleavage sites (CS)**.
Every position in a sequence is treated as a candidate site. The program scores the cis-elements
around that position and combines them into a single weighted score:

| Feature | Region (relative to candidate) | Method |
|---|---|---|
| NUE (near upstream element, AAUAAA-like) | −30 … −6 | best 6-mer PWM match |
| FUE (far upstream element) | −140 … −41 | best 10-mer PWM match |
| NUE / FUE / downstream composition | −50…−1, −150…−51, 0…+59 | pentamer probabilities |
| Cleavage site profile | −10 … +9 | 20-position PWM (built in) |
| NUE–CS distance | from NUE hit | binned distance prior |
| NUE hexamer weight | NUE motif | lookup table |

A candidate is reported if the combined score is **≥ 68.043**. At the end the program prints how many
sequences had a prediction inside the expected site window (positions 140–160) and how many had any prediction.

## Project layout

```
src/ppa_program/
  Main.java          pipeline: reads sequences, scores each candidate position, writes predictions
  MotifScoring.java  PWM motif search, pentamer scoring, CS profile, distance features
  DataLoader.java    readers for PWM matrices, pentamer and hexamer tables
sample/              small synthetic dataset for a test run (not real model parameters)
```

## Requirements

- Java JDK 8 or newer (`javac -version` to check)

## Build

```bash
mkdir -p out
javac -d out src/ppa_program/*.java
```

## Run

```bash
java -cp out ppa_program.Main <matrixDir> <input.fasta> <output.txt>
```

| Argument | Meaning | Default if omitted |
|---|---|---|
| `matrixDir` | folder holding the model files listed below | `C:\ppa_matrices` |
| `input.fasta` | sequences to scan | `D:\Loke_8k_train.seq` |
| `output.txt` | where predictions are written | `D:\pparf_output.txt` |

### Try it with the sample data

```bash
mkdir -p out
javac -d out src/ppa_program/*.java
java -cp out ppa_program.Main sample sample/sample_sequences.fa sample_output.txt
```

Expected console output:

```
 Analyzed : 5 NUE = 453.0 FUEshort = 0.0 FUElong = 965.0 CE = 0.0 NDE = 0.0 ...
count : 4 totalPos = 4
Processs Completed :
```

Four of the five sample sequences carry an inserted poly(A) signal and get predictions. The fifth
(`sample_5_random`) is a GC-balanced control and gets none. Start of `sample_output.txt`:

```
NUE     FUElong     NuePentaScore     DownPentaScore     FuePentaScore     CS     NUECSdis
>sample_1_with_signal
Pos: 151.0     Score: 89.1934089408385
Pos: 152.0     Score: 87.21226833409938
...
>sample_5_random
```

## Model files (in `matrixDir`)

| File | Contents |
|---|---|
| `nue_ji_final.txt` | NUE PWM (4 × 6) |
| `fue_long.res` | FUE PWM (4 × 10) |
| `pentamer_prob_cds.txt` | pentamer probabilities (downstream, NUE, FUE) |
| `weighted_hexamers_NUE_new.txt` | NUE hexamer counts and weights |
| `ce_m2.res`, `fue_short.res`, `nde.res` | loaded at start-up but not used in the score |

The files in `sample/` are synthetic and only show the format. Use the real trained files for actual predictions.

### File formats

**PWM files**: 2 header lines (skipped), then 4 rows (A, C, G, T), with one single-space-separated
value per motif position and no trailing spaces:

```
NUE (synthetic sample)
rows: A C G T
1.50 1.50 -0.50 1.50 1.50 1.50
-0.50 -0.50 -0.50 -0.50 -0.50 -0.50
-0.50 -0.50 -0.50 -0.50 -0.50 -0.50
-0.50 -0.50 1.50 -0.50 -0.50 -0.50
```

**Pentamer file**: the 5-mer, exactly 3 spaces, then the downstream, NUE and FUE values separated by single spaces:

```
AAAAA   0.9000 1.0000 0.9500
AAAAC   0.7200 0.8000 0.7600
```

**NUE hexamer weight file**: the 6-mer, exactly 5 spaces, the count, one space, then the weight:

```
AATAAA     1520 0.8731
AATGAA     410 0.2214
```

**Input sequences**: FASTA, DNA alphabet (T, not U). Each sequence should be at least ~212 nt long, with the
annotated site near position 150. Scanning starts at position 151 and stops 61 nt before the end.

```
>sample_1_with_signal
CAACTAAACACCGAGTAAAAAGGAGACCCGTATGTAACTTTAATAATATTGTACTGATAG
...
```

## Known limitations

- Every position above the threshold is listed. Nearby hits are not merged into one peak.
- File parsers expect the exact spacing shown above. Extra spaces or tabs cause a `NumberFormatException`.

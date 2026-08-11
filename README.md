# GRASP–SMPS almucantar diagnostics

Small MATLAB utility for inspecting AERONET Version 3 **raw almucantar** (`.alm`) scans used in the GRASP–SMPS comparison.

The immediate purpose is to test whether close good/bad GRASP–SMPS retrieval pairs differ in the radiance scan content or almucantar branch symmetry before introducing additional case-selection thresholds.

## What the script calculates

For the GRASP wavelengths **440, 675, 870 and 1020 nm**, the script calculates for each selected scan:

- number of valid raw radiance measurements;
- counts in the manuscript angular intervals `[2°,6°]`, `[6°,30°]`, `[30°,80°]`, and `>80°`;
- an additional non-overlapping angular partition for bookkeeping;
- solar zenith angle (SZA);
- maximum **actual scattering angle** present in the scan;
- number of matched measurements on the two almucantar branches for scattering angle `>6°`;
- mean, median and maximum branch difference;
- number of branch pairs exceeding a configurable 20% reference threshold;
- scan-level summaries across the four wavelengths.

Branch difference is calculated as

```text
100 * abs(I_plus - I_minus) / ((I_plus + I_minus)/2)
```

The script also writes the individual branch-pair values and can generate a four-panel figure of branch difference versus actual scattering angle.

## Important interpretation caveat

The code analyses the **raw AERONET almucantar file**. Therefore `NValidRaw` means valid measurements present in the raw scan.

If additional points were removed during preprocessing before the GRASP `SDAT` file was written, the raw counts should not automatically be described as the exact measurements *supplied to GRASP*. To establish that, compare the result with the `SDAT` file or with the exact preprocessing code.

## Input

The script accepts either:

1. an extracted AERONET `.alm` file, or
2. the `.zip` downloaded from AERONET containing the `.alm` file.

No special MATLAB toolboxes are required.

## Usage

Download `almucantar_scan_diagnostics.m`, place it in a working directory, and edit the input path near the top:

```matlab
inputFile = '20240710_20240815_Magurele_Inoe.zip';
```

The default target cases are:

- 13 July 2024 ~14:35 UTC — good GRASP–SMPS agreement;
- 13 July 2024 ~14:54 UTC — poor agreement;
- 12 August 2024 ~14:11 UTC — poor agreement;
- 12 August 2024 ~14:46 UTC — good agreement.

The times are approximate scan-centre times. The script locates the nearest complete four-wavelength scan within the configured tolerance.

To analyse all complete scans instead, set:

```matlab
targetTimes = [];
targetLabels = strings(0,1);
```

## Output files

The script creates:

- `almucantar_wavelength_diagnostics.csv` — wavelength-level diagnostics;
- `almucantar_scan_summary.csv` — one-row-per-scan summary;
- `almucantar_branch_pairs.csv` — individual matched branch pairs;
- `almucantar_branch_symmetry.png` — optional diagnostic plot.

## Cross-check values for the July/August test file

For `20240710_20240815_Magurele_Inoe.zip`, approximate values expected from the current implementation are:

### 13 July 2024, good scan (~14:35)

- 63 valid raw measurements at every wavelength;
- maximum actual scattering angle: about 111.5–113.0°;
- wavelength mean branch differences: approximately 7.78%, 6.77%, 7.23%, 5.08% for 1020, 870, 675, 440 nm respectively;
- mean of the four wavelength means: about 6.72%;
- largest individual branch difference: about 18.85%.

### 13 July 2024, poor scan (~14:54)

- 63 valid raw measurements at every wavelength;
- maximum actual scattering angle: about 118.2–119.7°;
- wavelength mean branch differences: approximately 7.10%, 7.69%, 8.24%, 5.73% for 1020, 870, 675, 440 nm respectively;
- mean of the four wavelength means: about 7.19%;
- largest individual branch difference: about 15.63%.

### 12 August 2024

- poor scan (~14:11): mean branch difference across wavelengths about 4.30%, maximum individual difference about 11.03%;
- good scan (~14:46): mean branch difference across wavelengths about 2.43%, maximum individual difference about 11.45%.

These values are intended as implementation cross-checks, not as universal QC thresholds.

## Angular-bin detail

The manuscript intervals overlap at 6° and 30°. Therefore their counts are **not additive**. For that reason, the script reports both:

- `Ncrit_*`: counts using the exact overlapping manuscript intervals;
- `Npart_*`: a non-overlapping partition whose sum equals `NValidRaw`.

## Scientific use

The diagnostic is intended to answer a narrow question: whether differences between apparently comparable GRASP–SMPS cases can be explained by obvious changes in raw almucantar scan content, angular coverage, or branch asymmetry.

It should be used together with other diagnostics such as GRASP residuals and retrieval errors, SMPS temporal variability, MLH/ABLH representativeness, and sensitivity to initialisation/regularisation.

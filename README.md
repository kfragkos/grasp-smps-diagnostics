# GRASP–SMPS diagnostics

Small MATLAB utilities supporting the GRASP–SMPS comparison.

The repository currently contains three diagnostics:

1. `almucantar_sdat_crosscheck.m` — reproduces the photometer-side selection used when writing the GRASP SDAT file and gives a simple cross-check of the exact azimuths, counts and geometry retained after the common four-wavelength mask.
2. `almucantar_scan_diagnostics.m` — inspects AERONET Version 3 **raw almucantar** (`.alm`) scans to test branch symmetry and angular/scattering-angle diagnostics. This is a separate diagnostic and does **not** reproduce SDAT construction.
3. `grasp_avp_normalization_check.m` — inspects GRASP classic inversion-output text files to check the normalization of the printed aerosol vertical profile (AVP).

---

## 1. SDAT almucantar cross-check

`almucantar_sdat_crosscheck.m` follows the photometer preprocessing described for the manuscript SDAT construction.

### What it does

For each selected almucantar scan it:

- requires the four GRASP photometer wavelengths: 440, 675, 870 and 1020 nm;
- identifies the 30 positive/outward raw almucantar azimuth positions from 2° to 180°;
- removes the repeated 6° entry in the raw AERONET header so that the initial set is exactly 30 azimuths;
- reports the number of valid radiances at each wavelength before the common mask;
- checks angular coverage;
- applies the same common four-wavelength mask used for SDAT: if an azimuth has an invalid/NaN radiance at **any** wavelength, that azimuth is removed from **all four** wavelengths;
- reports the retained azimuths and `RAA = azimuth + 180`;
- reports the SZA at each wavelength and `VZA = 180 - SZA`;
- writes an angle-level table showing exactly which wavelength caused each discarded azimuth.

### Important point

The SDAT cross-check **does not read or use scattering angle**, because scattering angle is not required when the SDAT is written. The scattering-angle/branch-symmetry analysis remains a separate diagnostic in `almucantar_scan_diagnostics.m`.

Raw radiances are included in the detailed CSV only to trace missing values. They are **not** normalized in this script by `E0` and Earth–Sun distance, so their absolute values should not be compared directly with the normalized SDAT radiances unless the identical normalization is subsequently applied.

### Usage

Place `almucantar_sdat_crosscheck.m` with the annual AERONET raw almucantar file and edit:

```matlab
inputFile = '20240101_20241231_Magurele_Inoe.alm';
```

A ZIP containing one `.alm` file can also be used.

The current script contains the 12 GRASP–SMPS case times discussed in the analysis. The target list can be edited if the case set changes.

Run:

```matlab
almucantar_sdat_crosscheck
```

### What Mariana should compare with the SDAT

The useful columns in `almucantar_sdat_crosscheck_summary.csv` are:

- `N_SDAT_common` — should equal `n`, the number of angular measurements written for each photometer wavelength in SDAT;
- `RetainedAzimuth_deg` — exact azimuth positions remaining after the common four-wavelength mask;
- `RetainedRAA_deg` — the corresponding `azimuth + 180` values written as RAA;
- `SZA440_deg`, `SZA675_deg`, `SZA870_deg`, `SZA1020_deg` — the four wavelength-specific SZA values;
- `VZA440_deg`, etc. — `180-SZA`, for direct comparison with SDAT geometry;
- `RemovedAzimuth_deg` — positions removed because at least one wavelength was invalid.

The detailed `almucantar_sdat_crosscheck_angles.csv` contains one row for each of the original 30 azimuths and flags `Valid440`, `Valid675`, `Valid870`, `Valid1020`, and `RetainedInSDAT`. This is the easiest file to use when an SDAT count does not match.

For example, for the 13 July 2024 14:35 scan, the raw positive branch contains 30 initial azimuths and the common mask removes 2° and 2.5°, leaving 28 positions. This matches the 28 photometer radiances shown in the corresponding GRASP inversion output.

---

## 2. Raw almucantar branch-symmetry diagnostics

### What the script calculates

For the GRASP wavelengths **440, 675, 870 and 1020 nm**, `almucantar_scan_diagnostics.m` calculates for each selected scan:

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

### Important interpretation caveat

The code analyses the **raw AERONET almucantar file**, including both branches. Therefore `NValidRaw` is not the number of measurements finally written to SDAT. Use `almucantar_sdat_crosscheck.m` for that comparison.

Branch symmetry is a diagnostic of possible scan asymmetry/horizontal inhomogeneity and is **not** proposed as a universal GRASP quality-control threshold.

---

## 3. GRASP AVP normalization diagnostic

`grasp_avp_normalization_check.m` reads GRASP classic `*_inversion_output.txt` files and tests the printed aerosol vertical profile (AVP), which GRASP labels in units of `1/m`.

For each file it calculates:

- integral of the retrieved AVP points only;
- contribution from the starred ground point printed by GRASP;
- integral of the complete printed AVP;
- missing fraction `1 - integral(AVP)`;
- uppermost printed AVP value;
- the extra vertical depth associated with the missing normalized area;
- `AOD1064 * integral(printed AVP)`, for comparison with AOD obtained by integrating the printed extinction profile.

Subsequent inspection of the GRASP source code and GRASPpac literature showed that, for the LUT profile used here, GRASP adds the ground and a top-of-atmosphere point at 40 km. Above the highest retrieved lidar level, extinction is represented by a linear connection from the uppermost retrieved value to zero at 40 km; below the lowest LUT level, the lowest value is held constant to the ground. Thus an integration that stops at the highest printed retrieval level does not include the complete GRASP column.

### Usage

Place the GRASP inversion-output text files and `grasp_avp_normalization_check.m` in the same folder. By default the script processes all files matching:

```matlab
*_inversion_output.txt
```

Run:

```matlab
grasp_avp_normalization_check
```

The script prints a summary table and writes:

```text
grasp_avp_normalization_summary.csv
```

---

## Scientific use

These scripts answer narrow diagnostic questions. They should be interpreted together with GRASP residuals and retrieval errors, SMPS averaging/uncertainty, boundary-layer representativeness, and targeted inversion-stability tests rather than used as universal QC filters.

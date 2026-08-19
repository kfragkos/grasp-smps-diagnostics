# GRASP–SMPS diagnostics

MATLAB utilities supporting diagnostics for GRASP retrievals that combine sun/sky photometer and ceilometer observations, with particular emphasis on GRASP–SMPS comparisons.

The repository currently contains three complementary tools:

1. `almucantar_sdat_crosscheck.m` — reproduces the photometer-side angular selection used when constructing a GRASP SDAT file and reports the retained azimuths, counts and geometry.
2. `almucantar_scan_diagnostics.m` — examines AERONET Version 3 raw almucantar (`.alm`) scans for branch-symmetry and scattering-angle diagnostics. This analysis is independent of SDAT construction.
3. `grasp_avp_normalization_check.m` — checks how the aerosol vertical profile (AVP) printed in GRASP classic inversion output relates to the full profile normalization used by GRASP.

---

## 1. SDAT almucantar cross-check

`almucantar_sdat_crosscheck.m` reproduces the photometer preprocessing used to select angular radiances before writing them to an SDAT file.

### What it does

For each selected almucantar scan it:

- requires measurements at 440, 675, 870 and 1020 nm;
- identifies the 30 positive/outward almucantar azimuth positions from 2° to 180°;
- removes the repeated 6° entry in the raw AERONET header so that the initial set contains 30 azimuths;
- reports the number of valid radiances at each wavelength before the common mask;
- checks angular-range coverage;
- applies a common four-wavelength mask: if an azimuth has an invalid/NaN radiance at any wavelength, that azimuth is removed from all four wavelengths;
- reports the retained azimuths and `RAA = azimuth + 180`;
- reports the wavelength-specific SZA values and `VZA = 180 - SZA`;
- writes an angle-level table showing which wavelength caused each discarded azimuth.

### Scattering angle

The SDAT cross-check does **not** read or use scattering angle, because scattering angle is not required for this SDAT construction step. Scattering-angle and branch-symmetry diagnostics are handled separately by `almucantar_scan_diagnostics.m`.

Raw radiances are included in the detailed CSV only to trace missing values. They are not normalized in this script by extraterrestrial irradiance (`E0`) and Earth–Sun distance, so their absolute values should only be compared with SDAT radiances after applying the same normalization.

### Usage

Place `almucantar_sdat_crosscheck.m` with an annual AERONET raw almucantar file and edit:

```matlab
inputFile = '20240101_20241231_Magurele_Inoe.alm';
```

A ZIP containing one `.alm` file can also be used.

The script contains an example list of selected retrieval times from the study. This list can be edited, or all complete scans can be analysed by setting:

```matlab
targetTimes = [];
targetLabels = strings(0,1);
```

Run:

```matlab
almucantar_sdat_crosscheck
```

### Output

The script writes:

- `almucantar_sdat_crosscheck_summary.csv`
- `almucantar_sdat_crosscheck_angles.csv`

Useful summary fields include:

- `N_SDAT_common` — number of azimuths retained after the common four-wavelength mask;
- `RetainedAzimuth_deg` — retained azimuth positions;
- `RetainedRAA_deg` — corresponding `azimuth + 180` values;
- `SZA440_deg`, `SZA675_deg`, `SZA870_deg`, `SZA1020_deg` — wavelength-specific SZA values;
- `VZA440_deg`, etc. — `180-SZA`;
- `RemovedAzimuth_deg` — azimuths removed because at least one wavelength was invalid.

The angle-level file contains one row for each initial azimuth and flags `Valid440`, `Valid675`, `Valid870`, `Valid1020`, and `RetainedInSDAT`. This allows any difference in retained angular counts to be traced to a specific wavelength and azimuth.

---

## 2. Raw almucantar branch-symmetry diagnostics

`almucantar_scan_diagnostics.m` is intentionally separate from the SDAT cross-check. It uses the full raw AERONET scan, including both almucantar branches and the scattering-angle fields, to test possible scan asymmetry or horizontal inhomogeneity.

It calculates, among other quantities:

- raw valid radiance counts;
- angular coverage;
- maximum actual scattering angle;
- matched plus/minus branch pairs;
- mean, median and maximum branch difference.

Branch difference is calculated as

```text
100 * abs(I_plus - I_minus) / ((I_plus + I_minus)/2)
```

The 20% line used in the diagnostic plots is a reference value only and is not proposed as a universal GRASP quality-control threshold.

`NValidRaw` refers to the full raw scan and should not be interpreted as the number of angular measurements written to an SDAT file. For the latter, use `N_SDAT_common` from `almucantar_sdat_crosscheck.m`.

---

## 3. GRASP AVP normalization diagnostic

`grasp_avp_normalization_check.m` reads GRASP classic `*_inversion_output.txt` files and evaluates the printed aerosol vertical profile (AVP), reported in units of `1/m`.

For a LUT aerosol profile, GRASP adds boundary points at the ground and at the top of the model atmosphere. The lower boundary is assigned the value of the lowest LUT level. Above the highest retrieved LUT altitude, the profile is connected linearly to zero at the model top. For the configuration considered here, the model top is 40 km.

The script calculates:

- the integral of the retrieved AVP points;
- the contribution from the additional ground point;
- the integral of the complete printed AVP up to the highest retrieved altitude;
- the missing normalized fraction;
- the upper-triangle contribution expected from the linear connection to zero at 40 km;
- the model-top altitude inferred from the printed AVP and its missing normalization;
- closure of the complete profile normalization;
- the corresponding 1064-nm AOD contributions.

This makes explicit why integrating only the printed profile up to the highest retrieved altitude yields a value below the full GRASP AOD: the upper linear extension to 40 km is not included in that truncated integration.

---

## Scientific use

These scripts address specific diagnostic questions and are not intended as universal quality-control filters. Their outputs should be interpreted together with retrieval residuals and uncertainties, measurement representativeness, and sensitivity tests appropriate to the application.

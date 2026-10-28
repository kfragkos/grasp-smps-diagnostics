# GRASP–SMPS diagnostics

MATLAB utilities supporting diagnostics for GRASP retrievals that combine sun/sky photometer and ceilometer observations, with particular emphasis on GRASP–SMPS comparisons.

The repository contains four complementary tools:

1. `almucantar_sdat_crosscheck.m` — reproduces the photometer-side angular selection used when constructing a GRASP SDAT file and reports the retained azimuths, counts and geometry.
2. `almucantar_scan_diagnostics.m` — examines AERONET Version 3 raw almucantar (`.alm`) scans for branch-symmetry and scattering-angle diagnostics. This analysis is independent of SDAT construction.
3. `grasp_avp_normalization_check.m` — checks how the aerosol vertical profile (AVP) printed in GRASP classic inversion output relates to the full profile normalization used by GRASP.
4. `almucantar_pair_comparison.m` — directly pairs both AERONET almucantar branches, applies the 20% branch-screening plus a common four-wavelength mask, and compares specifically requested case pairs.

---

## 1. SDAT almucantar cross-check

`almucantar_sdat_crosscheck.m` checks the positive-branch photometer angles before writing an SDAT file, **but its common missing-value mask does not implement the additional screening of both almucantar branches**. Therefore, its `N_SDAT_common` is a preliminary count, **not necessarily the actual SDAT count**. Use the new `almucantar_pair_comparison.m` for a cross-check that includes the 20% branch-difference screening.

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


---

## 4. Selected 2024 almucantar case-pair comparisons

The self-contained MATLAB script `almucantar_pair_comparison.m` compares two GRASP–SMPS pairs from **26 May 2024** (14:26:03, 16:00:28 UTC) and **23 September 2024** (07:19:17, 13:40:58 UTC).

It reads the full AERONET raw `.alm` scan, identifies the 30 outward azimuth angles (+2° to +180°, with one duplicated +6° entry removed), pairs each with the corresponding return-branch measurement, and calculates

\[
D = 200\,|I_+-I_-|/(I_++I_-)\quad [\%].
\]

An angular position is rejected if **either** branch is missing or nonpositive, or if **any one** of the four wavelengths (440, 675, 870, 1020 nm) has `D > 20%`. The resulting common angular mask is applied to all four wavelengths. This implements the branch comparison described for the photometer preprocessing; it is **not** a new general-purpose retrieval-quality threshold. The multiplicative radiance normalization does not affect `D`.

The script reports SZA, raw and retained maximum scattering angle, mean/maximum branch differences for **actual scattering angles >6°**, number of retained angles, and removed RAA values. Scattering-angle fill values (e.g., 999°) are excluded from these diagnostics. A 28 October 2026 code correction aligned the threshold with actual scattering angle rather than nominal azimuth; recompute results produced using earlier revisions before making exact numerical comparisons. **Scattering angle is a diagnostic only and does not enter the SDAT-screening mask.**

### Run

Raw almucantar data can be downloaded using the official [AERONET Version 3 raw-sky data service](https://aeronet.gsfc.nasa.gov/print_web_data_help_v3_raw_sky.html) (`ALM00=1`, `AVG=10`). For example, request [26 May 2024](https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_raw_sky_v3?site=Magurele_Inoe&year=2024&month=5&day=26&year2=2024&month2=5&day2=27&ALM00=1&AVG=10&if_no_html=1) or [23 September 2024](https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_raw_sky_v3?site=Magurele_Inoe&year=2024&month=9&day=23&year2=2024&month2=9&day2=24&ALM00=1&AVG=10&if_no_html=1), save the plain-text response with a `.alm` extension, and run the script on each day's file. The AERONET endpoint was not reachable from this analysis environment, so these day-specific downloads remain unverified here.

Download an annual AERONET Version 3 raw almucantar `.alm` file for Magurele_Inoe and place it alongside the MATLAB script. Then execute:

```matlab
[summary, angles] = almucantar_pair_comparison( ...
    '20240101_20241231_Magurele_Inoe.alm');
```

A ZIP containing a single `.alm` is also supported. Shorter time ranges may be used, in which case out-of-range requested scans are explicitly marked `MISSING`. Only scans matching the requested times within five minutes are accepted; a nearby scan (e.g. 07:38 instead of 07:19) must not be substituted.

The script writes `results/almucantar_case_pair_summary.csv`, `results/almucantar_case_pair_angles.csv`, and a per-day four-wavelength PNG for each date with matched scans. The detailed CSV gives all per-wavelength branch differences, the common mask, and retained RAA values for direct comparison with the independent GRASP SDAT. These output files are generated locally when running MATLAB; they are not automatically uploaded to GitHub.

### Available results

See [results/2024_case_pair_diagnostics.csv](results/2024_case_pair_diagnostics.csv). This file combines the selected 2024 GRASP–SMPS statistics with the **archived** earlier raw almucantar diagnostic results for the May pair. The two exact September raw scans have **not** been analyzed, because their raw `.alm` data were not available. The earlier September 07:38 scan is **not** the 07:19 retrieval and has not been substituted. Empty fields in the CSV mean unavailable measurements, **not zero**.

| Case (UTC) | VSD comparison | VSD NRMSB | Max raw scattering angle | Mean raw branch difference (>6°) | Mean after reported angular exclusion (>6°) |
|---|---|---:|---:|---:|---:|
| 26 May 14:26 | Satisfactory | 0.619 | 113.46° | 2.06% | 2.06% |
| 26 May 16:00 | Good | 0.182 | 146.88° | 3.04% | 2.22% |
| 23 Sep 07:19 | Good | 0.138 | not available | not available | not available |
| 23 Sep 13:40 | Satisfactory | 0.390 | not available | not available | not available |

For May 16:00, the raw 440-nm branch difference reaches 47.05% at nominal azimuth 180° (RAA 360°). That azimuth was reported as excluded from SDAT. Excluding it at all four wavelengths lowers the overall maximum retained branch difference (for paired angles >6°) to **14.31%**. This is consistent with the absence of a marked branch-asymmetry problem in the retained measurements of either May case. The better May VSD match has larger scattering-angle coverage; this association alone does **not** establish a causal effect of angular coverage on VSD retrieval accuracy.

### Cross-check of the 20% implementation

The branch-screening/common-mask procedure was independently checked with the available July–August 2024 raw AERONET archive (not included in this repository). It reproduced reported retained-angle counts and several reported removed RAA values: 14 July 14:35 (25 angles; RAA 250°, 260°, 360° plus two missing low-angle positions); 15 July 14:54 (17); 8 August 16:08 (25; RAA 196°, 320°, 360° plus two missing low-angle positions); and 11 August 14:12 (23; RAA 240°, 250°, 260°, 270°, 300° plus two missing low-angle positions). The code has not been executed in a MATLAB runtime here; these cross-checks used an independent implementation of the same numerical rules.

**Interpretation:** This is a limited diagnostic of the radiance input, not a criterion for selecting good GRASP–SMPS VSD comparisons. The contrast in VSD agreement cannot be straightforwardly attributed to almucantar branch symmetry, and reduced angular coverage remains an unverified possible contributor.

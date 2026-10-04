# WellLog Pro: User Guide (v1.0)

WellLog Pro is a Windows desktop app for quick-look well-log analysis. Load a LAS or CSV file, check data quality, view the logs, and calculate shale volume, porosity, water saturation, pay flags and net-to-gross, all with visible, editable parameters.

> **Engineering disclaimer:** Results are engineering estimates. They depend on data quality, environmental corrections, fluid properties, model assumptions and your parameter choices. Always review the results. They are not a substitute for professional interpretation.

---

## 1. Installation

**Option A: ready-made app.** Double-click `WellLogPro.exe`. Nothing else is needed.
If Windows shows *"Windows protected your PC"*, click **More info → Run anyway**. This is normal for unsigned software.

**Option B: build it yourself.**
1. Install Python 3.9+ from python.org and tick **"Add python.exe to PATH"**.
2. Put `app.py`, `petro.py`, `build.bat` and `requirements.txt` in one folder.
3. Double-click `build.bat`. When it finishes, your app is `dist\WellLogPro.exe`.

To run from source instead: `pip install numpy matplotlib` then `python app.py`.

---

## 2. Quick start (5 steps)

1. Click **Load LAS / CSV...** and choose a file (try `sample_synthetic_well_log.csv`).
2. Open the **QC** tab and read the data quality report.
3. Review the parameters on the left (especially **GR min/max**, **Rw**, **matrix density**).
4. Click **CALCULATE**.
5. Inspect the **Log Viewer**, **Crossplots** and **Zones / Summary** tabs, then export.

---

## 3. Input files

### LAS (version 2.0)
- The app reads the well name, curve names, units, depth and the null value (usually -999.25). Null values are treated as missing.
- The first curve must be depth.

### CSV
- First row = column headers. One column per curve.
- Depth column should be named `DEPT`, `DEPTH` or `MD`. Otherwise the first column is used as depth.
- Empty cells and non-numeric values are treated as missing.
- CSV files carry no unit information, so the depth axis is labelled "m". Interpret it according to your data.

### Automatic curve recognition

| Curve | Recognised names |
|---|---|
| Gamma ray (GR) | GR, GAM, GRC, CGR, SGR |
| Deep resistivity (RT) | RT, LLD, ILD, AT90, RDEP, RD, RLLD, RILD |
| Shallow resistivity (RS) | LLS, ILM, RS, AT30, RLLS |
| Micro/flushed zone (RXO) | MSFL, RXO, SFL, RMSFL |
| Bulk density (RHOB) | RHOB, DEN, RHOZ, RHOM, DENS |
| Neutron porosity (NPHI) | NPHI, TNPH, NEU, NPOR, CNL |
| Sonic (DT) | DT, DTC, AC, DTCO, DT4P |

The mapping used is printed at the bottom of the **QC** tab. If your curve has a different name, rename it in the file (for example to `GR` or `RHOB`) and reload. The app never assumes a curve that isn't there.

**Neutron units:** if NPHI values go above 1.5 the app treats them as percent and divides by 100.

---

## 4. Parameters panel

Every value is editable. Defaults are common starting points, not recommendations for your formation.

| Parameter | Default | Meaning |
|---|---|---|
| GR min / GR max | auto | Clean-sand and shale GR baselines. **Auto-filled on load** from the 5th and 95th percentiles. Check and adjust. |
| Matrix density | 2.65 | g/cc. Presets: Sandstone 2.65, Limestone 2.71, Dolomite 2.87 |
| Fluid density | 1.0 | g/cc |
| Rw | 0.08 | Formation water resistivity (ohm-m) |
| Rsh | 2.0 | Shale resistivity (ohm-m), used by Simandoux and Indonesia |
| Archie a, m, n | 1, 2, 2 | Tortuosity, cementation, saturation exponent |
| Vsh cutoff | 0.35 | Pay requires Vsh **below** this |
| Porosity cutoff | 0.10 | Pay requires PHI **above** this |
| Sw cutoff | 0.60 | Pay requires Sw **below** this |

Use the dropdowns to choose the **Vsh method**, **Porosity** method and **Sw model**.

---

## 5. Methods and equations

**Shale volume.** First the Gamma Ray Index: `IGR = (GR - GRmin) / (GRmax - GRmin)`, limited to 0–1.
- Linear: `Vsh = IGR`
- Larionov (tertiary): `Vsh = 0.083 (2^(3.7 IGR) - 1)`
- Larionov (older): `Vsh = 0.33 (2^(2 IGR) - 1)`
- Clavier: `Vsh = 1.7 - sqrt(3.38 - (IGR + 0.7)^2)`
- Steiber: `Vsh = IGR / (3 - 2 IGR)`

**Porosity.**
- Density: `PHI = (rho_ma - rho_b) / (rho_ma - rho_fl)`, limited to 0–0.6
- Neutron: NPHI from the log (v/v), limited to 0–0.6
- Density-Neutron average: `(PHI_D + PHI_N) / 2`

**Water saturation.** Sh = 1 - Sw. Results are limited to 0–1.
- Archie: `Sw = ((a Rw) / (Rt PHI^m))^(1/n)`
- Simandoux (modified): `Sw = (a Rw / (2 PHI^m)) [ sqrt((Vsh/Rsh)^2 + 4 PHI^m / (a Rw Rt)) - Vsh/Rsh ]`. **This form assumes n = 2, and the n box is not used.**
- Indonesia: `1/sqrt(Rt) = [ Vsh^(1 - Vsh/2) / sqrt(Rsh) + PHI^(m/2) / sqrt(a Rw) ] Sw^(n/2)`

**Pay flag.** A sample is pay when `Vsh < cutoff AND PHI > cutoff AND Sw < cutoff`. Samples with missing data are never flagged as pay.

**Zone statistics.** Gross = base - top. Net pay = number of pay samples × median depth step. NTG = Net / Gross. Mean, minimum and maximum are given for Vsh, PHI and Sw.

---

## 6. The tabs

**Log Viewer.** Tracks appear only for curves that exist: GR, Resistivity (log scale), RHOB/NPHI overlay, DT, and after CALCULATE a calculated track (Vsh, PHI, Sw with gold pay shading). Move the mouse over any track and a depth line follows across all tracks. Use the toolbar to zoom, pan, reset or save.

**Crossplots.** Left: density-neutron plot coloured by GR. Right: Pickett plot (Rt vs PHI, log-log, coloured by depth) with dashed Sw = 1, 0.75, 0.5 and 0.25 lines built from your a, m, n and Rw.

**QC.** Depth coverage, sample step, duplicate or non-increasing depths, percent valid per curve, and samples outside plausible ranges. Problems are **flagged, never deleted**.

**Zones / Summary.** The parameters used, plus statistics for each zone you defined.

---

## 7. Defining zones

In the **Zones** box on the left, enter one zone per line as `Name,top,base`:

```
Sand A,2458,2487
Sand B,2521,2545
```

On load, one zone covering the whole interval is added for you. Click **CALCULATE** again after editing zones.

---

## 8. Exporting

Run **CALCULATE** first.
- **Export CSV:** depth, all input curves, VSH, PHI, SW, SH and PAY (1/0). The first line is a comment containing every parameter used, so the result can be reproduced.
- **Export PNG:** saves the figure on the current tab (Log Viewer or Crossplots).
- **Export PDF report:** the QC and summary text plus the log and crossplot figures.

Original files are never modified.

---

## 9. Troubleshooting

| Problem | Fix |
|---|---|
| "Required curves not found" | The needed curve wasn't recognised. Rename it in the file (see section 3). Required: GR and deep resistivity, plus RHOB (density/average) or NPHI (neutron/average). |
| Everything shows Sw = 1 | Check that Rw is realistic, resistivity is in ohm-m, and porosity looks sensible. |
| Vsh is all 0 or all 1 | Adjust GR min and max. |
| Porosity looks wrong | Check matrix density, and whether NPHI is in v/v or percent. |
| "Load failed" | The LAS may be unusual (for example version 3.0, or wrapped data). Save it as LAS 2.0 or CSV. |
| Windows or antivirus warning | The app is unsigned. Use More info → Run anyway, or add an exception. |

---

## 10. Current limitations (v1.0)

- Inputs: LAS 2.0 and CSV only (no Excel or LAS 3.0).
- One well at a time. No multi-well correlation.
- No permeability estimate.
- No environmental or borehole corrections, and no uncertainty estimates.
- Pay flags are criteria-based screening, not a definitive reservoir call.
- Windows only.

For support, include your Windows version, the error message and, if possible, the first 30 lines of your input file.

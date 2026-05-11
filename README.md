# 📋 Project 8 — GNSS Data Quality Report

> **Completeness · SNR · Cycle Slips · Multipath · Gap Map · Auto-Generated Report | Auckland, NZ**

---

## 📌 Overview

Before any GNSS data is used for positioning, atmospheric research, or geodetic
applications, a **data quality check is mandatory**. This is the first thing a
GNSS data engineer does when receiving a new observation file.

This project produces a complete automated quality assessment of a RINEX 3
observation file — covering every satellite from every constellation.

| Metric | What It Measures |
|--------|-----------------|
| **Coverage** | % of total epochs with valid pseudorange per satellite |
| **Mean SNR** | Average signal strength [dB-Hz] |
| **Low-SNR fraction** | % of epochs below 30 dB-Hz |
| **Cycle slips** | Phase discontinuities via L4 geometry-free combination |
| **o/slps ratio** | Observations per cycle slip (higher = cleaner tracking) |
| **MP1 RMS** | Multipath proxy on L1 code [metres] |

---

## 🖼️ Output Files

### Plot 1 — 4-Panel Quality Dashboard
Four stacked bar charts for every satellite:
- Coverage [%] with 80% threshold line
- Mean SNR [dB-Hz] with 25 and 35 dB-Hz reference lines
- Cycle slips count (green < 2, orange 2–5, red > 5)
- MP1 RMS multipath proxy [m] (green < 0.3, orange 0.3–0.5, red > 0.5)

### Plot 2 — Data Gap Map
Binary heatmap: green = satellite tracked, dark = no observation.
Shows exactly when each satellite was visible and tracked over 24 hours.

### Plot 3 — Quality Scatter (Coverage vs SNR)
Each satellite as a dot, coloured by constellation, sized by cycle slip count.
Top-right = best satellites. Enables instant quality ranking.

### quality_report.txt — Auto-Generated Text Report
Plain-text file containing:
- File metadata (station, date, sampling rate)
- Constellation summary with per-system averages
- File-level quality flags (⚠ warnings / ✅ pass)
- Per-satellite detail table

---

## 📂 File Structure

```
project7-data-quality-report/
├── project7_data_quality_report.ipynb   ← Main notebook
├── requirements.txt                      ← Python dependencies
├── LICENSE                               ← MIT License
└── README.md                             ← This file
```

---

## ⚙️ How to Run

### 1. Install dependencies

```bash
pip install -r requirements.txt
```

### 2. Set your RINEX file path

Update **Step 2** of the notebook:

```python
obs_path = "/path/to/your/file.rnx"
```

### 3. Run all cells

```bash
jupyter notebook project7_data_quality_report.ipynb
```

Outputs: `plot1_quality_dashboard.png`, `plot2_gap_map.png`,
`plot3_quality_scatter.png`, and `quality_report.txt`

---

## 📐 Quality Metrics — Formulas

### Cycle slip detection (L4 method):
```
L4(t) = Φ₁(t) − Φ₂(t)     [metres]
ΔL4   = L4(t) − L4(t−1)    [per epoch]
Slip if |ΔL4| > 0.15 m
```

### Multipath proxy MP1:
```
MP1 = C1C − Φ₁ − 2·(f₂²/(f₁²−f₂²))·(Φ₁ − Φ₂)
    ≈ multipath on L1 + noise
```
- RMS(MP1) < 0.3 m → clean environment
- RMS(MP1) > 0.5 m → significant multipath (buildings, trees, water)

### Observations per slip:
```
o/slps = n_valid_epochs / (n_slips + 1)
```
Higher = fewer interruptions per observation → better data continuity.

---

## 🛠️ Dependencies

| Package | Purpose |
|---------|---------|
| `georinex` | Parse RINEX 3 observation files |
| `xarray` | N-dimensional labelled arrays |
| `pandas` | Time series and DataFrame operations |
| `numpy` | Numerical computations |
| `matplotlib` | Publication-quality plotting |

---

## 👤 Author

**Hakim El Azzouzi**  
MSc Global Navigation Satellite Systems  
Mohammed First University, Oujda, Morocco  
📧 elazzouzihakim10@gmail.com  
🔗 [linkedin.com/in/Hakim-El-Azzouzi](https://linkedin.com/in/Hakim-El-Azzouzi)  
📍 Luxembourg 🇱🇺

---

## 📜 License

MIT License — see [LICENSE](LICENSE) for details.

---

## 🔗 GNSS RINEX Analysis Series — Complete

| # | Project |
|---|---------|
| 1 | Single GPS Satellite — Pseudorange & SNR Heatmap |
| 2 | All GPS Satellites — Fleet Pseudorange & SNR Heatmap |
| 3 | Multi-Constellation GNSS — One Satellite per System |
| 4 | Pseudorange vs Carrier-Phase Comparison |
| 5 | Constellation Summary — Pie Chart & Histograms |
| 6 | Ionospheric Delay — Geometry-Free Combination |
| 7 | Multipath Analysis |
| **8** | **Data Quality Report** ← You are here |

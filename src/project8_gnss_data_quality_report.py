=============================================================================
🛰️ Project 8 — GNSS Data Quality Report
=============================================================================
-----------------------------------------------------------------------------
 Station  : AUCK00NZL  —  Auckland, New Zealand  (CDDIS / LINZ Network)
 File     : AUCK00NZL_R_20260010000_01D_30S_MO.rnx
 Receiver : TRIMBLE ALLOY
 Antenna  : TRM115000.00
 Format   : RINEX 3.05  |  Mixed Constellations  |  30-second sampling
 Date     : 2026-01-01  (Day-of-Year 001)
-----------------------------------------------------------------------------
 Description
 -----------
Before any GNSS data is used for positioning, atmospheric research, or geodetic
applications, a **data quality check** must be performed.

This project produces a **complete automated data quality report** from a RINEX 3.05
observation file, covering:

| Check | What It Detects |
|-------|-----------------|
| Observation completeness | Missing epochs, data gaps, coverage per satellite |
| Signal quality | SNR distribution, low-SNR epochs per satellite |
| Cycle slip detection | Phase discontinuities via L4 and Melbourne-Wübbena |
| Code-phase consistency | Multipath proxy MP1 per satellite |
| Constellation health | Availability, redundancy, DOP proxy |
| Text report | Auto-generated summary saved as `quality_report.txt` |

---
##  Quality Metrics Used

### 1. Observation completeness
```
Coverage [%] = (valid epochs / total epochs) × 100
```

### 2. Cycle slip detection — L4 method
```
L4(t) = Φ₁(t) − Φ₂(t)     [metres]
ΔL4(t) = L4(t) − L4(t−1)   [metres per epoch]

Cycle slip if |ΔL4| > threshold  (typically 0.15–0.20 m)
```

### 3. Multipath proxy — MP1
```
MP1 = C1C − Φ₁ − 2·(f₂²/(f₁²−f₂²))·(Φ₁−Φ₂)
    ≈ multipath on L1 code + noise
```
RMS(MP1) < 0.3 m → good environment
RMS(MP1) > 0.5 m → significant multipath

### 4. o/slps ratio (observations per slip)
```
o/slps = valid epochs / number of cycle slips
```
Higher = better (fewer interruptions relative to observations)
---
-----------------------------------------------------------------------------
 **About the projects**
 ----------------------
# Step1: Install & Import Libraries
# Step2: Load the RINEX File
# Step3: Compute Quality Metrics for Every Satellite
# Step4: Plot 1: Coverage & SNR Dashboard
# Step5: Plot 2: Data Gap Map
# Step6: Plot 3: Quality Score Summary Scatter
# Step7: Generate the Text Quality Report
=============================================================================
# ───────────────────────────────────
# Step 1 — Install & Import Libraries
# ───────────────────────────────────
# !pip install --upgrade georinex

import georinex as gr
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
import warnings
from datetime import datetime
import os
os.makedirs('../content/output', exist_ok=True)

warnings.filterwarnings('ignore')
plt.rcParams.update({'figure.dpi': 120})

# Physical constants
C      = 299_792_458.0
F1_GPS = 1_575.42e6
F2_GPS = 1_227.60e6
LAM1   = C / F1_GPS   # ≈ 0.19029 m
LAM2   = C / F2_GPS   # ≈ 0.24421 m

# Multipath coefficient for MP1
# MP1 = C1 - L1 - 2*(f2²/(f1²-f2²))*(L1-L2)
ALPHA_MP = 2 * F2_GPS**2 / (F1_GPS**2 - F2_GPS**2)   # ≈ 1.5457

# Cycle slip L4 threshold (metres)
CS_THRESHOLD = 0.15

# Constellation colours
CONST_COLORS = {
    'G': '#2196F3', 'R': '#F44336',
    'E': '#4CAF50', 'C': '#FF9800', 'J': '#9C27B0', 'S': '#00BCD4',
}
CONST_NAMES = {
    'G': 'GPS', 'R': 'GLONASS', 'E': 'Galileo',
    'C': 'BeiDou', 'J': 'QZSS', 'S': 'SBAS',
}

print('✅ Libraries loaded')
print(f'   Cycle slip L4 threshold : {CS_THRESHOLD} m')
print(f'   MP1 alpha coefficient   : {ALPHA_MP:.4f}')

# ─────────────────────────────────────────────────────────────────
# Step 2 — Load the RINEX File
# ─────────────────────────────────────────────────────────────────
# RINEX FILE PATH HERE
obs_path = "/AUCK00NZL_R_20260010000_01D_30S_MO.rnx"  # ← change this path
# Read the file header first (fast — no data loaded yet)
print("📋 FILE HEADER")
print("=" * 60)
header = gr.rinexheader(obs_path)

for key, value in header.items():
    print(f"{key:<25}: {value}")

print()

# Load all observation data (interval=30 means keep 30-sec rate)

print("⏳ Loading observation data (this may take 1–2 minutes)...")
obs = gr.load(obs_path, interval=30)
print()
print("✅ Data loaded!")
print(obs)

# ────────────────────────────────────────────────────
# Step 3 — Compute Quality Metrics for Every Satellite
# ────────────────────────────────────────────────────
all_sv       = obs.sv.values
time_index   = pd.to_datetime(obs.time.values)
total_epochs = len(time_index)

# Storage for quality metrics

results = []   # list of dicts, one per satellite

print("⏳ Computing quality metrics for all satellites...")
print()

for sat in sorted(all_sv):
    prefix = sat[0]
    rec = {'sat': sat, 'system': CONST_NAMES.get(prefix, prefix)}

    # ── Observation completeness (pseudorange) ──────────────
    pr_code = 'C1C' if 'C1C' in obs.data_vars else (
               'C1X' if 'C1X' in obs.data_vars else None)
    if pr_code:
        pr = obs[pr_code].sel(sv=sat).to_series()
        rec['n_obs']    = pr.notna().sum()
        rec['coverage'] = rec['n_obs'] / total_epochs * 100
    else:
        rec['n_obs']    = 0
        rec['coverage'] = 0.0

    # ── Mean SNR ────────────────────────────────────────────
    snr_val = np.nan
    for snr_code in ['S1C', 'S1X', 'S1P']:
        if snr_code in obs.data_vars:
            s = obs[snr_code].sel(sv=sat).to_series().dropna()
            if len(s) > 0:
                snr_val = s.mean()
                break
    rec['mean_snr'] = snr_val

    # ── Low-SNR fraction ────────────────────────────────────
    if not np.isnan(snr_val) and snr_code in obs.data_vars:
        s_all = obs[snr_code].sel(sv=sat).to_series().dropna()
        rec['pct_low_snr'] = (s_all < 30).sum() / len(s_all) * 100
    else:
        rec['pct_low_snr'] = np.nan

    # ── Cycle slip detection via L4 ─────────────────────────
    n_slips = 0
    if prefix == 'G' and 'L1C' in obs.data_vars and 'L2W' in obs.data_vars:
        phi1 = obs['L1C'].sel(sv=sat).to_series().dropna() * LAM1
        phi2 = obs['L2W'].sel(sv=sat).to_series().dropna() * LAM2
        common = phi1.index.intersection(phi2.index)
        if len(common) > 5:
            L4 = phi1[common] - phi2[common]
            dL4 = L4.diff().dropna()
            n_slips = (np.abs(dL4) > CS_THRESHOLD).sum()
    elif prefix == 'E' and 'L1X' in obs.data_vars and 'L5X' in obs.data_vars:
        lam1_e = C / 1_575.42e6
        lam5_e = C / 1_176.45e6
        phi1 = obs['L1X'].sel(sv=sat).to_series().dropna() * lam1_e
        phi2 = obs['L5X'].sel(sv=sat).to_series().dropna() * lam5_e
        common = phi1.index.intersection(phi2.index)
        if len(common) > 5:
            L4 = phi1[common] - phi2[common]
            dL4 = L4.diff().dropna()
            n_slips = (np.abs(dL4) > CS_THRESHOLD).sum()

    rec['n_slips'] = n_slips
    rec['o_slps']  = rec['n_obs'] / (n_slips + 1)   # +1 avoids division by zero

    # ── Multipath proxy MP1 (GPS only, needs C1C + L1C + L2W) ──
    mp1_rms = np.nan
    if prefix == 'G' and all(c in obs.data_vars for c in ['C1C', 'L1C', 'L2W']):
        try:
            code = obs['C1C'].sel(sv=sat).to_series().dropna()
            ph1  = obs['L1C'].sel(sv=sat).to_series().dropna() * LAM1
            ph2  = obs['L2W'].sel(sv=sat).to_series().dropna() * LAM2
            comm = code.index.intersection(ph1.index).intersection(ph2.index)
            if len(comm) > 10:
                # MP1 = C1 - L1 - 2*alpha*(L1 - L2)
                mp1 = code[comm] - ph1[comm] - ALPHA_MP * (ph1[comm] - ph2[comm])
                # Remove mean (contains ambiguity and ionosphere terms)
                mp1 = mp1 - mp1.mean()
                mp1_rms = np.sqrt(np.mean(mp1**2))
        except Exception:
            pass
    rec['mp1_rms'] = mp1_rms

    results.append(rec)

# Convert to DataFrame for easy handling
df = pd.DataFrame(results).set_index('sat')

print(f"✅ Quality metrics computed for {len(df)} satellites")
print()

# Quick summary
print(f"{'Sat':<6} {'System':<10} {'Coverage%':>10} {'MeanSNR':>9} {'LowSNR%':>9} {'Slips':>7} {'o/slps':>8} {'MP1_RMS':>9}")
print("-" * 75)
for sat, row in df.sort_values('coverage', ascending=False).iterrows():
    mp = f"{row['mp1_rms']:.3f} m" if not np.isnan(row['mp1_rms']) else "  N/A  "
    snr = f"{row['mean_snr']:.1f}" if not np.isnan(row['mean_snr']) else " N/A"
    low = f"{row['pct_low_snr']:.1f}" if not np.isnan(row['pct_low_snr']) else " N/A"
    print(f"  {sat:<4} {row['system']:<10} {row['coverage']:>9.1f}% {snr:>9} {low:>9}%"
          f" {int(row['n_slips']):>7} {row['o_slps']:>8.0f} {mp:>9}")

# ─────────────────────────────────────────
# Step 4 — Plot 1: Coverage & SNR Dashboard
# ─────────────────────────────────────────
# Sorting satellites by constellation then PRN

df_sorted = df.sort_values(['system', 'sat'])
sats      = df_sorted.index.tolist()
x         = np.arange(len(sats))

# Bar colours by constellation
bar_colors = [CONST_COLORS.get(s[0], '#888888') for s in sats]

fig = plt.figure(figsize=(18, 14), facecolor='#0d1117')
gs  = gridspec.GridSpec(4, 1, hspace=0.5, figure=fig)

fig.suptitle(
    'GNSS Data Quality Dashboard | AUCK00NZL | Auckland, New Zealand | 2026-01-01\n'
    'All constellations — sorted by system',
    fontsize=14, fontweight='bold', color='#ffffff', y=0.98
)

def style_ax(ax, title):
    ax.set_facecolor('#111827')
    ax.tick_params(colors='#aaaaaa', labelsize=7)
    ax.grid(True, axis='y', color='#222222', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')
    ax.set_title(title, color='white', fontsize=10, loc='left', pad=4)
    ax.set_xlim(-0.5, len(sats) - 0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(sats, rotation=70, ha='right', fontsize=6.5)
    for tick, color in zip(ax.get_xticklabels(), bar_colors):
        tick.set_color(color)

# ─── Panel 1: Coverage % ────────────────────────────────────
ax1 = fig.add_subplot(gs[0])
ax1.bar(x, df_sorted['coverage'], color=bar_colors, edgecolor='#0d1117', lw=0.5)
ax1.axhline(80, color='#FFEB3B', ls='--', lw=1.0, label='80% threshold')
ax1.set_ylabel('Coverage [%]', color='#aaaaaa', fontsize=9)
ax1.set_ylim(0, 110)
style_ax(ax1, '① Observation Coverage  (% of total epochs with valid pseudorange)')
ax1.legend(fontsize=8, framealpha=0.3, facecolor='#1a1a2e',
           edgecolor='#444444', labelcolor='white')

# ─── Panel 2: Mean SNR ──────────────────────────────────────
ax2 = fig.add_subplot(gs[1])
snr_vals = df_sorted['mean_snr'].fillna(0)
ax2.bar(x, snr_vals, color=bar_colors, edgecolor='#0d1117', lw=0.5)
ax2.axhline(35, color='#4CAF50', ls='--', lw=1.0, label='Good (35 dB-Hz)')
ax2.axhline(25, color='#F44336', ls='--', lw=1.0, label='Poor (25 dB-Hz)')
ax2.set_ylabel('Mean SNR [dB-Hz]', color='#aaaaaa', fontsize=9)
ax2.set_ylim(0, 55)
style_ax(ax2, '② Mean Signal-to-Noise Ratio L1  (higher = better signal quality)')
ax2.legend(fontsize=8, framealpha=0.3, facecolor='#1a1a2e',
           edgecolor='#444444', labelcolor='white')

# ─── Panel 3: Cycle slips ───────────────────────────────────
ax3 = fig.add_subplot(gs[2])
slip_vals = df_sorted['n_slips'].fillna(0)
slip_colors = ['#F44336' if v > 5 else '#FF9800' if v > 1 else '#4CAF50'
               for v in slip_vals]
ax3.bar(x, slip_vals, color=slip_colors, edgecolor='#0d1117', lw=0.5)
ax3.axhline(5, color='#F44336', ls='--', lw=1.0, label='>5 slips = concern')
ax3.set_ylabel('Cycle Slips (L4)', color='#aaaaaa', fontsize=9)
style_ax(ax3, '③ Cycle Slips Detected via L4  (GPS & Galileo only — needs dual-frequency)')
ax3.legend(fontsize=8, framealpha=0.3, facecolor='#1a1a2e',
           edgecolor='#444444', labelcolor='white')

# ─── Panel 4: MP1 RMS ───────────────────────────────────────
ax4 = fig.add_subplot(gs[3])
mp_mask   = df_sorted['mp1_rms'].notna()
mp_sats   = df_sorted[mp_mask]
mp_x      = [sats.index(s) for s in mp_sats.index]
mp_vals   = mp_sats['mp1_rms'].values
mp_colors_list = ['#F44336' if v > 0.5 else '#FF9800' if v > 0.3 else '#4CAF50'
                  for v in mp_vals]
ax4.bar(mp_x, mp_vals, color=mp_colors_list, edgecolor='#0d1117', lw=0.5, width=0.8)
ax4.axhline(0.3, color='#FF9800', ls='--', lw=1.0, label='Moderate MP (0.3 m)')
ax4.axhline(0.5, color='#F44336', ls='--', lw=1.0, label='High MP (0.5 m)')
ax4.set_ylabel('MP1 RMS [m]', color='#aaaaaa', fontsize=9)
style_ax(ax4, '④ Multipath Proxy MP1 RMS  (GPS only — lower = cleaner environment)')
ax4.legend(fontsize=8, framealpha=0.3, facecolor='#1a1a2e',
           edgecolor='#444444', labelcolor='white')

# Constellation legend (coloured patches on the right)
from matplotlib.patches import Patch
legend_patches = [Patch(color=v, label=CONST_NAMES.get(k, k))
                  for k, v in CONST_COLORS.items()
                  if any(s.startswith(k) for s in sats)]
fig.legend(
    handles=legend_patches,
    loc='lower center', ncol=len(legend_patches),
    fontsize=10, framealpha=0.2,
    facecolor='#1a1a2e', edgecolor='#444444', labelcolor='white',
    bbox_to_anchor=(0.5, 0.01)
)

plt.savefig('/../content/output/plot1_quality_dashboard.png', dpi=150,
            bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Plot saved: plot1_quality_dashboard.png')
print()
print('💡 Interpretation:')
print('   • Coverage: any satellite < 20% had very limited visibility or signal issues')
print('   • SNR:      bars below the red line (25 dB-Hz) indicate poor signal quality')
print('   • Slips:    red bars = satellites with frequent phase interruptions')
print('   • MP1:      red bars = high multipath (reflections from nearby objects)')

# ─────────────────────────────
# Step 5 — Plot 2: Data Gap Map
# ─────────────────────────────
# Only showing satellites with at least some data

active_sats = df_sorted[df_sorted['coverage'] > 1].index.tolist()
n_active    = len(active_sats)

# Determine best pseudorange code
pr_code = 'C1C' if 'C1C' in obs.data_vars else 'C1X'

# Build presence matrix: 1 = data present, 0 = gap
presence_matrix = np.zeros((n_active, total_epochs))

for i, sat in enumerate(active_sats):
    try:
        series = obs[pr_code].sel(sv=sat).to_series().reindex(time_index)
        presence_matrix[i, :] = series.notna().astype(int).values
    except Exception:
        pass

sat_colors_gap = [CONST_COLORS.get(s[0], '#888888') for s in active_sats]

# Custom colourmap: white = gap, coloured = data
# We use a simple binary: 0 → white, 1 → dark blue

row_height = 0.25
fig_height = max(6, n_active * row_height + 2.5)

fig, ax = plt.subplots(figsize=(16, fig_height), facecolor='#0d1117')
ax.set_facecolor('#0d1117')

# Custom binary colourmap: 0 = white/light grey (gap), 1 = dark (data)
gap_cmap = LinearSegmentedColormap.from_list(
    'gap_map', ['#1a1a2e', '#00e676'], N=2
)

ax.imshow(
    presence_matrix,
    aspect='auto',
    cmap=gap_cmap,
    vmin=0, vmax=1,
    extent=[
        mdates.date2num(time_index[0]),
        mdates.date2num(time_index[-1]),
        -0.5, n_active - 0.5
    ],
    origin='upper'
)
ax.xaxis_date()

ax.set_yticks(range(n_active))
ax.set_yticklabels(active_sats, fontsize=7.5)
for tick, color in zip(ax.get_yticklabels(), sat_colors_gap):
    tick.set_color(color)

ax.set_title(
    f'Observation Gap Map — {n_active} Active Satellites | AUCK00NZL | 2026-01-01\n'
    'Green = data present  |  Dark = no observation (satellite below horizon or signal lost)',
    fontsize=12, fontweight='bold', color='#ffffff'
)
ax.set_xlabel('UTC Time (HH:MM)', fontsize=11, color='#e0e0e0')

ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
plt.xticks(rotation=30, color='#aaaaaa')

for spine in ax.spines.values():
    spine.set_edgecolor('#333333')

plt.tight_layout()
plt.savefig('/../content/output/plot2_gap_map.png', dpi=150,
            bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Plot saved: plot2_gap_map.png')
print()
print('💡 Interpretation:')
print('   • Green rows = satellite continuously tracked during this window')
print('   • Fragmented rows = satellite near horizon, frequently entering/leaving mask angle')
print('   • Completely dark rows = satellite not visible from Auckland on this day')

# ──────────────────────────────────────────────
# Step 6 — Plot 3: Quality Score Summary Scatter
# ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 8), facecolor='#0d1117')
ax.set_facecolor('#111827')

ax.set_title(
    'Satellite Quality Score: Coverage vs Mean SNR | AUCK00NZL | 2026-01-01\n'
    'Top-right = best quality  |  Dot size = cycle slips detected (larger = more slips)',
    fontsize=12, fontweight='bold', color='#ffffff'
)

# Quality quadrant shading
ax.axhspan(35, 60, xmin=0.8, alpha=0.05, color='#00e676', label='Ideal zone')
ax.axhline(35, color='#4CAF50', ls='--', lw=0.8, alpha=0.6)
ax.axhline(25, color='#F44336', ls='--', lw=0.8, alpha=0.6)
ax.axvline(80, color='#4CAF50', ls='--', lw=0.8, alpha=0.6)

# Text labels for quality zones
ax.text(82, 36, 'Ideal\n(high coverage,\ngood signal)',
        color='#00e676', fontsize=8, alpha=0.7)
ax.text(82, 24, 'Poor SNR', color='#F44336', fontsize=8, alpha=0.7)
ax.text(2,  36, 'Low visibility', color='#FF9800', fontsize=8, alpha=0.7)

for prefix in ['G', 'R', 'E', 'C', 'J', 'S']:
    sub = df_sorted[df_sorted.index.str.startswith(prefix)]
    if sub.empty:
        continue

    valid = sub.dropna(subset=['mean_snr'])
    if valid.empty:
        continue

    # Dot size scales with cycle slips (larger = more slips)
    sizes = 80 + valid['n_slips'] * 30

    sc = ax.scatter(
        valid['coverage'],
        valid['mean_snr'],
        c=CONST_COLORS.get(prefix, '#888888'),
        s=sizes,
        alpha=0.85,
        edgecolors='white',
        linewidths=0.5,
        label=CONST_NAMES.get(prefix, prefix),
        zorder=5
    )

    # Label each satellite
    for sat, row in valid.iterrows():
        ax.annotate(
            sat,
            xy=(row['coverage'], row['mean_snr']),
            xytext=(3, 3), textcoords='offset points',
            fontsize=6.5, color=CONST_COLORS.get(prefix, 'white'), alpha=0.9
        )

ax.set_xlabel('Coverage [% of total epochs]', color='#aaaaaa', fontsize=11)
ax.set_ylabel('Mean SNR [dB-Hz]', color='#aaaaaa', fontsize=11)
ax.set_xlim(-2, 105)
ax.set_ylim(15, 55)

ax.tick_params(colors='#aaaaaa')
ax.grid(True, color='#222222', linewidth=0.5)
for spine in ax.spines.values():
    spine.set_edgecolor('#333333')

legend = ax.legend(
    fontsize=10, loc='lower right',
    framealpha=0.3, facecolor='#1a1a2e', edgecolor='#444444'
)
for t in legend.get_texts():
    t.set_color('white')

plt.tight_layout()
plt.savefig('/../content/output/plot3_quality_scatter.png', dpi=150,
            bbox_inches='tight', facecolor=fig.get_facecolor())
plt.show()

print('✅ Plot saved: plot3_quality_scatter.png')
print()
print('💡 Interpretation:')
print('   • Top-right = best performing satellites (long arc, strong signal)')
print('   • Bottom-left = satellites that barely rose above the horizon')
print('   • Larger dot = more cycle slips detected during the tracking arc')
print('   • Red dashed lines = quality thresholds (SNR 25/35 dB-Hz, coverage 80%)')

# ─────────────────────────────────────────
# Step 7 — Generate the Text Quality Report
# ─────────────────────────────────────────

# Compute file-level summary statistics

n_sats_total   = len(df)
n_sats_active  = (df['coverage'] > 5).sum()
mean_coverage  = df[df['coverage'] > 5]['coverage'].mean()
mean_snr_all   = df['mean_snr'].dropna().mean()
total_slips    = int(df['n_slips'].sum())
gps_mp_mean    = df[df.index.str.startswith('G')]['mp1_rms'].dropna().mean()

# Constellations present
prefixes_present = sorted(set(s[0] for s in df.index if df.loc[s, 'coverage'] > 5))
const_str = ', '.join(CONST_NAMES.get(p, p) for p in prefixes_present)

# Worst performing satellites (coverage < 20% or SNR < 28)
poor_cov  = df[(df['coverage'] > 1) & (df['coverage'] < 20)].index.tolist()
poor_snr  = df[df['mean_snr'] < 28].dropna(subset=['mean_snr']).index.tolist()
high_slip = df[df['n_slips'] > 5].index.tolist()

# Build report string

now = datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')

report_lines = [
    "=" * 70,
    "GNSS DATA QUALITY REPORT",
    "=" * 70,
    f"Generated        : {now}",
    f"RINEX file       : {obs_path}",
    f"Station          : AUCK00NZL — Auckland, New Zealand",
    f"Date             : 2026-01-01  (Day of Year 001)",
    f"Sampling rate    : 30 seconds",
    f"Total epochs     : {total_epochs:,}",
    f"Duration         : {total_epochs * 30 / 3600:.1f} hours",
    "",
    "=" * 70,
    "CONSTELLATION SUMMARY",
    "=" * 70,
    f"Systems present  : {const_str}",
    f"Total SVs in file: {n_sats_total}",
    f"Active SVs (>5%) : {n_sats_active}",
    "",
]

for prefix in prefixes_present:
    sub = df[df.index.str.startswith(prefix) & (df['coverage'] > 1)]
    report_lines += [
        f"  {CONST_NAMES.get(prefix, prefix):<10}: {len(sub):>2} satellites | "
        f"coverage {sub['coverage'].mean():.1f}% avg | "
        f"SNR {sub['mean_snr'].dropna().mean():.1f} dB-Hz avg"
    ]

report_lines += [
    "",
    "=" * 70,
    "SIGNAL QUALITY",
    "=" * 70,
    f"Mean SNR (all SVs)     : {mean_snr_all:.1f} dB-Hz",
    f"Mean coverage (active) : {mean_coverage:.1f}%",
    f"Total cycle slips (L4) : {total_slips}  (GPS + Galileo)",
    f"GPS MP1 RMS mean       : {gps_mp_mean:.3f} m" if not np.isnan(gps_mp_mean) else "GPS MP1 RMS mean  : N/A",
    "",
    "=" * 70,
    "QUALITY FLAGS",
    "=" * 70,
]

if poor_cov:
    report_lines.append(f"  ⚠ Low coverage (<20%)   : {', '.join(poor_cov)}")
else:
    report_lines.append("  ✅ No satellites with unusually low coverage")

if poor_snr:
    report_lines.append(f"  ⚠ Low SNR (<28 dB-Hz)   : {', '.join(poor_snr)}")
else:
    report_lines.append("  ✅ All satellites with good SNR (>28 dB-Hz)")

if high_slip:
    report_lines.append(f"  ⚠ High cycle slips (>5)  : {', '.join(high_slip)}")
else:
    report_lines.append("  ✅ No satellites with excessive cycle slips")

if not np.isnan(gps_mp_mean):
    mp_status = '✅ Clean' if gps_mp_mean < 0.3 else '⚠ Moderate' if gps_mp_mean < 0.5 else '❌ High'
    report_lines.append(f"  GPS multipath MP1       : {mp_status} ({gps_mp_mean:.3f} m RMS mean)")

report_lines += [
    "",
    "=" * 70,
    "PER-SATELLITE DETAIL",
    "=" * 70,
    f"{'Sat':<6} {'System':<10} {'Cov%':>7} {'SNR':>8} {'LowSNR%':>9} {'Slips':>7} {'o/slps':>8} {'MP1':>8}",
    "-" * 70,
]

for sat, row in df.sort_values('coverage', ascending=False).iterrows():
    mp = f"{row['mp1_rms']:.3f}" if not np.isnan(row['mp1_rms']) else "  N/A "
    snr = f"{row['mean_snr']:.1f}" if not np.isnan(row['mean_snr']) else " N/A"
    low = f"{row['pct_low_snr']:.1f}" if not np.isnan(row['pct_low_snr']) else " N/A"
    report_lines.append(
        f"  {sat:<4} {row['system']:<10} {row['coverage']:>6.1f}% {snr:>8}"
        f" {low:>9}% {int(row['n_slips']):>7} {row['o_slps']:>8.0f} {mp:>8}"
    )

report_lines += [
    "",
    "=" * 70,
    "END OF REPORT",
    "=" * 70,
]

report_text = "\n".join(report_lines)

# Save to file
with open('quality_report.txt', 'w') as f:
    f.write(report_text)

# Print to screen
print(report_text)


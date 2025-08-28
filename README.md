# 🌌 Star Party Planner

**Plan an observing night with confidence!**  
This tool generates observing lists tailored to your site, date, and time window. It prioritizes crowd-pleasers such as Saturn and Jupiter, filters by altitude and Moon separation, and outputs **night-vision-friendly HTML**, CSV tables, and console listings.

---

## ✨ Features

- Planets + Moon + **custom deep-sky catalog** (CSV, e.g. Messier + Caldwell)
- Filters by **minimum altitude**, **maximum magnitude**, and **Moon separation**
- Assigns an **interest score** to prioritize crowd favorites
- Finds each object’s **best observing time** (max altitude in your window)
- Outputs:
  - **Master list** (sorted by interestingness, then best time)
  - **Hourly “point your scope now” tables**
- Exports to:
  - `*_master.csv` – Master observing list
  - `*_hourly.csv` – Hourly shortlist
  - `*.html` – Red-on-black, night-vision-safe web page (accordion or tabbed UI)

---

## 📦 Installation

Requires **Python 3.10+**.

1. Clone or download this repository:

   ```bash
   git clone https://github.com/yourusername/star-party-planner.git
   cd star-party-planner
   ```

2. Install dependencies:

   ```bash
   pip install skyfield numpy pandas pytz python-dateutil
   ```

   Recommended versions:
   - [Skyfield](https://rhodesmill.org/skyfield/) `1.45+`
   - NumPy `1.20+`
   - Pandas `1.3+`

---

## 🪐 Ephemeris Data (`de440s.bsp`)

This script **requires a local planetary ephemeris file**.  
We recommend **`de440s.bsp`** (~31 MB, accurate 1550–2650 CE).

1. Create the data directory:

   ```bash
   mkdir -p ./skyfield_data
   ```

2. Download the file from NAIF (NASA/JPL):

   ```bash
   curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
   ```

3. Move it into place:

   ```bash
   mv de440s.bsp ./skyfield_data/
   ```

4. Verify you have:

   ```
   ./skyfield_data/de440s.bsp
   ```

> ⚠️ **Important**: This script is designed to **never fetch files from the internet automatically.** You must provide the BSP manually with `--bsp`.

---

## 📊 Usage

### Example: Eagle Lake Observatory, Minnesota

```bash
python starparty_planner.py \
  --lat 44.810265 --lon -93.939783 --elev 296 \
  --date 2025-08-30 --start 20:00 --end 01:00 --tz America/Chicago \
  --catalog messier_caldwell.csv --min_alt 20 --moon_sep_min 20 --max_mag 9 \
  --top_n_per_hour 16 --out_prefix starparty \
  --html starparty.html \
  --bsp ./skyfield_data/de440s.bsp \
  --html_ui tabs \
  --min_alt_planets 5 --min_alt_moon 0
```

### Outputs

- `starparty_master.csv` – Master list (sorted by interestingness)
- `starparty_hourly.csv` – Per-hour shortlist
- `starparty.html` – Red-on-black web UI  
  - **Accordion mode** (`--html_ui accordion`) – collapsible per-hour sections  
  - **Tabs mode** (`--html_ui tabs`) – Master list / By Hour tabs

### Defaults

- If `--date` is missing → defaults to **today**
- If `--start` is missing → defaults to **local sunset (rounded to nearest hour)**
- If `--end` is missing → defaults to **hour of local sunrise**
- If `--min_alt_planets` or `--min_alt_moon` are missing → uses safer low defaults (10° and 5° respectively)

---

## 📂 Input Catalogs

- Catalogs are **CSV files** with required columns:

  ```
  name,ra_deg,dec_deg,type,mag,notes
  ```

- Example row:

  ```
  M13,250.421,36.461,Globular cluster,5.8,"Great Hercules Cluster"
  ```

- This repository distributes a ready-to-use **`messier_caldwell.csv`** with all Messier and Caldwell objects.

---

## 🌙 HTML Output

- Fully **night-vision safe** (red on black; no white autofill)
- **Searchable tables** (instant filtering by name, type, or notes)
- **Sticky controls**: search box, hour links, and tabs remain visible
- **Responsive design**: mobile, tablet, and laptop friendly
- **Sortable Master List**: click column headers (Name, Type, Mag, Best Local Time)

---

## ⚡ Tips & Tricks

- Use `--top_n_per_hour` to control list length (e.g., 5 for a tight schedule).
- Adjust `--moon_sep_min` to relax/tighten Moon glare filtering.
- Planets and the Moon have **separate minimum altitude thresholds**, so Saturn doesn’t get excluded unfairly.
- On iOS/Android, save the HTML to your homescreen for quick night-vision-friendly access.
- CSVs can be imported into Excel, Google Sheets, or SkySafari for custom planning.

---

## 🔭 Example Workflow

1. Prepare your DSO catalog (`messier_caldwell.csv`).
2. Download ephemeris (`de440s.bsp`).
3. Run the script for your observing site.
4. Open `starparty.html` on your phone/tablet in night mode.
5. Enjoy your observing session!

---

## 🧑‍💻 Development Notes

- Written in **Python**
- No network access is attempted by Skyfield; all data is local
- Extendable:
  - Add new columns to catalogs
  - Adjust scoring system in `interest_score()`
  - Customize HTML/CSS for different UI themes

---

## 📜 License

MIT License.  
Created to help star parties run smoothly and maximize “wow” factor objects in the eyepiece.

---

🚀 *Clear skies and good seeing!*

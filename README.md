# 🪐 Star Party Planner

**Plan an observing night with confidence!**  
This tool generates observing lists tailored to your site, date, and time window. It prioritizes crowd-pleasers such as Saturn and Jupiter, filters by altitude and Moon separation, and outputs **night-vision-friendly HTML**, CSV tables, and console listings.

---

## Highlights

- Planets + Moon + **custom deep-sky catalog** (CSV; Messier + Caldwell included)
- Filters by **minimum altitude**, **maximum magnitude**, and **Moon separation**
- Assigns an **interest score** (crowd-pleasers bubble to the top)
- Finds each object’s **best observing time** (peak altitude in your window)
- Night‑vision HTML with **dark red theme** and mobile-first layout
- **Clickable rows** → modal with object details and a **DSS2 (red) preview**
- **Local image cache** to avoid re-downloading previews daily

Outputs:
- `*_master.csv` – Master observing list
- `*_hourly.csv` – Hour‑by‑hour shortlist
- `*.html` – Night‑vision web page (accordion or tabbed UI)

---

## Installation

Requires **Python 3.10+**.

```bash
git clone https://github.com/educationfutures/star-party-planner.git
cd star-party-planner

# Core deps
pip install skyfield numpy pandas pytz python-dateutil
```

Recommended versions:
- Skyfield `1.45+`
- NumPy `1.20+`
- Pandas `1.3+`

---

## Ephemeris (`de440s.bsp`)

The script uses a **local** JPL ephemeris for precise planet/Moon positions. We recommend **`de440s.bsp`** (~31 MB, valid 1550–2650 CE).

```bash
mkdir -p ./skyfield_data
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
mv de440s.bsp ./skyfield_data/
# You should now have: ./skyfield_data/de440s.bsp
```

> The planner never auto-downloads BSP files; you must provide them with `--bsp`.

---

## Catalogs

A CSV catalog is required with at least the following columns:

```
name,ra_deg,dec_deg,type,mag,notes
```

Example row:

```
M13,250.421,36.461,Globular cluster,5.8,"Great Hercules Cluster"
```

This repository ships with **`messier_caldwell.csv`** covering Messier + Caldwell objects.

---

## Usage

### Example (Eagle Lake Observatory, Minnesota)

```bash
python starparty_planner.py \
  --lat 44.810265 --lon -93.939783 --elev 296 \
  --date 2025-08-30 --start 20:00 --end 01:00 --tz America/Chicago \
  --catalog messier_caldwell.csv --min_alt 20 --moon_sep_min 20 --max_mag 9 \
  --top_n_per_hour 16 --out_prefix starparty \
  --html starparty.html \
  --bsp ./skyfield_data/de440s.bsp \
  --html_ui tabs \
  --min_alt_planets 5 --min_alt_moon 0 \
  --preview_cache_dir image_cache --preview_px 800 --preview_fov_deg 0.6
```

> ⚠️ This script may run slowly when downloading previews of images for the first time.

---

## Command-Line Arguments

### Core Location & Time
- `--lat` (float, required) – Latitude in decimal degrees (north positive)  
- `--lon` (float, required) – Longitude in decimal degrees (east positive)  
- `--elev` (float, default=0.0) – Elevation in meters
- `--date` (YYYY-MM-DD, default=today) – Observation date (local)  
- `--start` (HH:MM, default=sunset rounded) – Session start time  
- `--end` (HH:MM, default=sunrise hour) – Session end time  
- `--tz` (string, default=UTC) – IANA timezone, e.g. `America/Chicago`

### Catalog & Filtering
- `--catalog` (path, default=`messier_caldwell.csv`) – Input catalog of deep-sky objects  
- `--min_alt` (float, default=20.0) – Minimum altitude (deg) for DSOs  
- `--max_mag` (float, default=9.0) – Maximum magnitude (fainter objects are excluded)  
- `--moon_sep_min` (float, default=15.0) – Minimum separation from the Moon (deg)  
- `--hour_step` (int, default=1) – Step size in hours for hourly tables  
- `--top_n_per_hour` (int, default=16) – Max targets shown per hour slot  

### Outputs
- `--out_prefix` (string, default=`starparty`) – Prefix for CSV output files  
- `--html` (string, optional) – Path to save HTML output file  
- `--bsp` (path, default=`./skyfield_data/de440s.bsp`) – Planetary ephemeris file  
- `--html_ui` (`tabs` | `accordion`, default=`tabs`) – HTML layout style  

### Per-Type Altitude Thresholds
- `--min_alt_planets` (float, default=5.0) – Minimum altitude for planets  
- `--min_alt_moon` (float, default=0.0) – Minimum altitude for the Moon  

### Preview & Caching
- `--no_previews` – Disable image previews entirely  
- `--preview_cache_dir` (string, default=`image_cache`) – Directory for cached previews  
- `--preview_px` (int, default=800) – Pixel size for DSS2 images  
- `--preview_fov_deg` (float, default=0.6) – Field of view for DSS2 images (degrees)  
- `--refresh_previews` – Force re-download of all previews  
- `--clean_preview_cache` – Remove cached images not used in this run  

> **Note:** `--cache_dir` and `--clean_cache` apply to the *HTML preview images* (DSS2). Catalogs and BSPs are not affected.

---

## HTML Output & UI

- **Red-on-black** theme suitable for night use (avoid max screen brightness).
- **Search** input filters both Master and Hourly tables.
- Master & By‑Hour views via tabs (or per-hour accordions).
- **Non-sticky** table headers (fixes overlay on resize); sticky controls only for the top navbar.
- Columns emphasize what you need at the eyepiece: **Alt (°), Az (°), Dir**.

**Row Modals**  
Click any row to open a modal showing:

- Object metadata (type, magnitude, best time)
- **DSS2 Red** preview (via NASA SkyView) – rendered through a red filter
- Notes and coordinates (RA/Dec are hidden in the table but available to the modal)

> If SkyView is temporarily slow or blocked by the browser, the modal will still open with the object’s data. The image may appear after a delay or not at all depending on network and CORS conditions.

---

## Image Previews & Caching

The HTML output includes **clickable rows** that open a modal with details and a preview.

### Source
- **NASA SkyView (DSS2 Red)** for deep‑sky targets and planets (for a consistent, red‑friendly look).

### Night‑Vision Friendly
- Previews are displayed with a **red‑only filter**, keeping the interface dark-sky friendly.

### Local Cache
- Previews are **cached** in `--cache_dir` (default `./.cache`).  
- On subsequent runs, existing images are reused to keep the page fast and minimize external fetches.

### Cache Maintenance
- **Wipe & rebuild** the preview cache:
  ```bash
  python starparty_planner.py ... --clean_cache
  ```

> Tip: Keep your cache directory on persistent storage (e.g., not `/tmp`) so your cron job can reuse it daily.

---

## Example Workflow

1. Prepare your DSO catalog (`messier_caldwell.csv`).
2. Download ephemeris (`de440s.bsp`) into `./skyfield_data/`.
3. Run the script for your site with your preferred flags.
4. Open `starparty.html` on your phone/tablet at the telescope.
5. Enjoy your observing session!

💡 **Automation:** run this in a **daily cron job** and publish the HTML to your website so the plan is always fresh before sunset.

---

## Deploy to a Web Server (Daily at 7:00 AM)

Assume you want the output available at `/var/www/star-party/index.html` on a Linux host running Nginx/Apache.

1. Create a writable output directory and set ownership (replace `www-data:www-data` with your web user):
   ```bash
   sudo mkdir -p /var/www/star-party
   sudo chown -R $USER:www-data /var/www/star-party
   sudo chmod -R 775 /var/www/star-party
   ```

2. Create a small shell script (e.g., `~/bin/run_starparty.sh`):
   ```bash
   #!/usr/bin/env bash
   set -euo pipefail

   REPO="$HOME/star-party-planner"
   cd "$REPO"

   python starparty_planner.py \
     --lat 44.810265 --lon -93.939783 --elev 296 \
     --tz America/Chicago \
     --catalog messier_caldwell.csv --min_alt 20 --moon_sep_min 20 --max_mag 9 \
     --top_n_per_hour 16 --out_prefix starparty \
     --html /var/www/star-party/index.html \
     --bsp ./skyfield_data/de440s.bsp \
     --html_ui tabs \
     --min_alt_planets 5 --min_alt_moon 0 \
     --cache_dir /var/www/star-party/.cache \
     --clean_cache
   ```

   Make it executable:
   ```bash
   chmod +x ~/bin/run_starparty.sh
   ```

3. Add a cron entry to run daily at **07:00**:
   ```bash
   crontab -e
   ```

   Add this line:
   ```
   0 7 * * * /home/youruser/bin/run_starparty.sh >> /home/youruser/starparty_cron.log 2>&1
   ```

   > The script will compute **today’s date**, **sunset**, and the **sunrise hour** automatically when you omit `--date`, `--start`, and `--end`.

4. Visit `http(s)://your-server/star-party/` to see the latest plan.

---

## Troubleshooting

- **Blank/white preview**: Some browsers block cross‑origin iframes or image loads under certain privacy modes. Try a normal tab or allow the request in your content blocker. Cached images will still display if present.
- **No BSP found**: Ensure `--bsp` points to `./skyfield_data/de440s.bsp` and that the file exists.
- **Tables seem empty**: Relax filters (`--min_alt`, `--moon_sep_min`, `--max_mag`) or expand the time window.
- **Wrong times**: Confirm your `--tz` and system clock. The planner assumes an IANA timezone string.
- **Permissions**: Make sure the web server user can read the generated HTML and the cache directory.

---

## Development Notes

- Written in **Python**, using Skyfield for accurate ephemerides
- All heavy astronomy data (BSP) is **local**; the planner does not auto-fetch ephemerides
- The UI is plain HTML/CSS/JS with a night‑vision theme
- Extendable:
  - Change scoring in `interest_score()`
  - Add columns to your catalog(s)
  - Tweak UI (accordion vs tabs)

---

## License

MIT License

---

✨ *Clear skies!*
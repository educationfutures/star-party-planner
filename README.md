# 🪐 Star Party Planner

**Plan an observing night with confidence!**  
This tool generates observing lists tailored to your site, date, and time window. It prioritizes crowd‑pleasers (Saturn! Jupiter!), filters by altitude and Moon separation, and outputs **night‑vision friendly HTML**, CSV tables, and console listings.

A daily demo (Minnesota Astronomical Society — Eagle Lake Observatory) lives at https://party.knowmad.org

---

## Highlights

- Planets + Moon + **custom deep‑sky catalog** (CSV; Messier + Caldwell included)
- Smart filtering by **minimum altitude**, **maximum magnitude**, and **Moon separation**
- Assigns an **interest score** so crowd‑pleasers bubble to the top
- Supports **type_weights.csv** for customizable interest scoring per object type
- Finds each object’s **best observing time** (peak altitude in your window)
- **Night‑vision HTML**: dark red theme, mobile‑first, lightweight JS
- **Now view** that auto‑updates the HTML view in sync with the clock (every 4 minutes)
- **Directional/altitude “obstruction” filters** (up to 5 azimuth ranges)
- **Top‑N limiter** for the Now view (e.g., show top 16 targets only)
- **Clickable rows** → modal with object details and a preview image
- **Local image cache** (DSS2 Red via hips2fits; Wikipedia fallback) to avoid re‑downloading

**Outputs**
- `*_master.csv` – Master observing list
- `*_hourly.csv` – Hour‑by‑hour shortlist
- `*.html` – Night‑vision web page (tabbed UI with Now + Master; hourly tabs)

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

The planner uses a **local** JPL ephemeris for precise planet/Moon positions. We recommend **`de440s.bsp`** (~31 MB, valid 1550–2650 CE).

```bash
mkdir -p ./skyfield_data
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
mv de440s.bsp ./skyfield_data/
# You should now have: ./skyfield_data/de440s.bsp
```

> ⚠️ The planner never auto‑downloads BSP files; you must provide them with `--bsp`.

---

## Catalogs

Provide a CSV with at least the following columns:

```
name,ra_deg,dec_deg,type,mag,notes
```

Example row:

```
M13,250.421,36.461,Globular cluster,5.8,"Great Hercules Cluster"
```

This repo includes **`messier_caldwell.csv`** (Messier + Caldwell) and **`extended_targets.csv`** with additional crowd‑pleasers.

---

## Type Weights

The `type_weights.csv` file allows you to customize how different object types are scored in the prioritization process. Each row should contain a `type` and a corresponding `weight`, where higher weights increase the priority of those object types in the observing list.

The repository includes a default `type_weights.csv` tuned for a 16-inch telescope. You can copy and edit this file to suit your preferences or provide your own via the `--type_weights` command-line argument.

Example `type_weights.csv` snippet:

```
type,weight,diffuse,partial_diffuse
Globular cluster,45,FALSE,FALSE
Open cluster,33,FALSE,FALSE
Planetary nebula,40,FALSE,TRUE
Galaxy cluster,32,TRUE,FALSE
```

Weights act as relative multipliers to the base interest score and work together with other factors such as magnitude, altitude, and moonlight penalty to determine the final priority of each target.

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
  --preview_cache_dir image_cache --preview_px 800 --preview_fov_deg 0.6 \
  --type_weights type_weights.csv \
  --now_padding_min 30
```

> ⚠️ First run may spend a long time fetching previews to the local cache.

---

## Command‑Line Arguments

### Core Location & Time
- `--lat` (float, required) – Latitude in decimal degrees (north +)  
- `--lon` (float, required) – Longitude in decimal degrees (east +)  
- `--elev` (float, default=0.0) – Elevation in meters
- `--date` (YYYY‑MM‑DD, default=today) – Observation date (local)  
- `--start` (HH:MM, default=sunset rounded) – Session start time  
- `--end` (HH:MM, default=sunrise hour) – Session end time  
- `--tz` (string, default=UTC) – IANA timezone, e.g. `America/Chicago`

### Catalog & Filtering
- `--catalog` (path, default=`messier_caldwell.csv`) – DSO catalog  
- `--bsp` (path, default=`./skyfield_data/de440s.bsp`) – Planetary ephemeris  
- `--min_alt` (float, default=20.0) – Minimum altitude (deg) for DSOs  
- `--max_mag` (float, default=9.0) – Maximum magnitude (fainter excluded)  
- `--moon_sep_min` (float, default=15.0) – Minimum separation from the Moon (deg) for DSOs  
- `--hour_step` (int, default=1) – Step size in hours for hourly tables  
- `--top_n_per_hour` (int, default=16) – Max targets per hour slot
- `--moonlight_penalty_max` (float, default=18) – Max points subtracted from diffuse targets at full Moon when high in the sky

### Scoring & Weights
- `--type_weights` (path, optional) – CSV file specifying per-type interest weights to customize scoring. Overrides default scoring weights for object types.

### Outputs
- `--out_prefix` (string, default=`starparty`) – Prefix for CSV files  
- `--html` (string, optional) – Path to save HTML output  
- `--html_title` (string, default=`Star Party Planner`) – Heading/title for HTML output  

- `--html_ui` (`tabs` | `accordion`, default=`tabs`) – HTML layout style

### Per‑Type Altitude Thresholds
- `--min_alt_planets` (float, default=5.0) – Minimum altitude for planets  
- `--min_alt_moon` (float, default=0.0) – Minimum altitude for the Moon

### Preview & Caching (HTML previews)
- `--no_previews` – Disable image previews entirely  
- `--preview_cache_dir` (string, default=`image_cache`) – Directory for cached previews  
- `--preview_px` (int, default=800) – Pixel size for DSS2 preview images  
- `--preview_fov_deg` (float, default=0.6) – Field of view for DSS2 preview images (degrees)  
- `--refresh_previews` – Force re‑download of previews this run  
- `--clean_preview_cache` – After generation, remove cached images not referenced by this run

### Now View Timing
- `--now_padding_min` (int, default=0) – Number of minutes before sunset and after sunrise to keep the Now view active (pads the observing window).

> **Note:** These options apply to the **HTML preview images** only. Catalogs and BSPs are not affected.

---

## HTML Output & UI

- **Tabs** for **Now** and **Master List**; per‑hour tables are reachable via the hour “chips” in the navbar.
- **Search** filters the visible table(s) instantly.
- **Sticky navbar**, non‑sticky table headers (prevents overlay on resize).
- Emphasis on eyepiece‑friendly columns: **Alt (°), Az (°), Dir**.

### Now view (auto‑refreshing)
- Updates in step with the clock (every four minutes).  
- Renders the **top‑N** targets after filters—configurable in the Filters dialog.  
- The Now view window can be padded before sunset and after sunrise by specifying `--now_padding_min`.
- Scoring now incorporates customizable per-type weights from `type_weights.csv` if provided.
- Click any row to open the modal with details + preview.

### Filters dialog
- **Enable/disable** directional filtering.  
- Define up to **five azimuth ranges** with a **minimum clear altitude** each (great for trees/buildings).  
- Set the **Top‑N** limit for the Now view.  
- Settings are stored in `localStorage` and persist per browser.

### Row modal & previews
- Shows type, magnitude, best time, direction, coordinates, and notes.  
- Preview image:
  - **DSS2 Red** via **hips2fits** (CDS service), or  
  - **Wikipedia/Wikimedia** image (planets/Moon, or DSO fallback).  
- Previews are tinted for night vision.

---

## Image Previews & Caching

- During HTML generation the planner fetches preview tiles and writes them to `--preview_cache_dir` (default: `image_cache`).
- Subsequent runs reuse existing files so the page loads instantly in the field.
- If a preview is missing or outdated, run with `--refresh_previews`.
- To prune unused files after a run, add `--clean_preview_cache`.

---

## Example Workflow

1. Prepare your DSO catalog (`messier_caldwell.csv`).
2. Download the ephemeris (`de440s.bsp`) into `./skyfield_data/`.
3. Run the script for your site with your preferred flags.
4. Open `starparty.html` on your phone/tablet at the telescope.
5. Enjoy your observing session!

---

## Deploy to a Web Server (Example: Daily at 07:00)

**Automation:** run this via a **daily cron job** and publish the HTML so the plan is always fresh before sunset.

Assume you want the output at `/var/www/star-party/index.html` on a Linux host (Nginx/Apache).

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
     --type_weights type_weights.csv \
     --top_n_per_hour 16 --out_prefix starparty \
     --html /var/www/star-party/index.html \
     --bsp ./skyfield_data/de440s.bsp \
     --html_ui tabs \
     --min_alt_planets 5 --min_alt_moon 0 \
     --preview_cache_dir /var/www/star-party/images \
     --preview_px 800 \
     --now_padding_min 30 \
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

   > If you omit `--date`, `--start`, and `--end`, the script picks **today**, rounds **sunset** for the start, and uses the **sunrise hour** as the end.

4. Visit `http(s)://your-server/star-party/` to see the latest plan.

---

## Troubleshooting

- **Missing previews**: Re‑run with `--refresh_previews`. Ensure `--preview_cache_dir` is writable by the user running the planner.
- **No BSP found**: Point `--bsp` to `./skyfield_data/de440s.bsp` (or your chosen kernel) and verify the file exists.
- **Tables seem empty**: Relax filters (`--min_alt`, `--moon_sep_min`, `--max_mag`) or expand the time window.
- **Wrong times**: Confirm your `--tz` and system clock. The planner expects an IANA timezone.
- **Permissions**: Make sure the web server user can read the generated HTML and the preview cache directory.

---

## Development Notes

- Written in **Python**, using Skyfield for accurate ephemerides
- All heavy astronomy data (BSP) is **local**; the planner does not auto‑fetch ephemerides
- UI is plain HTML/CSS/JS, tinted for night vision
- Easy to extend:
  - Adjust scoring in `interest_score()`
  - Add/modify catalogs
  - Tweak the UI (tabs vs accordion)

---

## License

MIT License

---

✨ *Clear skies!*
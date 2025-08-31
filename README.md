🪐 Star Party Planner

Plan an observing night with confidence.

Star Party Planner generates tailored observing lists for your site, date, and time window. It prioritizes crowd‑pleasers (Saturn, Jupiter, the Moon), filters deep‑sky targets by altitude, magnitude, and Moon separation, and produces night‑vision‑friendly HTML, CSV tables, and console listings.

A public demo for Minnesota Astronomical Society’s Eagle Lake Observatory:
	•	https://party.knowmad.org

⸻

Highlights
	•	Planets + Moon + custom deep‑sky catalog (CSV; Messier + Caldwell provided)
	•	Filters by minimum altitude, maximum magnitude, and Moon separation
	•	Assigns an interest score so crowd‑pleasers bubble to the top
	•	Finds each object’s best local time (peak altitude in your window)
	•	Night‑vision HTML with a red-on-black theme and mobile‑first layout
	•	Optional directional sector filters to block parts of the sky obscured by trees/buildings
	•	Clickable rows → modal with object details and a preview image (DSS2 or Wikipedia)
	•	Local preview cache so images aren’t re‑downloaded every run
	•	“Now” view that updates on the nearest 2‑minute tick without re‑running the planner

Outputs
	•	*_master.csv — Master observing list (by overall interest)
	•	*_hourly.csv — Hour‑by‑hour shortlist (top N each hour)
	•	*.html — Night‑vision web page (tabs; per‑hour chips)

⸻

Installation

Requires Python 3.10+.

git clone https://github.com/educationfutures/star-party-planner.git
cd star-party-planner

# Core dependencies
pip install skyfield numpy pandas pytz python-dateutil requests

Recommended versions
	•	Skyfield 1.45+
	•	NumPy 1.20+
	•	Pandas 1.3+

The HTML preview feature fetches DSS2 cutouts and (optionally) Wikipedia images — hence requests.

⸻

Ephemeris (JPL de440s.bsp)

The planner uses a local JPL ephemeris for precise positions. Recommended: de440s.bsp (~31 MB; valid 1550–2650 CE).

mkdir -p ./skyfield_data
curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
mv de440s.bsp ./skyfield_data/

*Result: ./skyfield_data/de440s.bsp*

The planner does not auto‑download BSP files. You must provide one via --bsp.

⸻

Catalogs

Provide a CSV with at least these columns:

name,ra_deg,dec_deg,type,mag,notes

Example:

M13,250.421,36.461,Globular cluster,5.8,"Great Hercules Cluster"

This repo includes:
	•	messier_caldwell.csv — Messier + Caldwell
	•	extended_targets.csv — additional crowd‑pleasers in N/S skies

⸻

Usage

Example (Eagle Lake Observatory, Minnesota)

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

First run may take longer if many preview images need to be cached.

⸻

Command‑Line Arguments

Core Location & Time
	•	--lat (float, required) — Latitude in decimal degrees (north positive)
	•	--lon (float, required) — Longitude in decimal degrees (east positive)
	•	--elev (float, default=0.0) — Elevation in meters
	•	--date (YYYY‑MM‑DD, default=today) — Observation date (local)
	•	--start (HH:MM, default=sunset rounded) — Session start time
	•	--end (HH:MM, default=sunrise hour) — Session end time
	•	--tz (string, default=UTC) — IANA timezone, e.g., America/Chicago

Catalog & Filtering
	•	--catalog (path, default=messier_caldwell.csv) — Deep‑sky catalog CSV
	•	--min_alt (float, default=20.0) — Minimum altitude (deg) for DSOs
	•	--max_mag (float, default=9.0) — Maximum magnitude (fainter excluded)
	•	--moon_sep_min (float, default=15.0) — Minimum Moon separation (deg) for DSOs
	•	--hour_step (int, default=1) — Step size in hours for hourly tables
	•	--top_n_per_hour (int, default=16) — Max targets shown per hour slot
	•	--moonlight_penalty_max (float, default=18) — Max points subtracted from diffuse targets at full Moon when the Moon is high

Outputs
	•	--out_prefix (string, default=starparty) — Prefix for CSV outputs
	•	--html (string, optional) — Save path for HTML output
	•	--bsp (path, default=./skyfield_data/de440s.bsp) — Planetary ephemeris file
	•	--html_ui (tabs | accordion, default=tabs) — HTML layout style

Per‑Type Altitude Thresholds
	•	--min_alt_planets (float, default=5.0) — Minimum altitude for planets
	•	--min_alt_moon (float, default=0.0) — Minimum altitude for the Moon

Previews & Caching (HTML)
	•	--no_previews — Disable previews in the HTML entirely
	•	--preview_cache_dir (string, default=image_cache) — Directory for cached previews
	•	--preview_px (int, default=800) — Pixel size for DSS2 images
	•	--preview_fov_deg (float, default=0.6) — Field of view for DSS2 images (degrees)
	•	--refresh_previews — Force re‑download of previews
	•	--clean_preview_cache — Remove cached previews not used by this run

Note: These cache flags apply only to HTML preview images. Catalogs and BSPs are unaffected.

⸻

HTML Output & UI
	•	Theme: red‑on‑black, designed for night use (keep device brightness low).
	•	Navbar: search field, tabs, and per‑hour “chips” (links) stay accessible.
	•	Views:
	•	Now — auto‑updates every ~2 minutes (client‑side), applies your saved filters, and shows the top N by score.
	•	Master List — overall targets sorted by interest.
	•	By Hour — click an hour chip to load that hour’s table; sector filters can hide rows and show a “Showing X of Y” indicator.
	•	Hidden columns: RA/Dec are hidden in tables but shown in the modal.

Row Modals
	•	Click a row for details (type, magnitude, best time, alt/az, direction, notes) and a preview image.
	•	Images are cached locally and rendered with a red tint to preserve night vision.

Filters dialog
	•	Open via the Filters button (funnel icon) in the navbar.
	•	Options:
	•	Enable filters (master on/off)
	•	Limit count of objects in “Now” view (Top N)
	•	Directional sector filters (up to 5 ranges): Each row = Azimuth Start°, End°, and Min Alt°. If current azimuth falls within a sector, the object must clear that sector’s minimum altitude to be shown.
	•	Settings persist in localStorage per browser.

Why no explicit “Refresh” button?
	•	The Now view aligns to 2‑minute time slots and auto‑refreshes. When the tab regains focus, it re-syncs to the current slot.

⸻

Preview Images & Sources
	•	DSS2 (Red) cutouts via CDS HiPS hips2fits service for most DSOs.
	•	Wikipedia/Wikimedia thumbnails for planets and as fallback when DSS2 is unavailable.
	•	Each image is cached under --preview_cache_dir and reused on subsequent runs.

Network hiccups don’t break the page: the modal still opens with object data; the image appears when available.

⸻

Example Workflow
	1.	Prepare or edit your catalog (messier_caldwell.csv).
	2.	Download de440s.bsp to ./skyfield_data/.
	3.	Run the planner for your site with preferred flags.
	4.	Open the generated HTML on your phone/tablet at the telescope.
	5.	Observe, share, tweak filters as needed.

Want it always ready before sunset? Automate it.

⸻

Deploy to a Web Server (Daily at 07:00)

Assume you’ll publish to /var/www/star-party/index.html on a Linux host.
	1.	Create the web directory and set permissions (adjust user/group as needed):

sudo mkdir -p /var/www/star-party
sudo chown -R $USER:www-data /var/www/star-party
sudo chmod -R 775 /var/www/star-party

	2.	Create ~/bin/run_starparty.sh:

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
  --preview_cache_dir /var/www/star-party/images \
  --clean_preview_cache

Make it executable:

chmod +x ~/bin/run_starparty.sh

	3.	Add a cron entry (07:00 daily):

crontab -e

Add:

0 7 * * * /home/youruser/bin/run_starparty.sh >> /home/youruser/starparty_cron.log 2>&1

When --date, --start, and --end are omitted, the planner computes today’s date, sunset, and the sunrise hour automatically.

⸻

Troubleshooting
	•	Blank / slow previews — Try a normal (non‑private) tab, or allow the request in a content blocker. Cached images will still show when present.
	•	No BSP found — Ensure --bsp points to a real file, e.g., ./skyfield_data/de440s.bsp.
	•	Empty tables — Relax filters (--min_alt, --moon_sep_min, --max_mag) or expand the time window.
	•	Wrong times — Check your --tz (IANA name) and system clock.
	•	Permissions — Ensure the web server can read the HTML and preview cache directory.

⸻

Development Notes
	•	Python + Skyfield for accurate ephemerides
	•	All heavy data (BSP) is local; the planner never auto‑fetches ephemerides
	•	UI is plain HTML/CSS/JS; red‑tinted for night vision
	•	Extendable:
	•	Adjust scoring in interest_score()
	•	Add new columns to your catalogs
	•	Swap UI layout (--html_ui tabs|accordion)

⸻

⚖️ License

MIT License

⸻

✨ Clear skies!
# Star Party Planner — Quick Reference

A concise cheat sheet for operators at the scope. For full details, see `README.md`.

---

## Quick Start

```bash
python starparty_planner.py \
  --lat 44.810265 --lon -93.939783 --elev 296 \
  --tz America/Chicago --catalog messier_caldwell.csv \
  --bsp ./skyfield_data/de440s.bsp \
  --html starparty.html
```

- If `--date/--start/--end` omitted: date=today, start≈sunset (rounded), end≈sunrise (hour).
- HTML opens with red-on-black UI; tables are searchable; rows open modal previews.

---

## Full Example (all common flags)

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
  --preview_cache_dir image_cache --preview_fov_deg 0.6 --preview_px 800
```

---

## Flags Reference (TL;DR)

### Site & Time
- `--lat FLOAT` (req) · Latitude in decimal degrees (N+)
- `--lon FLOAT` (req) · Longitude in decimal degrees (E+)
- `--elev FLOAT` · Elevation meters
- `--date YYYY-MM-DD` · Local date (default=today)
- `--start HH:MM` · Local start (default≈sunset rounded)
- `--end HH:MM` · Local end (default=hour of sunrise)
- `--tz IANA` · Timezone (e.g., `America/Chicago`)

### Catalog & Filters
- `--catalog FILE.csv` · DSO catalog with `name,ra_deg,dec_deg,type,mag,notes`
- `--min_alt DEG` · Min altitude for DSOs (default 20)
- `--max_mag MAG` · Max magnitude (fainter=larger) (default 9)
- `--moon_sep_min DEG` · Min separation from Moon (default 15)
- `--hour_step INT` · Hour step for per-hour tables (default 1)
- `--top_n_per_hour INT` · Cap list length per hour (default 16)

### Output
- `--out_prefix STR` · Prefix for CSVs (default `starparty`)
- `--html FILE.html` · Emit night-vision HTML (optional)
- `--bsp PATH` · Path to local ephemeris BSP (default `./skyfield_data/de440s.bsp`)
- `--html_ui {accordion|tabs}` · HTML layout (default `tabs`)

### Per-Type Altitude Thresholds
- `--min_alt_planets DEG` · Planet min altitude (default 5)
- `--min_alt_moon DEG` · Moon min altitude (default 0)

### Preview & Caching (HTML modal images)
- `--no_previews` · Disable all previews
- `--preview_cache_dir DIR` · Cache dir for images (default `image_cache`)
- `--preview_px INT` · Image pixel size (default 800)
- `--preview_fov_deg DEG` · DSS2 field of view (default 0.6)
- `--refresh_previews` · Force re-download images
- `--clean_preview_cache` · Delete cached images **not** used by this run

> Previews:
> - Planets/Moon → Wikipedia image (cached)
> - DSOs → DSS2 Red via hips2fits (cached), fallback to Wikipedia
> - All previews displayed with night-vision **red filter**

---

## Maintenance Tips

- **Ephemeris** is local-only; download `de440s.bsp` once to `./skyfield_data/`.
- To keep cache fresh on daily rebuilds:
  - Use `--refresh_previews` weekly or when images change.
  - Use `--clean_preview_cache` in daily cron to trim unused files.

### Cron Example (daily at 09:00 local)

```cron
0 9 * * * cd /path/to/star-party-planner && \
  /usr/bin/python3 starparty_planner.py \
  --lat 44.810265 --lon -93.939783 --elev 296 \
  --tz America/Chicago --catalog messier_caldwell.csv \
  --bsp ./skyfield_data/de440s.bsp \
  --html /var/www/html/starparty.html \
  --out_prefix /var/www/html/starparty \
  --top_n_per_hour 16 --min_alt 20 --moon_sep_min 20 --max_mag 9 \
  --min_alt_planets 5 --min_alt_moon 0 \
  --preview_cache_dir /var/www/html/image_cache \
  --preview_fov_deg 0.6 --preview_px 800 \
  --clean_preview_cache
```

---

## Troubleshooting

- **No objects listed** → Relax filters (`--min_alt`, `--moon_sep_min`, `--max_mag`) or extend time window.
- **No images in modal** → Ensure `--preview_cache_dir` is writable; try `--refresh_previews`.
- **White/washed images** → This UI forces a red filter on previews; device brightness still matters.
- **Time zone confusion** → Confirm `--tz` is a valid IANA name.

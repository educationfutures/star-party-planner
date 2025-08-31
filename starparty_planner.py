#!/usr/bin/env python3
"""
Star Party Planner — hourly observing list for a given site & night.

Features
- Planets + Moon + custom deep-sky catalogue (CSV)
- Filters by altitude, magnitude, and Moon separation
- Picks a "best time" (max altitude within your window)
- Prints:
  (1) A master list (sorted by interestingness, then time)
  (2) Hourly "point your scope now" tables (20:00, 21:00, ...)
- Exports CSV files for both views
- Optional night-vision HTML output (red on black) with tabs/accordions
- Clickable rows -> modal with preview image (cached on disk)

Usage (examples)
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

Dependencies
  pip install skyfield numpy pandas pytz python-dateutil requests

Notes
  - Uses a LOCAL ephemeris only (default: ./skyfield_data/de440s.bsp).
    Download once from NAIF and place it there:
      mkdir -p ./skyfield_data
      curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
      mv de440s.bsp ./skyfield_data/
  - Deep-sky objects come from a CSV you provide (e.g., messier_caldwell.csv).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import List, Dict, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from dateutil import tz
import requests

from skyfield.api import Loader, wgs84, load_file, Star
from skyfield import almanac

# Friendly UA for Wikipedia/Wikimedia calls
WIKI_UA = (
    "StarPartyPlanner/1.0 (+https://example.com; contact: you@example.com) "
    "requests"
)

# Disambiguation-safe page titles for planets
PLANET_TITLE_OVERRIDES = {
    "Mercury": "Mercury (planet)",
    "Venus": "Venus",
    "Earth": "Earth",
    "Moon": "Moon",
    "Mars": "Mars",
    "Jupiter": "Jupiter",
    "Saturn": "Saturn",
    "Uranus": "Uranus",
    "Neptune": "Neptune",
}

# ---------------------------- Config helpers ----------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Hourly observing list generator")
    p.add_argument("--lat", type=float, required=True, help="Latitude in decimal degrees (N+)")
    p.add_argument("--lon", type=float, required=True, help="Longitude in decimal degrees (E+)")
    p.add_argument("--elev", type=float, default=0.0, help="Elevation meters (optional)")

    # Defaults: date=today, start=sunset (rounded), end=sunrise (hour-only)
    p.add_argument("--date", type=str, default=None, help="Local date YYYY-MM-DD (default=today)")
    p.add_argument("--start", type=str, default=None, help="Local start time HH:MM (default=sunset rounded)")
    p.add_argument("--end", type=str, default=None, help="Local end time HH:MM (default=sunrise hour)")
    p.add_argument("--tz", type=str, default="UTC", help="IANA timezone, e.g., America/Chicago")

    p.add_argument("--catalog", type=str, default="messier_caldwell.csv", help="CSV of DSOs")
    p.add_argument("--min_alt", type=float, default=20.0, help="Minimum altitude for DSOs (deg)")
    p.add_argument("--max_mag", type=float, default=9.0, help="Max magnitude (fainter=larger) for DSOs")
    p.add_argument("--moon_sep_min", type=float, default=15.0, help="Min separation from Moon (deg)")
    p.add_argument("--hour_step", type=int, default=1, help="Hour step for hourly tables")
    p.add_argument("--top_n_per_hour", type=int, default=16, help="Max targets per hour table")

    p.add_argument("--out_prefix", type=str, default="starparty", help="Output prefix (csv)")
    p.add_argument("--html", type=str, default="", help="Optional HTML output filepath")
    p.add_argument("--bsp", type=str, default="./skyfield_data/de440s.bsp",
                   help="Path to local planetary ephemeris BSP (e.g., de440s.bsp)")
    p.add_argument("--html_ui", type=str, default="tabs", choices=["accordion", "tabs"],
                   help="HTML layout")

    # Per-type altitude thresholds
    p.add_argument("--min_alt_planets", type=float, default=5.0, help="Minimum altitude for planets (deg)")
    p.add_argument("--min_alt_moon", type=float, default=0.0, help="Minimum altitude for the Moon (deg)")

    # ---- Preview/caching options ----
    p.add_argument("--no_previews", action="store_true", help="Disable preview image generation entirely")
    p.add_argument("--preview_cache_dir", type=str, default="image_cache", help="Directory for cached preview images")
    p.add_argument("--preview_px", type=int, default=800, help="Pixel size for DSS2 hips2fits tiles")
    p.add_argument("--preview_fov_deg", type=float, default=0.6, help="FOV for DSS2 hips2fits (degrees)")
    p.add_argument("--refresh_previews", action="store_true", help="Force re-download of previews")
    p.add_argument("--clean_preview_cache", action="store_true", help="Remove cached files not used by this run")

    # ---- Moonlight penalty control ----
    p.add_argument("--moonlight_penalty_max", type=float, default=18.0,
                   help="Max points subtracted from diffuse targets at full Moon when the Moon is high (default=18)")

    return p.parse_args()

# ---------------------------- Astronomy core ----------------------------

def load_ephemeris(bsp_path: str):
    local_bsp = Path(bsp_path)
    if not local_bsp.exists():
        raise RuntimeError(
            f"Missing ephemeris BSP: {local_bsp}\n"
            "Download de440s.bsp (~31 MB) and point --bsp to it."
        )
    load = Loader(str(local_bsp.parent), expire=False)
    ts = load.timescale()
    eph = load_file(str(local_bsp))
    return load, ts, eph

@dataclass
class Target:
    name: str
    ra_deg: float
    dec_deg: float
    type: str
    mag: Optional[float] = None
    notes: str = ""

def read_catalog(path: str) -> List[Target]:
    targets: List[Target] = []
    for enc in ["utf-8", "latin-1"]:
        try:
            with open(path, newline="", encoding=enc, errors="replace") as f:
                r = csv.DictReader(f)
                for row in r:
                    try:
                        targets.append(Target(
                            name=row["name"].strip(),
                            ra_deg=float(row["ra_deg"]),
                            dec_deg=float(row["dec_deg"]),
                            type=row.get("type","").strip(),
                            mag=float(row["mag"]) if row.get("mag","").strip() else None,
                            notes=row.get("notes","").strip(),
                        ))
                    except Exception as e:
                        print(f"Skipping row due to error: {e} -> {row}")
            return targets
        except UnicodeDecodeError:
            continue
    return targets

def round_nearest_hour(dt: datetime) -> datetime:
    if dt.minute >= 30:
        return dt.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1)
    else:
        return dt.replace(minute=0, second=0, microsecond=0)

def fill_time_defaults(args, ts, eph):
    T = tz.gettz(args.tz)
    if args.date is None:
        args.date = datetime.now(T).strftime("%Y-%m-%d")
    day = datetime.strptime(args.date, "%Y-%m-%d").replace(tzinfo=T)

    site = wgs84.latlon(args.lat, args.lon, args.elev)
    t0 = ts.from_datetime(day.replace(hour=0, minute=0))
    t1 = ts.from_datetime(day.replace(hour=23, minute=59))
    f = almanac.sunrise_sunset(eph, site)
    times, events = almanac.find_discrete(t0, t1, f)

    sunrise, sunset = None, None
    for t, e in zip(times, events):
        dt_local = t.utc_datetime().astimezone(T)
        if e == 1: sunrise = dt_local
        elif e == 0: sunset = dt_local

    if args.start is None:
        args.start = round_nearest_hour(sunset).strftime("%H:%M") if sunset else "20:00"
    if args.end is None:
        args.end = sunrise.strftime("%H:00") if sunrise else "06:00"
    return args

def hours_list(local_date: str, start: str, end: str, tzname: str) -> List[datetime]:
    T = tz.gettz(tzname)
    day = datetime.strptime(local_date, "%Y-%m-%d").replace(tzinfo=T)
    s_hour, s_min = map(int, start.split(":"))
    e_hour, e_min = map(int, end.split(":"))
    start_dt = day.replace(hour=s_hour, minute=s_min)
    end_dt = day.replace(hour=e_hour, minute=e_min)
    if end_dt <= start_dt:
        end_dt += timedelta(days=1)
    out = []
    t = start_dt
    while t <= end_dt:
        out.append(t.replace(minute=0, second=0, microsecond=0))
        t += timedelta(hours=1)
    return out

def cardinal_from_az(az_deg: float) -> str:
    dirs = ["N","NNE","NE","ENE","E","ESE","SE","SSE","S","SSW","SW","WSW","W","WNW","NW","NNW"]
    idx = int((az_deg % 360) / 22.5 + 0.5) % 16
    return dirs[idx]

def angular_sep(ra1_deg, dec1_deg, ra2_deg, dec2_deg) -> float:
    r1, d1, r2, d2 = map(np.deg2rad, [ra1_deg, dec1_deg, ra2_deg, dec2_deg])
    return np.rad2deg(np.arccos(np.sin(d1)*np.sin(d2) + np.cos(d1)*np.cos(d2)*np.cos(r1-r2)))

def get_apparent(ts, observer, earth, target):
    return (earth + observer).at(ts).observe(target).apparent()

def to_altaz(ts, observer, earth, ra_deg: float, dec_deg: float):
    star = Star(ra_hours=ra_deg/15.0, dec_degrees=dec_deg)
    app = get_apparent(ts, observer, earth, star)
    alt, az, _ = app.altaz()
    return alt.degrees, az.degrees

from skyfield.api import Star as SfStar  # just for type hints

def moon_ra_dec(eph, ts, observer, earth) -> Tuple[float,float]:
    app = get_apparent(ts, observer, earth, eph['moon'])
    ra, dec, _ = app.radec()
    return ra.hours*15.0, dec.degrees

def planet_altaz(eph, ts, observer, earth, name: str) -> Tuple[float,float,float]:
    candidates_by_name = {
        "Mercury": ["mercury", "mercury barycenter"],
        "Venus":   ["venus", "venus barycenter"],
        "Mars":    ["mars barycenter", "mars"],
        "Jupiter": ["jupiter barycenter", "jupiter"],
        "Saturn":  ["saturn barycenter", "saturn"],
        "Uranus":  ["uranus barycenter", "uranus"],
        "Neptune": ["neptune barycenter", "neptune"],
    }
    target = None
    for key in candidates_by_name[name]:
        try:
            target = eph[key]
            break
        except KeyError:
            continue
    if target is None:
        raise KeyError(f"No suitable target in BSP for {name}")
    app = (earth + observer).at(ts).observe(target).apparent()
    alt, az, _ = app.altaz()
    return alt.degrees, az.degrees, None

def moon_altaz_phase(eph, ts, observer, earth) -> Tuple[float,float,float]:
    app = get_apparent(ts, observer, earth, eph['moon'])
    alt, az, _ = app.altaz()
    phase = almanac.moon_phase(eph, ts).degrees
    return alt.degrees, az.degrees, phase

def best_time_in_window(ts_arr, alts) -> int:
    return int(np.nanargmax(alts))

def build_observer(load, ts, lat, lon, elev):
    return wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=elev)

# ---------------------------- Interest scoring ----------------------------

INTEREST_BASE = {
    "Saturn": 100, "Jupiter": 90, "Moon": 70, "Mars": 80, "Venus": 70,
    "Uranus": 55, "Neptune": 50,
}
TYPE_BONUS = {
    "Globular cluster": 45, "Open cluster": 35, "Planetary nebula": 40, "Galaxy cluster": 32,
    "Emission Nebula": 40, "Reflection Nebula": 35, "Dark nebula": 30, "Nebula with cluster": 45,
    "H II region nebula with cluster": 40, "Spiral galaxy": 35, "Peculiar galaxy": 30, 
    "Elliptical galaxy": 30, "Starburst galaxy": 30, "Galaxy": 30, "Lenticular galaxy": 32,
    "Supernova Remnant": 42, "Milky Way star cloud": 38, "Quasar": 25,
    "Asterism": 50, "Optical Double": 50, "Double star": 50, "Multiple star": 50,
}
CROWD_BONUS = {"Saturn": 30, "Jupiter": 25, "Moon": 20, "Mars": 12, "Venus": 10, "Cr 399": 15, "Albireo": 15}

def interest_score(name: str, typ: str, best_alt: float, alt_now: Optional[float] = None) -> float:
    base = INTEREST_BASE.get(name, 0)
    if base == 0:
        base = TYPE_BONUS.get(typ, 25)
    base += CROWD_BONUS.get(name, 0)
    alt_term = 0.10 * best_alt + 0.15 * (alt_now if alt_now is not None else best_alt)
    return min(200.0, base + alt_term)

# ---- Moonlight penalty helpers ----

# Targets that visibly wash out under moonlight
DIFFUSE_TYPES = {
    "Emission Nebula", "Reflection Nebula", "Dark nebula",
    "Milky Way star cloud", "Nebula with cluster", "H II region nebula with cluster",
    "Spiral galaxy", "Elliptical galaxy", "Lenticular galaxy", "Starburst galaxy",
    "Galaxy", "Peculiar galaxy", "Galaxy cluster", "Supernova Remnant"
}
# Compact nebulae: penalize half-strength
PARTIALLY_DIFFUSE_TYPES = {"Planetary nebula"}

def lunar_illum_fraction_from_phase_deg(phase_deg: float) -> float:
    """Convert Moon-Sun elongation (0°=New, 180°=Full) to illumination fraction [0,1]."""
    return 0.5 * (1.0 - math.cos(math.radians(phase_deg)))

def moon_alt_scale(alt_deg: Optional[float]) -> float:
    """
    Scale penalty by Moon altitude using sin(alt).
    - Below horizon -> 0
    - At horizon -> 0
    - 30° -> 0.5
    - 90° -> 1.0
    """
    if alt_deg is None:
        return 1.0
    return max(0.0, math.sin(math.radians(max(-5.0, float(alt_deg)))))  # small clamp

def moonlight_penalty_points(obj_type: str, illum_frac: float, max_penalty: float,
                             moon_alt_deg: Optional[float] = None) -> float:
    """
    Points to subtract based on Moon illumination and altitude.
    Illumination sets the ceiling; altitude scales it by visibility.
    """
    t = (obj_type or "").strip()
    if t in DIFFUSE_TYPES:
        return illum_frac * max_penalty * moon_alt_scale(moon_alt_deg)
    if t in PARTIALLY_DIFFUSE_TYPES:
        return 0.5 * illum_frac * max_penalty * moon_alt_scale(moon_alt_deg)
    return 0.0
def _round3(v):
    return "" if pd.isna(v) else round(float(v), 3)

# ---------------------------- Preview caching ----------------------------

def slugify(s: str) -> str:
    s = re.sub(r"[^\w\-]+", "_", s.strip(), flags=re.U)
    s = re.sub(r"_+", "_", s)
    return s.strip("_").lower() or "object"

def cache_write(path: Path, content: bytes, refresh: bool):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not refresh:
        return
    path.write_bytes(content)

def fetch_wikipedia_image(name: str, cache_dir: Path, refresh: bool,
                          title_override: Optional[str] = None) -> Optional[Path]:
    """
    Try to download a Wikipedia/Wikimedia image for `name`.
    Strategy:
      A) Wikipedia REST summary (thumbnail/originalimage)
      B) MediaWiki API 'pageimages'
    """
    slug = slugify(name) + "_wiki.jpg"
    out = cache_dir / slug
    if out.exists():
        try:
            if out.stat().st_size < 1024:
                refresh = True
        except Exception:
            refresh = True
    if out.exists() and not refresh:
        return out

    headers = {"accept": "application/json", "User-Agent": WIKI_UA}
    title = title_override or name

    try:
        summary_url = f"https://en.wikipedia.org/api/rest_v1/page/summary/{requests.utils.quote(title)}"
        r = requests.get(summary_url, headers=headers, timeout=20, allow_redirects=True)
        if r.ok:
            j = r.json()
            img = (j.get("originalimage") or {}).get("source") or (j.get("thumbnail") or {}).get("source")
            if img:
                ir = requests.get(img, headers={"User-Agent": WIKI_UA}, timeout=30, allow_redirects=True)
                if ir.ok and ir.content:
                    cache_write(out, ir.content, True)
                    return out
    except Exception:
        pass

    try:
        api_url = "https://en.wikipedia.org/w/api.php"
        params = {
            "action": "query", "format": "json", "prop": "pageimages",
            "piprop": "original|thumbnail", "pithumbsize": "1000", "redirects": "1",
            "titles": title,
        }
        r = requests.get(api_url, params=params, headers=headers, timeout=25, allow_redirects=True)
        if r.ok:
            j = r.json()
            pages = j.get("query", {}).get("pages", {})
            for _, page in pages.items():
                src = None
                if isinstance(page, dict):
                    if "original" in page and "source" in page["original"]:
                        src = page["original"]["source"]
                    elif "thumbnail" in page and "source" in page["thumbnail"]:
                        src = page["thumbnail"]["source"]
                if src:
                    ir = requests.get(src, headers={"User-Agent": WIKI_UA}, timeout=30, allow_redirects=True)
                    if ir.ok and ir.content:
                        cache_write(out, ir.content, True)
                        return out
    except Exception:
        pass
    return None

def fetch_dss2_hips2fits(ra_deg: float, dec_deg: float, fov_deg: float, px: int,
                         cache_dir: Path, name: str, refresh: bool) -> Optional[Path]:
    """Download DSS2 Red cutout via hips2fits and cache it."""
    if any(map(lambda v: v is None or (isinstance(v, float) and math.isnan(v)), [ra_deg, dec_deg])):
        return None
    slug = f"{slugify(name)}_dss2red_{fov_deg:.3f}deg_{px}px.jpg"
    out = cache_dir / slug
    if out.exists() and not refresh:
        return out
    base = "https://alasky.u-strasbg.fr/hips-image-services/hips2fits"
    params = {
        "hips": "CDS/P/DSS2/red",
        "ra": f"{ra_deg:.6f}", "dec": f"{dec_deg:.6f}",
        "fov": f"{fov_deg:.6f}", "width": str(px), "height": str(px),
        "format": "jpg", "projection": "TAN",
    }
    try:
        r = requests.get(base, params=params, timeout=60)
        if r.status_code != 200:
            return None
        cache_write(out, r.content, refresh)
        return out
    except Exception:
        return None

def build_preview_cache(master_df: pd.DataFrame, cache_dir: Path, refresh: bool,
                        fov_deg: float, px: int) -> Dict[str, Dict[str, str]]:
    """
    Returns map: name -> { "path": "image_cache/xxx.jpg", "kind": "dss2|wiki" }
    Policy:
      - Planets + Moon -> Wikipedia
      - DSOs -> DSS2 first, fallback to Wikipedia if DSS2 fails
    """
    preview: Dict[str, Dict[str, str]] = {}
    if master_df.empty:
        return preview

    planet_like = {"Mercury","Venus","Earth","Moon","Mars","Jupiter","Saturn","Uranus","Neptune"}

    for _, row in master_df.iterrows():
        name = str(row.get("Name", "")).strip()
        if not name or name in preview:
            continue
        ra = row.get("RA (deg)")
        dec = row.get("Dec (deg)")
        typ = str(row.get("Type", "")).strip().lower()

        chosen: Optional[Path] = None
        kind = None

        if (name in planet_like) or (typ in {"planet", "moon"}):
            title = PLANET_TITLE_OVERRIDES.get(name, name)
            chosen = fetch_wikipedia_image(name, cache_dir, refresh, title_override=title)
            kind = "wiki" if chosen else None
        else:
            chosen = fetch_dss2_hips2fits(
                float(ra) if pd.notna(ra) else None,
                float(dec) if pd.notna(dec) else None,
                fov_deg, px, cache_dir, name, refresh
            )
            if chosen:
                kind = "dss2"
            else:
                chosen = fetch_wikipedia_image(name, cache_dir, refresh)
                kind = "wiki" if chosen else None

        if chosen and kind:
            rel = str(Path(cache_dir.name) / chosen.name)
            preview[name] = {"path": rel, "kind": kind}

    return preview

def cleanup_preview_cache(cache_dir: Path, preview_map: Dict[str, Dict[str, str]]):
    """Delete cached files not referenced by this run's preview_map."""
    if not cache_dir.exists():
        return
    keep = {Path(v["path"]).name for v in preview_map.values()}
    for p in cache_dir.glob("*"):
        if p.is_file() and p.name not in keep:
            try: p.unlink()
            except Exception: pass

# ---------------------------- HTML rendering ----------------------------

def df_to_html_table(df: pd.DataFrame, id_attr: str = "") -> str:
    if df.empty:
        return "<p>No data.</p>"
    return f'<div class="table-wrap">{df.to_html(index=False, escape=True, border=0, table_id=id_attr)}</div>'

def _hour_anchor_label(hour_str: str) -> str:
    try:
        t = pd.to_datetime(hour_str)
        return t.strftime("%H:%M")
    except Exception:
        return hour_str

def write_html(output_path: str, site_lat: float, site_lon: float, tzname: str, date_str: str,
               start: str, end: str, master_df: pd.DataFrame, hourly_df: pd.DataFrame,
               ui_mode: str, preview_map: Dict[str, Dict[str,str]], now_data: List[Dict]):
    # ---------------- Build hourly lookup + chips ----------------
    hours_lookup = {}
    hour_links = []
    if not hourly_df.empty:
        for hour, sub in hourly_df.groupby("Hour"):
            anchor = f"hour-{hour.replace(' ','_').replace(':','')}"
            hour_label = _hour_anchor_label(hour)
            hour_links.append(f'<a href="#{anchor}">{hour_label}</a>')
            cols = ['Name','Type','Mag','Alt (°)','Az (°)','Dir','RA (deg)','Dec (deg)','Notes']
            present = [c for c in cols if c in sub.columns]
            tbl = df_to_html_table(sub[present], id_attr=f"tbl-{anchor}")
            hours_lookup[anchor] = f"<h3>{hour}</h3>{tbl}"
    master_html = df_to_html_table(master_df, id_attr="tbl-master") if not master_df.empty else "<p>No targets passed the filters. Try adjusting filters.</p>"

    # ---------------- Night-vision CSS ----------------
    css = """
    <style>
      :root { color-scheme: dark; }
      html, body { background: #000; color: #f33; }
      body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Ubuntu, Cantarell, 'Helvetica Neue', Arial, sans-serif; line-height: 1.35; }
      ::selection { background: #400; color: #fdd; }
      a { color: #f66; text-decoration: underline; }
      h1,h2,h3 { color: #f44; margin: 0.8rem 0 0.4rem; }
      .container { max-width: 1100px; margin: 0 auto; padding: 1rem; }

      .hr { border: 0; height: 1px; background: #400; margin: 1rem 0; }
      .small { font-size: 0.9rem; color: #f66; }
      .warn { color:#f77; font-style: italic; }

      .table-wrap { overflow-x: auto; -webkit-overflow-scrolling: touch; }
      table { width: 100%; border-collapse: collapse; margin: 0.5rem 0 1rem; }
      th, td { border: 1px solid #700; padding: 0.45rem 0.5rem; }
      th { background: #100; text-align: center; }
      tr:nth-child(even) { background: #070707; }
      .table-wrap table tbody tr { cursor: pointer; }
      .table-wrap table tbody tr:hover { background: #111; }
      th.hide, td.hide { display: none; }

      .meta-bar { position: static; padding: 0.5rem 0; display:flex; gap:0.6rem; align-items:center; flex-wrap:wrap; }
      .pill { border:1px solid #700; padding:0.25rem 0.5rem; border-radius:999px; color:#f66; }

      /* Sticky nav containing search, tabs, and hour chips */
      .navbar { position: sticky; top: 0; background: rgba(0,0,0,0.98); border-bottom: 1px solid #400; padding: 0.5rem; z-index: 5; display:flex; flex-wrap:wrap; gap:0.6rem; align-items:center; }
      .navbar .left { display:flex; gap:0.6rem; align-items:center; flex-wrap:wrap; }
      .navbar .right { flex: 1 1 100%; display:flex; gap:0.4rem; flex-wrap:wrap; align-items:center; margin-top: 0.35rem; }

      input[type="search"] { background: #160000; border: 1px solid #700; color: #f55; padding: 0.4rem 0.6rem; border-radius: 6px; min-width: 220px; caret-color: #f55; outline: none; }
      input[type="search"]::placeholder { color:#f66; }
      input[type="search"]:focus { box-shadow: 0 0 0 2px #500 inset; border-color:#900; }

      .tabs { display:flex; gap:0.4rem; }
      .tab { padding:0.35rem 0.7rem; border:1px solid #500; border-radius:8px; cursor:pointer; user-select:none; color:#f66; background:#0a0000; }
      .tab.active { background:#180000; border-color:#700; color:#f66; }

      .hours { display:flex; gap:0.35rem; flex-wrap:wrap; }
      .hours a { display:inline-block; padding: 0.25rem 0.55rem; border:1px solid #500; border-radius:6px; text-decoration:none; }

      .hidden { display:none; }

      .modal-backdrop { position: fixed; inset: 0; background: rgba(0,0,0,0.75); display: none; align-items: center; justify-content: center; z-index: 50; }
      .modal {
        width: min(900px, 92vw);
        max-height: 90vh; overflow: auto;
        background: #050505; border: 1px solid #600; border-radius: 12px; padding: 1rem;
        color: #f55; box-shadow: 0 0 0 2px #300;
      }
      .modal h3 { margin: 0 0 .5rem; color: #f66; }
      .modal .meta { display: grid; grid-template-columns: repeat(2, minmax(0,1fr)); gap: .5rem .75rem; margin-bottom: .75rem; }
      .modal .imgbox { text-align: center; margin: .5rem 0 0; }
      .modal .imgbox img { max-width: 100%; border: 1px solid #400; border-radius: 8px; }
      .modal .close {
        float: right; border:1px solid #700; background:#180000; color:#f66; border-radius:8px; padding:.25rem .6rem; cursor:pointer;
      }

      .night-red { filter: sepia(1) saturate(6) hue-rotate(-50deg) brightness(0.9) contrast(1.1); }

      /* ---- Filters modal ---- */
        .filters-grid {
        display: grid;
        grid-template-columns: 1fr 1fr 1fr;
        gap: .5rem .75rem;
        margin: .75rem 0;
        }
        .filters-grid label { font-size:.9rem; color:#f66; }
        .filters-row {
        display: grid;
        grid-template-columns: repeat(3, minmax(0,1fr));
        gap: .5rem .75rem;
        margin-bottom: .4rem;
        }
        .filters-row input[type="number"],
        .filters-grid input[type="number"] {
        background:#160000; border:1px solid #700; color:#f55;
        padding:.35rem .5rem; border-radius:6px; width:100%;
        }
        .filters-actions { display:flex; gap:.5rem; justify-content:flex-end; margin-top:.75rem; }
        .filters-actions button {
        border:1px solid #700; background:#180000; color:#f66;
        border-radius:8px; padding:.35rem .7rem; cursor:pointer;
        }
        .filters-switch { display:flex; align-items:center; gap:.5rem; margin:.25rem 0 .5rem; }
        .filters-note { font-size:.85rem; color:#f77; }
        .filters-dialog label {
        display: block;
        font-size: 0.85rem;
        margin-bottom: 0.2rem;
        color: #f77;
        }

        .filters-dialog input[type="number"],
        .filters-dialog input[type="checkbox"],
        .filters-dialog select {
        background: #100;
        border: 1px solid #700;
        color: #f55;
        font-size: 0.85rem;
        padding: 0.25rem 0.4rem;
        border-radius: 4px;
        width: 5rem;
        text-align: center;
        }

        .filters-dialog .sector-row {
        display: flex;
        gap: 0.5rem;
        margin-bottom: 0.35rem;
        }

        .filters-dialog .sector-row input {
        flex: 1 1 auto;
        max-width: 5rem;
        }

        .filters-dialog .topn-row {
        margin-bottom: 0.8rem;
        display: flex;
        align-items: center;
        gap: 0.6rem;
        }

        .filters-dialog .buttons {
        display: flex;
        justify-content: flex-end;
        gap: 0.6rem;
        margin-top: 0.6rem;
        }

        .filters-dialog button {
        border: 1px solid #700;
        background: #180000;
        color: #f66;
        border-radius: 6px;
        padding: 0.25rem 0.7rem;
        cursor: pointer;
        font-size: 0.85rem;
        }

      @media (max-width: 640px) {
        body { font-size: 15px; }
        th, td { padding: 0.35rem 0.45rem; }
        .pill { font-size: 0.85rem; }
        .tab { padding:0.3rem 0.55rem; }
        .modal .meta{ grid-template-columns: 1fr; }
      }
    </style>
    """

    # ---------------- JS payload ----------------
    preview_json = json.dumps(preview_map)
    now_json = json.dumps(now_data)
    hours_json = json.dumps(hours_lookup)

    js_template = """
    <script>
    /* ======= injected data ======= */
    const PREVIEW_MAP = __PREVIEW_JSON__;
    const NOW_DATA    = __NOW_JSON__;
    const HOURS_HTML  = __HOURS_JSON__;

    /* ==============================
    Local storage + sector logic
    ============================== */
    const LS_KEYS = {
    enabled: 'obs.enabled',
    sectors: 'obs.sectors',
    topN:    'now.topN',
    };

    function _clamp(v, lo, hi) { v = Number(v); return isFinite(v) ? Math.min(Math.max(v, lo), hi) : NaN; }

    function loadSettings(defaultTopN = 16) {
    const enabled = (localStorage.getItem(LS_KEYS.enabled) || 'false') === 'true';
    let sectors = [];
    try { sectors = JSON.parse(localStorage.getItem(LS_KEYS.sectors) || '[]'); } catch { sectors = []; }
    let topN = parseInt(localStorage.getItem(LS_KEYS.topN) || `${defaultTopN}`, 10);
    if (!Number.isFinite(topN) || topN <= 0) topN = defaultTopN;

    // sanitize
    const clean = [];
    for (const s of sectors) {
        const az1 = _clamp(s.az1, 0, 360);
        const az2 = _clamp(s.az2, 0, 360);
        const minAlt = _clamp(s.minAlt, 0, 90);
        if (Number.isFinite(az1) && Number.isFinite(az2) && Number.isFinite(minAlt)) clean.push({ az1, az2, minAlt });
    }
    return { enabled, sectors: clean.slice(0, 5), topN };
    }

    function saveSettings({ enabled, sectors, topN }) {
    localStorage.setItem(LS_KEYS.enabled, enabled ? 'true' : 'false');
    localStorage.setItem(LS_KEYS.sectors, JSON.stringify((sectors || []).slice(0, 5)));
    localStorage.setItem(LS_KEYS.topN, String(topN));
    }

    // ---------- Sector math ----------
    function normAz(x)    { x = Number(x); return ((x % 360) + 360) % 360; }
    function azInSector(az, a1, a2) {
    az = normAz(az); a1 = normAz(a1); a2 = normAz(a2);
    if (a1 === a2) return true;               // treat 0-length as full circle
    return (a1 <= a2) ? (az >= a1 && az <= a2)
                        : (az >= a1 || az <= a2); // wrap-around
    }
    function requiredAltForAz(sectors, az) {
    let req = 0;
    for (const s of sectors) if (azInSector(az, s.az1, s.az2)) req = Math.max(req, s.minAlt);
    return req;
    }
    function passesSectors(row, sectors) {
    const az = Number(row["Az (°)"]);
    const alt = Number(row["Alt (°)"]);
    if (!Number.isFinite(az) || !Number.isFinite(alt)) return true;
    const need = requiredAltForAz(sectors, az);
    return alt >= need;
    }

    /* ==============================
    UI helpers
    ============================== */
    function filterTables(q) {
    const needle = q.trim().toLowerCase();
    document.querySelectorAll('table').forEach(tbl => {
        const rows = tbl.tBodies[0]?.rows || [];
        for (let r of rows) {
        const txt = r.innerText.toLowerCase();
        r.style.display = (!needle || txt.includes(needle)) ? '' : 'none';
        }
    });
    }

    function activateMainTab(targetId) {
    document.querySelectorAll('.tabs .tab').forEach(el => el.classList.remove('active'));
    document.querySelectorAll('.tab-panel').forEach(el => el.classList.add('hidden'));
    const tabBtn = document.querySelector('.tab[data-tab="' + targetId + '"]');
    if (tabBtn) tabBtn.classList.add('active');
    const panel = document.getElementById(targetId);
    if (panel) panel.classList.remove('hidden');
    window.scrollTo({ top: 0, behavior: 'instant' });
    }

    function initTabsBehavior() {
    const tabs = document.querySelectorAll('[data-tab]');
    tabs.forEach(tab => {
        tab.addEventListener('click', () => {
        activateMainTab(tab.getAttribute('data-tab'));
        });
    });
    }

    function hideDataColumns() {
    const toHide = new Set(["RA (deg)", "Dec (deg)"]);
    document.querySelectorAll('table').forEach(tbl => {
        const head = tbl.tHead?.rows?.[0];
        if (!head) return;
        const idxs = [];
        [...head.cells].forEach((th, i) => {
        const label = (th.innerText || "").trim();
        if (toHide.has(label)) { th.classList.add('hide'); idxs.push(i); }
        });
        tbl.querySelectorAll('tbody tr').forEach(tr => {
        idxs.forEach(i => tr.cells[i]?.classList.add('hide'));
        });
    });
    }

    /* ==============================
    Preview modal
    ============================== */
    function buildMetaHTML(fields) {
    const items = [];
    const push = (label, val) => { if (val !== undefined && val !== "" && val != null && val !== "nan") items.push(`<div><strong>${label}:</strong> ${val}</div>`); };
    push("Type", fields.type); push("Magnitude", fields.mag);
    push("Best Time", fields.bestTime); push("Dir", fields.dir);
    push("Alt (°)", fields.alt); push("Az (°)", fields.az);
    push("RA (deg)", fields.ra); push("Dec (deg)", fields.dec);
    push("Notes", fields.notes);
    return items.join("");
    }

    function openModalForRow(tr) {
    const tbl = tr.closest('table');
    const headers = [...(tbl.tHead?.rows?.[0]?.cells || [])].map(th => (th.innerText || "").trim());
    const get = (label) => { const i = headers.indexOf(label); return (i >= 0 && tr.cells[i]) ? tr.cells[i].innerText.trim() : ""; };

    const data = {
        name: get("Name"),
        type: get("Type"),
        mag:  get("Mag"),
        bestTime: get("Best Time") || get("Best Local Time"),
        dir:  get("Dir"),
        alt:  get("Alt (°)"),
        az:   get("Az (°)"),
        ra:   get("RA (deg)"),
        dec:  get("Dec (deg)"),
        notes:get("Notes"),
    };

    const backdrop = document.getElementById('modal-backdrop');
    const title = document.getElementById('modal-title');
    const meta = document.getElementById('modal-meta');
    const imgbox = document.getElementById('modal-image');

    title.textContent = data.name || "Object";
    meta.innerHTML = buildMetaHTML(data);
    imgbox.innerHTML = `<p class="small">Loading preview…</p>`;
    backdrop.style.display = 'flex';
    backdrop.setAttribute('aria-hidden', 'false');

    const entry = PREVIEW_MAP[data.name] || null;
    if (!entry) { imgbox.innerHTML = `<p class="small">No preview available.</p>`; return; }
    imgbox.innerHTML =
        '<img class="night-red" src="' + entry.path + '" alt="Preview of ' + data.name + '" loading="lazy">' +
        '<div class="small" style="margin-top:.35rem;">' +
        (entry.kind === "dss2" ? "Image: DSS2 Red (CDS hips2fits)" : "Image: Wikipedia/Wikimedia — red-tinted") +
        "</div>";
    }

    function closeModal() {
    const backdrop = document.getElementById('modal-backdrop');
    backdrop.style.display = 'none';
    backdrop.setAttribute('aria-hidden', 'true');
    }

    function enableRowModals() {
    document.querySelectorAll('.table-wrap table').forEach(tbl => {
        tbl.addEventListener('click', (e) => {
        const tr = e.target.closest('tr');
        if (!tr || tr.parentElement.tagName !== 'TBODY') return;
        openModalForRow(tr);
        }, { passive: true });
    });
    const backdrop = document.getElementById('modal-backdrop');
    backdrop?.addEventListener('click', (e) => {
        if (e.target.id === 'modal-backdrop' || e.target.classList.contains('close')) closeModal();
    });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closeModal(); });
    }

    /* ==============================
    Hour chips -> single hour view
    ============================== */
    function handleHourNavClick(e) {
    const href = e.currentTarget.getAttribute('href') || '';
    if (!href.startsWith('#')) return;
    e.preventDefault();
    const id = href.slice(1);
    activateMainTab('panel-hour');
    const wrap = document.getElementById('hour-panel-wrap');
    if (wrap) {
        wrap.innerHTML = HOURS_HTML[id] || '<p class="small">No data for that hour.</p>';
        enableRowModals();
        hideDataColumns();
        applyHourFilters();   // apply sector filters + show counts
    }
    }

    function initHourLinks() {
    document.querySelectorAll('.navbar .hours a').forEach(a => {
        a.addEventListener('click', handleHourNavClick);
    });
    }

    // Hide hour rows that don't meet sector min-alt, and show counts
    function applyHourFilters() {
    const panel = document.getElementById('panel-hour');
    if (!panel || panel.classList.contains('hidden')) return; // only if hour panel visible
    const tbl = panel.querySelector('table');
    if (!tbl) return;

    const s = loadSettings(16);
    const rows = Array.from(tbl.tBodies?.[0]?.rows || []);
    let total = rows.length;
    let shown = 0;

    // headers once
    const headers = [...(tbl.tHead?.rows?.[0]?.cells || [])].map(th => (th.innerText || '').trim());
    const idxAlt = headers.indexOf("Alt (°)");
    const idxAz  = headers.indexOf("Az (°)");

    for (const tr of rows) {
        if (!s.enabled || !s.sectors.length) {
        tr.style.display = '';
        shown++;
        continue;
        }
        const rowObj = {
        "Az (°)": idxAz >= 0 ? tr.cells[idxAz]?.innerText?.trim() : "",
        "Alt (°)": idxAlt >= 0 ? tr.cells[idxAlt]?.innerText?.trim() : "",
        };
        const ok = passesSectors(rowObj, s.sectors);
        tr.style.display = ok ? '' : 'none';
        if (ok) shown++;
    }

    // Add/Update a tiny “Showing X of Y” indicator under the hour title
    let info = panel.querySelector('.hour-count-info');
    if (!info) {
        info = document.createElement('div');
        info.className = 'hour-count-info small';
        const before = panel.querySelector('.table-wrap');
        if (before) panel.insertBefore(info, before);
        else panel.appendChild(info);
    }
    info.textContent = s.enabled && s.sectors.length ? `Showing ${shown} of ${total}` : `Showing ${total} of ${total}`;
    }

    /* ==============================
    NOW view (2-minute granularity)
    ============================== */
    const NOW_INDEX = (() => {
    const m = new Map();
    for (const entry of NOW_DATA || []) {
        if (entry && entry.time && Array.isArray(entry.rows)) {
        m.set(entry.time, entry.rows);
        }
    }
    return m;
    })();
    const NOW_KEYS = Array.from(NOW_INDEX.keys()).sort();

    function _fmt2(n) { return n < 10 ? "0" + n : "" + n; }
    function _keyToDate(keyStr) { return new Date(keyStr.replace(" ", "T") + ":00"); }
    function _currentTwoMinuteKeyRaw() {
    const d = new Date();
    const y = d.getFullYear(), mo = _fmt2(d.getMonth()+1), da = _fmt2(d.getDate());
    const hh = _fmt2(d.getHours()); const mm = d.getMinutes() - (d.getMinutes()%2);
    return `${y}-${mo}-${da} ${hh}:${_fmt2(mm)}`;
    }

    function rowsToHTML(rows) {
    if (!rows || rows.length === 0) return '<p class="small">No targets above your thresholds right now.</p>';
    const head = `
        <table id="tbl-now">
        <thead>
            <tr>
            <th>Name</th><th>Type</th><th>Mag</th>
            <th>Alt (°)</th><th>Az (°)</th><th>Dir</th><th>Notes</th>
            </tr>
        </thead>
        <tbody>`;
    const body = rows.map(r => `
        <tr>
        <td>${r["Name"] ?? ""}</td>
        <td>${r["Type"] ?? ""}</td>
        <td>${r["Mag"] ?? ""}</td>
        <td>${r["Alt (°)"] ?? ""}</td>
        <td>${r["Az (°)"] ?? ""}</td>
        <td>${r["Dir"] ?? ""}</td>
        <td>${r["Notes"] ?? ""}</td>
        </tr>`).join("");
    return '<div class="table-wrap">' + head + body + "</tbody></table></div>";
    }

    function renderNow(forceKey = null) {
    const wrap = document.getElementById("now-table-wrap");
    if (!wrap) return;
    if (!NOW_KEYS.length) { wrap.innerHTML = '<p class="small">No current data in range.</p>'; return; }

    const firstKey = NOW_KEYS[0], lastKey = NOW_KEYS[NOW_KEYS.length - 1];
    const firstDt = _keyToDate(firstKey), lastDt = _keyToDate(lastKey);

    let key = forceKey || _currentTwoMinuteKeyRaw();
    const keyDt = _keyToDate(key);

    if (keyDt < firstDt || keyDt > lastDt) {
        wrap.innerHTML = '<p class="small warn">The current time is outside of the observation window.</p>';
        return;
    }

    if (!NOW_INDEX.has(key)) {
        let fallback = firstKey;
        for (const k of NOW_KEYS) { if (_keyToDate(k) <= keyDt) fallback = k; else break; }
        key = fallback;
    }

    const settings = loadSettings(16);
    let rows = NOW_INDEX.get(key) || [];

    // Apply sector filters first
    if (settings.enabled && settings.sectors.length) {
        rows = rows.filter(r => passesSectors(r, settings.sectors));
    }

    // Sort by _Score desc, then Name, then slice to TopN
    rows.sort((a, b) => (Number(b._Score||0) - Number(a._Score||0)) || String(a.Name).localeCompare(String(b.Name)));
    rows = rows.slice(0, settings.topN);

    wrap.innerHTML = rowsToHTML(rows);
    enableRowModals();
    hideDataColumns();
    }

    /* ==============================
    Filters dialog UI
    ============================== */
    function openFilters() {
    const b = document.getElementById('filters-backdrop');
    if (!b) return;
    const s = loadSettings(16);
    document.getElementById('obs-enabled').checked = !!s.enabled;
    const topNEl = document.getElementById('now-topn');
    if (topNEl) topNEl.value = s.topN;

    const rows = document.querySelectorAll('#sector-rows .filters-row');
    for (let i = 0; i < rows.length; i++) {
        const r = rows[i];
        const sct = s.sectors[i] || {};
        const az1 = r.querySelector('.az1'); const az2 = r.querySelector('.az2'); const ma = r.querySelector('.minAlt');
        if (az1) az1.value = (sct.az1 ?? "");
        if (az2) az2.value = (sct.az2 ?? "");
        if (ma)  ma.value  = (sct.minAlt ?? "");
    }

    b.style.display = 'flex';
    b.setAttribute('aria-hidden', 'false');
    }
    function closeFilters() {
    const b = document.getElementById('filters-backdrop');
    if (!b) return;
    b.style.display = 'none';
    b.setAttribute('aria-hidden', 'true');
    }
    function readSettingsFromForm() {
    const enabled = !!document.getElementById('obs-enabled')?.checked;
    const topN = _clamp(document.getElementById('now-topn')?.value, 1, 64) || 16;
    const out = [];
    const rows = document.querySelectorAll('#sector-rows .filters-row');
    rows.forEach(r => {
        const az1 = _clamp(r.querySelector('.az1')?.value, 0, 360);
        const az2 = _clamp(r.querySelector('.az2')?.value, 0, 360);
        const minAlt = _clamp(r.querySelector('.minAlt')?.value, 0, 90);
        if (Number.isFinite(az1) && Number.isFinite(az2) && Number.isFinite(minAlt)) out.push({ az1, az2, minAlt });
    });
    return { enabled, sectors: out.slice(0, 5), topN };
    }
    function initFiltersUI() {
    document.getElementById('btn-filters')?.addEventListener('click', openFilters);
    document.querySelector('#filters-backdrop .close')?.addEventListener('click', closeFilters);
    document.getElementById('filters-cancel')?.addEventListener('click', closeFilters);
    document.getElementById('filters-reset')?.addEventListener('click', () => {
        saveSettings({ enabled:false, sectors:[], topN:16 });
        closeFilters();
        renderNow();
        applyHourFilters();
    });
    document.getElementById('filters-save')?.addEventListener('click', () => {
        const s = readSettingsFromForm();
        saveSettings(s);
        closeFilters();
        renderNow();
        applyHourFilters();
    });
    document.getElementById('filters-backdrop')?.addEventListener('click', (e) => {
        if (e.target.id === 'filters-backdrop') closeFilters();
    });
    document.addEventListener('keydown', (e) => { if (e.key === 'Escape') closeFilters(); });
    }

    /* ==============================
    Boot
    ============================== */
    document.addEventListener('DOMContentLoaded', () => {
    const input = document.getElementById('q');
    if (input) input.addEventListener('input', () => filterTables(input.value));

    initTabsBehavior();
    initHourLinks();
    hideDataColumns();
    enableRowModals();
    initFiltersUI();

    activateMainTab('panel-now');  // land on Now
    renderNow();

    const btn = document.getElementById('btn-refresh-now');
    if (btn) btn.addEventListener('click', () => renderNow());

    setInterval(renderNow, 30000);

    applyHourFilters(); // safe no-op unless Hour panel is visible
    });
    </script>
    """

    js = (
        js_template
        .replace("__PREVIEW_JSON__", preview_json)
        .replace("__NOW_JSON__", now_json)
        .replace("__HOURS_JSON__", hours_json)
    )

    # ---------------- HTML skeleton ----------------
    now_str = datetime.now(tz.gettz(tzname)).strftime("%Y-%m-%d %H:%M %Z")
    hour_links_html = " &middot; ".join(hour_links) if hour_links else "No hourly targets"

    content = f"""
      <section id="panel-now" class="tab-panel">
        <div class="small" style="margin:.2rem 0 .6rem;">
          Targets currently in view. Updates every two minutes.
          <button id="btn-refresh-now" class="tab" style="margin-left:.5rem;">Refresh</button>
        </div>
        <div id="now-table-wrap" class="table-wrap"><p class="small">Loading current targets…</p></div>
      </section>

      <section id="panel-master" class="tab-panel hidden">
        {master_html}
      </section>

      <section id="panel-hour" class="tab-panel hidden">
        <div id="hour-panel-wrap" class="table-wrap">
          <p class="small">Pick an hour from the chips above to load its table.</p>
        </div>
      </section>
    """

    html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Star Party Planner – {date_str}</title>
{css}
<div class="container">
  <h1>Star Party Planner</h1>

  <div class="meta-bar">
    <span class="pill">Date: {date_str}</span>
    <span class="pill">Observation window: {start}–{end} {tzname}</span>
  </div>

  <div class="navbar">
    <div class="left">
      <input id="q" type="search" placeholder="Search targets… (name, type, notes)" aria-label="Search">
      <div class="tabs">
        <div class="tab active" data-tab="panel-now">Now</div>
        <div class="tab" data-tab="panel-master">Master List</div>
      </div>
      <button id="btn-filters" class="tab" type="button" style="margin-left:.25rem;">Filters</button>
    </div>
    <div class="right">
      <span class="hours">{hour_links_html}</span>
    </div>
  </div>

  {content}

  <!-- Object preview modal -->
  <div id="modal-backdrop" class="modal-backdrop" role="dialog" aria-modal="true" aria-hidden="true">
    <div class="modal" role="document">
      <button class="close" aria-label="Close">Close ✕</button>
      <h3 id="modal-title">Object</h3>
      <div class="meta" id="modal-meta"></div>
      <div class="imgbox" id="modal-image"><p class="small">Loading preview…</p></div>
    </div>
  </div>

  <!-- Filters dialog (compact layout) -->
  <div id="filters-backdrop" class="modal-backdrop" role="dialog" aria-modal="true" aria-hidden="true">
    <div class="modal filters-dialog" role="document" style="max-width:720px;">
      <button class="close" aria-label="Close">Close ✕</button>
      <h3>Directional + Altitude Filters</h3>

      <!-- Master enable -->
      <label style="display:flex;align-items:center;gap:.5rem;margin:.5rem 0 .75rem;">
        <input type="checkbox" id="obs-enabled">
        <span>Enable filters</span>
      </label>

      <!-- Top N for Now -->
      <div class="topn-row">
        <label for="now-topn">Top N for Now:</label>
        <input id="now-topn" type="number" min="1" max="64" step="1" value="16">
        <span class="small">Applies after sector filtering.</span>
      </div>

      <!-- Sectors -->
      <h4 style="margin:.5rem 0 .35rem; color:#f66;">Sectors (up to 5)</h4>
      <div class="small" style="margin-bottom:.35rem;">Each row: Az start°, Az end°, Min altitude°</div>

      <div id="sector-rows">
        <div class="filters-row">
          <input class="az1" type="number" placeholder="Az start (0–360)" min="0" max="360">
          <input class="az2" type="number" placeholder="Az end (0–360)" min="0" max="360">
          <input class="minAlt" type="number" placeholder="Min alt (0–90)" min="0" max="90">
        </div>
        <div class="filters-row">
          <input class="az1" type="number" placeholder="Az start (0–360)" min="0" max="360">
          <input class="az2" type="number" placeholder="Az end (0–360)" min="0" max="360">
          <input class="minAlt" type="number" placeholder="Min alt (0–90)" min="0" max="90">
          </div>
        <div class="filters-row">
          <input class="az1" type="number" placeholder="Az start (0–360)" min="0" max="360">
          <input class="az2" type="number" placeholder="Az end (0–360)" min="0" max="360">
          <input class="minAlt" type="number" placeholder="Min alt (0–90)" min="0" max="90">
        </div>
        <div class="filters-row">
          <input class="az1" type="number" placeholder="Az start (0–360)" min="0" max="360">
          <input class="az2" type="number" placeholder="Az end (0–360)" min="0" max="360">
          <input class="minAlt" type="number" placeholder="Min alt (0–90)" min="0" max="90">
        </div>
        <div class="filters-row">
          <input class="az1" type="number" placeholder="Az start (0–360)" min="0" max="360">
          <input class="az2" type="number" placeholder="Az end (0–360)" min="0" max="360">
          <input class="minAlt" type="number" placeholder="Min alt (0–90)" min="0" max="90">
        </div>
      </div>

      <!-- Actions -->
      <div class="buttons" style="margin-top:.75rem;">
        <button type="button" id="filters-reset">Reset</button>
        <button type="button" id="filters-cancel">Cancel</button>
        <button type="button" id="filters-save">Save</button>
      </div>
    </div>
  </div>

  <p class="small">Location: {site_lat:.6f}, {site_lon:.6f}. Generated: {now_str}</p>
</div>
{js}
</html>
"""
    Path(output_path).write_text(html, encoding="utf-8")

    Path(output_path).write_text(html, encoding="utf-8")

# ---------------------------- Planning & tables ----------------------------

def plan_for_site(args):
    load, ts, eph = load_ephemeris(args.bsp)
    args = fill_time_defaults(args, ts, eph)

    earth = eph['earth']
    site = build_observer(load, ts, args.lat, args.lon, args.elev)
    hours = hours_list(args.date, args.start, args.end, args.tz)

    start_dt = hours[0]
    end_dt = hours[-1]
    if end_dt <= start_dt:
        end_dt += timedelta(days=1)

    minute_times_local = []
    t_cursor = start_dt
    while t_cursor <= end_dt:
        minute_times_local.append(t_cursor)
        t_cursor += timedelta(minutes=2)
    ts_minute = ts.from_datetimes(minute_times_local)

    # Moon RA/Dec over grid for separation tests
    moon_ras, moon_decs = [], []
    for t in ts_minute:
        ra_m, dec_m = moon_ra_dec(eph, t, site, earth)
        moon_ras.append(ra_m); moon_decs.append(dec_m)
    moon_ras = np.array(moon_ras); moon_decs = np.array(moon_decs)

    # Night-wide illumination fraction (mid-window)
    mid_idx = len(ts_minute) // 2
    phase_mid_deg = almanac.moon_phase(eph, ts_minute[mid_idx]).degrees
    illum_frac_night = lunar_illum_fraction_from_phase_deg(phase_mid_deg)

    planet_names = ["Mercury","Venus","Mars","Jupiter","Saturn","Uranus","Neptune"]
    targets: List[Dict] = []

    dso_list = read_catalog(args.catalog)

    def eval_target(name, ra_deg, dec_deg, typ, mag=None, notes=""):
        alts, azs = [], []
        for t in ts_minute:
            alt, az = to_altaz(t, site, earth, ra_deg, dec_deg)
            alts.append(alt); azs.append(az)
        alts = np.array(alts); azs = np.array(azs)

        seps = angular_sep(ra_deg, dec_deg, moon_ras, moon_decs)
        visible_mask = (alts >= args.min_alt) & (seps >= args.moon_sep_min)
        if not np.any(visible_mask):
            return

        idx_best = best_time_in_window(ts_minute, alts)
        best_local = minute_times_local[idx_best]
        best_alt = float(alts[idx_best]); best_az = float(azs[idx_best])
        best_dir = cardinal_from_az(best_az)

        # Moon altitude at the object's best time
        alt_moon_best, _, _ = moon_altaz_phase(eph, ts.from_datetime(best_local), site, earth)

        hourly_rows = []
        for h in hours:
            t_h = ts.from_datetime(h)
            alt_h, az_h = to_altaz(t_h, site, earth, ra_deg, dec_deg)
            if alt_h >= args.min_alt:
                dir_h = cardinal_from_az(az_h)
                prio = interest_score(name, typ, best_alt, alt_now=alt_h)
                # Moon altitude at this hour
                alt_moon_h, _, _ = moon_altaz_phase(eph, t_h, site, earth)
                prio -= moonlight_penalty_points(typ, illum_frac_night, args.moonlight_penalty_max,
                                                 moon_alt_deg=alt_moon_h)
                hourly_rows.append((h, float(alt_h), float(az_h), dir_h, max(0.0, prio)))

        if not hourly_rows:
            return

        master_prio = interest_score(name, typ, best_alt)
        master_prio -= moonlight_penalty_points(typ, illum_frac_night, args.moonlight_penalty_max,
                                                moon_alt_deg=alt_moon_best)

        targets.append({
            "name": name, "type": typ, "mag": mag, "notes": notes,
            "best_time": best_local, "best_alt": round(best_alt,1),
            "best_dir": best_dir, "best_az": round(best_az,1),
            "ra_deg": ra_deg, "dec_deg": dec_deg,
            "priority": max(0.0, master_prio),
            "hourly": hourly_rows
        })

    # DSOs
    for d in dso_list:
        if d.mag is not None and d.mag > args.max_mag:
            continue
        eval_target(d.name, d.ra_deg, d.dec_deg, d.type, d.mag, d.notes)

    # Planets
    for name in planet_names:
        alts, azs = [], []
        for t in ts_minute:
            alt, az, _ = planet_altaz(eph, t, site, earth, name)
            alts.append(alt); azs.append(az)
        alts = np.array(alts); azs = np.array(azs)
        if not np.any(alts >= args.min_alt_planets):
            continue
        idx_best = best_time_in_window(ts_minute, alts)
        best_local = minute_times_local[idx_best]
        best_alt = float(alts[idx_best]); best_az = float(azs[idx_best])
        best_dir = cardinal_from_az(best_az)

        hourly_rows = []
        for h in hours:
            t_h = ts.from_datetime(h)
            alt_h, az_h, _ = planet_altaz(eph, t_h, site, earth, name)
            if alt_h >= args.min_alt_planets:
                dir_h = cardinal_from_az(az_h)
                prio = interest_score(name, "Planet", best_alt, alt_now=alt_h)
                hourly_rows.append((h, float(alt_h), float(az_h), dir_h, prio))

        if hourly_rows:
            targets.append({
                "name": name, "type": "Planet", "mag": None, "notes": "",
                "best_time": best_local, "best_alt": round(best_alt,1),
                "best_dir": best_dir, "best_az": round(best_az,1),
                "ra_deg": np.nan, "dec_deg": np.nan,
                "priority": interest_score(name, "Planet", best_alt),
                "hourly": hourly_rows
            })

    # Moon
    alts, azs, phases = [], [], []
    for t in ts_minute:
        alt, az, phase = moon_altaz_phase(eph, t, site, earth)
        alts.append(alt); azs.append(az); phases.append(phase)
    alts = np.array(alts); azs = np.array(azs)
    if np.any(alts >= args.min_alt_moon):
        idx_best = best_time_in_window(ts_minute, alts)
        best_local = minute_times_local[idx_best]
        best_alt = float(alts[idx_best])
        hourly_rows = []
        for h in hours:
            t_h = ts.from_datetime(h)
            alt_h, az_h, phase_h = moon_altaz_phase(eph, t_h, site, earth)
            if alt_h >= args.min_alt_moon:
                dir_h = cardinal_from_az(az_h)
                prio = interest_score("Moon", "Moon", best_alt, alt_now=alt_h)
                hourly_rows.append((h, float(alt_h), float(az_h), dir_h, prio, float(phase_h)))
        if hourly_rows:
            targets.append({
                "name": "Moon", "type": "Moon", "mag": None, "notes": "",
                "best_time": best_local, "best_alt": round(best_alt,1),
                "best_dir": cardinal_from_az(float(azs[idx_best])),
                "best_az": round(float(azs[idx_best]),1),
                "ra_deg": np.nan, "dec_deg": np.nan,
                "priority": interest_score("Moon", "Moon", best_alt),
                "hourly": hourly_rows
            })

    # Master table
    master_rows = []
    for t in targets:
        master_rows.append({
            "Name": t["name"],
            "Type": t["type"],
            "Mag": t["mag"] if t["mag"] is not None else "",
            "Best Time": t["best_time"].strftime("%H:%M"),
            "Dir": t["best_dir"],
            "Alt (°)": t["best_alt"],
            "Az (°)": t["best_az"],
            "Notes": t["notes"],
            "RA (deg)": round(t["ra_deg"], 3) if not pd.isna(t["ra_deg"]) else "",
            "Dec (deg)": round(t["dec_deg"], 3) if not pd.isna(t["dec_deg"]) else "",
            "_Priority": t["priority"],
            "_BestDT": t["best_time"],
        })
    master_df = pd.DataFrame(master_rows)
    if not master_df.empty:
        master_df.sort_values(by=["_Priority","_BestDT"], ascending=[False, True], inplace=True)
        master_df.drop(columns=["_Priority","_BestDT"], inplace=True)

    # Hourly tables
    hour_tables = []
    for t in targets:
        for row in t["hourly"]:
            h, alt_h, az_h, dir_h, prio, *rest = row
            out = {
                "Hour": h.strftime("%Y-%m-%d %H:%M"),
                "Name": t["name"],
                "Type": t["type"],
                "Mag": t["mag"] if t["mag"] is not None else "",
                "Dir": dir_h,
                "Alt (°)": round(alt_h, 1),
                "Az (°)": round(az_h, 1),
                "RA (deg)": round(t["ra_deg"], 3) if not pd.isna(t["ra_deg"]) else "",
                "Dec (deg)": round(t["dec_deg"], 3) if not pd.isna(t["dec_deg"]) else "",
                "Notes": t.get("notes", ""),
                "_Priority": prio,
            }
            if t["name"] == "Moon" and rest:
                out["Moon Phase (°)"] = round(rest[-1], 1)
            hour_tables.append(out)

    hourly_df = pd.DataFrame(hour_tables)
    if not hourly_df.empty:
        hourly_df["Hour_dt"] = pd.to_datetime(hourly_df["Hour"])
        hourly_df.sort_values(["Hour_dt","_Priority"], ascending=[True,False], inplace=True)
        hourly_df["rank"] = hourly_df.groupby("Hour")["_Priority"].rank(method="first", ascending=False)
        hourly_df = hourly_df[hourly_df["rank"] <= args.top_n_per_hour]
        hourly_df.drop(columns=["rank","Hour_dt","_Priority"], inplace=True, errors="ignore")

    # -------- NOW view (2-minute granularity, client-side) --------
    # Build an array of dicts: one entry per 2-minute slot with top N objects
    # Uses same scoring rules (incl. moonlight penalty) as hourly
    now_data = []
    top_n_now = args.top_n_per_hour  # reuse the same cap

    # Pre-create light-weight target descriptors we can iterate quickly
    # We will re-compute alt/az at each 2-min tick (fast enough for typical lists)
    simple_targets = []
    for t in targets:
        simple_targets.append({
            "name": t["name"],
            "type": t["type"],
            "mag": t["mag"],
            "ra_deg": t["ra_deg"],
            "dec_deg": t["dec_deg"],
            "min_alt": (
                args.min_alt_planets if t["type"].lower() == "planet"
                else (args.min_alt_moon if t["type"].lower() == "moon" else args.min_alt)
            ),
        })

    # We will reuse illum_frac_night (night-wide), but scale by the *instant* Moon altitude
    def _altaz_for_target_at(ts_obj, targ):
        typ = (targ["type"] or "").lower()
        if typ == "planet":
            alt, az, _ = planet_altaz(eph, ts_obj, site, earth, targ["name"])
        elif typ == "moon":
            alt, az, _ = moon_altaz_phase(eph, ts_obj, site, earth)
        else:
            alt, az = to_altaz(ts_obj, site, earth, targ["ra_deg"], targ["dec_deg"])
        return float(alt), float(az)

    for dt_local in minute_times_local:
        ts_tick = ts.from_datetime(dt_local)

        # Moon altitude for penalty scaling at this instant
        moon_alt_tick, _, moon_phase_tick = moon_altaz_phase(eph, ts_tick, site, earth)
        # We keep illumination fraction fixed for the night; if you prefer per-tick,
        # replace illum_frac_night with lunar_illum_fraction_from_phase_deg(moon_phase_tick).
        rows = []
        for targ in simple_targets:
            alt, az = _altaz_for_target_at(ts_tick, targ)
            if alt < targ["min_alt"]:
                continue

            # Dir, prio with penalties as needed
            dir_txt = cardinal_from_az(az)

            # Base interest uses 'best_alt' vs 'alt_now'; at "Now" we only have current alt
            # We can approximate by using alt as both best_alt and alt_now (works well in practice)
            base_prio = interest_score(targ["name"], targ["type"], best_alt=alt, alt_now=alt)

            # Moonlight penalty for diffuse DSOs only (planets/Moon unaffected)
            prio = base_prio
            if targ["type"].lower() not in {"planet", "moon"}:
                prio -= moonlight_penalty_points(
                    targ["type"], illum_frac_night, args.moonlight_penalty_max,
                    moon_alt_deg=moon_alt_tick
                )
            if prio <= 0:
                continue

            rows.append({
                "Name": targ["name"],
                "Type": targ["type"],
                "Mag": targ["mag"] if targ["mag"] is not None else "",
                "Alt (°)": round(alt, 1),
                "Az (°)": round(az, 1),
                "Dir": dir_txt,
                "RA (deg)": _round3(targ["ra_deg"]),
                "Dec (deg)": _round3(targ["dec_deg"]),
                "_Priority": prio
            })

        if rows:
            rows.sort(key=lambda r: (-r["_Priority"], r["Name"]))
            # Preserve score for client-side sorting/filtering; do NOT slice here
            for r in rows:
                r["_Score"] = float(r.pop("_Priority", 0.0))

            now_data.append({
                "time": dt_local.strftime("%Y-%m-%d %H:%M"),
                "rows": rows
            })

    # Save CSVs
    master_csv = f"{args.out_prefix}_master.csv"
    hourly_csv = f"{args.out_prefix}_hourly.csv"
    master_df.to_csv(master_csv, index=False); hourly_df.to_csv(hourly_csv, index=False)

    # ---- Build preview cache & cleanup if requested ----
    preview_map: Dict[str, Dict[str,str]] = {}
    cache_dir = Path(args.preview_cache_dir)
    if not args.no_previews:
        preview_map = build_preview_cache(master_df, cache_dir, args.refresh_previews,
                                          args.preview_fov_deg, args.preview_px)
        if args.clean_preview_cache:
            cleanup_preview_cache(cache_dir, preview_map)

    # HTML output
    if args.html:
        write_html(args.html, args.lat, args.lon, args.tz,
                   args.date, args.start, args.end,
                   master_df, hourly_df, args.html_ui, preview_map,
                   now_data=now_data)

    # Console
    print("\n=== MASTER LIST (by interestingness) ===")
    if not master_df.empty:
        print(master_df[["Name","Type","Alt (°)","Az (°)","Dir","Best Time","Mag","Notes"]].to_string(index=False))
    else:
        print("No targets passed the filters. Try adjusting filters or time window.")

    print("\n=== HOURLY 'POINT YOUR SCOPE' TABLES (by interest) ===")
    if not hourly_df.empty:
        for hour, sub in hourly_df.groupby("Hour"):
            print(f"\n# {hour}")
            cols = [c for c in ["Name","Type","Mag","Alt (°)","Az (°)","Dir"] if c in sub.columns]
            print(sub[cols].to_string(index=False))
    else:
        print("No hourly targets above your altitude threshold.")

if __name__ == "__main__":
    args = parse_args()
    plan_for_site(args)
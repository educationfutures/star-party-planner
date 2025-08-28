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

Usage (examples)
  python starparty_planner.py --lat 44.98 --lon -93.26 --date 2025-08-28 \
      --start 20:00 --end 01:00 --tz America/Chicago --catalog objects_sample.csv

Dependencies
  pip install skyfield numpy pandas pytz python-dateutil
  # Skyfield 1.45+ recommended

Notes
  - Uses a LOCAL ephemeris only (default: ./skyfield_data/de440s.bsp).
    Download once from NAIF and place it there:
      mkdir -p ./skyfield_data
      curl -LO https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
      mv de440s.bsp ./skyfield_data/
  - Deep-sky objects come from a CSV you provide (see objects_sample.csv).
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import List, Dict, Tuple, Optional
from pathlib import Path

import numpy as np
import pandas as pd
from dateutil import tz

from skyfield.api import Loader, wgs84, load_file, Star
from skyfield.almanac import dark_twilight_day
from skyfield import almanac

# ---------------------------- Config helpers ----------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Hourly observing list generator")
    p.add_argument("--lat", type=float, required=True, help="Latitude in decimal degrees (N+)")
    p.add_argument("--lon", type=float, required=True, help="Longitude in decimal degrees (E+)")
    p.add_argument("--elev", type=float, default=0.0, help="Elevation meters (optional)")
    p.add_argument("--date", type=str, required=True, help="Local date YYYY-MM-DD")
    p.add_argument("--start", type=str, default="20:00", help="Local start time HH:MM (default 20:00)")
    p.add_argument("--end", type=str, default="01:00", help="Local end time HH:MM (default 01:00)")
    p.add_argument("--tz", type=str, default="UTC", help="IANA timezone, e.g., America/Chicago")
    p.add_argument("--catalog", type=str, default="objects_sample.csv", help="CSV of DSOs")
    p.add_argument("--min_alt", type=float, default=20.0, help="Minimum altitude for DSOs (deg)")
    p.add_argument("--max_mag", type=float, default=9.0, help="Max magnitude (fainter=larger) for DSOs")
    p.add_argument("--moon_sep_min", type=float, default=15.0, help="Min separation from Moon (deg)")
    p.add_argument("--hour_step", type=int, default=1, help="Hour step for hourly tables")
    p.add_argument("--top_n_per_hour", type=int, default=10, help="Max targets per hour table")
    p.add_argument("--out_prefix", type=str, default="starparty", help="Output prefix")
    p.add_argument("--html", type=str, default="", help="Optional HTML output filepath (night-vision themed)")
    p.add_argument("--bsp", type=str, default="./skyfield_data/de440s.bsp",
                   help="Path to local planetary ephemeris BSP (e.g., de440s.bsp)")
    p.add_argument("--html_ui", type=str, default="accordion", choices=["accordion", "tabs"],
                   help="HTML layout: accordion (per-hour collapsible) or tabs (Master/By Hour)")
    # Per-type altitude thresholds
    p.add_argument("--min_alt_planets", type=float, default=10.0, help="Minimum altitude for planets (deg)")
    p.add_argument("--min_alt_moon", type=float, default=5.0, help="Minimum altitude for the Moon (deg)")
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
    phase = almanac.moon_phase(eph, ts).degrees  # 0=new, 180=full
    return alt.degrees, az.degrees, phase

def best_time_in_window(ts_arr, alts) -> int:
    return int(np.nanargmax(alts))

def build_observer(load, ts, lat, lon, elev):
    return wgs84.latlon(latitude_degrees=lat, longitude_degrees=lon, elevation_m=elev)

def compute_dark_window(eph, ts, tzname: str, local_date: str, lat: float, lon: float) -> Tuple[datetime, datetime]:
    """Astronomical night (Sun < -18 deg). Not used directly, but kept for future."""
    T = tz.gettz(tzname)
    day = datetime.strptime(local_date, "%Y-%m-%d").replace(tzinfo=T)
    t0 = day.replace(hour=12, minute=0)
    t1 = t0 + timedelta(days=2)
    ts0 = ts.from_datetime(t0)
    ts1 = ts.from_datetime(t1)
    f = dark_twilight_day(eph, wgs84.latlon(lat, lon), horizon=-18.0)
    times, events = almanac.find_discrete(ts0, ts1, f)
    night_starts = []
    for i in range(len(events)-1):
        if events[i] == 3 and events[i+1] == 4:
            start = times[i].utc_datetime().astimezone(T)
            end = times[i+1].utc_datetime().astimezone(T)
            if start.date() <= day.date() <= end.date() or start.date()==day.date():
                night_starts.append((start, end))
    if night_starts:
        return night_starts[0]
    else:
        return day.replace(hour=21, minute=0), (day+timedelta(days=1)).replace(hour=4, minute=0)


# ---------------------------- Interest scoring ----------------------------

INTEREST_BASE = {
    "Saturn": 100, "Jupiter": 95, "Moon": 90, "Mars": 80, "Venus": 70,
    "Uranus": 55, "Neptune": 50,
}
TYPE_BONUS = {
    "Globular cluster": 45, "Open cluster": 35, "Planetary nebula": 40,
    "Emission Nebula": 40, "Reflection Nebula": 35, "Nebula with cluster": 40,
    "H II region nebula with cluster": 40, "Spiral galaxy": 35,
    "Elliptical galaxy": 30, "Starburst galaxy": 34, "Galaxy": 32,
    "Supernova Remnant": 42, "Milky Way star cloud": 38,
    "Asterism": 28, "Optical Double": 20,
}
# Extra bonus for sure-fire crowd-pleasers
CROWD_BONUS = {"Saturn": 30, "Jupiter": 25, "Moon": 20, "Mars": 12, "Venus": 10}

def interest_score(name: str, typ: str, best_alt: float, alt_now: Optional[float] = None) -> float:
    base = INTEREST_BASE.get(name, 0)
    if base == 0:
        base = TYPE_BONUS.get(typ, 25)
    base += CROWD_BONUS.get(name, 0)
    alt_term = 0.10 * best_alt + 0.15 * (alt_now if alt_now is not None else best_alt)
    return min(200.0, base + alt_term)


# ---------------------------- HTML rendering ----------------------------

def df_to_html_table(df: pd.DataFrame, id_attr: str = "") -> str:
    if df.empty:
        return "<p>No data.</p>"
    # Wrap in a scroller for mobile
    return f'<div class="table-wrap">{df.to_html(index=False, escape=True, border=0, table_id=id_attr)}</div>'

def _hour_anchor_label(hour_str: str) -> str:
    try:
        t = pd.to_datetime(hour_str)
        return t.strftime("%H:%M")
    except Exception:
        return hour_str

def write_html(output_path: str, site_lat: float, site_lon: float, tzname: str, date_str: str,
               start: str, end: str, master_df: pd.DataFrame, hourly_df: pd.DataFrame,
               ui_mode: str = "accordion"):
    # Group hourly sections
    hourly_sections = []
    hour_links = []
    if not hourly_df.empty:
        for hour, sub in hourly_df.groupby("Hour"):
            anchor = f"hour-{hour.replace(' ','_').replace(':','')}"
            hour_label = _hour_anchor_label(hour)
            hour_links.append(f'<a href="#{anchor}">{hour_label}</a>')
            tbl = df_to_html_table(sub[['Name','Type','Alt (°)','Dir','Az (°)']], id_attr=f"tbl-{anchor}")
            if ui_mode == "accordion":
                hourly_sections.append(f"""
                <details id="{anchor}" class="acc">
                  <summary><span class="acc-time">{hour_label}</span>
                    <span class="acc-count">{len(sub)} targets</span>
                  </summary>
                  {tbl}
                </details>
                """)
            else:
                hourly_sections.append(f"""
                <section id="{anchor}" class="hour-tab-panel hidden">
                  <h3>{hour}</h3>
                  {tbl}
                </section>
                """)
    hourly_html = "\n".join(hourly_sections) if hourly_sections else "<p>No hourly targets above your altitude threshold.</p>"

    master_html = df_to_html_table(master_df, id_attr="tbl-master") if not master_df.empty else "<p>No targets passed the filters. Try adjusting filters.</p>"

    # Night-vision CSS + mobile fixes (no white anywhere; disable sticky header on phones)
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

      /* Tables */
      .table-wrap { overflow-x: auto; -webkit-overflow-scrolling: touch; }
      table { width: 100%; border-collapse: collapse; margin: 0.5rem 0 1rem; }
      th, td { border: 1px solid #700; padding: 0.45rem 0.5rem; }
      th { background: #100; position: sticky; top: 48px; z-index: 1; }
      tr:nth-child(even) { background: #070707; }
      tr:hover { background: #111; }
      th.sortable { cursor: pointer; }
      th.sortable:after { content: " ⇅"; color:#f66; font-weight: normal; }

      /* Toolbar */
      .toolbar { position: sticky; top: 0; background: rgba(0,0,0,0.98); border-bottom: 1px solid #400; padding: 0.5rem; z-index: 5; display:flex; gap:0.6rem; align-items:center; flex-wrap:wrap;}
      .pill { border:1px solid #700; padding:0.25rem 0.5rem; border-radius:999px; color:#f66; }
      input[type="search"] { background: #160000; border: 1px solid #700; color: #f55; padding: 0.4rem 0.6rem; border-radius: 6px; min-width: 220px; caret-color: #f55; outline: none; }
      input[type="search"]::placeholder { color:#f66; }
      input[type="search"]:focus { box-shadow: 0 0 0 2px #500 inset; border-color:#900; }
      * { -webkit-tap-highlight-color: rgba(255, 0, 0, 0.2); }

      /* Prevent white autofill on iOS/Safari/Chrome */
      input:-webkit-autofill,
      input:-webkit-autofill:hover,
      input:-webkit-autofill:focus {
        -webkit-text-fill-color: #f55;
        transition: background-color 5000s ease-in-out 0s;
        box-shadow: 0 0 0px 1000px #160000 inset;
        border: 1px solid #700;
      }

      /* Hour links */
      .hours { display:flex; gap:0.35rem; flex-wrap:wrap; max-height: 2.3rem; overflow:auto; }
      .hours a { display:inline-block; padding: 0.25rem 0.55rem; border:1px solid #500; border-radius:6px; text-decoration:none; }
      .hours a:focus { outline:none; box-shadow: 0 0 0 2px #500 inset; }

      /* Accordions */
      .acc { border: 1px solid #400; border-radius: 8px; margin: 0.5rem 0; background:#050505; }
      .acc summary { cursor: pointer; padding: 0.5rem 0.7rem; list-style:none; display:flex; justify-content:space-between; align-items:center; }
      .acc summary::-webkit-details-marker { display:none; }
      .acc summary:hover { background:#0a0a0a; }
      .acc-time { font-weight:600; color:#f55; }
      .acc-count { font-size:0.9rem; color:#f77; }

      /* Tabs */
      .tabs { display:flex; gap:0.4rem; margin: 0.6rem 0 0.8rem; }
      .tab { padding:0.35rem 0.7rem; border:1px solid #500; border-radius:8px; cursor:pointer; user-select:none; color:#f66; background:#0a0000; }
      .tab.active { background:#180000; border-color:#700; color:#f66; } /* no white */

      .hidden { display:none; }

      /* Phone tweaks */
      @media (max-width: 640px) {
        body { font-size: 15px; }
        th, td { padding: 0.35rem 0.45rem; }
        .pill { font-size: 0.85rem; }
        .tabs { gap:0.3rem; }
        .tab { padding:0.3rem 0.55rem; }
        .toolbar { gap:0.5rem; }
        /* Disable sticky headers on phones to avoid overlay issues */
        th { position: static !important; top: auto !important; z-index: auto !important; }
      }
    </style>
    """

    # JS: search filter + tab/hour navigation + sortable master table
    js = """
    <script>
    function filterTables(q) {
      const needle = q.trim().toLowerCase();
      const tables = document.querySelectorAll('table');
      tables.forEach(tbl => {
        const rows = tbl.tBodies[0]?.rows || [];
        for (let r of rows) {
          const txt = r.innerText.toLowerCase();
          r.style.display = (!needle || txt.includes(needle)) ? '' : 'none';
        }
      });
    }

    function activateMainTab(targetId) {
      document.querySelectorAll('.tab').forEach(el => el.classList.remove('active'));
      document.querySelectorAll('.tab-panel').forEach(el => el.classList.add('hidden'));
      const tabBtn = document.querySelector('.tab[data-tab="' + targetId + '"]');
      if (tabBtn) tabBtn.classList.add('active');
      const panel = document.getElementById(targetId);
      if (panel) panel.classList.remove('hidden');
    }

    function showHourPanel(hourId) {
      activateMainTab('panel-hourly');
      const all = document.querySelectorAll('.hour-tab-panel');
      all.forEach(el => el.classList.add('hidden'));
      const target = document.getElementById(hourId);
      if (target) {
        target.classList.remove('hidden');
        const y = target.getBoundingClientRect().top + window.scrollY - 70;
        window.scrollTo({ top: y, behavior: 'instant' });
      }
    }

    function openAccordion(hourId) {
      const det = document.getElementById(hourId);
      if (det && det.tagName.toLowerCase() === 'details') {
        det.open = true;
        const y = det.getBoundingClientRect().top + window.scrollY - 70;
        window.scrollTo({ top: y, behavior: 'instant' });
      }
    }

    function handleHourNavClick(e) {
      const href = e.currentTarget.getAttribute('href') || '';
      if (!href.startsWith('#')) return;
      e.preventDefault();
      const id = href.slice(1);
      if (document.getElementById('panel-hourly')) {
        showHourPanel(id);
      } else {
        openAccordion(id);
      }
      history.replaceState(null, '', href);
    }

    function initTabsBehavior() {
      const tabs = document.querySelectorAll('[data-tab]');
      tabs.forEach(tab => {
        tab.addEventListener('click', () => {
          const target = tab.getAttribute('data-tab');
          activateMainTab(target);
          if (target === 'panel-hourly') {
            const first = document.querySelector('.hour-tab-panel');
            if (first) {
              document.querySelectorAll('.hour-tab-panel').forEach(el => el.classList.add('hidden'));
              first.classList.remove('hidden');
            }
          }
        });
      });
    }

    function initHourLinks() {
      document.querySelectorAll('.toolbar .hours a').forEach(a => {
        a.addEventListener('click', handleHourNavClick);
      });
    }

    function handleHashOnLoad() {
      const hash = window.location.hash;
      if (!hash) return;
      const id = hash.slice(1);
      if (document.getElementById('panel-hourly')) {
        showHourPanel(id);
      } else {
        openAccordion(id);
      }
    }

    // ---- Sortable Master Table ----
    function parseHHMM(s) {
      // expects "HH:MM"
      const m = /^\\s*(\\d{1,2}):(\\d{2})\\s*$/.exec(s || "");
      if (!m) return NaN;
      return parseInt(m[1], 10) * 60 + parseInt(m[2], 10);
    }

    function sortTableByColumn(table, colIndex, numeric=false, time=false, asc=true) {
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.querySelectorAll('tr'));
      const dir = asc ? 1 : -1;

      rows.sort((a, b) => {
        let A = a.children[colIndex]?.innerText || "";
        let B = b.children[colIndex]?.innerText || "";
        if (time) {
          A = parseHHMM(A); B = parseHHMM(B);
          if (isNaN(A)) A = -Infinity; if (isNaN(B)) B = -Infinity;
        } else if (numeric) {
          A = parseFloat(A); B = parseFloat(B);
          if (isNaN(A)) A = -Infinity; if (isNaN(B)) B = -Infinity;
        } else {
          A = A.toLowerCase(); B = B.toLowerCase();
        }
        if (A < B) return -1 * dir;
        if (A > B) return  1 * dir;
        return 0;
      });

      rows.forEach(r => tbody.appendChild(r));
    }

    function makeMasterTableSortable() {
      const tbl = document.getElementById('tbl-master');
      if (!tbl) return;

      // Find the columns we want clickable
      const headers = Array.from(tbl.tHead?.rows[0]?.cells || []);
      const colMap = {}; // name -> (idx, opts)
      headers.forEach((th, i) => {
        const label = (th.innerText || "").trim().toLowerCase();
        if (["name","type"].includes(label)) {
          colMap[label] = { idx: i, numeric: false, time: false };
        } else if (label.startsWith("mag")) {
          colMap["mag"] = { idx: i, numeric: true, time: false };
        } else if (label.startsWith("best local time")) {
          colMap["best local time"] = { idx: i, numeric: false, time: true };
        }
      });

      // Mark sortable headers and attach listeners
      Object.entries(colMap).forEach(([key, cfg]) => {
        const th = headers[cfg.idx];
        th.classList.add('sortable');
        let asc = true;
        th.addEventListener('click', () => {
          sortTableByColumn(tbl, cfg.idx, cfg.numeric, cfg.time, asc);
          asc = !asc;
        });
      });
    }

    document.addEventListener('DOMContentLoaded', () => {
      const input = document.getElementById('q');
      if (input) input.addEventListener('input', () => filterTables(input.value));
      initTabsBehavior();
      initHourLinks();
      makeMasterTableSortable();
      handleHashOnLoad();
    });
    </script>
    """

    now = datetime.now(tz.gettz(tzname)).strftime("%Y-%m-%d %H:%M %Z")
    hour_links_html = " &middot; ".join(hour_links) if hour_links else "No hourly targets"

    if ui_mode == "tabs":
        content = f"""
          <div class="tabs">
            <div class="tab active" data-tab="panel-master">Master List</div>
            <div class="tab" data-tab="panel-hourly">By Hour</div>
          </div>
          <section id="panel-master" class="tab-panel">
            {master_html}
          </section>
          <section id="panel-hourly" class="tab-panel hidden">
            <p class="small">Top targets each hour, prioritized by interest score. Directions are compass points.</p>
            {hourly_html}
          </section>
        """
    else:
        content = f"""
          <h2>Master List (by interest)</h2>
          {master_html}
          <hr class="hr">
          <h2>Hourly “Point Your Scope Now”</h2>
          <p class="small">Top targets each hour, prioritized by interest score. Directions are compass points.</p>
          {hourly_html}
        """

    html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Star Party Planner – {date_str} ({tzname})</title>
{css}
<div class="container">
  <h1>Star Party Planner</h1>

  <div class="toolbar">
    <span class="pill">Location: {site_lat:.6f}, {site_lon:.6f}</span>
    <span class="pill">Date: {date_str}</span>
    <span class="pill">Window: {start}–{end} {tzname}</span>
    <span class="pill">Generated: {now}</span>
    <input id="q" type="search" placeholder="Search targets… (name, type, notes)" aria-label="Search">
    <span class="hours">{hour_links_html}</span>
  </div>

  {content}

  <p class="small warn">Night-vision mode: red on black. Keep device brightness low at the scope.</p>
</div>
{js}
</html>
"""
    Path(output_path).write_text(html, encoding="utf-8")


# ---------------------------- Planning & tables ----------------------------

def plan_for_site(args):
    load, ts, eph = load_ephemeris(args.bsp)
    earth = eph['earth']
    site = build_observer(load, ts, args.lat, args.lon, args.elev)

    hours = hours_list(args.date, args.start, args.end, args.tz)
    T = tz.gettz(args.tz)

    # Minute sampling grid to find "best time"
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

    # Moon RA/Dec over grid
    moon_ras = []
    moon_decs = []
    for t in ts_minute:
        ra_m, dec_m = moon_ra_dec(eph, t, site, earth)
        moon_ras.append(ra_m)
        moon_decs.append(dec_m)
    moon_ras = np.array(moon_ras)
    moon_decs = np.array(moon_decs)

    planet_names = ["Mercury","Venus","Mars","Jupiter","Saturn","Uranus","Neptune"]
    targets: List[Dict] = []

    dso_list = read_catalog(args.catalog)

    # Helper to evaluate a fixed RA/Dec over the minute grid
    def eval_target(name, ra_deg, dec_deg, typ, mag=None, notes=""):
        alts = []
        azs = []
        for t in ts_minute:
            alt, az = to_altaz(t, site, earth, ra_deg, dec_deg)
            alts.append(alt); azs.append(az)
        alts = np.array(alts); azs = np.array(azs)

        # Moon separation (deg)
        seps = angular_sep(ra_deg, dec_deg, moon_ras, moon_decs)

        # Filters (DSOs use --min_alt and Moon separation)
        visible_mask = (alts >= args.min_alt) & (seps >= args.moon_sep_min)
        if not np.any(visible_mask):
            return

        idx_best = best_time_in_window(ts_minute, alts)
        best_local = minute_times_local[idx_best]
        best_alt = float(alts[idx_best])
        best_az = float(azs[idx_best])
        best_dir = cardinal_from_az(best_az)

        # Hourly snapshots
        hourly_rows = []
        for h in hours:
            t_h = ts.from_datetime(h)
            alt_h, az_h = to_altaz(t_h, site, earth, ra_deg, dec_deg)
            if alt_h >= args.min_alt:
                dir_h = cardinal_from_az(az_h)
                prio = interest_score(name, typ, best_alt, alt_now=alt_h)
                hourly_rows.append((h, float(alt_h), float(az_h), dir_h, prio))

        if not hourly_rows:
            return

        targets.append({
            "name": name,
            "type": typ,
            "mag": mag,
            "notes": notes,
            "best_time": best_local,
            "best_alt": round(best_alt,1),
            "best_dir": best_dir,
            "best_az": round(best_az,1),
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
            "priority": interest_score(name, typ, best_alt),
            "hourly": hourly_rows
        })

    # Evaluate DSOs
    for d in dso_list:
        if d.mag is not None and d.mag > args.max_mag:
            continue
        eval_target(d.name, d.ra_deg, d.dec_deg, d.type, d.mag, d.notes)

    # Evaluate planets (use planets-specific min altitude)
    for name in planet_names:
        alts = []
        azs = []
        for t in ts_minute:
            alt, az, _ = planet_altaz(eph, t, site, earth, name)
            alts.append(alt); azs.append(az)
        alts = np.array(alts); azs = np.array(azs)

        visible_mask = (alts >= args.min_alt_planets)
        if not np.any(visible_mask):
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
                "name": name,
                "type": "Planet",
                "mag": None,
                "notes": "",
                "best_time": best_local,
                "best_alt": round(best_alt,1),
                "best_dir": best_dir,
                "best_az": round(best_az,1),
                "ra_deg": np.nan,
                "dec_deg": np.nan,
                "priority": interest_score(name, "Planet", best_alt),
                "hourly": hourly_rows
            })

    # Add Moon (context & crowd-pleaser) using moon-specific min altitude
    alts = []; azs = []; phases = []
    for t in ts_minute:
        alt, az, phase = moon_altaz_phase(eph, t, site, earth)
        alts.append(alt); azs.append(az); phases.append(phase)
    alts = np.array(alts); azs = np.array(azs)
    if np.any(alts >= args.min_alt_moon):
        idx_best = best_time_in_window(ts_minute, alts)
        best_alt = float(alts[idx_best])

        hourly_rows = []
        for h in hours:
            t_h = ts.from_datetime(h)
            alt_h, az_h, _ = moon_altaz_phase(eph, t_h, site, earth)
            if alt_h >= args.min_alt_moon:
                dir_h = cardinal_from_az(az_h)
                prio = interest_score("Moon", "Moon", best_alt, alt_now=alt_h)
                hourly_rows.append((h, float(alt_h), float(az_h), dir_h, prio))

        if hourly_rows:
            targets.append({
                "name": "Moon",
                "type": "Moon",
                "mag": None,
                "notes": "",
                "best_time": minute_times_local[idx_best],
                "best_alt": round(best_alt,1),
                "best_dir": cardinal_from_az(float(azs[idx_best])),
                "best_az": round(float(azs[idx_best]),1),
                "ra_deg": np.nan,
                "dec_deg": np.nan,
                "priority": interest_score("Moon", "Moon", best_alt),
                "hourly": hourly_rows
            })

    # Master table (sort by priority desc, then best-time)
    master_rows = []
    for t in targets:
        master_rows.append({
            "Name": t["name"],
            "Type": t["type"],
            "Mag": t["mag"] if t["mag"] is not None else "",
            "Best Local Time": t["best_time"].strftime("%H:%M"),  # time only
            "Best Alt (°)": t["best_alt"],
            "Best Dir": t["best_dir"],
            "Best Az (°)": t["best_az"],
            "Notes": t["notes"],
            "_Priority": t["priority"],
            "_BestDT": t["best_time"],
        })
    master_df = pd.DataFrame(master_rows)
    if not master_df.empty:
        master_df.sort_values(by=["_Priority","_BestDT"], ascending=[False, True], inplace=True)
        master_df.drop(columns=["_Priority","_BestDT"], inplace=True)

    # Hourly tables (rank by interest within hour)
    hour_tables = []
    for t in targets:
        for row in t["hourly"]:
            h, alt_h, az_h, dir_h, prio, *rest = row
            out = {
                "Hour": h.strftime("%Y-%m-%d %H:%M"),
                "Name": t["name"],
                "Type": t["type"],
                "Alt (°)": round(alt_h,1),
                "Dir": dir_h,
                "Az (°)": round(az_h,1),
                "_Priority": prio,
            }
            hour_tables.append(out)

    hourly_df = pd.DataFrame(hour_tables)
    if not hourly_df.empty:
        hourly_df["Hour_dt"] = pd.to_datetime(hourly_df["Hour"])
        hourly_df.sort_values(["Hour_dt","_Priority"], ascending=[True,False], inplace=True)
        # keep top-N by interest per hour
        hourly_df["rank"] = hourly_df.groupby("Hour")["_Priority"].rank(method="first", ascending=False)
        hourly_df = hourly_df[hourly_df["rank"] <= args.top_n_per_hour]
        hourly_df.drop(columns=["rank","Hour_dt","_Priority"], inplace=True, errors="ignore")

    # Save CSVs
    master_csv = f"{args.out_prefix}_master.csv"
    hourly_csv = f"{args.out_prefix}_hourly.csv"
    master_df.to_csv(master_csv, index=False)
    hourly_df.to_csv(hourly_csv, index=False)

    # HTML output
    if args.html:
        write_html(args.html, args.lat, args.lon, args.tz,
                   args.date, args.start, args.end,
                   master_df, hourly_df,
                   ui_mode=args.html_ui)

    # Console
    print("\n=== MASTER LIST (by interestingness) ===")
    if not master_df.empty:
        print(master_df.to_string(index=False))
    else:
        print("No targets passed the filters. Try lowering min_alt or max_mag, or extend the time window.")

    print("\n=== HOURLY 'POINT YOUR SCOPE' TABLES (by interest) ===")
    if not hourly_df.empty:
        for hour, sub in hourly_df.groupby("Hour"):
            print(f"\n# {hour}")
            print(sub[["Name","Type","Alt (°)","Dir","Az (°)"]].to_string(index=False))
    else:
        print("No hourly targets above your altitude threshold.")

if __name__ == "__main__":
    args = parse_args()
    plan_for_site(args)
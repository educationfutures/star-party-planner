Star Party Planner — quick start

Install once:
  pip install skyfield numpy pandas pytz python-dateutil

Example (Eagle Lake Observatory):
  python starparty_planner.py     --lat 44.810265 --lon -93.939783 --elev 296     --date 2025-08-28 --start 20:00 --end 01:00 --tz America/Chicago     --catalog objects_sample.csv --min_alt 20 --moon_sep_min 20 --max_mag 9     --top_n_per_hour 12 --out_prefix starparty     --html starparty.html

Outputs:
  - starparty_master.csv   (best-time list)
  - starparty_hourly.csv   (per-hour shortlist)
  - starparty.html         (night-vision: red text on black)

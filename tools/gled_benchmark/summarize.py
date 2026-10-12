"""Summaries of results.tsv: per-target table and per-method scores."""
import csv
import math
import sys


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def load(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stats(rows, key_est, two_sided_only=True):
    """Scores of the point estimate in `key_est` over rows with a window."""
    abs_log, signed, abs_yr, inwin = [], [], [], 0
    n = 0
    for r in rows:
        lo, hi, est = f(r["lo"]), f(r["hi"]), f(r[key_est])
        if None in (lo, hi, est) or est <= 0:
            continue
        mid = f(r["center"]) if f(r.get("center")) else math.sqrt(lo * hi)
        n += 1
        e = math.log(est / mid)
        abs_log.append(abs(e)); signed.append(e)
        abs_yr.append(abs(est - mid))
        inwin += lo <= est <= hi
    if not n:
        return None
    srt = sorted(abs_log)
    return dict(n=n, med_abs_log=srt[len(srt) // 2], mean_signed_log=sum(signed) / n,
                mae_yr=sum(abs_yr) / n, rmse_yr=math.sqrt(sum(x * x for x in abs_yr) / n),
                in_window=inwin)


def interval_stats(rows):
    n = cov = cont = 0
    wid, iscore = [], []
    for r in rows:
        lo, hi, a, b = f(r["lo"]), f(r["hi"]), f(r["hpd_lo"]), f(r["hpd_hi"])
        if None in (lo, hi, a, b):
            continue
        n += 1
        cov += a <= hi and b >= lo
        cont += a <= lo and b >= hi
        wid.append(b - a)
        sc = f(r["iscore"])
        if sc is not None:
            iscore.append(sc)
    return dict(n=n, overlap=cov, contains=cont,
                med_width=sorted(wid)[len(wid) // 2] if wid else None,
                mean_iscore=sum(iscore) / len(iscore) if iscore else None)


def main(path):
    rows = load(path)
    two = [r for r in rows if r["use"] in ("interval", "point")]
    for v in sorted({r["variant"] for r in rows}):
        rv = [r for r in two if r["variant"] == v]
        print(f"## clepsydra, {v}: {stats(rv, 'median')}  {interval_stats(rv)}")
    ids = {}
    for r in two:
        ids.setdefault(r["id"], r)
    print("## ASJP (Holman 2011 constants):", stats(list(ids.values()), "asjp_date"))
    for r in rows:
        pass


if __name__ == "__main__":
    main(sys.argv[1])

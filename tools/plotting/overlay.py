"""
Prior-only overlay: each node's 90% HPD under the prior alone against its
posterior, from the `prior_*` columns of one `date-tree` run.

The prior here is the full prior of that run (tree prior, calibrations, clock
prior) with the cognate likelihood switched off. Where the two bars coincide,
the data did not move the node; its age is whatever the calibrations and the
age prior say.

Writes an SVG (stdlib only) and prints a Markdown table.

Usage:
  python3 tools/plotting/overlay.py examples/romance_tree/output/posterior_nodes.tsv \\
      --svg examples/romance_tree/prior_overlay.svg --title "Romance"
"""
import argparse
import sys

ROW_H = 30
LEFT = 170
RIGHT = 30
TOP = 62
BOTTOM = 48
WIDTH = 760

PRIOR = "#9aa3ad"
POST = "#1f6fb4"
INK = "#1d232a"
GRID = "#e3e7eb"


def read_nodes(path):
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        rows = []
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            r = dict(zip(header, parts))
            rows.append({
                "node": r["node"],
                "med": float(r["median_bp"]),
                "lo": float(r["hpd90_lo_bp"]),
                "hi": float(r["hpd90_hi_bp"]),
                "pmed": float(r["prior_median_bp"]),
                "plo": float(r["prior_hpd90_lo_bp"]),
                "phi": float(r["prior_hpd90_hi_bp"]),
                "overlap": float(r["prior_overlap"]),
            })
    return rows


def nice_step(span):
    raw = span / 6.0
    for e in range(0, 7):
        for m in (1, 2, 2.5, 5):
            step = m * 10 ** e
            if step >= raw:
                return step
    return raw


def svg(rows, title):
    xmax = max(max(r["hi"], r["phi"]) for r in rows)
    step = nice_step(xmax)
    xmax = step * (int(xmax / step) + 1)
    plot_w = WIDTH - LEFT - RIGHT
    height = TOP + ROW_H * len(rows) + BOTTOM

    def x(bp):  # older to the left, as on a time tree
        return LEFT + plot_w * (1.0 - bp / xmax)

    out = []
    w = out.append
    w('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
      'viewBox="0 0 %d %d" font-family="Helvetica, Arial, sans-serif" '
      'font-size="12">' % (WIDTH, height, WIDTH, height))
    w('<rect width="100%" height="100%" fill="#ffffff"/>')
    w('<text x="%d" y="20" font-size="14" font-weight="bold" fill="%s">%s</text>'
      % (LEFT, INK, title))
    # legend
    lx = LEFT
    w('<rect x="%d" y="31" width="24" height="7" fill="%s"/>' % (lx, PRIOR))
    w('<text x="%d" y="39" fill="%s">prior only</text>' % (lx + 30, INK))
    w('<rect x="%d" y="31" width="24" height="7" fill="%s"/>' % (lx + 110, POST))
    w('<text x="%d" y="39" fill="%s">posterior</text>' % (lx + 140, INK))
    w('<text x="%d" y="39" fill="%s">bars: 90%% HPD; ticks: median</text>'
      % (lx + 230, INK))
    # grid and axis
    y_axis = TOP + ROW_H * len(rows) + 6
    t = 0.0
    while t <= xmax + 1e-9:
        xt = x(t)
        w('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s"/>'
          % (xt, TOP - 4, xt, y_axis, GRID))
        w('<text x="%.1f" y="%d" text-anchor="middle" fill="%s">%s</text>'
          % (xt, y_axis + 16, INK, "{:,}".format(int(t))))
        t += step
    w('<text x="%.1f" y="%d" text-anchor="middle" fill="%s">years BP</text>'
      % (LEFT + plot_w / 2.0, y_axis + 34, INK))
    for i, r in enumerate(rows):
        yc = TOP + ROW_H * i + ROW_H / 2.0
        w('<text x="%d" y="%.1f" text-anchor="end" fill="%s">%s</text>'
          % (LEFT - 10, yc + 4, INK, r["node"]))
        w('<rect x="%.1f" y="%.1f" width="%.1f" height="7" fill="%s"/>'
          % (x(r["phi"]), yc - 9, max(x(r["plo"]) - x(r["phi"]), 1.0), PRIOR))
        w('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
          'stroke-width="2"/>' % (x(r["pmed"]), yc - 11, x(r["pmed"]), yc,
                                  INK))
        w('<rect x="%.1f" y="%.1f" width="%.1f" height="7" fill="%s"/>'
          % (x(r["hi"]), yc + 2, max(x(r["lo"]) - x(r["hi"]), 1.0), POST))
        w('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" '
          'stroke-width="2"/>' % (x(r["med"]), yc, x(r["med"]), yc + 11, INK))
    w('</svg>')
    return "\n".join(out) + "\n"


def table(rows):
    lines = ["| node | prior median (90% HPD) | posterior median (90% HPD) "
             "| overlap |", "|---|---|---|---|"]
    for r in rows:
        lines.append("| %s | %s (%s–%s) | %s (%s–%s) | %.2f |" % (
            r["node"],
            "{:,.0f}".format(r["pmed"]), "{:,.0f}".format(r["plo"]),
            "{:,.0f}".format(r["phi"]),
            "{:,.0f}".format(r["med"]), "{:,.0f}".format(r["lo"]),
            "{:,.0f}".format(r["hi"]), r["overlap"]))
    return "\n".join(lines) + "\n"


def main(argv):
    ap = argparse.ArgumentParser(prog="overlay")
    ap.add_argument("nodes", help="posterior_nodes.tsv from date-tree")
    ap.add_argument("--svg", required=True)
    ap.add_argument("--title", default="Prior only against posterior")
    args = ap.parse_args(argv[1:])
    rows = read_nodes(args.nodes)
    with open(args.svg, "w", encoding="utf-8") as f:
        f.write(svg(rows, args.title))
    sys.stdout.write(table(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

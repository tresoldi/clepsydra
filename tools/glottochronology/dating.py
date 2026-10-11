"""
Classical lexicostatistical dating of a fixed tree, for comparison with the
Bayesian engine: Swadesh (1955), Lees (1953), Starostin (2000) and a
self-calibrated rate.

Every method turns a proportion C of shared cognates into a time. For a
node, C is pooled over every pair of tips on either side of it (matches
summed over compared concepts), the classical practice of dating a split by
its cross-clade comparisons. Intervals come from a bootstrap over concepts,
which keeps the dependence between pairs that share a concept; they cover
sampling error in C only, not error in the rate constant.

Dated tips. A tip of age tau has been separate for t - tau of the node's t
years. With total separation S = (t - tau_a) + (t - tau_b):
  Swadesh/Lees   C = r^S                       t = (ln C / ln r + tau_a + tau_b) / 2
  Starostin      C = exp(-lam C (S/2)^2 * 2)   t = sqrt(ln C / (-2 lam C)) + (tau_a + tau_b) / 2
Starostin's form is nonlinear in time, so the tip correction uses the mean
separation S/2; it is exact for contemporaneous tips. lam = 0.05 per
millennium per lineage (Starostin 2000, as formalised by Vasilyev 2009); the
widely quoted single-lambda form sqrt(ln C / (-0.05 C)) is reported too.

Stdlib only.

  python3 tools/glottochronology/dating.py --forms F --tree T --tip-age a=131 ...
"""
import argparse
import math
import random
import sys

# Swadesh (1955) 100-item list and the Lees (1953) / Swadesh 200-item list,
# as concept IDs of the Arawan forms (C_<GLOSS>). Only exact or standard
# equivalents are mapped (FINGERNAIL = claw, FLESHORMEAT = flesh, FAT = grease,
# SOIL = earth, HILL = mountain, SPEAK = say, LIEREST = lie, STRIKEORBEAT =
# hit, EARTHWORM = worm, INSIDE = in); GO is left out (walk/come).
SWADESH_100 = """I YOU ALL MANY ONE TWO BIG LONG WOMAN MAN FISH DOG LOUSE TREE SEED
LEAF ROOT BARK SKIN FLESHORMEAT BLOOD BONE FAT EGG TAIL HAIR HEAD EAR NOSE
MOUTH TOOTH TONGUE FINGERNAIL FOOT KNEE HAND BELLY NECK BREAST HEART LIVER
DRINK EAT BITE SEE HEAR KNOW SLEEP DIE KILL LIEREST SIT STAND GIVE SPEAK SUN
MOON STAR WATER RAIN STONE SAND SOIL CLOUD FIRE BURN PATH HILL RED YELLOW
WHITE BLACK NIGHT COLD FULL NEW GOOD DRY NAME""".split()
SWADESH_200_EXTRA = """CHILD WIFE HUSBAND MOTHER FATHER ANIMAL SNAKE EARTHWORM
STICK FRUIT FLOWER LEG BACK SUCK SPIT VOMIT BLOW BREATHE LAUGH SMELL
STRIKEORBEAT CUT DIG FALL SQUEEZE WASH PUSH TIE RIVER LAKE WIND OLD BAD
ROTTEN DIRTY WET NEAR HERE THREE HEAVY INSIDE""".split()
SWADESH_200 = SWADESH_100 + SWADESH_200_EXTRA

R_SWADESH = 0.86     # per millennium, 100-item list (Swadesh 1955)
R_LEES = 0.805       # per millennium, 200/215-item list (Lees 1953)
LAM_STAROSTIN = 0.05 # per millennium per lineage (Starostin 2000)


def read_forms(path):
    """language -> concept -> set of cognate ids (synonyms kept)."""
    out = {}
    with open(path, encoding="utf-8") as f:
        h = f.readline().rstrip("\n").split("\t")
        il, ic, ik = h.index("language_id"), h.index("concept_id"), h.index("cognate_id")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= ik or p[ik] in ("", "."):
                continue
            out.setdefault(p[il], {}).setdefault(p[ic], set()).add(p[ik])
    return out


def parse_newick(s):
    """Minimal labelled-Newick parser: returns (children, root)."""
    s = s.strip().rstrip(";")
    pos = 0
    children = {}
    anon = [0]

    def label():
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        lab = s[start:pos].split(":")[0]
        if pos < len(s) and s[pos] == ":":
            pos += 1
            while pos < len(s) and s[pos] not in ",();":
                pos += 1
        return lab

    def node():
        nonlocal pos
        kids = []
        if s[pos] == "(":
            pos += 1
            kids.append(node())
            while s[pos] == ",":
                pos += 1
                kids.append(node())
            pos += 1  # ')'
        lab = label()
        if not lab:
            anon[0] += 1
            lab = "node%d" % anon[0]
        children[lab] = kids
        return lab

    root = node()
    return children, root


def tips_under(children, v):
    if not children[v]:
        return [v]
    out = []
    for c in children[v]:
        out += tips_under(children, c)
    return out


def match_counts(forms, a, b, concepts):
    """(matches, compared): a concept matches when the two share a cognate."""
    m = n = 0
    for c in concepts:
        sa, sb = forms.get(a, {}).get(c), forms.get(b, {}).get(c)
        if sa and sb:
            n += 1
            m += 1 if sa & sb else 0
    return m, n


def date_swadesh(C, r, tau_a, tau_b):
    if C <= 0.0:
        return float("inf")
    if C >= 1.0:
        return (tau_a + tau_b) / 2.0
    return (math.log(C) / math.log(r) * 1000.0 + tau_a + tau_b) / 2.0


def date_starostin(C, lam, tau_a, tau_b, pair_factor=2.0):
    if C <= 0.0:
        return float("inf")
    if C >= 1.0:
        return (tau_a + tau_b) / 2.0
    return math.sqrt(math.log(C) / (-pair_factor * lam * C)) * 1000.0 + (tau_a + tau_b) / 2.0


def node_date(forms, children, v, concepts, method, tip_age):
    """Pool C over cross pairs, date with the pairs' mean tip ages."""
    kids = children[v]
    if len(kids) < 2:
        return None
    m = n = 0
    taus = []
    for i in range(len(kids)):
        for j in range(i + 1, len(kids)):
            for a in tips_under(children, kids[i]):
                for b in tips_under(children, kids[j]):
                    mm, nn = match_counts(forms, a, b, concepts)
                    if nn:
                        m += mm; n += nn
                        taus.append((tip_age.get(a, 0.0), tip_age.get(b, 0.0), nn))
    if n == 0:
        return None
    C = m / n
    w = sum(t[2] for t in taus)
    ta = sum(t[0] * t[2] for t in taus) / w
    tb = sum(t[1] * t[2] for t in taus) / w
    return C, n, method(C, ta, tb)


METHODS = {
    "swadesh100": (SWADESH_100, lambda C, a, b: date_swadesh(C, R_SWADESH, a, b)),
    "lees200": (SWADESH_200, lambda C, a, b: date_swadesh(C, R_LEES, a, b)),
    "starostin100": (SWADESH_100, lambda C, a, b: date_starostin(C, LAM_STAROSTIN, a, b, 2.0)),
    "starostin100_1lam": (SWADESH_100, lambda C, a, b: date_starostin(C, LAM_STAROSTIN, a, b, 1.0)),
}


def bootstrap(forms, children, v, concepts, method, tip_age, reps, rng):
    out = []
    for _ in range(reps):
        sample = [rng.choice(concepts) for _ in concepts]
        r = node_date(forms, children, v, sample, method, tip_age)
        if r is not None and math.isfinite(r[2]):
            out.append(r[2])
    out.sort()
    if not out:
        return float("nan"), float("nan")
    return out[int(0.05 * len(out))], out[min(len(out) - 1, int(0.95 * len(out)))]


def upgma(forms, langs, concepts, method, tip_age):
    """UPGMA on pairwise method dates: (newick, root date). The method's own
    tree, as glottochronology would build it without an external topology."""
    d = {}
    for i, a in enumerate(langs):
        for b in langs[i + 1:]:
            m, n = match_counts(forms, a, b, concepts)
            t = method(m / n, tip_age.get(a, 0.0), tip_age.get(b, 0.0)) if n else float("inf")
            d[(a, b)] = d[(b, a)] = t
    cl = {a: ([a], a) for a in langs}
    top = 0.0
    while len(cl) > 1:
        ks = list(cl)
        best = None
        for i, a in enumerate(ks):
            for b in ks[i + 1:]:
                x = sum(d[(p, q)] for p in cl[a][0] for q in cl[b][0])
                x /= len(cl[a][0]) * len(cl[b][0])
                if best is None or x < best[0]:
                    best = (x, a, b)
        x, a, b = best
        (ma, na), (mb, nb) = cl.pop(a), cl.pop(b)
        key = "(%s,%s):%.0f" % (na, nb, x)
        cl[key] = (ma + mb, key)
        top = x
    return list(cl.values())[0][1], top


def main(argv):
    ap = argparse.ArgumentParser(prog="glottochronology.dating")
    ap.add_argument("--forms", required=True)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--tip-age", action="append", default=[],
                    help="TIP=YEARS for a dated tip (repeatable)")
    ap.add_argument("--calibrate", nargs=3, metavar=("TIP_A", "TIP_B", "YEARS"),
                    help="also fit the Swadesh-form rate on one pair of known separation")
    ap.add_argument("--upgma", default="",
                    help="comma-separated tips: also build each method's own UPGMA tree")
    ap.add_argument("--all-concepts", action="store_true",
                    help="every method uses every concept (for data whose concept IDs"
                         " are not Arawan glosses; a list-size mismatch for Swadesh-100)")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv[1:])

    forms = read_forms(args.forms)
    children, root = parse_newick(open(args.tree, encoding="utf-8").read())
    tip_age = {}
    for t in args.tip_age:
        k, v = t.split("=")
        tip_age[k] = float(v)
    all_concepts = sorted({c for l in forms.values() for c in l})
    present = set(all_concepts)

    methods = dict(METHODS)
    if args.calibrate:
        a, b, years = args.calibrate[0], args.calibrate[1], float(args.calibrate[2])
        m, n = match_counts(forms, a, b, all_concepts)
        C = m / n
        # Total separation of the pair: (years - tau_a) + (years - tau_b).
        sep = (2.0 * years - tip_age.get(a, 0.0) - tip_age.get(b, 0.0)) / 1000.0
        r_cal = C ** (1.0 / sep) if C < 1.0 else 1.0
        print("# self-calibration on %s-%s: C = %d/%d = %.3f over %.0f years of separation"
              " -> r = %.4f per millennium" % (a, b, m, n, C, sep * 1000, r_cal))
        methods["selfcal_all"] = (all_concepts,
                                  lambda C, x, y, r=r_cal: date_swadesh(C, r, x, y))

    if args.all_concepts:
        methods = {k: (all_concepts, fn) for k, (_, fn) in methods.items()}
    rng = random.Random(args.seed)
    internal = [v for v in children if children[v]]
    print("method\tnode\tn_concepts_list\tn_compared\tC\tdate_bp\tboot90_lo\tboot90_hi")
    for name, (clist, fn) in methods.items():
        concepts = [c if c in present else "C_" + c for c in clist]
        concepts = [c for c in concepts if c in present]
        for v in internal:
            r = node_date(forms, children, v, concepts, fn, tip_age)
            if r is None:
                continue
            lo, hi = bootstrap(forms, children, v, concepts, fn, tip_age, args.reps, rng)
            print("%s\t%s\t%d\t%d\t%.4f\t%.0f\t%.0f\t%.0f"
                  % (name, v, len(concepts), r[1], r[0], r[2], lo, hi))
    if args.upgma:
        langs = args.upgma.split(",")
        print("\nmethod\tupgma_root_bp\tboot90_lo\tboot90_hi\tupgma_tree")
        for name, (clist, fn) in methods.items():
            concepts = [c if c in present else "C_" + c for c in clist]
            concepts = [c for c in concepts if c in present]
            tree, top = upgma(forms, langs, concepts, fn, tip_age)
            boots = []
            for _ in range(args.reps):
                sample = [rng.choice(concepts) for _ in concepts]
                t = upgma(forms, langs, sample, fn, tip_age)[1]
                if math.isfinite(t):
                    boots.append(t)
            boots.sort()
            print("%s\t%.0f\t%.0f\t%.0f\t%s;" % (name, top, boots[int(0.05 * len(boots))],
                  boots[min(len(boots) - 1, int(0.95 * len(boots)))], tree))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

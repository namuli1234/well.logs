"""Core petrophysics engine: LAS/CSV parsing, curve mapping, QC, calculations.
Every equation is documented; nothing is hidden or assumed silently."""
import re, csv
import numpy as np

ALIASES = {
    'GR': ['GR', 'GAM', 'GRC', 'CGR', 'SGR', 'GR_1'],
    'RT': ['RT', 'LLD', 'ILD', 'AT90', 'RDEP', 'RD', 'RLLD', 'RILD'],
    'RS': ['LLS', 'ILM', 'RS', 'AT30', 'RLLS'],
    'RXO': ['MSFL', 'RXO', 'SFL', 'RMSFL'],
    'RHOB': ['RHOB', 'DEN', 'RHOZ', 'RHOM', 'DENS'],
    'NPHI': ['NPHI', 'TNPH', 'NEU', 'NPOR', 'CNL'],
    'DT': ['DT', 'DTC', 'AC', 'DTCO', 'DT4P'],
    'PEF': ['PEF', 'PE', 'PEFZ'],
    'CALI': ['CALI', 'CAL', 'HCAL', 'CALS'],
    'SP': ['SP', 'SPR'],
}
MATRIX = {'Sandstone': 2.65, 'Limestone': 2.71, 'Dolomite': 2.87}


def map_curves(names, overrides=None):
    """Return {standard: actual_mnemonic}. Overrides allow manual mapping."""
    out, up = {}, {n.upper(): n for n in names}
    for std, al in ALIASES.items():
        for a in al:
            if a in up:
                out[std] = up[a]
                break
    out.update(overrides or {})
    return out


def parse_las(path):
    sec, meta, curves, tokens = None, {}, [], []
    with open(path, errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            if s.startswith('~'):
                sec = s[1:2].upper()
                continue
            if sec == 'A':
                tokens.extend(s.split())
                continue
            m = re.match(r'^([^.\s:]+)\s*\.([^\s:]*)\s*(.*?)\s*:\s*(.*)$', s)
            if not m:
                continue
            mn, un, val, desc = [x.strip() for x in m.groups()]
            if sec == 'C':
                curves.append((mn.upper(), un, desc))
            elif sec in ('V', 'W', 'P'):
                meta[mn.upper()] = (val, un)
    if not curves or not tokens:
        raise ValueError('No curve definitions or ASCII data found in LAS file.')
    n = len(curves)
    arr = np.array(tokens[: len(tokens) // n * n], dtype=float).reshape(-1, n)
    try:
        null = float(meta.get('NULL', ('-999.25',))[0])
    except ValueError:
        null = -999.25
    arr[np.isclose(arr, null)] = np.nan
    names, seen = [], {}
    for c in curves:
        k = seen.get(c[0], 0)
        seen[c[0]] = k + 1
        names.append(c[0] if k == 0 else f'{c[0]}_{k}')
    return dict(meta=meta, names=names[1:], units={nm: c[1] for nm, c in zip(names, curves)},
                depth=arr[:, 0], depth_unit=curves[0][1] or 'm',
                data={nm: arr[:, i] for i, nm in enumerate(names) if i > 0},
                well=meta.get('WELL', ('UNKNOWN',))[0])


def parse_csv(path):
    with open(path, newline='', errors='ignore') as f:
        rows = list(csv.reader(f))
    head, body = [h.strip().upper() for h in rows[0]], rows[1:]

    def fl(x):
        try:
            return float(x)
        except ValueError:
            return np.nan
    arr = np.array([[fl(x) for x in r[: len(head)]] + [np.nan] * (len(head) - len(r)) for r in body])
    di = next((i for i, h in enumerate(head) if h in ('DEPT', 'DEPTH', 'MD')), 0)
    names = [h for i, h in enumerate(head) if i != di]
    return dict(meta={}, names=names, units={}, depth=arr[:, di], depth_unit='m',
                data={h: arr[:, i] for i, h in enumerate(head) if i != di}, well=path.split('/')[-1])


def norm_nphi(x):
    """Neutron porosity as v/v. If values look like percent (max > 1.5), divide by 100."""
    x = np.asarray(x, float)
    return x / 100.0 if np.nanmax(x) > 1.5 else x


RANGES = {'GR': (0, 500), 'RT': (0.01, 1e5), 'RHOB': (1.0, 3.3), 'NPHI': (-0.15, 1.0), 'DT': (40, 250)}


def qc_report(w, mapped):
    d = w['depth']
    ok = ~np.isnan(d)
    L = ['DATA QUALITY REPORT', f"Well: {w['well']}",
         f"Depth coverage: {np.nanmin(d):.2f} - {np.nanmax(d):.2f} {w['depth_unit']}  ({ok.sum()} samples)"]
    dd = np.diff(d[ok])
    if len(dd):
        L.append(f"Median step: {np.median(dd):.4f}")
        if (dd < 0).any(): L.append('FLAG: depth is NOT monotonic increasing')
        if (dd == 0).any(): L.append(f'FLAG: {(dd == 0).sum()} duplicate depth values')
        if np.ptp(np.round(dd, 4)) > 1e-3: L.append('FLAG: irregular depth sampling')
    L.append('\nCurve validity (flags only - nothing is deleted):')
    for nm in w['names']:
        v = w['data'][nm]
        pct = 100 * np.mean(~np.isnan(v))
        line = f'  {nm:8s} {pct:5.1f}% valid'
        std = next((k for k, a in mapped.items() if a == nm), None)
        if std in RANGES and pct > 0:
            lo, hi = RANGES[std]
            x = norm_nphi(v) if std == 'NPHI' else v
            bad = np.nansum((x < lo) | (x > hi))
            if bad: line += f'  | {int(bad)} samples outside plausible range [{lo}, {hi}]'
        L.append(line)
    L.append('\nCurve mapping: ' + (', '.join(f'{k}<-{v}' for k, v in mapped.items()) or 'none recognised'))
    return '\n'.join(L)


# ---------------- Calculations ----------------
def vsh(gr, grmin, grmax, method='Linear'):
    """IGR=(GR-GRmin)/(GRmax-GRmin). Larionov (tertiary) 0.083(2^(3.7 IGR)-1);
    Clavier 1.7-sqrt(3.38-(IGR+0.7)^2); Steiber IGR/(3-2 IGR)."""
    igr = np.clip((gr - grmin) / (grmax - grmin), 0, 1)
    if method == 'Larionov (tertiary)': v = 0.083 * (2 ** (3.7 * igr) - 1)
    elif method == 'Larionov (older)': v = 0.33 * (2 ** (2 * igr) - 1)
    elif method == 'Clavier': v = 1.7 - np.sqrt(3.38 - (igr + 0.7) ** 2)
    elif method == 'Steiber': v = igr / (3 - 2 * igr)
    else: v = igr
    return np.clip(v, 0, 1)


def phi_density(rhob, rma, rfl):
    """PHI_D = (rho_ma - rho_b)/(rho_ma - rho_fl), clipped to [0, 0.6]."""
    return np.clip((rma - rhob) / (rma - rfl), 0, 0.6)


def sw_archie(rt, phi, rw, a, m, n):
    """Sw = ((a*Rw)/(Rt*PHI^m))^(1/n), clipped to [0,1]."""
    with np.errstate(all='ignore'):
        rt = np.where(rt > 0, rt, np.nan); phi = np.where(phi > 0, phi, np.nan)
        return np.clip((a * rw / (rt * phi ** m)) ** (1.0 / n), 0, 1)


def sw_simandoux(rt, phi, v, rw, rsh, a, m, n=2.0):
    """Modified Simandoux (n=2): Sw = (a*Rw/(2*PHI^m)) * [sqrt((Vsh/Rsh)^2 + 4*PHI^m/(a*Rw*Rt)) - Vsh/Rsh]."""
    with np.errstate(all='ignore'):
        rt = np.where(rt > 0, rt, np.nan); phi = np.where(phi > 0, phi, np.nan)
        pm = phi ** m
        s = (a * rw / (2 * pm)) * (np.sqrt((v / rsh) ** 2 + 4 * pm / (a * rw * rt)) - v / rsh)
        return np.clip(s, 0, 1)


def sw_indonesia(rt, phi, v, rw, rsh, a, m, n):
    """Indonesia: 1/sqrt(Rt) = [Vsh^(1-Vsh/2)/sqrt(Rsh) + PHI^(m/2)/sqrt(a*Rw)] * Sw^(n/2)."""
    with np.errstate(all='ignore'):
        rt = np.where(rt > 0, rt, np.nan); phi = np.where(phi > 0, phi, np.nan)
        den = v ** (1 - v / 2) / np.sqrt(rsh) + phi ** (m / 2) / np.sqrt(a * rw)
        return np.clip(((1 / np.sqrt(rt)) / den) ** (2.0 / n), 0, 1)


def pay_flag(v, phi, sw, cut):
    """PAY = (Vsh<cut) AND (PHI>cut) AND (Sw<cut). NaN inputs -> False."""
    with np.errstate(invalid='ignore'):
        return (v < cut['vsh']) & (phi > cut['phi']) & (sw < cut['sw'])


def zone_stats(depth, res, top, base):
    m = (depth >= top) & (depth <= base)
    step = np.nanmedian(np.diff(depth))
    gross = base - top
    net = res['PAY'][m].sum() * step
    def st(k):
        x = res[k][m]; x = x[~np.isnan(x)]
        return (x.mean(), x.min(), x.max()) if len(x) else (np.nan,) * 3
    return dict(gross=gross, net=net, ntg=net / gross if gross > 0 else np.nan,
                vsh=st('VSH'), phi=st('PHI'), sw=st('SW'))

#!/usr/bin/env python3
"""Mesure de la précision de tildes.py sur des observations marégraphiques.

Outil de développement (numpy). Les PM et BM « observées » sont obtenues sans
recourir à aucun modèle : autour de chaque maximum (minimum) discret des
mesures au pas de 10 min, un polynôme de degré 6 est ajusté sur ± 2 h et
son extremum donne l'heure et la hauteur observées. Elles sont appariées aux
PM et BM calculées par tildes.extrema_periode (le code réellement utilisé),
puis on mesure les écarts d'heure et de hauteur.

La hauteur observée contient la surcote météorologique ; on donne aussi
l'écart « hors surcote », où la surcote est estimée par la moyenne glissante
sur 25 h du résidu observé − prédit.
"""

import math
import os
import sys
from datetime import timedelta

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import harmonique as hm  # noqa: E402
import tildes  # noqa: E402

DEMI_FENETRE = 120.0 / 1440.0  # jours
DEGRE = 6


def extrema_observes(t, h, approx):
    """[(jour, type, hauteur)] observés près des instants approximatifs `approx`.

    approx : [(jour J2000, 'PM'|'BM')], p. ex. les extrema d'un modèle ; ils ne
    servent qu'à localiser la recherche (± 3 h), le résultat ne dépend que des
    mesures.
    """
    sortie = []
    for ja, typ in approx:
        i0, i1 = np.searchsorted(t, [ja - 0.125, ja + 0.125])
        if i1 - i0 < 25:
            continue
        seg = h[i0:i1]
        k = i0 + (np.argmax(seg) if typ == "PM" else np.argmin(seg))
        tc = t[k]
        for _ in range(2):
            j0, j1 = np.searchsorted(t, [tc - DEMI_FENETRE, tc + DEMI_FENETRE])
            x = (t[j0:j1] - tc) / DEMI_FENETRE
            if j1 - j0 < 15 or np.max(np.diff(t[j0:j1])) > 15 / 1440:  # trou
                tc = None
                break
            c = np.polynomial.polynomial.polyfit(x, h[j0:j1], DEGRE)
            xs = np.linspace(-0.8, 0.8, 1601)
            v = np.polynomial.polynomial.polyval(xs, c)
            m = np.argmax(v) if typ == "PM" else np.argmin(v)
            if m in (0, len(xs) - 1):
                tc = None
                break
            tc, hc = tc + xs[m] * DEMI_FENETRE, v[m]
        if tc is not None:
            sortie.append((tc, typ, hc))
    return sortie


def extrema_modele(port, j0, j1):
    """[(jour, type, hauteur, coeff)] calculés par tildes entre deux jours J2000."""
    debut = hm.datetime_depuis_jours(j0)
    fin = hm.datetime_depuis_jours(j1)
    res, a = [], debut
    while a < fin:  # par mois (mémoire / coût des listes Brest)
        b = min(fin, a + timedelta(days=31))
        for e in tildes.extrema_periode(a, b, port):
            j = (e["heure"] - tildes.J2000).total_seconds() / 86400.0
            res.append((j, e["type"], e["hauteur"], e["coeff"]))
        a = b
    return res


def surcote(t, h, port):
    """Résidu observé − prédit lissé sur 25 h (m), aux instants t."""
    r = h - hm.predire(t, port.z0, port.constantes)
    n = 151  # 25 h au pas de 10 min
    noyau = np.ones(n) / n
    lisse = np.convolve(r, noyau, mode="same")
    poids = np.convolve(np.ones_like(r), noyau, mode="same")
    return lisse / poids, r


def comparer(port, t, h, j0, j1, obs=None, mod=None):
    """Statistiques d'écart modèle − observation sur [j0, j1]."""
    sel = (t >= j0) & (t < j1)
    t, h = t[sel], h[sel]
    if mod is None:
        mod = extrema_modele(port, j0, j1)
    if obs is None:
        obs = extrema_observes(t, h, [(j, typ) for j, typ, _, _ in mod])
    s, r = surcote(t, h, port)
    r_hf = r - s
    # appariement : extremum calculé de même type le plus proche (< 2 h)
    jm = np.array([m[0] for m in mod])
    lignes = []
    coeff_courant = None
    coeffs = {}
    for j, typ, hh, c in mod:
        if typ == "PM":
            coeff_courant = c
        coeffs[j] = coeff_courant if c is None else c
    for jo, typ, ho in obs:
        i = int(np.argmin(np.abs(jm - jo)))
        jmod, tmod, hmod, _ = mod[i]
        if tmod != typ or abs(jmod - jo) > 2 / 24:
            continue
        k = min(np.searchsorted(t, jo), len(t) - 1)
        lignes.append((typ, (jmod - jo) * 1440.0, hmod - ho, hmod - (ho - s[k]), coeffs[jmod] or 0))
    return lignes, float(np.std(r_hf)), obs, mod


def resume(lignes, titre=""):
    def st(v):
        v = np.asarray(v)
        return v.mean(), math.sqrt(np.mean(v * v)), np.max(np.abs(v)) if len(v) else 0.0
    out = {}
    for typ in ("PM", "BM"):
        L = [l for l in lignes if l[0] == typ]
        dt = [l[1] for l in L]
        dh = [l[3] for l in L]
        out[typ] = (len(L), st(dt), st(dh), st([l[2] for l in L]),
                    np.mean(np.abs(dt) <= 5) * 100, np.mean(np.abs(dt) <= 10) * 100)
    return out


def afficher(lignes, sigma, titre):
    r = resume(lignes)
    print(f"--- {titre}  (σ résidu hors surcote : {sigma * 100:.2f} cm)")
    for typ in ("PM", "BM"):
        n, (bt, rt, mt), (bh, rh, mh), (bb, rb, _), p5, p10 = r[typ]
        print(f"  {typ} n={n:5d}  heure biais {bt:+5.2f} RMS {rt:5.2f} max {mt:5.1f} min"
              f" | ±5 min {p5:4.1f} % ±10 min {p10:4.1f} %"
              f" | hauteur hors surcote biais {bh * 100:+5.1f} RMS {rh * 100:4.1f} cm"
              f" (brute RMS {rb * 100:4.1f})")
    for nom, f in (("vives-eaux ≥ 75", lambda c: c >= 75), ("mortes-eaux < 55", lambda c: c < 55)):
        L = [l for l in lignes if f(l[4])]
        for typ in ("PM", "BM"):
            dt = np.array([l[1] for l in L if l[0] == typ])
            if len(dt):
                print(f"    {nom:17s} {typ} n={len(dt):4d} heure biais {dt.mean():+5.2f} RMS {math.sqrt(np.mean(dt * dt)):5.2f} min")


def port_depuis_json(chemin, base=tildes.CHERBOURG):
    """Port candidat enregistré par ajuster.py --json (définitions + constantes)."""
    import json
    r = json.load(open(chemin))
    for o, (d, regle) in r["definitions"].items():
        if o not in tildes.ONDES:
            if isinstance(regle, list):
                regle = tuple((m, k) for m, k in regle)
            tildes.ONDES[o] = (tuple(d), regle)
    cst = {o: tuple(v) for o, v in r["constantes"].items()}
    return tildes.Port(base.nom, base.latitude, base.longitude, r["z0"], cst,
                       inferer=False, exclues=frozenset())


def main():
    import argparse
    from datetime import datetime, timezone
    p = argparse.ArgumentParser(description="Écart de tildes.py aux observations REFMAR.")
    p.add_argument("debut", help="AAAA-MM-JJ")
    p.add_argument("fin", help="AAAA-MM-JJ (exclu)")
    p.add_argument("--station", type=int, default=13)
    p.add_argument("--candidat", action="append", default=[],
                   help="résultat de ajuster.py --json à comparer (répétable)")
    a = p.parse_args()
    j0, j1 = (hm.jours_depuis_j2000([datetime.fromisoformat(x).replace(tzinfo=timezone.utc)])[0]
              for x in (a.debut, a.fin))
    ans = range(int(a.debut[:4]), int(a.fin[:4]) + 1)
    t, h = hm.charger(a.station, 3, ans)
    t2, h2 = hm.charger(a.station, 2, ans)  # brutes, là où les validées manquent
    sel = t2 > (t[-1] if len(t) else -1e9)
    t, h = np.concatenate([t, t2[sel]]), np.concatenate([h, h2[sel]])
    port = {13: tildes.CHERBOURG, 3: tildes.BREST}[a.station]
    mod = extrema_modele(port, j0, j1)
    obs = extrema_observes(t[(t >= j0) & (t < j1)], h[(t >= j0) & (t < j1)], [(j, ty) for j, ty, _, _ in mod])
    lignes, sigma, _, _ = comparer(port, t, h, j0, j1, obs=obs, mod=mod)
    afficher(lignes, sigma, f"tildes.py ({port!r})")
    sel = (t >= j0) & (t < j1)
    print(f"  niveau moyen observé − Z0 : {np.mean(h[sel] - hm.predire(t[sel], port.z0, port.constantes)) * 100:+.1f} cm")
    for c in a.candidat:
        cand = port_depuis_json(c, port)
        lignes, sigma, _, _ = comparer(cand, t, h, j0, j1, obs=obs)
        afficher(lignes, sigma, f"{os.path.basename(c)} ({cand!r})")
        print(f"  niveau moyen observé − Z0 : {np.mean(h[sel] - hm.predire(t[sel], cand.z0, cand.constantes)) * 100:+.1f} cm")


if __name__ == "__main__":
    main()

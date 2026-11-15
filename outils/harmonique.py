#!/usr/bin/env python3
"""Analyse harmonique par moindres carrés d'observations marégraphiques.

Outil de développement (nécessite numpy) ; tildes.py n'en dépend pas.

Le modèle ajusté est exactement celui de tildes.py :

    h(t) = Z0 + b·(t − t_ref) + Σ f(t)·H·cos(V(t) + u(t) − G)

avec les mêmes arguments astronomiques V (Meeus, nombres de Doodson de
tildes.ONDES) et les mêmes corrections nodales f, u (OHI). Les constantes
(H, G) obtenues sont donc directement utilisables par tildes.Port.
"""

import math
import os
import sys
from datetime import datetime, timezone

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import tildes  # noqa: E402

J2000 = tildes.J2000


def jours_depuis_j2000(instants):
    """Tableau numpy de jours (float64) depuis J2000 pour des datetimes."""
    return np.array([(tildes.en_utc(t) - J2000).total_seconds() / 86400.0 for t in instants])


def arguments(jours):
    """Arguments (τ, s, h, p, N', p', 90°) en degrés, forme (7, n), et (N, p).

    Mêmes polynômes que tildes.arguments_astronomiques, vectorisés.
    """
    T = jours / 36525.0
    T2, T3, T4 = T * T, T ** 3, T ** 4
    s = 218.3164477 + 481267.88123421 * T - 0.0015786 * T2 + T3 / 538841.0 - T4 / 65194000.0
    h = 280.46646 + 36000.76983 * T + 0.0003032 * T2
    p = 83.3532465 + 4069.0137287 * T - 0.0103200 * T2 - T3 / 80053.0 + T4 / 18999000.0
    N = 125.04452 - 1934.136261 * T + 0.0020708 * T2 + T3 / 450000.0
    pp = 282.93735 + 1.71946 * T + 0.00046 * T2
    tau = np.mod(jours, 1.0) * 360.0 + h - s
    return np.vstack([tau, s, h, p, -N, pp, np.full_like(T, 90.0)]), np.mod(N, 360.0), np.mod(p, 360.0)


def nodal_par_jour(jours, ondes):
    """f et u (rad) de chaque onde, évalués au jour (midi) de chaque instant.

    Les corrections nodales varient sur 4,4 à 18,6 ans : une valeur par jour
    est exacte à mieux que 1e-4 près.
    """
    jour_entier = np.floor(jours + 0.5)  # jours J2000 entiers (midi)
    uniques, inverse = np.unique(jour_entier, return_inverse=True)
    _, N, p = arguments(uniques)
    F = np.empty((len(uniques), len(ondes)))
    U = np.empty((len(uniques), len(ondes)))
    for i in range(len(uniques)):
        nod = tildes.corrections_nodales(float(N[i]), float(p[i]))
        for j, o in enumerate(ondes):
            F[i, j], U[i, j] = nod[o]
    return F[inverse], np.radians(U[inverse])


def phases(jours, ondes):
    """V + u (rad) et f pour chaque instant et chaque onde : formes (n, k)."""
    args, _, _ = arguments(jours)
    D = np.array([tildes.ONDES[o][0] for o in ondes], dtype=float)  # (k, 7)
    V = np.radians(np.mod(D @ args, 360.0)).T  # (n, k)
    F, U = nodal_par_jour(jours, ondes)
    return V + U, F


def vitesses(ondes):
    return np.array([tildes.vitesse(o) for o in ondes])


def verifier_separation(ondes, duree_heures):
    """Lève une erreur si deux ondes ne sont pas séparables (critère de Rayleigh)."""
    w = vitesses(ondes)
    ordre = np.argsort(w)
    for a, b in zip(ordre, ordre[1:]):
        dw = abs(w[b] - w[a])
        if dw * duree_heures < 360.0:
            raise ValueError(f"{ondes[a]} et {ondes[b]} inséparables sur {duree_heures / 8766:.1f} ans "
                             f"(écart {dw:.6f} °/h)")


def analyser(jours, hauteurs, ondes, tendance=True, jour_ref=None, bloc=40000):
    """Ajuste Z0 (+ tendance) et (H, G) de chaque onde par moindres carrés.

    Retourne (z0 à jour_ref, tendance en m/an, {onde: (H m, G °)}, rms résiduel).
    Équations normales accumulées par blocs pour tenir en mémoire.
    """
    ondes = list(ondes)
    verifier_separation(ondes, (jours.max() - jours.min()) * 24.0)
    if jour_ref is None:
        jour_ref = 0.5 * (jours.min() + jours.max())
    k = len(ondes)
    ncol = 1 + int(tendance) + 2 * k
    AtA = np.zeros((ncol, ncol))
    Atb = np.zeros(ncol)
    for i in range(0, len(jours), bloc):
        A, y = _bloc(jours[i:i + bloc], hauteurs[i:i + bloc], ondes, tendance, jour_ref)
        AtA += A.T @ A
        Atb += A.T @ y
    x = np.linalg.solve(AtA, Atb)
    # rms résiduel
    sse, n = 0.0, 0
    for i in range(0, len(jours), bloc):
        A, y = _bloc(jours[i:i + bloc], hauteurs[i:i + bloc], ondes, tendance, jour_ref)
        r = y - A @ x
        sse += float(r @ r)
        n += len(r)
    z0 = x[0]
    pente = x[1] * 365.25 if tendance else 0.0
    off = 1 + int(tendance)
    constantes = {}
    for j, o in enumerate(ondes):
        a, b = x[off + 2 * j], x[off + 2 * j + 1]
        constantes[o] = (math.hypot(a, b), math.degrees(math.atan2(b, a)) % 360.0)
    return z0, pente, constantes, math.sqrt(sse / n)


def _bloc(jours, hauteurs, ondes, tendance, jour_ref):
    ph, F = phases(jours, ondes)
    cols = [np.ones_like(jours)]
    if tendance:
        cols.append(jours - jour_ref)
    A = np.empty((len(jours), len(cols) + 2 * len(ondes)))
    for c, v in enumerate(cols):
        A[:, c] = v
    off = len(cols)
    # f·H·cos(V+u−G) = (H cos G)·f cos(V+u) + (H sin G)·f sin(V+u)
    A[:, off::2] = F * np.cos(ph)
    A[:, off + 1::2] = F * np.sin(ph)
    return A, hauteurs


def predire(jours, z0, constantes, bloc=100000):
    """Hauteurs prédites (m) aux instants donnés (jours J2000), sans tendance."""
    ondes = list(constantes)
    H = np.array([constantes[o][0] for o in ondes])
    G = np.radians([constantes[o][1] for o in ondes])
    sortie = np.empty_like(jours)
    for i in range(0, len(jours), bloc):
        ph, F = phases(jours[i:i + bloc], ondes)
        sortie[i:i + bloc] = z0 + (F * H * np.cos(ph - G)).sum(axis=1)
    return sortie


def charger(station, source, annees):
    """(jours J2000, hauteurs) des fichiers téléchargés par refmar.py."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import refmar
    t, h = [], []
    for a in annees:
        chemin = os.path.join(refmar.RACINE, f"{station}_s{source}_{a}.csv")
        if not os.path.exists(chemin):
            continue
        brut = np.loadtxt(chemin, delimiter=",", dtype=str, ndmin=2)
        if brut.size == 0:
            continue
        instants = np.char.rstrip(brut[:, 0], "Z").astype("datetime64[s]")
        t.append((instants - np.datetime64("2000-01-01T12:00:00")).astype(float) / 86400.0)
        h.append(brut[:, 1].astype(float))
    if not t:
        return np.empty(0), np.empty(0)
    t, h = np.concatenate(t), np.concatenate(h)
    ordre = np.argsort(t, kind="stable")
    t, h = t[ordre], h[ordre]
    garde = np.concatenate([[True], np.diff(t) > 1e-6])  # doublons
    return t[garde], h[garde]


def datetime_depuis_jours(j):
    return datetime.fromtimestamp((J2000.timestamp() + j * 86400.0), tz=timezone.utc)

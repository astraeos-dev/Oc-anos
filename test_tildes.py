#!/usr/bin/env python3
"""Tests de cohérence du calculateur (bibliothèque standard : python3 -m unittest)."""

import unittest
from datetime import date, datetime, timedelta, timezone
from zoneinfo import ZoneInfo

import tildes

PARIS = ZoneInfo("Europe/Paris")


class TestOndes(unittest.TestCase):
    def test_vitesses_de_reference(self):
        # Vitesses de la liste des constituantes de l'OHI (°/h).
        attendues = {
            "M2": 28.9841042, "S2": 30.0, "N2": 28.4397295, "K2": 30.0821373,
            "K1": 15.0410686, "O1": 13.9430356, "P1": 14.9589314, "Q1": 13.3986609,
            "M4": 57.9682085, "MS4": 58.9841042, "M6": 86.9523127, "SA": 0.0410686,
        }
        for onde, v in attendues.items():
            self.assertAlmostEqual(tildes.vitesse(onde), v, places=5, msg=onde)

    def test_ondes_composees_vitesse_somme_des_meres(self):
        for onde, (_, regle) in tildes.ONDES.items():
            if isinstance(regle, tuple) and onde in tildes.INFERENCES:
                somme = sum(k * tildes.vitesse(m) for m, k in regle)
                self.assertAlmostEqual(tildes.vitesse(onde), somme, places=5, msg=onde)

    def test_corrections_nodales_bornees(self):
        for N in range(0, 360, 15):
            nod = tildes.corrections_nodales(float(N), 37.0)
            self.assertTrue(0.96 < nod["M2"][0] < 1.04)
            self.assertTrue(0.88 < nod["K1"][0] < 1.12)
            self.assertTrue(0.80 < nod["O1"][0] < 1.20)
            self.assertTrue(abs(nod["M2"][1]) < 2.2)

    def test_inference_ne_remplace_pas_les_observations(self):
        port = tildes.CHERBOURG
        for onde, valeur in port.observees.items():
            if onde in tildes.ONDES_NON_RETENUES:
                self.assertNotIn(onde, port.constantes)
            else:
                self.assertEqual(port.constantes[onde], valeur)
        self.assertIn("MK4", port.constantes)
        self.assertNotIn("MK4", port.observees)


class TestExtrema(unittest.TestCase):
    def verifier_journee(self, jour):
        ex = tildes.extrema_journee(jour)
        self.assertIn(len(ex), (3, 4), f"{jour} : {len(ex)} extrema")
        types = [e["type"] for e in ex]
        for a, b in zip(types, types[1:]):
            self.assertNotEqual(a, b, f"{jour} : deux {a} consécutives")
        utc = [e["heure"].astimezone(timezone.utc) for e in ex]
        self.assertEqual(utc, sorted(utc))
        for e in ex:
            self.assertEqual(e["heure"].date(), jour)
        return ex

    def test_jours_de_changement_d_heure(self):
        for jour in (date(2026, 3, 29), date(2026, 10, 25), date(2027, 3, 28), date(2027, 10, 31)):
            self.verifier_journee(jour)

    def test_continuite_jour_apres_jour(self):
        debut = date(2027, 1, 1)
        par_jour = []
        for i in range(60):
            par_jour += self.verifier_journee(debut + timedelta(days=i))
        d = datetime(2027, 1, 1, tzinfo=PARIS)
        bloc = tildes.extrema_periode(d, d + timedelta(days=60))
        self.assertEqual(len(par_jour), len(bloc))
        for a, b in zip(par_jour, bloc):
            self.assertEqual(a["type"], b["type"])
            self.assertLess(abs((a["heure"] - b["heure"]).total_seconds()), 2.0)
        pm = [e["heure"] for e in par_jour if e["type"] == "PM"]
        for a, b in zip(pm, pm[1:]):
            self.assertTrue(timedelta(hours=11) < b - a < timedelta(hours=14))

    def test_extremum_reel(self):
        for e in tildes.extrema_journee(date(2027, 6, 14)):
            t = e["heure"]
            avant = tildes.hauteur_maree(t - timedelta(minutes=5))
            apres = tildes.hauteur_maree(t + timedelta(minutes=5))
            if e["type"] == "PM":
                self.assertGreater(e["hauteur"], max(avant, apres))
            else:
                self.assertLess(e["hauteur"], min(avant, apres))
            self.assertAlmostEqual(tildes.hauteur_maree(t), e["hauteur"], places=4)

    def test_courbe_coherente_avec_les_extrema(self):
        jour = date(2027, 6, 14)
        ex = tildes.extrema_journee(jour)
        d = datetime(2027, 6, 14, tzinfo=PARIS)
        h = [v for _, v in tildes.courbe(d, d + timedelta(days=1), pas_minutes=1)]
        self.assertAlmostEqual(max(h), max(e["hauteur"] for e in ex), delta=0.005)
        self.assertAlmostEqual(min(h), min(e["hauteur"] for e in ex), delta=0.005)


class TestCoefficient(unittest.TestCase):
    def test_bornes_et_vive_eau_d_equinoxe(self):
        coeffs = []
        for i in range(30):
            coeffs += [e["coeff"] for e in tildes.extrema_journee(date(2027, 9, 1) + timedelta(days=i))
                       if e["type"] == "PM"]
        self.assertTrue(all(20 <= c <= 120 for c in coeffs))
        self.assertGreaterEqual(max(coeffs), 100)  # vives-eaux d'équinoxe
        self.assertLessEqual(min(coeffs), 40)      # mortes-eaux

    def test_amplitude_unite(self):
        self.assertEqual(tildes.coefficient_depuis_amplitude(tildes.UNITE_HAUTEUR_BREST), 100)
        self.assertEqual(tildes.coefficient_depuis_amplitude(0.1), 20)
        self.assertEqual(tildes.coefficient_depuis_amplitude(9.0), 120)


class TestCompatibilite(unittest.TestCase):
    def test_datetime_naif_lu_en_utc(self):
        naif = datetime(2027, 2, 1, 12, 0)
        self.assertEqual(tildes.hauteur_maree(naif),
                         tildes.hauteur_maree(naif.replace(tzinfo=timezone.utc)))

    def test_ancienne_fonction(self):
        jour = date(2027, 2, 1)
        self.assertEqual(len(tildes.trouver_extrema_journee(jour)), len(tildes.extrema_journee(jour)))


if __name__ == "__main__":
    unittest.main()

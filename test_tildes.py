#!/usr/bin/env python3
"""Tests de cohérence du calculateur (bibliothèque standard : python3 -m unittest)."""

import unittest
from datetime import date, datetime, timedelta, timezone
from unittest import mock
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import tildes

PARIS = tildes.fuseau("Europe/Paris")


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
        port = tildes.BREST_TICON
        for onde, valeur in port.observees.items():
            if onde in tildes.ONDES_NON_RETENUES:
                self.assertNotIn(onde, port.constantes)
            else:
                self.assertEqual(port.constantes[onde], valeur)
        self.assertIn("MK4", port.constantes)
        self.assertNotIn("MK4", port.observees)

    def test_satellites_nodaux(self):
        # Même onde que la mère, à ±1 sur N' près, sans correction nodale.
        for onde, (coeffs, regle) in tildes.SATELLITES_NODAUX.items():
            mere, signe = onde[:-2], 1 if onde.endswith("+N") else -1
            attendu = list(tildes.ONDES[mere][0])
            attendu[4] += signe
            self.assertEqual(tuple(attendu), coeffs, onde)
            self.assertIsNone(regle)

    def test_raies_empiriques_composees(self):
        # La règle nodale d'une raie est une combinaison d'ondes mères de même
        # fréquence (τ, s, h, p).
        for onde, (coeffs, regle) in tildes.RAIES_EMPIRIQUES.items():
            if regle is None:
                continue
            somme = [sum(k * tildes.ONDES[m][0][i] for m, k in regle) for i in range(4)]
            self.assertEqual(tuple(somme), coeffs[:4], onde)

    def test_ondes_separables(self):
        # Les ondes de chaque port sont séparables sur la durée des
        # observations analysées (critère de Rayleigh, ≥ 26 ans).
        for port in (tildes.CHERBOURG, tildes.BREST):
            w = sorted(tildes.vitesse(o) for o in port.constantes)
            ecart = min(b - a for a, b in zip(w, w[1:]))
            self.assertGreater(ecart * 26 * 8766, 360.0, port.nom)

    def test_ondes_definies(self):
        for port in (tildes.CHERBOURG, tildes.BREST, tildes.BREST_TICON):
            self.assertTrue(set(port.constantes) <= set(tildes.ONDES))
            self.assertTrue(3.5 < port.z0 < 4.5)


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

    def test_decalage_brest(self):
        self.assertAlmostEqual(tildes.decalage_brest(tildes.BREST), 0.0, delta=0.1)
        self.assertAlmostEqual(tildes.decalage_brest(tildes.CHERBOURG), 4.1, delta=0.2)

    def test_coefficient_identique_a_brest_et_cherbourg(self):
        # Une PM de Cherbourg reçoit le coefficient de la PM de Brest de la
        # même onde de marée, environ 4 h plus tôt.
        brest = [e for e in tildes.extrema_journee(date(2027, 3, 10), tildes.BREST) if e["type"] == "PM"]
        for e in brest:
            t = e["heure"] + timedelta(hours=tildes.decalage_brest(tildes.CHERBOURG))
            proche = min((x for x in tildes.extrema_periode(t - timedelta(hours=3), t + timedelta(hours=3))
                          if x["type"] == "PM"), key=lambda x: abs(x["heure"] - t))
            self.assertEqual(proche["coeff"], e["coeff"])

    def test_amplitude_unite(self):
        self.assertEqual(tildes.coefficient_depuis_amplitude(tildes.UNITE_HAUTEUR_BREST), 100)
        self.assertEqual(tildes.coefficient_depuis_amplitude(0.1), 20)
        self.assertEqual(tildes.coefficient_depuis_amplitude(9.0), 120)


class TestHeureLegale(unittest.TestCase):
    """Repli sans base IANA (Python sous Windows sans le paquet tzdata)."""

    F = tildes.HEURE_LEGALE_FRANCE

    def test_bascules_2026(self):
        u = lambda *a: datetime(*a, tzinfo=timezone.utc)
        self.assertEqual(u(2026, 3, 29, 0, 59).astimezone(self.F).utcoffset(), timedelta(hours=1))
        self.assertEqual(u(2026, 3, 29, 1, 0).astimezone(self.F).utcoffset(), timedelta(hours=2))
        self.assertEqual(u(2026, 10, 25, 0, 59).astimezone(self.F).utcoffset(), timedelta(hours=2))
        self.assertEqual(u(2026, 10, 25, 1, 0).astimezone(self.F).utcoffset(), timedelta(hours=1))
        self.assertEqual(u(2026, 7, 1).astimezone(self.F).tzname(), "CEST")

    def test_identique_a_zoneinfo(self):
        try:
            iana = ZoneInfo("Europe/Paris")
        except ZoneInfoNotFoundError:
            self.skipTest("base IANA absente")
        t = datetime(2020, 1, 1, tzinfo=timezone.utc)
        while t < datetime(2031, 1, 1, tzinfo=timezone.utc):
            a, b = t.astimezone(self.F), t.astimezone(iana)
            self.assertEqual((a.replace(tzinfo=None), a.fold, a.utcoffset()),
                             (b.replace(tzinfo=None), b.fold, b.utcoffset()), t)
            t += timedelta(minutes=30)
        for jour in (date(2026, 3, 29), date(2026, 10, 25), date(2027, 3, 28), date(2027, 10, 31)):
            for minute in range(0, 24 * 60, 15):
                for pli in (0, 1):
                    loc = datetime.combine(jour, datetime.min.time()) + timedelta(minutes=minute)
                    self.assertEqual(loc.replace(tzinfo=self.F, fold=pli).utcoffset(),
                                     loc.replace(tzinfo=iana, fold=pli).utcoffset(), (loc, pli))

    def test_marees_sans_base_iana(self):
        with mock.patch.object(tildes, "ZoneInfo", side_effect=ZoneInfoNotFoundError):
            self.assertIs(tildes.fuseau(), self.F)
            with self.assertRaises(ZoneInfoNotFoundError):
                tildes.fuseau("Europe/London")
            for jour in (date(2026, 3, 29), date(2026, 10, 25), date(2027, 6, 14)):
                sans = tildes.extrema_journee(jour)
                self.assertTrue(all(e["heure"].tzinfo is self.F for e in sans))
                with mock.patch.object(tildes, "fuseau", return_value=PARIS):
                    avec = tildes.extrema_journee(jour)
                self.assertEqual([(e["type"], e["heure"].astimezone(timezone.utc), e["coeff"]) for e in sans],
                                 [(e["type"], e["heure"].astimezone(timezone.utc), e["coeff"]) for e in avec])
                self.assertEqual([e["heure"].strftime("%Hh%M") for e in sans],
                                 [e["heure"].strftime("%Hh%M") for e in avec])


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

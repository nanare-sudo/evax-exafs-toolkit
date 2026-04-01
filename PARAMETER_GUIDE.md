# EvAX Parameter Guide

Alle Parameter, die in der YAML-Konfiguration und in `parameters.dat` vorkommen,
mit Erklaerung, empfohlenen Werten und Auswirkungen.

---

## Vergleichsraum

| Parameter | Werte | Empfehlung | Erklaerung |
|-----------|-------|------------|------------|
| `Space` | `k`, `R`, `w` | **`w`** (Wavelet) | In welchem Raum Theorie und Experiment verglichen werden. `k` = direkt im k-Raum (schnell, aber wenig sensitiv). `R` = Fourier-Transformierte (sensitiver, Standard in konventioneller EXAFS-Analyse). `w` = Wavelet-Raum (2D in k und R gleichzeitig, **am sensitivsten**, aber am rechenintensivsten). |

## k-Bereich und Gewichtung

| Parameter | Typ | Empfehlung | Erklaerung |
|-----------|-----|------------|------------|
| `k_min` | Float (1/Angstrom) | **3** | Untere k-Grenze. Unter ~3 dominieren Hintergrundeffekte und die Normierung ist unsicher. `k_min=2` kann bei sehr sauberen Daten Vorteile bringen. |
| `k_max` | Float (1/Angstrom) | **11-14** | Obere k-Grenze. Abhaengig vom Signal-Rausch-Verhaeltnis. Schwere Elemente (Pt) haben staerkeres Signal bei hohem k und erlauben hoehere Werte. Kann pro Kante verschieden sein: `[11.5, 11.5, 14]`. |
| `k_power` | 1, 2 oder 3 | **2** | k-Gewichtung: chi(k) wird mit k^n multipliziert. `1` betont niedrige k (leichte Streuer). `2` ist ein Kompromiss. `3` betont hohe k (schwere Streuer). |

## Mehrfachstreuung

| Parameter | Werte | Empfehlung | Erklaerung |
|-----------|-------|------------|------------|
| `N_legs` | 2 oder 4 | **4** | Max. Anzahl der Streuarme pro Pfad. `2` = nur Einfachstreuung (Absorber → Nachbar → zurueck). `4` = inkl. Dreifach- und Vierfachstreuung ueber Zwischenatome. Deutlich genauer, aber langsamere FEFF-Rechnung. |
| `R_max_for_FEFF` | Float (Angstrom) | **6** | Radius des Clusters um jeden Absorber fuer die FEFF-Rechnung. `4` = nur naechste Nachbarn (schnell). `6` = bis 3. Nachbarschale (empfohlen). `8` = sehr weit, meist unnoetig. Kann pro Kante verschieden sein. |

## Optimierungsparameter

| Parameter | Typ | Empfehlung | Erklaerung |
|-----------|-----|------------|------------|
| `Maximal_step_length` | Float (Angstrom) | **0.005** | Max. Verschiebung eines Atoms pro Iteration. Kleinere Werte = feinere Optimierung, aber langsamer. `0.001` ist sehr fein, `0.01` ist grob. |
| `Maximal_displacement` | Float (Angstrom) | **0.4** | Max. kumulative Verschiebung eines Atoms vom Startgitter. Begrenzt, wie weit sich die Struktur veraendern kann. `0.2` = konservativ (Atome bleiben nah am Gitterplatz). `0.4` = Standard. `0.8` = sehr flexibel (Vorsicht: kann zu unphysikalischen Strukturen fuehren). |
| `Initial_displacement` | Float (Angstrom) | **0.01** | Zufaellige Anfangsauslenkung aller Atome. Bricht die perfekte Gittersymmetrie, damit die Optimierung starten kann. |
| `Number_of_states` | Int | **16** | Anzahl paralleler Strukturkopien im evolutionaeren Algorithmus. Mehr States = bessere Exploration, aber mehr RAM. 16 ist Standard. |

## Konvergenz und Abbruchkriterien

Es gibt **zwei unabhaengige** Abbruchbedingungen (siehe EvAX Manual Kap. 6.11):

| Parameter | Typ | Empfehlung | Erklaerung |
|-----------|-----|------------|------------|
| `Stop_after` | Int | 100-10000 | **Harter Stop**: Run endet nach exakt dieser Anzahl Iterationen, unabhaengig von `Froze_in`. Fuer kurze Tests: 50-100. Fuer Produktion: 1000-5000. Auf `-1` setzen um nur `Froze_in` als Abbruch zu nutzen. |
| `Froze_in` | Int | 500-3000 | **Konvergenzbasierter Stop**: Ab dieser Iteration werden nur noch Verschiebungen akzeptiert, die den Fit verbessern (Simulated Annealing wird eingefroren). Dann wird alle `time2` Iterationen geprueft, ob sich der Residual ξ noch verbessert. Falls nicht → Stop. |
| `time2` | Int | 100 | Pruefintervall nach `Froze_in`: alle `time2` Iterationen wird geprueft, ob ξ sich verbessert hat. Standard: 100. |

**Wichtig**: `Stop_after` und `Froze_in` sind **unabhaengig**. Wenn `Stop_after=100` gesetzt ist, endet der Run nach 100 Iterationen — egal ob `Froze_in=150` oder `Froze_in=5000`.

Typische Konfigurationen:
- **Nur harter Stop**: `Stop_after=5000`, `Froze_in` auf sehr hohen Wert → Run endet nach genau 5000 Iterationen
- **Nur Konvergenz**: `Stop_after=-1`, `Froze_in=1500`, `time2=100` → Run laeuft bis 1500 Iterationen ohne Verbesserung (empfohlen fuer Produktion)
- **Kombination**: `Stop_after=5000`, `Froze_in=1500`, `time2=100` → Run endet bei dem Kriterium, das zuerst eintritt

## FEFF-Steuerung

| Parameter | Typ | Erklaerung |
|-----------|-----|------------|
| `FEFF_path` | Pfad | Pfad zum FEFF-Binary (wird automatisch gesetzt). |
| `Update_basis_every` | Int | Nach wie vielen Iterationen die gesamte Streubasis aktualisiert wird. Standard: automatisch. |
| `Clustering_precision` | Float | Schwelle, ab der eine veraenderte lokale Umgebung als neuer Clustertyp erkannt wird. Hoeherer Wert = seltener neue FEFF-Rechnung = schneller, aber ungenauer. |

## Kantenspezifische Parameter

| Parameter | Typ | Erklaerung |
|-----------|-----|------------|
| `Edge` | Liste | Absorberelemente, z.B. `Co Ni Cu` oder `Co Pt`. |
| `Edge_type` | Liste | Kantentyp pro Element. 3d-Metalle: `K`. 5d-Metalle (Pt, Au): `L3` (K-Kante bei >70 keV zu hoch). |
| `S02` | Float pro Kante | Amplituden-Reduktionsfaktor (0.7-1.0). Wird mit `-autoES` automatisch bestimmt. |
| `dE0` | Float pro Kante (eV) | Energieverschiebung zwischen Theorie und Experiment. Wird mit `-autoES` bestimmt. |

## Zeitabschaetzung

Die Rechenzeit haengt hauptsaechlich von der Anzahl der Atome, FEFF-Aufrufe und Iterationen ab:

| Setting | Typische Dauer pro Iteration | Erklaerung |
|---------|------------------------------|------------|
| 256 Atome, 3 Kanten, N_legs=4, R=6 | ~3-4 min | Voller Produktionslauf |
| 256 Atome, 2 Kanten, N_legs=4, R=6 | ~2 min | Weniger Kanten = schneller |
| 256 Atome, 3 Kanten, N_legs=2, R=4 | ~30 sec | Schnell aber ungenau |

**Gesamtzeit** = Anzahl Iterationen × Zeit/Iteration

Beispiel: Stop_after=5000, 128 Atome, 2 Kanten, N_legs=4, R=6 → 5000 × ~2 min = **~7 Tage**

## Empfohlene Parameter-Kombinationen

### Schneller Test (Funktionstest, ~30 min)
```yaml
Stop_after: 50
Froze_in: -1       # nur harten Stop nutzen
N_legs: 2
R_max_for_FEFF: [4, 4, 4]
```

### Parameter-Scan (~15-30 min pro Kombination)
```yaml
Stop_after: 100
Froze_in: -1        # nur harten Stop nutzen
```

### Produktion (konvergenzbasiert)
```yaml
Stop_after: -1      # kein harter Stop
Froze_in: 1500
time2: 100
N_legs: 4
R_max_for_FEFF: [6, 6, 6]
Space: w
k_power: 2
Maximal_step_length: 0.005
```

### Langzeit-Produktion (bestmoegliches Ergebnis)
```yaml
Stop_after: -1
Froze_in: 3000
time2: 100
```

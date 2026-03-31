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

## Konvergenz

| Parameter | Typ | Empfehlung | Erklaerung |
|-----------|-----|------------|------------|
| `Stop_after` | Int | 50-300 | Iterationen in der **Screening-Phase**. In dieser Phase wird die FEFF-Streubasis haeufiger aktualisiert (= genauer, aber langsamer). Fuer Parameter-Scans reichen 50. Fuer Produktion: 100-300. |
| `Froze_in` | Int | 50-500 | **Abbruchkriterium**: Der Lauf endet, wenn sich der Residual ueber diese Anzahl an Iterationen nicht mehr verbessert hat. Fuer Scans: 50. Fuer Produktion: 150-500. |

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

**Gesamtzeit** = (Stop_after + erwartete Iterationen nach Screening) × Zeit/Iteration

Beispiel: Stop_after=100, Froze_in=150 → max. ~250 Iterationen × 3 min = **~12h**

## Empfohlene Parameter-Kombinationen

### Schneller Test (zum Ausprobieren, ~1h)
```yaml
Stop_after: 30
Froze_in: 30
N_legs: 2
R_max_for_FEFF: [4, 4, 4]
```

### Parameter-Scan (~15 min pro Kombination)
```yaml
Stop_after: 50
Froze_in: 50
```

### Produktion (~7h)
```yaml
Stop_after: 100
Froze_in: 150
N_legs: 4
R_max_for_FEFF: [6, 6, 6]
Space: w
k_power: 2
Maximal_step_length: 0.005
```

### Langzeit-Produktion (~24-48h, bestmoegliches Ergebnis)
```yaml
Stop_after: 300
Froze_in: 500
```

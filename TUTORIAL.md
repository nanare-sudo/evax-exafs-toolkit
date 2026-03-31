# Tutorial: EvAX EXAFS-Analyse Schritt fuer Schritt

Dieses Tutorial beschreibt die komplette Analyse einer EXAFS-Probe mit EvAX:
- Probe konfigurieren
- Parameterscan durchfuehren
- Verschiedene Kristallstrukturen vergleichen
- Ergebnisse interpretieren

---

## Voraussetzungen

### Software installieren

```bash
# Python-Abhaengigkeiten
pip install numpy matplotlib pyyaml

# Optional (fuer chi(k)-Extraktion aus Athena-Projekten):
pip install xraylarch
```

### Binaries einrichten

EvAX und FEFF sind nicht im Repository enthalten (zu gross fuer Git).

```bash
# EvAX-Binary:
mkdir -p software
cp /pfad/zum/EvAX-6.16_linux.exe software/
chmod +x software/EvAX-6.16_linux.exe

# FEFF-Binary:
mkdir -p feff85exafs/legacy/mod/MONO
cp /pfad/zum/feff85L feff85exafs/legacy/mod/MONO/
chmod +x feff85exafs/legacy/mod/MONO/feff85L
```

**Test:**
```bash
./software/EvAX-6.16_linux.exe --help     # Sollte Usage-Info zeigen
./feff85exafs/legacy/mod/MONO/feff85L     # Sollte "No input file" melden
```

---

## Teil 1: Bestehende Probe analysieren (MEA CoNiCu)

Die MEA-Probe ist komplett vorkonfiguriert und eignet sich als Einstieg.

### 1.1 Trockenlauf — pruefen ob alles stimmt

```bash
python run.py configs/MEA_CoNiCu.yaml --dry
```

Du solltest sehen:
```
============================================================
  Probe: MEA_CoNiCu
  Medium Entropy Alloy CoNiCu — 3 Kanten gleichzeitig
  Kanten: Co (K), Ni (K), Cu (K)
============================================================

Modus: Einzellauf (Produktion)
  ...
(Trockenlauf — nichts wird ausgefuehrt)
```

Falls Fehler kommen ("Binary nicht gefunden", "chi-Datei fehlt"),
pruefe die Pfade in `configs/MEA_CoNiCu.yaml`.

### 1.2 Parameterscan starten

Der Parameterscan testet systematisch verschiedene Einstellungen und zeigt,
welche Kombination den besten Fit liefert.

```bash
# Erst anschauen was passieren wuerde:
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --dry

# Dann wirklich starten (dauert einige Stunden):
python run.py configs/MEA_CoNiCu.yaml --parameter-scan
```

**Was passiert im Hintergrund:**
- Fuer jede Parameterkombination wird ein eigenes Verzeichnis unter
  `output/MEA_CoNiCu/parameter_scan/` erstellt
- EvAX laeuft mit reduzierten Iterationen (50 statt 150+)
- Bis zu 4 Instanzen laufen parallel

**Tipp:** Fuer einen schnelleren ersten Test mit weniger Kombinationen
lassen sich in der YAML-Config die Listen kuerzen:
```yaml
parameter_scan:
  phase1:
    k_power: [2]        # Statt [1, 2, 3] — nur k² testen
    Space: [w, R]       # Statt [w, R, k] — k-Raum weglassen
```

### 1.3 Ergebnisse des Scans ansehen

```bash
# Tabelle im Terminal:
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --results

# Plots generieren:
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --plot
```

Die Tabelle ist nach Residual sortiert. **Niedrigerer Residual = besserer Fit.**
Der beste Run steht oben.

### 1.4 Strukturvergleich

Um die Kristallstruktur der Probe zu bestimmen (FCC, BCC oder HCP):

```bash
python run.py configs/MEA_CoNiCu.yaml --structure-compare
```

Das laeuft 3 EvAX-Instanzen parallel (eine pro Struktur).
Dauert ca. 5-7 Stunden.

```bash
# Ergebnis:
python run.py configs/MEA_CoNiCu.yaml --structure-compare --results
```

Die Struktur mit dem **niedrigsten Residual** passt am besten zu den Daten.

```bash
# Vergleichsplot:
python run.py configs/MEA_CoNiCu.yaml --structure-compare --plot
```

### 1.5 Ergebnisse interpretieren

Nach einem erfolgreichen Lauf befinden sich in `output/MEA_CoNiCu/` folgende Dateien:

| Datei | Was drin ist |
|-------|-------------|
| `output.dat` | Konvergenzverlauf (Iteration vs. Residual) |
| `Co/EXAFS.dat` | Berechnetes vs. experimentelles chi(k) fuer Co |
| `Co/RDF_Co_Ni.dat` | Radiale Verteilungsfunktion Co-Ni |
| `final.xyz` | Optimierte 3D-Atomstruktur |

**Physikalische Groessen, die sich ableiten lassen:**
- **Kristallstruktur?** → Struktur mit niedrigstem Residual
- **Koordinationszahlen?** → Aus den RDF-Dateien (Integration des 1. Peaks)
- **Bindungslaengen?** → Peak-Position in der RDF
- **Chemische Nahordnung?** → Warren-Cowley-Parameter aus final.xyz

---

## Teil 2: Eigene Probe einrichten

### 2.1 Config-Datei erstellen

```bash
cp configs/MEA_CoNiCu.yaml configs/MeineProbe.yaml
```

Die neue Datei oeffnen und folgende Felder anpassen:

#### Probenname und Beschreibung
```yaml
name: MeineProbe
description: "Kurze Beschreibung"
```

#### Elemente und Kanten
```yaml
edges: [Fe, Co]
edge_types: [K, K]
```

**Welcher Kantentyp?**
| Element | Kantentyp | Kantenenergie | Grund |
|---------|-----------|---------------|-------|
| Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn | **K** | 5-10 keV | 3d-Metalle, K-Kante im harten Roentgenbereich |
| Zr, Mo, Ru, Rh, Pd, Ag | **K** | 18-26 keV | 4d-Metalle |
| W, Re, Os, Ir, Pt, Au | **L3** | 10-14 keV | 5d-Metalle, K-Kante viel zu hoch (>60 keV) |
| La, Ce, Pr, Nd (Seltenerd) | **L3** | 5-7 keV | 4f-Elemente |

#### chi(k)-Dateien zuordnen
```yaml
chi_files:
  Fe: data/chi_extracted/meine_Fe_chi.dat
  Co: data/chi_extracted/meine_Co_chi.dat
```

Lege die Dateien dort ab:
```bash
cp /pfad/zu/meinen/chi-dateien/*.dat data/chi_extracted/
```

**Format:** Zweispaltig, ASCII, Leerzeichen-getrennt:
```
# k (1/Angstrom)    chi(k)
0.000000    0.123456
0.050000    0.234567
...
```

#### Startstruktur
```yaml
structure_file: structures/MeineProbe.p1
```

POSCAR-Datei erstellen oder `scripts/generate_structures.py` als Vorlage verwenden.
Die wichtigsten Punkte:
- **~256 Atome** (guter Kompromiss zwischen Genauigkeit und Rechenzeit)
- **Richtige Stoechometrie** (z.B. 128 Fe + 128 Co fuer 50/50)
- **Richtige Gitterkonstante** (aus Literatur oder XRD)

#### Parameter anpassen (optional)

Die Standardwerte funktionieren fuer die meisten Proben. Aendere nur wenn noetig:

```yaml
parameters:
  k_max: [12, 14]   # Pt hat staerkeres Signal, also hoehere k_max moeglich
```

Siehe [PARAMETER_GUIDE.md](PARAMETER_GUIDE.md) fuer Details.

### 2.2 Test und Start

```bash
# Alles korrekt konfiguriert?
python run.py configs/MeineProbe.yaml --dry

# Schneller Testlauf (30 min):
python run.py configs/MeineProbe.yaml --screening 30 --froze-in 30

# Wenn der Test funktioniert — Parameterscan:
python run.py configs/MeineProbe.yaml --parameter-scan
```

---

## Teil 3: Parameter verstehen und tunen

### Was beeinflusst was?

```
Genauigkeit ↑    |    Rechenzeit ↑
─────────────────┼──────────────────
N_legs=4         |    N_legs=4
R_max_for_FEFF=6 |    R_max_for_FEFF=6
Space=w          |    Space=w
Stop_after=300   |    Stop_after=300
Froze_in=500     |    Froze_in=500
```

Es gibt einen Trade-off: Genauere Einstellungen dauern laenger.
Der Parameterscan hilft, die optimale Kombination zu finden.

### Typische Probleme und Loesungen

**Problem: Run dauert ewig (>24h)**
```yaml
# Reduziere Iterationen:
Stop_after: 50
Froze_in: 100
# Oder reduziere FEFF-Genauigkeit:
N_legs: 2
R_max_for_FEFF: [4, 4, 4]
```

**Problem: Schlechter Fit (hoher Residual)**
- Pruefe ob die chi(k)-Daten korrekt normiert sind
- Probiere verschiedene k_min/k_max Werte
- Pruefe die Startstruktur (richtige Gitterkonstante?)
- Erhoehe Iterationen (Stop_after, Froze_in)

**Problem: "Binary nicht gefunden"**
- Pruefe die Pfade in der YAML-Config unter `software:`
- Pruefe Ausfuehrrechte: `chmod +x software/EvAX-6.16_linux.exe`

**Problem: FEFF-Fehler im Log**
- Meistens: zu viele Atome im Cluster oder ungueltige Koordinaten
- Pruefe die POSCAR-Datei (Gittervektoren muessen positiv sein)
- Reduziere R_max_for_FEFF

### Reihenfolge fuer eine neue Analyse

1. **Schnelltest** (30 min): `--screening 30 --froze-in 30`
   → Funktioniert alles? Gibt es offensichtliche Fehler?

2. **Parameterscan** (einige Stunden): `--parameter-scan`
   → Welche Parameter sind optimal?

3. **Beste Parameter in Config uebernehmen**
   → Aus der Ergebnistabelle ablesen und in die YAML schreiben

4. **Produktionslauf** (7+ Stunden): `python run.py configs/...yaml`
   → Finaler Lauf mit optimalen Parametern

5. **Strukturvergleich** (optional): `--structure-compare`
   → Nur wenn unklar ist, welche Kristallstruktur vorliegt

---

## Teil 4: Uebersicht der Befehle

| Befehl | Was es macht |
|--------|-------------|
| `python run.py CONFIG` | Einzelner Produktionslauf |
| `python run.py CONFIG --dry` | Zeigt was passieren wuerde |
| `python run.py CONFIG --parameter-scan` | Parameterstudie |
| `python run.py CONFIG --parameter-scan --dry` | Zeigt alle Kombinationen |
| `python run.py CONFIG --structure-compare` | FCC vs BCC vs HCP |
| `python run.py CONFIG --results` | Ergebnistabelle anzeigen |
| `python run.py CONFIG --parameter-scan --results` | Scan-Ergebnisse |
| `python run.py CONFIG --plot` | Plots generieren |
| `python run.py CONFIG --parallel 6` | 6 parallele Instanzen |
| `python run.py CONFIG --screening 30 --froze-in 30` | Kurzer Testlauf |

`CONFIG` ist immer der Pfad zur YAML-Datei, z.B. `configs/MEA_CoNiCu.yaml`.

---

## Haeufige Fragen

**Kann ich den Lauf abbrechen und spaeter fortsetzen?**
Ja, Ctrl+C bricht ab. Bereits fertige Teil-Runs werden beim naechsten Start
uebersprungen. Nur der gerade laufende Run geht verloren.

**Wie viel RAM braucht ein Lauf?**
Ca. 200-500 MB pro EvAX-Instanz. Bei 4 parallelen Instanzen also ~2 GB.

**Laeuft das auf einem Cluster?**
Ja. Repository dorthin kopieren, Binaries einrichten, und starten mit
`nohup python run.py ... &` oder per Batch-System (SLURM etc.).

**Was bedeutet `-autoES`?**
EvAX optimiert automatisch S0² (Amplitudenfaktor) und dE0 (Energieverschiebung)
pro Kante, bevor die eigentliche Strukturoptimierung beginnt.
Das ist fast immer die richtige Wahl.

**Was machen die Skripte in `scripts/` einzeln?**
`scripts/evax_engine.py` ist die zentrale Engine, die von `run.py` aufgerufen wird.
Die uebrigen Skripte (`parameter_scan.py`, `structure_comparison.py`, etc.) sind
eigenstaendige Implementierungen fuer spezifische Aufgaben und koennen als
Referenz fuer eigene Anpassungen dienen.

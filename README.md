# EvAX EXAFS-Analyse Toolkit

Automatisierte Multi-Edge EXAFS-Analyse mit [EvAX](http://www.baas.uni-due.de/) (Evolutionary Algorithm + Reverse Monte Carlo) und FEFF 8.5.

Entwickelt im Rahmen des Kurses Experimentelle Roentgenphysik (ERP), Universitaet Wuppertal.

---

## Schnellstart

```bash
# 1. Repo klonen
git clone https://github.com/nanare-sudo/evax-exafs-toolkit.git
cd evax-exafs-toolkit

# 2. Abhaengigkeiten installieren
pip install numpy matplotlib pyyaml

# 3. Binaries einrichten (siehe unten)

# 4. Daten in data/chi_extracted/ ablegen (siehe unten)

# 5. Analyse starten (ein Befehl!)
python run.py configs/MEA_CoNiCu.yaml
```

---

## Einrichtung

### Verzeichnisstruktur

```
ERP_Hausarbeit/
├── run.py                      # Hauptskript — einziger Einstiegspunkt
├── configs/                    # Proben-Konfigurationen (YAML)
│   ├── MEA_CoNiCu.yaml        # Beispiel: CoNiCu Medium Entropy Alloy
│   └── CoPt_nanoparticles.yaml # Beispiel: CoPt Nanopartikel
├── scripts/                    # Python-Skripte (Engine, Plots, etc.)
├── structures/                 # Startstrukturen (POSCAR-Format)
├── data/                       # Experimentelle Daten (NICHT im Git)
│   └── chi_extracted/          # chi(k)-Dateien (2 Spalten: k, chi)
├── software/                   # EvAX-Binary (NICHT im Git)
├── feff85exafs/                # FEFF-Binary (NICHT im Git)
├── output/                     # Ergebnisse (automatisch erstellt)
├── figures/                    # Plots (automatisch erstellt)
├── PARAMETER_GUIDE.md          # Alle EvAX-Parameter erklaert
└── HOWDOESITWORK.txt           # Erklaerung der Methode
```

### 1. Binaries besorgen

Die Binaries sind zu gross fuer Git und muessen separat eingerichtet werden:

```bash
# EvAX (Version 6.16):
mkdir -p software
cp /pfad/zu/EvAX-6.16_linux.exe software/
chmod +x software/EvAX-6.16_linux.exe

# FEFF 8.5:
mkdir -p feff85exafs/legacy/mod/MONO
cp /pfad/zu/feff85L feff85exafs/legacy/mod/MONO/
chmod +x feff85exafs/legacy/mod/MONO/feff85L
```

> Die Pfade muessen mit den Angaben in der YAML-Config uebereinstimmen.
> Standard: `software/EvAX-6.16_linux.exe` und `feff85exafs/legacy/mod/MONO/feff85L`

### 2. Experimentelle Daten vorbereiten

EvAX erwartet pro Absorptionskante eine **zweispaltige ASCII-Datei**:

```
# k (1/Angstrom)    chi(k) (dimensionslos)
0.00000000e+00    1.23117930e+00
5.00000000e-02    1.14055848e+00
...
```

Die Dateien werden in `data/chi_extracted/` abgelegt. Der Dateiname ist frei waehlbar,
muss aber in der YAML-Config angegeben werden.

**Woher kommen die chi(k)-Dateien?**
- Aus Athena/Demeter: Spektrum auswaehlen → "Save chi(k)" als ASCII
- Aus larch (Python): `autobk()` → `np.savetxt('chi.dat', [k, chi])`
- Aus anderen XAS-Programmen (z.B. Xanda): Export als zweispaltige Datei

> **Wichtig bei Multi-Edge-Messungen** (z.B. CoNiCu in einem Scan):
> Die Normierung muss fuer jede Kante separat gemacht werden, mit
> begrenztem Energiefenster, damit der Pre/Post-Edge der naechsten Kante
> nicht stoert. Siehe `scripts/extract_chi.py` als Beispiel.

### 3. Startstruktur erstellen

EvAX braucht ein 3D-Atommodell als Startpunkt im **POSCAR-Format** (aus VASP).
Beispiele liegen in `structures/`.

**Fuer eine neue Probe:**
1. Kristallstruktur und Gitterkonstante aus der Literatur oder XRD bestimmen
2. Superzelle mit ~256 Atomen erzeugen (z.B. 4x4x4 FCC = 256)
3. Elemente zufaellig auf Gitterplaetze verteilen
4. Als `.p1`-Datei in `structures/` ablegen

> `scripts/generate_structures.py` zeigt wie man BCC/HCP-Strukturen automatisch
> mit Python erzeugt. Fuer einfache Strukturen kann man auch VESTA oder ASE nutzen.

---

## Verwendung

Alles laeuft ueber `run.py` mit einer YAML-Konfiguration:

### Einzelner Produktionslauf

```bash
python run.py configs/MEA_CoNiCu.yaml
```

Startet EvAX mit den in der Config definierten Parametern. Ergebnisse landen in
`output/MEA_CoNiCu/production/`.

### Parameterstudie

```bash
# Zeigt erst was laufen wuerde (Trockenlauf):
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --dry

# Startet den Scan (laeuft mehrere Stunden):
python run.py configs/MEA_CoNiCu.yaml --parameter-scan

# Mit mehr parallelen Instanzen (Standard: 4):
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --parallel 6
```

Testet systematisch verschiedene Parameterkombinationen (Wavelet vs R-Raum,
k-Gewichtung, Mehrfachstreuung, etc.) und zeigt welche am besten funktionieren.

Ergebnisse in `output/MEA_CoNiCu/parameter_scan/`.

### Strukturvergleich (FCC vs BCC vs HCP)

```bash
python run.py configs/MEA_CoNiCu.yaml --structure-compare
```

Laeuft parallel fuer alle in der Config definierten Strukturen und erstellt
Vergleichsplots. Die Struktur mit dem niedrigsten Residual passt am besten.

### Nur Plots oder Ergebnisse anzeigen

```bash
# Ergebnistabelle im Terminal:
python run.py configs/MEA_CoNiCu.yaml --results
python run.py configs/MEA_CoNiCu.yaml --parameter-scan --results

# Plots generieren (ohne EvAX-Lauf):
python run.py configs/MEA_CoNiCu.yaml --plot
python run.py configs/MEA_CoNiCu.yaml --structure-compare --plot
```

### Iterationen anpassen (fuer schnelle Tests)

```bash
# Kurzer Testlauf (~30 min):
python run.py configs/MEA_CoNiCu.yaml --screening 30 --froze-in 30

# Langer Produktionslauf (~24h):
python run.py configs/MEA_CoNiCu.yaml --screening 300 --froze-in 500
```

---

## Eigene Probe einrichten

### 1. Config-Datei kopieren und anpassen

```bash
cp configs/MEA_CoNiCu.yaml configs/MeineProbe.yaml
```

Dann in der neuen Datei aendern:

```yaml
name: MeineProbe
description: "Beschreibung der Probe"

# Elemente anpassen:
edges: [Fe, Ni]           # Welche Elemente?
edge_types: [K, K]         # K fuer 3d-Metalle, L3 fuer 5d-Metalle (Pt, Au, ...)

# chi(k)-Dateien:
chi_files:
  Fe: data/chi_extracted/meine_Fe_chi.dat
  Ni: data/chi_extracted/meine_Ni_chi.dat

# Startstruktur:
structure_file: structures/MeineProbe.p1
```

### 2. Daten ablegen

```bash
cp /pfad/zu/meinen/daten/*_chi.dat data/chi_extracted/
cp /pfad/zu/meiner/struktur.p1 structures/
```

### 3. Ausfuehren

```bash
# Erst Trockenlauf um zu pruefen ob alles stimmt:
python run.py configs/MeineProbe.yaml --dry

# Dann Parameterscan:
python run.py configs/MeineProbe.yaml --parameter-scan
```

---

## Was die Parameter bedeuten (Kurzversion)

| Parameter | Was es macht | Empfehlung |
|-----------|-------------|------------|
| `Space` | Vergleichsraum (k/R/Wavelet) | `w` (Wavelet) |
| `k_min` / `k_max` | k-Bereich fuer den Fit | 3 / 11.5 |
| `k_power` | k-Gewichtung (k^n) | 2 |
| `N_legs` | Mehrfachstreuung (2=einfach, 4+=inkl. MS, Default: 8) | 4 |
| `R_max_for_FEFF` | FEFF-Clusterradius (Angstrom) | 6 |
| `Maximal_step_length` | Schrittweite pro Iteration | 0.005 |
| `Maximal_displacement` | Max. Gesamtverschiebung | 0.4 |
| `Stop_after` | Harter Stop nach N Iterationen (-1 = aus) | 5000 (Prod.) / -1 (konvergenzbasiert) |
| `Froze_in` | SA einfrieren ab Iteration N, dann Konvergenzpruefung | 1500 (Prod.) / 50 (Scan) |

**Ausfuehrliche Erklaerung**: Siehe [PARAMETER_GUIDE.md](PARAMETER_GUIDE.md)

---

## Wie EvAX funktioniert (Kurzversion)

1. **Startmodell**: 3D-Superzelle (z.B. 256 Atome auf FCC-Gitter)
2. **FEFF-Rechnung**: Berechnet theoretische EXAFS-Spektren ab initio
3. **Atom verschieben**: Zufaelliges Atom wird leicht verschoben
4. **Vergleich**: Neues Spektrum wird mit Experiment verglichen (Residual)
5. **Akzeptieren/Verwerfen**: Metropolis-Kriterium (wie Simulated Annealing)
6. **Evolutionaere Selektion**: 32 parallele Strukturen konkurrieren (Default)
7. **Konvergenz**: Wenn sich nichts mehr verbessert → fertig

Ergebnis: Optimierte 3D-Struktur, aus der man Koordinationszahlen,
Bindungslaengen und chemische Nahordnung (Warren-Cowley) ablesen kann.

**Ausfuehrliche Erklaerung**: Siehe [HOWDOESITWORK.txt](HOWDOESITWORK.txt)

---

## Vorhandene Proben

### MEA CoNiCu (`configs/MEA_CoNiCu.yaml`)
- Medium Entropy Alloy, 3 Kanten (Co-K, Ni-K, Cu-K)
- Multi-Edge-Messung in einem Scan
- FCC-Startstruktur, 256 Atome
- Optimale Parameter aus Parameterscan: Wavelet, k^2, N_legs=4

### CoPt Nanopartikel (`configs/CoPt_nanoparticles.yaml`)
- 2 Kanten: Co-K + Pt-L3
- Pt verwendet die L3-Kante (11564 eV), da die K-Kante bei ~78 keV fuer die meisten Beamlines zu hoch liegt

---

## Tipps

- **Schneller Test**: `Stop_after` und `Froze_in` auf 30 setzen
  (oder `--screening 30 --froze-in 30`), um die Konfiguration zu pruefen.
- **Parameter veraendern**: Einfach in der YAML-Datei aendern und `run.py`
  neu starten. Kein Code aendern noetig.
- **Parallele Laeufe**: `--parallel N` steuert wie viele EvAX-Instanzen
  gleichzeitig laufen. Abhaengig von CPU-Kernen und RAM.
- **Ergebnisse ansehen**: `--results` zeigt eine sortierte Tabelle im Terminal.
  `--plot` generiert Vergleichsplots.
- **Lauf abbrechen**: Ctrl+C oder `kill <PID>`. Bereits fertige Teillaeufe
  bleiben erhalten (werden beim naechsten Start uebersprungen).
- **Einzelskripte**: Die Skripte in `scripts/` koennen auch unabhaengig von
  `run.py` als Referenz oder fuer eigene Anpassungen genutzt werden.

---

## Abhaengigkeiten

- Python 3.8+
- numpy
- matplotlib
- pyyaml
- EvAX 6.16 (Binary, nicht im Repo)
- FEFF 8.5 (Binary, nicht im Repo)
- Optional: larch (fuer chi(k)-Extraktion aus Athena-Projekten)

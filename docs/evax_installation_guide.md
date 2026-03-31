# EvAX + FEFF8L Installationsanleitung (Ubuntu Linux)

## Überblick

EvAX ist ein RMC/EA-basierter Code zur Analyse von EXAFS-Spektren. Er benötigt FEFF8L als externen Rechenkern für Streuphasen und -amplituden.

**Referenz:** J. Timoshenko, A. Kuzmin, J. Purans, *J. Phys.: Condens. Matter* **26** (2014) 055401

---

## 1. FEFF8L kompilieren (mit erhöhtem natx)

EvAX benötigt eine modifizierte Version von FEFF8L mit `natx = 10000` (Standard: 1000).

### 1.1 Voraussetzungen installieren

```bash
sudo apt update
sudo apt install -y gfortran git
```

### 1.2 Quellcode herunterladen

```bash
git clone https://github.com/xraypy/feff85exafs.git
cd feff85exafs/legacy/mod/MONO/
```

Die Datei `feff85L.f` ist der monolithische Fortran-Quellcode (~63.000 Zeilen).

### 1.3 natx und nattx auf 10000 erhöhen

```bash
# Backup erstellen
cp feff85L.f feff85L_original.f

# Alle natx=1000 und nattx=1000 durch 10000 ersetzen
sed -i 's/parameter (natx =1000)/parameter (natx =10000)/g' feff85L.f
sed -i 's/parameter (nattx =1000)/parameter (nattx =10000)/g' feff85L.f

# Überprüfen (sollte jeweils 0 ergeben):
grep -c "natx =1000)" feff85L.f
grep -c "nattx =1000)" feff85L.f
```

### 1.4 Kompilieren

```bash
gfortran -o feff85L feff85L.f -O2 -w -std=legacy -fallow-argument-mismatch -fd-lines-as-comments -ffixed-form
```

**Erklärung der Compiler-Flags:**
- `-O2`: Optimierung
- `-w`: Warnungen unterdrücken (der Code hat viele harmlose Warnungen)
- `-std=legacy`: Erlaubt alten Fortran-Stil
- `-fallow-argument-mismatch`: Nötig für ältere Fortran-Aufrufkonventionen
- `-fd-lines-as-comments`: Behandelt `d`-Zeilen als Kommentare (Debug-Zeilen)
- `-ffixed-form`: Fixed-Form Fortran-Format

### 1.5 Testen

```bash
./feff85L
# Erwartete Ausgabe:
#  Feff 8.50
#  # read error: cannot find file
#  # >> file name = feff.inp
```

### 1.6 Installieren

```bash
# In einen sinnvollen Pfad kopieren, z.B.:
sudo cp feff85L /usr/local/bin/feff85L
# oder im Home-Verzeichnis:
mkdir -p ~/bin
cp feff85L ~/bin/feff85L
```

**Merke dir den Pfad** – er wird in der EvAX-Eingabedatei benötigt!

---

## 2. EvAX installieren

### 2.1 Download

- **Windows:** Direkt von https://www.dragon.lv/evax/downloads.html
- **Linux:** Per E-Mail anfragen bei janis.timoshenko@gmail.com

### 2.2 Konfiguration

In der EvAX-Eingabedatei muss der Pfad zum FEFF-Programm angegeben werden, z.B.:

```
FEFF_PATH /usr/local/bin/feff85L
```

(Genaues Format siehe EvAX-Manual)

---

## 3. Zusätzliche Software

### 3.1 VESTA (3D-Strukturvisualisierung)

```bash
# Download von: https://jp-minerals.org/vesta/en/
# Oder via Flatpak/Snap, falls verfügbar
```

VESTA wird genutzt für:
- Visualisierung der .xyz Ausgabedateien von EvAX
- Erstellung der Anfangsgeometrie (.p1 Dateien als Input für EvAX)

### 3.2 GNUplot (Ergebnis-Plots)

```bash
sudo apt install -y gnuplot
```

EvAX liefert .plt-Skripte für die grundlegende Visualisierung der Ergebnisse.

### 3.3 Python (alternativ für Plots)

```bash
sudo apt install -y python3 python3-matplotlib python3-numpy
```

---

## 4. Beispiel-Daten

Auf Zenodo gibt es ein vollständiges Beispiel für eine fcc High-Entropy Alloy (AuCuNiPdPt), das als Vorlage dienen kann:

**DOI:** https://doi.org/10.5281/zenodo.15349402

```bash
wget https://zenodo.org/records/15349402/files/AuCuNiPdPt_273K_fcc_RMC_example.zip
unzip AuCuNiPdPt_273K_fcc_RMC_example.zip
```

---

## 5. Workflow-Übersicht

1. **Anfangsstruktur** definieren (z.B. mit VESTA → .p1 Datei exportieren)
2. **Experimentelle EXAFS-Daten** vorbereiten (χ(k) Spektren)
3. **EvAX-Eingabedatei** erstellen (Pfad zu FEFF, Strukturdatei, experimentelle Daten)
4. **EvAX ausführen** (RMC/EA-Simulation)
5. **Ergebnisse analysieren**:
   - ASCII-Dateien → GNUplot/Python
   - .xyz Dateien → VESTA

---

## Zusammenfassung der Befehle

```bash
# Alles auf einmal:
sudo apt update && sudo apt install -y gfortran git gnuplot

git clone https://github.com/xraypy/feff85exafs.git
cd feff85exafs/legacy/mod/MONO/

cp feff85L.f feff85L_original.f
sed -i 's/parameter (natx =1000)/parameter (natx =10000)/g' feff85L.f
sed -i 's/parameter (nattx =1000)/parameter (nattx =10000)/g' feff85L.f

gfortran -o feff85L feff85L.f -O2 -w -std=legacy -fallow-argument-mismatch -fd-lines-as-comments -ffixed-form

mkdir -p ~/bin && cp feff85L ~/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
```

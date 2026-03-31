# Hausarbeit: Experimentelle Röntgenphysik — EvAX-Analyse

## Überblick

LaTeX-Projekt für die Hausarbeit zur Vorlesung "Experimentelle Röntgenphysik".
Thema: Lokale Strukturanalyse von CoNiCu-MEA und CoPt-Nanopartikeln mittels EXAFS und EvAX.

## Kompilieren

### Overleaf
1. ZIP-Datei bei Overleaf hochladen ("New Project" → "Upload Project")
2. `main.tex` als Hauptdokument auswählen
3. Compiler auf **pdfLaTeX** und Bibliographie auf **Biber** stellen
4. Kompilieren

### Lokal
```bash
pdflatex main.tex
biber main
pdflatex main.tex
pdflatex main.tex
```

## Projektstruktur

```
EvAX_Hausarbeit/
├── main.tex              # Hauptdatei (Präambel + Kapitel-Includes)
├── references.bib        # Literaturverzeichnis (BibTeX)
├── README.md             # Diese Datei
├── chapters/
│   ├── deckblatt.tex     # Titelseite
│   ├── abstract.tex      # Zusammenfassung
│   ├── einleitung.tex    # Kap. 1: Einleitung
│   ├── theorie.tex       # Kap. 2: Theoretische Grundlagen
│   ├── evax.tex          # Kap. 3: Die Software EvAX
│   ├── methodik.tex      # Kap. 4: Methodik
│   ├── auswertung.tex    # Kap. 5: Ergebnisse
│   ├── diskussion.tex    # Kap. 6: Diskussion
│   └── fazit.tex         # Kap. 7: Fazit und Ausblick
└── figures/
    └── xas-principle.tex # TikZ-Abbildung: XAS-Messprinzip
```

## Platzhalter

Alle Stellen, die noch ergänzt werden müssen, sind mit `\TODO{...}` markiert
und erscheinen in roter Schrift im PDF. Suche nach `TODO` um alle zu finden.

### Wichtigste Ergänzungen:
- **Deckblatt**: Autorennamen, Matrikelnummern, Betreuer, Datum
- **Abbildungen**: Plots aus der EvAX-Analyse in `figures/` ablegen und
  `\includegraphics`-Befehle einkommentieren
- **Ergebnistabellen**: Residuen, Koordinationszahlen, SRO-Parameter
- **CoPt-Probe**: Messbedingungen und Ergebnisse

## Abbildungen einfügen

Die Plots aus der Analyse (PNG/PDF) in den `figures/`-Ordner kopieren:
- `fig8_parameter_scan.png`
- `fig9_parameter_heatmap.png`
- `fig10_structure_comparison.png`
- `struct_fcc_exafs_ft.png`, `struct_fcc_rdf.png`, etc.

Dann in `chapters/auswertung.tex` die `\includegraphics`-Zeilen einkommentieren.

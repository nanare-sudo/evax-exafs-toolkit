#!/usr/bin/env python3
"""
EvAX EXAFS Analysis — Hauptskript
==================================

Ein einziger Einstiegspunkt fuer alle Analysen.
Konfiguration erfolgt ueber YAML-Dateien in configs/.

Verwendung:
    python run.py configs/MEA_CoNiCu.yaml                     # Einzelner Produktionslauf
    python run.py configs/MEA_CoNiCu.yaml --parameter-scan     # Parameterstudie
    python run.py configs/MEA_CoNiCu.yaml --structure-compare  # Strukturvergleich (FCC/BCC/HCP)
    python run.py configs/MEA_CoNiCu.yaml --plot               # Nur Plots generieren
    python run.py configs/MEA_CoNiCu.yaml --results            # Nur Ergebnistabelle
    python run.py configs/MEA_CoNiCu.yaml --dry                # Zeigt was laufen wuerde

Optionen:
    --parallel N        Anzahl paralleler EvAX-Instanzen (Standard: aus Config)
    --screening N       Screening-Iterationen ueberschreiben
    --froze-in N        Froze_in ueberschreiben
"""

import argparse
import sys
from pathlib import Path

# Add scripts directory to path
sys.path.insert(0, str(Path(__file__).parent / 'scripts'))

from evax_engine import (
    load_config,
    run_single,
    run_parameter_scan,
    run_structure_comparison,
    show_results,
    generate_plots,
)


def main():
    parser = argparse.ArgumentParser(
        description='EvAX EXAFS Analysis — Konfigurationsbasierte Ausfuehrung',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Beispiele:
  python run.py configs/MEA_CoNiCu.yaml                     Einzelner Run
  python run.py configs/MEA_CoNiCu.yaml --parameter-scan    Alle Parameter testen
  python run.py configs/MEA_CoNiCu.yaml --structure-compare FCC vs BCC vs HCP
  python run.py configs/MEA_CoNiCu.yaml --plot              Nur Plots erstellen
  python run.py configs/MEA_CoNiCu.yaml --results           Nur Ergebnisse anzeigen
        """,
    )
    parser.add_argument('config', help='Pfad zur YAML-Konfigurationsdatei')
    parser.add_argument('--parameter-scan', action='store_true',
                        help='Parameterstudie durchfuehren')
    parser.add_argument('--structure-compare', action='store_true',
                        help='Strukturvergleich (FCC/BCC/HCP) durchfuehren')
    parser.add_argument('--plot', action='store_true',
                        help='Nur Plots generieren (keine EvAX-Laeufe)')
    parser.add_argument('--results', action='store_true',
                        help='Nur Ergebnistabelle anzeigen')
    parser.add_argument('--dry', action='store_true',
                        help='Nur anzeigen, was ausgefuehrt wuerde')
    parser.add_argument('--parallel', type=int, default=None,
                        help='Anzahl paralleler EvAX-Instanzen')
    parser.add_argument('--screening', type=int, default=None,
                        help='Screening-Iterationen ueberschreiben')
    parser.add_argument('--froze-in', type=int, default=None,
                        help='Froze_in ueberschreiben')

    args = parser.parse_args()

    # Load configuration
    config_path = Path(args.config)
    if not config_path.exists():
        print(f"FEHLER: Konfigurationsdatei nicht gefunden: {config_path}")
        print(f"\nVerfuegbare Konfigurationen:")
        for f in sorted(Path('configs').glob('*.yaml')):
            print(f"  {f}")
        sys.exit(1)

    cfg = load_config(config_path)

    # Apply command-line overrides
    if args.screening is not None:
        cfg['parameters']['Stop_after'] = args.screening
    if args.froze_in is not None:
        cfg['parameters']['Froze_in'] = args.froze_in

    # Validate that required files exist
    _validate_config(cfg)

    print(f"{'='*60}")
    print(f"  Probe: {cfg['name']}")
    print(f"  {cfg['description']}")
    print(f"  Kanten: {', '.join(f'{e} ({t})' for e, t in zip(cfg['edges'], cfg['edge_types']))}")
    print(f"{'='*60}\n")

    # Dispatch to the right mode
    if args.results:
        if args.parameter_scan:
            show_results(cfg, mode='parameter_scan')
        elif args.structure_compare:
            show_results(cfg, mode='structure_compare')
        else:
            show_results(cfg, mode='single')

    elif args.plot:
        if args.parameter_scan:
            generate_plots(cfg, mode='parameter_scan')
        elif args.structure_compare:
            generate_plots(cfg, mode='structure_compare')
        else:
            generate_plots(cfg, mode='single')

    elif args.parameter_scan:
        n_parallel = args.parallel or cfg.get('parameter_scan', {}).get('max_parallel', 4)
        run_parameter_scan(cfg, n_parallel=n_parallel, dry_run=args.dry)

    elif args.structure_compare:
        n_parallel = args.parallel or cfg.get('structure_compare', {}).get('max_parallel', 3)
        run_structure_comparison(cfg, n_parallel=n_parallel, dry_run=args.dry)

    else:
        run_single(cfg, dry_run=args.dry)


def _validate_config(cfg):
    """Check that required files and binaries exist."""
    base = cfg['_base_dir']
    errors = []

    # Check software
    evax = base / cfg['software']['evax']
    if not evax.exists():
        errors.append(f"EvAX-Binary nicht gefunden: {evax}")
    feff = base / cfg['software']['feff']
    if not feff.exists():
        errors.append(f"FEFF-Binary nicht gefunden: {feff}")

    # Check chi files
    for edge, path in cfg['chi_files'].items():
        full = base / path
        if not full.exists():
            errors.append(f"chi(k)-Datei fuer {edge} nicht gefunden: {full}")

    # Check structure file
    struct = base / cfg['structure_file']
    if not struct.exists():
        errors.append(f"Strukturdatei nicht gefunden: {struct}")

    if errors:
        print("FEHLER in der Konfiguration:\n")
        for e in errors:
            print(f"  - {e}")
        print(f"\nBitte die Pfade in {cfg['_config_path']} pruefen.")
        print("Siehe README.md fuer die Einrichtung.")
        sys.exit(1)


if __name__ == '__main__':
    main()

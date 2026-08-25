"""
CLI module for Enzyme Kinetics.
"""
import argparse
import csv
import json
import sys
from typing import List

from .engine import (
    michaelis_menten,
    michaelis_menten_batch,
    lineweaver_burk,
    eadie_hofstee,
    hanes_woolf,
    calculate_kcat,
    calculate_catalytic_efficiency,
    competitive_inhibition,
    uncompetitive_inhibition,
    noncompetitive_inhibition,
    mixed_inhibition,
    substrate_inhibition,
    substrate_inhibition_optimal_conc,
    calculate_ic50,
    ic50_to_ki,
    fit_michaelis_menten,
    fit_with_kcat,
    full_kinetic_analysis,
)


def _load_data(filepath: str):
    """Load substrate concentrations and velocities from CSV or JSON."""
    if filepath.endswith('.json'):
        with open(filepath, 'r') as f:
            data = json.load(f)
        if isinstance(data, dict):
            return data['concentrations'], data['velocities']
        return [d['concentration'] for d in data], [d['velocity'] for d in data]
    else:
        concentrations = []
        velocities = []
        with open(filepath, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                concentrations.append(float(row.get('concentration', row.get('s', 0))))
                velocities.append(float(row.get('velocity', row.get('v', 0))))
        return concentrations, velocities


def cmd_fit(args):
    """Fit Michaelis-Menten parameters."""
    concentrations, velocities = _load_data(args.input)
    
    if args.enzyme_conc:
        params = fit_with_kcat(concentrations, velocities, args.enzyme_conc)
    else:
        params = fit_michaelis_menten(concentrations, velocities)
    
    result = {
        'vmax': round(params.vmax, 6),
        'km': round(params.km, 6),
        'r_squared': round(params.r_squared, 6) if params.r_squared else None,
    }
    if params.kcat is not None:
        result['kcat'] = round(params.kcat, 6)
    if params.catalytic_efficiency is not None:
        result['catalytic_efficiency'] = round(params.catalytic_efficiency, 6)
    
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print("Michaelis-Menten Fit Results")
        print("=" * 40)
        print(f"  Vmax:    {params.vmax:.6f}")
        print(f"  Km:      {params.km:.6f}")
        if params.r_squared:
            print(f"  R²:      {params.r_squared:.6f}")
        if params.kcat is not None:
            print(f"  kcat:    {params.kcat:.6f}")
        if params.catalytic_efficiency is not None:
            print(f"  kcat/Km: {params.catalytic_efficiency:.6f}")
    
    return 0


def cmd_linearize(args):
    """Run all linearization methods."""
    concentrations, velocities = _load_data(args.input)
    result = full_kinetic_analysis(concentrations, velocities)
    
    output = {
        'lineweaver_burk': {
            'vmax': round(result.lineweaver_burk['vmax'], 6),
            'km': round(result.lineweaver_burk['km'], 6),
            'r_squared': round(result.lineweaver_burk['r_squared'], 6),
        },
        'eadie_hofstee': {
            'vmax': round(result.eadie_hofstee['vmax'], 6),
            'km': round(result.eadie_hofstee['km'], 6),
            'r_squared': round(result.eadie_hofstee['r_squared'], 6),
        },
        'hanes_woolf': {
            'vmax': round(result.hanes_woolf['vmax'], 6),
            'km': round(result.hanes_woolf['km'], 6),
            'r_squared': round(result.hanes_woolf['r_squared'], 6),
        },
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("Linearization Results")
        print("=" * 50)
        for method, data in output.items():
            print(f"\n{method.replace('_', ' ').title()}:")
            print(f"  Vmax: {data['vmax']:.6f}")
            print(f"  Km:   {data['km']:.6f}")
            print(f"  R²:   {data['r_squared']:.6f}")
    
    return 0


def cmd_inhibition(args):
    """Calculate inhibition kinetics."""
    s_range = [i * args.km / 10 for i in range(1, 51)]
    
    output = {'inhibitor_concentration': args.inhibitor, 'ki': args.ki}
    
    v_no_inhib = [michaelis_menten(s, args.vmax, args.km) for s in s_range]
    
    if args.type == 'competitive':
        v_inhib = [competitive_inhibition(s, args.vmax, args.km, args.inhibitor, args.ki) for s in s_range]
    elif args.type == 'uncompetitive':
        v_inhib = [uncompetitive_inhibition(s, args.vmax, args.km, args.inhibitor, args.ki) for s in s_range]
    elif args.type == 'noncompetitive':
        v_inhib = [noncompetitive_inhibition(s, args.vmax, args.km, args.inhibitor, args.ki) for s in s_range]
    elif args.type == 'mixed':
        if not args.ki_prime:
            print("Error: --ki-prime required for mixed inhibition", file=sys.stderr)
            return 1
        v_inhib = [mixed_inhibition(s, args.vmax, args.km, args.inhibitor, args.ki, args.ki_prime) for s in s_range]
    else:
        print(f"Unknown inhibition type: {args.type}", file=sys.stderr)
        return 1
    
    output['substrate_concentrations'] = [round(s, 4) for s in s_range]
    output['v_no_inhibitor'] = [round(v, 6) for v in v_no_inhib]
    output['v_with_inhibitor'] = [round(v, 6) for v in v_inhib]
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print(f"Inhibition Analysis: {args.type}")
        print(f"  [I] = {args.inhibitor}, Ki = {args.ki}")
        print(f"{'[S]':>10} {'v (no I)':>12} {'v (with I)':>12} {'% Inhib':>10}")
        print("-" * 48)
        for s, v0, vi in zip(s_range[:15], v_no_inhib[:15], v_inhib[:15]):
            pct = (1 - vi / v0) * 100 if v0 > 0 else 0
            print(f"{s:>10.3f} {v0:>12.6f} {vi:>12.6f} {pct:>9.1f}%")
    
    return 0


def cmd_ic50(args):
    """Calculate IC50 using Cheng-Prusoff equation."""
    ic50 = calculate_ic50(args.vmax, args.km, args.substrate, args.ki, args.type)
    
    output = {
        'ic50': round(ic50, 6),
        'ki': args.ki,
        'substrate_concentration': args.substrate,
        'km': args.km,
        'inhibition_type': args.type,
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print(f"IC50 Calculation (Cheng-Prusoff)")
        print(f"  Inhibition type: {args.type}")
        print(f"  Ki:              {args.ki}")
        print(f"  [S]:             {args.substrate}")
        print(f"  Km:              {args.km}")
        print(f"  IC50:            {ic50:.6f}")
    
    return 0


def cmd_substrate_inhibition(args):
    """Analyze substrate inhibition."""
    s_opt = substrate_inhibition_optimal_conc(args.km, args.ki_substrate)
    v_opt = substrate_inhibition(s_opt, args.vmax, args.km, args.ki_substrate)
    
    output = {
        'vmax': args.vmax,
        'km': args.km,
        'ki_substrate': args.ki_substrate,
        'optimal_concentration': round(s_opt, 6),
        'max_velocity': round(v_opt, 6),
    }
    
    if args.json:
        print(json.dumps(output, indent=2))
    else:
        print("Substrate Inhibition Analysis")
        print("=" * 40)
        print(f"  Vmax:           {args.vmax}")
        print(f"  Km:             {args.km}")
        print(f"  Ki (substrate): {args.ki_substrate}")
        print(f"  Optimal [S]:    {s_opt:.4f}")
        print(f"  Max velocity:   {v_opt:.6f}")
    
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='enzyme-kinetics-michaelis-menten',
        description='Enzyme Kinetics: Michaelis-Menten analysis, linearization, inhibition, IC50',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    # Fit
    p_fit = subparsers.add_parser('fit', help='Fit Michaelis-Menten parameters')
    p_fit.add_argument('-i', '--input', required=True, help='CSV or JSON data file')
    p_fit.add_argument('--enzyme-conc', type=float, help='Enzyme concentration for kcat')
    p_fit.add_argument('--json', action='store_true', help='JSON output')
    
    # Linearize
    p_lin = subparsers.add_parser('linearize', help='All linearization methods')
    p_lin.add_argument('-i', '--input', required=True, help='CSV or JSON data file')
    p_lin.add_argument('--json', action='store_true', help='JSON output')
    
    # Inhibition
    p_inhib = subparsers.add_parser('inhibition', help='Inhibition analysis')
    p_inhib.add_argument('--vmax', type=float, required=True)
    p_inhib.add_argument('--km', type=float, required=True)
    p_inhib.add_argument('--type', required=True, choices=['competitive', 'uncompetitive', 'noncompetitive', 'mixed'])
    p_inhib.add_argument('--inhibitor', type=float, required=True, help='[I] concentration')
    p_inhib.add_argument('--ki', type=float, required=True, help='Ki value')
    p_inhib.add_argument('--ki-prime', type=float, help="Ki' for mixed inhibition")
    p_inhib.add_argument('--json', action='store_true', help='JSON output')
    
    # IC50
    p_ic50 = subparsers.add_parser('ic50', help='IC50 calculation (Cheng-Prusoff)')
    p_ic50.add_argument('--vmax', type=float, required=True)
    p_ic50.add_argument('--km', type=float, required=True)
    p_ic50.add_argument('--substrate', type=float, required=True, help='[S] at IC50 measurement')
    p_ic50.add_argument('--ki', type=float, required=True)
    p_ic50.add_argument('--type', default='competitive', choices=['competitive', 'uncompetitive', 'noncompetitive'])
    p_ic50.add_argument('--json', action='store_true', help='JSON output')
    
    # Substrate inhibition
    p_sub = subparsers.add_parser('substrate-inhibition', help='Substrate inhibition analysis')
    p_sub.add_argument('--vmax', type=float, required=True)
    p_sub.add_argument('--km', type=float, required=True)
    p_sub.add_argument('--ki-substrate', type=float, required=True)
    p_sub.add_argument('--json', action='store_true', help='JSON output')
    
    args = parser.parse_args(argv)
    
    commands = {
        'fit': cmd_fit,
        'linearize': cmd_linearize,
        'inhibition': cmd_inhibition,
        'ic50': cmd_ic50,
        'substrate-inhibition': cmd_substrate_inhibition,
    }
    
    return commands[args.command](args)


if __name__ == '__main__':
    sys.exit(main())

#!/usr/bin/env python3
"""Repair ordering of already ACE/NME-capped protein PDBs.

Requires only Python 3. Optional --validate requires OpenMM.
Assumes residue numbering (including insertion codes) reflects peptide order
within each chain, including its caps. Geometry checks reject broken ordering.
Does not build caps, add atoms/hydrogens, or repair molecular geometry.
Rejects alternate locations, multiple models, uncapped fragments, duplicate
residue identities/atom names, and non-protein residues rather than guessing.
Original files are never overwritten. Output filenames are preserved.
"""
import argparse
from collections import Counter, OrderedDict
from math import dist, isfinite
from pathlib import Path
import sys
import traceback

AMINO_ACIDS = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL HID HIE HIP HSD HSE HSP ASH GLH LYN CYM CYX'.split())


def repair(text, minimum=1.0, maximum=1.8):
    residues = OrderedDict()
    atoms_by_serial = {}
    anisou = {}
    conect = []
    headers = []
    for line in text.splitlines():
        record = line[:6].strip()
        if record in {'MODEL', 'ENDMDL'}:
            raise ValueError('MODEL/ENDMDL records: supply a single-model PDB without model wrappers.')
        if record in {'ATOM', 'HETATM'}:
            if line[16:17].strip():
                raise ValueError('Alternate locations present; select one conformation first.')
            name = line[17:20].strip()
            if name not in AMINO_ACIDS | {'ACE', 'NME'}:
                raise ValueError(f'Unsupported residue {name}; input must contain only capped protein.')
            key = (line[21:22], int(line[22:26]), line[26:27])
            residue = residues.setdefault(key, {'name': name, 'atoms': OrderedDict()})
            atom = line[12:16].strip()
            if residue['name'] != name or atom in residue['atoms']:
                raise ValueError(f'Duplicate/ambiguous residue or atom at {key}, {atom}.')
            serial = int(line[6:11])
            if serial in atoms_by_serial:
                raise ValueError(f'Duplicate atom serial {serial}.')
            xyz = tuple(float(line[i:i+8]) for i in (30, 38, 46))
            if not all(isfinite(v) for v in xyz):
                raise ValueError(f'Nonfinite coordinates at atom {serial}.')
            residue['atoms'][atom] = (line, xyz)
            atoms_by_serial[serial] = line
        elif record == 'ANISOU':
            serial = int(line[6:11])
            if serial in anisou:
                raise ValueError(f'Duplicate ANISOU serial {serial}.')
            anisou[serial] = line
        elif record == 'CONECT':
            conect.append([int(line[i:i+5]) for i in range(6, len(line), 5) if line[i:i+5].strip()])
        elif record in {'TER', 'END', 'MASTER', ''}:
            pass  # Rebuild terminators; MASTER counts would become stale.
        elif record in {'SIGATM', 'SIGUIJ'}:
            raise ValueError(f'Unsupported atom-indexed record {record}.')
        else:
            headers.append(line)
    if not residues:
        raise ValueError('No atoms found.')
    for serial, line in anisou.items():
        if serial not in atoms_by_serial or line[12:27] != atoms_by_serial[serial][12:27]:
            raise ValueError(f'Unmatched ANISOU record {serial}.')

    chains = OrderedDict()
    for key in residues:
        chains.setdefault(key[0], []).append(key)
    fragments = []
    for keys in chains.values():
        active = []
        for key in sorted(keys, key=lambda k: (k[1], k[2])):
            name = residues[key]['name']
            if name == 'ACE':
                if active:
                    raise ValueError(f'ACE {key} encountered before preceding NME.')
                active = [key]
            elif not active:
                raise ValueError(f'{name} {key} has no preceding ACE in residue-number order.')
            else:
                active.append(key)
            if name == 'NME':
                if len(active) < 3:
                    raise ValueError(f'No amino acids between caps ending at {key}.')
                fragments.append(active)
                active = []
        if active:
            raise ValueError(f'Fragment beginning at {active[0]} has no final NME.')

    for fragment in fragments:
        for left, right in zip(fragment, fragment[1:]):
            try:
                distance = dist(residues[left]['atoms']['C'][1], residues[right]['atoms']['N'][1])
            except KeyError:
                raise ValueError(f'Missing backbone C or N between {left} and {right}.') from None
            if not minimum <= distance <= maximum:
                raise ValueError(f'C-N distance {distance:.3f} A between {left} and {right}; expected {minimum}-{maximum} A. Check numbering/cap placement.')

    output = list(headers)
    serial_map = {}
    serial = 1
    for fragment in fragments:
        for key in fragment:
            for line, _ in residues[key]['atoms'].values():
                old = int(line[6:11])
                serial_map[old] = serial
                output.append(line[:6] + f'{serial:5d}' + line[11:])
                if old in anisou:
                    output.append(anisou[old][:6] + f'{serial:5d}' + anisou[old][11:])
                serial += 1
        last = next(reversed(residues[fragment[-1]]['atoms'].values()))[0]
        output.append(f'TER   {serial:5d}      {last[17:20]} {last[21]}{last[22:27]}')
        serial += 1
    if serial - 1 > 99999:
        raise ValueError('Atom/TER serials exceed PDB format limit.')
    for row in conect:
        if not row or any(old not in serial_map for old in row):
            raise ValueError('CONECT references a missing atom.')
        output.append('CONECT' + ''.join(f'{serial_map[old]:5d}' for old in row))
    output.append('END')
    before = Counter(line[:6] + line[11:] for line in atoms_by_serial.values())
    after = Counter(line[:6] + line[11:] for line in output if line[:6].strip() in {'ATOM', 'HETATM'})
    # if before != after:
    #     raise ValueError('Internal error: atom preservation check failed.')
    return '\n'.join(output) + '\n', len(fragments), len(atoms_by_serial)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('input', type=Path, help='PDB file or directory')
    parser.add_argument('-o', '--output-dir', required=True, type=Path)
    parser.add_argument('--recursive', action='store_true', help='Search subdirectories and preserve relative paths')
    parser.add_argument('--validate', action='store_true', help='Require OpenMM ff19SB system creation before saving')
    parser.add_argument('--forcefield', default='amber19/protein.ff19SB.xml', help='Protein XML used with --validate')
    args = parser.parse_args()
    source, destination = args.input.resolve(), args.output_dir.resolve()
    if not source.exists():
        parser.error(f'Input does not exist: {source}')
    root = source if source.is_dir() else source.parent
    if destination == root or destination in root.parents or root in destination.parents:
        parser.error('Use a separate output directory outside the input directory tree.')
    files = ([source] if source.is_file() else sorted(p for p in (source.rglob('*') if args.recursive else source.iterdir()) if p.is_file() and p.name.lower().endswith('_h.pdb')))
    if not files:
        parser.error('No PDB files found.')
    ff = None
    if args.validate:
        try:
            from openmm.app import PDBFile, ForceField, NoCutoff
            from io import StringIO
            ff = ForceField(args.forcefield)
        except (ImportError, ValueError) as exc:
            parser.error(f'Cannot initialize OpenMM validation: {exc}')
    failed = 0
    for path in files:
        target = destination / path.relative_to(root)
        try:
            if target.exists():
                raise ValueError(f'Output exists: {target}; use a new/empty output directory.')
            fixed, count, atoms = repair(path.read_text())
            if ff is not None:
                pdb = PDBFile(StringIO(fixed))
                ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('x') as handle:
                handle.write(fixed)
            print(f'OK   {path.name}: {count} fragments, {atoms} atoms -> {target}')
        except Exception as exc:
            failed += 1
            print(f'FAIL {path}: {exc}', file=sys.stderr)
            traceback.print_exc()
    print(f'Completed: {len(files)-failed} succeeded, {failed} failed.')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())

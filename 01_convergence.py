"""
SET BEFORE RUNNING
Section 1: input_dir, output_dir, reference_cutoff, reference_mesh,
mesh_cutoff and selected_lreal. Cutoffs are in eV. number_of_atoms=702
must agree with the source structure.

Compare numerical convergence at the fixed A702 geometry.
Inputs: inputs/01_convergence/results.csv and the corresponding OUTCAR.gz files.
Outputs: results/01_convergence/cutoff_comparison.tsv and mesh_comparison.tsv.
No earlier analysis step is required.

A finite plane-wave basis and k mesh approximate the electronic problem.
Total energies are therefore compared while changing one numerical setting
at a time. The CSV contains settings, completion status, free energies,
energies per atom and maximum forces. Each status is checked against the
electronic convergence message, force block and normal completion in OUTCAR.
The recorded total energy and energy per atom are also checked. A missing
output may belong to an incomplete row, but supplies no energy observation.

The comparison uses the 12 completed LREAL=Auto cases. The cutoff series
keeps the mesh fixed and uses 800 eV as its reference. The mesh series keeps
the cutoff at 500 eV and uses 3x1x1 as its reference. Subtracting the reference
and multiplying eV/atom by 1000 gives meV/atom. Both tables retain the recorded
maximum force; the energy differences also appear in the terminal.

The compared quantity is VASP's free energy TOTEN. The extrapolated
zero-smearing energy is a separate quantity. Atom count, geometry and
potentials must agree across the series. Small energy differences alone
do not establish converged forces or band edges; this unrelaxed structure
can still have large forces. No completed LREAL=False comparison is available.
"""

from pathlib import Path
import csv
import gzip
import re

# 1. Paths and parameters
# ------------------------------------------------------------
package_dir = Path(__file__).resolve().parents[1]
input_dir = package_dir / 'inputs/01_convergence'
output_dir = package_dir / 'results/01_convergence'
reference_cutoff = 800
reference_mesh = '3x1x1'
mesh_cutoff = 500
selected_lreal = 'AUTO'
number_of_atoms = 702

# 2. Read the recorded table and check completion against OUTCAR
# ------------------------------------------------------------
with (input_dir / 'results.csv').open() as file:
    all_rows = list(csv.DictReader(file))
rows = []
for row in all_rows:
    relative = row['folder'].removeprefix('calculations/')
    output_path = input_dir / relative / 'OUTCAR.gz'
    text = ''
    if output_path.exists():
        with gzip.open(output_path, 'rt') as file:
            text = file.read()
    # All three pieces of evidence are required; a filename or CSV label alone
    # cannot establish that a static energy was computed successfully.
    complete = ('General timing and accounting' in text
                and 'aborting loop because EDIFF is reached' in text
                and 'TOTAL-FORCE (eV/Angst)' in text)
    assert complete == (row['status'] == 'completed'), relative
    if complete:
        # The final TOTEN is the free energy for this calculation. Checking the
        # stored per-atom value also catches use of a different atom count.
        raw_energy = float(re.findall(r'free  energy\s+TOTEN\s*=\s*([-\d.]+)', text)[-1])
        assert abs(raw_energy - float(row['energy'])) < 1e-7
        assert abs(raw_energy / number_of_atoms - float(row['energy_per_atom'])) < 1e-9
        if row['lreal'] == selected_lreal:
            rows.append(row)
assert len(rows) == 12

# 3. Compare cutoffs at a fixed k mesh
# ------------------------------------------------------------
# Keep the k mesh fixed so this comparison isolates the plane-wave cutoff.
cutoff_rows = [row for row in rows if row['kpoints'] == reference_mesh]
reference = next(row for row in cutoff_rows if int(row['encut']) == reference_cutoff)
reference_energy = float(reference['energy_per_atom'])
output_dir.mkdir(parents=True, exist_ok=True)
with (output_dir / 'cutoff_comparison.tsv').open('w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['ENCUT_eV', 'free_energy_difference_meV_per_atom', 'Fmax_eV_per_A'])
    for row in sorted(cutoff_rows, key=lambda row: int(row['encut'])):
        difference = 1000 * (float(row['energy_per_atom']) - reference_energy)
        writer.writerow([row['encut'], f'{difference:.6f}', row['max_force']])
        print(f"ENCUT {row['encut']} eV: {difference:.6f} meV/atom relative to {reference_cutoff} eV")

# 4. Compare k meshes at a fixed cutoff
# ------------------------------------------------------------
# Hold ENCUT fixed and use the same 3x1x1 reference mesh.
mesh_rows = [row for row in rows if int(row['encut']) == mesh_cutoff]
reference = next(row for row in mesh_rows if row['kpoints'] == reference_mesh)
with (output_dir / 'mesh_comparison.tsv').open('w') as file:
    writer = csv.writer(file, delimiter='\t')
    writer.writerow(['mesh', 'free_energy_difference_meV_per_atom', 'Fmax_eV_per_A'])
    for row in mesh_rows:
        difference = 1000 * (float(row['energy_per_atom']) - float(reference['energy_per_atom']))
        writer.writerow([row['kpoints'], f'{difference:.6f}', row['max_force']])
        print(f"Mesh {row['kpoints']}: {difference:.6f} meV/atom relative to {reference_mesh}")
print('The 500 eV choice remains a working compromise; LREAL=False has no completed comparison.')

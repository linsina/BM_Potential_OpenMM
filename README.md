# BM_Potential_OpenMM
This repository provides a utility function for applying a Born–Mayer potential to nonbonded interactions involving selected ions in OpenMM, while all other interactions are modeled using standard Lennard‑Jones 12‑6 or 12‑6‑4 potentials.

The implementation allows hybrid nonbonded models where specific ion interactions use the Born–Mayer form, while the rest of the system continues to use conventional Lennard‑Jones potentials.

# Usage
After creating the OpenMM System object, call the function add_126_4_or_126_4_bm_opt in your simulation script.

```python
nonbMethod = 'PME'
nonbCutoff = 1.0*nanometer
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbCutoff, constraints=HBonds)
add_126_4_or_126_4_bm_opt(top_file, system, nonbMethod, nonbCutoff, use_ig=True, lj_1264=False, lj_1264_bm_opt=True)
```

top_file must be an AMBER .prmtop file that includes additional sections containing Born–Mayer parameters for the ions. The expected format of these additional sections is defined in the script. Please inspect the implementation to see the exact structure required for the topology file.

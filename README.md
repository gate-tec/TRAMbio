![PyPI - Python Version](https://img.shields.io/pypi/pyversions/TRAMbio)
![GitHub](https://img.shields.io/github/license/gate-tec/trambio)
[![Build & PyTests](https://github.com/gate-tec/TRAMbio/actions/workflows/tests.yml/badge.svg?branch=develop)](https://github.com/gate-tec/TRAMbio/actions/workflows/tests.yml)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](.github/CODE_OF_CONDUCT.md)

# TRAMbio

*Topological Rigidity Analysis in Molecular Biology* (TRAMbio)
is a package based on and centered on the pebble game. It provides
functionality for general applications in graph theory including testing
for $(k,l)$-sparsity as well as determining the $(k,l)$-spanning
subgraphs, i.e., the $(k,l)$-rigid components. With regard to
molecular data, in particular proteins, TRAMbio provides tools for the
rigidity analysis on an atom or residue-level with further functionality
towards specialized tasks like (a) simulated protein unfolding ([Rader et al. 2002])
on single-state protein data and (b) time-based tracking of
rigid component changes in molecular dynamics (MD) trajectories.

[Rader et al. 2002]: https://doi.org/10.1073/pnas.062492699

- Source: https://github.com/gate-tec/TRAMbio
- Bug reports: https://github.com/gate-tec/TRAMbio/issues

<details><summary>Table of contents</summary>

- [Installation](#installation)
- [Examples](#examples)
- [CLI commands](#cli)
- [Envoronment Variables](#environment-variables)
- [Data Formats](#data-formats)
  - [Component-XML Files](#component-xml-files)
  - [Edge Data Files](#edge-data-files)
  - [PyMol Command Files](#pymol-command-files)
  - [Residue-level Component-JSON Files](#residue-level-component-json-files)
- [Licenses](#licenses)
- [Citation](#citation)
</details>

## Installation

### System requirements
- python (>=3.8)

> [!NOTE]
> If you intend to use TRAMbio for analysis of protein trajectories, please note that building Cython wheels for
> [MDAnalysis](https://pypi.org/project/MDAnalysis/) for Python 3.8 on macOS may fail.
> In that case, users are encouraged to use [mdtraj](https://pypi.org/project/mdtraj/) as an alternative.

### Install from PyPI
The default installation only carries dependencies for the [`tram-pebble`](#cli) command as well as the API of the general (k,l)-Component Pebble Game:
```bash
pip install TRAMbio
```
Additional installation profiles are as follows:
- `molecule` is required for calculations involving protein files or data in [PDB v.3 format](https://www.wwpdb.org/documentation/file-format).
- `trajectory` extends the previous profile with dependencies required for analysis of [MD Simulation trajectories](#cli). As an alternative, you can manually install [biopandas](https://pypi.org/project/biopandas) and either [MDAnalysis](https://pypi.org/project/MDAnalysis/) or [mdtraj](https://pypi.org/project/mdtraj/).
- `xml` enables utility for verifying TRAM's [XML output](#component-xml-files). As this depends on the non-python [libxml2](https://github.com/GNOME/libxml2) library, the dependency is made optional.
- `all` combines all the above profiles.

An example installation command for protein analysis and XML verification but without support for trajectories would be:
```bash
pip install TRAMbio[molecule,xml]
```

Alternatively, you can manually download TRAMbio from [GitHub](https://github.com/gate-tec/TRAMbio/releases).
To install one of these versions, unpack it and run the following from the top-level source directory using the Terminal:
```bash
pip install .[all]
```

## Examples

Execution of the general (k,l)-Pebble Game via API:
```pycon
>>> import networkx as nx
>>> from TRAMbio.pebble_game.pebble_game import run_pebble_game
>>> graph = nx.complete_graph(5)
>>> rho, _ , graph_type = run_pebble_game(graph, k=2, l=3)
>>> print(f"The graph is {graph_type} with {rho} redundant edges.")
'The graph is over-constrained with 3 redundant edges.'
```

Command-line execution of protein analysis:
```bash
# Fetch and analyse
tram-pdb --pdb 1l2y
# Convert to PyMol script
tram-pymol --pdb 1l2y --xml 1l2y_components.xml
# Visualize
pymol 1l2y_components.pml
```

## CLI

- `tram-pdb`: command for analyzing a single-state protein from a PDB file.

  The notable arguments are:
  - `-p`, `--pdb`: Protein input file in PDB v3 format.
  - `-o`, `--out-dir`: Directory for output files. (default: next to input file)
  - `-n`, `--name`: Alternate name for protein (used as output prefix).  
  If not specified, name is derived from input file name.
- `tram-xtc`: command for analyzing an MD Simulation trajectory.

  The notable arguments are:
  - `-x`, `--xtc`: Trajectory file in XTC format.
  - `-p`, `--pdb`: Protein input file in PDB v3 format.
  - `-o`, `--out-dir`: Directory for output files. (default: next to input file)
  - `-n`, `--name`: Alternate name for protein (used as output prefix).  
  If not specified, name is derived from input file name.
  - `-s`, `--stride`: Only processes every stride-th frame. (default: 50)  
  Negative values result in only the first frame being processed.
  - `-m`, `--module`: Base module for trajectory loading.  
  Requires either module [MDAnalysis](https://pypi.org/project/MDAnalysis/) (current default) or [mdtraj](https://pypi.org/project/mdtraj/) to load trajectories.
- `tram-pymol`: command for creating a PyMol visualization from the output of above commands.

  The notable arguments are:
  - `-p`, `--pdb`: Protein input file in PDB v3 format.
  - `-x`, `--xml`: Path to components XML resulting from the component analysis.
  - `--xtc`: Trajectory file in XTC format. Required for visualizing results from `tram-trajectory`.
  - `-o`, `--out-dir`: Directory for output files. (default: next to input file)
  - `-n`, `--name`: Alternate name for protein (used as output prefix).  
  If not specified, name is derived from input file name.
  - `-b`, `--bnd-file`: Optionally, the `.bnd` file resulting from the component analysis can be provided, if hydrogen bonds should be included in the visualization. (Not recommended for trajectories)
  - `-m`, `--module`: Base module for trajectory loading. Only used when `--xtc` is present. (default: MDAnalysis)  
  Requires either module [MDAnalysis](https://pypi.org/project/MDAnalysis/) (current default) or [mdtraj](https://pypi.org/project/mdtraj/) to load trajectories.
- `tram-residue`: command for converting (atom-level) component results to residue-level.

  The notable arguments are:
  - `-x`, `--xml`: Path to components XML resulting from the component analysis.
  - `-o`, `--out-dir`: Directory for output files. (default: next to input file)
  - `-n`, `--name`: Alternate name for protein (used as output prefix).  
  If not specified, name is derived from input file name.
- `tram-pebble`: command for applying the general (k,l)-Component Pebble Game to provided graphs.

  The notable arguments are:
  - `-g`, `--graph`: Path to the input file in [GRAPHML format](http://graphml.graphdrawing.org/).
  - `-k`, `--k-param`: Parameter `k` of the Pebble Game. Needs to be a positive integer. (default: 2)
  - `-l`, `--l-param`: Parameter `l` of the Pebble Game. Nedds to be in interval `[0,2k)` (default: 3)
  - `-o`, `--out-dir`: Directory for output files. (default: next to input file)
  - `-n`, `--name`: Alternate name for protein (used as output prefix).  
  If not specified, name is derived from input file name.

## Environment Variables

Multiple features of the CLI commands can be customized via environment variables.

<details><summary>Environment variables for controlling atomic interactions:</summary>

<table>
<thead>
<tr><th>Name</th><th>Description</th><th>Data Type</th><th>Default</th></tr>
</thead>
<tbody>
<tr><td><code>TRAM_HYDROGEN_INCLUDE</code></td><td>Whether to insert hydrogen bonds</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_HYDROGEN_MINIMAL_LENGTH</code></td><td>Minimum distance between atoms to consider for a hydrogen bond (in Angstroms)</td><td><code>float</code></td><td><code>2.6</code></td></tr>
<tr><td><code>TRAM_HYDROGEN_ENERGY_THRESHOLD</code></td><td>Energy threshold for inclusion of hydrogen bonds. All bonds with energy lower or equal to this threshold are included</td><td><code>float</code></td><td><code>-0.1</code></td></tr>
<tr><td><code>TRAM_HYDROGEN_CUTOFF_DISTANCE</code></td><td>Maximum distance between atoms to consider for a hydrogen bond (in Angstroms)</td><td><code>float</code></td><td><code>3.0</code></td></tr>
<tr><td><code>TRAM_HYDROGEN_STRONG_ENERGY_THRESHOLD</code></td><td>Separate energy threshold for 'strong' hydrogen bonds. Every bond with energies greater than this will be considered 'weak' and modeled with a gradually fewer number of bars</td><td><code>float</code></td><td><code>0.0</code></td></tr>
<tr><td><code>TRAM_HYDROGEN_BAR_COUNT</code></td><td>Number of bars modeling a hydrogen bond in the pebble game</td><td><code>integer</code></td><td><code>5</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_INCLUDE</code></td><td>Whether to insert hydrophobic interactions</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_MINIMAL_LENGTH</code></td><td>Whether to only keep the hydrophobic interaction with the shortest length for each atom</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_USE_POTENTIAL</code></td><td>Whether to use hydrophobic potential calculation instead of surface calculation</td><td><code>bool</code></td><td><code>false</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_SURFACE_CUTOFF_DISTANCE</code></td><td>Maximum distance between the VDW surfaces of atoms in hydrophobic interactions (Surface only, in Angstroms)</td><td><code>float</code></td><td><code>0.25</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_POTENTIAL_CUTOFF_DISTANCE</code></td><td>Maximum distance between atoms to consider for hydrophobic interactions (Potential only, in Angstroms)</td><td><code>float</code></td><td><code>9.0</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_SCALE_14</code></td><td>Scaling factor for energy potential in hydrophobic interaction between 3rd-degree (1-4) covalent neighbors (Potential only)</td><td><code>float</code></td><td><code>0.5</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_SCALE_15</code></td><td>Scaling factor for energy potential in hydrophobic interaction between 4th-degree (1-5) covalent neighbors (Potential only)</td><td><code>float</code></td><td><code>1.0</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_SCALE_UNBOUNDED</code></td><td>Scaling factor for energy potential in hydrophobic interaction between 5th or higher degree covalent neighbors (Potential only)</td><td><code>float</code></td><td><code>0.5</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_ENERGY_THRESHOLD</code></td><td>Energy threshold for inclusion of hydrophobic interactions. All interactions with energy lower or equal to this threshold are included (Potential only)</td><td><code>float</code></td><td><code>-0.1</code></td></tr>
<tr><td><code>TRAM_HYDROPHOBIC_BAR_COUNT</code></td><td>Number of bars modeling a hydrophobic interaction in the pebble game</td><td><code>integer</code></td><td><code>3</code></td></tr>
<tr><td><code>TRAM_DISULPHIDE_INCLUDE</code></td><td>Whether to insert disulphide bridges</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_DISULPHIDE_CUTOFF_DISTANCE</code></td><td>Maximum distance between sulphur atoms to consider for a disulphide bridge(in Angstroms)</td><td><code>float</code></td><td><code>3.0</code></td></tr>
<tr><td><code>TRAM_CATION_PI_INCLUDE</code></td><td>Whether to insert cation-pi interactions</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_CATION_PI_CUTOFF_DISTANCE</code></td><td>Maximum distance between atoms to consider for a cation-pi interaction (in Angstroms)</td><td><code>float</code></td><td><code>3.0</code></td></tr>
<tr><td><code>TRAM_CATION_PI_BAR_COUNT</code></td><td>Number of bars modeling a cation-pi interaction in the pebble game</td><td><code>integer</code></td><td><code>3</code></td></tr>
<tr><td><code>TRAM_AROMATIC_INCLUDE</code></td><td>Whether to insert aromatic interactions</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_AROMATIC_CUTOFF_DISTANCE_PI</code></td><td>Maximum distance between aromatic centers to consider for a pi-stacking (in Angstroms)</td><td><code>float</code></td><td><code>7.0</code></td></tr>
<tr><td><code>TRAM_AROMATIC_CUTOFF_DISTANCE_T</code></td><td>Maximum distance between aromatic centers to consider for a t-stacking (in Angstroms)</td><td><code>float</code></td><td><code>5.0</code></td></tr>
<tr><td><code>TRAM_AROMATIC_ANGLE_VARIANCE</code></td><td>Maximum allowed variance of the interaction angle in aromatic interactions (in degrees)</td><td><code>float</code></td><td><code>5.0</code></td></tr>
<tr><td><code>TRAM_AROMATIC_BAR_COUNT</code></td><td>Number of bars modeling an aromatic interaction in the pebble game</td><td><code>integer</code></td><td><code>3</code></td></tr>
<tr><td><code>TRAM_SSBOND_INCLUDE</code></td><td>Whether to insert SSBOND records from PDB data</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_LINK_INCLUDE</code></td><td>Whether to insert LINK records from PDB data</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_CONECT_INCLUDE</code></td><td>Whether to insert CONECT records from PDB data</td><td><code>bool</code></td><td><code>true</code></td></tr>
</tbody>
</table>
</details>

<details><summary>Environment variables for loading PDB data:</summary>

<table>
<thead>
<tr><th>Name</th><th>Description</th><th>Data Type</th><th>Default</th></tr>
</thead>
<tbody>
<tr><td><code>TRAM_PDB_UNIQUE_BONDS</code></td><td>Whether to limit annotations for atomic edges to a single, unique label</td><td><code>bool</code></td><td><code>false</code></td></tr>
<tr><td><code>TRAM_PDB_KEEP_HETS</code></td><td>Whether to include <code>HETATM</code> records from PDB data</td><td><code>bool</code></td><td><code>true</code></td></tr>
</tbody>
</table>
</details>

<details><summary>Environment variables for controlling the workflows:</summary>

<table>
<thead>
<tr><th>Name</th><th>Description</th><th>Data Type</th><th>Default</th></tr>
</thead>
<tbody>
<tr><td><code>TRAM_VERBOSE</code></td><td>Verbosity switch for various outputs like progress bars</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_PEBBLE_GAME_K</code></td><td>Parameter k for the pebble game (Overwritten from command line)</td><td><code>integer</code></td><td><code>2</code></td></tr>
<tr><td><code>TRAM_PEBBLE_GAME_L</code></td><td>Parameter l for the pebble game (Overwritten from command line)</td><td><code>integer</code></td><td><code>3</code></td></tr>
<tr><td><code>TRAM_PEBBLE_GAME_THREADS</code></td><td>Number of threads for multiprocessing. A value of 1 specifies no multiprocessing (Overwritten from command line)</td><td><code>integer</code></td><td><code>1</code></td></tr>
<tr><td><code>TRAM_XTC_MODULE</code></td><td>The selected third-party module for loading XTC files (Overwritten from command line)</td><td><code>MDAnalysis</code> or <code>mdtraj</code></td><td><code>MDAnalysis</code></td></tr>
<tr><td><code>TRAM_XTC_STRIDE</code></td><td>Selected stride for XTC frames (Overwritten from command line)</td><td><code>integer</code></td><td><code>50</code></td></tr>
<tr><td><code>TRAM_XTC_DYNAMIC_SCALING</code></td><td>Whether to activate dynamic core allocation during multiprocessing setup. If False, 2/3 of the available cores will be used for graph construction and 1/3 for the pebble game runs</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_RESIDUE_MIN_KEY</code></td><td>Minimum value for state keys. Either a float for single frame PDBs, indicating the minimum strength for present hydrogen bonds, or a starting frame number for trajectories (Overwritten from command line)</td><td><code>string</code></td><td>No limit</td></tr>
<tr><td><code>TRAM_RESIDUE_MAX_STATES</code></td><td>Maximum number of states to visualize. Values lower than 1 indicate no limit. (Overwritten from command line)</td><td><code>integer</code></td><td><code>0</code></td></tr>
<tr><td><code>TRAM_RESIDUE_THRESHOLD</code></td><td>Percentage of required residue atoms within a component to count it as present on the residue level. Requires PDB data to be provided</td><td><code>float</code></td><td><code>0.8</code></td></tr>
<tr><td><code>TRAM_RESIDUE_USE_MAIN_CHAIN</code></td><td>Whether to consider a residue present in the residue-level component if the N-CA-C main chain atoms are contained</td><td><code>bool</code></td><td><code>true</code></td></tr>
<tr><td><code>TRAM_PYMOL_ALL_WEIGHTED_BONDS</code></td><td>Whether to include all weighted interactions (from a provided bond file) in the PyMol frames or only hydrogen bonds and salt-bridges (PDB to PyMol only)</td><td><code>bool</code></td><td><code>false</code></td></tr>
</tbody>
</table>
</details>

## Data Formats

### Component-XML Files

The structure analysis results are written as an XML file with a root `<graph>` tag, containing:
- `<components>` list of "atomic" components, i.e. components that are never subdivided, but instead fully present (possibly joined together).
  `size` attribute describes the number of base components.
  - `<omponent>` individual base component with a unique `id` and a denoted `size` (including halo). Some components have an optional `structure` attribute, if they form a known configuration.
    - `<nodes>` list of nodes, that form this component
      - `<node>` tag, containing the unique node-id
    - `<halo>` list of nodes that should be considered part of this component (cf. hinge properties in mechanical model)
      - `<node>` tag, containing the unique node-id
- `<states>` list of all states, being either configurations (dilution analysis) or individual frames (trajectory analysis)
  - `<state>` list of individual components for the specific `key` attribute, being either a minimal energy threshold (dilution) or the frame number (trajectory)
    - `<component>` individual component for this state
      - `<halo>` [see above]
      - `<components>` list of atomic components (sub-components) that build up this component
        - `<component/>` tag, referencing the `id` of the respective atomic component

### Edge Data Files

Data on interactions collected during analysis is written to a tab-separated `.bnd` bond-list file with the columns:
- `node1`: first node
- `node2`: second node
- `type`: the main interaction type between these nodes
- `type_set`: the full set of possible interactions (if not filtered)
- `key`: currently only bonding energy for hydrogen bonds
- `extra`: list of additional information for certain interaction types
  - hydrogen bonds: up to three values for known angles in order of theta (donor-hydrogen-acceptor), phi (hydrogen-acceptor-base), and gamma (angle between normals of the two sp²-planes)
  - cation-pi interactions: interaction angle between cation and aromatic normal vector
  - aromatic interactions: interaction angle between the two aromatic normal vectors followed by the two angles between the aromatic distance vector and the respective aromatic normal vector

Additionally, only for processed trajectories, the file documents interaction changes between frames, indicated by the following columns:
- `frame`: the specific frame where this interaction changed
- `_merge`: indicating the type of change,
  that this interaction was either added (`right_only`) or removed (`left_only`) in this frame

### PyMol Command Files

The `tram-pymol` command generates a PDB and a `.pml` file. The latter can be directly run in a PyMol instance from its directory.

### Residue-level Component-JSON Files

The residue-level component structure resulting from `tram-residue` is stored in JSON file with the following structure:
```json
{
  "KEY-1": [
    [
      "AMINO-ACID-1",
      "AMINO-ACID-2"
    ],
    [
      "AMINO-ACID-X",
      "AMINO-ACID-Y",
      ...
    ],
    ...
  ],
  "KEY-2": [...]
}
```
The `KEY-X` entries are the keys of the corresponding `<state>` in the components XML. Each of them maps to a list of residue-level components.

## Licenses

- [NetworkX](https://github.com/networkx/networkx) is Revised BSD licensed.
- [BioPandas](https://github.com/BioPandas/biopandas) is Revised BSD licensed.
- TRAMbio itself is MIT licensed.

## Citation

If you use TRAMbio as part of your workflow in a scientific publication, please consider citing the TRAMbio repository as follows:

```bibtex
@article{handke2025trambio,
  title =     {\texttt{TRAMbio}: A Flexible Python Package for Graph Rigidity Analysis of Macromolecules},
  author =    {Handke, Nicolas and Gatter, Thomas and Reinhardt, Franziska and Stadler, Peter F.},
  year =      {2025},
  journal =   {BMC Bioinformatics},
  version =   {26},
  number =    {266},
  doi =       {10.1186/s12859-025-06300-3}
}
```

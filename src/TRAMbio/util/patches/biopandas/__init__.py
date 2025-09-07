"""
PandasPdb handles `charge` as float and since the column entry is empty for the most cases,
i.e., no charge (given), the **entire** column gets replaced by NaN values.
This patch changes the handling of `charge` entries to the exact format of PDB v3.3:

"Columns 79 - 80 indicate any charge on the atom, e.g., 2+, 1-. In most cases, these are blank."
"""
import re
import biopandas.pdb.engines as pdb_engines

__all__ = []

charge_handler = {
    "id": "charge",
    "line": [78, 80],
    "type": str,
    "strf": lambda x: (
        str(int(re.sub(r"[+-]", "", x)))[-1] + ("-" if "-" in x else "+") if len(x.strip()) > 0 else ""
    ),
}

pdb_engines.pdb_atomdict[-1] = charge_handler
pdb_engines.pdb_anisoudict[-1] = charge_handler

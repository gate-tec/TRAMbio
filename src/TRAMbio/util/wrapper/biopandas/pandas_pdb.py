import numpy as np
import pandas as pd
import re
import warnings
from TRAMbio.util.wrapper.pandas.warning import WarningWrapper
from looseversion import LooseVersion

from biopandas.pdb.engines import pdb_df_columns, pdb_records

from biopandas.pdb import PandasPdb
from io import StringIO

import sys
from urllib.request import urlopen

pd_version = LooseVersion(pd.__version__)


class CustomPandasPdb(PandasPdb):

    def get_model_start_end(self) -> pd.DataFrame:
        """Get the start and end of the models contained in the PDB file.

        Extracts model start and end line indexes based
          on lines labelled 'OTHERS' during parsing.

        Returns
        ---------
        pandas.DataFrame : Pandas DataFrame object containing
          the start and end line indexes of the models.
        """

        other_records = self.df["OTHERS"]

        idxs = other_records.loc[other_records["record_name"] == "MODEL", :].copy()
        ends = other_records.loc[other_records["record_name"] == "ENDMDL", :].copy()
        idxs.columns = ["record_name", "model_idx", "start_idx"]
        idxs.loc[:, "end_idx"] = ends.line_idx.values
        # Actual change: reset index after copying from OTHERS df
        idxs.reset_index(drop=True, inplace=True)
        # If structure only contains 1 model, create a dummy df mapping all lines to model_idx 1
        if len(idxs) == 0:
            n_lines = len(self.pdb_text.splitlines())
            idxs = pd.DataFrame(
                [
                    {
                        "record_name": "MODEL",
                        "model_idx": 1,
                        "start_idx": 0,
                        "end_idx": n_lines,
                    }
                ]
            )

        return idxs

    # Manual fix for incompatibility towards OTHERS entries in current function "to_pdb_stream"
    # Adapts function PandasPdb.to_pdb from biopandas.pdb.pandas_pd which handles OTHERS entries correctly
    def to_pdb_stream(self, records=None, append_newline=True) -> "StringIO":
        """Write record DataFrames to a PDB-stream.

        Parameters
        ----------
        records : iterable, default: None
            A list of PDB record sections in
            {'ATOM', 'HETATM', 'ANISOU', 'OTHERS'} that are to be written.
            Writes all lines to PDB-stream if `records=None`.

        append_newline : bool, default: True
            Appends a new line at the end of the PDB-stream if True

        """
        if not records:
            records = self.df.keys()

        dfs = {r: self.df[r].copy() for r in records if not self.df[r].empty}

        for r in dfs:
            column_order = []
            for col in pdb_records[r]:
                column_order.append(col["id"])
                if col["id"] not in dfs[r].columns:
                    default = 0 if col["type"] in [int, float] else " "
                    dfs[r][col["id"]] = pd.Series(default, index=dfs[r].index).apply(col["strf"])
                else:
                    dfs[r][col["id"]] = dfs[r][col["id"]].apply(col["strf"])
            column_order.append("line_idx")
            dfs[r] = dfs[r][column_order]
            dfs[r]["OUT"] = pd.Series("", index=dfs[r].index)

            for c in dfs[r].columns:
                # fix issue where coordinates with four or more digits would
                # cause issues because the columns become too wide
                if c in {"x_coord", "y_coord", "z_coord"}:
                    for idx in range(dfs[r][c].values.shape[0]):
                        if len(dfs[r][c].values[idx]) > 8:
                            dfs[r][c].values[idx] = str(dfs[r][c].values[idx]).strip()
                if c in {"line_idx", "OUT"}:
                    pass
                elif r in {"ATOM", "HETATM"} and c not in pdb_df_columns:
                    pass
                else:
                    dfs[r]["OUT"] = dfs[r]["OUT"] + dfs[r][c]

        if pd_version < LooseVersion("0.23.0"):
            df = pd.concat(dfs)
        else:
            df = pd.concat(dfs, sort=False)

        df.sort_values(by="line_idx", inplace=True)

        output = StringIO()

        s = df["OUT"].tolist()
        for idx in range(len(s)):
            if len(s[idx]) < 80:
                s[idx] = f"{s[idx]}{' ' * (80 - len(s[idx]))}"
        to_write = "\n".join(s)
        output.write(to_write)
        if append_newline:
            output.write("\n")

        output.seek(0)
        return output

    def export_first_model(self):
        with WarningWrapper(FutureWarning):
            model_start_end = self.get_model_start_end()
            first_model = int(model_start_end.loc[0, 'model_idx'])
            first_model_start_idx = int(model_start_end.loc[0, 'start_idx'])
            first_model_end_idx = int(model_start_end.loc[0, 'end_idx'])
            last_model_end_idx = int(model_start_end.loc[len(model_start_end) - 1, 'end_idx'])
            raw_df = self.get_model(first_model).df

        for record in ['ATOM', 'HETATM']:
            if record in raw_df.keys() and len(raw_df[record]) > 0:
                # left pad residue names to 3 characters
                raw_df[record].loc[:, 'residue_name'] = raw_df[record].residue_name.apply(
                    lambda x: str(x).rjust(3, ' ')[:3]
                )

                # correct malformed hydrogen atom name patterns
                pattern = re.compile(r"^[0-9]H[A-Za-z][0-9]$")
                raw_df[record].loc[:, 'atom_name'] = raw_df[record].atom_name.apply(
                    lambda x: x[1:] + x[0] if pattern.match(x) else x
                )

        if "OTHERS" in raw_df.keys():
            # slice out model specific OTHERS records
            raw_df["OTHERS"] = raw_df["OTHERS"].loc[
                raw_df["OTHERS"]["line_idx"].map(
                    lambda x: x < first_model_start_idx or first_model_start_idx <= x <= first_model_end_idx or x > last_model_end_idx,
                    na_action='ignore'
                ), :
            ].reset_index(drop=True)

        return raw_df

    @staticmethod
    def _fetch_pdb(pdb_code):
        """
        Load PDB file from rcsb.org.

        Overwrites existing method in order to propagate raised HTTPErrors and URLErrors.
        """
        url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
        response = urlopen(url)
        txt = response.read()
        txt = (
            txt.decode("utf-8") if sys.version_info[0] >= 3 else txt.encode("ascii")
        )
        return url, txt

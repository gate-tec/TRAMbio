# TRAMbio
# Author: Nicolas Handke <nicolas@bioinf.uni-leipzig.de>
# License: MIT
# Code Repository: https://github.com/gate-tec/TRAMbio
"""
TRAMbio
=======

*Topological Rigidity Analysis in Molecular Biology* (TRAMbio)
is a package based on and centered on the pebble game. It provides
functionality for general applications in graph theory including testing
for $(k,l)$-sparsity as well as determining the $(k,l)$-spanning
subgraphs, i.e., the $(k,l)$-rigid components. With regard to
molecular data, in particular proteins, TRAMbio provides tools for the
rigidity analysis on an atom or residue-level with further functionality
towards specialized tasks like (a) simulated protein unfolding
(`Rader et al. 2002`_)
on single-state protein data and (b) time-based tracking of
rigid component changes in molecular dynamics (MD) trajectories.

 .. _Rader et al. 2002:
    https://doi.org/10.1073/pnas.062492699
"""

__author__ = "Nicolas Handke <nicolas@bioinf.uni-leipzig.de>"
__version__ = "0.1.2"
__docformat__ = "numpy"

import sys
from typing import Optional

from loguru import logger
from tqdm import tqdm

# Apply patches
import TRAMbio.util.patches

def set_log_level(log_level: Optional[str]):
    if log_level is None or log_level == 'NONE':
        logger.disable('TRAMbio')
    else:
        try:
            logger.remove(None)
        except ValueError:
            pass
        try:
            logger.remove(0)
        except ValueError:
            pass
        logger.add(sys.stdout,
                   format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level}</level> | {message}",
                   filter=lambda record: "task" not in record["extra"].keys() or record["extra"]["task"] != "tqdm",
                   level=log_level)
        logger.add(lambda msg: tqdm.write(msg, end=""), colorize=True,
                   format=" {message}",
                   filter=lambda record: "task" in record["extra"].keys() and record["extra"]["task"] == "tqdm",
                   enqueue=True,
                   level=log_level)
        logger.enable('TRAMbio')

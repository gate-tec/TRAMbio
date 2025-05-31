from enum import Enum
import textwrap
from nicegui import ui


class ClassName(Enum):
    VIEWER = "tram-viewer"
    HIDDEN = "hidden"


def init_header():
    ui.add_head_html(textwrap.dedent(f"""\
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol.min.js"></script>
        <style>
            .{ClassName.VIEWER.value} {{
                position: relative;
                width: 100%;
                min-height: 650px;
            }}
        </style>
        <script type="text/javascript">const viewers={{}};const models={{}};</script>
    """))

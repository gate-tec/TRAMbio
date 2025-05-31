import textwrap
from typing import List, Union, TextIO
from enum import Enum

from nicegui import ui

from TRAMbio.gui.viewer.config import ClassName


class ViewerStyle(Enum):
    DEFAULT = "DEFAULT"
    COLORED_COMPONENTS = "COLORED_COMPONENTS"


class Viewer:

    def __init__(self, multi_model: bool = False, style: ViewerStyle = ViewerStyle.DEFAULT):
        self._current_model_idx = 1
        self.target_model_idx = 1
        self._max_models = 0

        self.multi_model = multi_model
        self.viewer_style = style

        self.viewer = None
        self._viewer_id = None
        self._pagination = None
        self._slider = None

    @property
    def has_model(self) -> bool:
        return self._max_models > 0

    def create_viewer(self):
        self.viewer = ui.element("div").classes(ClassName.VIEWER.value)
        self._viewer_id = self.viewer.html_id

    def _initialize_viewer(self):
        multi_model_string = "true" if self.multi_model else "false"

        ui.run_javascript(textwrap.dedent(f"""\
        viewers.{self._viewer_id} = $3Dmol.createViewer('{self._viewer_id}', {{onemol: true, multimodel: {multi_model_string}}});
        viewers.{self._viewer_id}.setBackgroundColor(0x000000);
        models.{self._viewer_id} = [];"""))

    def create_slider(self):
        if self._slider is not None:
            raise ValueError("Unable to create multiple sliders for one viewer.")
        self._slider = ui.slider(min=1, max=0).bind_value(self, 'target_model_idx')

    def create_pagination(self):
        if self._pagination is not None:
            raise ValueError("Unable to create multiple paginations for one viewer.")
        self._pagination = ui.pagination(
            min=1, max=0, on_change=lambda x: self.change_model(),
            direction_links=True
        ).bind_value(self, 'target_model_idx')

    def add_models(self, files: List[TextIO]):
        first_models = False
        if self._max_models == 0:
            if not self.multi_model and len(files) > 1:
                raise ValueError
            self._initialize_viewer()
            first_models = True
        elif not self.multi_model:
            # clear models
            ui.run_javascript(textwrap.dedent(f"""\
            viewers.{self._viewer_id}.removeAllModels()
            viewers.{self._viewer_id}.clear();
            models.{self._viewer_id} = [];
            viewers.{self._viewer_id}.stopAnimate();
            """))
            self._max_models = 0
            first_models = True

        for file in files:
            data = []
            num_models = 0
            for line in file.readlines():
                if line.startswith("MODEL"):
                    num_models += 1
                data.append(line.strip() + "\\n")
            self._max_models += (1 if num_models == 0 else num_models)

            data = "".join(data)

            if self.multi_model:
                ui.run_javascript(textwrap.dedent(f"""\
                            let new_models = viewers.{self._viewer_id}.addModels("{data}", 'pdb');

                            new_models.forEach((model) => {{
                                model.hide();
                                models.{self._viewer_id}.push(model);
                            }});"""))
            else:
                ui.run_javascript(textwrap.dedent(f"""\
                            let new_model = viewers.{self._viewer_id}.addModel("{data}", 'pdb');

                            models.{self._viewer_id}.push(new_model);"""))

        if self.viewer_style == ViewerStyle.DEFAULT:
            ui.run_javascript(f"viewers.{self._viewer_id}.setStyle({{stick:{{}}}});")
        elif self.viewer_style == ViewerStyle.COLORED_COMPONENTS:
            ui.run_javascript(textwrap.dedent(f"""\
            viewers.{self._viewer_id}.setStyle({{b: 0}},{{stick:{{color:'white'}}}});
            viewers.{self._viewer_id}.setStyle({{not: {{b: 0}}}},{{stick:{{color:'white'}},sphere:{{scale:0.5,opacity:0.75,colorscheme:{{prop:'b',gradient: new $3Dmol.Gradient.Sinebow($3Dmol.getPropertyRange(viewers.{self._viewer_id}.selectedAtoms(),'b'))}}}}}});
            """))

        if first_models:
            ui.run_javascript(textwrap.dedent(f"""\
            models.{self._viewer_id}[0].show();
            viewers.{self._viewer_id}.zoomTo();
            viewers.{self._viewer_id}.render();"""))

        if self._pagination is not None:
            self._pagination.max = self._max_models
            self._pagination.update()
        if self._slider is not None:
            self._slider.max = self._max_models
            self._slider.props['max'] = self._max_models
            self._slider.update()

    def change_model(self):
        if self._current_model_idx == self.target_model_idx:
            return
        # - hide current model
        # - show target model
        # - re-render
        ui.run_javascript(textwrap.dedent(f"""\
        models.{self._viewer_id}[{self._current_model_idx - 1}].hide();
        models.{self._viewer_id}[{self.target_model_idx - 1}].show();
        viewers.{self._viewer_id}.render();
        """))

        self._current_model_idx = self.target_model_idx

from io import StringIO


class ToggledStringIO(StringIO):

    __toggle = False

    @property
    def toggle(self):
        return self.__toggle

    @toggle.setter
    def toggle(self, var_toggle: bool):
        self.__toggle = var_toggle

    def close(self):
        if not self.__toggle:
            super().close()

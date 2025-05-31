from nicegui import app, ui


class ShutdownDialog(ui.dialog):
    def __init__(self):
        super().__init__()

        with self, ui.card():
            with ui.row().classes('w-full justify-end items-center'):
                ui.icon('info', color='primary').classes('text-5xl')
                ui.label("Shutdown confirmation").classes('text-2xl')
            ui.label("Are you sure you want to shutdown?")
            with ui.row().classes('w-full justify-end'):
                ui.button('Confirm Shutdown', on_click=lambda x: self.submit(True))
                ui.button('Cancel', on_click=self.close).props('outline')


async def shutdown() -> None:
    result = await ShutdownDialog()
    if result is not None and result:
        app.shutdown()


# https://stackoverflow.com/questions/73257387/how-to-un-pyinstaller-converted-python-app-with-shiny-for-python
# install pyinstaller
# create spec-file: pyinstaller -F file2.py
# manually edit the spec file (include pathes to shiny)
# create exe: pyinstaller file2.spec


from shiny import App, render, ui
print("Hello world")
app_ui = ui.page_fluid(
    ui.h2("Hello Shiny!"),
    ui.input_slider("n", "N", 0, 100, 20),
    ui.output_text_verbatim("txt"),
)
def server(input, output, session):
    @output
    @render.text
    def txt():
        return f"n*2 is {input.n() * 2}"
app = App(app_ui, server)
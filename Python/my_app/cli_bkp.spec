# -*- mode: python ; coding: utf-8 -*-

block_cipher = None
import os
                       # use OS equivalent for path below
shiny = os.path.abspath("../venv/Lib/site-packages/shiny")
shinywidgets = os.path.abspath("../venv/Lib/site-packages/shinywidgets")
ipywidgets = os.path.abspath("../venv/Lib/site-packages/ipywidgets")
networkx = os.path.abspath("../venv/Lib/site-packages/networkx")
plotly = os.path.abspath("../venv/Lib/site-packages/plotly")
_plotly_utils = os.path.abspath("../venv/Lib/site-packages/_plotly_utils")
fpdf = os.path.abspath("../venv/Lib/site-packages/fpdf")
fontTools = os.path.abspath("../venv/Lib/site-packages/fontTools")
tenacity =  os.path.abspath("../venv/Lib/site-packages/tenacity")
kaleido =   os.path.abspath("../venv/Lib/site-packages/kaleido")

a = Analysis(
    ['cli.py'],
    pathex=[],
    binaries=[],
    datas=[('app.py', './'), ('../myplotting.py', './'), ('../mypdf.py', './'),
            (shiny,'./shiny'), (shinywidgets,'./shinywidgets'), (ipywidgets,'./ipywidgets'),
            (networkx,'./networkx'),
            (plotly,'./plotly'), (_plotly_utils,'./_plotly_utils'),
            (fpdf,'./fpdf'), (fontTools,'./fontTools'), (tenacity,'./tenacity'), (kaleido,'./kaleido')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='cli',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

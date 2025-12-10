# -*- mode: python ; coding: utf-8 -*-

block_cipher = None
import os

shiny = os.path.abspath("./venv/Lib/site-packages/shiny")
shinywidgets = os.path.abspath("./venv/Lib/site-packages/shinywidgets")
ipywidgets = os.path.abspath("./venv/Lib/site-packages/ipywidgets")
jupyter_core = os.path.abspath("./venv/Lib/site-packages/jupyter_core")
dateutil = os.path.abspath("./venv/Lib/site-packages/dateutil")

a = Analysis(
    ['cli.py'],
    pathex=[],
    binaries=[('ADVENS_Fluid.exe', './')],
    datas=[('app.py', './'), ('../dataBase/DataFluids.txt', 'dataBase/'),
            (shiny,'./shiny'), (shinywidgets,'./shinywidgets'), (ipywidgets,'./ipywidgets'),
            (jupyter_core,'./jupyter_core'), (dateutil,'./dateutil')],
    hiddenimports=['comm', 'platformdirs', 'plotly', '_plotly_utils'],
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
    [],
    exclude_binaries=True,
    name='myspecfile',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='myspecfile',
)


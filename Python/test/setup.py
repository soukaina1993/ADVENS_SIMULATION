from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need
# fine tuning.
build_options = {'packages': ['shiny', 'uvicorn', 'anyio'], 'excludes': [], 'include_files': ['file1.py']}

base = 'console'

executables = [
    Executable('file2.py', base=base, target_name = 'Test.exe')
]

setup(name='Test',
      version = '1',
      description = 'First test',
      options = {'build_exe': build_options},
      executables = executables)

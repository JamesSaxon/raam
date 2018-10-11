import platform, distutils.core, distutils.extension, Cython.Build

import sys

## Macs require this extra build option.
ouff_mac = []
if sys.platform == "darwin":
  ouff_mac = ['-mmacosx-version-min=10.9']

EXTENSION = distutils.extension.Extension(
    name = 'pyraam', language = 'c++',
    sources = ['pyraam.pyx'],
    extra_compile_args = ['-Wno-unused-function', 
                          '-std=c++11', '-Wall'] + ouff_mac,
    undef_macros       = ["NDEBUG"],
    extra_link_args    = ouff_mac,
    )

EXT_MODULES=Cython.Build.cythonize([EXTENSION], language='c++')

distutils.core.setup(name = 'simple', ext_modules=EXT_MODULES,)


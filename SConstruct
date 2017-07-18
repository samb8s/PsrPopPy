# Starter SConstruct for enscons

import sys, os
from distutils import sysconfig
import pytoml as toml
import enscons

metadata = dict(toml.load(open('pyproject.toml')))['tool']['enscons']

# most specific binary, non-manylinux1 tag should be at the top of this list
import wheel.pep425tags
full_tag = '-'.join(next(tag for tag in wheel.pep425tags.get_supported() if not 'manylinux' in tag))

# full_tag = py2.py3-none-any # pure Python packages compatible with 2+3

env = Environment(tools=['default', 'packaging', enscons.generate],
                  PACKAGE_METADATA=metadata,
                  WHEEL_TAG=full_tag)

# set the compiler
env['CC'] = 'gfortran'

py_source = (Glob('lib/python/*.py') + ['lib/__init__.py'])

libpath = os.path.join('lib', 'fortran')

env.Append(CPPFLAGS=['-O2', '-fPIC', '-fno-second-underscore', '-c', '-std=legacy'])
env.Append(CPPPATH=[libpath])

# compile NE2001 library
ne2001_libname = 'libne2001'
ne2001_lib = os.path.join(libpath, ne2001_libname)
ne2001_files = ['ne2001.f', 'dm.f', 'psr_ne.f', 'dist.f', 'calc_xyz.f', 'density.f', 'glun.f']
ne2001_source = [os.path.join(libpath, srcfile) for srcfile in ne2001_files]

sharedlib = env.SharedLibrary(target=ne2001_lib, source=ne2001_source)
#staticlib = env.StaticLibrary(target=ne2001_lib, source=ne2001_source)

py_source += sharedlib

# compile ykarea library
ykarea_libname = 'libykarea'
ykarea_lib =  os.path.join(libpath, ykarea_libname)
ykarea_files = ['ykarea.f', 'psrran.f']
ykarea_source = [os.path.join(libpath, srcfile) for srcfile in ykarea_files]

sharedlib = env.SharedLibrary(target=ykarea_lib, source=ykarea_source)
#staticlib = env.StaticLibrary(target=ykarea_lib, source=ykarea_source)

py_source += sharedlib

# compile sla library
sla_libname = 'libsla'
sla_lib = os.path.join(libpath, sla_libname)
sla_files = ['galtfeq.f', 'sla.f']
sla_source = [os.path.join(libpath, srcfile) for srcfile in sla_files]

sharedlib = env.SharedLibrary(target=sla_lib, source=sla_source)
#staticlib = env.StaticLibrary(target=sla_lib, source=sla_source)

py_source += sharedlib

# compile vyyz library
vxyz_libname = 'libvxyz'
vxyz_lib = os.path.join(libpath, vxyz_libname)
vxyz_files = ['vxyz.f', 'rkqc.f', 'rk4.f']
vxyz_source = [os.path.join(libpath, srcfile) for srcfile in vxyz_files]

sharedlib = env.SharedLibrary(target=vxyz_lib, source=vxyz_source)
#staticlib = env.StaticLibrary(target=vxyz_lib, source=vxyz_source)

py_source += sharedlib

# compile gamma library
gamma_libname = 'libgamma'
gamma_lib = os.path.join(libpath, gamma_libname)
gamma_files = ['gamma.f']
gamma_source = [os.path.join(libpath, srcfile) for srcfile in gamma_files]

sharedlib = env.SharedLibrary(target=gamma_lib, source=gamma_source)
#staticlib = env.StaticLibrary(target=gamma_lib, source=gamma_source)

py_source += sharedlib
py_source += Glob('lib/fortran/*.*') + Glob('lib/fortran/lookuptables/*') + Glob('lib/python/models/*') + Glob('lib/surveys/*')

platlib = env.Whl('platlib', py_source, root='')
whl = env.WhlFile(source=platlib)

# Add automatic source files, plus any other needed files.
sdist_source=list(set(FindSourceFiles() +
                  ['PKG-INFO', 'setup.py'] +
                  Glob('lib/fortran/*.*') + Glob('lib/fortran/lookuptables/*') + Glob('lib/python/models/*') + Glob('lib/surveys/*')))

sdist_source += py_source

sdist = env.SDist(source=sdist_source)
env.Alias('sdist', sdist)

install = env.Command("#DUMMY", whl, ' '.join([sys.executable, '-m', 'pip', 'install', 'psrpoppy', '--no-deps', '--no-index', '--find-links', '$SOURCE']))
env.Alias('install', install)
env.AlwaysBuild(install)

env.Default(sdist)


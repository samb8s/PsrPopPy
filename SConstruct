# Starter SConstruct for enscons

import sys, os
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

# set the compiler to gfortran
env['CC'] = 'gfortran'

py_source = (Glob('psrpoppy/python/*.py') + ['psrpoppy/__init__.py'])

libpath = os.path.join('psrpoppy', 'fortran')

# check whether installing on a Mac
if 'darwin' in os.environ['OSTYPE']:
    env.Append(CFLAGS=['-m 32'])
    env.Append(CPPFLAGS=['-dynamiclib', '-O2', '-fPIC', '-fno-second-underscore', '-c', '-std=legacy'])
else:
    env.Append(CPPFLAGS=['-O2', '-fPIC', '-fno-second-underscore', '-c', '-std=legacy'])

env.Append(CPPPATH=[libpath])

# dictionary of libraries and files needed for library
LIBDIC = {}
LIBDIC['libne2001']  = ['ne2001.f', 'dm.f', 'psr_ne.f', 'dist.f', 'calc_xyz.f', 'density.f', 'glun.f']
LIBDIC['libykarea']  = ['ykarea.f', 'psrran.f']
LIBDIC['libsla']     = ['galtfeq.f', 'sla.f']
LIBDIC['libvxyz']    = ['vxyz.f', 'rkqc.f', 'rk4.f']
LIBDIC['libgamma']   = ['gamma.f']
LIBDIC['libgetseed'] = ['getseed.f', 'clock.f']

libs = []

# compile libraries
for libname in LIBDIC:
    lib = os.path.join(libpath, libname)
    libsources = [os.path.join(libpath, srcfile) for srcfile in LIBDIC[libname]]
    
    sharedlib = env.SharedLibrary(target=lib, source=libsources)
    #staticlib = env.StaticLibrary(target=lib, source=libsources)

    libs += sharedlib

otherfiles = Glob('psrpoppy/fortran/*.so') + Glob('psrpoppy/fortran/lookuptables/*') + Glob('psrpoppy/python/models/*') + Glob('psrpoppy/surveys/*')

platlib = env.Whl('platlib', py_source + libs + otherfiles, root='')
whl = env.WhlFile(source=platlib)

# Add automatic source files, plus any other needed files.
sdist_source=list(set(FindSourceFiles() +
                  ['PKG-INFO', 'setup.py'] +
                  Glob('psrpoppy/fortran/*.*') + Glob('psrpoppy/fortran/lookuptables/*') + Glob('psrpoppy/python/models/*') + Glob('psrpoppy/surveys/*')))

sdist_source += py_source

sdist = env.SDist(source=sdist_source)
env.Alias('sdist', sdist)

install = env.Command("#DUMMY", whl, ' '.join([sys.executable, '-m', 'pip', 'install', '--no-deps', '$SOURCE']))
env.Alias('install', install)
env.AlwaysBuild(install)

env.Default(sdist)


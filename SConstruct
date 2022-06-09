# Starter SConstruct for enscons

import sys
import os
import pytoml as toml
import enscons

# check if prefix is set
AddOption('--prefix', dest='prefix', type='string', nargs=1,
          action='store', metavar='DIR', default='/usr/local',
          help='Installation Prefix')

AddOption('--user', dest='user', action='store_true', default=False,
          help='Install as pip-like "--user" (this overrides "--prefix")')

metadata = dict(toml.load(open('pyproject.toml')))['tool']['enscons']

full_tag = enscons.get_binary_tag()

# full_tag = py2.py3-none-any # pure Python packages compatible with 2+3

env = Environment(tools=['default', 'packaging', enscons.generate],
                  PACKAGE_METADATA=metadata,
                  WHEEL_TAG=full_tag)

# set the compiler to gfortran
env['CC'] = 'gfortran'

py_source = Glob('psrpoppy/*.py')

libpath = os.path.join('psrpoppy', 'fortran')

# check whether installing on a Mac
if 'Darwin' in os.uname()[0]:
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

# install prefix
if not GetOption('user'):
    installprefix=GetOption('prefix')
else: # install in --user location
    installprefix=os.path.join(os.environ['HOME'], '.local')
pyprefix='psrpoppy'
executables = ['dosurvey', 'evolve', 'populate']

# install executables
insbins = env.InstallAs(target=[os.path.join(installprefix, 'bin', ex) for ex in executables],
               source=[os.path.join(pyprefix, ex+'.py') for ex in executables])

otherfiles = Glob('psrpoppy/fortran/*.so') + Glob('psrpoppy/fortran/lookuptables/*') + Glob('psrpoppy/models/*') + Glob('psrpoppy/surveys/*')

platlib = env.Whl('platlib', py_source + libs + otherfiles, root='')
whl = env.WhlFile(source=platlib)

# Add automatic source files, plus any other needed files.
sdist_source=list(set(FindSourceFiles() +
                  ['PKG-INFO', 'setup.py'] +
                  Glob('psrpoppy/fortran/*.f') + Glob('psrpoppy/fortran/*.inc') + Glob('psrpoppy/fortran/lookuptables/*') + Glob('psrpoppy/models/*') + Glob('psrpoppy/surveys/*')))

sdist_source += py_source

sdist = env.SDist(source=sdist_source)
env.Alias('sdist', sdist)

if GetOption('user'):
    install = env.Command("#DUMMY", whl, ' '.join([sys.executable, '-m', 'pip', 'install', '--no-deps', '--user', '$SOURCE']))
else:
    install = env.Command("#DUMMY", whl, ' '.join(['PYTHONUSERBASE={}'.format(installprefix), sys.executable, '-m', 'pip', 'install', '--no-deps', '--user', '$SOURCE']))
env.Alias('install', install + insbins)
env.AlwaysBuild(install + insbins)

env.Default(sdist)


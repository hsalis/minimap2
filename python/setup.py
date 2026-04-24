from pathlib import Path
import os
import platform

os.environ.setdefault("SETUPTOOLS_USE_DISTUTILS", "stdlib")
if platform.system() == "Darwin" and platform.machine() in ("aarch64", "arm64"):
    os.environ.setdefault("ARCHFLAGS", "-arch arm64")

from setuptools import Extension, setup

try:
    from Cython.Build import cythonize
except ImportError as exc:
    raise SystemExit("Cython is required to build mappy from python/setup.py") from exc

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

extra_compile_args = ['-DHAVE_KALLOC']
include_dirs = [str(ROOT)]

if platform.machine() in ("aarch64", "arm64"):
    include_dirs.append(str(ROOT / 'sse2neon'))
    extra_compile_args.extend(['-ftree-vectorize', '-DKSW_SSE2_ONLY', '-D__SSE2__'])
else:
    extra_compile_args.append('-msse4.1')

sources = [
    str(HERE / 'mappy.pyx'),
    str(ROOT / 'align.c'), str(ROOT / 'bseq.c'), str(ROOT / 'lchain.c'),
    str(ROOT / 'seed.c'), str(ROOT / 'format.c'), str(ROOT / 'hit.c'),
    str(ROOT / 'index.c'), str(ROOT / 'pe.c'), str(ROOT / 'jump.c'),
    str(ROOT / 'options.c'), str(ROOT / 'ksw2_extd2_sse.c'),
    str(ROOT / 'ksw2_exts2_sse.c'), str(ROOT / 'ksw2_extz2_sse.c'),
    str(ROOT / 'ksw2_ll_sse.c'), str(ROOT / 'kalloc.c'),
    str(ROOT / 'kthread.c'), str(ROOT / 'map.c'), str(ROOT / 'misc.c'),
    str(ROOT / 'sdust.c'), str(ROOT / 'sketch.c'), str(ROOT / 'esterr.c'),
    str(ROOT / 'splitidx.c'),
]

depends = [
    str(ROOT / 'minimap.h'), str(ROOT / 'bseq.h'), str(ROOT / 'kalloc.h'),
    str(ROOT / 'kdq.h'), str(ROOT / 'khash.h'), str(ROOT / 'kseq.h'),
    str(ROOT / 'ksort.h'), str(ROOT / 'ksw2.h'), str(ROOT / 'kthread.h'),
    str(ROOT / 'kvec.h'), str(ROOT / 'mmpriv.h'), str(ROOT / 'sdust.h'),
    str(HERE / 'cmappy.h'), str(HERE / 'cmappy.pxd'),
]

extension = Extension(
    'mappy',
    sources=sources,
    depends=depends,
    extra_compile_args=extra_compile_args,
    include_dirs=include_dirs,
    libraries=['z', 'm', 'pthread'],
)

setup(
    name='mappy',
    version='2.30',
    url='https://github.com/lh3/minimap2',
    description='Minimap2 python binding',
    long_description=(HERE / 'README.rst').read_text(),
    author='Heng Li',
    author_email='lh3@me.com',
    license='MIT',
    keywords='sequence-alignment',
    scripts=[str(HERE / 'minimap2.py')],
    ext_modules=cythonize([extension], compiler_directives={'language_level': '3'}),
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX',
        'Programming Language :: C',
        'Programming Language :: Cython',
        'Programming Language :: Python :: 3',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)

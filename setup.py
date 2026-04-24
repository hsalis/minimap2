from pathlib import Path
import os
import platform

os.environ.setdefault("SETUPTOOLS_USE_DISTUTILS", "stdlib")
if platform.system() == "Darwin" and platform.machine() in ("aarch64", "arm64"):
    os.environ.setdefault("ARCHFLAGS", "-arch arm64")

from setuptools import Extension, setup

ROOT = Path(__file__).resolve().parent
PYTHON_DIR = ROOT / "python"
README = ROOT / "README.md"

use_cython = False
module_source = "python/mappy.c"
try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None
else:
    if (PYTHON_DIR / "mappy.pyx").exists():
        use_cython = True
        module_source = "python/mappy.pyx"

extra_compile_args = ["-DHAVE_KALLOC"]
include_dirs = [".", "python"]

if platform.machine() in ("aarch64", "arm64"):
    include_dirs.append("sse2neon")
    extra_compile_args.extend(["-ftree-vectorize", "-DKSW_SSE2_ONLY", "-D__SSE2__"])
else:
    extra_compile_args.append("-msse4.1")

sources = [
    module_source,
    "align.c", "bseq.c", "lchain.c", "seed.c", "format.c", "hit.c",
    "index.c", "pe.c", "jump.c", "options.c", "ksw2_extd2_sse.c",
    "ksw2_exts2_sse.c", "ksw2_extz2_sse.c", "ksw2_ll_sse.c", "kalloc.c",
    "kthread.c", "map.c", "misc.c", "sdust.c", "sketch.c", "esterr.c",
    "splitidx.c",
]

depends = [
    "minimap.h", "bseq.h", "kalloc.h", "kdq.h", "khash.h", "kseq.h",
    "ksort.h", "ksw2.h", "kthread.h", "kvec.h", "mmpriv.h", "sdust.h",
    "python/cmappy.h", "python/cmappy.pxd",
]
if use_cython:
    depends.append("python/mappy.pyx")
else:
    depends.append("python/mappy.c")

extension = Extension(
    "mappy",
    sources=sources,
    depends=depends,
    extra_compile_args=extra_compile_args,
    include_dirs=include_dirs,
    libraries=["z", "m", "pthread"],
)

ext_modules = [extension]
if use_cython:
    ext_modules = cythonize(ext_modules, compiler_directives={"language_level": "3"})

setup(
    name="mappy",
    version="2.30",
    url="https://github.com/lh3/minimap2",
    description="Minimap2 Python bindings and CLI",
    long_description=README.read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    author="Heng Li",
    author_email="lh3@me.com",
    license="MIT",
    keywords="sequence-alignment bioinformatics minimap2",
    py_modules=["minimap2"],
    package_dir={"": "python"},
    entry_points={
        "console_scripts": [
            "minimap2py=minimap2:main",
        ]
    },
    ext_modules=ext_modules,
    install_requires=["numpy"],
    python_requires=">=3.9",
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX",
        "Programming Language :: C",
        "Programming Language :: Python :: 3",
        "Programming Language :: Cython",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)

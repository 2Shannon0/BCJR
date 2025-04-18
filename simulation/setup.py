from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "bcjr_decoder",
        ["bcjr_decoder.cpp"],
        include_dirs=[pybind11.get_include()],
        language="c++",
        extra_compile_args=["-std=c++17"],  # Добавили поддержку C++17
    ),
]

setup(
    name="bcjr_decoder",
    ext_modules=ext_modules,
)

# python3 setup.py build_ext --inplace

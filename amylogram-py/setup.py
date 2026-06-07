from setuptools import find_packages, setup


setup(
    name="amylogram-py",
    version="0.1.0",
    description="Fast Python-compatible AmyloGram predictor implementation",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    author="AmyloGram-Py contributors",
    license="GPL-3.0-or-later",
    python_requires=">=3.10",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    entry_points={"console_scripts": ["amylogram-py=amylogram_py.cli:main"]},
)

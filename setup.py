import setuptools

with open("requirements.txt", "r") as fileobj:
    requirements = fileobj.readlines()

with open("README.md", "r") as fileobj:
    readme = fileobj.read()

setuptools.setup(
    name="ase_cantera_microkinetics",
    version="0.1.0",
    url="https://github.com/raffaelecheula/ase_cantera_microkinetics.git",
    author="Raffaele Cheula",
    author_email="cheula.raffaele@gmail.com",
    description=(
        "Microkinetic modeling with ASE and Cantera."
    ),
    long_description=readme,
    license="GPL-3.0",
    packages=[
        "ase_cantera_microkinetics",
    ],
    package_dir={
        "ase_cantera_microkinetics": "ase_cantera_microkinetics"
    },
    install_requires=requirements,
    python_requires=">=3.5, <4",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
    ],
)

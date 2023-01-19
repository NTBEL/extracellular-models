import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="extmodels",
    version="0.2.0",
    python_requires='>=3.9',
    install_requires=['numpy', 'scipy', 'matplotlib', 'neuron'],
    author="Blake A. Wilson",
    author_email="blake.wilson@utdallas.edu",
    description="Reaction-diffusion models to simulate the dynamics of fluorescent dyes and peptides in the brain extracellular space implemented using the extracellular simulator in NEURON.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NTBEL/extracellular-models",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
)

from setuptools import setup, find_packages

setup(
    name="molforge",
    version="2.0.0",
    packages=find_packages(),
    install_requires=[
        'pandas>=2.2.0',      # Data manipulation
        'numpy>=1.26.3',      # Numerical operations (explicit dependency)
        'scipy>=1.16.3',      # Statistical analysis (distributions actor)
        'rdkit>=2025.9.1',    # Cheminformatics toolkit
        'requests>=2.31.0',   # HTTP requests (ChEMBL database downloader)
        'tqdm>=4.66.2',       # Progress bars (database downloads)
        'joblib>=1.3.2',      # Parallel processing (RDKit backend uses cloudpickle)
        'tabulate>=0.9.0',    # Formatted table output
    ],
    extras_require={
        'molbox': ['molbox'],                    # Molbox integration
        'cli': [
            'pyyaml>=5.4',                      # YAML config file support
            'questionary>=2.0.0',               # Interactive CLI prompts
        ],
        'openeye': ['openeye-toolkits'],        # OpenEye OMEGA backend (optional)
        'dev': [                                # Development tools
            'pytest>=7.0.0',
            'pytest-mock>=3.6.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'molforge=molforge.cli.main:main',
        ],
    },
    python_requires='>=3.12.1',
)
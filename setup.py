from setuptools import setup, find_packages

VERSION = '0.0.1'
DESCRIPTION = 'A package to perform Two-sample Mendelian Randomization'
LONG_DESCRIPTION = 'A package to perform Two-sample Mendelian Randomization using Torch'

setup(
    name="pyTWMR",
    version=VERSION,
    author="Sergey Oreshkov",
    author_email="Sergey.Oreshkov@unil.ch",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=['numpy', 'torch' , 'pandas'],
    entry_points={
        'console_scripts': [
            'TWMR = pyTWMR.__main__:main',
            'RevTWMR = pyRevTWMR.cli:main'
        ]
    },
    keywords=['python', 'MR', 'Mendelian Randomization', 'Torch', 'GPU'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
    ]
)
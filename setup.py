#!/usr/bin/env python

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'argparse',
    'h5py',
    'numpy',
    'pandas',
    'pyfaidx',
    'pysam',
    'scikit-allel',
    'zarr',
]

test_requirements = [
    'pip==19.2.3',
    'bump2version==0.5.11',
    'wheel==0.33.6',
    'watchdog==0.9.0',
    'flake8==3.7.8',
    'tox==3.14.0',
    'coverage==4.5.4',
    'Sphinx==1.8.5',
    'twine==1.14.0',
]

setup(
    author="Nace Kranjc",
    author_email='n.kranjc@imperial.ac.uk',
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Guido is a fashionable gRNA designer, specialised in gene drives.",
    entry_points={
        'console_scripts': [
            'guido=guido.guido:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='guido',
    name='guido',
    packages=find_packages(include=['guido', 'guido.*']),
    setup_requires=requirements,
    test_suite='tests',
    tests_require=requirements,
    url='https://github.com/nkran/guido',
    version='0.1.0',
    zip_safe=False,
)

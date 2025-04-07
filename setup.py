from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    'jsonschema>=4.23.0',
    'numpy>=2.2.4',
    'netCDF4>=1.7.2',
    'pandas>=2.2.3',
    'xarray>=2023.1.0'
]

TESTS_REQUIRE = [
    'pytest',
    'setuptools_scm',
]

EXTRAS_REQUIRE = {
    'testing': TESTS_REQUIRE
}

PACKAGE_DATA = {
    'aodntools.ncwriter': ['*.json'],
    'aodntools.timeseries_products': ['*.json']
}

PACKAGE_EXCLUDES = ['test_aodntools.*', 'test_aodntools']
PACKAGE_NAME = 'aodntools'

setup(
    name=PACKAGE_NAME,
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=find_packages(exclude=PACKAGE_EXCLUDES),
    package_data=PACKAGE_DATA,
    url='https://github.com/aodn',
    license='GPLv3',
    author='AODN',
    author_email='projectofficers@emii.org.au',
    description='AODN data tools library',
    zip_safe=False,
    python_requires='>=3.11, <3.12',
    install_requires=INSTALL_REQUIRES,
    tests_require=TESTS_REQUIRE,
    extras_require=EXTRAS_REQUIRE,
    test_suite='test_aodntools',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: Implementation :: CPython',
    ]
)

from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    'jsonschema>=2.6.0,<3.0.0',
    'numpy>=1.13.0',
    'netCDF4>=1.5.3',
    'pandas>=0.24.2',
    'xarray>=0.11.3'
]

TESTS_REQUIRE = [
    'pytest'
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
    version='0.0.0',
    packages=find_packages(exclude=PACKAGE_EXCLUDES),
    package_data=PACKAGE_DATA,
    url='https://github.com/aodn',
    license='GPLv3',
    author='AODN',
    author_email='projectofficers@emii.org.au',
    description='AODN data tools library',
    zip_safe=False,
    python_requires='>=3.5',
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
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython',
    ]
)

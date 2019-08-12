from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    'jsonschema==2.6.0',
    'numpy>=1.13.0',
    'netCDF4'
]

PACKAGE_DATA = {
    'ncwriter': ['*.json']
}

PACKAGE_EXCLUDES = ['ncwriter.*', 'ncwriter']
PACKAGE_NAME = 'aodntools'

setup(
    name=PACKAGE_NAME,
    version='0.2.1',
    packages=find_packages(exclude=PACKAGE_EXCLUDES),
    package_data=PACKAGE_DATA,
    url='https://github.com/aodn',
    license='GPLv3',
    author='AODN',
    author_email='projectofficers@emii.org.au',
    description='AODN data tools library',
    zip_safe=False,
    install_requires=INSTALL_REQUIRES,
    test_suite='test_aodntools'
)

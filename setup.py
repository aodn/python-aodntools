from setuptools import setup, find_packages

import versioneer

INSTALL_REQUIRES = [
    'jsonschema==2.6.0',
    'numpy>=1.13.0',
    'netCDF4'
]

# TODO: add this when we have JSON schema files and templates bundled
# PACKAGE_DATA = {
#     'ncwriter': [
#     ]
# }

PACKAGE_EXCLUDES = ['test_ncwrwiter.*', 'test_ncwriter']
PACKAGE_NAME = 'ncwriter'

setup(
    name=PACKAGE_NAME,
    version=versioneer.get_version(),
    packages=find_packages(exclude=PACKAGE_EXCLUDES),
    # package_data=PACKAGE_DATA,
    url='https://github.com/aodn',
    license='GPLv3',
    author='AODN',
    author_email='projectofficers@emii.org.au',
    description='AODN netCDF tools library',
    zip_safe=False,
    install_requires=INSTALL_REQUIRES,
    test_suite='test_ncwriter'
)

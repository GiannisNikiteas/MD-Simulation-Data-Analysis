from setuptools import setup, find_packages


with open('README.md', 'r') as fh:
    LONG_DESC = fh.read()


setup(
    name='md',
    version='0.0.1',
    packages=['md'],
    url='https://github.com/GiannisNikiteas/MD-Simulation-Data-Analysis',
    license='',
    author='Giannis Nikiteas',
    author_email='giannis_nikiteas@hotmail.com',
    description='A package responsible for Molecular Dynamics data analysis',
    long_description=LONG_DESC,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)




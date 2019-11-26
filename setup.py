from setuptools import setup   # find_packages


with open('README.md', 'r') as fh:
    LONG_DESC = fh.read()


setup(
    name='mdtools',
    version='0.1.0',
    packages=['mdtools'],
    url='https://github.com/GNikit/md-tools',
    license='',
    author='GNikit',
    author_email='giannis.nikiteas@gmail.com',
    description='A package responsible for Molecular Dynamics data analysis',
    long_description=LONG_DESC,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    install_requires=['matplotlib', 'scipy', 'numpy']
)

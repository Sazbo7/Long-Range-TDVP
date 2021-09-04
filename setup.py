from setuptools import setup, find_packages

setup(
    name='Long Range TDVP',
    version='1.0.0',
    url='https://github.com/mypackage.git',
    author='Author Name',
    author_email='jszabo94@gmail.com',
    description='''Simple examples codes extended from TeNpy mpo_exponential_decay (https://tenpy.readthedocs.io/en/latest/examples/advanced/mpo_exponential_decay.html). Evaluate real time evolution following quench from initial MPS states using hybrid 2-site and 1-site TDVP evolution. Calculates entanglement entropy, and 1 and 2 site observables/correlations.''',
    packages=find_packages(),
    install_requires=['numpy >= 1.11.1', 'tenpy>= 0.8.4'],
)

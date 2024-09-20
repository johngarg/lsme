from setuptools import setup, find_packages

setup(
    name='lsme',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[],
    author='John Gargalionis',
    author_email='john.gargalionis@gmail.com',
    description='One-loop matching data for the linear SM extensions',
    url='https://github.com/johngarg/lsme',
    license='GPLv3',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)

from setuptools import setup, find_packages

setup(
    name="IIBacFinder",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'IIBacFinder=scripts.predict:main',
        ],
    },
    author="Zhang Dengwei",
    description="Class II bacteriocin finder",
    python_requires=">=3.6",
)

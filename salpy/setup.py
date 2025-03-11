from setuptools import setup, find_packages

setup(
    name="salpy",
    version="0.0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "geopandas",
        "rasterio",
        "scipy",
        "scikit-learn",
        "pyproj",
        "shapely",
        "dtw-python",
        "deap",  # For genetic algorithms (equivalent to GA in R)
    ],
    author="Narumasa Tsutsumida",
    author_email="rsnaru.jp@gmail.com",
    description="Satellite LiDAR Point Adjustment for Python",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/naru-T/salpy",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)

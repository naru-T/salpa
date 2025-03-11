"""
Minimize Loss Between LiDAR and Reference Data

This module provides functionality to minimize the loss (difference) between
LiDAR point measurements and reference raster data using genetic algorithms.
"""

import random

import numpy as np
from deap import algorithms, base, creator, tools
from geopandas import GeoDataFrame

from .get_loss import get_loss


def minimize_loss(
    lidar_footprints: GeoDataFrame,
    input_rast,
    raster_transform,
    minimizing_method: str = "euclidean",
    target_variable: str = "mean",
    buf: float = 12.5,
    crs_code: int = None,
    lidar_value: str = None,
    lower_bounds: list = None,
    upper_bounds: list = None,
    pop_size: int = 50,
    max_iter: int = 100,
    parallel: bool = False
) -> dict:
    """
    Minimizes the loss between LiDAR point measurements and reference raster data
    using genetic algorithms.

    Parameters
    ----------
    lidar_footprints : GeoDataFrame
        Input satellite LiDAR footprints as a GeoDataFrame
    input_rast : numpy.ndarray
        Input raster data as a numpy array
    raster_transform : affine.Affine
        Affine transformation that maps pixel coordinates to coordinates in the CRS
    minimizing_method : str, optional
        Method for calculating the loss, by default "euclidean"
    target_variable : str, optional
        Statistic to compute for buffered areas, by default "mean"
    buf : float, optional
        Buffer size around each point, by default 12.5
    crs_code : int, optional
        Coordinate reference system code, by default None
    lidar_value : str, optional
        Name of the column containing LiDAR measurement values, by default None
    lower_bounds : list, optional
        Lower bounds for optimization parameters [x, y], by default [-30, -30]
    upper_bounds : list, optional
        Upper bounds for optimization parameters [x, y], by default [30, 30]
    pop_size : int, optional
        Population size for genetic algorithm, by default 50
    max_iter : int, optional
        Maximum number of iterations for genetic algorithm, by default 100
    parallel : bool, optional
        Whether to use parallel processing, by default False

    Returns
    -------
    dict
        Dictionary containing optimization results

    Examples
    --------
    >>> import geopandas as gpd
    >>> import rasterio
    >>> gdf = gpd.read_file('lidar_points.shp')
    >>> with rasterio.open('dem.tif') as src:
    ...     raster = src.read(1)
    ...     transform = src.transform
    >>> result = minimize_loss(
    ...     lidar_footprints=gdf,
    ...     input_rast=raster,
    ...     raster_transform=transform,
    ...     minimizing_method="euclidean",
    ...     target_variable="mean",
    ...     buf=12.5,
    ...     crs_code=3857,
    ...     lidar_value="elevation",
    ...     lower_bounds=[-30, -30],
    ...     upper_bounds=[30, 30],
    ...     pop_size=50,
    ...     max_iter=100,
    ...     parallel=False
    ... )
    """
    # Check if the input GeoDataFrame is empty
    if len(lidar_footprints) == 0:
        print("Skipping because of no high quality data found")
        return None

    # Set default bounds if not provided
    if lower_bounds is None:
        lower_bounds = [-30, -30]
    if upper_bounds is None:
        upper_bounds = [30, 30]

    # Define the objective function
    def objective_function(individual):
        x, y = individual
        result = get_loss(
            lidar_footprints=lidar_footprints,
            add_x=x,
            add_y=y,
            buf=buf,
            input_rast=input_rast,
            raster_transform=raster_transform,
            minimizing_method=minimizing_method,
            target_variable=target_variable,
            crs_code=crs_code,
            lidar_value=lidar_value
        )
        loss_value = result["loss_value"]

        if np.isnan(loss_value):
            print(f"NA values produced at x: {x} y: {y}")
            return 1000,

        return loss_value,

    # Set up the genetic algorithm
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()

    # Register genetic operators
    toolbox.register("attr_float_x", random.uniform, lower_bounds[0], upper_bounds[0])
    toolbox.register("attr_float_y", random.uniform, lower_bounds[1], upper_bounds[1])
    toolbox.register("individual", tools.initCycle, creator.Individual,
                     (toolbox.attr_float_x, toolbox.attr_float_y), n=1)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    # Register the evaluation function
    toolbox.register("evaluate", objective_function)

    # Register genetic operators
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # Set up parallel evaluation if requested
    if parallel:
        import multiprocessing
        pool = multiprocessing.Pool()
        toolbox.register("map", pool.map)

    # Create initial population
    population = toolbox.population(n=pop_size)

    # Set random seed for reproducibility
    random.seed(1118)

    # Run the genetic algorithm
    try:
        algorithms.eaSimple(
            population,
            toolbox,
            cxpb=0.7,  # Crossover probability
            mutpb=0.2,  # Mutation probability
            ngen=max_iter,  # Number of generations
            verbose=False
        )
    except Exception as e:
        print(f"Error during optimization: {e}")
        return None
    finally:
        if parallel:
            pool.close()

    # Get the best individual
    best_ind = tools.selBest(population, 1)[0]
    best_x, best_y = best_ind
    best_value = best_ind.fitness.values[0]

    # Return optimization results
    optim_result = {
        "best_x": best_x,
        "best_y": best_y,
        "best_value": best_value
    }

    return optim_result

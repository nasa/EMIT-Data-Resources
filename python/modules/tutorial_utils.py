"""
This module contains various functions used within the Exploring EMIT
L2B CH4 plume complexes notebook for searching, wrangling data and visualizations.

Author: Erik Bolch, ebolch@contractor.usgs.gov
Last Updated: 2024-02-16
"""

# Imports
from typing import List, Union
import re
import pandas as pd
import geopandas as gpd
import shapely
import earthaccess
import asyncio
import aiohttp


def convert_bounds(bbox, invert_y=False):
    """
    Helper method for changing bounding box representation to leaflet notation

    ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))``
    """
    x1, y1, x2, y2 = bbox
    if invert_y:
        y1, y2 = y2, y1
    return ((y1, x1), (y2, x2))


def flattent_column_names(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [
        re.sub("([A-Z]+)", r"_\1", col.split(".")[-1]).lower() for col in df.columns
    ]
    return df


def get_shapely_object(
    result: earthaccess.search.DataGranule,
) -> Union[shapely.geometry.base.BaseGeometry, None]:
    """
    Retrieve the coordinates from the umm metadata and convert to a shapely geometry.
    """
    shape = None
    try:
        geo = result["umm"]["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]
        keys = geo.keys()
        if "BoundingRectangles" in keys:
            bounding_rectangle = geo["BoundingRectangles"][0]
            # Create bbox tuple
            bbox_coords = (
                bounding_rectangle["WestBoundingCoordinate"],
                bounding_rectangle["SouthBoundingCoordinate"],
                bounding_rectangle["EastBoundingCoordinate"],
                bounding_rectangle["NorthBoundingCoordinate"],
            )
            # Create shapely geometry from bbox
            shape = shapely.geometry.box(*bbox_coords, ccw=True)
        elif "GPolygons" in keys:
            points = geo["GPolygons"][0]["Boundary"]["Points"]
            # Create shapely geometry from polygons
            shape = shapely.geometry.Polygon(
                [[p["Longitude"], p["Latitude"]] for p in points]
            )
        else:
            raise ValueError(
                "Provided result does not contain bounding boxes/polygons or is incompatible."
            )
    except Exception as e:
        pass
    return shape


def list_metadata_fields(results: List[earthaccess.search.DataGranule]) -> List[str]:
    metadata_fields = list(flattent_column_names(pd.json_normalize(results)).columns)
    return metadata_fields


def results_to_geopandas(
    results: List[earthaccess.search.DataGranule],
    fields: List[str] = [],
) -> gpd.GeoDataFrame:
    """
    Convert the results of an earthaccess search into a geodataframe using some default fields.
    Add additional ones with the fields kwarg.
    """
    default_fields = [
        "size",
        "concept_id",
        "dataset-id",
        "native-id",
        "provider-id",
        "_related_urls",
        "_single_date_time",
        "_beginning_date_time",
        "_ending_date_time",
        "geometry",
    ]
    results_df = pd.json_normalize(list(results), errors="ignore")

    results_df = flattent_column_names(results_df)
    if len(fields) == 0:
        fields = default_fields
    else:
        fields = list(set(fields + default_fields))

    results_df = results_df.drop(
        columns=[col for col in results_df.columns if col not in fields]
    )

    results_df["_related_urls"] = results_df["_related_urls"].apply(
        lambda links: [
            link
            for link in links
            if link["Type"]
            in [
                "GET DATA",
                "GET DATA VIA DIRECT ACCESS",
                "GET RELATED VISUALIZATION",
            ]
        ]
    )
    # Create shapely polygons for result
    geometries = [
        get_shapely_object(results[index]) for index in results_df.index.to_list()
    ]
    # Convert to GeoDataframe
    gdf = gpd.GeoDataFrame(results_df, geometry=geometries, crs="EPSG:4326")
    return gdf

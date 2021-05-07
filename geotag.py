#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" Fast geotagging of big data (with a progress bar!)

See README.md or geotag.py --help for more information.
"""

import dataclasses
import typing

import geopandas
import pandas
import shapely.geometry
import rtree
import tqdm

tqdm.tqdm.pandas()


@dataclasses.dataclass
class GeotagInput:
    input_file: str
    input_column: str
    output_column: str


class Geotagger:
    """ Class encapsulating fast lookups of unique identifiers from spatial
    data using an R-tree index.

    Attributes:
        id_strs: A dict that maps the integer representation of each unique ID
            to its original string representation - rtree only supports
            integer keys, so this allows us to get back things like leading
            zeroes that become lost in the type conversion.
        shapes: A dict that maps the integer representation of each unique ID
            to its corresponding shape.
        index: The R-tree index.
    """

    def __init__(self,
                 gdf: geopandas.geodataframe.GeoDataFrame,
                 id_column: str = "GEOID",
                 verbose: bool = False
                 ) -> None:
        """ Initialize Geotagger object.

        Args:
            gdf: A GeoDataFrame containing polygons to use for geotagging.
            id_column: The column of the GeoDataFrame (e.g. field of a
                shapefile) to pull unique IDs from (must be convertable to int).
            verbose: If True, print what is being done.
        """

        self.id_strs = {int(id_): id_ for id_ in gdf[id_column]}
        self.shapes = gdf.set_index(id_column)["geometry"].to_dict()
        self.index = rtree.index.Index()

        iterable = self.shapes.items()
        if verbose:
            iterable = tqdm.tqdm(
                self.shapes.items(),
                "Creating rtree index",
                unit=" indexed"
            )

        for id_, shape in iterable:
            self.index.insert(int(id_), shape.bounds)

    def lookup(self, x: float, y: float) -> typing.Optional[int]:
        """ Look up a coordinate pair's unique ID.

        Args:
            x: The longitude, as a float.
            y: The latitude, as a float.

        Returns:
            The unique ID, if any.
        """
        results = list(self.index.intersection((x, y, x, y)))

        # single result: return it
        if len(results) == 1:
            return self.id_strs[results[0]]

        # multiple results: check which polygon contains the point
        else:
            point = shapely.geometry.Point(x, y)
            for id_ in results:
                id_str = self.id_strs[id_]
                shape = self.shapes[id_str]
                if shape.contains(point):
                    return id_str


def parse_geotag_input(geotag_input: str) -> GeotagInput:
    """ Parse a geotag operation instructions string.

    Instruction strings should have the following format:

        input_file$input_column>output_column

    Args:
        geotag_input: A string containing instructions for a geotag operation.

    Returns: A GeotagInput containing the parsed geotag_input.
    """

    (input_file, other_fields) = geotag_input.split("$")
    (input_column, output_column) = other_fields.split(">")

    return GeotagInput(
        input_file.strip(), input_column.strip(), output_column.strip()
    )


def dummy_function(*_, **__) -> None:
    """ A function that does nothing. """
    pass


if __name__ == "__main__":
    import argparse
    import copy

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "geotag", nargs="+",
        help="A list of geotag operation instructions. These should be in the format \"input_file$input_column>output_column\". This is passed directly into geopandas.read_file(). To directly read compressed shapefile archives, use \"zip://path/to/shapefile.zip\". NOTE: Be careful about bash! Be sure to enclose in single quotes (lack of quotes will register \">\" as an output redirect; double quotes will still register \"$\" as a variable prefix."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="The path to the input file. This is passed directly into pandas.read_csv()."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="The path to the output file. This is passed directly into pandas.core.frame.DataFrame.write_csv(). For compression, simply append \".gz\", etc."
    )
    parser.add_argument(
        "-l", "--longitude", required=True, metavar="LONGITUDE_FIELD",
        help="The name of the field containing longitude data in the input file."
    )
    parser.add_argument(
        "-L", "--latitude", required=True, metavar="LATITUDE_FIELD",
        help="The name of the field containing latitude data in the input file."
    )
    parser.add_argument(
        "-s", "--subset", metavar="SUBSET_COLUMNS",
        help="Optional. Mutually exclusive with -r/--rownames-only. A comma-separated list of fields to subset the input file to. This is passed to pandas.read_csv(), so this can be useful for limiting memory usage."
    )
    parser.add_argument(
        "-r", "--rownames-only", action="store_true", default=False,
        help="Optional. Mutually exclusive with -s/--subset. Creates a new column containing the data frame rownames and drops all other columns from the input. Rownames will be R-style, starting at 1, rather than pandas-style, starting at 0."
    )
    parser.add_argument(
        "-f", "--force-overwrite", action="store_true", default=False,
        help="Optional. Allows overwriting of existing columns in the data frame post-subsetting. Will not allow for overwriting of the longitude or latitude columns."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", default=False,
        help="Optional. Causes geotag.py to print out everything it is doing (otherwise the script will run without any output)."
    )

    args = parser.parse_args()

    # verbosity settings
    display = print
    if not args.verbose:
        display = dummy_function

    # input validation
    geotag_inputs = []
    for geotag_input in args.geotag:
        try:
            geotag_inputs.append(parse_geotag_input(geotag_input))
        except:
            raise Exception("could not parse geotag instructions {}".format(geotag_input))

    original_subset_columns = None
    subset_columns = None
    if args.subset:
        original_subset_columns = args.subset.split(",")
        subset_columns = copy.copy(original_subset_columns)

        # need to force inclusion of longitude and latitude columns; remove later
        if args.longitude not in subset_columns:
            display(
                "Forcing inclusion of \"{}\" pandas.read_csv (will remove later)"
                .format(args.longitude)
            )
            subset_columns.append(args.longitude)
        if args.latitude not in subset_columns:
            display(
                "Forcing inclusion of \"{}\" in pandas.read_csv (will remove later)"
                .format(args.latitude)
            )
            subset_columns.append(args.latitude)

    if args.subset and args.rownames_only:
        raise Exception("-s/--subset and -r/--rownames-only are mutually-exclusive")

    # check for duplicate output columns
    output_columns = [
        geotag_input.output_column
        for geotag_input in geotag_inputs
    ]
    if len(output_columns) != len(set(output_columns)):
        for column in set(output_columns):
            # inefficient, but probably not an issue due to small size
            output_columns.remove(column)
        raise Exception(
            "The following output columns are duplicated: {}"
            .format(", ".join(output_columns))
        )
    if args.longitude in output_columns:
        raise Exception(
            "Longitude column \"{}\" overlaps with an output column"
            .format(args.longitude)
        )
    if args.latitude in output_columns:
        raise Exception(
            "Latitude column \"{}\" overlaps with an output column"
            .format(args.latitude)
        )

    # read input
    if args.subset:
        display(
            "Reading input file (subsetting to {} as per -s/--subset): {}"
            .format(subset_columns, args.input)
        )
        df = pandas.read_csv(args.input, usecols=subset_columns)
    elif args.rownames_only:
        required_columns = [args.longitude, args.latitude]
        display(
            "Reading input file (subsetting to {} as per -r/--rownames-only): {}"
                .format(required_columns, args.input)
        )
        df = pandas.read_csv(args.input, usecols=required_columns)
    else:
        display("Reading input file: {}".format(args.input))
        df = pandas.read_csv(args.input)

    # drop null coordinates
    n_original_rows = len(df)
    df = df.dropna(subset=[args.longitude, args.latitude])
    n_dropped_rows = n_original_rows - len(df)
    display(
        "Dropped {}/{} columns with missing coordinates ({:0.2f}%)".format(
            n_dropped_rows, n_original_rows, n_dropped_rows/n_original_rows*100
        )
    )

    # generate rownames
    if args.rownames_only:
        display("Generating rownames")
        df["rowname"] = pandas.Series(df.index).apply(lambda rowname: str(rowname + 1)).astype(str)


    # check for overwriting of existing columns
    overwritten = set(df.columns).intersection(set(output_columns))
    if (len(overwritten) > 0) and (not args.force_overwrite):
        raise Exception(
            "Refusing to overwrite the following input columns (use -f to force): {}"
            .format(list(overwritten))
        )

    for i, geotag_input in enumerate(geotag_inputs):
        display(
            "\n({}/{}) Operation: {}${}>{}".format(
                i+1, len(geotag_inputs),
                geotag_input.input_file,
                geotag_input.input_column,
                geotag_input.output_column
            )
        )

        display("Reading: {}".format(geotag_input.input_file))
        gdf = geopandas.read_file(geotag_input.input_file)

        display("Checking validity of geometries")
        gdf["_invalid_geometry"] = gdf.geometry.progress_apply(
            lambda geometry: not geometry.is_valid
        )

        if gdf["_invalid_geometry"].any():
            display(
                "Correcting {} invalid geometries via geometry.buffer(0)"
                .format(sum(gdf["_invalid_geometry"]))
            )
            gdf.geometry = gdf.progress_apply(
                lambda row:
                    row.geometry.buffer(0)
                    if row["_invalid_geometry"]
                    else row.geometry,
                axis=1
            )

        tagger = Geotagger(gdf, geotag_input.input_column, verbose=args.verbose)

        display(
            "Geotagging: creating column \"{}\" <- {}${}"
            .format(
                geotag_input.output_column,
                geotag_input.input_file,
                geotag_input.input_column
            )
        )

        # TODO: is there a more elegant way to do this?
        if args.verbose:
            df[geotag_input.output_column] = df[[args.longitude, args.latitude]]\
                .progress_apply(
                    lambda xy: tagger.lookup(*xy),
                    axis=1
                )
        else:
            df[geotag_input.output_column] = df[[args.longitude, args.latitude]] \
                .apply(
                lambda xy: tagger.lookup(*xy),
                axis=1
            )
    display("")

    if args.rownames_only:
        extra_columns = [args.longitude, args.latitude]
        display("Dropping {} (as per -r/--rownames-only)".format(extra_columns))
        df = df.drop(extra_columns, axis=1)
    elif args.subset:
        extra_columns = list(set(subset_columns) - set(original_subset_columns))
        display("Dropping extra columns: {} (as per -s/--subset)".format(extra_columns))
        df = df.drop(extra_columns, axis=1)

    display("Writing to: {}".format(args.output))
    df.to_csv(args.output, index=False)

# geotag.py

Geocode flat CSV data using various shapefiles

![screenshot](screenshot.png)

### About

`geotag.py` implements the Geotagger class which can do very fast geotagging of big data in real time with live feedback through the use of R-tree spatial indices (via [rtree](https://toblerity.org/rtree/)) and [tqdm](https://tqdm.github.io/) (a progress bar library).

This functionality can be useful by allowing big data to be geotagged prior to the analysis process to create crosswalks that can then later be joined, rather than having to do a costly spatial join in the middle of the analysis process.

Other useful features include:

* Automatically validates input geometries and corrects invalid geometries via Shapely's `buffer(0)` routine (more info [here](https://shapely.readthedocs.io/en/stable/manual.html#constructive-methods)).
  * Note that it is still preferred to do this manually, if possible, as this may lead to [unexpected results](https://github.com/Toblerity/Shapely/issues/462#issuecomment-277179669)
* Automatically drops null coordinates.
* Supports globbed inputs, e.g. `tl_2020*bg.shp$GEOID>geoid_bg` - all matched paths will be separately read and concatenated into a single GeoDataFrame before proceeding.

### Usage

`geotag.py` can be used either embedded in a script or as a command-line application.

For information on command-line usage, see `geotag.py --help`. An example can be seen below using the example data found in this repository:

```bash
$ python3 geotag.py \
>     --input example.csv \
>     --longitude lon \
>     --latitude lat \
>     --output example-geotagged.csv \
>     --rownames-only \
>     --verbose \
>     'zip://tl_2010_25_zcta510.zip$GEOID10>zcta_10'
Reading input file: example.csv
Generating rownames and dropping other columns

(1/1) Operation: zip://tl_2010_25_zcta510.zip$GEOID10>zcta_10
Reading: zip://tl_2010_25_zcta510.zip
Creating rtree index: 100%|███████████| 538/538 [00:00<00:00, 9753.82 indexed/s]
Geotagging: creating column "zcta_10" <- zip://tl_2010_25_zcta510.zip$GEOID10
100%|████████████████████████████████████████████| 2/2 [00:00<00:00, 640.30it/s]

Dropping ['lon', 'lat'] (as per -r/--rownames-only)
Writing to: example-geotagged.csv
```

Embedding `geotag.py` is also straightforward:

```python
>>> import pandas, geotag, geopandas
>>> df = pandas.read_csv("example.csv")
>>> gdf = geopandas.read_file("zip://tl_2010_25_zcta510.zip")
>>> tagger = Geotagger(gdf, "GEOID10")
>>> df["zcta_10"] = df[["lon", "lat"]].apply(
...     lambda xy: tagger.lookup(*xy), axis = 1
... )
>>> df.to_csv("example-geotagged.csv", index=False)
```

### Limitations / future directions

Specifying the coordinate reference system (CRS) of the input CSV file is currently not supported; inputs are assumed to be in WGS 84 (EPSG:4326). Support for alternative CRSs would be valuable e.g. for survey points which may use a different CRS.

CRS transformations are not currently supported. Manual and automatic transformations may be added in the future, though this may be beyond the scope of this project.

The routines in the Geotagger object could theoretically be parallelized through a MapReduce-like programming pattern for faster geotagging, but the memory usage would increase proportinally with the number of cores being used. For this reason, parallelization has not yet been implemented.

Currently, only CSV input is supported. Support for other geospatial data files e.g. shapefiles or GeoJSONs as inputs may be added in the future, though this could complicate the Geotagger class

In the future, the current libraries may also be replaced with more performant alternatives, or CPU-intensive parts reimplemented in a faster language like C++.

Contributions welcome!
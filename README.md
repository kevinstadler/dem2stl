# dem2stl

`dem2stl` provides functions to convert digital elevation map (dem) files to `.stl` (and other) 3d-printable model files. It is different from other packages in that:

- it supports the generation of arbitrary (i.e. non-rectangular) shapes, such as countries or continents
- the generated models exhibit the correct curvature of the respective landmass on the Earth spheroid (this feature is only perceptible for larger geographic areas)

## Installation

    devtools::install_github("kevinstadler/dem2stl")

`dem2stl` relies on R's `rgl` package which requires OpenGL libraries to be installed. On Fedora, these are available through the `mesa-libGLU-devel` package

## Usage

```bash
wget https://research.cip.cgiar.org/gis/downloads/dataserver/_msk_alt/NLD_msk_alt.zip
unzip NLD_msk_alt.zip
```

```R
nl <- raster::raster("NLD_msk_alt.grd")

# aggregate to 2x2 cells because otherwise West-Friesland is clipped
# (the Noordzeekanal, marked with NA values, is so wide that it would split the
# Netherlands into two separated landmasses, with only the biggest rendered)
nl <- raster::aggregate(nl, 2, fun=max)

# render to 10x10cm max, with a thin base
m <- dem2stl::dem2mesh(nl, size=100, thicknessratio=0.001)
dem2stl::mesh2stl("nederland.stl", m)
```

Full documentation is available at <http://kevinstadler.github.io/dem2stl/>
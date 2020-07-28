desmon
====

Functions for design and monitoring of clinical trials

Performs various calculations related to design of studies with either failure
time or binary endpoints, and also functions to facilitate monitoring and
analysis of such studies, including group sequential studies with failure time
endpoints and two-stage and randomized phase II studies.

Installing
====

The desmon R package (current and older builds) can be found on the BCB Unix
servers under `~gray/R` and `~rredd/r.packages/desmon`

### On any platform

In R,
```r
# install.packages('devtools')
devtools::install_github('raredd/desmon', build_vignettes = TRUE)
```

### Windows

In R,
```r
install.packages('path/to/desmon.zip', repos = NULL, type = 'source')
```

### OSX

#### R 2.8

If you have version 2.8 of R installed, make sure that you also have the
`gfortran.pkg` installed (it can be found on the CRAN web site or in the
original disk image of the R 2.8 download from CRAN, in the folder labeled
"Packages"). Once the gfortran library is installed, the only remaining step
is to install the desmon library. In R, choose "Package Installer" from the
"Packages & Data" menu and under "Package Respository" choose "Local Binary
Package". Click Install All and navigate to the `desmon_0.9-12.tgz` file.
(Note: if you are installing the source file `desmon_0.9-12.tar.gz` then
choose "Local Source Package" instead). You probably want to install the
package at the default "System Level (in R framework)". That will put the
desmon package in `/Library/Frameworks/R.framework/Versions/2.8/Resources/library/`
and after this last step you should be all set.

#### R 2.9

If you have version 2.9 of R installed, make sure that you also have the
`gfortran.pkg` installed (it can be found on the CRAN web site or in the
original disk image of the R 2.9 download from CRAN, in the folder labeled
"Packages"). Once the gfortran library is installed, you need to copy this
directory: `/Library/Frameworks/R.framework/Versions/2.9` to this one:
`/Library/Frameworks/R.framework/Versions/2.8` (just use the Finder)
because the desmon library expects a "2.8" path to `libgfortran.2.dylib`
and `libR.dylib` (and maybe other libraries as well under this tree).

Now you can install the desmon library and everything should be set up so it
can find the files it needs. In R, choose "Package Installer" from the
"Packages & Data" menu and under "Package Respository" choose "Local Binary
Package". Click Install All and navigate to the `desmon_0.9-12.tgz` file.
(Note: if you are installing the source file `desmon_0.9-12.tar.gz` then
choose "Local Source Package" instead). You probably want to install the
package at the default "System Level (in R framework)". That will put the
desmon package in
`/Library/Frameworks/R.framework/Versions/2.9/Resources/library/` and now you
should be all set.

#### R > 3.0

Make a copy of your current `/Library/Frameworks/R.framework/Versions/x.x`
folder in the *same* directory and rename it to `.../2.8` as described above.

In the terminal:

```sh
which gcc
which gfortran
```
should return paths to `gcc` and `gfortran`, respectively.

To resolve the error
`make: gfortran-4.8: No such file or directory`

In the terminal:

```sh
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

Note that there is a copy of `gfortran` in the `/inst` folder of this package.

Loading the library
-------------------

Under either version of R, you can load the library and check out its contents
with the following two commands:

```r
library('desmon')
ls('package:desmon')
```

Vignettes
-------------------

Depending on the build options, binary and source packages may include the
vignettes. They may be viewed or accessed with the following commands:

```r
## list all
vignette(package = 'desmon')

## view topic
vignette('desmon-intro', package = 'desmon')
```

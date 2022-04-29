desmon
====

Functions for design and monitoring of clinical trials

Performs various calculations related to design of studies with either failure
time or binary endpoints, and also functions to facilitate monitoring and
analysis of such studies, including group sequential studies with failure time
endpoints and two-stage and randomized phase II studies.

# Installing

## via binaries

The `desmon` R package (current and older binaries) can be found on the DS Unix
servers under `~gray/R` and `~rredd/r.packages/desmon`

## from source

In R,

```r
# install.packages('devtools')
devtools::install_github('raredd/desmon', build_vignettes = TRUE)
```

Note that building the package from source may require some additional
build tools to compile FORTRAN code. If the above throws errors regarding
compilation, try the steps below for the relevant platform.

### Windows

Users will require [RTools](https://cran.r-project.org/bin/windows/Rtools)
for building R packages from source on Windows. This includes GNU
[gcc](https://gcc.gnu.org) as well as a full build system and package manager
to build, install, and distribute external c/c++/fortran libraries needed by
R packages ([source](https://cran.r-project.org/bin/windows/Rtools).

### macOS

Users will need to follow [instructions](https://mac.r-project.org/tools) to
install tools required for package compilation. These include Xcode command
line tools provided by [Apple](https://developer.apple.com/xcode/resources)
as well as a GNU FORTRAN compiler as this is not included with Xcode.

Running the following lines in a terminal should locate each program in the
user's path:

```sh
which gcc
which gfortran
```

# Loading the library

Under either version of R, you can load the library and check out its contents
with the following two commands:

```r
library('desmon')
ls('package:desmon')
```

# Vignettes

Depending on the build options, binary and source packages may include
vignettes. These may be viewed or accessed with the following commands:

```r
## list all
vignette(package = 'desmon')

## view topic
vignette('desmon-intro', package = 'desmon')
```

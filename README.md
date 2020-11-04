# adcirc-dgswem-coupling

## Table of contents
* [Introduction](#introduction)
* [Status](#status)
* [Setup](#setup)
    * [Dependencies](#dependencies)
    * [Compilation](#compilation)
    * [Usage](#usage)
* [Methodology](#methodology)
* [License](#license)

## Introduction
This is a Fortran software for weakly coupling [ADCIRC](https://github.com/adcirc/adcirc-cg) and [DG-SWEM](https://github.com/foci/dgswem). You need access to the GitHub repositories of both ADCIRC and DG-SWEM in order to use this software.

## Status
This project is *in progress*.

## Setup
### Dependencies
- ADCIRC - [link to website](http://adcirc.org/), [link to GitHub repository](https://github.com/adcirc/adcirc-cg)
- DG-SWEM - [link to website](https://users.oden.utexas.edu/~michoski/dgswem_doc/index.html), [link to GitHub repository](https://github.com/foci/dgswem)

### Compilation
To be determined.

### Usage
To be determined.

## Methodology
We plan to couple ADCIRC and DG-SWEM over models that are non-overlapping, but may other wise share a boundary. For now, Gauss-Seidel-style weak coupling will be implemented, with the possibility of exploring Gauss-Jacobi-style weak coupling in the future.

## License
This software is licensed under the terms of BSD 3-Clause "New" or "Revised" License; see [the LICENSE file](LICENSE). Note that some or all of the dependencies of this software are **not** open source, and these should not be distributed with this software.

## Resubmission
This is a resubmission. In this version we have:

* Created executable examples for all exported functions demonstrating their use.
    - Created a 'geco_small_result' dataset containing a pre-run record linkage result object for use in example code.
* Changed use of T or F to TRUE or FALSE.
* Removed use of unexported functions from examples.
* Bumped version number from 1.0.0 to 1.0.2

We have not included a reference describing the methods in this package but we plan to do this in a future update. This reference is in review and is not yet available online. When it is available we will update the package with the reference.

## R CMD check results
Checks performed using rhub::check_for_cran() and devtools::check_win_devel().

There were no ERRORs or WARNINGs on any platform.

There were 2 NOTEs on various platforms:

```
* checking CRAN incoming feasibility ... NOTE

New submission
Maintainer: 'Ian Taylor <ian.taylor@colostate.edu>'
```

This note appeared on all platforms. This is my first submission to CRAN.

```
* checking HTML version of manual ... NOTE
Skipping checking HTML validation: no command 'tidy' found
```

This note appeared on Fedora devel only.

## Downstream dependencies
This package has no downstream dependencies.

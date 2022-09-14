# TBProfiler-to-GARC
Covert the TBProfiler's TBDB to GARC for use with `piezo`


## Assumptions
Within the `tbdb.other_annotations.csv` there are some rows which appear to use an undocumented format for `iHGVS`. These rows contain several mutations, separated by the `|` character (which is reserved within `iHGVS` for denoting DNA methylation).

At first this seems similar to the `multi-mutation`'s use of `&` introduced in `GARC` to deal with splitting lines of allelic calls. However, due to the usual usage of the `|` character to denote a logical `OR`, it is assumed that these lines are denoting several similar mutations which share the rest of the line.
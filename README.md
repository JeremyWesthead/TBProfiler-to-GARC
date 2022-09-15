# TBProfiler-to-GARC
Covert the TBProfiler's TBDB to GARC for use with `piezo`. Some default rules for resistance genes are added for completeness.

Due to the nature of this catalogue, it is highly specific, and as such is significantly larger than other GARC catalogues.


## Assumptions
Within the `tbdb.other_annotations.csv` there are some rows which appear to use an undocumented format for `iHGVS`. These rows contain several mutations, separated by the `|` character (which is reserved within `iHGVS` for denoting DNA methylation).

At first this seems similar to the `multi-mutation`'s use of `&` introduced in `GARC` to deal with splitting lines of allelic calls. However, due to the usual usage of the `|` character to denote a logical `OR`, it is assumed that these lines are denoting several similar mutations which share the rest of the line.

## Questionable lines
Some rows of the `tbdb` have questionable mutations listed:
| Row number | Row | Reason |
| ---------- | --- | ------ |
| 69 | ethA,c.-1242_166del,ethionamide,resistance,,https://www.who.int/publications/i/item/9789240028173,Assoc w R - interim | This implies a deletion of 1409 bases (far beyond clockwork limits) + not in WHO |
| 46 | ethA,c.-1058_968del,ethionamide,resistance,,https://www.who.int/publications/i/item/9789240028173,Assoc w R - interim | This implies a deletion of 2027 bases + not in WHO |
| 118 | ethA,c.-2501_1459del,ethionamide,resistance,,https://www.who.int/publications/i/item/9789240028173,Assoc w R - interim | This implies a deletion of 3961 bases + not in WHO |

There are several similar cases. In the case of row 118, this deletes the entire ethA gene (assuming it has a promoter of 2501 bases) - and as shown in the table, does not infact exist within the WHO catalogue.
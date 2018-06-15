
## Clinical data download

2018-06-15

-----

This simple script downloads the pan-cancer clinical data resource from
Supplementary Table 1 of [this
paper](https://www.cell.com/cell/fulltext/S0092-8674\(18\)30229-0).

``` r
library(here)
download.file("https://www.cell.com/cms/attachment/2119196605/2090459781/mmc1.xlsx", here("data/tcga-liu2018_clinData.xlsx"))
```

This piece now reads the first sheet which contains the 11,160 clinical
records and exports to an RData (RDS)
object.

``` r
dat <- readxl::read_xlsx(here("data", "tcga-liu2018_clinData.xlsx"), na=c("", "#N/A"), guess_max=2000)
dat <- dat[, -1] # First column is just row number
saveRDS(dat, file=here("data", "tcga-liu2018_clinData.rds"))
readr::write_csv(dat, here("data", "tcga-liu2018_clinData.csv"))
```

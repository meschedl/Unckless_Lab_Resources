# qPCR-check

libraries

``` r
library(ggplot2)
```

Load in data

``` r
df <- read.csv("/Users/maggieschedl/Desktop/Github/Unckless_Lab_Resources/qPCR_analysis/well-check/qPCR-well-check.csv")
```

plot count of well by primer?

``` r
ggplot(df, aes(x = well, fill = primer)) +
  geom_bar() + facet_wrap(~primer) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](qPCR-well-check_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
ggplot(df, aes(x = well, fill = data.from)) +
  geom_bar() + facet_wrap(~primer) + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
```

![](qPCR-well-check_files/figure-commonmark/unnamed-chunk-4-1.png)

I don’t see any pattern. I think the only reason some wells are missing
is because I never put a control in those wells yet.

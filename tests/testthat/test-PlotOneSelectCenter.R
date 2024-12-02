library(RegionalST)

data("example_sce")
example_sce <- mySpatialPreprocess(example_sce, platform="Visium")
example_sce$"array_col" <- example_sce$col
example_sce$"array_row" <- example_sce$row
example_sce$"pxl_col_in_fullres" <- example_sce$imagecol
example_sce$"pxl_row_in_fullres" <- example_sce$imagerow
PlotOneSelectedCenter(example_sce, ploti = 1)

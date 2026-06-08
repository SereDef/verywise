# renv::snapshot()

# tools::package_dependencies("scales", reverse = TRUE, which = "Imports")

dep = "scales"

pkg_deps <- renv::dependencies()  # what your project uses
all_deps <- tools::package_dependencies(pkg_deps$Package, recursive = TRUE, which = c("Imports", "Depends"))

# Which package in your dependency tree imports scales?
names(Filter(function(x) dep %in% x, all_deps))
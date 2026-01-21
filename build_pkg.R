# 1. Generate documentation -----
# This creates the .Rd help files (man pages) from Roxygen comments and updates the NAMESPACE file.
devtools::document()

# 2. Build the website -----
# This uses the _pkgdown.yml to build the documentation site.
pkgdown::build_site()


# 3. Install the package -----
# This installs the latest version of the package into the local library to use.
devtools::install(build_vignettes = TRUE)


# 4. Full check -----
# Make sure everything is CRAN-ready, run a full check.
devtools::check()


# 5. Push to GitHub

# """
# git add .
# git commit -m “Added CACTI-S module and updated documentation”
# git push origin main
# """

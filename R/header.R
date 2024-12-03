.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

   packageStartupMessage("\n## BASIC SPACE SCALING PACKAGE")
   packageStartupMessage("## 2009 - ", this.year)
   packageStartupMessage("## Keith Poole, Howard Rosenthal, Jeffrey Lewis, James Lo, Royce Carroll, and Christopher Hare")
   packageStartupMessage("## Support provided by the U.S. National Science Foundation")
   packageStartupMessage("## NSF Grant SES-0611974\n")

   Sys.setenv("_R_BUILD_COMPACT_VIGNETTES_"="qpdf")

}

.onUnload <- function(libpath) {
    library.dynam.unload("basicspace", libpath)
}

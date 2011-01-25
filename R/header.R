.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

   cat("\n## BASIC SPACE SCALING PACKAGE \n")
   cat("## Copyright 2009 -", this.year, "\n")
   cat("## Keith Poole, Howard Rosenthal, Jeffrey Lewis, James Lo, and Royce Carroll\n")
   cat("## Support provided by the U.S. National Science Foundation\n")
   cat("## NSF Grant SES-0611974\n\n")

}

.onUnload <- function(libpath) {
    library.dynam.unload("basicspace", libpath)
}

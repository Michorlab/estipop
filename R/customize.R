###################################################################################----
#'  Customization code for time-dependent rates
#'  taken with minor adaptions from:
#'
#'  McDonald, T. O., & Michor, F. (2017).
#'  SIApopr: a computational method to simulate evolutionary branching trees for
#'  analysis of tumor clonal evolution. Bioinformatics, 33(14), 2221-2223.
#'
#'  Original code author: T. O. McDonald (mcdonald(at)jimmy.harvard.edu)
#'


##------------------------------------------------------------------------
#' create_timedep_template
#'
#' Creates a time-dependent rate ".cpp" template file at the given location
#' which can be modified before running `compile_timedep`.
#'
#' @param cppfile file name of .cpp file to create
#'
#' @export
#' @examples
#' \dontrun{
#' create_timedep_template(cppfile = "custom_rate_plugin.cpp")
#' }
create_timedep_template <- function(cppfile = "custom_rate_plugin.cpp"){
  cppfile <- unlist(strsplit(cppfile, ".", fixed = T))
  if(cppfile[length(cppfile)] != "cpp")
  {
    hfile <- c(cppfile, "h")
    cppfile <- c(cppfile, "cpp")
  }
  else
  {
    hfile <- c(cppfile[1:(length(cppfile)-1)], "h")
  }

  cppfile <- paste(cppfile, collapse = ".")
  hfile <- paste(hfile, collapse = ".")


  cpp_location <- paste(.libPaths()[1], "/estipop/extras/timedependent_template.cpp", sep = "")
  h_location <- paste(.libPaths()[1], "/estipop/extras/timedependent_template.h", sep = "")

  file.copy(cpp_location, cppfile)
  file.copy(h_location, hfile)

  # Update cpp file's include statement
  cpp_con <- file(cppfile)
  cpp_lines <- readLines(cpp_con)
  include_header <- paste("#include \"", hfile, "\"", sep = "")
  cpp_lines <- c(include_header, cpp_lines)
  write(cpp_lines, cppfile)
  close(cpp_con)
}



##------------------------------------------------------------------------
#' compile_timedep
#'
#' Allows the user to create a time-dependent rate function in C++ for use in
#' ESTIpop. Creates the appropriate header file and compiles the source file then
#' links into a shared object. The created shared object is used as a plugin.
#'
#' @param cppfile a character string giving the path name of a cpp file
#'
#' @export
#' @examples
#' \dontrun{
#' compile_timedep(cppfile = "./plugin.cpp")
#' }
compile_timedep <- function(cppfile){
  cpproot <- .pop_off(cppfile, ".", fixed = T)
  cppsuffix <- .pop(cppfile, ".", fixed = T)

  concopy <- paste(cppfile, ".backup", sep = "")
  file.copy(cppfile, concopy)
  if(.Platform$OS.type == "windows")
  {
    file.copy(paste(.libPaths()[1], "/estipop/extras/Makevars.win", sep = ""), "Makevars.win")

  }
  else if(.Platform$OS.type == "unix")
  {
    file.copy(paste(.libPaths()[1], "/estipop/extras/Makevars", sep = ""), "Makevars")
  }
  else
  {
    stop("Customization only support on Windows, OS X, and Linux systems")
  }

  # Need to make windows version
  shlib <- paste("R CMD SHLIB", cppfile)
  system(shlib)
}



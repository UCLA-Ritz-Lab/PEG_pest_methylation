meffil.report <- function (ewas.summary, output.file = "ewas-report.html", author = "Analyst", 
                           study = "Illumina methylation data", ...) 
{
  meffil:::msg("Writing report as html file to", output.file)
  report.path <- system.file("reports", package = "meffil")
  require(knitr)
  require(gridExtra)
  opts <- opts_chunk$get()
  on.exit(opts_chunk$set(opts))
  opts_chunk$set(warning = FALSE, echo = FALSE, message = FALSE, 
                 results = "asis", fig.width = 6, fig.height = 6)
  knit.report(file.path(report.path, "ewas-report.rmd"), output.file, 
              ...)
}

knit.report <- function (input.filename, output.filename, ...) 
{
  input.filename <- normalizePath(input.filename)
  output.dir <- dirname(output.filename)
  if (!file.exists(output.dir)) 
    dir.create(output.dir, recursive = T)
  output.dir <- normalizePath(output.dir)
  output.filename <- file.path(output.dir, basename(output.filename))
  name <- gsub("\\.[^.]+$", "", basename(output.filename))
  suffix <- gsub(".*\\.([^.]+)$", "\\1", output.filename)
  is.html <- tolower(suffix) %in% c("htm", "html")
  if (is.html) 
    md.filename <- file.path(output.dir, paste(name, "md", 
                                               sep = "."))
  else md.filename <- output.filename
  current.dir <- opts_knit$get("output.dir")
  on.exit(opts_knit$set(output.dir = current.dir))
  opts_knit$set(output.dir = output.dir)
  cwd <- getwd()
  on.exit(setwd(cwd), add = TRUE)
  if (is.null(options("knitr.in.progress")[[1]])) {
    setwd(output.dir)
    knit(input.filename, output = md.filename, envir = parent.frame(), 
         ...)
  }
  else {
    opts_knit$set(progress = FALSE)
    out <- knit_child(input.filename, envir = parent.frame(), 
                      quiet = T)
    writeLines(out, con = md.filename)
  }
  lines <- readLines(md.filename)
  lines <- gsub("![plot", "\n\n![plot", lines, fixed = T)
  writeLines(lines, md.filename)
  if (is.html) 
    markdownToHTML(md.filename, output.filename)
}
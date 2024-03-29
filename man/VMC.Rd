% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VMC.R
\name{VMC}
\alias{VMC}
\title{Vertex Map Correction}
\usage{
VMC(
  data,
  Tree_DBH = "DBH",
  X_col = "X.m.",
  Y_col = "Y.m.",
  Alt_col = "ALTITUDE",
  Base_ID = "Base",
  Tree_ID = "TreeID",
  limit = 10,
  ref_num = 2,
  plot = FALSE,
  z_corr = FALSE,
  randomize = FALSE,
  seed = 11111
)
}
\arguments{
\item{data}{data.frame including coordinates, Tree-IDs and "Base-number"}

\item{Base_ID}{name for base ID column}

\item{Tree_ID}{name for tree ID column}

\item{limit}{number of iterations without change before loop stops (default = 10)}

\item{ref_num}{minimum number of overlapping trees (default = 2)}

\item{plot}{logical, if true plots two maps for comparison (requires ggplot2 and ggrepel)}

\item{z_corr}{logical, if true correction of z-column will be included}

\item{randomize}{logical, if true randomizes order for selecting base-pairs}

\item{seed}{seed for random selection of base-pairs}

\item{TreeDBH}{name for DBH column (only needed for plotting)}

\item{X_Col}{name for X-coordinates column}

\item{Y_Col}{name for Y-coordinates column}

\item{Alt_Col}{name for Altitude column (only needed for Z correction)}
}
\value{
data.frame with corrected X and Y values
}
\description{
This function tries to minimize the offset introduced when mapping an area using a vertex laser geo 360.
The function checks for reference trees that where measured from several positions (bases) and simply
loops through base-pairs with overlapping information. The offset is reduced by moving all data-points
measured from one base towards the other base. Since the data cannot fit perfectly it repeats the step with
different base-pairs as long as the sum of the squares of all offsets is further reduced.
}

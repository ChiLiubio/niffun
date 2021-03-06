% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_niffun.R
\name{trans_niffun}
\alias{trans_niffun}
\title{Create an R6 object for functional prediction of nitrogen fixers based on nifH sequences.}
\description{
This class is a wrapper for a series of functional analysis on diazotrophic communities based on the nifH sequences.
}
\examples{

## ------------------------------------------------
## Method `trans_niffun$new`
## ------------------------------------------------

\donttest{
data(dataset_nifH)
t1 <- trans_niffun$new(dataset = dataset)
}

## ------------------------------------------------
## Method `trans_niffun$cal_blast`
## ------------------------------------------------

\donttest{
t1$cal_blast()
}

## ------------------------------------------------
## Method `trans_niffun$cal_pathway`
## ------------------------------------------------

\donttest{
t1$cal_pathway()
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{trans_niffun$new()}}
\item \href{#method-cal_blast}{\code{trans_niffun$cal_blast()}}
\item \href{#method-cal_pathway}{\code{trans_niffun$cal_pathway()}}
\item \href{#method-clone}{\code{trans_niffun$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create the trans_niffun object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_niffun$new(dataset = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dataset}}{default NULL; the object of microtable Class.}
}
\if{html}{\out{</div>}}
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
data(dataset_nifH)
t1 <- trans_niffun$new(dataset = dataset)
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cal_blast"></a>}}
\if{latex}{\out{\hypertarget{method-cal_blast}{}}}
\subsection{Method \code{cal_blast()}}{
Taxonomic assignment with blast.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_niffun$cal_blast(
  blast_tool_path = NULL,
  path_to_temp_folder = NULL,
  evalue = 1e-05,
  max_target_seqs = 50,
  num_threads = 1,
  other_blast_parameters = ""
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{blast_tool_path}}{default NULL; the folder path, e.g. ncbi-blast-2.5.0+/bin ; blast tools folder downloaded from 
"ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+"  ; e.g. ncbi-blast-2.5.0+-x64-win64.tar.gz  for windows system; 
if blast_tool_path is NULL, search the tools in the environmental path variable.}

\item{\code{path_to_temp_folder}}{default NULL; The temporary folder to store the logfile, intermediate file and result files; if NULL, 
use the default temporary in the computer.}

\item{\code{evalue}}{default 1e-5; the E value threshold used in blast.}

\item{\code{max_target_seqs}}{default 50; Maximum number of aligned sequences to keep.}

\item{\code{num_threads}}{default 1; Number of threads (CPUs) to use in the BLAST search.}

\item{\code{other_blast_parameters}}{default ""; other parameters provided to BLAST, such as "-perc_identity 80".}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_tax_table stored in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
t1$cal_blast()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-cal_pathway"></a>}}
\if{latex}{\out{\hypertarget{method-cal_pathway}{}}}
\subsection{Method \code{cal_pathway()}}{
Predict KEGG orthology and pathway abundance.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_niffun$cal_pathway(min_identity_to_reference = 95)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{min_identity_to_reference}}{default 95; the sequences identity threshold used for finding the nearest species genome.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_KEGG_KO and res_KEGG_pathway stored in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
t1$cal_pathway()
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_niffun$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}

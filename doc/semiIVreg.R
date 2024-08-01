## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>", echo = TRUE, fig.retina = 2#,#fig.height=4, fig.width=5, fig.align='center'
)

## ----setup--------------------------------------------------------------------
library(semiIVreg)

## ----installation-rforge, eval=FALSE------------------------------------------
#  remotes::install_github("cbruneelzupanc/semiIVreg")

## ----eval=FALSE---------------------------------------------------------------
#  semiivreg(y~d|w0|w1|x, data)

## ----mte, fig.height=3, fig.width=7, fig.align='center'-----------------------
library(semiIVreg)
data(roydata) # load the data from a simulated Roy model

# semi-IV regression
semiiv = semiivreg(y~d|w0|w1, data=roydata) 

## ----mtr, fig.height=4, fig.width=5, fig.align='center'-----------------------
semiiv$plot$mtr


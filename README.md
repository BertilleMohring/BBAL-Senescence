# Code for: Environmental variability shapes life-history trade-offs within and between populations of a long-lived seabird

Bertille Mohring, Jonathan R. Potts, Alastair J. Wilson, Denis Réale, Richard A. Phillips, Henri Weimerskirch, Christophe Barbraud, Ashley Bennison, Karine Delord, Andrew G. Wood, Samuel Peroteau, Etienne Rouby, Francesco Ventura, Samantha C. Patrick

## Overview

This repository contains the code and data used in the following manuscript: Mohring B, Potts JR, Wilson AJ, Réale D, Phillips RA, Weimerskirch H, Barbraud C, Bennison A, Delord , Wood AG, Peroteau S, Rouby E, Ventura F, Patrick SC (under review in Ecology Letters) Environmental variability shapes life-history trade-offs within and between populations of a long-lived seabird

## Data details

The file data\_BBAL.csv contains the data used in this manuscript. Details of the data are presented below:

* population : population (Bird Island or Kerguelen)
* ReproS : Reproductive success (0: failed breeding attempt; 1: successful breeding attempt)
* Age : Individual age
* Year\_last\_seen : breeding season during which the individual was seen for the last time
* First\_breeding\_attempt : categorical variable for primiparity/multiparity of the breeding attempt of a given individual ("yes": first breeding attempt of an individual; "no": not the first breeding attempt)
* YearPop : breeding seasons and population corresponding to the breeding attemps
* cohortPop : cohort (birth year) and population
* lifetime\_repro\_output : Total number of succesful breeding attempts of an individual during the monitoring period (i.e., during an individual's lifetime)
* nb\_breeding\_attempts : Total number of breeding attempts of an individual during the monitoring period (i.e., during an individual's lifetime)
* ID : Individual identity

## Code details

This repository contains the code used in the manuscript and supporting informations. The details of the different R scripts is presented below:

* *1\_GLMM\_senescence\_run\_and\_save\_model\_normalPrior.R*: Run the main GLMM looking at age-related variation in reproductive success in black-browed albatrosses.
* *2\_Extract\_parameters\_of\_interest\_from\_model\_output.R*: Extract the parameters of interest estimated by the GLMM (notably: age at onset of senescence, senescence rate, early-life probability of successful reproduction and probability of successful reproduction at the onset of senescence).
* *3\_Figure2.R*: Code used to produce Figure 2.
* *4\_Get\_variance.R*: Extract the variance in parameters of interest. This code is used to produce Table 2.
* *5\_Figure3.R*: Code used to produce Figure 3.
* *6\_link\_with\_fitness.R*: Code used to investigate the link between lifetime reproductive success, onset and rate of senescence.
* *SupportingInformation\_S1.R*: Code used to produce Supporting information S1.
* *SupportingInformation\_S3.R*: Code used to produce Supporting information S3.
* *SupportingInformation\_S4.R*: Code used to produce Supporting information S4.
* *SupportingInformation\_S5.R*: Code used to produce Supporting information S5.
* *SupportingInformation\_S6.R*: Code used to produce Supporting information S6.
* *SupportingInformation\_S7.R*: Code used to produce Supporting information S7.


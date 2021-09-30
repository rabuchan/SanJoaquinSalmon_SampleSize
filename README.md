# SanJoaquinSalmon_SampleSize
R code to perform sample size simulations for acoustic telemetry studies of juvenile salmonid migration survival through Sacramento-San Joaquin River Delta

Author: Rebecca Buchanan, School of Aquatic and Fishery Sciences, University of Washington/n
Source: Buchanan, R. 2015. Sample Size for 2016 Chinook Tagging Study:  Supplemental Release. Prepared for USFWS, Stockton, CA.

Objective: Determine size of release group necessary to be able to estimate:
(1) survival from release to Chipps Island
(2) route-specific survival to Chipps Island
(3) route selection at the head of Old River and at Turner Cut

Estimation criteria considered:
C1: Parameter is estimable in specified minimum percentage of simulations (e.g., 95%)
C2: Probability estimate is not greater than specified maximum >1 in specified minimum percentage of simulations (e.g., not > 1.1 in 95% of simulations)
C3: Standard error on probability estimate is not greater than specified value (e.g., 0.05)
C4: Bias in parameter estimate is less than specified maximum (e.g., 0.05)

Release sites:
Primary: Durham Ferry (DF)
Secondary: Between Lathrop and Turner Cut (e.g., Stockton) (STK)

Detection sites (telemetry stations):
San Joaquin River at Lathrop: A1 (dual array)
San Joaquin River at MacDonald Island: A2 (dual array)
Old River at head of Old River: B1 (dual array)
Turner Cut: F1 (dual array)
Chipps Island: G2 (dual array)

Parameters estimated:
sR = overall survival from Durham Ferry to Chipps Island
sRO = survival from Durham Ferry to head of Old River
sA = survival from head of Old River to Chipps Island in the San Joaquin River route
sB = survival from head of Old River to Chipps Island in Old River route
sA1 = survival from head of Old River to Turner Cut Junction in the San Joaquin River
sA2 = survival from Turner Cut Junction to Chipps Island in the San Joaquin River route
sF2 = survival from Turner Cut Junction to Chipps Island in the Turner Cut route
psiA1 = probability of selecting San Joaquin River route at the head of Old River
psiA2 = probability of selecting San Joaquin River route at Turner Cut Junction
pA1 = conditional detection probability at site A1 (Lathrop)
pA2 = conditional detection probability at site A2 (MacDonald Island)
pB1 = conditional detection probability at site B1 (Old River at its head)
pF1 = conditional detection probability at site F1 (Turner Cut)
pG2 = conditional detection probability at site G2 (Chipps Island)
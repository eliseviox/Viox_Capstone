---
title: "Viox_Elise_Capstone"
author: "Elise Viox"
date: "4/28/2020"
output: html_document
---

1) Provide a brief background and significance about a specific research problem that interests you. It could be project you’re involved with now, or a rotation project, or something you’d like to work on. The reader will need to understand enough background to make sense of the experiment you propose below. Keep it brief. In one short paragraph.

**The main barrier to HIV-1 cure is the persistence of a reservoir of latently infected CD4+ T cells harboring replication-competent proviruses. One promising strategy for eliminating HIV-1 infected cells following viral reactivation is through antibody-dependent cellular cytotoxicity (ADCC). HIV-specific antibodies with ADCC effector functions are elicited by the majority of HIV-infected individuals and preferentially target the “open” CD4-bound conformation of HIV-1 envelope glycoproteins (Env). However, due to the “closed” conformation of most primary HIV-1 isolate Envs, these antibodies are ineffective in ADCC-mediated clearance of infected cells. In order to utilize these naturally elicited, ADCC-mediating antibodies, small molecule CD4 mimetics (CD4mc) have been designed to bind to the closed conformation of HIV-1 Env, force the envelope into an “open” CD4-bound conformation, and expose epitopes that are targeted by ADCC-mediating antibodies. In this proposal, we aim to evaluate the potential of a small molecule CD4mc compound to reduce the HIV-1 viral reservoir in anti-retroviral therapy (ART)-treated, Simian/Human Immunodeficiency Virus (SHIV)-infected rhesus macaques (RMs).**

2) Briefly state something that is unknown about this system that can be discovered through, and leads to, an experiment.  For example, "It is not known whether....."
**The CD4mc compound BNM-III-170 targets the conserved HIV-1 gp120 Phe43 cavity on the surface of infected cells and has been shown to sensitize SHIV-infected cells to ADCC killing mediated by sera from HIV-1-infected individuals in vitro. However, the biological activity and virologic impact of BNMIII-170 remains to be demonstrated in in vivo.**

3) Make an “if” “then” prediction that is related to item #2. It should be of the general form, “if X is true, then Y should happen”.
**If SHIV-infected RMs are treated with the CD4mc BNM-III-170 at ART interuption, then the SHIV-reservoir (aka the frequency of latently infected CD4+ T cells) will be reduced in these animals.**

4) What dependent variable will be observed to test this prediction in item #3? What predictor variable will be used to manipulate the system experimentally? Define the inherent properties of these variables (eg, are they sorted, ordered or measured).
**The dependent variable is the frequency of latently infected CD4+ T cells and it is a discrete, measured outcome variable. The predictor variable that will be used to manipulate the system experimentally is the treatment groups: untreated (control) and CD4mc BNM-III-170 treatment. The predictor variable is a nominal, categorical variable.**

5) Write a statistical hypothesis.  There should be a null and alternate. These should be explicitly consistent with the prediction in item #3 and the response variable in #4. In other words, make sure the statistical hypotheses that you write here serves as a test of the prediction made in item #3.

**The magnitude of the SHIV reservoir will be measured by Tat/rev Induced Limiting Dilution Assay (TILDA). TILDA measures the frequency of CD4+ T cells with inducible multiply-spliced SHIV RNA (msRNA).**

**Null hypothesis: The size of the SHIV reservoir (measured by the frequency of CD4+ T cells with inducible msRNA) in SHIV-infected RMs treated with CD4mc BNM-III-170 will be greater than or equal to the the size of the SHIV reservoir (measured by the frequency of CD4+ T cells with inducible msRNA) in untreated control SHIV-infected RMs.**

$H0 = \mu_{CD4mcTreated}\geq\mu_{Untreated}$

**Alternative hypothesis: SHIV-infected RMs treated with CD4mc BNM-III-170 will have a smaller SHIV reservoir (lower frequency of CD4+ t cells with inducible msRNA ) than SHIV-infected RMs in the untreated control group.**

$Ha = \mu_{CD4mcTreated} < \mu_{Untreated}$

**The alternative hypothesis is based on previous studies that show that BNM-III-170 results in conformational changes in HIV-1 Env expressed on the surface of HIV-1 infected target cells and exposes key epitopes such as the cluster A epitope that are commonly targeted by non-neutralizing, ADCC-mediating antibodies commonly elicited during natural infection. We expect that these ADCC-mediating antibodies may recruit Fc-receptor bearing monocytes and NK cells to actively lyse infected target cells and ultimately lead to a reduction in the size of the latent reservoir.**

6) What is the statistical test you would use to test the hypothesis in item #5? Briefly defend what makes this appropriate for the hypothesis and the experimental variables. If there are alternatives, why is this approach chosen instead? Points will not be awarded if the justification involves something like "because everybody does it this way".
**I would use the Mann Whitney Rank Sum Test for 2 independent groups to test the hypothesis in item #5. I selected the Mann Whitney Rank Sum Test for 2 independent groups because it is used to compare two groups that receive either of 2 levels of a predictor value. In the case of my experiment, one group of 6 independent replicates (6 RMs) will be exposed to control (no treatment) and a second group of 6 independent replicates (6 RMs) will be exposed to treatment with the small molecule CD4 mimetic BNM-III-170.**

7) List the procedures and decision rules you have for executing and interpreting the experiment. These procedures range from selection of experimental units, to randomization to primary endpoint to threshold decisions. Define (and defend) what you believe will be the independent replicate.

**Each rhesus macaque will count as a biological replicate. Rhesus macaques are incredibly expensive to maintain and the CD4mc BNM-III-170 compound is also very expensive to produce. To determine the minimum number of rhesus macaques that are necessary to include in each group to give my experiment the severe test it deserves, I will conduct a power analysis (see below). Based on this power analysis, I need to include a minimum of 4 RMs per treatment group in my study (8 RMs total). However, given that some SHIV-infected RMs virally control naturally during infection and may do so before we get the chance to initiate ART (and would therefore be removed from the study), I will include 2 additional animals per group (6 RMs per treatment group, 12 RMs in total). All RMs (n=12) will be intravenously infected with 200 TCID50 of the pathogenic CCR5 (R5)-tropic SHIV molecular clone SHIV AD8-EO. Assuming no animals virally control before 6 weeks post infection (p.i.), all 12 RMs will be placed on an optimized daily ART regimen (Tenofovir, TDF; Emtricitabine, FTC; and Dolutegravir, DTG) at 6 weeks p.i. RMs will then be randomly assigned to a control group (n=6 if no animals virally control before ART) and a CD4mc treatment intervention group (n=6 if no animals virally control before ART).**

**Beginning at week 18p.i., the animals in the intervention group will receive a total of 15 36mg/kg doses of BNM-III-170 administered subcutaneously every three days. 24 hours after the last dose of BNM-III-170 is administered to RMs in the treatment group, blood will be collected from animals in the treated and untreated groups. CD4+ T cells will be isolated from each blood sample and the frequency of CD4+ T cells with inducible, multiply-spliced SHIV RNA will be measured via the Tat-inducible limiting dilution assay (TILDA).**

**For the power analysis, I am declaring that a 25% difference in the size of the latent SHIV reservoir between the treated and untreated groups would be scientifically meaningful based on previous HIV cure intervention studies. This 25% difference would translate to a delta value (aka difference of means between the two groups) of approximately 50 msRNA+ CD4+ cells/ 10^6 CD4+ cells assuming that the mean msRNA of untreated, SHIV-infected RMs on ART is 200 msRNA+ CD4+ cells/ 10^6 CD4+ cells. I will also be using 5% for type1 error and 20% for type2 error (80% power) as error tolerance thresholds.**

8) Produce a graph of a simulation for the expected results. Create a dataMaker-like function in R to create and plot the data. Label and scale any axis. The graph should illustrate the magnitude of the expected response, or the level of response that you expect to see and would be minimally scientifically relevant. Be sure to illustrate any variation that is expected.

```{r}
#Load libraries
library(tidyverse)
library(ez)
library(viridis)

#The magnitude of the SHIV reservoir will be measured by Tat/rev Induced Limiting Dilution Assay (TILDA). TILDA measures the frequency of CD4+ cells with inducible multiply-spliced HIV RNA (msRNA+ CD4+ cells/million CD4+ cells)

#Size of the reservoir in SHIV-infected, ART-treated RMs (based on existing literature)
msRNA_AtBaseline <- 200

#Percent reduction in the reservoir with CD4mc BNM-III-170 treatment
msRNA_PercentReduction_CD4mcTreated <- 0.75

#Expected standard deviation (based on commonly observed std in TILDA assays in SHIV-infected RMs studies)
std <- 25

#Number of rhesus macaques in each group (control untreated and CD4mc-treated)
replicates <- 6

sims <- 200

resultsimulation <- function(replicates, msRNA_AtBaseline, msRNA_PercentReduction_CD4mcTreated, std) {
 Untreated <- rnorm(replicates, msRNA_AtBaseline, std)
 CD4mcTreatment <- rnorm(replicates, (msRNA_AtBaseline*msRNA_PercentReduction_CD4mcTreated), std)
  outcome <- c(Untreated, CD4mcTreatment)
  predictor <- c(rep(c("Untreated", "CD4mctreatment"), each = replicates))
  ID <- as.factor(c(1:length(predictor)))
  dataframe <- data.frame(ID, predictor, outcome)
}

CD4mcStudy_Data <- resultsimulation(replicates, msRNA_AtBaseline, msRNA_PercentReduction_CD4mcTreated, std)
CD4mcStudy_Data

ggplot(CD4mcStudy_Data, aes(predictor, outcome, color=predictor)) +
  geom_boxplot(outlier.colour="black", outlier.size=2) +
  geom_jitter(width=0.1, size=2, alpha=0.5) +
  labs(x = "Treatment", y = "msRNA+ CD4+ T cells / 10^6 CD4+ T cells") +
  ggtitle("Effect of CD4mc BNM-III-170 on the Magnitude of the Latent SHIV Reservoir") +
  theme_classic() +
  theme(legend.position="none")
```

9) Write and perform a Monte Carlo analysis to calculate a sample size necessary to test the hypothesis. This Monte Carlo must test the primary endpoint.
```{r}
#Power T test to determine number of RMs needed per group 
##A 25% difference in the size of the latent SHIV reservoir between the treated and untreated groups would be scientifically meaningful based on previous HIV cure intervention studies
###This 25% difference would translate to a delta value (aka difference of means between the two groups) of approximately 50 msRNA+ CD4+ cells/ 10^6 CD4+ cells assuming that the mean msRNA of untreated, SHIV-infected RMs on ART is 200 msRNA+ CD4+ cells/ 10^6 CD4+ cells. 
###I will be using 5% for type1 error and 20% for type2 error (80% power) as error tolerance thresholds.
PowerTTest <- power.t.test(n = NULL, power = 0.8, delta = 50, sd = 25, sig.level = 0.05, 
             type = "two.sample", alternative = "one.sided") 
PowerTTest
##For Power T Test, n=3.987 RMs; Since we can't have 0.987 of an RM, we will assign 4 RMs per group (although we would likely assign more RMs in case any of them virally control prior to ART initiation)
RMpergroup <- 4


#Typical frequency of msRNA+ CD4+ T cells in SHIV-infected RMs on ART is 200 msRNA+ CD4+ T cells/ 10^6 CD4+ T cells
Untreated_mean_msRNA <- 200
#Expected standard deviation (based on commonly observed std in TILDA assays in SHIV-infected RMs studies)
Untreated_std <- 25

#A 25% difference in the size of the latent SHIV reservoir between the treated and untreated groups would be scientifically meaningful based on previous HIV cure intervention studies
###If the mean msRNA+ CD4+ T cells/10^6 CD4+ T cells in untreated RMs is 200 then a 25% difference would translate to  150 msRNA+ CD4+ cells/ 10^6 CD4+ cells for CD4mc treated, SHIV-infected RMs.
CD4mcTreated_mean_msRNA <- 150
#Expected standard deviation (based on commonly observed std in TILDA assays in SHIV-infected RMs studies)
CD4mcTreated_std <- 25

#I will be using 5% for type1 error as an error tolerance threshold.
alpha <- 0.05

simulations <- 200

montecarlo_data1 <-numeric(simulations) 

for(i in 1:sims){
  Untreated <- rnorm(n = RMpergroup, mean = Untreated_mean_msRNA, sd = Untreated_std)
  CD4mcTreated <- rnorm(n = RMpergroup, mean = CD4mcTreated_mean_msRNA, sd = CD4mcTreated_std)
  data_Ttest <- t.test(Untreated, CD4mcTreated, paired = FALSE, alternative="greater")
  montecarlo_data1[i]<-data_Ttest$p.value
}

montecarlo_data2 <- length(which(montecarlo_data1 < alpha))
montecarlo_data2

power <- montecarlo_data2/simulations
power

ggplot(data.frame(montecarlo_data1)) + 
  geom_histogram(aes(montecarlo_data1), color = "#002878") + 
  labs(x = "p-value")



```


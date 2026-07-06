# Combining ability, heritability and stability for grain yield and phosphorus-use efficiency in a maize line × tester population evaluated under contrasting phosphorus supply in Mozambique

Stélio Boaventura Nuvunga^1,2^, Fausto Rodríguez-Zapata^2^, Rubén Rellán-Álvarez^2^, Joseph L. Gage^3^

^1^Instituto de Investigação Agrária de Moçambique (IIAM), Maputo, Mozambique

^2^Department of Molecular and Structural Biochemistry, North Carolina State University, Raleigh, NC, USA

^3^Department of Crop and Soil Sciences, North Carolina State University, Raleigh, NC, USA

*Corresponding author:* [email] · [ORCID]

*(ORCIDs, corresponding author and Stélio's affiliation footnotes to be confirmed. Note: this
draft is formatted with the supplied MWEN/IMWA template styles; substitute the Theoretical and
Applied Genetics template before submission.)*

# Abstract

Low soil-phosphorus availability limits maize productivity across much of sub-Saharan Africa,
and breeding phosphorus-efficient hybrids requires knowing how much of the variation in grain
yield (GY) and phosphorus-use efficiency (PUE) is heritable and combinable. We evaluated a
7 × 7 line × tester (North Carolina Design II) population — 49 F₁ hybrids from seven CIMMYT
donor lines crossed with seven IIAM tester lines — across five environments in Mozambique
(four locations, with the Chokwe site grown at both optimal and phosphorus-stress levels),
in a randomized complete block design with three replications. Eleven traits were analysed:
GY and PUE (selection), harvest index, root:shoot ratio, root, shoot and total dry mass
(physiological), and plant height, days to anthesis, days to silking and the
anthesis–silking interval (auxiliary). Combined line × tester analysis of variance,
REML variance components (`sommer`), genetic parameters (additive and dominance variance,
broad- and narrow-sense heritability, Baker's predictability ratio), trait correlations and
AMMI/GGE stability were estimated. Male-tester GCA was effectively null for GY and its
components; narrow-sense heritability was low for GY (h² = 0.12) and PUE (h² = 0.03) but high
for flowering time (h² ≈ 0.60, BR ≈ 0.98). PUE and GY were driven by harvest index and shoot
biomass (R² = 0.81 and 0.96). The low additive variance for yield is attributed to the narrow,
yield-selected parental base. Indirect selection on harvest index and biomass, and exploitation
of specific combining ability, are proposed for improving PUE.

**Keywords** phosphorus-use efficiency · line × tester · general combining ability ·
heritability · genotype × environment · maize

# Introduction

Phosphorus (P) is among the most limiting nutrients for maize (*Zea mays* L.) production on
the weathered, P-fixing soils that dominate sub-Saharan Africa, where most smallholders cannot
afford corrective fertilization. Breeding hybrids that yield well under low-P supply — that is,
with high phosphorus-use efficiency (PUE) — is therefore a strategic objective for national
programmes such as Mozambique's Instituto de Investigação Agrária de Moçambique (IIAM).

Realizing genetic gain for a complex trait such as PUE depends on the amount and nature of the
genetic variance available: additive (general combining ability, GCA) variance is directly
selectable and transmissible, whereas non-additive (specific combining ability, SCA) variance
is exploited through specific hybrid combinations. The line × tester (North Carolina Design II)
mating design partitions among-cross variation into GCA of each parental pool and SCA, and,
when evaluated across environments, into the corresponding genotype × environment (G×E)
interactions.

Here we analyse a 7 × 7 line × tester maize population — CIMMYT low-P donor lines crossed with
IIAM testers — evaluated across five Mozambican environments that contrast in P supply. Our
objectives were to (i) estimate combining ability, variance components, heritability and gene
action for grain yield, PUE and nine physiological and phenological traits; (ii) identify the
traits that drive PUE and could serve as indirect-selection criteria; (iii) assess hybrid
stability across environments; and (iv) interpret the observed genetic architecture in light of
the parental germplasm base.

# Materials and Methods

## Plant material

The population is a 7 × 7 factorial (North Carolina Design II) of 49 F₁ hybrids, with no selfs
or reciprocals. The **female (line) pool** comprised seven CIMMYT tropical donor inbreds
selected for low-P tolerance: CML364, CML366, CML434, CML435, CML439, CML530 and CML532. The
**male (tester) pool** comprised seven IIAM tester lines, designated by their IIAM donor
identifiers EN17, EN20, EN21, EN25, EN31, EN32 and EN64 (corresponding to the legacy NCDII
identifiers MOZL3, MOZL5, MOZL6, MOZL7, MOZL8, MOZL9 and MOZL10). The testers derive from three
CIMMYT-based pedigree backgrounds — (ZM421 × CML491), (ZM421 × CLRCY034) and (MATUBASG × CML539)
— and were retained as the top seven for yield from a broader set of elite candidates. A
released commercial hybrid, NAMULONGUE, was included as a repeated check.

## Locations, phosphorus treatments and experimental design

Hybrids were evaluated across **four locations** in Mozambique. At Chokwe the trial was grown
under **two phosphorus regimes** — optimal (50 kg P ha⁻¹) and stress (10 kg P ha⁻¹) — while
Maniquenique was grown under optimal supply and Nhacoongo and Sussundenga under P-stress,
giving **five environment × treatment combinations** (hereafter environments): CHOKWE OPT,
CHOKWE STS, MANIQUENIQ OPT, NHACOONGO STS and SUSSUNDENGA STS. Because both P levels were
applied to the same genotypes only at Chokwe, the phosphorus effect was estimated within that
paired site to avoid confounding with location.

Each environment was laid out as a **randomized complete block design (RCBD)** with three
replications. Within each replication, the 50 entries (49 hybrids + check) were sown in a
5-column × 10-row **serpentine** arrangement; field row and column were recovered from the
planting order and fitted as nested spatial covariates. The complete dataset comprised **750
plots** (49 hybrids × 5 environments × 3 replications, plus the check), a balanced 7 × 7
factorial. Planting dates ranged from early July to early August; days to flowering were
computed relative to each environment's planting date.

## Traits

Eleven traits were analysed, in three groups. **Selection traits:** grain yield (GY, t ha⁻¹)
and phosphorus-use efficiency (PUE = grain yield / applied P). **Physiological traits:**
harvest index (HI = grain dry weight / total aboveground dry biomass), root:shoot ratio
(RSR = root dry weight / shoot dry weight), root weight (RW), shoot dry weight (SW) and total
dry mass (TDM = SW + RW). **Auxiliary traits:** plant height (PH, cm), days to anthesis (DTA),
days to silking (DTS) and the anthesis–silking interval (ASI = DTS − DTA).

## Statistical analyses

All analyses were performed in R. Data quality was screened for missing values and outliers
(Tukey inter-quartile-range fences), and distributional form assessed per trait with the
Shapiro–Wilk test, Q–Q plots and histograms. The **phosphorus effect** at Chokwe was estimated
per trait with a mixed model (treatment fixed; replication, genotype and genotype × treatment
random), using Kenward–Roger F-tests (`lmerTest`) and estimated marginal means (`emmeans`).

**Combining ability** was analysed with the combined line × tester model across environments.
General combining ability of the lines (GCA_L) and testers (GCA_T) and specific combining
ability (SCA) were tested against their respective environment-interaction mean squares (the
appropriate error strata for a multi-environment trial), and GCA and SCA effects were estimated
from the balanced cross means. **Variance components** were estimated by REML with
`sommer::mmer` (environment fixed; female, male, cross and their environment interactions, and
the spatial terms, random). From the variance components we derived additive variance
(AV = 2[σ²_GCA(L) + σ²_GCA(T)]), dominance variance (DV = 4 σ²_SCA), entry-mean broad-sense
(H²) and narrow-sense (h²) heritabilities, and **Baker's predictability ratio**
(BR = 2σ²_GCA / [2σ²_GCA + σ²_SCA]; BR → 1 indicates additive control). Superior hybrids were
ranked by BLUP genotypic value.

Genetic (among-hybrid) **correlations** among traits were computed on the 49 cross means, and
standardized **multiple regressions** of PUE and GY on the physiological and auxiliary traits
identified the traits most associated with each (excluding the exact identities TDM = SW + RW
and ASI = DTS − DTA, and GY from the PUE model). **Stability** of GY and PUE was assessed with
AMMI and GGE biplots and the WAAS/WAASY indices (`metan`).

# Results

## Trait variation and data quality

The 49 hybrids were fully replicated (750 plots, no missing values in the analysed traits). The
grand-mean grain yield was 6.5 t ha⁻¹, ranging among environments from 5.25 t ha⁻¹ (Nhacoongo)
to 7.41 t ha⁻¹ (Sussundenga). Ratio and interval traits (PUE, RSR, HI, ASI) departed most from
normality (Shapiro–Wilk), flagging them as transformation candidates.

## Effect of phosphorus supply (Chokwe)

At Chokwe, phosphorus stress significantly reduced plant height (162 → 145 cm, P < 0.001),
delayed anthesis (73.0 → 75.4 d, P < 0.001) and silking (74.9 → 78.1 d, P < 0.001), and widened
the anthesis–silking interval (1.9 → 2.7 d, P = 0.03) — the classical phenological signatures of
abiotic stress. Grain yield did not differ significantly between P levels (6.49 vs 6.83 t ha⁻¹,
P = 0.46), nor did the biomass traits; the apparent PUE difference is a mechanical consequence
of the applied-P denominator. Genotype × treatment interaction was significant for plant height
and days to anthesis, indicating that hybrids differed in their phenological response to P.

## Combining ability

In the combined line × tester ANOVA, tester GCA (GCA_T) was significant only for maturity (days
to anthesis and silking) and for root and total dry mass; it was non-significant for grain yield
and its components. Line GCA (GCA_L) was significant for harvest index, root weight, plant
height and maturity. Specific combining ability was significant for grain yield, root weight and
total dry mass, and the G×E interaction dominated PUE and grain yield. Thus the male testers
contributed essentially no additive variation for yield, while the female donor lines and
specific combinations carried the yield-relevant genetic signal.

## Variance components, heritability and gene action

Narrow-sense heritability was highest for the flowering traits (DTA h² = 0.60, DTS h² = 0.61)
with Baker's ratio ≈ 0.98, indicating almost purely additive, predictable control. Grain yield
had moderate broad-sense but low narrow-sense heritability (H² = 0.39, h² = 0.12, BR = 0.47),
and PUE was lowest (H² = 0.25, h² = 0.03, BR = 0.20), reflecting a predominantly non-additive
and G×E-driven architecture. Plant height and the biomass traits were intermediate (h² = 0.19–
0.27, BR = 0.52–0.72). Male-tester GCA variance for grain yield and its components was estimated
at the zero boundary (singular), confirming negligible additive differentiation among testers
for yield.

## Traits associated with PUE

Genetically, PUE was most strongly correlated with grain yield (r = 0.87, by definition) and
then with total dry mass (0.67), shoot dry weight (0.66), root weight (0.58) and plant height
(0.41). Standardized multiple regression identified **shoot dry weight and harvest index** as
the dominant drivers of PUE (R² = 0.81) and of grain yield (harvest index and shoot dry weight,
R² = 0.96). Because these driver traits are moderately heritable whereas PUE is weakly heritable
directly, they are candidate indirect-selection criteria.

## Stability

AMMI analysis of grain yield showed a significant genotype main effect but a non-significant
genotype × environment interaction, and the first two interaction principal components captured
≈ 72 % of the (limited) G×E, indicating broadly stable hybrid rankings. The GGE which-won-where
biplot resolved the environments into few mega-environments. Jointly ranking mean performance
and stability (WAASY), the most stable high-yielding hybrids were CML434 × EN20, CML434 × EN17
and CML530 × EN31.

# Discussion

## Genetic architecture of yield and PUE

Grain yield and, especially, PUE showed low narrow-sense heritability and low Baker's ratios,
with the yield-relevant variance concentrated in specific combining ability and genotype ×
environment interaction rather than in transmissible additive (GCA) effects. In contrast,
flowering time was highly heritable and almost purely additive. This pattern — heritable
phenology but weakly heritable, non-additive yield — has direct consequences for the breeding
strategy: direct selection on yield or PUE across environments will be inefficient, whereas the
more heritable, correlated traits offer a more reliable route.

## A narrow, yield-selected parental base explains the low additive variance for yield

The most striking result is the near-null tester GCA for grain yield: the seven male testers do
not differ in their average (additive) contribution to hybrid yield, so among-tester additive
variance for yield is effectively zero. We interpret this as a consequence of the **narrow,
double-ascertained genetic base of the parental lines**, rather than of any defect of the trial.

First, the testers were **twice filtered on yield**: they were drawn from an already-elite,
advanced breeding pool and then retained as the top seven for yield. Truncation selection on a
trait mechanically compresses the among-line variance for that trait in the selected set — the
survivors are, by construction, those most alike at the top of the yield distribution.

Second, the elite germplasm from which they derive is independently documented as narrow.
Masuka et al. (2017), fingerprinting the parents of the best CIMMYT eastern- and southern-African
hybrids, reported a low mean genetic distance (≈ 0.29) and found that four lines accounted for
29–58 % of the released hybrids, explicitly attributing the narrowing to "selection pressure for
defined traits." Fato et al. (2025) reported that seven years of IIAM breeding produced **no
significant genetic gain for grain yield** (a yield plateau), concluding that the parental lines
may lack superior alleles for yield beyond the current plateau — precisely the condition that
drives tester GCA for yield to zero. Crucially, the same studies show the base is **not**
genome-wide monomorphic: the pool did respond to selection for non-yield traits (ear height,
ears per plant, flowering), and Afonso (2013) showed that Mozambican **landraces** remain highly
diverse (expected heterozygosity ≈ 0.67). The narrowness is therefore a property of the elite,
yield-selected slice — not of the national germplasm.

Third, the residual among-tester differentiation is confined to exactly the traits selection did
**not** target. Tester GCA remained significant for maturity (days to anthesis and silking) and
for biomass, even though it vanished for yield. This is the fingerprint of **ascertainment**: the
testers were homogenized where they were selected (yield and its correlates) but still carry the
phenological and biomass differences that come free with diverse elite lines. Consistent with a
genuine boundary-zero rather than an estimation artefact, the yield tester-variance estimates sat
on the zero boundary (singular fits) under full replication.

Finally, the specific combining ability that does exist should be read cautiously: the SCA and
SCA × environment estimates were imprecise (Wald Z ≈ 1.3; the arithmetic SCA spread was of the
same order as the sampling error of a cross mean), so the non-additive signal for yield, while
present, is weakly determined. In short, the low additive variance for yield reflects a
narrow, yield-selected parental base; it does not imply that the germplasm is uniform for all
traits, and it does not preclude gains through complementary specific combinations or through
broadening the base with more diverse donors.

## Implications for selection

Because PUE and grain yield are governed largely by **harvest index and shoot biomass** — traits
of moderate, usable heritability — indirect selection on these secondary traits is expected to be
more effective than direct selection on PUE, which is weakly heritable and, across a P gradient,
partly definitional. The heritable phenological traits provide additional, well-characterized
selection handles and are relevant to escape/adaptation under stress. Broadening the tester base
with more genetically distinct donors would be expected to restore additive variance for yield
and increase the response to selection.

# Conclusions

In this maize line × tester population, grain yield and phosphorus-use efficiency were weakly
heritable and dominated by non-additive and genotype × environment variance, with the male
testers contributing essentially no additive variation for yield. This low additive variance is
best explained by the narrow, yield-selected base of the elite parental lines, corroborated by
independent diversity and genetic-trend studies, rather than by a lack of variation in the wider
germplasm. Harvest index and shoot biomass emerged as the heritable drivers of PUE and are
proposed as indirect-selection criteria, alongside exploitation of specific high-performing,
stable hybrids and, in the longer term, broadening of the parental base.

# References

Afonso AV (2013) Genetic diversity of local maize (*Zea mays* L.) germplasm from eight
agro-ecological zones in Mozambique. MSc Thesis, Swedish University of Agricultural Sciences,
Alnarp

Baker RJ (1978) Issues in diallel analysis. Crop Sci 18:533–536.
https://doi.org/10.2135/cropsci1978.0011183X001800040001x

Comstock RE, Robinson HF (1948) The components of genetic variance in populations of biparental
progenies and their use in estimating the average degree of dominance. Biometrics 4:254–266.
https://doi.org/10.2307/3001412

Fato P, Chaúque P, Senete C, Nhamucho E, Sneller C, Mutiga S, Musundire L, Wegary D, Das B,
Prasanna BM (2025) Genetic trends in seven years of maize breeding at Mozambique's Institute of
Agricultural Research. Agronomy 15(2):449. https://doi.org/10.3390/agronomy15020449

Isik F, Holland J, Maltecca C (2017) Genetic data analysis for plant and animal breeding.
Springer International Publishing. https://doi.org/10.1007/978-3-319-55177-7

Masuka BP, van Biljon A, Cairns JE, Das B, Labuschagne M, MacRobert J, Makumbi D, Magorokosho C,
Zaman-Allah M, Ogugo V, Olsen M, Prasanna BM, Tarekegne A, Semagn K (2017) Genetic diversity
among selected elite CIMMYT maize hybrids in East and Southern Africa. Crop Sci 57:2395–2404.
https://doi.org/10.2135/cropsci2016.09.0754

Olivoto T, Lúcio ADC (2020) metan: an R package for multi-environment trial analysis. Methods
Ecol Evol 11:783–789. https://doi.org/10.1111/2041-210X.13384

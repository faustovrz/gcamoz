# Zero Male GCA, Singular Fits, and Where the Genetic Variance Goes

**Trait:** Husked ear yield (HEY), grain yield `GY` (t/ha)
**Design:** Multilocation 7 × 7 line × tester (NC-II) — 7 CIMMYT/CML donor lines (females) × 7 IIAM `EN` testers (males; the legacy `MOZL` IDs), **5 environments** × 3 replications (four sites, with Chokwe grown at optimal *and* P-stress levels), alpha-lattice / RCBD layout
**Source notebooks:** [`multilocation_alpha_lattice.html`](multilocation_alpha_lattice.html); the trait-by-trait refit proposed in §3.5 is carried out in [`male_gca_trait_sweep.html`](male_gca_trait_sweep.html)

This note documents three connected questions that came up while interpreting the
combining-ability model: (1) what `isSingular()` means and why the male term triggers
it, (2) how the yield variance components compare with a published sweet × waxy maize
study, and (3) the **germplasm hypothesis** behind the zero male GCA — that the seven
yield-selected IIAM (`EN`) testers carry so little additive diversity they behave as
"essentially the same male" — including whether that extends to all phenotypes or only
to yield.

---

## 1. What `isSingular() == TRUE` means

From `?lme4::isSingular`:

> Evaluates whether a fitted mixed model is (almost / near) **singular**, i.e., the
> parameters are on the **boundary of the feasible parameter space**: variances of one
> or more linear combinations of effects are (close to) zero.

A variance **cannot be negative**. REML/ML searches for the variance that best fits the
data but is constrained to be ≥ 0. When the *unconstrained* optimum would be negative,
the optimizer pins the estimate to the boundary (0). That pinned-at-zero fit is exactly
what `isSingular()` flags.

This is the male (tester) term here:

- Fixed-effects ANOVA: male MS = 1.48 **< residual MS = 3.03** (p = 0.82).
- Because the between-male mean square is *below* residual, the unconstrained REML
  variance is **negative** → clamped to 0 → singular.
- The arithmetic SCA "spread" from family means (≈ 0.21) is now of the same order as the
  sampling error of a 15-plot cross mean (residual / 15 ≈ 0.16); the REML SCA variance
  (0.086) is correspondingly small and not significant (Wald Z ≈ 1.3).

### Common misconception

"Singular" is often read as *"not enough samples to compute the variance,"* like a
variance-0 case with a single observation. That conflates two different causes that both
produce a zero variance:

| Cause | What's happening | This study? |
|---|---|---|
| **Genuine zero signal** | Plenty of data; the variance is truly ≈ 0 (or below noise). The estimate is well-determined and it is 0. | **✅ Yes.** Full replication; the testers simply don't differ. |
| **Unidentifiable / too few levels** | Not enough levels or replication to *separate* the variance from others — it isn't estimable. | ❌ No. |

It is **not** "we couldn't compute it." It *was* computed, with adequate data, and the
answer is zero. REML's boundary clamping also hides a slightly-negative unconstrained
value, which is why the male term doesn't just shrink — it disappears.

### Practical consequence

While a fit is singular, Wald statistics, likelihood-ratio tests, and profile
confidence intervals are unreliable (per the help page). The clean fix is to **drop the
term pinned to zero** (the male random effect) and keep the female GCA and female GCA×E
terms, which carry the real signal.

---

## 2. Comparison with a published study (PMC10819934)

**Reference:** *The Combining Ability and Heterosis Analysis of Sweet–Waxy Corn Hybrids
for Yield-Related Traits and Carotenoids* — maize, **husked ear yield (HEY)**, NC Design
II (8 super-sweet females × 3 waxy males), 24 F₁ hybrids, 2 seasons.
<https://pmc.ncbi.nlm.nih.gov/articles/PMC10819934>

### Same units — with a trap

Both studies measure HEY in **t/ha (ton/ha)**, so every variance component is in
**(t/ha)²** and is directly comparable. **But** the large numbers in the paper's table
(GCA♂ 30.35, GCA♀ 21.23, SCA 29.38, …) are **mean squares, not variances** — they live
on a different scale (a mean square ≈ error + replication-coefficient × variance
component, inflated by how many plots were averaged). Compare *variance components to
variance components*: the paper's are σ²ₐ = 0.01 and σ²_d = 0.99.

### Apples-to-apples (variance components, (t/ha)²)

| Quantity | This study | PMC10819934 |
|---|---|---|
| GCA female (additive, ♀) | 0.0132 | — |
| GCA male (additive, ♂) | **0.0000** (singular) | — |
| **Additive σ²ₐ total** | **≈ 0.013** | **0.01** |
| **Non-additive / SCA (σ²_d)** | **0.0856** | **0.99** |
| Total G×E | 0.179 | (significant; not split into VC) |
| Residual variance | 2.347 | error MS 0.39 |
| Narrow-sense h² | **0.045** | **0.01** |
| Broad-sense H² | 0.339 | — |

> The paper's mean-square columns (GCA♂ 30.35, etc.) are *not* on this scale — do not
> compare them to the variance components above.

### Interpretation — different genetic architecture, not "too small"

- **The published yield is almost pure non-additive:** σ²_d = 0.99 dwarfs σ²ₐ = 0.01, so
  narrow-sense h² ≈ 0.01. Yield comes from specific ♀×♂ combinations (heterosis) — as
  expected for a sweet × waxy hybrid program built to exploit dominance.
- **This study's yield is mostly non-additive too,** but for a different reason: the
  among-cross genetic variance is dominated by **SCA (0.086)** over a small additive part
  (**0.013**, all of it female — male GCA is zero), so narrow-sense h² is only **0.045**.
  Across five environments the donor lines barely separate on their *across-environment
  average* — most of their action shows up as **female × environment** (see §3.2), and the
  testers add nothing additive.
- On the trait scale (√variance): female-GCA SD ≈ 0.11 t/ha and SCA SD ≈ 0.29 t/ha on a
  residual SD ≈ 1.53 t/ha. The genetic SDs are small against ~1.5 t/ha plot noise, which
  is why across-environment heritability is low even after averaging over
  5 environments × 3 reps.

---

## 3. Hypothesis: low genetic diversity in the male testers ("essentially the same male")

The zero male GCA in Section 1 is not just a statistical curiosity — it has a concrete
**germplasm explanation**. The working hypothesis is:

> The seven IIAM (`EN`) testers (males) carry very little additive genetic diversity for yield —
> they behave, for combining-ability purposes, as if they were *essentially the same
> male*. The breedable genetic variation comes from the introduced CIMMYT/CML donor lines
> (females) and from how specific crosses behave across environments.

This section sets out (a) *why the breeding history predicts this*, (b) *what the data
show*, (c) *whether the prediction extends to all phenotypes* — the key question — and
(d) *a falsifiable test*.

### 3.1 Rationale — why we expect near-zero male diversity

**Double ascertainment of the testers.** Stelio screened **70+ elite lines** and kept the
**top 7 for yield under phosphorus sufficiency**. Two filters stack here:

1. The lines were already *elite* — i.e., drawn from a **narrow, advanced breeding pool**,
   not from the broad landrace base.
2. They were then **truncation-selected on yield itself**. Selecting the extreme upper
   tail of a trait mechanically compresses the among-line variance *for that trait* in the
   selected set: the survivors are, by construction, the ones that are most alike at the
   top.

**The elite Mozambican pool is itself narrow and yield-plateaued.** Two local sources
document this:

- **Fato et al. 2025** (*Agronomy* 15(2):449 — genetic trends in 7 years of IIAM maize
  breeding) found **no significant genetic gain for grain yield** over 2014–2020
  (p = 0.96; yield even declined within the AVT and VCU streams), against heavy
  G×E (≈ 29% of variance) and error (≈ 66%). The authors state plainly: *"it is also
  possible that the parental lines used in the program do not have superior alleles for
  developing superior hybrids with grain yield beyond the current plateau."* That is a
  direct, independent statement that the elite line pool has **little exploitable additive
  variance for yield** — exactly the condition that drives male GCA to zero here.
- The IIAM elite base traces to **a handful of founders**: the first hybrids came from one
  South African inbred line plus two inbred lines from Ghana, and from IITA-sourced
  populations, with later material from the CIMMYT regional network; the program's
  germplasm is supplied by CIMMYT and IITA (Fato et al. 2025; Afonso 2013, §1.4). A small
  founder set means elite lines can be closely related — low *general* (additive)
  divergence among them.
- **The upstream CIMMYT pool is independently measured as low-diversity — Masuka et al.
  2017.** GBS-fingerprinting (258,038 SNPs) of the **55 parents** behind 52 of the 67 best
  CIMMYT Eastern & Southern Africa hybrids released 2000–2010 gave a **mean pairwise genetic
  distance of only 0.294** (range 0.004–0.4005; *every* pair < 0.45). Their verdict, opening
  and closing the paper: ***"Diversity was low."*** Crucially, **just four lines — CML444,
  CML395, CML312, CML442 — were each used in 15–30 of the 52 hybrids (29–58%)**, which
  *"gave 29 to 58% of the hybrids a **narrow genetic base**."* And the stated cause is the
  same mechanism as Stelio's selection: ***"selection pressure for defined traits can result
  in the narrowing down of the genetic base"*** (all lines "improved for yield and adaptation
  to the mid-altitude environment"). The IIAM (`EN`) testers are drawn from exactly this
  yield-and-adaptation-selected, four-line-dominated pool — so near-zero additive divergence
  for yield among seven of them is the expected outcome, not a surprise.
- **Founder base, and the "Ecuador" recollection.** CIMMYT's tropical heterotic groups rest
  on a few founder populations: **Group A** is largely **Tuxpeño**-derived (within it the
  CMLs are "closely related," max genetic distance ≈ 0.34), while **Group B** is built on
  **ETO, *Ecuador 573*, Lancaster, Mo17 and Southern Cross**; African drought-tolerant lines
  also trace to Tuxpeño Sequía / La Posta Sequía. So the half-remembered "couple of CML lines
  from Ecuador" is most likely **Ecuador 573, a heterotic-group-B founder population** — a
  founder *source*, not specific CMLs. The four dominant ESA lines span different heterotic
  groups (pairwise 0.266–0.282), which is why they serve as complementary testers even though
  the *overall* released-hybrid base is narrow.

**The narrowness is a property of the *elite* pool, not of Mozambican maize.** This is the
important contrast that keeps the hypothesis honest. **Afonso (2013)** (SLU MSc thesis;
`data/afonso_a_130426.pdf`) genotyped 27 Mozambican **landrace** accessions at 11 SSRs and
found them **highly diverse** — 84 alleles (7.6/locus), expected heterozygosity He ≈ 0.67,
with **88.3% of variation residing *within* accessions** and no clean structure by
agro-ecological zone. So the raw national germplasm is rich; the elite testers are a
deliberately **narrow, selected slice** of it. That is precisely why the experiment brings
in **CIMMYT/CML donor lines as the females** — to inject diversity the local testers lack.
The design is essentially asking *"do the introduced donors add additive value that the
local elite testers cannot?"* — and the answer (Section 2) is **yes**.

### 3.2 What the data show — variance routing

Treat the hypothesis as a set of predictions about where variance lands, then read off the
fitted components:

| Statement | Channel | This study (t/ha²) |
|---|---|---|
| Males carry little additive diversity | σ²_GCA(male) → 0 | **0.0000** ✓ |
| Introduced female lines differ | σ²_GCA(female) | 0.0132 |
| Contribution "depends on the female" (♀×♂) | σ²_SCA | 0.0856 |
| Male response "depends on environment" (♂×E) | σ²_GCA(male)×E | 0.0308 |
| ♀ × environment | σ²_GCA(female)×E | **0.1484** (largest) |
| ♀×♂ × environment | σ²_SCA×E | 0.0000 (singular) |

**The main-effect prediction holds.** Vg = GCA_F + GCA_M + SCA = 0.0132 + 0 + 0.0856 =
0.099 — the male contributes nothing, so the among-cross genetic variance sits in
**female + SCA**, though even the female main effect is modest. In the ×E layer the male
interaction is small (male×E ≈ 0.031), **female×E (0.148) is now the single largest
genetic component**, and SCA×E has itself collapsed to zero (singular). The male additive
channel is dead; the live genetic signal is the donor lines' *environment-specific*
performance.

**Caveat (routing of "environment").** Much of the donor-line signal does **not** fold into
the female *main* effect — it routes into **female×E**: the CML lines re-rank across the
five environments, most sharply between the optimal sites and the Chokwe P-stress site.
That female×E is the largest genetic term is literally "which donor line is best depends on
the environment," and it is why the across-environment narrow-sense h² (0.045) is so low —
the additive signal that *would* transmit is largely environment-specific.

### 3.3 Does "essentially the same male" extend to *all* phenotypes?

This is the crux of your reasoning: *if the top seven are essentially the same line, the
zero-male-variance result should hold not just for yield but for every trait.* **That
conclusion does not follow automatically — it depends on which of two mechanisms produced
the zero.**

| Mechanism | What it compresses | Holds for all traits? |
|---|---|---|
| **A. Ascertainment / truncation on yield** | among-tester variance for **yield and traits genetically correlated with yield** | **No** — traits uncorrelated with yield can still vary |
| **B. Narrow genetic base / near-identity-by-descent** | among-tester variance for **every trait** (the lines are genome-wide alike) | **Yes** |

Your "all phenotypes" claim requires **Mechanism B**. The evidence says the truth is a
*mix that leans toward A*, not pure B:

- **Against pure B:** Fato et al. (2025) show the *same* elite pool **did** respond to
  selection for several non-yield traits over 7 years — significant genetic gains for
  plant height (+0.67%/yr), ear height (+1.74%/yr), ears per plant (+1.31%/yr), ear
  position (+1.22%/yr), husk cover, anthesis date and ASI (−4.9%/yr). If the elite lines
  were genome-wide near-identical, *none* of those traits could have moved. So the pool
  retains additive variation for non-yield traits — it is **not** monomorphic.
- **For A:** the plateau and the "no superior alleles … for grain yield" statement are
  **yield-specific**. The strongest compression is exactly on the trait the testers were
  selected on.

So the defensible version of your hypothesis is **trait-specific, not universal**: the
testers are "essentially the same male" *for yield* (and for whatever is tightly
correlated with yield), but they are **not** guaranteed to be interchangeable for traits
that selection never targeted — flowering time, height, disease scores, etc. For those,
male GCA could well be non-zero.

### 3.4 A conceptual guardrail — zero male GCA ≠ genetically identical testers

Even for yield, read the zero carefully. GCA(male) = 0 means only that the testers don't
differ in their *average* (general, additive) effect. SCA is a different thing — whether a
*specific* ♀×♂ pairing beats the prediction from the two parents' averages. If the males
were *literally* interchangeable, the male identity wouldn't matter, **SCA would also
collapse to ≈ 0**, and everything would pile onto female GCA *alone*. The non-zero SCA
(0.086) says the males are **not** identical: they retain modest *specific-combination*
diversity (some ♀×♂ pairings beat their parents' average) even though their
*additive/average* diversity is nil. Its environment interaction, SCA×E, is itself ≈ 0 —
the specific combinations that exist are stable across environments.

> Sharpened premise: the testers have little **general (additive)** diversity for yield,
> but retain modest **specific-combination (SCA)** diversity — and, per §3.3, likely retain
> additive diversity for traits other than yield.

### 3.5 A falsifiable test

The hypothesis makes a clean, testable prediction: **refit the identical model
trait-by-trait** (flowering, plant/ear height, ears per plant, husk cover, disease scores,
…) and inspect σ²_GCA(male) and its `isSingular()` flag for each.

- If male GCA is singular/zero **only for yield** (and yield-correlated traits) but
  **non-zero** for, say, flowering or height → **Mechanism A** confirmed: the zero is an
  *ascertainment artifact of selecting on yield*, and the result does **not** generalize.
- If male GCA is singular/zero for **essentially every trait** → **Mechanism B** supported:
  the seven really are near-identical-by-descent, and your "holds for all phenotypes"
  extrapolation stands.

Fato et al. (2025) makes the first outcome the more likely one a priori (non-yield traits
respond to selection in this germplasm), but it is worth confirming directly in *these*
seven testers.

**This refit has now been run** ([`male_gca_trait_sweep.html`](male_gca_trait_sweep.html)):
across all 16 plot-level phenotypes the male term is singular/null for grain yield and every
yield component, for biomass, and for the efficiency/index ratios, and differs (boundary-
corrected LRT p < 0.05) **only for maturity** — days to silking and days to anthesis. That
is exactly the **Mechanism A** signature: the testers are homogenised where they were
selected (yield and its correlates) but still carry the phenology differences selection
never targeted. The "essentially the same male" claim holds *for yield*, not genome-wide.

### Honest statistical caveat — is SCA real?

The SCA point estimate is imprecise — Wald Z ≈ 1.3 (under 2), and the arithmetic SCA spread
(≈ 0.21) is of the same order as the sampling error of a 15-plot family mean (≈ 0.16); SCA×E
is estimated at the boundary (≈ 0). So whether σ²_SCA is *real signal* or *sampling noise*
is genuinely uncertain. (Wald Z on variance components near a boundary is itself unreliable;
firmer support comes from LRT/ANOVA.) By the fixed-effects ANOVA the **female** main effect
still differs (p = 0.008), and **female×E** is the strongest random genetic term (Z ≈ 1.8,
p ≈ 0.03) — though the female *main-effect* REML variance is small precisely because most of
the donor-line signal has routed into that female×E term. The male is null on every test.

**Bottom line:** the genetic variance is **female-driven, but largely through
female×environment** rather than a stable female main effect; the **male is null for
yield** on every test; and SCA / SCA×E are small (SCA×E ≈ 0). The germplasm history explains
the null male cleanly; the open question — whether the testers are "the same male" *for
yield specifically* or *genome-wide* — is settled by the trait-by-trait refit (§3.5), which
now supports the yield-specific (Mechanism A) reading.

> *On the "Ecuador" recollection (now resolved):* it most plausibly refers to **Ecuador
> 573**, one of the **founder populations of CIMMYT heterotic group B** (alongside ETO,
> Lancaster, Mo17 and Southern Cross); group A is largely Tuxpeño-derived. So the source is a
> *founder population*, not a couple of Ecuadorian CMLs — see §3.1 and Masuka et al. (2017).

---

## Take-home

1. `isSingular() == TRUE` here is **correct and informative**: the male variance is a
   genuine boundary zero (signal below noise), not an estimation failure.
2. Yield variances are in the **same units (t/ha²)** as PMC10819934 — but compare
   *variance components*, never the paper's *mean squares*.
3. The genetic architecture differs: that study is **dominance-led** (h² ≈ 0.01); this one
   has a **null male tester**, a small additive part carried only by the females, and most
   of its genetic action in **SCA and female×environment** — so across-environment
   narrow-sense h² is low (≈ 0.05).
4. The null male tester has a **germplasm cause**: the seven were double-ascertained
   (an already-narrow elite pool, then truncation-selected on yield). That pool is
   independently documented as low-diversity (CIMMYT ESA mean genetic distance 0.294, four
   lines in 29–58% of hybrids — Masuka et al. 2017) and yield-plateaued with "no superior
   alleles" for yield (Fato et al. 2025) — even though Mozambican *landraces* are diverse
   (Afonso 2013). The "Ecuador" founder is **Ecuador 573**, a CIMMYT heterotic-group-B
   founder population.
5. "Essentially the same male" is well-supported **for yield**, but **not** as a
   genome-wide claim: the same elite pool responded to selection for height, ASI, ears
   per plant, etc., so male GCA may be non-zero for non-yield traits. The decisive test is
   a **trait-by-trait refit** (§3.5).

---

## Sources

- **Fato, P., Chaúque, P., Senete, C., Nhamucho, E., Sneller, C., Mutiga, S., Musundire, L.,
  Wegary, D., Das, B. & Prasanna, B.M. (2025).** *Genetic Trends in Seven Years of Maize
  Breeding at Mozambique's Institute of Agricultural Research.* **Agronomy** 15(2), 449.
  https://doi.org/10.3390/agronomy15020449 — local copy: `data/FATO2025.txt`. *Used for:*
  the yield genetic-gain plateau (p = 0.96), the "no superior alleles for grain yield"
  statement, significant non-yield gains (height, ASI, EPP, husk cover), and the
  CIMMYT/IITA founder history.
- **Masuka, B.P., van Biljon, A., Cairns, J.E., Das, B., Labuschagne, M., MacRobert, J.,
  Makumbi, D., Magorokosho, C., Zaman-Allah, M., Ogugo, V., Olsen, M., Prasanna, B.M.,
  Tarekegne, A. & Semagn, K. (2017).** *Genetic Diversity among Selected Elite CIMMYT Maize
  Hybrids in East and Southern Africa.* **Crop Science** 57, 2395–2404.
  https://doi.org/10.2135/cropsci2016.09.0754 (open access). *Used for:* low diversity of the
  CIMMYT ESA elite parent pool (mean GD 0.294; "diversity was low"), the four lines
  (CML444/CML395/CML312/CML442) dominating 29–58% of hybrids, and the heterotic-group founder
  populations (Tuxpeño for A; ETO/Ecuador 573/Lancaster/Mo17/SC for B).
- **Afonso, A.V. (2013).** *Genetic diversity of local maize (Zea mays L.) germplasm from
  eight agro-ecological zones in Mozambique.* MSc thesis, Swedish University of Agricultural
  Sciences (SLU), Dept. of Plant Breeding, Alnarp — local copy: `data/afonso_a_130426.pdf`.
  *Used for:* SSR diversity of Mozambican landraces (84 alleles, He ≈ 0.67, 88.3% variation
  within accessions) and the IIAM breeding-strategy / founder background (§1.4).
- *The Combining Ability and Heterosis Analysis of Sweet–Waxy Corn Hybrids for
  Yield-Related Traits and Carotenoids* (maize HEY, NC Design II), PMC10819934 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10819934 (Section 2 comparison; authors/year not
  captured here).
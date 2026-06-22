# Zero Male GCA, Singular Fits, and Where the Genetic Variance Goes

**Trait:** Husked ear yield (HEY), `Yield t/ha`
**Design:** Multilocation 7 × 7 line × tester (NC-II) — 7 CIMMYT/CML donor lines (females) × 7 MOZL testers (males), 4 environments × 3 replications, alpha-lattice / RCBD layout
**Source notebook:** [`multilocation_alpha_lattice.html`](multilocation_alpha_lattice.html)

This note documents three connected questions that came up while interpreting the
combining-ability model: (1) what `isSingular()` means and why the male term triggers
it, (2) how the yield variance components compare with a published sweet × waxy maize
study, and (3) whether "all genetic variance loads onto the female and SCA channels" is
the right prediction when the male testers carry little diversity.

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

- Fixed-effects ANOVA: male MS = 1.55 **< residual MS = 2.23** (p = 0.65).
- Because the between-male mean square is *below* residual, the unconstrained REML
  variance is **negative** → clamped to 0 → singular.
- The apparent SCA "spread" from family means (≈ 0.18) is fully accounted for by the
  sampling error of a 12-plot mean (residual / 12 ≈ 0.19).

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
| GCA female (additive, ♀) | 0.0715 | — |
| GCA male (additive, ♂) | **0.0000** (singular) | — |
| **Additive σ²ₐ total** | **≈ 0.072** | **0.01** |
| **Non-additive / SCA (σ²_d)** | **0.0716** | **0.99** |
| Total G×E | 0.195 | (significant; not split into VC) |
| Residual variance | 1.530 | error MS 0.39 |
| Narrow-sense h² | **0.224** | **0.01** |
| Broad-sense H² | 0.448 | — |

> The paper's mean-square columns (GCA♂ 30.35, etc.) are *not* on this scale — do not
> compare them to the variance components above.

### Interpretation — different genetic architecture, not "too small"

- **The published yield is almost pure non-additive:** σ²_d = 0.99 dwarfs σ²ₐ = 0.01, so
  narrow-sense h² ≈ 0.01. Yield comes from specific ♀×♂ combinations (heterosis) — as
  expected for a sweet × waxy hybrid program built to exploit dominance.
- **This study's yield is roughly half additive, half non-additive,** and the additive
  half is entirely the **female (CML donor) lines** — male/tester GCA is zero. That gives
  a real narrow-sense h² of 0.22: the donor lines carry breedable additive value; the
  testers don't differentiate.
- On the trait scale (√variance): female-GCA SD ≈ 0.27 t/ha and SCA SD ≈ 0.27 t/ha on a
  residual SD ≈ 1.24 t/ha. "Small" is relative — modest genetic SDs against ~1.2 t/ha
  noise, which is why heritability is only moderate, and only after averaging over
  4 environments × 3 reps.

---

## 3. "All genetic variance loads onto female + SCA" — assessment

**Premise:** the male testers have little genetic diversity, and their contribution
depends on the female (cross) and the environment.
**Prediction:** all genetic variance goes into the female and SCA channels.

**Verdict: mostly right, with one refinement and one conceptual correction.**

### Routing rules and the data

| Statement | Channel | This study (t/ha²) |
|---|---|---|
| Males have little genetic diversity | σ²_GCA(male) → 0 | **0.0000** ✓ |
| Female lines differ | σ²_GCA(female) | 0.0715 |
| Contribution "depends on the female" (♀×♂) | σ²_SCA | 0.0716 |
| "depends on environment" (♂×E) | σ²_GCA(male)×E | 0.0043 |
| ♀ × environment | σ²_GCA(female)×E | 0.0728 |
| ♀×♂ × environment | σ²_SCA×E | **0.1176** (largest) |

**The main-effect prediction is exactly right.** Vg = GCA_F + GCA_M + SCA =
0.0715 + 0 + 0.0716 = 0.143 — the male contributes nothing, so all genetic variance sits
in **female + SCA**. The same pattern repeats in the ×E layer: male×E ≈ 0.004
(negligible), while female×E + SCA×E carry the G×E.

### Refinement — "depends on environment" is its own channel

Environment-dependence does **not** fold into the female or SCA *main* effects; it routes
into the **×E** terms. Here those happen to be female-side (0.073) and SCA-side (0.118),
so the prediction holds — but only if the ×E versions are explicitly included. SCA×E
being the *largest* genetic component is precisely "the cross contribution depends on the
environment."

### Conceptual correction — zero male GCA ≠ genetically uniform males

GCA(male) = 0 means only that the testers don't differ in their *average* (general,
additive) effect. SCA measures something different — whether a *specific* ♀×♂ pairing
beats the prediction from the two parents' averages.

The trap: if the males were *truly* genetically uniform (interchangeable), the male
identity wouldn't matter, so **SCA would also collapse to ≈ 0** and *everything* would
pile onto female GCA *alone* — not female + SCA. The nonzero SCA (0.0716) and SCA×E
(0.118) are positive evidence the males are **not** uniform: they carry *specific
/interaction* diversity (they combine and respond to environments differently) even
though they lack *additive/average* diversity. Sharpen the premise to:

> Males have little **general (additive)** diversity, but retain **specific
> (interaction)** diversity.

These are independent; the data show the first is zero and the second is not.

### Honest caveat — is SCA real?

The SCA and SCA×E point estimates are imprecise — Wald Z-ratios ≈ 1.1 (well under 2), and
the arithmetic SCA spread (≈ 0.18) is about equal to the sampling error of a 12-plot
family mean (≈ 0.19). So whether σ²_SCA is *real signal* or *sampling noise* is genuinely
uncertain. (Wald Z on variance components near a boundary is itself unreliable; the firm
support comes from LRT/ANOVA.) The only genetic terms with solid LRT support are **female
GCA (p = 0.0002)** and **female GCA×E (p = 0.004)**.

**Bullet-proof conclusion:** the genetic variance is **female-driven (GCA and GCA×E)**,
the **male is null**, and SCA/SCA×E are where any male-by-cross-by-environment signal
*would* land — but in these data that signal is too small to distinguish from noise. If
the SCA is in fact noise, the situation reduces to the clean end-member of the original
prediction: *males uniform → all variance in female GCA.*

---

## Take-home

1. `isSingular() == TRUE` here is **correct and informative**: the male variance is a
   genuine boundary zero (signal below noise), not an estimation failure.
2. Yield variances are in the **same units (t/ha²)** as PMC10819934 — but compare
   *variance components*, never the paper's *mean squares*.
3. The genetic architecture differs: that study is **dominance-led** (h² ≈ 0.01); this
   one is **female-additive-led** (h² ≈ 0.22) with a null male tester.
4. "Variance loads onto female + SCA" is right for the main effects; remember the **×E**
   channels for environment-dependence, and don't read zero male GCA as genetically
   identical testers.
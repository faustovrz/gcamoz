# Delete candidates — scripts/qtl/

Proposed only. Nothing here has been deleted. FVRZ decides.

Already removed (written today, never run, superseded by `build_map.R`):
`qtl_build_map_airmine.R`, `qtl_map_teonam_filters.R`.

## Dead ends — the CML530 pileup route

Superseded by `cml530_alleles_direct.R`, which reads the assembly base directly.
The pileup route failed twice over: every site was ~DP=1 where `bcftools call -m`
defaults to REF (2:1 ref/alt skew, 53.4% phase consistency), and "matches B73" is
not the ref/alt question at all — DArT's ref allele equals the B73 base at only
~76% of markers.

| file | why it can go |
|---|---|
`cml530_genotypes_pileup.sh` | the DP=1 pileup pipeline |
`cml530_pileup_consensus.R` | raw mpileup retry; same wrong question |
`cml530_attrition.R` | attrition accounting for the dead route |
`test_pileup_coords.R` | diagnosed the dead route (proved coords were fine, 97.6% delta==0) |

Keep the *finding* in the notebook; the scripts are not needed to reproduce it.

## Dead ends — my filters and my map attempts

| file | why it can go |
|---|---|
`make_abh_framework_only.R` | encoded only the 1,384 `qc_pass` markers. Superseded by `encode_abh.R` (3,728) |
`build_map_with_my_filters.R` | produced the 880 cM / 878-marker map. Includes my `LINK_MIN` filter |
`phase_markers.R` | phasing with my `LINK_MIN` marker drop bolted on; phasing now lives in `encode_abh.R`'s successor path |
`map_errorprob_sweep.R` | written, never run |
`detect_inversions.R` | the reversal test was mathematically vacuous — reversing a contiguous block leaves internal adjacency unchanged, so every gain was identically 0 |
`diagnose_chr4_chr9.R` | chased chr4/chr9 as bad markers; the Marey plot showed order was fine |
`diagnose_phase.R` | one-off; its conclusion (assembly read fixed phase, 53.4% → 93.3%) is in the notebook |
`marker_accounting.R` | one-off funnel accounting |
`audits.R` | the four inline blocks I should have written as files in the first place |

## Keep — the working pipeline

```
phenotype_dictionary_unify.R     data/phenotype_dictionary.csv
data_audit.R                     inventory, marker + sample QC
structure_check.R                single-F2 verification, no parent genotyped
liftover_v4_to_v5.R              AGPv4 -> NAM v5, tag FASTA
install_minibwa_from_source.sh   minibwa WITH the GPL low-memory BWT (-l)
index_genomes.sh                 CML530 + B73 v5 minibwa indexes (~45 min)
cml530_alleles_direct.R          A = CML530 from the assembly (3,728 markers)
encode_abh.R                     A/H/B, all 3,728, no filters
build_map.R                      FVRZ QC: distortion -> est.map -> find_quirky -> est.map
```

## Keep — supporting diagnostics

```
rqtl_format_probe.R              proved csvsr takes TWO leading columns
minibwa_index_bench.sh           index size/time
minibwa_index_threads_bench.sh   -l costs 8.9x; thread scaling untestable in this build
index_wait_verify.sh             artefact-based index validation
plot_map.R                       ideogram + Marey plots
pick_trait_by_h2.R               diallel h2 ranking (TDM 0.266 highest F2-mappable)
```

## Note

Several dead-end scripts are referenced by `qtl_marker_pipeline.qmd`, which is
committed (`87f0b87`). Deleting them requires editing the notebook in the same
pass or its citations break:
`cml530_genotypes_pileup.sh`, `cml530_attrition.R`, `test_pileup_coords.R`,
`make_abh_framework_only.R`.

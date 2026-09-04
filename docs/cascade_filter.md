# Cascade Filtering

The cascade filter is RNAConSnake's multi-stage screening that identifies statistically significant RNA structure predictions. It combines three independent quality metrics applied in sequence to real and null alignments.

## Metrics

### 1. RNAz Probability
[More text about RNAz...]

The RNAz class probability filter selects alignments with high structural conservation. Candidates failing this stage never reach downstream filters.

- **Threshold**: `calibration.rnaz_prob_threshold` (default: 0.9)
- **Stage-one threshold**: `calibration.stage1_rnaz_prob` (default: 0.5) — AlifoldZ is computed only on stage-one survivors when using two-stage calibration

### 2. AlifoldZ Z-Score
The AlifoldZ filter screens for thermodynamic stability of consensus structure. A negative z-score indicates the consensus structure is less stable than random structure for the sequence length and composition.

- **Threshold**: `calibration.alifoldz_threshold` (default: -2.0)

### 3. R-scape Covariation Statistics
R-scape detects covariation in aligned sequences — an evolutionary signal that indicates structural constraint. Starting with the full set of base pair covariations it observes, the RNAConSnake cascade now evaluates three independent R-scape statistics:

#### Minimum Covarying Pairs
The count of base pairs observed to covary in the alignment. This count must be at least the configured minimum to pass the filter.

- **Threshold**: `calibration.rscape_min_pairs` (default: 1)
- **Context**: A count of 0 indicates no base pair covariations were detected (a real, reported result, not a missing value)

#### Average Confidence
The average conditional probability for observed covariations. Confidence measures the probability of covariation given the current base pair assignment and is stable across different sequence compositions.

- **Threshold**: `calibration.rscape_min_confidence` (default: 0.5, range [0, 1])
- **Interpretation**: Values near 1.0 indicate high confidence; values near 0.5 indicate borderline confidence

#### Mutual Information
The mutual information content (in bits) of observed covariations. Mutual information quantifies the dependency between paired bases and is sensitive to the number of sequences in the alignment.

- **Threshold**: `calibration.rscape_min_mutual_info` (default: 0.1, range [0, ∞])
- **Context**: A value near 0 indicates weak covariation signal; values >0.5 typically indicate strong covariation

## Cascade Behavior

The cascade applies thresholds in sequence on the de-replicated locus set:

1. **RNAz filter**: Loci with `rnazprob >= rnaz_prob_threshold` survive
2. **AlifoldZ filter**: Survivors from step 1 with `alifoldz <= alifoldz_threshold` survive
3. **R-scape filter**: Survivors from step 2 with ALL of the following survive:
   - `rscape_covary_count >= rscape_min_pairs`
   - `rscape_avg_confidence >= rscape_min_confidence`
   - `rscape_mutual_info >= rscape_min_mutual_info`

## Disabling R-scape in the Cascade

When R-scape was not computed (e.g., `do_rscape: false`), the cascade automatically drops the R-scape stage, reporting counts only for RNAz and AlifoldZ filters. This ensures that run configuration does not artificially inflate or deflate headline counts.

## Null-Model Calibration

The cascade filter is evaluated on both the real alignment and parallel null-model runs (generated on shuffled or randomized alignments). The empirical false discovery rate (FDR) for cascade survivors is the mean number of null survivors divided by the number of real survivors:

```
FDR(cascade) = mean_null_survivors / real_survivors
```

Per-stage FDR is also computed for intermediate filters. The q-value for each locus is derived from the monotone envelope of FDR across all threshold combinations.

## Configuration Example

```yaml
calibration:
  rnaz_prob_threshold: 0.9
  alifoldz_threshold: -2.0
  stage1_rnaz_prob: 0.5
  rscape_min_pairs: 1
  rscape_min_confidence: 0.5
  rscape_min_mutual_info: 0.1
```

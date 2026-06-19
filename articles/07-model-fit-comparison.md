# Model fit and model comparisons

Running a vertex-wise analysis with `verywise` means fitting a model
hundreds of thousands of times. It better be a good model. So in this
article we go over some considerations around estimation method,
convergence, and model fit that are particularly important to ensure the
best model.

We talk about:

- what convergence warnings mean and when to act on them,
- the choice between REML and ML estimation,
- how to interpret per-vertex model fit statistics,
- how to compare model specifications vertex-wise

### Singular fits

A singular model is one where one or more random-effect variance
components are estimated at (or near) zero.

This may happen at vertices with very low signal, and are not
necessarily cause for alarm. However, if singular fits affect many of
your vertices (e.g. \> 30%), the random structure you specified is
likely too complex for the data and should be simplified (e.g. remove a
random slope or switch from correlated to uncorrelated random effects).

Singular fit are also sometimes caused by too few random-effect levels
(e.g. \< 5 groups; see more info
[here](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#singular-fits)),

Singular fit counts per vertex are stored in the
`<hemi>.<measure>.singular_fits.mgh` output file for post-hoc
inspection.

### Other converge warnings

Any other warnings emitted by `lme4` during a vertex-wise run are logged
in the `<hemi>.<measure>.issues.log`. The most common ones are failed
convergence warnings: when the optimizer did not reach a satisfactory
solution. They come in a variety of flavours, but you can resolve most
of them using `lmm_control`.

> Note: do not dismiss failed convergence warnings, parameter estimates
> are unreliable when the models do not converge.

In the [performance tips
article](https://seredef.github.io/verywise/articles/04-run-slurm-array.md)
we noted that setting `lmm_control = lmerControl(calc.derivs = FALSE)`
skips the finite-difference gradient and Hessian calculations that are
performed post-optimization to verify convergence. This speeds up
fitting but also silences some convergence diagnostics, so it is better
suited to final runs with a model you are already confident in.

### Model comparison

Vertex-wise AIC values (average in case of multiple imputations) are
stored in the `<hemi>.<measure>.aic.mgh` output file.

These values are useful for *comparing model specifications*, for
example, when evaluating whether to include a fixed or random term in
the model. They are generally not interpretable as standalone global fit
metrics.

More on AIC in the context of LMMs
[here](https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#can-i-use-aic-for-mixed-models-how-do-i-count-the-number-of-degrees-of-freedom-for-a-random-effect).

#### REML vs. ML

`verywise` takes a `REML` argument (***Restricted maximum likelihood***)
that controls the estimation criterion used. This is set to `TRUE` by
default, which is recommended in order produce unbiased estimates of
variance components and estimate random effects correctly.

However, when comparing models with \***different fixed-effects terms**
via AIC, you should fit them using ML (`REML = FALSE`). Models differing
only in their random-effects structure can be compared using REML fits.

#### Example model comparison pipeline

\[TODO\]

### ICC, marginal and conditional R²

`verywise` also outputs three complementary model fit statistics that
together help you get a more complete picture of the model.

- The intra-class correlation coefficient (ICC; stored in
  `<hemi>.<measure>.icc.mgh`) indicates the degree of clustering in each
  vertex, i.e. how much variance in the brain metric is attributable to
  the grouping factor (e.g. site, subject). A high ICC confirms that the
  random effect captures meaningful clustering and that the mixed-model
  structure is warranted.
- The marginal R² (stored in `<hemi>.<measure>.r2_marginal.mgh`)
  captures the variance explained by the fixed effects alone (relative
  to the total variance)
- The conditional R² (stored in `<hemi>.<measure>.r2_conditional.mgh`)
  captures the variance explained by both the fixed and random effects
  combined (i.e. the full model)

For example, in a longitudinal study, a large ICC with a small marginal
R² suggests that most of the variance comes from inter-individual brain
variability while fixed-effects predictors do not explain much.

Using the ([plotting functions in
`verywise`](https://seredef.github.io/verywise/articles/05-visualize-results.md))
you can visualize spatial patterns in marginal R², revealing, for
example, which cortical regions are most strongly associated with the
fixed-effect predictors in the model.

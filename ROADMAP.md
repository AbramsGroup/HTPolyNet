# Roadmap

Ideas worth doing that we haven't done yet. This is a living list, not a
commitment or a schedule. When something here gets done, it moves to
`CHANGELOG.md` and comes off this page.

Rough ordering within each section is by value, not by effort.

## Container and deployment

- **The container's Gromacs is generic `AVX2_256`, on both tags.** The
  `:cuda` image fixes the GPU half of this problem and not the SIMD half:
  conda-forge builds for a portable baseline, so on an AVX-512 host both
  images leave single-core throughput on the table relative to a natively
  built or module-provided Gromacs. This is inherent to installing Gromacs
  from conda-forge and cannot be fixed by choosing a different build string
  -- it would take building Gromacs in the image, which trades the weekly
  rebuild's freshness for a long build and a host-specific artifact.

  Confirmed on hardware 2026-08-25: the `:cuda` image on Picotte's gpu001
  reports `SIMD instructions: AVX2_256` on a node whose CPUs support
  AVX-512, and Gromacs itself prints the hint that "AVX_512 ...
  instructions will perform best on this hardware". Still unquantified, and
  that is what would decide whether it matters: nobody has run the same
  system against Picotte's own `abramsGrp-gromacs/2021.2/cpu-gpu` module on
  the same node. Note Gromacs also observes that AVX2 is often the better
  choice for runs that offload to a GPU anyway, so the penalty may be small
  for exactly the case the `:cuda` image serves.

- **Releases before v2.6.1 have no version tag in the registry.** From
  v2.6.1 on, `docker.yml` pushes `:v<version>` and `:<version>` (and the
  `cuda-` variants) alongside the commit sha, but the earlier releases were
  pushed as `:latest` and a sha only. So `pull ...:v2.5.0` returns `manifest
  unknown`, which reads like a failed build rather than a differently-named
  tag. The docs now say how to resolve an old release to its sha, which is
  the cheap half. Retagging the existing images would be better and is a
  registry operation, not a code change: `docker buildx imagetools create -t
  ghcr.io/cameronabrams/htpolynet:v2.5.0 ghcr.io/cameronabrams/htpolynet:<sha>`
  for each past release, which needs a token with `write:packages`. Worth
  doing if anyone reports hitting it a second time.

- **The image's Gromacs and AmberTools are unpinned, so two images with
  identical htpolynet code can compute different numbers.** `docker/Dockerfile`
  installs `ambertools`, `gromacs`, `parmed` and `rdkit` from conda-forge with
  no version constraints, and the weekly rebuild exists precisely to pick up
  whatever is newest. So the commit stamp now tells you which htpolynet
  produced a build, and still does not tell you which force-field tool chain
  did. Observed live: the image built 2026-08-25 carries AmberTools 26.0 and
  Gromacs 2026.3, while panacea's native environment is on Gromacs
  2025.4 -- a different major version against the same htpolynet.

  Unpinning is deliberate and mostly right: pinning would freeze the image on
  old Gromacs and defeat the point of a weekly rebuild. The gap is that
  nothing *records* what a given image resolved to, so the versions are
  discoverable only by running `htpolynet info` inside it and writing the
  answer down by hand. The better of two candidates is an image
  **label** carrying a `conda list` export, written at build time: `docker
  inspect` then answers the question without running the image, so the record
  survives reaching someone who cannot run the container at all -- a reviewer,
  an archive, a future reader holding only the digest. The alternative, making
  `htpolynet info` machine-readable so a build can capture it, requires the
  image to still be runnable, which is the weaker guarantee. This is the same
  requirement as the build manifest below, one layer further down -- what
  produced this build, all the way to the compilers.

  Raised by the calibration study, which found the 2025.4/2026.3 split only
  because it went looking after an unrelated prompt.

- **Retire `ghcr.io/abramsgroup/htpolynet`.** Superseded by the
  `cameronabrams` package; still public and still serving a June image to
  anyone with an old link.

## Testing and CI

Coverage as of the last measurement: **38.8%** overall.

- **`repair/`'s driver has no test.** `test_cap_placement.py` and
  `test_repair_conversion.py` now cover the placement search and the
  reported statistics, but `triazine_to_cyanate_cap` itself and all of
  `topology_surgery.py` (125 statements) are still untouched. This is the
  highest-value gap: the postcure repair stage makes the strongest
  correctness claim in the project ("atom conservation is exact"), and
  right now the only thing checking it in-tree is reading a residue census
  at the end of a multi-hour build. (An external audit of 54 builds at
  v2.6.2 found the accounting exact everywhere -- `3*TAZ_final + CYN ==
  720` in all 54, every surviving triazine with exactly 3 aryl-ether bonds
  across 7,501 examined, zero bare -OH -- and a six-defect negative control
  was caught by 2-10 checks each. That is real evidence, but it is not a
  test and it does not run on a PR.) It is pure topology manipulation, so
  it can be tested deterministically in milliseconds against a synthetic
  `TopoCoord` carrying triazines at k=0,1,2,3 — assert atom counts, the
  residue census, cap placement, and that no unreacted bridge -OH
  survives.
- **An end-to-end example in CI.** A deliberately tiny build (a
  20-molecule, few-ps variant of example 0) run inside the container
  would cover `core/runtime.py` and `cure/curecontroller.py` — 1,056
  statements, both at 0% — in the only way that is honest, since faking
  the whole AmberTools/Gromacs tool chain to unit-test the orchestration
  is a large effort for less confidence.
- **Remaining zero-coverage modules**: `analysis/postsim.py` (203),
  `cli.py` (165), `analysis/analyze.py` (137), `utils/vmd_viz.py` (75),
  `utils/checkpoint.py` (54).
- **`analysis/plot.py` is at 34%** after the smoke-test pass. The
  diagnostics-log parsers (`diagnostics_graphs`, `_token_match`,
  `_parse_data`) are the part most likely to rot silently — they already
  broke once when modules were renamed — and they are testable against a
  small captured log fixture.
- **Coverage reporting in CI**, so the number is visible on a PR rather
  than something we remember to measure by hand.

## Release and distribution

- **The release preflight cannot tell whether the *previous* release
  actually shipped to conda-forge.** `scripts/check-conda-sync.py` compares
  `pyproject.toml`'s runtime deps against the feedstock recipe, which
  catches dependency drift the autotick bot cannot handle. It says nothing
  about whether the bot's last PR ever merged. That gap ran for four
  releases: `pip check` in the recipe's test section started failing when
  `ambertools` began pulling in distributions with unsatisfiable metadata,
  so the bot PRs for 2.2.0, 2.3.0, 2.3.1 and 2.4.0 all sat red and unmerged
  while conda-forge served 2.1.0 from June. Every one of those releases
  passed the preflight, because the deps genuinely did match. Nobody
  noticed until a bot email got read.

  Fixed on 2026-08-26 by dropping `pip check` from the recipe's test section
  (it fails on `ambertools`' bundled distributions, never on ours) and
  merging the 2.4.0 bump, which closed the other three. **The damage is
  permanent and visible**: conda-forge's version list for this package now
  reads 1.0.9 -> 2.1.0 -> 2.4.0, because 2.2.0, 2.3.0 and 2.3.1 were never
  built there and never will be. That is worth knowing as a diagnostic in
  its own right -- a published version list that skips releases the project
  actually made is the retroactive signature of this failure, in any
  package, without needing to have run any check at the time.

  The check is cheap: query `api.anaconda.org/package/conda-forge/htpolynet`
  for `latest_version` and compare it against the version being superseded.
  A mismatch does not have to block the release -- the fix is usually on the
  feedstock, not here -- but it must be loud, because the failure mode is
  silence. Consider also listing open PRs on the feedstock, since a red bot
  PR is the specific thing to look at.

  The deeper point is that publishing to conda-forge is the one leg of the
  release that completes *after* `release.sh` exits and on someone else's
  infrastructure, so it is the only one that can fail without anything here
  noticing.

- **External services still keyed to the old repo identity.** The Aug 2026
  transfer from `AbramsGroup/HTPolyNet` to `cameronabrams/htpolynet` moved
  the code but left every integration pointing at the old owner. Three
  broke and were fixed during the 2.2.0 release: the GHCR package path
  (docs referenced a package that had never been published under the new
  owner), PyPI trusted publishing (`invalid-publisher` — the claim no
  longer matched, so the v2.2.0 upload failed until the publisher was
  re-registered), and Read the Docs, whose project `repository.url` is
  **still** `https://github.com/AbramsGroup/htpolynet`. Builds succeed
  anyway because GitHub redirects the clone, but tag versions never sync —
  `/en/v2.2.0/` 404s and the API reports "No Version matches the given
  query". Fix the RTD project URL, then activate the tagged version. Worth
  keeping this list as the checklist if the repo ever moves again.

  A fourth instance turned up on 2026-08-26 and is fixed: the conda-forge
  recipe's `about:` block still gave `AbramsGroup/HTPolyNet` for `home` and
  `dev_url`, and a `doc_url` of `abramsgroup.github.io/HTPolyNet` that
  returns 404. Corrected in the same feedstock PR that unblocked the version
  bumps.


- **Mint a software DOI.** Enable the Zenodo GitHub integration for
  `cameronabrams/htpolynet`, then the next `scripts/release.sh` run
  archives the release automatically. Afterwards, add the concept DOI (not
  the version DOI) as a README badge and an `identifiers` entry in
  `CITATION.cff`. Note that enabling is not retroactive: releases before
  the toggle are not archived. Also note that adding a `.zenodo.json`
  would make Zenodo ignore `CITATION.cff` entirely — only worth doing if
  we need Zenodo-specific fields such as `grants` for funder linkage.
- **Per-minor Python classifiers.** `pyproject.toml` declares only
  `Programming Language :: Python :: 3`, so the PyPI Python badge reads an
  uninformative `python: 3`. Adding `:: 3.10` through `:: 3.13` would make
  it read `3.10 | 3.11 | 3.12 | 3.13`, which CI already verifies at both
  ends.

## Cure and repair

- **`cap_min_clearance` has never been calibrated against a working metric.**
  The 0.150 nm default was chosen in 2.6.0 against a clearance that was pinned
  at 0.136 nm by a bug, so it fired on 100 % of caps for a reason that had
  nothing to do with crowding, and the 0.22 nm value rejected before it was
  rejected on the same broken measurement. The metric is real as of the fix,
  but the only two data points on a real box -- the accept260 builds on
  picotte, `/ifs/groups/abramsGrp/cfa22/htpolynet/cyanate-bridge-series/accept260/{unbiased,biased}/`
  -- were taken with the bug present and are uninformative. The geometry now
  bounds the achievable clearance at 0.272 nm (cap carbon to the attachment
  oxygen's aryl carbon, antiparallel), with the chemically sensible ~120
  degree placement at 0.236 nm -- but the aryl carbon has since been excluded
  from the metric precisely so it does not saturate, so that bound no longer
  applies and there is no analytic ceiling to reason from at all. What is
  unknown was the fire rate in a box at polymer density, and that has now been
  measured (see the real-box block at the end of this entry). What is left is
  the recalibration itself: set the default so it flags the tail rather than
  the bulk, using `blind_min_clearance_nm`, `blind_median_clearance_nm` and
  `n_below_target`. Do **not** calibrate against `min_clearance_nm` -- it
  saturates at whatever target is set. A threshold that fires on everything is
  not a threshold. The `n_preferred_out_of_angle` escape hatch is closed: it
  came back at 6 %, so the 90--150 degree window is not biting real O-H
  vectors and is not the thing to change.

  One thing the synthetic sweep did settle, so nobody re-derives it: the
  clipping the aryl-carbon exclusion removes is density-dependent, and at melt
  density it was mild. Median clearance against random neighbourhoods at 5,
  15, 33 and 60 heavy atoms per nm^3 ran 0.316 / 0.230 / 0.193 / 0.177 nm with
  the aryl carbon excluded and 0.233 / 0.222 / 0.191 / 0.174 nm with it in.
  It pins hard at 0.233 nm -- the C-O-C geometry at 118 degrees -- only in
  sparse surroundings, where the caps are comfortable anyway; by 33/nm^3 real
  neighbours are usually closer than the geometric bound and the two agree to
  about 2 %. So the median was censored from above rather than constant, and
  the exclusion buys an uncensored statistic rather than rescuing a dead one.
  Worth knowing before treating a change in that median as physics.

  For the cyanate-ester system specifically, the study session measured the
  box: heavy-atom count is conversion-independent at 7560 (14,040 atoms less
  6,480 H, and the cure removes only hydrogen), post-repair mass 100,192
  g/mol, and measured densities of 1.171--1.198 g/cm^3 across the as-cured,
  ladder-cold-end and monomer-melt states bracket **~53 heavy atoms per
  nm^3**. Repair runs before the postcure anneal so the placement-time value
  is not directly measured, but the spread between those three is smaller than
  the uncertainty over which stage to attribute. That is inside the swept
  range, so nothing here extrapolates. Rerunning the sweep at 53/nm^3 gives
  median blind clearance 0.143 nm (p10 0.075, p90 0.212), median searched
  clearance 0.179 nm, and about 45 % of caps keeping the O-H direction, with
  the censored and uncensored medians differing by under 1 %. Three things
  follow.

  First, and this needs stating at commit granularity rather than as "the
  clearance fix", because the two commits move the numbers by wildly different
  amounts in the same release: 5a58d3d, dropping the attachment oxygen, turns
  `median_clearance_nm` from the constant 0.136 into a distribution, so every
  placement field moves enormously between the accept260 runs and whatever
  comes next. That is the released bug being fixed and is expected. 82ddd3e,
  dropping the aryl carbon, is the <1 % above. Read as one change, the large
  shift becomes evidence about the aryl-carbon exclusion, which it is not.

  Second, a `blind_median_clearance_nm` anywhere near 0.233 in a box believed
  to be at this density is unambiguous evidence the box is wrong, since that
  value is only legitimate in the sparse regime.

  Third, the default sits *above* the blind median: 0.150 against 0.143, which
  is the same fact as only ~45 % of caps keeping the O-H vector. So at melt
  density roughly half of all cap sites cannot reach `cap_min_clearance` along
  the preferred direction and the search is doing real work rather than
  rubber-stamping. That makes 0.150 demanding but not absurd -- which is more
  than it was when it was unreachable by construction -- and it suggests the
  recalibration may move it *down* rather than up. Not the direction either
  session guessed, and worth not being surprised by.

  Independent check on the scale, by the study session: uncorrelated Poisson
  at 53/nm^3 has a median nearest-neighbour distance of 0.1462 nm from a probe
  point, dropping to 0.1236 when minimized over both atoms of a -C#N group at
  0.116 nm separation. The swept 0.143 sits inside [0.124, 0.146], near the
  top, which is where the excluded hole around the attachment oxygen should
  put it. Script at
  `~/devtests/htpolynet/bridge-series/verify_blind_clearance_scale.py`.

  A trap in that check, flagged by the study session and worth disarming
  before someone rederives it: the disjoint-sphere floor,
  `(ln 2 / (8/3 pi rho))^(1/3)`, evaluates to 0.116 nm at rho = 53, which is
  also `cn_len`, the C-N bond length. That is arithmetic coincidence at this
  one density and nothing else. The floor tracks number density -- 0.127 at
  40/nm^3, 0.116 at 53, 0.111 at 60, 0.101 at 80 -- while the bond length is
  fixed. Anyone who writes down "the floor is the C-N bond length" has made a
  claim that is invisible here and wrong everywhere else, and has attributed a
  number-density result to cap geometry.

  Two caveats on those figures. The sweep places uncorrelated Poisson points
  with a hole around the oxygen, which is not a melt -- real packing is
  correlated and has excluded volume between the neighbours too -- so it gives
  a scale, not a prediction, and in particular its `n_below_target` of zero is
  certainly optimistic. And replicate scatter in this system is topological
  rather than numerical: the study session measures density sd rising from
  0.71 kg/m^3 at chi_OCN 0 to 3.25 at 0.746, a factor of 4.6, because at zero
  conversion every build is the same molecular liquid. Cap sites are network
  sites, so expect the scatter in `blind_median_clearance_nm` to grow with
  conversion the same way, and size any series meant to calibrate the default
  accordingly rather than assuming the high-conversion points are as tight as
  the low ones.

  **Measured on real boxes, 2026-08-27** (study session, job 22107390): 14
  independent cured BPA-cyanate-ester boxes at v2.6.1, chi_bond 0.740--0.901,
  1955 caps, ~53 heavy atoms/nm^3 -- the first placement numbers taken without
  the clearance bug present. `blind_median_clearance_nm` 0.1198 +/- 0.0061
  against the synthetic sweep's 0.143: real boxes are ~16 % *tighter*, in the
  direction the sweep's own caveat predicted, so the sweep's absolute values
  should not be used to set the default. `blind_min_clearance_nm` 0.006--0.048.
  722 of 1955 caps (37 %) would have been placed inside 0.10 nm blind.
  `n_preferred_out_of_angle` 109 (6 %). `n_below_target` 1 of 1955, against
  the sweep's 0 that was flagged "certainly optimistic".

  And a correction to this entry's own advice, which the numbers force:
  `min_clearance_nm` is **not** a tail statistic and must not be calibrated
  against. Across all 14 boxes it came in at 0.1503 +/- 0.0006 against a 0.150
  target -- the target seen from above, not a measurement. The search exits at
  the first direction reaching target, so whenever it succeeds for every cap
  the worst-placed cap is one that only just cleared, and the minimum is
  pinned to the threshold by construction. The one box below it (cb0773,
  0.1488) is the one box with `n_below_target` = 1. This is the searched-median
  defect one statistic further along; it was missed because the median's
  version of it was the one being written up. Docs corrected.

  So the recalibration is unblocked and wants the `blind_*` columns. Note that
  0.150 sits well above the measured blind median of 0.120, not marginally
  above 0.143 as the sweep suggested -- so the search runs for well over half
  of all caps -- yet it succeeds for 1954 of 1955. Demanding and reachable at
  once, which is a different situation from the one this entry was written
  against, and it weakens the earlier guess that recalibration would move the
  default down.

  **The placement change is density-neutral at high conversion, 2026-08-28**
  (study session; picotte 22107390 complete, 14/14, exit 0:0). v2.6.1 anchors
  at chi_OCN 0.7438 give 1197.59 kg/m^3 (n = 2) against v2.3.0's 1197.94 +/-
  3.25 at 0.7458 (n = 8): offset -0.34 on a combined SE of 2.38, i.e. 0.1
  sigma, bounding |offset| < 4.7 kg/m^3 at 2 sigma. So the net of 5a58d3d and
  82ddd3e does not move bulk density where it has been checked. Two limits.
  The bound is small against the 21 kg/m^3 rise across the series, so the two
  versions pool for curve shape -- but it is *not* small against the 3.0 kg/m^3
  basin depth, so it licenses nothing about the basin. And it does not
  calibrate `cap_min_clearance`: a placement change being invisible in bulk
  density is a much weaker statement than the threshold being right, and this
  item stays open exactly as written above.

- **The cap direction search stops at the first adequate direction, not the
  best one.** `_choose_cap_placement` breaks out as soon as a direction
  reaches `cap_min_clearance`, so a cap that could have had 0.25 nm of room
  may be left with 0.15. Two consequences. The placement is worse than it
  needs to be, for compute that is genuinely negligible -- the angle window
  admits about 43 % of a 48-direction Fibonacci spiral, so a full scan is
  ~21 vectorized distance evaluations against a neighbour list of order 100.
  And `median_clearance_nm` is partly a readout of the threshold rather than
  of the box, since the distribution is truncated from below at the target --
  as is `min_clearance_nm`, which is pinned *at* the target rather than merely
  pulled toward it, measured at 0.1503 +/- 0.0006 across 14 real boxes.
  Removing the early exit would give both fields their meaning back.

  Deliberately not changed for 2.6.1: taking the best direction moves every
  cap that needed a search, and doing that on the same release as the
  clearance fix means the next real build cannot attribute a change in the
  numbers to either one. That precondition is now satisfied -- the 14-box
  measurement of 2026-08-27 recorded above is the clean pre-change baseline,
  so this is ready to do. Note that changing it invalidates the achieved
  columns of anything built before it while leaving the `blind_*` columns
  comparable, which is another reason those are the ones to correlate on.

- **The shortfall from the cube law between `f` and many iterations is
  unexplained, and is now known to be narrow.** Below `f` iterations
  crosslinker conversion is exactly zero, a counting constraint; the four runs
  of the nine-iteration cohort land +0.6 to +3.5 % above the cube (attested at
  nine for two of the four; the other two are the same cohort at the same
  target with the count unrecorded). In between, measured ratios
  to the cube law are 0.06 and 0.46--0.50 at three iterations and 0.86 at
  four -- an eightfold spread at an *identical* iteration count, so the
  iteration count does not determine it. The band above that is now measured
  and shows no shortfall at all: 14 runs at `n`/`f` 1.67--2.67 sit at
  1.020 +/- 0.029 of the cube, with 3 of the 14 *below* it to -3.7 % against
  replicate scatter of 0.8--1.8 %. So the unexplained region is only
  `f` to about 1.7`f`, and across 1.7--2.7`f` the cube is a two-sided estimate
  rather than the floor the docs used to call it -- corrected there.

  What that does **not** settle is 3`f`, and the temptation to settle it by
  arithmetic should be resisted: applying the 1.7--2.7`f` scatter to the
  nine-iteration cohort gives P(below cube) = 0.245 and makes four-of-four
  above a one-in-three outcome, but that presupposes the deviation is
  regime-independent, which is the claim in question. The cohort on its own is
  +2.64 % sd 1.39 % (n = 4), giving P = 0.03 -- an order of magnitude apart --
  and at n = 4 the sd's 95 % CI is [0.79, 5.18], so the two cannot be
  distinguished. 3`f` is untested, not disproved. Settling it wants more runs
  at 3`f`, not a wider distribution borrowed from a lower one. Whatever
  mechanism explains the shortfall has to switch off within a factor of two
  of `f`.

  The one-bond-per-residue-per-iteration rule was simulated against this and
  moves the n=3 prediction only from 30.0 to 25.9 against 15 observed: a fifth
  of the gap. The leading untested guess is spatial anti-correlation --
  bonds forming preferentially on crosslinkers that already have room, leaving
  the rest to compete. This matters because it is the regime where a user's
  reported conversion is most wrong, and because a warning better than the
  current `iterations < f` one needs a mechanism to threshold on. Widening
  that warning to `n < 2f` was considered and rejected: no threshold on `n`
  alone can work.

  That conclusion survives a correction to its own evidence, and comes out
  stronger. The 8x pair was **confounded**: conv40 r3 and conv50 r1 are both
  n = 3 but sit at chi_bond 0.40 and 0.50, so they never isolated `n` either.
  The new low-conversion builds split it -- at fixed chi_bond the residual
  spread is 1.36x, not 8x, and the chi_bond 0.40 ratio-to-cube reproduces to
  three digits across v2.3.0 and v2.6.1 five months apart. So at fixed `n` the
  ratio moves 7--10x purely with chi_bond, which is a sharper refutation of an
  `n` threshold than the confounded pair was. What does *not* survive is the
  explanation: the "product of per-iteration bond fractions" account is
  superseded by plain chi_bond, which is simpler and reproduces across
  versions. Docs corrected.

  **Reinstated 2026-08-29**, and worth recording because it was withdrawn and
  reinstated inside two days. The "product of per-iteration bond fractions"
  account was not superseded -- it was under-evidenced, and the iteration logs
  now make it the best-evidenced thing in this entry. Bonds per iteration decay
  steeply (bpa-fl0500: 156, 145, 59; bpa-cb0900: 162, 148, 128, 92, 46, 29, 17,
  26), so `720 x chi_bond / n` is an average no build actually realises. What
  matters is the LAST iteration, the one that has to supply a triazine's third
  bond: bpa-fl0400 spent 9 bonds there against an average of 96 and came in at
  ratio 0.064, and conv40 r3 -- same conversion, v2.3.0, five months earlier --
  spent 10 on its own last iteration and also came in at 0.064.

  So chi_bond is a **proxy for the distribution, not a replacement for it**.
  The algebra above still holds for *average* bonds-per-iteration, which is
  exactly collinear with chi_bond at fixed n (r = 1.0000). The *last-iteration*
  count is a third variable at r = 0.980 with chi_bond -- separable in
  principle, not on four points. The withdrawal on 08-28 was made on the
  strength of a correlation, which is the move this entry already warns about
  one paragraph earlier.

  **Measured below the band, 2026-08-28** (study session; picotte 22132281,
  12 BPA builds at v2.6.1, cures final, ladders still running): at bond
  conversions of 0.40--0.73 the deviation from the cube is one-sided and
  large -- 11 of 12 below, mean -19.9 %, worsening monotonically as bond
  conversion falls, -93.6 % at 0.40. The single build inside the documented
  band (0.7306, +1.5 %) agrees with it, so nothing here disturbs the
  1.7--2.7`f` scoping.

  The open question this creates is which variable is doing the work. The
  documented band is *both* 1.7--2.7`f` and bond conversion 0.74--0.90; this
  array varies bond conversion and its iteration counts are not yet in hand
  (they live in `diagnostics.log`, returned only when a task ends -- expected
  early 2026-08-29). So the docs state the degradation against bond
  conversion, which is what was measured, and do not attribute it to `n`/`f`.
  **Resolved 2026-08-29**, counts returned with the array. Four builds sit at
  exactly n = 3 spanning 15x in ratio-to-cube (0.064 / 0.633 / 0.795 / 0.937 at
  chi_bond 0.40 / 0.50 / 0.55 / 0.58), monotone in chi_bond. So "iteration
  count does not determine the shortfall" is established directly, on four
  same-n runs, rather than inferred from the confounded pair. Over all 26
  builds r(ratio, chi_bond) = +0.813 against r(ratio, n/f) = +0.567, with
  r(chi_bond, n/f) = +0.905; counts run 3 to 8, n/f 1.00 to 2.67. Docs cite
  the four-run evidence now.

  **The deviation changes sign, 2026-09-03** (study session, 54 builds
  audited at v2.6.2): near `chi_bond` 0.900 the observed crosslinker count
  runs *above* the cube -- 180.60 against 175.39 predicted over 30 builds,
  +2.17 pp at t = +10.8 against a replicate sd of 2.63 triazines -- while at
  `chi_bond` <= 0.64 it runs below, 32.17 against 41.28 over 6 builds, -3.80
  pp. The compact way to say it is that the effective exponent
  `ln(chi_OCN)/ln(chi_bond)` is not 3 and is not constant: ~3.5 near
  `chi_bond` 0.55, ~3.0 near 0.73, ~2.6--2.7 near 0.90. No single power law
  fits, so the cube is a local estimate around 0.9 and nothing more. This is
  consistent with the 1.020 +/- 0.029 above, and puts a significance on it.

- **A config-time version of that warning, before any compute is spent.**
  Largely superseded by the completed-crosslinker check that shipped -- nothing
  needs a proxy for a quantity that is exactly known by the time the cure
  ends -- and worth keeping only for
  what the post-hoc check cannot do: warn *before* a multi-hour build rather
  than after it. That is a weaker claim on much worse evidence, so it is a
  separate feature and not a substitute. The reasoning that motivated it:
  chi_bond predicts the ratio-to-cube to within 1.36x, and unlike a bond
  distribution -- which does not exist until the cure runs -- it is knowable
  before the cure starts. So
  htpolynet could warn at config time that a run will yield far less
  crosslinker conversion than the cube law implies, which is the number a user
  is actually reasoning with. The existing `iterations < f` warning fires only
  after the fact and only in the exactly zero regime.

  **Do not key it on `desired_conversion`, and do not key it on chi_bond
  either.** chi_bond is a config-time proxy, not the governing quantity, and
  the collinearity is not an accident of sampling that more builds would fix.
  Within a fixed iteration count `bonds = 720 x chi_bond` for this system, so
  bonds-per-iteration `= 720 x chi_bond / n` is *exactly* proportional to
  chi_bond. No dataset at fixed `n` can separate them, however large. Breaking
  the proportionality needs `max_conversion_per_iteration` or
  `min_bonds_per_iteration` varied AT FIXED `desired_conversion` -- which is
  a deliberate experiment nobody has run, not something a wider sweep of
  conversions will deliver.

  There is a further reason to keep it qualitative: the quantity that appears
  to govern the shortfall is the number of bonds formed in the *last*
  iteration, and that is not knowable at config time at all -- it is an
  outcome of the cure, not a setting. Any config-time warning is necessarily
  keyed on a proxy for it.

  Those same two directives are also why `desired_conversion` is the wrong
  key in ordinary use: they set bonds-per-iteration independently of it, which
  is exactly what the warning at `cure/curecontroller.py:674` already tells a
  user to do to spend more iterations at the same conversion, and what this
  project plans in order to reach 9+ iterations. A threshold keyed on
  `desired_conversion` alone would misfire on precisely the runs that took
  htpolynet's own advice.

  So keep it qualitative until something separates the two. A warning saying
  "below about 0.7 the cube law overstates badly, and worsens as you go lower"
  needs no attribution and could ship soon. One that quotes a number needs both
  the separation and a sense of how much of the curve is BADCy-specific -- it
  is one chemistry, one force field, one cure protocol, n = 1 per target apart
  from a single duplicate.

- **`completion_bias` biases the `B` side only, and that is a convention,
  not a law.** The new `CURE.controls.completion_bias` ranks bond candidates
  by how many bonds their `B`-side residue already carries, because
  htpolynet's A2+B3 idiom puts the multifunctional crosslinker in the `B`
  position -- as example 6 does. A user who declares the crosslinker as `A`
  gets the bias pointed at their *bridges* instead: not a no-op, a different
  claim about which partly-reacted species is a reactive intermediate. That
  case is now warned about at the first cure iteration, by comparing the
  initial functionality of the two sides, so it is visible rather than silent
  -- but it is still not *served*. The generalization is small (rank on the
  `A` side, or on the sum of both) and each variant is a different physical
  claim, none of them run. Do it when someone has a chemistry that needs it,
  and make them say which side rather than inferring it from functionality:
  the warning's heuristic is good enough to flag a probable mistake and not
  good enough to silently redirect the ranking on.

- **Nobody has run example 6 with `completion_bias` on.** The ranking is
  unit-tested -- ordering, tie-breaks, missing residues, the fallback when a
  `.grx` predates the `nreactions` attribute -- but the acceptance criteria
  that matter are system-scale and need gmx and a multi-hour build:
  incomplete triazines should collapse from tens to ~1 at
  `desired_conversion: 0.90`, crosslinker conversion should rise from ~0.76
  to ~0.90, `TAZ + CYN/3` must still equal the initial triazine count
  exactly, and atom conservation must still be exact. If the incomplete
  count does *not* collapse, the sort key is not surviving as far as the
  truncation and that is where to look. Until someone runs it, the directive
  is documented as untested at scale.

- **Cap placement is greedy and never backtracks.** Transferred `-C#N`
  fragments are now placed against the local neighbourhood rather than blindly
  along the old O-H vector, which removes the catastrophic overlaps: in a
  synthetic box at polymer density, blind placement put 91 of 158 caps within
  0.10 nm of a neighbour and the worst at 0.008 nm, while the neighbourhood
  search puts none below 0.128 nm. What it does not do is relax. Each cap is
  placed once, in the clearest direction available *at that moment*, and later
  caps must work around it; in a genuinely over-packed box a residue of
  sub-target placements remains (17 of 158, 33 of 378 in the same synthetic
  test) and is reported rather than fixed. The next lever is minimizing
  incrementally as caps land, or placing in order of how constrained each site
  is rather than in match order. Do it if the reported clearance warnings turn
  out to correlate with builds that still die -- the instrumentation to decide
  that now exists, and did not before. One thing the instrumentation still
  does not report is how far a transferred cap actually travelled:
  `_greedy_match` can reach the globally nearest free oxygen after ten radius
  doublings, and one audited build transferred 72 of 177 caps, but
  `repair-summary.yaml` carries clearance statistics only, so an external
  audit could not bound the transfer distance. Report it alongside the
  clearances.

- **`bdf.loc[:abs_max]` takes one bond more than the limit.** In
  `curecontroller.py::_searchbonds`, the truncation that applies
  `max_conversion_per_iteration` uses `.loc` with a slice, which is
  inclusive of its endpoint, so an iteration limited to `n` bonds forms
  `n + 1`. Harmless in practice -- the limit is a throttle, not a
  correctness bound -- but it is off by one, and fixing it changes the
  trajectory of every existing config by one bond per throttled iteration.
  Worth doing at a version boundary where a small reproducibility break is
  already expected, not before.

- **An algebraic self-check that no build currently runs.** For an A2 monomer
  of `A` atoms in the 360-bisphenol/240-triazine example, the total atom count
  is `360*A + 2160 - 1440` regardless of conversion -- verified at 12600
  (bpa, A = 33), 9720 (bpo, A = 25) and 14760 (tmb, A = 39). It is
  conversion-independent and closed-form, so any atom dropped or duplicated
  anywhere in cure or repair breaks it. Cheap to assert at the end of repair.

- **`_distance_attenuation` runs an open-loop stage ladder.** In
  `cure/curecontroller.py:606`, `this_nstages = int(maxL/d['increment'])` is
  computed once from the initial maximum bond length, and the loop then runs
  that many stages whatever happens. `maxL` *is* recomputed inside the loop,
  but only to log it -- it never extends, shortens or aborts the ladder -- and
  `restore_bond_parameters(saveT)` runs each stage whether or not the bonds
  actually closed. So a stage schedule that fails to bring a bond in has no
  way to say so. Robustness only; it was ruled out as the cause of the
  halogen failure below (quadrupling the stage count made that *worse*).
  The density-convergence entry below is the same shape of problem -- a
  number computed once and trusted thereafter -- and if both are done, they
  are one change to the same loop.

- **`CURE.relax` should gate on measured density convergence, not on a fixed
  step count.** Raised by Cameron 2026-09-05, from the observation that
  pestifer has a gate system htpolynet does not.

  Step zero is instrumentation, and it is missing: `_do_relax` delegates to
  `_distance_attenuation`, which never calls `TopoCoord.equilibrate()` --
  the only method that runs `gmx_energy_trace(..., ['Density'])`. So the
  relax stages **do not observe density at all**, even though each one
  already writes an `.edr` that contains it. Density is traced only in
  `_do_equilibrate`, and that stage defaults to 300 K, below `Tg`, where the
  box cannot densify.

  The window is small. The default relax equilibration sequence ends in an
  NPT segment of 2000 steps, and `relax-npt.mdp` runs `dt = 0.001`, so that
  is 2 ps per stage; at the default `nstages: 6` over ~10 iterations it is
  roughly **120 ps of above-`Tg` constant-pressure time across an entire
  cure**, against a 39.2 ns production ladder.

  Measured (study session, 4 replicates per chemistry, slope over the second
  half of cure where the topology is nearly final): density is still climbing
  when cure terminates -- bpa +1.29 +/- 0.62 and bpf +1.48 +/- 0.28 kg/m3 per
  iteration, 2.1 and 5.3 sigma. Cure stops on a conversion criterion with no
  reference to whether the box settled. The honest caveat is that a rising
  density is *partly expected*, since crosslinking genuinely densifies, and
  this slope cannot separate real densification from incomplete relaxation.
  That inability is the argument: without an observed criterion there is no
  way to tell which builds settled.

  Why it is more than an equilibration nicety: if the box is under-relaxed
  while bonds are forming, bonds form between whichever atoms happen to be
  adjacent in a too-open configuration, and **the topology is then
  permanent** -- a later remelt relaxes coordinates but cannot rewire bonds.
  That one mechanism ties together three otherwise-awkward results from the
  study's 64-build series: post-cure annealing is refuted as a cause of the
  density deficit (the ladder opens 2 ns at `Tg`+112 K, erasing pre-ladder
  history), the uncured monomer melt matches experimental dilatometry to
  0.19 % while the cured network is 2.25 % low, and cure shrinkage is
  strongly bridge-dependent (+0.033 to -0.004 ml/g) and anti-correlates with
  monomer van der Waals volume -- which is what a *fixed* relaxation window
  does when different chemistries have different relaxation times. A gate
  fixes the bridge dependence for free, because bulkier networks simply get
  more time.

  **Scope it to relax.** The production ladder does not need this: two
  structures with opposite histories held 5 ns at 480 K converge from
  opposite sides and stall 4.25 kg/m3 apart with the ladder's own hold inside
  that bracket; within-hold drift on the glassy branch is -0.35 to +0.53
  kg/m3 per 667 ps across four replicates; and replicate scatter is
  topological rather than equilibrative (sd 3.2 kg/m3 at high conversion
  against 0.7 at `chi_OCN` 0).

  What to port from `pestifer/util/density_convergence.py` is only the
  criterion -- an **autocorrelation-corrected SEM**, `sigma/sqrt(N/tau_int)`,
  because NPT cell density is autocorrelated over hundreds of steps and a
  naive block-means SEM is optimistic by ~1.5x; this also makes the gate
  size-aware for free, since `sigma/mean ~ 1/sqrt(N_atoms)` while tau is
  roughly size-independent -- plus an explicit **ceiling outcome**, so a
  build that never settled says so instead of silently reporting a density.
  Do **not** port the chunking (`next_chunk_steps`, `is_patch_grid_crash`):
  that exists because NAMD fixes its patch/PME grid at the start of each
  `run`, and GROMACS rescales the box within one `mdrun` without that
  failure mode. The chunking is most of pestifer's complexity and none of
  its value here.

  Minimal design, in shippable order: (1) read Density from each relax NPT
  `.edr` and log it -- worth shipping on its own, since it makes the question
  answerable for every existing user; (2) test the trailing window after each
  relax stage and extend or repeat the last stage until it converges or hits
  a ceiling; (3) record converged/ceiling plus residual drift in
  `diagnostics.log` and the repair summary yaml. One optional block under
  `CURE.relax`, defaulting to off, so existing results stay reproducible.

  **Do not build the gate until the crude arm reports.** The study has a run
  testing `CURE.relax` NPT `nsteps` 2000 -> 32000 (~120 ps -> ~1.9 ns of
  above-`Tg` relaxation) on bpa, which reaches 47 % of experimental cure
  shrinkage, and bpf, which already matches experiment and is therefore the
  control that should *not* move. If bpa gains and bpf does not, the
  mechanism above is confirmed. If both move, or bpf moves more, this whole
  account is wrong and the gate would be solving the wrong problem.

  **Steps 2-3 have a published competitor, and it may be the better
  design.** Amended 2026-09-05 by the study session, which found it after
  sending the note above. Schichtel & Chattopadhyay 2020, *Comput. Mater.
  Sci.*, gate *bonding* on an Arrhenius reaction probability built from cure
  temperature, cutoff distance and an activation energy -- kinetics-gated,
  not density-gated. That attacks the cause, a bonding rate outrunning
  relaxation, where a density gate monitors a symptom; it is also the
  physically motivated knob, since the thing a real cure has and this one
  does not is a reaction rate. Their sec 2.4 reportedly *demonstrates* the
  failure mode described above rather than merely flagging it -- a density
  trajectory they call "not physical", ending in a metastable configuration,
  caused by accelerated reaction kinetics. That is the same mechanism as the
  permanent-topology argument, arrived at independently and published, which
  is worth more than the argument on its own.

  Nobody here has read the paper yet; this is the study session's report of
  it, and the citation should be checked before it is repeated anywhere
  public. Read it before designing either gate -- the two are not exclusive,
  and a kinetics gate would change what a density gate is even for. Step 1,
  observe and log, is unaffected either way and still worth shipping alone.

## Simulation defaults

- **The halogen constraint failure is fixed but not explained.** v2.7.0
  switched the per-iteration equilibration to `constraints = h-bonds` with
  `lincs_order = 8`, which took a fluorinated bisphenol from 0 of 7 builds to
  4 of 4, and the eight-bridge series from 7 failures to 32 of 32. What is
  still not understood is the *timing*. The constraint artifact is present
  from the very first NPT step -- a melt of pristine, uncrosslinked BAF
  monomers at production density shows step-0 pressure of -1.34e5 bar, with
  no cure involved at all -- yet the old builds survived densification,
  precure and 4-6 cure iterations before dying. So it was a predisposing
  cause plus a threshold in the network's tolerance for it, not a single
  trigger, and nobody knows what sets the threshold. The failing structures
  were not preserved, so the LINCS atom indices were never mapped to
  residues; anyone chasing this should trap and keep them first. The
  diagnostic that identified it is worth reusing: halving `dt` made the
  artifact ~4x *worse*, and 1/dt^2 scaling is the signature of a one-shot
  constraint start-up projection rather than a physical clash.

- **`mdp_to_dict` crashes on any mdp line containing two `=` signs.** In
  `external/gromacs.py:239` it does `k,v = l.split('=')`, which raises
  `ValueError: too many values to unpack` rather than reporting anything
  useful. No packaged template trips it today -- checked across all eleven --
  but any comment containing an equals sign does, including the obvious one
  someone would write next to a timestep. It also does not strip `;`
  comments, so the comment text becomes part of the value and is re-emitted
  attached to it; GROMACS tolerates that, but `mdp_get` on a commented
  numeric key would not survive a `float()`. Split on the first `=` only,
  and strip comments.

- **`postsim` inherits `npt.mdp` wholesale, so cure defaults govern the
  production measurement.** `analysis/postsim.py:81` checks out `npt` and
  renames it; `build_mdp()` overrides `ref_t`, `ref_p`, `nsteps`, velocity
  generation and `tcoupl` (to `v-rescale`) but never touches `dt`,
  `constraints`, `lincs_*` or `pcoupl`. So whatever is chosen above also sets
  the timestep and constraint scheme of every postsim run, and the barostat
  there is still **Berendsen**, which does not sample a correct NPT ensemble
  and is deprecated in modern GROMACS. Changing the barostat is a physics
  change to published numbers; changing it silently would be worse than
  leaving it, so it wants a release note and probably a config knob.

## Usability

- **The atom-serial column in `final.gro` wraps modulo 10000, and nothing
  says so.** That is the GROMACS `.gro` format, not an htpolynet bug -- atom
  10000 prints as `0` -- but every htpolynet system large enough to matter
  crosses it, and the natural analysis script joins `final.gro` to `final.top`
  on that column. In one 12,600-atom build 2,601 atoms (21%) have a `.gro`
  serial that differs from their `.top` index, so such a script is silently
  wrong for a fifth of the box with no error anywhere. Worth a sentence in the
  analysis docs saying to join on line order, not on the serial. Same family
  as the existing trap that `min_clearance_nm` is a substring of
  `blind_min_clearance_nm` while PyYAML emits keys alphabetically, so a naive
  grep of `repair-summary.yaml` picks the wrong one.

- **`gen-slurm-script` doesn't stage to scratch.** The emitted script
  runs in the submit directory. A cure run does heavy small-file I/O
  every iteration, so on a cluster whose home and group storage are NFS
  that is the wrong place. The cluster-correct pattern is to run in
  node-local or parallel scratch and copy results back, with a `trap` so
  partial results survive a timeout. It also cannot emit a job array,
  which is the right shape for sweeping a set of configs.
- **`-restart` is documented as "EXPERIMENTAL: broken at the cure
  stage".** A build that dies mid-cure currently has to start over — the
  worst possible time to lose work, since cure is the longest stage.
- **Make `input-check` a real config linter.** It currently reports only
  the initial atom count. Everything a new config gets wrong is checkable
  cheaply and statically, before the user spends hours discovering it:
  (a) that each `symmetry_equivalent_atoms` group really is topologically
  equivalent — RDKit canonical ranks with `breakTies=False` settle it in
  about ten lines, and a wrong group silently generates reaction templates
  the cure stage will never match; (b) that the A2+B3 site counts actually
  balance at the given monomer `count`s, and that `desired_conversion` is
  reachable given them; (c) that no `reactive_atoms` name collides with the
  bare element names `_reset_names_to_element` assigns to every unmapped
  atom, since the repair stage looks atoms up by name; (d) that the
  requested `charge_method` matches the provenance record of any cached
  parameterization of the same name -- the build itself now rejects a
  mismatch, but reporting it in the pre-flight is cheaper than discovering
  mid-run that half the molecules need re-parameterizing. Building the
  cyanate-ester bridge series meant doing (a) and (b) by hand in a throwaway
  RDKit script, which is exactly the work a user should not have to
  reinvent.
- **Generated topologies cannot be compared byte-wise, because ParmEd
  stamps them.** Every `.top` htpolynet writes opens with a ParmEd header
  recording the invoking user, the host, and the date:

      ;   File TAZ.top  was generated
      ;   By user: cfa (1000)
      ;   On host: panacea.chemeng.drexel.edu
      ;   At date: Fri. May  5 15:55:38 2026

  So two physically identical parameterizations never hash the same, and the
  obvious check -- "did these two builds produce the same parameters for the
  shared monomers?" -- returns a false difference. The calibration study hit
  this on its bridge series and nearly reported the series confounded on the
  strength of an md5; the fix on their side was to strip comment lines before
  comparing, which every user will have to reinvent. Note the header is
  itself provenance, but the wrong kind: it records *when and where* rather
  than *what directives*, which is what the `parm` record now covers, and it
  actively prevents the comparison you would want. The `.itp` files carry no
  such header and the `parm` records are JSON with sorted keys, so both
  already compare cleanly. Worth either normalizing the `.top` header away or
  shipping a comparison helper; this belongs with the manifest entry below,
  since both are about being able to answer "is this the same build?".

- **`-lib` is read-only, and nothing says so.** `pfs.checkout()` and
  `pfs.exists()` consult the `-lib` user library first, then the user cache,
  then the system library -- but `pfs.checkin()` writes unconditionally to the
  user cache, and `UserLibrary` has no `checkin()` method at all. So pointing
  a run at `-lib somewhere` redirects where it *reads* parameterizations from
  while leaving where it *writes* them untouched. That asymmetry is invisible
  from the flag, from `--help`, and from the docs, and it is a natural thing
  to get wrong: a user who passes `-lib` for provenance reasonably concludes
  their products are being contained there. Either `checkin()` should prefer
  the user library when one is configured, or the flag and docs should say
  plainly that it governs lookup only. Verified empirically 2026-08-23: with
  `userlibrary` set, `pfs.checkin()` still landed the file in the user cache.

  What makes this more than a naming problem is that the writes are silent
  and cumulative. Nothing logs them, nothing warns, and a run that is
  quietly depositing molecules into a shared library looks identical to one
  that is not -- the only way to discover it is to audit mtimes, which
  nobody does unprompted. The calibration study ran three sweeps believing
  `-lib` contained it, and would have added roughly seventy entries to
  `~/.htpolynet` across its remaining bridges; it was caught by reading
  `checkin()`, not by anything the tool said. Whatever the fix, a run should
  be able to say where its parameterizations went -- which is the same
  requirement as the build manifest below, approached from the write side
  rather than the read side. The two may well be one item: *where did this
  come from, and what produced it*.

- **No way to run without writing to the user library.** `pfs.checkin()`
  declines to *replace* an entry unless `--force-checkin` is given, but it
  always *writes* one the library does not yet hold, so any run adds every
  molecule name it produces to `~/.htpolynet`. The flag's help text said
  "force check-in of generated parameter files to the system library", which
  reads as though check-in happens only with the flag, and it named the wrong
  library; that wording is fixed, but the missing capability is real. The
  calibration study wanted exactly this: a way to develop against a config
  without its intermediate products accumulating in a library shared with
  other work. Today the only lever is pointing `HTPOLYNET_CACHE` elsewhere,
  which is a blunt instrument because it also hides the entries you *do* want
  to reuse. A `--no-checkin` flag threaded through
  `Runtime._checkin_parameterization()` would cover it in a few lines.

- **Two different defaults for `charge_method`.** `AMBERTOOLS_DEFAULTS` in
  `external/ambertools.py` says `bcc`, which is what a direct
  `GAFFParameterize()` call with no directives gets; `Runtime.runtime_defaults`
  says `gas`, which is what every actual build gets, because
  `_apply_runtime_defaults()` fills it in before AmberTools is ever reached.
  Both were already there and the provenance record is consistent either way
  -- whichever default applies is the one recorded -- so this is a
  readability trap rather than a live bug. It is still worth collapsing to
  one value, and the answer is probably `gas`, since changing what builds
  default to would silently change everyone's charges, which is the exact
  class of harm the record was added to prevent.

- **Nothing durable records that a build reused parameterizations of
  unverified provenance.** The parameterization stage now warns per molecule
  and again in a block at the end of the stage, but both live only in the
  log, and a log is the first thing discarded. Someone reading a result six
  months later -- or a reviewer asking what a published network was actually
  parameterized with -- has no artifact to check. htpolynet writes no build
  manifest today; `profile.json` holds timings only, and the diagnostic log
  is the log. A small `build-manifest.json` in the project directory listing
  each molecule, its origin (`newly parameterized` / `previously
  parameterized`), and its provenance record where one exists would carry
  that, and would be worth more than this one flag: it is also the natural
  home for the config hash, the htpolynet and AmberTools versions, and the
  seed once seeds exist. Raised by the calibration study, which caught the
  original cache bug from a wall-clock anomaly rather than from any log line.
  See also the `-lib` entry above: "a run should be able to say where its
  parameterizations went" and "a result should be able to say what produced
  it" are the same requirement from two directions, and a manifest that
  recorded the check-in destination would answer both.

- **A cached parameterization is not checked against the input structure
  it was built from.** The provenance record added alongside the cache now
  covers `charge_method`, `net_charge` and `atom_type` -- everything in the
  AmberTools invocation -- but not the structure antechamber consumed. So
  editing `lib/molecules/inputs/TAZ.mol2` and re-running still reuses the
  parameterization of the *old* geometry, silently, exactly as the charge
  method used to. Hashing the input would close it, and for a monomer that
  is easy: `Molecule.parameterize()` has the input file in the working
  directory at the moment it runs. Two things stopped it going in with the
  rest: (a) a molecule built by a reaction has no stable input to hash --
  `generate()` writes its mol2 from the merged reactant TopoCoord, whose
  coordinates vary run to run, so hashing it would make every generated
  molecule a permanent cache miss; the hash would have to be recorded only
  for `origin == 'unparameterized'` monomers and compared only when both
  sides carry one. (b) `Runtime._cached_parameterization_mismatch()` runs
  before `generate()` checks the input structure out, and `pfs` has no way
  to resolve a library file to an absolute path without copying it into the
  working directory -- `checkout()` always copies. A `pfs.locate(filename)`
  returning the resolved source path across user library, user cache and
  system library is the missing piece, and is worth having on its own.

- **No seed control anywhere, so a build cannot be reproduced.** Three
  independent sources of randomness are all unseeded: `random.sample` for
  conformer selection (`core/runtime.py`), `np.random.random()` for the
  per-bond probability test (`cure/curecontroller.py`), and every `.mdp`
  sets `gen-vel = yes` without `gen-seed`, so Gromacs picks a pseudo-random
  one per run. Two runs of one config therefore diverge — example 6 gave 57
  vs 60 incomplete triazines on two machines at the same 0.90 conversion.
  That is convenient in one direction (independent replicas of a
  quenched-disorder ensemble come free by re-running) but it means a build
  reported in a paper cannot be reproduced exactly, and a failure seen once
  may not reappear. The fix is a top-level `seed:` that feeds all three:
  seed `random` and `np.random` at startup and write `gen-seed` into every
  generated mdp, with replicas then requested by varying it rather than by
  relying on entropy.
- **Drop the pre-3.5 matplotlib fallback** in `analysis/plot.py`'s
  `_get_cmap()` once `matplotlib>=3.6` is a safe floor; the
  `matplotlib.colormaps` registry is then always present.

## Example depot

- **Example 6's postcure anneal peaks too close to *T*:sub:`g` to relax the
  network, and the fix is not simply "run it longer".** Measured on four
  independent BPA builds: the postcure NPT plateau gives 1.1712 ± 0.0029
  g/cm³ where the same systems melted and slowly re-cooled give 1.1983 ±
  0.0037 g/cm³ -- 2.31% apart, with the plateau 2.6-2.8% below experiment
  and the re-cooled value within 1%. The plateau is an under-relaxed
  structure.

  The mechanism, and the reason the obvious fix is the wrong one: the
  example anneals at a 500 K peak, and the measured *T*:sub:`g` for this
  system is 487.9 K. That is 12.1 K above the glass transition, for 80 ps.
  The Tg ladder that produced the relaxed value spent 15,785 ps above
  *T*:sub:`g`, peaking 112 K above it -- 197x the time, in the regime where
  a crosslinked network actually moves. Lengthening an anneal that sits 12 K
  above *T*:sub:`g` extends a process that barely moves anything, so the
  lever is the peak temperature, not the duration.

  **But raising the peak is unlikely to close the gap entirely.** Two
  structures of the same network were held for 5 ns at a 480 K setpoint
  that thermostatted to 477.1 K -- about 11 K *below* the measured
  *T*:sub:`g`, and far longer than the 80 ps the example spends near it.
  Relaxation there is fast at first and then stops: about 87% of an
  initial 34.8 kg/m³ density gap -- the gap measured over the first 200 ps
  -- closes within 3 ns, and the remaining 4.25 kg/m³ persists at
  7.2 sigma, with a last-half trend of +0.030 kg/m³/ns, i.e. not decaying.
  (State the baseline and the window whenever this number is quoted: the
  same data give 83.5-87.8% across the defensible choices of each, so a
  bare percentage is not a fact.)  So heat and time buy most of the
  relaxation quickly and then buy nothing, and a residual difference
  survives that a longer anneal at this scale does not appear able to
  remove.  Whatever replaces the current anneal should
  therefore be argued as a large improvement, not as a fix; an entry
  promising equilibration would oversell what was measured.

  Note what that hold also says about *T*:sub:`g` itself: **relaxation does
  not switch off at the glass transition**, it slows continuously through
  it.  A transition 20-34 K wide by experiment means 11 K below
  *T*:sub:`g` is inside the transition, not deep in the glass, so fast
  early relaxation there is expected rather than anomalous.  This does not
  soften the argument above -- the 197x figure is unchanged and a peak well
  clear of the transition is still the lever -- but it does mean *T*:sub:`g`
  is a soft kinetic boundary for this material, which any protocol
  expressed relative to it has to accommodate.

  **Untested.** Nobody has yet built with a hotter anneal and compared it
  against the melt-and-recool value, which is the experiment that would
  settle it; the argument above is mechanism, not measurement. That is why
  the example is documented rather than changed. Raising the peak also costs
  wall-clock in a tutorial that should stay small, so the test needs to
  report both whether it works and what it costs.

  The general problem is worth more than the instance: **any protocol whose
  relaxation depends on being above *T*:sub:`g` is fragile when *T*:sub:`g`
  is unknown at configuration time**, which it always is -- you cannot know
  it before building the thing. A fixed absolute peak temperature is a guess
  that happens to be right or wrong per chemistry, silently. Something
  expressed relative to an estimated or measured *T*:sub:`g`, or a postcure
  that measures *T*:sub:`g` and then anneals accordingly, would be robust
  where an absolute number is not. That is a real design change, not a
  config edit.

  Raised by the calibration study, which had already published the
  pessimistic density before catching it -- the failure mode is reporting a
  protocol artifact as a force-field result.

- **Example 5 is still not laptop-scale.** The `68d81fb` retune cut the
  monomer pool to make HTPB/IPDI "fit in a reasonable wall time on a
  laptop", but a validation run on 16 CPU cores took **10h17m** — about
  ten times any other example (the next longest, PACM/DGEBA, was 2h15m).
  The cost is structural: densification starts at 10 kg/m³ in a 27 nm box
  and needs all 50 NPT repeats to compress geometrically (a steady 2.9%
  box-side reduction per repeat) up to ~860 kg/m³, and the cure stage then
  does drag/relax on long chains. Worth revisiting whether the low
  starting density is really necessary, or whether the example can start
  denser with a shorter densification.


- **More cyanate-ester variants.** Example 6 covers bisphenol-A
  dicyanate. The same reaction and repair machinery carries over to other
  bisphenol bridges with only a SMILES swap in `constituents` — bisphenol
  F (methylene), bisphenol E (methylethylidene), hexafluorobisphenol A,
  thioether, sulfone, dicyclopentadiene. A homologous series would
  exercise the repair stage across chemistries and give the tutorials a
  structure-property story.

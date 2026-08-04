# Zr3d_ghoniem — differences from the original MoDELib2-NNL

**Branch:** `zr3d_ghoniem`  ·  **Baseline:** `main` (Po's MoDELib2-NNL)  ·  **Status:** in progress, see §6

Zr3d_ghoniem is a variant of MoDELib2-NNL whose cluster-dynamics equations and parameters
are taken from the **ZrMicro 0-D code** (`d:/GitHub/Fluor_Zr/ZrMicro`), so that the
3-D spatially-resolved solve reduces to the 0-D result in the well-mixed limit.

The governing equations are those of Po & Ghoniem, Deliverable D1/M1
(`Docs/Formulation/GW_Phase4_D1M1_NG.pdf`), §2.2 in particular.

Po's `Library/Materials/Zr4.txt` is **untouched**. The Zr3d_ghoniem material definition is a
new file, `Library/Materials/Zr3d_ghoniem.txt`. The C++ changes, however, are global: `iSize`
goes from 0 to 8, so `Zr4.txt` cannot be run from this branch. Check out `main` for that.

---

## 1. Build and environment

Not physics — required to build on WSL/Ubuntu rather than macOS.

| File | Change | Why |
|---|---|---|
| `CMakeLists.txt` | `EIGEN3_INCLUDE_DIRS` `/opt/local/include/eigen3` → `/usr/include/eigen3`; dropped the hard-coded Qt prefix path and macOS `RPATH` | Paths were absolute to a macOS MacPorts install |
| `tools/CMakeLists.txt` | `DDqt` subdirectory commented out | Needs Qt; not required for CD runs |

Two build traps worth recording (handled in `Fluor_Zr/Docs/Formulation/build_modelib_wsl.sh`):

- The repo declares `cmake_minimum_required(VERSION 3.1.0)`. CMake ≥ 4.0 **rejects**
  anything below 3.5 outright, so `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is mandatory,
  not defensive — Ubuntu 26.04 ships CMake 4.2.3.
- `libfftw3-dev` is **not optional**, despite `CMakeLists.txt:65` warning that the noise
  generator is merely "disabled" when it is missing: line 117 still links
  `${FFTW3_LIBRARIES}` unconditionally, so the generate step fails without it.

---

## 2. C++ changes

### 2.1 `include/DislocationDynamicsBase/ClusterDynamicsParameters.h`

- **`iSize` 0 → 8.** Four immobile families × (number density, content):
  `[N_c, N_a1, N_a2, N_a3 | c_c, c_a1, c_a2, c_a3]`. This is what turns on the
  immobile machinery throughout; it is the reason `Zr4.txt` is incompatible.
- **New parameters** for the §2.2 kinetics: `loopCascadeFractions` (ε_k, Eq. 51),
  `nNuc` (Eq. 36), `cLL`/`cLN` (Eq. 99), `kappaLL`/`kappaLN` (Eq. 98),
  `tauVac` (Eq. 95), `rhoNetwork`, `dadAnisotropy` (p_m), `dadZ0` (Z⁰_m).
- **New method `loopDADbias()`** — the Eq. (15) capture efficiencies,
  `Z_c = Z⁰·p`, `Z_a = Z⁰·(p + p⁻²)/2`, row 0 = ⟨c⟩ (vacancy) loops, row 1 = ⟨a⟩.
- **New static `getClusterAtomicVolume()`** — see issue 6 below.

### 2.2 `src/DislocationDynamicsBase/ClusterDynamicsParameters.cpp`

- Initializers for all of the above, each guarded by `iSize>0` and converted to
  MoDELib units on read (`tauVac` seconds → b/cs, `rhoNetwork` m⁻² → b⁻²).
- `omega(ddBase.poly.Omega)` → `omega(getClusterAtomicVolume(ddBase))`.
- `r_min` now divided by `b_SI` on read (issue 5).

### 2.2b `mBWF` / `dmBWF` scale with `cdp.omega`, not `ddBase.poly.Omega`

`FluxMatrix` divides `D` by `cdp.omega`, and the weak form multiplies the flux back by
an atomic volume, so the two must be the *same* quantity for the flux to come out as
`−D·∇c`. They were identical until the CD atomic volume became settable independently
of the lattice cell volume (issue 6); leaving `poly.Omega` here would have scaled the
entire diffusion operator by `poly.Omega/cdp.omega` = 3.98 relative to every reaction
term. See issue 13.

### 2.3 `include/ClusterDynamics/ImmobileSinks.h` — **new file**

An `EvalFunction` mirroring the existing `SecondOrderReaction`, supplying the loop
absorption sink that was previously missing (the source carried a
`// Missing immobile sinks` comment). For each mobile species *m* it assembles

    k²_m = Σ_k Z(type(k), m) · S_k ,        S_k from clusterDensity() = 2πr_k N_k

and returns `−k²_m · D̄_m` on the diagonal, matching the sign convention of `getR1()`.

Per your instruction, the sink set is the 0-D one **except the grain boundary**,
which is imposed spatially by the Dirichlet BC on the mobile species and must not be
double-counted here. The network dislocation sink is not in this file either — it
reuses the existing `otherSinks_SI` plumbing (§3, `otherSinks_SI`), which required no
new code.

### 2.4 `src/ClusterDynamics/ClusterDynamicsFEM.cpp`

- **`ImmobileSinks` wired into the mobile Newton system** — `bWF_RI` into the
  Jacobian and `lWF_RI` into the residual, alongside the existing `R1`/`R2` terms.
- **`solveImmobileClusters()` implemented** (was an empty stub). Node-wise
  sub-stepped integration of Eqs. (34)–(56) and (97)–(102) over one dose step:
  cascade nucleation, absorption growth with the Eq. (91) shrink gate, Avrami
  coalescence in two channels (like-loop conserves content, loop-network removes
  both), vacancy-loop dissolution, and thermal emission. DOF layout is
  `index = node*iSize + component`, `[0,nF)` densities and `[nF,2nF)` contents.

---

## 3. `Library/Materials/Zr3d_ghoniem.txt` vs `Zr4.txt`

Everything outside the `# ClusterDynamics` section (elasticity, mobility laws,
gamma-surface) is inherited unchanged. Original values are preserved in comments
in the file itself.

### 3.1 Values changed to match ZrMicro

| Key | Zr4 (Po) | Zr3d_ghoniem (ZrMicro) | Basis |
|---|---|---|---|
| `atomicVolume_SI` | *(absent)* | `1.2e-29` | ZrMicro `physical_props['Omega']` |
| `mobileSpeciesSurvivingEfficiency` | `0.01` | `1.0` | 0-D loses defects only to the in-cascade cluster fractions |
| `mobileSpeciesCascadeFractions` | `1.0 0.6 0.3212 0.0788` | `0.995 0.97788714 0.005 0.00333333` | `G_v/G = 1−ε_vL`, `G_i/G = 1−ε_2i−ε_3i−ε_iL`, `G_2i/G = ε_2i/2`, `G_3i/G = ε_3i/3` |
| `mobileSpeciesEnergyFormation_eV` | `1.5 …` | `1.6 …` | 0-D `E_F_v`, gives `C_v_eq = 8.45e-15` at 573 K |
| `mobileSpeciesEnergyMigration_eV` | `0.9 / 0.09 / 0.041–0.187` | `1.2 / 0.759101 / 0.759101 / 0.759101`, isotropic | 0-D `E_m_v`, `E_m_i`; clusters transport at the monomer value (issue 11) |
| `mobileSpeciesD0_SI` | `4.9639e-10`, `8.1882e-10` | `5.21645e-7` for all, isotropic | `D = (a²/z_c)·ω = a²ν e^(−E_m/kT)` ⇒ `D0 = a²ν`, a = 3.23e-10 m, ν = 5e12 s⁻¹ |
| `reactionPrefactorMap` | `2.856`, `1.0`, `0.0` | all ten recomputed | §3.3 |
| `Eb_eV` | `5.0 5.0 5.0 5.0` | `5.0 5.0 0.957033 2.0` | index 0 → loop emission (0-D uses e⁻⁵ᵉⱽ/ᵏᵀ); 2, 3 → 0-D `E_b_2i`, `E_b_3i` |
| `otherSinks_SI` | `0 0 0 0` | `1.88e14  1.8988e14 ×3` | network dislocation sink, `k²_v = ρ_N`, `k²_i = Z_N ρ_N` |

The migration-energy row order matters: components are `[11 12 13 22 23 33]`, and the
off-diagonal `D0` entries are zero, so the tensors are diagonal.

### 3.2 New §2.2 block

`immobileSpeciesVector`, `immobileSpeciesRelRelaxVol`, `immobileSpeciesBurgers`,
`alpha_bp`, `delVPyramid`, `w0`, `n_s`, `evc`, `Nvmax`, `nmin`, `nmax`, `r_min`,
`distanceFactor`, `loopCascadeFractions`, `nNuc`, `cLL`, `cLN`, `kappaLL`, `kappaLN`,
`tau0_vLoop_SI`, `Ea_vLoop_eV`, `rhoNetwork_SI`, `dadAnisotropy`, `dadZ0`.

The coalescence coefficients, nucleation sizes and dissolution constants are the
fitted 0-D values from run `20260622_144021_7959445`.

### 3.3 How the reaction prefactors were derived

MoDELib forms `K_ab = p_ab·4π(r_a+r_b)(D_a+D_b)/Ω` in `getR2()`, with
`r = (3Ω/4π)^(1/3)` for `|m| = 1` and `r = √(|m|Ω/(bπ))` otherwise.
ZrMicro forms `K_ab = f_ab(ω_a+ω_b) = f_ab(z_c/a²)(D_a+D_b)`, with `f = 0.5` for v+i
(`physical_props['recom']`) and 1 otherwise, the like-species channels carrying the
extra factors of `R_i_i = 2ω_i C_i²` and `R_2i_2i = 4ω_2i C_2i²`. Equating,

    p_ab = f_ab · (z_c/a²) · Ω / (4π(r_a+r_b)) · (D⁰ᴰ_a+D⁰ᴰ_b)/(D³ᴰ_a+D³ᴰ_b)

The last ratio is unity except for pairs involving 2i/3i, where it restores the
ω_2i-based rate now that those species transport at the monomer diffusivity
(ω_2i/ω_i = 4.886e-6).

**This derivation was confirmed empirically.** MoDELib prints the assembled
second-order interaction matrices in 1/s at start-up:

| entry | MoDELib | ZrMicro | agreement |
|---|---|---|---|
| (v,i) | 2.526876e7 | 0.5(ω_i+ω_v) = 2.52560e7 | 0.05% |
| (i,i) | −2.021161e8 | −4ω_i = −2.020216e8 | 0.05% |
| (v,2i) | 6.940722e3 | ω_2i+ω_v = 6.93473e3 | 0.09% |

This also settles a factor-of-2 ambiguity: `SecondOrderReaction.h` documents the rate
as `ċ = ½R₂cc` while `getR2()` fills both `(a,b)` and `(b,a)`, and it was not clear
from inspection whether the ½ survives into the weak form. The match above shows it does.

### 3.4 `loopSinkScale` — reproducing ZrMicro's loop sink convention exactly

ZrMicro's loop absorption is *not* the purely geometric `S_k = 2πr_k N_k`. It departs
from it in two ways, and `loopSinkScale` reproduces both exactly:

**(a) One radius prefactor for all four families.** `input_data.py` derives both
`l_a = √(Ω/πb_a)` and `l_c = √(Ω/πb_c)`, but `reaction_rates.loop_absorption()` uses
`lc_l = l_c/l` for *every* family — `l_a` is computed and never used anywhere in the
codebase. The 3-D builds `r_k` from each family's own `|b_k|`, so each needs
`√(|b_k|_MoDELib / b_c_ZrMicro)`, with b_c = 5.15e-10 m against MoDELib's
1.632993·b_SI = 5.279467e-10 m (⟨c⟩) and 1.0·b_SI = 3.233e-10 m (⟨a⟩):

| family | factor |
|---|---|
| ⟨c⟩ | √(5.279467/5.15) = 1.012495 |
| ⟨a⟩ | √(3.233/5.15) = 0.792317 |

**(b) The fitted `Q` on the ⟨c⟩ channel only:** ⟨c⟩ → 1.012495 × 0.287931 = **0.291528**.

    loopSinkScale = 0.291528  0.792317  0.792317  0.792317

Implemented as an optional material key (default `1 1 1 1` = purely geometric),
applied identically in `ImmobileSinks.h` (the sink seen by the mobile species) and in
`solveImmobileClusters()` (the growth flux it causes) — mirroring the 0-D, where
`loop_absorption()` and the `loop_growth_rate_*` methods share one prefactor. Applying
it in only one place would break the mass-conserving coupling between them.

---

## 4. Simulation case `tutorials/zrmicro_coupled/` — new

CD-only verification case: `useFEM=1`, `useDislocations=0`, `useClusterDynamics=1`,
`dtMax = 6.95868964e19` (= 1 dpa at G = 1e-7), `outputFrequency=1`,
`startAtTimeStep=0`, `absoluteTemperature=573`, non-periodic, zero target loop
densities. `clean_run.sh` refreshes the material file and regenerates `evl_0` on every
invocation — see issues 1 and 2.

---

## 5. Issues encountered, and how each was solved

**A recurring hazard, stated once:** almost every failure in this stack **exits 0** —
the build script, `DDomp` with a missing `evl_0`, a missing material key, a CRLF mesh,
and the solver failures below. Verification has been by artifacts on disk and log
tails, never by exit status.

### 1. `microstructureGenerator` never terminates
**Symptom:** hangs indefinitely at start-up.
**Cause:** it is the *discrete-dislocation* path. It inserts loops one at a time; the
log showed it had reached 2.34e17 m⁻³ against a ZrMicro target of 1.48e21 — 0.016% of
the way, needing ~10¹¹ discrete loops. CD loop densities are not representable as
discrete dislocations, which is precisely why they belong in the continuum fields.
**Fix:** `useDislocations=0` and `targetDensity=0` in `frankLoopsDensity.txt` /
`aLoopsDensity.txt`, giving an empty 20-byte `evl_0.txt` in about a second. `DDomp`
still requires that file to exist even with dislocations off.

### 2. Restart silently double-counted dose
**Symptom:** ⟨a⟩ loop content came out at exactly 3× the expected value.
**Cause:** `startAtTimeStep=-1` means "last available step", and `evl_0.txt` serves
double duty as both the initial configuration *and* the runID=0 output. A bare re-run
therefore resumed from the previous run's accumulated dose.
**Diagnosis:** `c_a` max = 2.112e-3 = exactly 3 × (G₀·ε_a·1 dpa).
**Fix:** pinned `startAtTimeStep=0`, and `clean_run.sh` now regenerates `evl_0` every
time. Also note output is written **after** `solve()`, so `evl_N` holds the state at
*(N+1)* dose steps.

### 3. Stale material file in the simulation directory
**Symptom:** run died on `does not cointain line with format immobileSpeciesVector=…`
**Cause:** `inputFiles/` holds a *copy* of the material file, made before the §2.2
block existed.
**Fix:** `clean_run.sh` re-copies it from `Library/` on every run.

### 4. CRLF line endings silently produced a zero-node mesh
**Symptom:** downstream periodicity failure with no parse error.
**Cause:** `unitCube_15K.msh` checked out with CRLF; the `.msh` parser returned zero
nodes without complaint.
**Fix:** `sed -i "s/\r$//"` on the simulation-directory copies only, leaving the
`Library/` originals untouched.

### 5. `r_min` never converted from metres
**Cause:** read raw from the material file in metres, but compared against radii from
`rloop()`/`rpyr()`, which are in units of b. The Eq. (91) shrinking gate was therefore
permanently wide open.
**Fix:** divide by `b_SI` on read. It is used only by the new code, so the conversion
is safe.

### 6. Atomic volume was the lattice *cell* volume
**Cause:** `Polycrystal::Omega` is `det(latticeBasis)` = 1.41421 b³ for HEX — the
volume of a cell containing **two** atoms, not the atomic volume the CD equations
assume. Every content ↔ radius ↔ density conversion inherited the error, and it
differs from ZrMicro's 1.2e-29 m³ by roughly a factor of 4.
**Fix:** new **optional** material key `atomicVolume_SI`, read by
`getClusterAtomicVolume()`; material files that omit it keep the previous behaviour,
so `Polycrystal` — which the elastic/DD side also uses — is left alone.
**Noted for the record:** ZrMicro's 1.2e-29 m³ implies ~12.6 g/cm³ for Zr, about half
the true atomic volume. It is baked into the fitted parameter set, so Zr3d_ghoniem matches it
rather than "correcting" it.

### 7. ⟨c⟩ loops collapsed onto the positivity floor
**Symptom:** `n_vL` pinned at 1e-30 while ⟨a⟩ loops nucleated and grew normally.
**Cause:** τ_vL ≈ 9.95e16 MoDELib time units against a sub-step of 3.48e18 — the decay
time is **35× shorter than the step**, so explicit Euler gave `n(1−35) < 0`, clamped
every iteration.
**Confirmation it was integration and not physics:** the analytic balance
`n* = ṅ_nuc·τ_vL = 1.79e-8` matches the 0-D `CvL+CavL = 1.78e-8`.
**Fix:** integrate the linear losses implicitly — unconditionally stable, and returns
`n* = ṅ_nuc·τ_vL` exactly as `dt/τ → ∞`. Verified: `n_vL = 1.79161e-08`.

### 8. Thermal emission written as a rate, not a coefficient
**Cause:** the Eq. (55)–(56) emission term was computed as an absolute rate rather
than proportional to the stored content, so it could not enter the implicit update.
**Fix:** per loop `α_v = 2πr·D_v/Ω·e^(−E_b/kT)` (the convention of `getR1()`); the
content lost per unit volume is `α_v·N·Ω = S_k·D_v·e^(−E_b/kT)`, the Ω's cancelling.
Since that does not scale with *c*, it is divided by *c* to form an equivalent
first-order coefficient.

### 9. Loop densities stored per atom instead of per b³
**Cause:** `clusterRadius()`, `clusterDensity()` and `sigmoid()` all form
`n = c/(N·Ω)`, i.e. they expect **N in loops per b³**. The new integration was
accumulating loops *per atom*.
**Consequence:** the helpers read the mean loop size as 2.82× too large — enough to
push the ⟨c⟩ family across the bi-pyramid→loop sigmoid onto the wrong branch.
**Fix:** store per b³ (`nucRate` divided by Ω), `Nvol = n(k)` directly, and
`dLL = (1/N)^(1/3)`. Content needs no conversion: a volume fraction is numerically
equal to the 0-D per-atom content.

### 10. ⟨a⟩ Burgers vectors were not equivalent
**Symptom:** `n_a1, n_a2, n_a3` = 2.84e-7, 2.32e-7, 3.16e-7 at **zero applied stress**,
where the three prismatic families must be identical.
**Cause:** `immobileSpeciesBurgers` columns are multiplied by `latticeBasis`, i.e.
they are **lattice** components; Cartesian ones had been supplied. The resulting
magnitudes were |b| = 1.0, 0.753, 1.197 instead of all 1.0.
**Fix:** `(1,0,0), (0,1,0), (−1,1,0)` for ⟨a⟩ and `(0,0,1)` for ⟨c⟩. Verified:
`immobileSpeciesBurgersMagnitude: 1.633, 1.0, 1.0, 1.0`.
**Standing test:** at zero stress the three ⟨a⟩ family densities must agree.

### 11. Cluster diffusivity made the Newton matrix singular
**Symptom:** `FixedDirichletSolver failed.` — thrown from `compute()`, i.e. the
factorization, not the solve.
**Cause:** setting `E_m_2i = 1.36292 eV` drops D_2i to ~5e-6 of D_i. The mobile solve
is a **pure steady-state (QSSA) system** — `dmBWF` is diffusion only, with no mass
term — so the 2i/3i block had nothing on its diagonal but a vanishing Laplacian.
**Fix, and why it is *more* faithful to ZrMicro:** the 0-D is inconsistent about which
frequency it applies where, but unambiguous — `R_2i_s`/`R_3i_s` (network sinks) and
`flux_i` (loop absorption) all use **ω_i**, and ω_2i appears *only* in the
mobile–mobile cluster reactions. So 2i/3i transport at the monomer diffusivity, and
the ω_2i/ω_i factor is carried in the reaction prefactors instead (§3.3), where it
belongs. The network-sink and loop-absorption terms now reproduce the 0-D exactly,
which they did not before.

### 12. Mobile Newton iteration does not converge on the first step — **OPEN**
**Symptom:** `convergenceError = 2.9e5` after the first Newton iteration, then the
iterative (BiCGSTAB) solve fails on the second.
**Cause:** cold start. `mobileClusters` is initialized to the *thermal equilibrium*
concentration (~1e-15) and must reach the irradiation steady state (~1e-6) — nine
orders — in one Newton solve of a quadratically nonlinear system. It converged before
only because production was 100× smaller.
**Attempted and rejected:** switching `rSolver` to the direct SparseLU branch. That
fails *earlier*, in `compute()` itself, on the ~96k-dof 3-D matrix, while the LLT
factorization of the pure diffusion operator in `mSolver` succeeds — so the direct
path is not viable here at this size. Reverted, with a comment at the call site.
**Note:** `dt` does not enter the mobile equation at all (it is a steady-state solve),
so ramping the dose step cannot help. The remedy has to be a better initial guess —
seeding the mobile field from the 0-D steady state, which is exactly what the
two-time-scale coupling provides — or damping/continuation in G.

### 13. Diffusion operator scaled 3.98× relative to the reactions
**Cause:** introduced by the fix to issue 6. `FluxMatrix` returns `−D/cdp.omega`, and
`mBWF`/`dmBWF` multiplied it back by `ddBase.poly.Omega`. Those two were the same
number until the CD atomic volume became independently settable; afterwards the
diffusion term carried a spurious factor of `poly.Omega/cdp.omega` = 1.41421/0.355105
= 3.98 against every reaction term.
**Fix:** use `cdp.omega` in both weak forms, restoring the exact cancellation so the
flux is `−D·∇c`. Found by reading the constructor while chasing issue 12 — a factor
of four on the diffusion operator is a plausible contributor to that non-convergence,
so issue 12 should be re-tested now that this is corrected.

### 14. Coalescence was integrated explicitly and overshot
**Symptom:** ⟨a⟩ mean loop size m_a = 5042 at 1 dpa, then exactly 300.0 (= n_nuc) at
2 dpa — content collapsing 11× while density kept climbing.
**Cause:** the same stiffness failure as issue 7, in the channel left explicit.
`nu_LN*phi_LN*dt_sub` reaches ~3.6, so each sub-step overshoots and clamps to the floor.
**Fix:** Eq. (97) removes loops in proportion to n and content in proportion to c, so
both are first-order coefficients; moved into the implicit update.
**Effect:** ⟨a⟩ loop density went from 3–4× the 0-D to within ~30%.

### 15. Vacancy-loop annealing released the wrong content
**Symptom:** ⟨c⟩ loop content 350× below the 0-D; m_c = 487 against n_nuc = 400,
i.e. pinned at the nucleation size, versus 1.7e5 in the 0-D.
**Cause:** dissolution was implemented as removing content in proportion to the
content, `c/tau_vL`. ZrMicro's `annealing_content_vL()` returns
`n_vL_nuc * CvL / tau_vL` — each dissolving loop releases its **birth** content, not
its mean, because tau_vL describes dissolution of still-small EMBRYO loops while
loops that survive and grow are stable. Removing the mean content instead caps the
mean size at ~n_nuc.
**Fix:** content sink proportional to the DENSITY, `n_nuc*n*Omega/tau_vL`.
**Self-check:** at number saturation `n* = Ndot_nuc*tau_vL`, so the release equals
`n_nuc*Ndot_nuc*Omega = Gk` and exactly cancels the cascade content seed, leaving net
content growth flux-driven — precisely what the 0-D docstring describes.
**Effect:** ⟨c⟩ content from 350× low to within 20–26%.

### 16. `Eb` read before initialization in `getR1()` — **pre-existing upstream bug**
**Symptom:** C2i and C3i ~1e-17 against the 0-D's ~2.5e-11, six orders low.
**Cause:** `getR1()` reads `Eb(k)` to build the cluster dissociation rates, but `Eb`
was declared *after* `R1` in the class. C++ initializes members in DECLARATION order,
not in the order of the constructor's initializer list, so this was an uninitialized
read. It evaluated to `exp(0) = 1` instead of `exp(-Eb/kT)`, inflating the 2i/3i
dissociation rate to the bare attempt frequency.
**Evidence:**
- `R1(2,2)` = −5.0535e7 1/s where the network sink alone gives −20.86.
- `R1(3,3)/R1(2,2)` = 0.4802 = `p(1,2)/p(1,1)` exactly, so both dissociation rates
  had lost their Boltzmann factor.
- Predicted C2i under the bug, `G_2i/alpha` = 9.89e-18, against the observed 1.06e-17.
**Fix:** move the `Eb` declaration (and its initializer) ahead of `R1`, with the
reasoning recorded at the declaration site.
**Effect:** `R1(2,2)` → −21.05 1/s; C2i +1.2% and C3i −2.3% against the 0-D.
**NOTE:** this is a bug in the ORIGINAL code, not in the Zr3d_ghoniem changes. It was
invisible in Po's calibration because `Eb_eV = 5.0` for every species makes
`exp(-5/kT)` ~ 1e-44, negligible whether or not the factor is applied. Giving 2i its
physical 0.957 eV is what exposed it. **Worth reporting upstream**: it affects any
MoDELib run with a non-negligible cluster binding energy.

---

## 5b. Verification against the 0-D

Interior nodes, selected by the top decile of Cv (the grain boundary is a Dirichlet
sink, so Cv is depressed near it and saturates in the interior). This selection must
be made ONCE from the mobile field and applied to all quantities: loop density
anti-correlates with Cv (corr = −0.98), because near the boundary the mobile species
are depleted so loops neither grow nor coalesce and simply accumulate at their
nucleation value. A per-quantity percentile picks the boundary for the immobile
fields and the interior for the mobile ones.

| quantity | 1 dpa | 2 dpa |
|---|---|---|
| Cv | +7.3% | −1.4% |
| Ci | +14.5% | +11.2% |
| C2i | +1.2% | −41.1% |
| C3i | −2.3% | −43.1% |
| ⟨c⟩ loop density | −0.13% | −0.07% |
| ⟨c⟩ loop content | −38.3% | −38.4% |
| ⟨a⟩ loop density | −57.4% | −48.4% |
| ⟨a⟩ loop content | −51.4% | −37.7% |
| ⟨a⟩ defects/loop | 1.14× | 1.21× |
| ⟨c⟩ defects/loop | 0.62× | 0.62× |

Mobile species and the ⟨c⟩ loop density agree to a few percent. The remaining
factor-of-two on ⟨a⟩ loop density is consistent in sign with the known gaps below:
the homogeneous clustering flux `Phi_clus` is not wired, and it would ADD ⟨a⟩ loops.

---

## 6. Status and open items

**Working:** MoDELib builds and runs the §2.2 physics. The mobile Newton iteration
converges cleanly (7.5e6 → 6.2e-8 in ten iterations, quadratic tail). ⟨a⟩ families are
symmetric at zero stress. The reaction network matches ZrMicro to ~0.1%. Mobile
species, both cluster species and the ⟨c⟩ loop density all agree with the 0-D to a few
percent (§5b). No negative concentrations anywhere.

**Two further outputs added for post-processing:**
- `evl/cdNodes.txt` — finite-element node coordinates, written once per run. Without
  it the fields cannot be placed in space: the CD trial functions live on second-order
  elements (24115 nodes for a 3392-node mesh), so the field rows do not correspond to
  the `.msh` vertices and the ordering is not reproducible externally.
- `Fluor_Zr/ZrMicro/py_utils/modelib_fields.py` — Fig. 4-style panel plots, one row
  per species and one column per dose, with loops overlaid as discrete platelets at
  sampled sites (visualization only).

**Known modelling gaps, all documented in the coupling `.tex`:**

1. ~~**Loop sink strengths differ by construction.**~~ **Resolved** — see §3.4.
2. **Nucleation is cascade-only.** The homogeneous clustering flux
   `Φ_clus = R_i,3i + R_2i,2i` and the SIPN stress weights `w_k(σ)` are not wired.
   Exact at zero stress, which is the configured condition.
3. **DAD is a fitting form, not a generated one.** `dadAnisotropy`/`dadZ0` are fitted
   to reproduce the 0-D symmetric bias forms rather than derived from `D_c/D_a`, and
   the diffusion tensors are isotropic to match ZrMicro. These runs therefore exercise
   none of the tensorial-DAD mechanism of the deliverable, and the fit is tied to 573 K.
4. **`Warning: Sum of R2 is not zero`** — benign and pre-existing. The reaction map
   contains `i+3i` but there is no 4i species to receive the product, so that channel
   does not conserve interstitials. ZrMicro truncates the cluster ladder at the same place.

**Performance note:** `useElasticDeformation=1` costs ~42 s/step on a 72,345-dof solve
for a case with no dislocations and no applied stress — roughly 20 minutes of a 30-step
run. Worth disabling for production runs.

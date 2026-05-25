# Reader's guide to `conductivity_paper.py`

The file itself is **read-only by convention** (sibling dependency from a published paper — do not modify). This document explains what each function does, what equation it implements, and why.

If you're learning this code for the first time, read this in order:
1. The physical picture (below)
2. `_shedlovsky` — what it computes
3. `variant_shedlovsky` — the public single-salt entry point
4. `msa_transport` — the multi-salt path (binary/ternary mixtures)

---

## 1. The physical picture

We measure **conductivity** (μS/cm) in the lab. We need **concentration** (mM) for the diafiltration model. These two quantities are related but not identical:

```
   σ (conductivity)   =   c (concentration)  ×  Λ (equivalent conductivity)
                                                 ^^^^^^^^^^^^^^^^^^^^^^^^^
                                                 NOT constant — depends on c
```

In a dilute solution, Λ is roughly constant ("limiting equivalent conductivity" Λ₀). But as concentration rises, ions start interacting electrostatically and Λ DROPS. The Shedlovsky model gives us a formula for how Λ depends on c, so we can do the inversion `σ → c` correctly.

There are two regimes the codebase supports:

| Regime | Function to use | Why |
|---|---|---|
| Single salt (NaCl, CaCl₂, LaCl₃, ...) | `variant_shedlovsky` | Closed-form Shedlovsky expression, fast |
| Binary or ternary salt mixtures | `msa_transport` | MSA = Mean Spherical Approximation; handles ion-ion coupling that Shedlovsky can't |

---

## 2. `_shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion)`

### What it computes

The **equivalent conductivity** Λ as a function of salt concentration `c`:

```
Λ(c)  =  Λ₀  −  (B₁·Λ₀ + B₂) · √I  /  (1 + a·B·√I)
         ^^      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         what    correction term — gets BIGGER as c (and √I) rise.
         Λ would Captures how ion-ion interaction slows charge transport.
         be at
         c → 0
```

### What the constants mean

| Symbol | Meaning | Units |
|---|---|---|
| Λ₀ | limiting equivalent conductivity (at c → 0) | cm²·S/equiv |
| I | ionic strength = ½ Σ cᵢ zᵢ² | M |
| B | Debye–Hückel parameter for activity | 1/cm |
| B₁ | relaxation correction coefficient | dimensionless |
| B₂ | hydrodynamic (electrophoretic) correction | cm²·S/equiv |
| a | ion size parameter (distance of closest approach) | cm |
| ε | dielectric constant of solvent (water) | dimensionless |
| η | viscosity of solvent | poise |
| T | temperature | K |

### The temperature-dependent prefactors

```
B   =  50.29  ×  10⁸  ×  (ε·T)⁻¹/²
B₁  =  2.801  ×  10⁶  · |z₁·z₂| · q  /  [ (ε·T)^(3/2) · (1 + √q) ]
B₂  =  41.25  · (|z₁| + |z₂|)  /  [ η · √(ε·T) ]
```

where `q` is a dimensionless function of valencies and limiting conductivities of the cation and anion.

**The connection to the EC25 compensation:** Notice that B, B₁, B₂ all depend on T through the (ε·T) factor. So if we feed `temp` as the actual lab temperature (e.g., 294 K = 21 °C), Shedlovsky gives us Λ at 21 °C — which the inversion `σ → c` would then over-correct compared to a 25 °C-anchored table. That's why our loader does the EC25 compensation upstream and then forces `temp = 298.15 K` (25 °C) before calling Shedlovsky.

### The ionic strength loop

```python
I[i] = 0.5 * (cation_conc[i] * z₁² + Cl_conc[i] * z₂²)
```

For a single salt MₓCl₍ᵧ₎ that fully dissociates:
- cation_conc = c (whatever you passed in)
- Cl_conc = |z₁| · c (charge balance: |z₁| anions per cation)

So I scales linearly with c. The √I in the correction term is what makes the model **non-linear** in concentration.

### Return value

A list of Λ values, same length as `conc`, each in **cm²·S/equiv**.

---

## 3. `variant_shedlovsky(...)` — the single-salt public entry point

### What it computes

The **specific conductivity** κ (what your meter actually reads) in **mS/cm**:

```
κ  =  c_equiv · Λ  ×  conversion-to-mS-per-cm
```

where `c_equiv = c · |z₁|` is the equivalent concentration (one mole of charge per equivalent).

### How it differs from `_shedlovsky`

`_shedlovsky` returns Λ in cm²·S/equiv. `variant_shedlovsky` calls `_shedlovsky` and then converts Λ → κ by multiplying by `c_equiv / 1000` (M → eq/cm³) and then by 1000 (S → mS).

### Why the public function is called "variant"

It's the variant that returns **specific conductivity** (κ) rather than **equivalent conductivity** (Λ). The loader uses this version because κ is what the meter reads.

### Inverse use in our loader

The loader actually does the OPPOSITE of `variant_shedlovsky`: it has κ (measured) and wants c (unknown). It does a **bisection inversion**: pick a c, compute κ_predicted = variant_shedlovsky(c, ...), compare to κ_measured, adjust c. Repeats until they match.

```
   measured κ  ──►  bisect over c  ──►  c that produces matching κ
```

This bisection is implemented as `_invert_shedlovsky` (or equivalent) in `refactored_ucb_library.py`.

---

## 4. `msa_transport(...)` — the multi-salt path

### When you need this

The single-salt Shedlovsky model assumes one cation + one anion. For mixtures (e.g., NaCl + CaCl₂, or NaCl + CaCl₂ + LaCl₃), ions of different valencies interact in ways Shedlovsky can't capture. The Mean Spherical Approximation (MSA) is a thermodynamic model that handles this.

### The setup

For an N-species system (cations + one common anion, usually Cl⁻):

| Input | Meaning |
|---|---|
| `valency` | list of integer charges, e.g. `[1, 2, -1]` for Na⁺, Ca²⁺, Cl⁻ |
| `diameters` | hard-sphere diameters of each ion (m), used in pair-correlation terms |
| `diff_coeff` | infinite-dilution diffusion coefficients (m²/s) |
| `lambda_0` | limiting molar conductivities (S·m²/mol) |
| `salt_*_conc` | per-time-step concentrations of each salt (mM) |

### The math, in plain words

For each timestep, the MSA computes:

1. **Number densities** of each ion species (concentration × Avogadro).
2. **Debye length** κ (screening length set by total ionic strength).
3. **Bare mobility** ω of each ion from its diffusion coefficient: ω = D / (k_B · T).
4. **Mean mobility** ω̄ as a weighted average across species.
5. **Transport numbers** t_j = fraction of current carried by species j.
6. **Hydrodynamic correction** (Δv/v): how much the bulk fluid drag slows each ion at finite concentration.
7. **Relaxation correction** (Δκ/κ): how much the ionic atmosphere lags behind a moving ion.
8. Then it assembles the **per-species conductivity**:
   ```
   κ_species = (e² · n · D · z²) · (1 + Δv/v) · (1 + Δκ/κ) / (k_B · T)
   ```
9. Sum across species → total bulk conductivity in S/m, then convert to mS/cm.

### The bisection inside MSA

Step 4 above requires solving an implicit equation:

```
α · Σⱼ  t_j / [ (ω_j/ω̄)² − (α/ω̄)² ]  =  0
```

for the unknown α (a kind of effective relaxation rate). `bisect()` from `scipy.optimize` does this — there's one bisection root call per ion species per timestep.

### Why this matters for our DATA3 single-salt sheets

The single-salt sheets (which we're focusing on first) use `variant_shedlovsky`, not `msa_transport`. The MSA path is dormant for those sheets. It's there for the multi-salt sheets in the wider 26-sheet NF270 campaign that we may extend to later.

---

## 5. Quick reference: which function does the loader call?

```
loader (refactored_ucb_library.py)
   │
   ├─ EC25 compensate raw conductivity to 25 °C-equivalent
   │  (using the alpha table at the top of the library file)
   │
   ├─ Call into conductivity_paper.py:
   │     single salt   →  variant_shedlovsky (then bisect to invert σ → c)
   │     mixture       →  msa_transport      (then bisect to invert σ → c)
   │
   └─ Store result as cF_exp (mM) in data_stru
```

Everything downstream of this — the Pyomo model, the figures, the contour panels — consumes `cF_exp` in mM and never sees conductivity again.

---

## 6. References

- Robinson & Stokes, *Electrolyte Solutions* (Dover, 2002): the canonical Shedlovsky derivation and the limiting-mobility tables we use for Λ₀ and αᵢ.
- Shedlovsky, T., *J. Am. Chem. Soc.* **54** (1932): the original Λ(c) expression for single salts.
- Bernard & Blum, "Mean Spherical Approximation for transport in electrolytes" (1996+): the MSA expressions implemented in `msa_transport`.
- APHA Standard Methods 2510 B: the EC25 temperature compensation we apply *before* calling Shedlovsky.

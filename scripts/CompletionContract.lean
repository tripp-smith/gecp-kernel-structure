import GECPKernelStructure.Theorems

/-!
Downstream type-level regression for the original SPEC's delivered interfaces.
Only public re-exports are imported. These examples consume existing proofs;
they neither add a research result nor replace review of their dependencies.
-/

open GECPKernelStructure
open GECPKernelStructure.GECP
open GECPKernelStructure.Fermionic

section Algebra

variable {X Y F : Type*} [Field F] {K : Kernel X Y F}

example (run : Run K) {x : X} (hx : x ∈ run.rows) (y : Y) :
    run.finalResidual x y = 0 :=
  Run.gecp_interpolates_selected_rows run hx y

example (run : Run K) {y : Y} (hy : y ∈ run.columns) (x : X) :
    run.finalResidual x y = 0 :=
  Run.gecp_interpolates_selected_cols run hy x

example (run : Run K) : run.finSelectedCore.det = run.pivots.prod :=
  gecp_core_det_eq_prod_pivots run

example (run : Run K) : run.finSelectedCore.det ≠ 0 :=
  gecp_core_nonsingular run

end Algebra

section Fermionic

example (t omega : ℝ) : fermionicKernel (1 - t) (-omega) = fermionicKernel t omega :=
  fermionicKernel_reflection t omega

example (t omega : ℝ) :
    fermionicKernel t omega = Real.exp (-(t - (1 / 2 : ℝ)) * omega) /
      (2 * Real.cosh (omega / 2)) :=
  fermionicKernel_centered t omega

example : Continuous (fun z : ℝ × ℝ => fermionicKernel z.1 z.2) :=
  fermionicKernel_continuous

example (Λ epsilon : ℝ) (hΛ : 1 ≤ Λ) (he : 0 < epsilon) (he1 : epsilon < 1) :
    ∃ terms : List SeparatedTerm,
      (terms.length : ℝ) ≤ 128 * Real.log (1 + Λ) * Real.log (128 / epsilon) ∧
      ∀ t omega : ℝ, 0 ≤ t → t ≤ 1 → -Λ ≤ omega → omega ≤ Λ →
        |fermionicKernel t omega - evalSeparated terms t omega| ≤ epsilon :=
  fermionicKernel_separatedApprox_real Λ epsilon hΛ he he1

example (Λ : ℝ) (hΛ : 0 ≤ Λ) (run : Run fermionicKernel) :
    ∃ p q : ℝ, CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -Λ ≤ omega ∧ omega ≤ Λ) run.finalResidual p q :=
  fermionicKernel_finiteRun_exists_completePivot Λ hΛ run

variable (Λ : ℝ) (hΛ : 1 ≤ Λ) (run : Run fermionicKernel)
  (complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
    (fun omega : ℝ => -Λ ≤ omega ∧ omega ≤ Λ))

example (t omega : ℝ) (ht : 0 ≤ t ∧ t ≤ 1) (hw : -Λ ≤ omega ∧ omega ≤ Λ) :
    |run.finalResidual t omega| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ)) :=
  fermionicKernel_finiteRun_error_le_exp Λ hΛ run complete t omega ht hw

example (t omega : ℝ) (ht : 0 ≤ t ∧ t ≤ 1) (hw : -Λ ≤ omega ∧ omega ≤ Λ) :
    |fermionicKernel t omega - run.crossApproximation t omega| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ)) :=
  fermionicKernel_finiteRun_cross_error_le_exp Λ hΛ run complete t omega ht hw

example (μ : MeasureTheory.Measure ℝ) (rho : ℝ → ℝ)
    (hrho : MeasureTheory.Integrable rho (μ.restrict (Set.Icc (-Λ) Λ)))
    (t : ℝ) (ht : 0 ≤ t ∧ t ≤ 1) :
    |(∫ y, rho y * fermionicKernel t y ∂μ.restrict (Set.Icc (-Λ) Λ)) -
        (∫ y, rho y * run.crossApproximation t y ∂μ.restrict (Set.Icc (-Λ) Λ))| ≤
      (1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (1 + Λ))) *
        ∫ y, |rho y| ∂μ.restrict (Set.Icc (-Λ) Λ) :=
  fermionicKernel_finiteRun_greenError_le_exp_interval Λ hΛ run complete μ rho hrho t ht

end Fermionic

section ParameterBoundary

example (Λ : ℝ) (hΛ : 0 < Λ) (run : Run fermionicKernel)
    (complete : run.CompleteOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
      (fun omega : ℝ => -Λ ≤ omega ∧ omega ≤ Λ))
    (t omega : ℝ) (ht : 0 ≤ t ∧ t ≤ 1) (hw : -Λ ≤ omega ∧ omega ≤ Λ) :
    |run.finalResidual t omega| ≤
      1024 * Real.exp (-((Real.log 2) ^ 2 / 256) * run.rows.length / Real.log (2 + Λ)) :=
  fermionicKernel_finiteRun_error_le_allCutoff_exp Λ hΛ run complete t omega ht hw

example (C c : ℝ) (hC : 0 < C) (hc : 0 < c) :
    ∃ Λ : ℝ, 0 < Λ ∧ Λ ≤ 1 ∧
      CompletePivotOn (fun t : ℝ => 0 ≤ t ∧ t ≤ 1)
        (fun omega : ℝ => -Λ ≤ omega ∧ omega ≤ Λ) fermionicKernel 0 Λ ∧
      C * Real.exp (-c / Real.log (1 + Λ)) <
        |residualUpdate fermionicKernel 0 Λ (fermionicKernel_pos 0 Λ).ne' 1 (-Λ)| :=
  fermionicKernel_smallCutoff_logarithmic_obstruction C c hC hc

example {R : Kernel ℝ ℝ ℝ} {D C : ℝ → Prop} {p q : ℝ}
    (complete : CompletePivotOn D C R p q) (hzero : R p q = 0)
    (x y : ℝ) (hx : D x) (hy : C y) : R x y = 0 :=
  residual_eq_zero_of_complete_zero complete hzero x y hx hy

end ParameterBoundary

section Baseline

open GECPKernelStructure.PositiveDefinite

variable {X : Type*} {K : Kernel X X ℝ} (run : Run K) (hK : Matrix.PosDef K)

example : run.CompleteOn (fun _ => True) (fun _ => True) ↔ DiagonalMaxRun run :=
  posDef_completeRun_iff_pivotedCholesky run hK

example (hc : run.CompleteOn (fun _ => True) (fun _ => True)) (x y : X) :
    |run.finalResidual x y| ^ 2 ≤ run.finalResidual x x * run.finalResidual y y :=
  posDef_completeRun_residual_abs_le_diag run hK hc x y

example [PseudoMetricSpace X] (hc : run.CompleteOn (fun _ => True) (fun _ => True))
    {L : NNReal} {h : ℝ} (hL : ∀ y, LipschitzWith L (fun x => K x y))
    (hfill : ∀ x, ∃ s ∈ run.rows, dist x s ≤ h) :
    ∀ x y, |run.finalResidual x y| ≤ 2 * (L : ℝ) * h :=
  posDef_completeRun_error_le_original_fillDistance run hK hc hL hfill

end Baseline

section Discretization

open GECPKernelStructure.Discretization

example {X : Type*} [PseudoMetricSpace X] {f : X → ℝ} {grid : Set X}
    {L : NNReal} {h b : ℝ} (hL : LipschitzWith L f)
    (hb : ∀ x ∈ grid, f x ≤ b) (hcover : ∀ x, ∃ z ∈ grid, dist x z ≤ h) :
    ∀ x, f x ≤ b + L * h :=
  grid_sup_le_max_add_lipschitz hL hb hcover

example {X : Type*} {f : X → ℝ} {candidate upper eta : ℝ}
    (heta : 0 ≤ eta) (hupper : ∀ x, f x ≤ upper) (hgap : eta * upper ≤ candidate) :
    ∀ x, eta * f x ≤ candidate :=
  approxPivot_of_grid_certificate heta hupper hgap

end Discretization

example : GECP.exactTwoPivotRealRun.CompleteOn (fun _ => True) (fun _ => True) :=
  GECP.exactTwoPivotRealRun_complete

example : GECP.exactTwoPivotRealRun.pivots = [2, 1 / 2] :=
  GECP.exactTwoPivotRealRun_pivots

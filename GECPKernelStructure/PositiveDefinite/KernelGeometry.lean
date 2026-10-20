import GECPKernelStructure.PositiveDefinite.ResidualPSD
import GECPKernelStructure.GECP.Interpolation
import Mathlib.Topology.Order.Compact

namespace GECPKernelStructure.PositiveDefinite

open Matrix GECP

universe u
variable {X : Type u} {K : Kernel X X ℝ}

/-- Positivity on every finite sample, including repeated nodes. -/
def KernelPosSemidef (K : Kernel X X ℝ) : Prop :=
  ∀ (n : ℕ) (nodes : Fin n → X),
    Matrix.PosSemidef (fun i j => K (nodes i) (nodes j))

namespace KernelPosSemidef

theorem of_matrix (hK : Matrix.PosSemidef K) : KernelPosSemidef K :=
  fun _ nodes => hK.submatrix nodes

theorem symmetric (hK : KernelPosSemidef K) (x y : X) : K y x = K x y := by
  simpa using (hK 2 ![x, y]).isHermitian.apply 0 1

theorem diag_nonneg (hK : KernelPosSemidef K) (x : X) : 0 ≤ K x x :=
  (hK 1 (fun _ => x)).diag_nonneg (i := 0)

theorem abs_sq_le (hK : KernelPosSemidef K) (x y : X) :
    |K x y| ^ 2 ≤ K x x * K y y := by
  simpa using posDef_gecp_residual_abs_le_diag _ (hK 2 ![x, y]) 0 1

/-- A bound on the diagonal is exactly an absolute bound on all entries. -/
theorem diagonal_bound_iff (hK : KernelPosSemidef K) (b : ℝ) :
    (∀ x, K x x ≤ b) ↔ ∀ x y, |K x y| ≤ b := by
  constructor
  · intro hb x y
    have h := (posDef_gecp_residual_sup_eq_diag_sup _ (hK 2 ![x, y]) b).mp
      (fun i => hb (![x, y] i)) 0 1
    exact h
  · intro hb x
    simpa [abs_of_nonneg (hK.diag_nonneg x)] using hb x x

theorem energy_nonneg (hK : KernelPosSemidef K) (x s : X) :
    0 ≤ K x x - 2 * K x s + K s s := by
  have h := (hK 2 ![x, s]).dotProduct_mulVec_nonneg ![1, -1]
  simp [dotProduct, mulVec, Fin.sum_univ_two, hK.symmetric x s] at h
  linarith

/-- Append the pivot to a finite sample and restrict its finite Schur update. -/
theorem residualUpdate (hK : KernelPosSemidef K) (p : X) (hp : K p p ≠ 0) :
    KernelPosSemidef (GECP.residualUpdate K p p hp) := by
  intro n nodes
  let enlarged : Fin (n + 1) → X := Fin.cons p nodes
  have hpos : 0 < K p p := lt_of_le_of_ne (hK.diag_nonneg p) (Ne.symm hp)
  have h := posDef_gecp_residual_posSemidefinite _ (hK (n + 1) enlarged) 0 hpos
  exact h.submatrix Fin.succ

/-- A positive diagonal update removes a PSD rank-one kernel. -/
theorem removedUpdate (hK : KernelPosSemidef K) (p : X) (hp : K p p ≠ 0) :
    KernelPosSemidef (fun x y => K x y - GECP.residualUpdate K p p hp x y) := by
  intro n nodes
  have h := (Matrix.posSemidef_vecMulVec_self_star (fun i => K (nodes i) p)).smul
    (inv_nonneg.mpr (hK.diag_nonneg p))
  have heq : (fun i j => K (nodes i) (nodes j) -
      GECP.residualUpdate K p p hp (nodes i) (nodes j)) =
      ((K p p)⁻¹ • Matrix.vecMulVec (fun i => K (nodes i) p)
        (star (fun i => K (nodes i) p)) : Matrix (Fin n) (Fin n) ℝ) := by
    ext i j
    simp [GECP.residualUpdate, Matrix.vecMulVec, hK.symmetric p (nodes j), div_eq_mul_inv]
    ring
  rw [heq]
  exact h

theorem add {H : Kernel X X ℝ} (hK : KernelPosSemidef K)
    (hH : KernelPosSemidef H) : KernelPosSemidef (fun x y => K x y + H x y) :=
  fun n nodes => (hK n nodes).add (hH n nodes)

end KernelPosSemidef

/-- Every successful pivot is diagonal; no maximization is required here. -/
def DiagonalRun : {K : Kernel X X ℝ} → Run K → Prop
  | _, .nil _ => True
  | _, .step p q _ tail => p = q ∧ DiagonalRun tail

/-- Finite diagonal elimination preserves both residual positivity and Loewner order. -/
theorem diagonalRun_residual_and_removed_posSemidef (run : Run K)
    (hdiag : DiagonalRun run) (hK : KernelPosSemidef K) :
    KernelPosSemidef run.finalResidual ∧
      KernelPosSemidef (fun x y => K x y - run.finalResidual x y) := by
  induction run with
  | nil K =>
    refine ⟨hK, ?_⟩
    intro n nodes
    have heq : (fun i j => K (nodes i) (nodes j) -
        (Run.nil K).finalResidual (nodes i) (nodes j)) =
        (0 : Matrix (Fin n) (Fin n) ℝ) := by ext; simp [Run.finalResidual]
    rw [heq]
    exact Matrix.PosSemidef.zero
  | @step K p q hp tail ih =>
    obtain ⟨rfl, htail⟩ := hdiag
    obtain ⟨hR, hremoved⟩ := ih htail (hK.residualUpdate p hp)
    refine ⟨hR, ?_⟩
    have hsum := (hK.removedUpdate p hp).add hremoved
    convert hsum using 1
    ext x y
    change K x y - tail.finalResidual x y = _
    ring

/-- Power-function energy comparison uses original-kernel geometry. -/
theorem diagonalRun_power_le_original_energy (run : Run K)
    (hdiag : DiagonalRun run) (hK : KernelPosSemidef K)
    {s : X} (hs : s ∈ run.rows) (x : X) :
    run.finalResidual x x ≤ K x x - 2 * K x s + K s s := by
  obtain ⟨hR, hremoved⟩ := diagonalRun_residual_and_removed_posSemidef run hdiag hK
  have h := hremoved.energy_nonneg x s
  have hrow := run.gecp_interpolates_selected_rows hs
  have hcol : run.finalResidual x s = 0 := by
    rw [← hR.symmetric x s, hrow x]
  simp only [hcol, hrow s, sub_zero] at h
  linarith

/-- A successful finite run preserves joint continuity on arbitrary topological domains. -/
theorem continuous_finalResidual_kernel [TopologicalSpace X] (run : Run K)
    (hK : Continuous (fun z : X × X => K z.1 z.2)) :
    Continuous (fun z : X × X => run.finalResidual z.1 z.2) := by
  induction run with
  | nil => exact hK
  | @step K p q hp tail ih =>
    apply ih
    exact hK.sub (((hK.comp (continuous_fst.prodMk continuous_const)).mul
      (hK.comp (continuous_const.prodMk continuous_snd))).div_const (K p q))

/-- Compact continuous PSD residuals have a diagonal complete maximum, possibly zero. -/
theorem diagonalRun_exists_diagonal_completePivot [TopologicalSpace X]
    [CompactSpace X] [Nonempty X] (run : Run K) (hdiag : DiagonalRun run)
    (hK : KernelPosSemidef K) (hcont : Continuous (fun z : X × X => K z.1 z.2)) :
    ∃ p : X, ∀ x y, |run.finalResidual x y| ≤ run.finalResidual p p := by
  have hR := (diagonalRun_residual_and_removed_posSemidef run hdiag hK).1
  have hc := (continuous_finalResidual_kernel run hcont).comp
    (continuous_id.prodMk continuous_id)
  obtain ⟨p, _, hp⟩ := isCompact_univ.exists_isMaxOn Set.univ_nonempty hc.continuousOn
  exact ⟨p, (hR.diagonal_bound_iff _).mp (fun x => hp (Set.mem_univ x))⟩

end GECPKernelStructure.PositiveDefinite

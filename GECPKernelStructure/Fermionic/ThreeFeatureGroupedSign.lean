import Mathlib.Tactic

namespace GECPKernelStructure
namespace Fermionic

/-- The first three positive Taylor features, evaluated at one row and one
column node. -/
noncomputable def positiveTaylorThreeKernel (x y : ℝ) : ℝ :=
  1 + x * y + (x * y) ^ 2 / 2

/-- The left bordered coordinate of degree zero after selecting `(x0,y0)`. -/
noncomputable def positiveTaylorThreeLeftZero (x0 x1 y0 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 - positiveTaylorThreeKernel x1 y0

/-- The left bordered coordinate of degree one. -/
noncomputable def positiveTaylorThreeLeftOne (x0 x1 y0 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 * x1 -
    x0 * positiveTaylorThreeKernel x1 y0

/-- The left bordered coordinate of degree two, including the Taylor
factorial `2!`. -/
noncomputable def positiveTaylorThreeLeftTwo (x0 x1 y0 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 * (x1 ^ 2 / 2) -
    (x0 ^ 2 / 2) * positiveTaylorThreeKernel x1 y0

/-- The right bordered coordinate of degree zero after selecting `(x0,y0)`. -/
noncomputable def positiveTaylorThreeRightZero (x0 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 - positiveTaylorThreeKernel x0 y1

/-- The right bordered coordinate of degree one. -/
noncomputable def positiveTaylorThreeRightOne (x0 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 * y1 -
    positiveTaylorThreeKernel x0 y1 * y0

/-- The right bordered coordinate of degree two. -/
noncomputable def positiveTaylorThreeRightTwo (x0 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeKernel x0 y0 * y1 ^ 2 -
    positiveTaylorThreeKernel x0 y1 * y0 ^ 2

noncomputable def positiveTaylorThreeProductZero
    (x0 x1 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeLeftZero x0 x1 y0 *
    positiveTaylorThreeRightZero x0 y0 y1

noncomputable def positiveTaylorThreeProductOne
    (x0 x1 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeLeftOne x0 x1 y0 *
    positiveTaylorThreeRightOne x0 y0 y1

noncomputable def positiveTaylorThreeProductTwo
    (x0 x1 y0 y1 : ℝ) : ℝ :=
  positiveTaylorThreeLeftTwo x0 x1 y0 *
    positiveTaylorThreeRightTwo x0 y0 y1

theorem positiveTaylorThreeProductZero_factor
    (x0 x1 y0 y1 : ℝ) :
    positiveTaylorThreeProductZero x0 x1 y0 y1 =
      (x0 - x1) * (y0 - y1) / 4 *
        (x0 * y0 * (x0 * y0 + x0 * y1 + 2) *
          (x0 * y0 + x1 * y0 + 2)) := by
  unfold positiveTaylorThreeProductZero positiveTaylorThreeLeftZero
    positiveTaylorThreeRightZero positiveTaylorThreeKernel
  ring

theorem positiveTaylorThreeProductOne_factor
    (x0 x1 y0 y1 : ℝ) :
    positiveTaylorThreeProductOne x0 x1 y0 y1 =
      (x0 - x1) * (y0 - y1) / 4 *
        ((x0 * x1 * y0 ^ 2 - 2) * (x0 ^ 2 * y0 * y1 - 2)) := by
  unfold positiveTaylorThreeProductOne positiveTaylorThreeLeftOne
    positiveTaylorThreeRightOne positiveTaylorThreeKernel
  ring

theorem positiveTaylorThreeProductTwo_factor
    (x0 x1 y0 y1 : ℝ) :
    positiveTaylorThreeProductTwo x0 x1 y0 y1 =
      (x0 - x1) * (y0 - y1) / 2 *
        ((x0 * x1 * y0 + x0 + x1) *
          (x0 * y0 * y1 + y0 + y1)) := by
  unfold positiveTaylorThreeProductTwo positiveTaylorThreeLeftTwo
    positiveTaylorThreeRightTwo positiveTaylorThreeKernel
  ring

/-- The first adjacent coordinate group factors into the two ordered gaps and
a polynomial with positive coefficients. -/
theorem positiveTaylorThree_prefix_factor
    (x0 x1 y0 y1 : ℝ) :
    positiveTaylorThreeProductZero x0 x1 y0 y1 +
        positiveTaylorThreeProductOne x0 x1 y0 y1 =
      (x0 - x1) * (y0 - y1) / 4 *
        (x0 ^ 3 * x1 * y0 ^ 3 * y1 + x0 ^ 3 * y0 ^ 3 +
          x0 ^ 3 * y0 ^ 2 * y1 + x0 ^ 2 * x1 * y0 ^ 3 +
          x0 ^ 2 * x1 * y0 ^ 2 * y1 + 4 * x0 ^ 2 * y0 ^ 2 +
          4 * x0 * y0 + 4) := by
  rw [positiveTaylorThreeProductZero_factor,
    positiveTaylorThreeProductOne_factor]
  ring

/-- The second adjacent coordinate group has the same ordered-gap factor and
another positive-coefficient polynomial. -/
theorem positiveTaylorThree_suffix_factor
    (x0 x1 y0 y1 : ℝ) :
    positiveTaylorThreeProductOne x0 x1 y0 y1 +
        positiveTaylorThreeProductTwo x0 x1 y0 y1 =
      (x0 - x1) * (y0 - y1) / 4 *
        (x0 ^ 3 * x1 * y0 ^ 3 * y1 +
          2 * x0 ^ 2 * x1 * y0 ^ 2 * y1 +
          4 * x0 * x1 * y0 * y1 + 2 * x0 * y0 + 2 * x0 * y1 +
          2 * x1 * y0 + 2 * x1 * y1 + 4) := by
  rw [positiveTaylorThreeProductOne_factor,
    positiveTaylorThreeProductTwo_factor]
  ring

/-- The first adjacent group is strictly positive at arbitrary strictly
ordered nonnegative row and column nodes. -/
theorem positiveTaylorThree_prefix_pos
    {x0 x1 y0 y1 : ℝ}
    (x0_nonneg : 0 ≤ x0) (x_order : x0 < x1)
    (y0_nonneg : 0 ≤ y0) (y_order : y0 < y1) :
    0 < positiveTaylorThreeProductZero x0 x1 y0 y1 +
      positiveTaylorThreeProductOne x0 x1 y0 y1 := by
  have x1_nonneg : 0 ≤ x1 := le_trans x0_nonneg x_order.le
  have y1_nonneg : 0 ≤ y1 := le_trans y0_nonneg y_order.le
  have gaps_pos : 0 < (x0 - x1) * (y0 - y1) :=
    mul_pos_of_neg_of_neg (sub_neg.mpr x_order) (sub_neg.mpr y_order)
  rw [positiveTaylorThree_prefix_factor]
  positivity

/-- The second adjacent group is strictly positive on the same full domain. -/
theorem positiveTaylorThree_suffix_pos
    {x0 x1 y0 y1 : ℝ}
    (x0_nonneg : 0 ≤ x0) (x_order : x0 < x1)
    (y0_nonneg : 0 ≤ y0) (y_order : y0 < y1) :
    0 < positiveTaylorThreeProductOne x0 x1 y0 y1 +
      positiveTaylorThreeProductTwo x0 x1 y0 y1 := by
  have x1_nonneg : 0 ≤ x1 := le_trans x0_nonneg x_order.le
  have y1_nonneg : 0 ≤ y1 := le_trans y0_nonneg y_order.le
  have gaps_pos : 0 < (x0 - x1) * (y0 - y1) :=
    mul_pos_of_neg_of_neg (sub_neg.mpr x_order) (sub_neg.mpr y_order)
  rw [positiveTaylorThree_suffix_factor]
  positivity

/-- The degree-zero coordinate is positive when the selected nodes themselves
are positive. -/
theorem positiveTaylorThree_productZero_pos
    {x0 x1 y0 y1 : ℝ}
    (x0_pos : 0 < x0) (x_order : x0 < x1)
    (y0_pos : 0 < y0) (y_order : y0 < y1) :
    0 < positiveTaylorThreeProductZero x0 x1 y0 y1 := by
  have x1_pos : 0 < x1 := lt_trans x0_pos x_order
  have y1_pos : 0 < y1 := lt_trans y0_pos y_order
  have gaps_pos : 0 < (x0 - x1) * (y0 - y1) :=
    mul_pos_of_neg_of_neg (sub_neg.mpr x_order) (sub_neg.mpr y_order)
  rw [positiveTaylorThreeProductZero_factor]
  positivity

/-- The degree-two coordinate is positive even when an initial node is zero. -/
theorem positiveTaylorThree_productTwo_pos
    {x0 x1 y0 y1 : ℝ}
    (x0_nonneg : 0 ≤ x0) (x_order : x0 < x1)
    (y0_nonneg : 0 ≤ y0) (y_order : y0 < y1) :
    0 < positiveTaylorThreeProductTwo x0 x1 y0 y1 := by
  have x1_pos : 0 < x1 := lt_of_le_of_lt x0_nonneg x_order
  have y1_pos : 0 < y1 := lt_of_le_of_lt y0_nonneg y_order
  have gaps_pos : 0 < (x0 - x1) * (y0 - y1) :=
    mul_pos_of_neg_of_neg (sub_neg.mpr x_order) (sub_neg.mpr y_order)
  rw [positiveTaylorThreeProductTwo_factor]
  positivity

/-- The first nontrivial contiguous step is positive for every nonnegative
weight at positive selected nodes. -/
theorem positiveTaylorThree_firstStep_pos
    {x0 x1 y0 y1 weight : ℝ}
    (x0_pos : 0 < x0) (x_order : x0 < x1)
    (y0_pos : 0 < y0) (y_order : y0 < y1)
    (weight_nonneg : 0 ≤ weight) :
    0 < positiveTaylorThreeProductZero x0 x1 y0 y1 +
      weight * (positiveTaylorThreeProductOne x0 x1 y0 y1 +
        positiveTaylorThreeProductTwo x0 x1 y0 y1) := by
  have head_pos := positiveTaylorThree_productZero_pos x0_pos x_order y0_pos y_order
  have tail_pos := positiveTaylorThree_suffix_pos x0_pos.le x_order y0_pos.le y_order
  positivity

/-- The second nontrivial contiguous step is positive for every nonnegative
weight, including at zero boundary nodes. -/
theorem positiveTaylorThree_secondStep_pos
    {x0 x1 y0 y1 weight : ℝ}
    (x0_nonneg : 0 ≤ x0) (x_order : x0 < x1)
    (y0_nonneg : 0 ≤ y0) (y_order : y0 < y1)
    (weight_nonneg : 0 ≤ weight) :
    0 < positiveTaylorThreeProductZero x0 x1 y0 y1 +
        positiveTaylorThreeProductOne x0 x1 y0 y1 +
      weight * positiveTaylorThreeProductTwo x0 x1 y0 y1 := by
  have head_pos := positiveTaylorThree_prefix_pos x0_nonneg x_order y0_nonneg y_order
  have tail_pos := positiveTaylorThree_productTwo_pos x0_nonneg x_order y0_nonneg y_order
  positivity

end Fermionic
end GECPKernelStructure

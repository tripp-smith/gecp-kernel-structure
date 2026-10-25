//! Optional reference time-design accelerator. No mathematical certificate is claimed.
//! Every MP operation uses 128 bits; input normalization and exports are binary64.
use gmp_mpfr_sys::{gmp, mpfr};
use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyBytes, PyDict, PyFloat, PyInt, PyList, PyString};
use rug::Float;
use std::collections::HashSet;
use std::ffi::CStr;

const BITS: u32 = 128;
const MAX_RANK: usize = 128;
const MAX_ROWS: usize = 16_384;
const MAX_CELLS: usize = 1 << 20;
const MAX_TOKEN: usize = 96;
const MAX_FACTOR_EXP: i32 = 65_536;
const REQUIRED_EXP_RANGE: i64 = 1 << 28;

/// MPFR contexts are thread-local. Do not alter exponent limits or caller flags.
struct FlagsGuard {
    saved: mpfr::flags_t,
}
impl FlagsGuard {
    fn begin() -> Result<Self, String> {
        // SAFETY: MPFR is initialized by its linked library; scalar queries have
        // no pointer arguments. Require TLS before detaching from Python.
        unsafe {
            if mpfr::buildopt_tls_p() == 0 {
                return Err("native backend requires a thread-safe MPFR build".into());
            }
            if (mpfr::get_emin() as i64) > -REQUIRED_EXP_RANGE
                || (mpfr::get_emax() as i64) < REQUIRED_EXP_RANGE
            {
                return Err(
                    "caller MPFR exponent range is too narrow for native evaluation".into(),
                );
            }
            let saved = mpfr::flags_save();
            mpfr::clear_flags();
            Ok(Self { saved })
        }
    }
    fn check(&self) -> Result<(), String> {
        // SAFETY: scalar, thread-local MPFR flag queries.
        let bad = unsafe {
            mpfr::underflow_p() != 0
                || mpfr::overflow_p() != 0
                || mpfr::nanflag_p() != 0
                || mpfr::divby0_p() != 0
                || mpfr::erangeflag_p() != 0
        };
        if bad {
            Err("native MPFR arithmetic raised an underflow, overflow, invalid, division or range flag".into())
        } else {
            Ok(())
        }
    }
    fn clear_after_export(&self) {
        // Final binary64 underflow-to-zero follows the existing export policy.
        // It cannot erase an arithmetic flag: check() is required before export.
        unsafe { mpfr::clear_flags() }
    }
}
impl Drop for FlagsGuard {
    fn drop(&mut self) {
        // SAFETY: restore exactly the thread's actual prior flags on all exits.
        unsafe {
            mpfr::clear_flags();
            mpfr::flags_set(self.saved);
        }
    }
}

fn list<'py>(value: &'py Bound<'py, PyAny>, name: &str) -> PyResult<&'py Bound<'py, PyList>> {
    if !value.is_exact_instance_of::<PyList>() {
        return Err(PyTypeError::new_err(format!(
            "{name} must be a builtin list"
        )));
    }
    value
        .cast::<PyList>()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be a builtin list")))
}
fn real(value: &Bound<'_, PyAny>, name: &str) -> PyResult<f64> {
    if value.is_instance_of::<PyBool>() {
        return Err(PyTypeError::new_err(format!("{name} must not be boolean")));
    }
    // Arbitrary __float__ hooks could mutate an input list during validation,
    // invalidating its checked rank before the owned computation starts.
    if !value.is_exact_instance_of::<PyFloat>() && !value.is_exact_instance_of::<PyInt>() {
        return Err(PyTypeError::new_err(format!(
            "{name} must be a builtin float or int"
        )));
    }
    let number = value
        .extract::<f64>()
        .map_err(|_| PyTypeError::new_err(format!("{name} must be real")))?;
    if !number.is_finite() {
        return Err(PyValueError::new_err(format!("{name} must be finite")));
    }
    Ok(number)
}
fn factor(value: &Bound<'_, PyAny>) -> PyResult<Float> {
    if !value.is_exact_instance_of::<PyString>() {
        return Err(PyTypeError::new_err(
            "LU factors must be builtin decimal strings",
        ));
    }
    let token = value
        .cast::<PyString>()
        .map_err(|_| PyTypeError::new_err("LU factors must be decimal strings"))?;
    // Check the exact string's character count before requesting UTF-8; the
    // latter may otherwise allocate an unbounded copy of a non-ASCII input.
    if token.len()? > MAX_TOKEN {
        return Err(PyValueError::new_err("LU factor token exceeds its budget"));
    }
    let text = token.to_str()?;
    if text.is_empty()
        || text.len() > MAX_TOKEN
        || !text
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'.' | b'+' | b'-' | b'e' | b'E'))
    {
        return Err(PyValueError::new_err(
            "LU factor must be a bounded finite decimal token",
        ));
    }
    if let Some(offset) = text.find(['e', 'E']) {
        let exponent = text[offset + 1..]
            .parse::<i32>()
            .map_err(|_| PyValueError::new_err("invalid decimal factor exponent"))?;
        if !(-MAX_FACTOR_EXP..=MAX_FACTOR_EXP).contains(&exponent) {
            return Err(PyValueError::new_err(
                "decimal factor exponent exceeds native range",
            ));
        }
    }
    let parsed =
        Float::parse(text).map_err(|_| PyValueError::new_err("invalid decimal LU factor"))?;
    let number = Float::with_val(BITS, parsed);
    if !number.is_finite()
        || number
            .get_exp()
            .is_some_and(|e| !(-MAX_FACTOR_EXP..=MAX_FACTOR_EXP).contains(&e))
    {
        return Err(PyValueError::new_err(
            "LU factor is nonfinite or outside native exponent range",
        ));
    }
    Ok(number)
}
fn matrix(value: &Bound<'_, PyAny>, rank: usize, name: &str) -> PyResult<Vec<Vec<Float>>> {
    let rows = list(value, name)?;
    if rows.len() != rank {
        return Err(PyValueError::new_err(format!("{name} must have rank rows")));
    }
    let mut result = Vec::with_capacity(rank);
    for row in rows.iter() {
        let columns = list(&row, name)?;
        if columns.len() != rank {
            return Err(PyValueError::new_err(format!(
                "{name} must be square without ragged rows"
            )));
        }
        result.push(
            columns
                .iter()
                .map(|v| factor(&v))
                .collect::<PyResult<Vec<_>>>()?,
        );
    }
    Ok(result)
}

#[pyclass(frozen, module = "kernelgecp._native")]
struct ReferenceCardinal128 {
    omegas: Vec<Float>,
    order: Vec<usize>,
    lower: Vec<Vec<Float>>,
    upper: Vec<Vec<Float>>,
}
impl ReferenceCardinal128 {
    fn evaluate(&self, times: &[f64]) -> Result<Vec<u8>, String> {
        let guard = FlagsGuard::begin()?;
        let rank = self.omegas.len();
        let mut bytes = Vec::with_capacity(times.len() * rank * 8);
        for time in times {
            let t = Float::with_val(BITS, *time);
            let rhs: Vec<Float> = self
                .omegas
                .iter()
                .map(|w| {
                    let exponent = if w >= &0 {
                        -Float::with_val(BITS, &t * w)
                    } else {
                        Float::with_val(BITS, Float::with_val(BITS, 1 - &t) * w)
                    };
                    let numerator = exponent.exp();
                    let denominator_exponent = if w >= &0 { -w.clone() } else { w.clone() };
                    let denominator = Float::with_val(BITS, 1 + denominator_exponent.exp());
                    Float::with_val(BITS, numerator / denominator)
                })
                .collect();
            let mut forward = vec![Float::new(BITS); rank];
            let mut solution = vec![Float::new(BITS); rank];
            for i in 0..rank {
                let products: Vec<Float> = (0..i)
                    .map(|j| Float::with_val(BITS, &self.lower[i][j] * &forward[j]))
                    .collect();
                let subtotal = Float::with_val(BITS, Float::sum(products.iter()));
                let residual = Float::with_val(BITS, &rhs[self.order[i]] - subtotal);
                forward[i] = Float::with_val(BITS, residual / &self.lower[i][i]);
            }
            for i in (0..rank).rev() {
                let products: Vec<Float> = ((i + 1)..rank)
                    .map(|j| Float::with_val(BITS, &self.upper[i][j] * &solution[j]))
                    .collect();
                let subtotal = Float::with_val(BITS, Float::sum(products.iter()));
                let residual = Float::with_val(BITS, &forward[i] - subtotal);
                solution[i] = Float::with_val(BITS, residual / &self.upper[i][i]);
            }
            guard.check()?;
            for number in solution {
                if !number.is_finite() {
                    return Err("native solve produced a nonfinite MPFR value".into());
                }
                let exported = number.to_f64();
                if !exported.is_finite() {
                    return Err("native solve exceeds finite binary64 export range".into());
                }
                bytes.extend_from_slice(&exported.to_le_bytes());
            }
            guard.clear_after_export();
        }
        Ok(bytes)
    }
}
#[pymethods]
impl ReferenceCardinal128 {
    #[new]
    fn new(
        precision_bits: &Bound<'_, PyAny>,
        omega_nodes: &Bound<'_, PyAny>,
        order: &Bound<'_, PyAny>,
        lower: &Bound<'_, PyAny>,
        upper: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        if precision_bits.is_instance_of::<PyBool>()
            || !precision_bits.is_exact_instance_of::<PyInt>()
            || precision_bits.extract::<u32>().ok() != Some(BITS)
        {
            return Err(PyValueError::new_err(
                "native reference backend requires exactly 128 integer precision bits",
            ));
        }
        let raw_omegas = list(omega_nodes, "omega_nodes")?;
        let rank = raw_omegas.len();
        if !(1..=MAX_RANK).contains(&rank) {
            return Err(PyValueError::new_err("native rank must lie in [1,128]"));
        }
        let mut frequencies = Vec::with_capacity(rank);
        let mut seen = HashSet::with_capacity(rank);
        for value in raw_omegas.iter() {
            let w = real(&value, "omega node")?;
            if w.abs() > 16_384.0 || !seen.insert(if w == 0.0 { 0 } else { w.to_bits() }) {
                return Err(PyValueError::new_err(
                    "omega nodes must be distinct with absolute value <=16384",
                ));
            }
            frequencies.push(w);
        }
        let raw_order = list(order, "order")?;
        if raw_order.len() != rank {
            return Err(PyValueError::new_err(
                "row order must be a rank-sized permutation",
            ));
        }
        let mut ordering = Vec::with_capacity(rank);
        let mut used = vec![false; rank];
        for value in raw_order.iter() {
            if value.is_instance_of::<PyBool>() || !value.is_exact_instance_of::<PyInt>() {
                return Err(PyTypeError::new_err(
                    "row order must contain nonboolean integers",
                ));
            }
            let i = value
                .extract::<usize>()
                .map_err(|_| PyValueError::new_err("invalid row permutation index"))?;
            if i >= rank || used[i] {
                return Err(PyValueError::new_err("row order must be a permutation"));
            }
            used[i] = true;
            ordering.push(i);
        }
        let guard = FlagsGuard::begin().map_err(PyRuntimeError::new_err)?;
        let lo = matrix(lower, rank, "lower")?;
        let up = matrix(upper, rank, "upper")?;
        for i in 0..rank {
            if lo[i][i] != 1 || up[i][i] == 0 {
                return Err(PyValueError::new_err(
                    "LU requires a unit lower diagonal and nonzero upper pivots",
                ));
            }
            for j in 0..rank {
                if (j > i && lo[i][j] != 0) || (j < i && up[i][j] != 0) {
                    return Err(PyValueError::new_err("LU factors must be triangular"));
                }
            }
        }
        let omegas = frequencies
            .iter()
            .map(|w| Float::with_val(BITS, *w))
            .collect();
        guard.check().map_err(PyRuntimeError::new_err)?;
        Ok(Self {
            omegas,
            order: ordering,
            lower: lo,
            upper: up,
        })
    }
    #[getter]
    fn rank(&self) -> usize {
        self.omegas.len()
    }
    /// Evaluate at already-normalized binary64 t=tau/beta. Output is little-endian f64.
    fn design_tau<'py>(
        &self,
        py: Python<'py>,
        normalized_times: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyBytes>> {
        let raw = list(normalized_times, "normalized_times")?;
        let count = raw.len();
        if count > MAX_ROWS || count > MAX_CELLS / self.omegas.len() {
            return Err(PyValueError::new_err(
                "native query exceeds row/work budget",
            ));
        }
        let mut times = Vec::with_capacity(count);
        for value in raw.iter() {
            let time = real(&value, "normalized time")?;
            if !(0.0..=1.0).contains(&time) {
                return Err(PyValueError::new_err("normalized time must lie in [0,1]"));
            }
            times.push(time);
        }
        let output = py
            .detach(|| self.evaluate(&times))
            .map_err(PyRuntimeError::new_err)?;
        Ok(PyBytes::new(py, &output))
    }
}
#[pyfunction]
fn build_info<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let info = PyDict::new(py);
    info.set_item("protocol_version", 1)?;
    info.set_item("working_precision_bits", BITS)?;
    info.set_item("package_version", env!("CARGO_PKG_VERSION"))?;
    info.set_item("build_target", env!("KERNELGECP_BUILD_TARGET"))?;
    info.set_item("source_revision", env!("KERNELGECP_BUILD_SOURCE_REVISION"))?;
    info.set_item("max_rank", MAX_RANK)?;
    info.set_item("max_rows_per_call", MAX_ROWS)?;
    info.set_item("max_cells_per_call", MAX_CELLS)?;
    info.set_item("rug_version", "1.30.0")?;
    info.set_item("pyo3_version", "0.27.2")?;
    info.set_item(
        "export",
        "binary64_little_endian; final underflow-to-zero allowed",
    )?;
    // SAFETY: version strings are static, null-terminated strings from linked libraries.
    unsafe {
        info.set_item(
            "mpfr_version",
            CStr::from_ptr(mpfr::get_version())
                .to_string_lossy()
                .as_ref(),
        )?;
        info.set_item(
            "gmp_version",
            CStr::from_ptr(gmp::version).to_string_lossy().as_ref(),
        )?;
        info.set_item("mpfr_tls", mpfr::buildopt_tls_p() != 0)?;
    }
    Ok(info)
}
#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("NATIVE_PROTOCOL_VERSION", 1)?;
    m.add("WORKING_PRECISION_BITS", BITS)?;
    m.add_class::<ReferenceCardinal128>()?;
    m.add_function(wrap_pyfunction!(build_info, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn scalar(omega: f64, pivot: Float) -> ReferenceCardinal128 {
        ReferenceCardinal128 {
            omegas: vec![Float::with_val(BITS, omega)],
            order: vec![0],
            lower: vec![vec![Float::with_val(BITS, 1)]],
            upper: vec![vec![pivot]],
        }
    }
    #[test]
    fn scalar_kernel_and_little_endian_export() {
        let handle = scalar(0.0, Float::with_val(BITS, 1));
        let output = handle.evaluate(&[0.0, 0.5, 1.0]).expect("valid scalar");
        let expected: Vec<u8> = [0.5_f64, 0.5, 0.5]
            .iter()
            .flat_map(|x| x.to_le_bytes())
            .collect();
        assert_eq!(output, expected);
        assert!(handle.evaluate(&[]).expect("empty query").is_empty());
    }
    #[test]
    fn prior_flags_restored_on_success_and_export_failure() {
        unsafe {
            let prior = mpfr::flags_save();
            mpfr::clear_flags();
            mpfr::set_underflow();
            mpfr::set_inexflag();
            let expected = mpfr::flags_save();
            let handle = scalar(0.0, Float::with_val(BITS, 1));
            assert!(handle.evaluate(&[0.3]).is_ok());
            assert_eq!(mpfr::flags_save(), expected);
            let tiny = Float::with_val(BITS, Float::parse("1e-1000").expect("bounded token"));
            // Export overflows binary64, but the MPFR computation stays finite.
            assert!(scalar(0.0, tiny).evaluate(&[0.3]).is_err());
            assert_eq!(mpfr::flags_save(), expected);
            mpfr::clear_flags();
            mpfr::flags_set(prior);
        }
    }
    #[test]
    fn narrow_caller_exponents_rejected_without_mutation() {
        unsafe {
            let low = mpfr::get_emin();
            let high = mpfr::get_emax();
            let flags = mpfr::flags_save();
            assert_eq!(mpfr::set_emin(-1024), 0);
            assert!(FlagsGuard::begin().is_err());
            assert_eq!(mpfr::get_emin(), -1024);
            assert_eq!(mpfr::get_emax(), high);
            assert_eq!(mpfr::flags_save(), flags);
            assert_eq!(mpfr::set_emin(low), 0);
        }
    }
    #[test]
    fn arithmetic_flags_fail_closed_then_restore() {
        unsafe {
            let prior = mpfr::flags_save();
            {
                let guard = FlagsGuard::begin().expect("normal context");
                mpfr::set_underflow();
                assert!(guard.check().is_err());
            }
            assert_eq!(mpfr::flags_save(), prior);
        }
    }
    #[test]
    fn binary64_underflow_export_allowed_and_does_not_leak_flags() {
        unsafe {
            let pivot = Float::with_val(BITS, Float::parse("1e1000").expect("bounded token"));
            let prior = mpfr::flags_save();
            let output = scalar(0.0, pivot)
                .evaluate(&[0.3])
                .expect("finite MP arithmetic");
            assert_eq!(output, 0.0_f64.to_le_bytes());
            assert_eq!(mpfr::flags_save(), prior);
        }
    }
}

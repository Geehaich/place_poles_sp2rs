use nalgebra::allocator::{Allocator, Reallocator};
use nalgebra::{DefaultAllocator, Matrix, Unit};
use nalgebra::{ComplexField,RealField};
use nalgebra::dimension::{ Dim, DimMin, Const, DimMinimum};
use nalgebra::base::{OVector, OMatrix , DMatrix};
use nalgebra::storage::{Storage, StorageMut};
use nalgebra::linalg::householder;
use nalgebra::geometry::Reflection;
use nalgebra::constraint::{SameNumberOfRows, ShapeConstraint};
use num_traits::Zero;
use std::mem::MaybeUninit;



/// The QR decomposition of a general matrix. Same code as nalgebra's but ignores sparseness when generating Q and R.
#[cfg_attr(feature = "serde-serialize-no-std", derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "serde-serialize-no-std",
    serde(bound(serialize = "DefaultAllocator: Allocator<T, R, C> +
                           Allocator<T, DimMinimum<R, C>>,
         OMatrix<T, R, C>: Serialize,
         OVector<T, DimMinimum<R, C>>: Serialize"))
)]
#[cfg_attr(
    feature = "serde-serialize-no-std",
    serde(bound(deserialize = "DefaultAllocator: Allocator<T, R, C> +
                           Allocator<T, DimMinimum<R, C>>,
         OMatrix<T, R, C>: Deserialize<'de>,
         OVector<T, DimMinimum<R, C>>: Deserialize<'de>"))
)]
#[derive(Clone, Debug)]
pub struct full_QR<T: ComplexField, R: DimMin<C>, C: Dim>
where
    DefaultAllocator: Allocator<T, R, C> + Allocator<T, DimMinimum<R, C>>,
{
    pub qr: OMatrix<T, R, C>,
    pub diag: OVector<T, DimMinimum<R, C>>,
}


impl<T: ComplexField, R: DimMin<C>, C: Dim> full_QR<T, R, C>
where
    DefaultAllocator: Allocator<T, R, C> + Allocator<T, R> + Allocator<T, DimMinimum<R, C>>,
{
    /// Computes the QR decomposition using householder reflections.
    pub fn new(mut matrix: OMatrix<T, R, C>) -> Self {
        let (nrows, ncols) = matrix.shape_generic();
        let min_nrows_ncols = nrows.min(ncols);

        if min_nrows_ncols.value() == 0 {
            return full_QR {
                qr: matrix,
                diag: Matrix::zeros_generic(min_nrows_ncols, Const::<1>),
            };
        }

        let mut diag = Matrix::uninit(min_nrows_ncols, Const::<1>);

        for i in 0..min_nrows_ncols.value() {
            diag[i] =
                MaybeUninit::new(householder::clear_column_unchecked(&mut matrix, i, 0, None));
}

// Safety: diag is now fully initialized.
let diag = unsafe { diag.assume_init() };
full_QR { qr: matrix, diag }
}

/// Retrieves the upper trapezoidal submatrix `R` of this decomposition.
#[inline]
#[must_use]
pub fn r(&self) -> OMatrix<T, DimMinimum<R, C>, C>
where
DefaultAllocator: Allocator<T, DimMinimum<R, C>, C>,
{
let (nrows, ncols) = self.qr.shape_generic();
let mut res = self.qr.rows_generic(0, nrows.min(C::from_usize(usize::MAX))).upper_triangle();
res.set_partial_diagonal(self.diag.iter().map(|e| T::from_real(e.clone().modulus())));
res
}

/// Retrieves the upper trapezoidal submatrix `R` of this decomposition.
///
/// This is usually faster than `r` but consumes `self`.
#[inline]
pub fn unpack_r(self) -> OMatrix<T, DimMinimum<R, C>, C>
where
DefaultAllocator: Reallocator<T, R, C, DimMinimum<R, C>, C>,
{
let (nrows, ncols) = self.qr.shape_generic();
let mut res = self.qr.resize_generic(nrows.min(C::from_usize(usize::MAX)), ncols, T::zero());
res.fill_lower_triangle(T::zero(), 1);
res.set_partial_diagonal(self.diag.iter().map(|e| T::from_real(e.clone().modulus())));
res
}

/// Computes the orthogonal matrix `Q` of this decomposition.
#[must_use]
pub fn q(&self) -> DMatrix<T>
where
DefaultAllocator: Allocator<T, R, DimMinimum<R, C>>,
{
let (nrows, ncols) = self.qr.shape();

let mut res = DMatrix::identity(nrows, ncols.max(nrows));
let dim = self.diag.len();

for i in (0..dim).rev() {
    let axis = self.qr.view_range(i.., i);
    // TODO: sometimes, the axis might have a zero magnitude.
    let refl = Reflection::new(Unit::new_unchecked(axis), T::zero());

    let mut res_rows = res.view_range_mut(i.., i..);
    refl.reflect_with_sign(&mut res_rows, self.diag[i].clone().signum());
}

res.to_owned()
}

/// Unpacks this decomposition into its two matrix factors.
pub fn unpack(
self
) -> (
DMatrix<T>,
OMatrix<T, DimMinimum<R, C>, C>,
)
where
DimMinimum<R, C>: DimMin<C, Output = DimMinimum<R, C>>,
DefaultAllocator:
    Allocator<T, R, DimMinimum<R, C>> + Reallocator<T, R, C, DimMinimum<R, C>, C>,
{
(self.q(), self.unpack_r())
}

#[doc(hidden)]
pub fn qr_internal(&self) -> &OMatrix<T, R, C> {
&self.qr
}


#[must_use]
pub(crate) fn diag_internal(&self) -> &OVector<T, DimMinimum<R, C>> {
&self.diag
}

/// Multiplies the provided matrix by the transpose of the `Q` matrix of this decomposition.
pub fn q_tr_mul<R2: Dim, C2: Dim, S2>(&self, rhs: &mut Matrix<T, R2, C2, S2>)
// TODO: do we need a static constraint on the number of rows of rhs?
where
S2: StorageMut<T, R2, C2>,
{
let dim = self.diag.len();

for i in 0..dim {
    let axis = self.qr.view_range(i.., i);
    let refl = Reflection::new(Unit::new_unchecked(axis), T::zero());

    let mut rhs_rows = rhs.rows_range_mut(i..);
    refl.reflect_with_sign(&mut rhs_rows, self.diag[i].clone().signum().conjugate());
}
}

}

impl<T: ComplexField, D: DimMin<D, Output = D>> full_QR<T, D, D>
where
    DefaultAllocator: Allocator<T, D, D> + Allocator<T, D>,
{
    /// Solves the linear system `self * x = b`, where `x` is the unknown to be determined.
    ///
    /// Returns `None` if `self` is not invertible.
    #[must_use = "Did you mean to use solve_mut()?"]
    pub fn solve<R2: Dim, C2: Dim, S2>(
        &self,
        b: &Matrix<T, R2, C2, S2>,
    ) -> Option<OMatrix<T, R2, C2>>
    where
        S2: Storage<T, R2, C2>,
        ShapeConstraint: SameNumberOfRows<R2, D>,
        DefaultAllocator: Allocator<T, R2, C2>,
    {
        let mut res = b.clone_owned();

        if self.solve_mut(&mut res) {
            Some(res)
        } else {
            None
        }
    }

    /// Solves the linear system `self * x = b`, where `x` is the unknown to be determined.
    ///
    /// If the decomposed matrix is not invertible, this returns `false` and its input `b` is
    /// overwritten with garbage.
    pub fn solve_mut<R2: Dim, C2: Dim, S2>(&self, b: &mut Matrix<T, R2, C2, S2>) -> bool
    where
        S2: StorageMut<T, R2, C2>,
        ShapeConstraint: SameNumberOfRows<R2, D>,
    {
        assert_eq!(
            self.qr.nrows(),
            b.nrows(),
            "QR solve matrix dimension mismatch."
        );
        assert!(
            self.qr.is_square(),
            "QR solve: unable to solve a non-square system."
        );

        self.q_tr_mul(b);
        self.solve_upper_triangular_mut(b)
    }

    // TODO: duplicate code from the `solve` module.
    fn solve_upper_triangular_mut<R2: Dim, C2: Dim, S2>(
        &self,
        b: &mut Matrix<T, R2, C2, S2>,
    ) -> bool
    where
        S2: StorageMut<T, R2, C2>,
        ShapeConstraint: SameNumberOfRows<R2, D>,
    {
        let dim = self.qr.nrows();

        for k in 0..b.ncols() {
            let mut b = b.column_mut(k);
            for i in (0..dim).rev() {
                let coeff;

                unsafe {
                    let diag = self.diag.vget_unchecked(i).clone().modulus();

                    if diag.is_zero() {
                        return false;
                    }

                    coeff = b.vget_unchecked(i).clone().unscale(diag);
                    *b.vget_unchecked_mut(i) = coeff.clone();
                }

                b.rows_range_mut(..i)
                    .axpy(-coeff, &self.qr.view_range(..i, i), T::one());
            }
        }

        true
    }

    /// Computes the inverse of the decomposed matrix.
    ///
    /// Returns `None` if the decomposed matrix is not invertible.
    #[must_use]
    pub fn try_inverse(&self) -> Option<OMatrix<T, D, D>> {
        assert!(
            self.qr.is_square(),
            "QR inverse: unable to compute the inverse of a non-square matrix."
        );

        // TODO: is there a less naive method ?
        let (nrows, ncols) = self.qr.shape_generic();
        let mut res = OMatrix::identity_generic(nrows, ncols);

        if self.solve_mut(&mut res) {
            Some(res)
        } else {
            None
        }
    }

    /// Indicates if the decomposed matrix is invertible.
    #[must_use]
    pub fn is_invertible(&self) -> bool {
        assert!(
            self.qr.is_square(),
            "QR: unable to test the invertibility of a non-square matrix."
        );

        for i in 0..self.diag.len() {
            if self.diag[i].is_zero() {
                return false;
            }
        }

        true
    }

    // /// Computes the determinant of the decomposed matrix.
    // pub fn determinant(&self) -> T {
    //     let dim = self.qr.nrows();
    //     assert!(self.qr.is_square(), "QR determinant: unable to compute the determinant of a non-square matrix.");

    //     let mut res = T::one();
    //     for i in 0 .. dim {
    //         res *= unsafe { *self.diag.vget_unchecked(i) };
    //     }

    //     res self.q_determinant()
    // }
}
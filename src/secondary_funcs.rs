use core::panic;


use na::{Complex, ComplexField, DMatrix, DVector, RealField };
use num_traits::real;

use super::EPoleMethod;

pub fn real_to_cpx<T:RealField + Copy>(r_mat: &DMatrix<T>) -> DMatrix<Complex<T>>
{
    let mut c_mat = DMatrix::<Complex<T>>::zeros(r_mat.shape().0,r_mat.shape().1);
    let zipped = std::iter::zip(c_mat.iter_mut(),r_mat.iter());
    
    for (c,r) in zipped
    {
        (*c).re = *r;
    }
    
    c_mat
}

pub fn cpx_to_real<T : RealField + Copy>(cp_mat : &DMatrix<Complex<T>>) -> DMatrix<T>
{
    let mut real_part = DMatrix::<T>::zeros(cp_mat.shape().0,cp_mat.shape().1);
    let zipped = std::iter::zip(real_part.iter_mut(),cp_mat.iter());
    
    for (r,cp) in zipped
    {
        *r = (*cp).re;
    }
    
    real_part
}

//check validity of inputs for the place_poles function. less checks than in Python since Rust less dynamic.
pub fn valid_inputs(mat_a : &DMatrix<f64> , poles : &DVector<Complex<f64>>, rtol :f64, maxiter: u32 , method : &EPoleMethod)
-> DVector<Complex<f64>>
{
    if !mat_a.is_square() {panic!("A must be square");}
    if poles.shape().0 != mat_a.shape().0 {panic!( "Dimension mismatch : A has {} poles, {} provided",mat_a.shape().0,poles.shape().0);}
    if maxiter < 1 {        panic!("maxiter must be at least equal to 1") ;}
    if rtol > 1.0 {panic!("rtol can't be greater than 1");}
    
    if !poles.iter().all(|x| x.im==0.0) && matches!(method,EPoleMethod::KNV0) {panic!("KNV0 Method only works with Real poles")};
    
    order_complex_poles(poles)
}

pub fn order_complex_poles(poles : &DVector::<Complex<f64>>) -> DVector<Complex<f64>>
{
    
    let mut real_poles = Vec::new();
    let mut cpx_poles = Vec::new();
    for i in 0..poles.len()
    {
        if poles[i].im == 0.0 { real_poles.push(poles[i])}
        else if poles[i].im <0.0 {cpx_poles.push(poles[i]);}
    }
    real_poles.sort_by(|x,y| x.re.partial_cmp(&y.re).unwrap());
    cpx_poles.sort_by(|x,y| x.re.partial_cmp(&y.re).unwrap());

    let mut valid_cpx_poles = Vec::with_capacity(cpx_poles.len());
    for i in 0..cpx_poles.len()
    {
        for j in i..cpx_poles.len()
        {
            if cpx_poles[i] == cpx_poles[j].conjugate() //add pole if conjugate present in list
            {
                if cpx_poles[i].im < 0.0 {valid_cpx_poles.push(cpx_poles[i]); valid_cpx_poles.push(cpx_poles[j]);} //respect cpx lexicographic sorting order.
                else {valid_cpx_poles.push(cpx_poles[j]); valid_cpx_poles.push(cpx_poles[i]);}
            }
        }
    }
    
    let ordered_poles = [real_poles,valid_cpx_poles].concat();
    if ordered_poles.len() != poles.len() {panic! ("Complex poles must come with their conjugates");}
    
    DVector::from_vec(ordered_poles)
}


pub fn is_close_cpxcst(mat : &DMatrix<Complex<f64>>, cst : Complex<f64>,tol : f64) ->bool
{
    let mut res = true;
    for elem in mat.iter()
    {
        res &= (*elem-cst).modulus() < tol;
        if res {break;}
    }
    res
}

#[macro_export]
macro_rules! vstack {
    ( $( $matrix:expr ),* ) => {
        {
            let mut ref_row_vec = Vec::new();
            let mut cols : usize = 0;
            let mut total_rows = 0;
            $(
                if cols==0 { cols = $matrix.shape().1; }
                else {assert_eq!($matrix.shape().1,cols, "Dimension mismatch, one array has an incorrect # of rows");}
                ref_row_vec.push($matrix);
                total_rows += $matrix.shape().0;
            )*
            
            let mut stacked = DMatrix::zeros(total_rows,cols);
            let mut rows_ind = 0;
            for i in 0..ref_row_vec.len()
            {
                stacked.rows_mut(rows_ind,ref_row_vec[i].shape().0).copy_from(&ref_row_vec[i]);
                rows_ind += ref_row_vec[i].shape().0;
            }
            stacked
            
        }
    };
}


#[macro_export]
macro_rules! hstack {
    ( $( $matrix:expr ),* ) => {
        {
            let mut ref_vec = Vec::new();
            let mut rows : usize = 0;
            let mut total_cols = 0;
            $(
                if rows==0 { rows = $matrix.shape().0; }
                else {assert_eq!($matrix.shape().0,rows, "Dimension mismatch, one array has an incorrect # of columns");}
                ref_vec.push(&$matrix);
                total_cols += $matrix.shape().1;
            )*
            
            let mut stacked = DMatrix::zeros(rows,total_cols);
            let mut cols_ind = 0;
            for i in 0..ref_vec.len()
            {
                stacked.columns_mut(cols_ind,ref_vec[i].shape().1).copy_from(ref_vec[i]);
                cols_ind += ref_vec[i].shape().1;
            }
            stacked
            
        }
    }
}



use core::panic;
use num_traits::{Num, Zero};


use na::{Complex, ComplexField, DMatrix, DVector, RealField };

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

pub fn real_to_cpx_vec<T:RealField + Copy>(r_mat: &DVector<T>) -> DVector<Complex<T>>
{
    let mut c_mat = DVector::<Complex<T>>::zeros(r_mat.shape().0);
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

//get real and imaginary parts of a complex array
pub fn  cpx_decompose<T:RealField+Num+Zero+Copy>(cp_mat : &DMatrix<Complex<T>>) ->(DMatrix<T>,DMatrix<T>)
{
    let mut real_part = DMatrix::<T>::zeros(cp_mat.shape().0,cp_mat.shape().1);
    let mut imag_part = DMatrix::<T>::zeros(cp_mat.shape().0,cp_mat.shape().1);

    let zipped  = std::iter::zip(std::iter::zip(real_part.iter_mut(),imag_part.iter_mut()),cp_mat.iter());

    for ((r,i),c) in zipped
    {
        *r = (*c).re;
        *i = (*c).im;
    } 
    (real_part,imag_part)
}

//get real and imaginary parts but each array is complex
pub fn  cpx_decompose_as_cpx<T:RealField+Num+Zero+Copy>(cp_mat : &DMatrix<Complex<T>>) ->(DMatrix<Complex<T>>,DMatrix<Complex<T>>)
{
    let mut real_part = DMatrix::<Complex<T>>::zeros(cp_mat.shape().0,cp_mat.shape().1);
    let mut imag_part = DMatrix::<Complex<T>>::zeros(cp_mat.shape().0,cp_mat.shape().1);

    let zipped  = std::iter::zip(std::iter::zip(real_part.iter_mut(),imag_part.iter_mut()),cp_mat.iter());

    for ((r,i),c) in zipped
    {
        (*r).re = (*c).re;
        (*i).re = (*c).im;
    } 
    (real_part,imag_part)
}

//check validity of inputs for the place_poles function. less checks than in Python since Rust less dynamic.
pub fn valid_inputs<T:ComplexField+PartialOrd>(mat_a : &DMatrix<T> , poles : &DVector<Complex<T>>, rtol : T, maxiter: u32 , method : &EPoleMethod)
-> (bool, DVector<Complex<T>>)
{
    if !mat_a.is_square() {panic!("A must be square");}
    if poles.shape().0 != mat_a.shape().0 {panic!( "Dimension mismatch : A has {} poles, {} provided",mat_a.shape().0,poles.shape().0);}
    if maxiter < 1 {        panic!("maxiter must be at least equal to 1") ;}
    if rtol > T::one() {panic!("rtol can't be greater than 1");}

    if !poles.iter().all(|x| x.im.is_zero()) && matches!(method,EPoleMethod::KNV0) {panic!("KNV0 Method only works with Real poles")};

    order_complex_poles(poles)
}

pub fn order_complex_poles<T:ComplexField+PartialOrd>(poles : &DVector::<Complex<T>>) -> (bool, DVector<Complex<T>>)
{

    let mut actually_cpx = false;
    let mut ordered_poles = poles.clone();


    ordered_poles.as_mut_slice().sort_by(|a, b| a.re.partial_cmp(&b.re).unwrap());

    for _i in 0..ordered_poles.shape().0
    {

        if !ordered_poles[_i].im.is_zero(){ actually_cpx = true;}

        let mut conj_present = false;
        for _j in 0.. ordered_poles.shape().0
        {
            if ordered_poles[_j] == ordered_poles[_i].conj()
            {
                conj_present = true;
                break;
            }
        }
        if conj_present == false {panic!("complex poles must come in conjugate pairs");}
    }

    (actually_cpx , ordered_poles)
}




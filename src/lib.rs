
extern crate approx; // For the macro relative_eq!
extern crate nalgebra as na;


use num_traits::Zero;
use na::{ zero, Complex, DMatrix, DVector };

#[cfg(test)]
use na::{dmatrix,dvector};

mod secondary_funcs;
mod gipsy_kings;
mod knv0;
mod yangtits;
use secondary_funcs::*;


#[allow(non_snake_case)]
#[allow(dead_code)]
pub struct FullStateFeedBack
{
    pub gain_matrix : DMatrix<f64>,
    pub computed_poles : DVector<Complex<f64>>,
    pub requested_poles : DVector<Complex<f64>>,
    pub X : DMatrix<Complex<f64>>,
    pub rtol : f64,
    pub nb_iter : u32,
    pub err : f64
}

impl FullStateFeedBack
{
    pub fn new(size : usize) ->FullStateFeedBack
    {
        FullStateFeedBack
        {
            gain_matrix : DMatrix::<f64>::zeros(size,size),
            computed_poles : DVector::<Complex<f64>>::zeros(size),
            requested_poles : DVector::<Complex<f64>>::zeros(size),
            X : DMatrix::<Complex<f64>>::zeros(size,size),
            rtol : f64::zero(),
            nb_iter : 0,
            err : 0.0
            
        }
    }
}

impl std::fmt::Display for FullStateFeedBack
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "(gain matrix: {:.3}, computed_poles: {:.3},rtol :{},nb iterations: {})", self.gain_matrix, self.computed_poles,self.rtol,self.nb_iter)
    }
}







pub enum EPoleMethod {YT, KNV0}

#[allow(non_snake_case)]
pub fn place_poles_f64(
    mat_a : &DMatrix<f64> , 
    mat_b : &DMatrix<f64>, 
    poles : &DVector<na::Complex<f64>>, 
    method : EPoleMethod , 
    rtol : f64, 
    maxiter: u32,
    initial_transfer_matrix : Option<DMatrix<Complex<f64>>>,
) -> FullStateFeedBack {
        
        let (_is_cpx,ordered_poles) = valid_inputs(mat_a, poles, rtol, maxiter,&method);
        let cpx_a = real_to_cpx(&mat_a);
        
        
        
        let mut cur_rtol = f64::zero();
        let mut nb_iter : u32 = 0;
        let mut gain_matrix : DMatrix<Complex<f64>>;
        let mut transfer_matrix : DMatrix<Complex<f64>>;
        
        let (u,z) = gipsy_kings::FullQr::new(mat_b.clone()).unpack();
        
        
        let u = real_to_cpx(&u);
        let z = real_to_cpx(&z);
        let rank_b = mat_b.clone().rank(f64::EPSILON);
        let u0 = u.columns(0,rank_b).to_owned();
        
        let u1 = u.columns(rank_b, u.shape().1-rank_b).to_owned();
        let z = z.rows(0,rank_b).to_owned();
        
        //println!("{} {} {}",u0,u1,z);
        
        //If we can use the identity matrix as X the solution is obvious
        if mat_b.shape().0 == rank_b
        {
            let mut diag_poles = DMatrix::<f64>::zeros( mat_a.shape().0,mat_a.shape().1); 
            let mut idx = 0;
            
            while idx < ordered_poles.shape().0
            {
                diag_poles[(idx,idx)] = ordered_poles[idx].re;
                if !ordered_poles[idx].im.is_zero()
                {
                    diag_poles[(idx,idx+1)] = -ordered_poles[idx].im;
                    diag_poles[(idx+1,idx+1)] = ordered_poles[idx].re;
                    diag_poles[(idx+1,idx)] = ordered_poles[idx].im;
                    idx+=1; // skip conjugate
                }
                
                idx+=1;
                
            }
            
            gain_matrix = real_to_cpx(&gipsy_kings::FullQr::new(mat_b.clone()).solve(&(diag_poles-mat_a)).unwrap());
            transfer_matrix = DMatrix::identity(mat_a.shape().0,mat_a.shape().0);
        }
        
        else
        {
            let mut ker_pole : Vec<DMatrix<Complex<f64>>> = Vec::new();
            let mut skip_conjugate = false;
            
            match &initial_transfer_matrix
            {
                Some(X) => {transfer_matrix = X.clone();}
                None => {transfer_matrix = DMatrix::from_element(mat_a.shape().0,1,Complex::zero());}
            }
            
            
            for j in 0..mat_b.shape().0
            {
                if skip_conjugate
                {
                    skip_conjugate = false;
                    continue;
                }
                
                let pole_eye = DMatrix::<Complex<f64>>::from_diagonal_element(mat_b.shape().0,mat_b.shape().0,ordered_poles[j]);
                let pole_space_j = (u1.transpose() * (&cpx_a - &pole_eye) ).transpose();
                
                // println!("----{}------ {} {}",j, pole_space_j,(&cpx_a - &pole_eye));
                
                
                let Q  = gipsy_kings::FullQr::new(pole_space_j.to_owned()).q();
                
                
                
                let ker_pole_j = Q.columns(pole_space_j.shape().1,Q.shape().1-pole_space_j.shape().1);
                
                let mut transfer_matrix_j = DVector::zeros(mat_b.shape().0);
                
                //sum columns of ker_pol_j
                if initial_transfer_matrix.is_none()
                {
                    transfer_matrix_j = ker_pole_j.column_sum();
                    transfer_matrix_j/= Complex::new( ker_pole_j.norm(),zero()) ;
                }
                
                
                
                
                
                if ordered_poles[j].im != zero() //complex pole
                {
                    let mut real = DMatrix::<Complex<f64>>::from_fn(ker_pole_j.shape().0,ker_pole_j.shape().1, |i,j| Complex::new(ker_pole_j[(i,j)].re,f64::zero()));
                    let imag = DMatrix::<Complex<f64>>::from_fn(ker_pole_j.shape().0,ker_pole_j.shape().1, |i,j| Complex::new(ker_pole_j[(i,j)].im,f64::zero()));
                    //transfer_matrix_j = hstack::<Complex<f64>>(vec![real,imag]);
                    real.extend(imag.iter().cloned());
                    
                    if initial_transfer_matrix.is_none()
                    {transfer_matrix_j.copy_from(&real);}
                    
                    ker_pole.extend([ker_pole_j.into_owned(),ker_pole_j.into_owned()]);    
                    
                    //skip next pole as it's conjugate of this one
                    skip_conjugate = true;
                }
                
                else
                {
                    ker_pole.push(ker_pole_j.into_owned());
                }
                
                if initial_transfer_matrix.is_none()
                {
                    if j ==0 { transfer_matrix.copy_from(&transfer_matrix_j); }
                    else
                    {
                        transfer_matrix.extend(transfer_matrix_j.iter().cloned());
                    }
                }
                
            }
            
            let stop;
            if rank_b >1
            {
                match method
                {
                    EPoleMethod::YT => {  (stop,cur_rtol,nb_iter) = yangtits::YT_loop(&ker_pole,&mut transfer_matrix,poles,maxiter,rtol);},
                    EPoleMethod::KNV0 => 
                    {
                        
                        (stop,cur_rtol,nb_iter) = knv0::knv0_loop(&ker_pole,&mut transfer_matrix,mat_b,maxiter,rtol);
                    }
                };
                
                if !stop && rtol> f64::EPSILON
                {
                    println!("Convergence wasn't reached in {} iterations. asked for tolerance {}, got {}",nb_iter,rtol,cur_rtol);
                }
            }
            
            let mut idx = 0;
            
            while idx < ordered_poles.shape().0
            {
                if !(ordered_poles[idx].im ==0.0)
                {
                    let rel = transfer_matrix.column_mut(idx).clone_owned();
                    let img = transfer_matrix.column_mut(idx+1).clone_owned()*Complex::i();
                    
                    transfer_matrix.column_mut(idx).copy_from( &(&rel- &img));
                    transfer_matrix.column_mut(idx+1).copy_from( &(&rel + &img));
                    idx+=1;
                }
                idx+=1;
            }
            
            
            
            let diag_poles = DMatrix::from_diagonal(&ordered_poles);
            let m = transfer_matrix.transpose().qr().solve(&(&diag_poles*transfer_matrix.transpose())).unwrap();
            let m = m.transpose();
            
            gain_matrix = z.qr().solve(&(&u0.transpose()*(&m-&cpx_a))).unwrap();
            
            
            
        }
        
        gain_matrix.scale_mut(-1.0);
        
        let gain_matrix_re = cpx_to_real(&gain_matrix);
        
        
        let mut comp_poles = (mat_a - mat_b*&gain_matrix_re).complex_eigenvalues(); //sort computed poles by 
        comp_poles = order_complex_poles(&comp_poles).1;

        let _err = (&comp_poles-&ordered_poles).norm();
        
        
        FullStateFeedBack
        {
            gain_matrix: gain_matrix_re,
            computed_poles : comp_poles,
            X : transfer_matrix,
            requested_poles : ordered_poles.clone(),
            rtol : cur_rtol,
            nb_iter : nb_iter,
            err : _err
            
        }
        
    }
    

#[test]
fn svdcheck()
{
    let mat = dmatrix![ 
        0.16392308, -1.63535915, -0.84471624,  1.90584585,  1.22223669;
        1.33882279,  0.6401222 ,  0.27215501,  0.33321721, -0.07478135;
        0.32941376, -0.45005802,  0.74089024, -0.59396737,  1.23153441;
        -1.65792534, -0.11407759,  0.54935731, -0.29308824,  1.04328475;
        1.03143172,  0.77403473, -0.42167852, -0.26050874,  0.29815126];
    let eigs = &mat.symmetric_eigen();
    let eigvals : &[f64] = eigs.eigenvalues.as_slice();

    let mut indices = Vec::<usize>::with_capacity(5); for i in 0..5 {indices.push(i);}
    indices.sort_by(|x,y| eigvals[*x].partial_cmp(&eigvals[*y]).unwrap());

    println!("{:?} {:?}",eigvals,indices);

}


#[test]
fn basic_qr_test()
{
    let mat = dmatrix![ 
        0.16392308, -1.63535915, -0.84471624,  1.90584585,  1.22223669;
        1.33882279,  0.6401222 ,  0.27215501,  0.33321721, -0.07478135;
        0.32941376, -0.45005802,  0.74089024, -0.59396737,  1.23153441;
        -1.65792534, -0.11407759,  0.54935731, -0.29308824,  1.04328475;
        1.03143172,  0.77403473, -0.42167852, -0.26050874,  0.29815126];
    
    let mat_b : DMatrix<f64> = dmatrix! [
    0.0,      5.679 ;
    1.136,  1.136 ;
    -0.0,      0.0    ;
    -3.146,  0.0     ];
    

    let mat_big = gipsy_kings::FullQr::new(mat_b.clone());
    let (q,r) = mat_big.unpack();
    println!("{} {} {}",&q,&r, &q*&r);
    
    let (q,r) = mat_b.clone().qr().unpack();
    println!("{} {} ",&q,&r);
}

#[test]
fn scipy_test() 
{
    let mut mat_a = dmatrix![
        1.380,  -0.4,  6.715, -5.676;
        -0.5814, -4.290,   35.0,      0.6750;
        1.18,   4.273,  -6.9,  5.893;
        0.0480,  -4.273,   10.343, -2.104;
        ];
    
    let mat_b = dmatrix![
        90.0,  6.679;
        1.136, 1.6;
        -10.0, 0.0;
        -3.146, 0.0;
        ];
    let p = dvector![
        Complex::new(-0.2,0.0), 
        Complex::new(-0.5,0.0), 
        Complex::new(-5.0566,0.0), 
        Complex::new(-8.6659,0.0),
        ];
    
    let mut fsf0 = place_poles_f64(&mat_a, &mat_b, &p, EPoleMethod::YT, 1e-3, 20,None);
    println!("{}",fsf0);
    //let mut Tmat = fsf0.X;
    // for i in 0..100
    // {
    //     mat_a[(2,2)] = -4.290 + ((i as f64)*0.017/2.0).cos();
    //     mat_a[(0,0)] = 1.380 + (( (i+20) as f64)*0.017/2.0).sin();
        
    //     fsf0 = place_poles_f64(&mat_a, &mat_b, &p, EPoleMethod::KNV0, 1e-3, 20, None);
    //     let fsfupd = place_poles_f64(&mat_a, &mat_b, &p, EPoleMethod::KNV0, 1e-3, 20, Some(Tmat));
    //     Tmat = fsfupd.X;
    //     println!("{} {} {} {}",fsfupd.nb_iter , fsf0.nb_iter, fsfupd.err , fsf0.err);
    // }
}


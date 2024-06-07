use na::{ Complex, ComplexField, DMatrix, MatrixView, RealField , Dyn,Dim};

use super::gipsy_kings;



fn all_close_cst_cpx<T : Copy + ComplexField + RealField,R: Dim,C: Dim>( mat_a : &MatrixView<Complex<T>,R,C>, b: T, eps : T) ->bool
{
    let mut max_dist:T = na::zero();

    for elem in mat_a.iter()
    {
        max_dist = max_dist.max((*elem - b).modulus());
    }
    max_dist<eps
}


/// Single pass of KNV0 algorithm. Implemented as an independent function since YT uses this method for some poles.
pub fn knv0_pass(ker_pole : &Vec<DMatrix<Complex<f64>>>, 
    transfer_matrix : &mut DMatrix::<Complex<f64>>,
    j : usize)
    {
        let transfer_matrix_not_j = transfer_matrix.clone().remove_column(j as usize);
        let mut q = gipsy_kings::FullQr::new(transfer_matrix_not_j).q();

        let mat_ker_pj = ker_pole[j].clone() * ker_pole[j].transpose();

        let mut yj = mat_ker_pj.view((0,0),mat_ker_pj.shape())* q.columns_mut(q.shape().1-1,1);        

        if !all_close_cst_cpx::<f64,Dyn,Dyn>(&yj.view((0,0),yj.shape()),0.0, f64::EPSILON)
        {
            yj.scale_mut(1.0/yj.norm());

            transfer_matrix.column_mut(j as usize).copy_from(&yj);
        }

    }

///Loop over poles and use the KNV0 algorithm to fit the transfer matrix
pub fn knv0_loop(ker_pole : &Vec<DMatrix<Complex<f64>>>, 
    transfer_matrix : &mut DMatrix::<Complex<f64>>,        
    mat_b : &DMatrix<f64>,
    maxiter: u32,
    rtol : f64) -> (bool, f64,u32)
    {
        let mut stop = false;
        let mut nb_try = 0;

        let mut cur_rtol = rtol;
        
        while nb_try < maxiter && !stop
        {
            let det_transfer_matrixb = transfer_matrix.determinant().abs();
            for j in 0..mat_b.shape().0
            {
                knv0_pass(ker_pole, transfer_matrix, j);
            }


            let det_transfer_matrix  = (transfer_matrix.determinant().abs()).max(f64::EPSILON.sqrt()) ;
            cur_rtol  = ( (det_transfer_matrix - det_transfer_matrixb) /   det_transfer_matrix ).abs();
            if cur_rtol < rtol && det_transfer_matrix > f64::EPSILON.sqrt()
            {
                // Convergence test from YT page 21
                stop = true ;
            }

            nb_try += 1;
        }

        (stop, cur_rtol, nb_try)

    }

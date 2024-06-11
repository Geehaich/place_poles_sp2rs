use core::panic;

use num_traits::Zero;


use na::{Complex, ComplexField,  DMatrix, DVector};

use crate::{hstack,vstack, secondary_funcs::*};

use super::gipsy_kings;



#[allow(non_snake_case)]
#[allow(dead_code)]
fn YT_add_orders( update_order:&mut(DVector<isize>,DVector<isize>), to_add : &DVector<isize>,add_ones:bool) ->()
{
    if add_ones
    {
        update_order.0.extend([1]);
        update_order.1.extend([1]);
    }
    update_order.0.extend(to_add.iter().cloned());
    update_order.1.extend(to_add.add_scalar(1).iter().cloned());
}

#[allow(non_snake_case)]
#[allow(dead_code)]
pub fn YT_loop(ker_pole : &Vec<DMatrix<Complex<f64>>>, transfer_matrix : &mut DMatrix::<Complex<f64>>, poles : &DVector<Complex<f64>>,
    maxiter: u32, rtol : f64) -> (bool, f64,u32)
    {

        let nb_real :isize = poles.iter().filter(|&x| x.imaginary().is_zero()).count() as isize;
        let hnb :isize = nb_real/2;
        
        let mut update_order = (DVector::<isize>::from_vec(vec![]),DVector::<isize>::from_vec(vec![]));
        if nb_real >0 //update the biggest real pole with the smallest one
        {
            update_order.0.extend([nb_real]);
            update_order.1.extend([1]);
        }
        
        //1.a
        
        let r_comp = (nb_real+1..poles.len() as isize +1).step_by(2).into_iter();
        let r_comp = DVector::from_iterator(r_comp.size_hint().0,r_comp);
        
        let r_p = (1..hnb+nb_real%2).into_iter();
        let r_p   = DVector::<isize>::from_iterator(r_p.size_hint().0 as usize,r_p);

        
        
        update_order.0.extend((2*&r_p).iter().cloned());
        update_order.1.extend((2*&r_p).add_scalar(1).iter().cloned());
        
        //1.b
        
        update_order.0.extend(r_comp.iter().cloned());
        update_order.1.extend(r_comp.add_scalar(1).iter().cloned());
        
        
        // step 1.c
        let r_p = (1..hnb+1).into_iter();
        let r_p = DVector::from_iterator(r_p.size_hint().0,r_p);
        
        update_order.0.extend((2*&r_p).add_scalar(-1).iter().cloned());
        update_order.1.extend((2*&r_p).iter().cloned());
        // step 1.d
        
        YT_add_orders(&mut update_order, &r_comp, hnb == 0 && poles[0].imaginary().is_zero());
        // step 2.a
        
        let r_j = 2 .. hnb+nb_real % 2;
        
        for j in r_j
        {
            for i in 1..hnb+1
            {
                update_order.0.extend([1]);
                update_order.1.extend([i+j]);
            }
        }
        // step 2.b
        
        
        YT_add_orders(&mut update_order, &r_comp, hnb == 0 && poles[0].imaginary().is_zero());
        
        
        // step 2.c
        
        
        for j in 2..hnb+nb_real%2
        {
            for i in hnb+1..nb_real+1
            {
                let mut idx_1 = i+j;
                if idx_1 >nb_real
                {
                    idx_1 = i+j-nb_real;
                }
                update_order.0.extend([i]);
                update_order.1.extend([idx_1]);
            }
        }
        
        // step 2.d
        
        YT_add_orders(&mut update_order, &r_comp, hnb == 0 && poles[0].imaginary().is_zero());
        
        
        // step 3.a
        for i in 1.. hnb+1
        {
            update_order.0.extend([i]);
            update_order.1.extend([i+hnb]);
        }
        // step 3.b
        
        YT_add_orders(&mut update_order, &r_comp, hnb == 0 && poles[0].imaginary().is_zero());
        
        let mut stop = false;
        let mut nb_try : u32 = 0;
        let mut cur_rtol = f64::NAN;

        update_order.0.add_scalar_mut(-1);
        update_order.1.add_scalar_mut(-1);
        


        while nb_try< maxiter && !stop
        {
            
            let det_transfermatrixb = transfer_matrix.determinant().abs();      
            for order in 0..update_order.0.len()
            {
                let (i,j) = (update_order.0[order] as usize, update_order.1[order] as usize);
                
                if i==j
                {
                    assert!(i==0,"i!=0 for KNV call in YT");
                    assert!(poles[i as usize].imaginary().is_zero(), "calling KNV on a complex pole");
                    crate::knv0::knv0_pass(&ker_pole,transfer_matrix,j as usize);
                }
                
                else
                {
                    let transfer_matrix_not_i_j =  transfer_matrix.clone().remove_column(i.max(j)).remove_column(i.min(j));
                    let Q = gipsy_kings::FullQr::new(transfer_matrix_not_i_j.clone()).q();

                    //panic!("{} {}",&transfer_matrix_not_i_j.clone(),Q);
                    
                    if poles[i as usize].im == 0.0
                    {
                        assert!(poles[j as usize].im==0.0, "mixing real and complex in YT_real");
                        _yt_real(&ker_pole,Q,transfer_matrix,i ,j);
                    }
                    else 
                    {
                        assert!(poles[j as usize].im!=0.0, "mixing real and complex in YT_real");
                        _yt_complex(&ker_pole, Q, transfer_matrix, i as usize, j as usize);
                        
                    }
                    
                    
                }
                
                let det_transfer_matrix = transfer_matrix.determinant().abs().max(f64::EPSILON.sqrt());
                cur_rtol =(det_transfer_matrix-det_transfermatrixb).abs()/det_transfermatrixb;
                
                if cur_rtol < rtol && det_transfer_matrix > f64::EPSILON.sqrt()
                {
                    stop = true;
                }
                
            }
            nb_try+=1;
        }
        
        
        (stop,cur_rtol,nb_try)
        
        
    }
    
    //    Applies algorithm from YT section 6.1 page 19 related to real pairs
    fn _yt_real(ker_pole : &Vec<DMatrix<Complex<f64>>>, q_mat : DMatrix<Complex<f64>> , transfer_matrix : &mut DMatrix::<Complex<f64>>, i:usize, j:usize)
    {
        let u = q_mat.column(q_mat.shape().1-2);
        let v = q_mat.column(q_mat.shape().1-1);
        
        let (rows, _) = transfer_matrix.shape();

        let m = &ker_pole[i].transpose()* ( u*v.transpose() - v*u.transpose() ) * &ker_pole[j];
        let svd_m = m.clone().svd(true,true);
        let (um,vm) = (svd_m.u.unwrap().transpose(),svd_m.v_t.unwrap().transpose());
        
        let mu1 = um.column(0);
        let mu2 = um.column(1);
        
        let nu1 = vm.column(0);
        let nu2 = vm.column(1);
        
        let mut ker_pole_mu_nu : DMatrix<Complex<f64>>;
        
        let (ci,cj) = (transfer_matrix.column(i),transfer_matrix.column(j));
        let transfer_matrix_j_mo_transfer_matrix_j = crate::vstack![ci,cj];


        if (svd_m.singular_values[0]-svd_m.singular_values[1]).abs() > 1e-8
        {
            let ker_pole_imo_mu1 = &ker_pole[i]*mu1;
            let ker_pole_i_nu1 =  &ker_pole[j]*nu1;
            ker_pole_mu_nu  = crate::vstack![&ker_pole_imo_mu1,&ker_pole_i_nu1];
        }
        
        else 
        {
            let _z = DMatrix::<Complex<f64>>::zeros(ker_pole[i].nrows(),ker_pole[i].ncols());
            let ker_istack = crate::hstack![ker_pole[i],_z];
            let _z = DMatrix::<Complex<f64>>::zeros(ker_pole[j].nrows(),ker_pole[j].ncols());
            let ker_jstack = crate::hstack![ker_pole[j],_z];
            let ker_pole_ij = crate::vstack![&ker_istack,&ker_jstack];
            
            let mumu = crate::hstack![mu1,mu2];
            let nunu = crate::hstack![nu1,nu2];
            let mu_nu_matrix = crate::vstack![&mumu,&nunu];
            
            ker_pole_mu_nu = ker_pole_ij*mu_nu_matrix;
        }
        

        let mut transfer_matrix_ij = (&ker_pole_mu_nu*&ker_pole_mu_nu.transpose())*&transfer_matrix_j_mo_transfer_matrix_j;

        
        if is_close_cpxcst(&transfer_matrix_ij, Complex::zero(), 1e-7) == false
        {
            transfer_matrix_ij.scale_mut(2.0.sqrt()/transfer_matrix_ij.norm());
            transfer_matrix.column_mut(i).copy_from(&transfer_matrix_ij.view((0,0),(rows,1)));
            transfer_matrix.column_mut(j).copy_from(&transfer_matrix_ij.view((rows,0),(rows,1)));
            
        }
        
        else
        {
            transfer_matrix.column_mut(i).copy_from(&ker_pole_mu_nu.view((0,0),(rows,1)));
            transfer_matrix.column_mut(j).copy_from(&ker_pole_mu_nu.view((rows,0),(rows,1)));
        }
        
    }
    
    fn _yt_complex(ker_pole : &Vec<DMatrix<Complex<f64>>>, q_mat : DMatrix<Complex<f64>> , transfer_matrix : &mut DMatrix::<Complex<f64>>, i:usize, j:usize)
    {
        let ur = q_mat.column(q_mat.shape().1-2).scale(2.0.sqrt());
        let ui = q_mat.column(q_mat.shape().1-1).scale(2.0.sqrt());
        let mut u = DMatrix::<Complex<f64>>::zeros(q_mat.shape().0,1);
        for i in 0..q_mat.shape().0
        { u[i] = ur[i] + ui[i]*Complex::i();}
        
        
        let ker_pole_ij = &ker_pole[i];
        let m = ker_pole_ij.transpose().conjugate() * (&u*&u.conjugate().transpose() - &u.conjugate()*&u.transpose())*ker_pole_ij;
        
        let s_eig_m = m.symmetric_eigen();
        let (e_val, e_vec) = (s_eig_m.eigenvalues, s_eig_m.eigenvectors);
        
        let n_eigs = e_val.nrows();
        let mut e_val_idx = Vec::<usize>::with_capacity(n_eigs); //ordering which would sort e_val by modulus
        for i in 0..e_val.len() { e_val_idx.push(i);}
        e_val_idx.sort_by(|i,j| e_val[*i].modulus_squared().partial_cmp(&e_val[*j].modulus_squared()).unwrap());
        
        let mu1 = e_vec.columns(e_val_idx[n_eigs-1],1);
        let mu2 = e_vec.columns(e_val_idx[n_eigs-2],1);
        
        let mut transfer_matrix_j_mo_transfer_matrix_j = DMatrix::<Complex<f64>>::zeros(transfer_matrix.nrows(),1);
        for row in 0..transfer_matrix.nrows()
        {
            transfer_matrix_j_mo_transfer_matrix_j[row] = transfer_matrix[(row,i)]+transfer_matrix[(row,j)]*Complex::i();
        }

        
        let ker_pole_mu : DMatrix<Complex<f64>>;
        let mut transfer_matrix_i_j : DMatrix<Complex<f64>> = DMatrix::zeros(0,0);
        
        if (e_val[e_val_idx[n_eigs-1]].abs() - e_val[e_val_idx[n_eigs-2]].abs()).abs() >= 1e-7
        {
            ker_pole_mu = ker_pole_ij*mu1;
        }
        
        else
        {    
            let mu1_mu2_matrix = hstack![mu1,mu2];

            ker_pole_mu = ker_pole_ij*mu1_mu2_matrix;
            transfer_matrix_i_j = (&ker_pole_mu*&ker_pole_mu.transpose().conjugate())*transfer_matrix_j_mo_transfer_matrix_j;
        }
        
        if is_close_cpxcst(&transfer_matrix_i_j, Complex::zero(), 1e-7) == false
        {
            transfer_matrix_i_j.scale_mut(2.0.sqrt()/transfer_matrix_i_j.norm());
            for row in 0..transfer_matrix.nrows()
            {
                transfer_matrix[(row,i)].re = transfer_matrix_i_j[(row,0)].re;
                transfer_matrix[(row,i)].im = 0.0;
                transfer_matrix[(row,j)].im = 0.0;
                transfer_matrix[(row,j)].re = transfer_matrix_i_j[(row,0)].im;
            }
            
        }        
        else
        {
            for row in 0..transfer_matrix.nrows()
            {
                transfer_matrix[(row,i)].re = ker_pole_mu[(row,0)].re;
                transfer_matrix[(row,i)].im = 0.0;
                transfer_matrix[(row,j)].im = 0.0;
                transfer_matrix[(row,j)].re = ker_pole_mu[(row,0)].im;
            }
        }
        
        
    }

use num_traits::Zero;


use na::{Complex, ComplexField,  DMatrix, DVector};

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
fn YT_loop(ker_pole : Vec<DMatrix<Complex<f64>>>, transfer_matrix : &mut DMatrix::<Complex<f64>>, poles : &DVector<f64>,
    B : &DMatrix<f64>,rtol : f32, maxiter: u8) -> (bool, f32, u8)
    {
        let mut nb_real :isize = poles.iter().filter(|&x| x.imaginary().is_zero()).count() as isize;
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
        let r_p   = DVector::<isize>::from_iterator(r_p.size_hint().0,r_p);

        update_order.1.extend((2*&r_p).iter().cloned());
        update_order.0.extend((2*&r_p).add_scalar(1).iter().cloned());

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
        let mut nb_try = 0;
        
        while nb_try< maxiter && !stop
        {
         let det_transfermatrix = transfer_matrix.determinant().abs();
        }

        for order in 0..update_order.0.len()
        {
            let (i,j) = (update_order.0[order], update_order.1[order]);
            
            if i==j
            {
                assert!(i==0,"i!=0 for KNV call in YT");
                assert!(poles[i as usize].imaginary().is_zero(), "calling KNV on a complex pole");
                crate::knv0::knv0_pass(&ker_pole,transfer_matrix,j as usize);
            }

            else
            {
                let transfer_matrix_not_i_j =  transfer_matrix.clone().remove_column(j as usize).remove_row(i as usize);
                let (Q,_) = gipsy_kings::full_qr_c64(&transfer_matrix_not_i_j);

                if poles[i].im == 0.0
                {
                    assert!(poles[j].im==0.0, "mixing real and complex in YT_real");
                    _YT_real(ker_pole,Q,transfer_matrix,i,j);
                }
                else 
                {
                    assert!(poles[i].im==0.0, "mixing real and complex in YT_real");
                    _YT_complex(ker_pole, Q, transfer_matrix, i, j);
                }

                
            }
            
        }
        

        (stop,0.0,nb_try)


    }

//    Applies algorithm from YT section 6.1 page 19 related to real pairs
    fn _YT_real(ker_pole : Vec<DMatrix<Complex<64>>>, q_mat : DMatrix<Complex<F64>> , transfer_matrix : &mut DMatrix::<Complex<f64>>, i:usize, j:usize)
    {
        let u = Q.co
    }




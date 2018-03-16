program fort_subroutines

  !Empty program, only want individual subroutines to be used by python
  !
  !This file contains 4 subroutines:
  !1. calc_dr_dv which calculates dv and dr for every possible pair of stars
  !It works whether or not there are velocity errors via an optional argument
  !
  !2 and 3 are sort_no_errors and sort_errors
  !These both sort the pairs into dr bins, calculate the average dv in each bin,
  !and the error on that.
  !The first is to be used when there are errors on the velocities, and the latter if not.
  !The maths of how the final errors are calculated is different enough in these
  !two cases to warrent these being two different subroutines.
  !
  !4. calc_av_dv calculates average dv for all pairs of stars

end program fort_subroutines




!This subroutine calculates dr and dv for every pair of stars.
!It also returns the largest distance between a pair in this cluster, max_dr
!rdim and vdim are inputs, they are the number of dimensions the r and v data is in
!It is necessary to supply this information because fortran is a compliled language and so
!needs to know the size of the arrays before it makes them.
subroutine calc_dr_dv(n_stars, n_pairs, bin_width, rdim, r, vdim, v, error_flag, dv_default, verr_in, dr_dv, max_dr)

  !Delare variables.
  integer, intent(in) :: n_stars, n_pairs, rdim, vdim
  integer :: star_a, star_b, count, dim
  real, intent(in) :: bin_width
  real, intent(in) :: v(vdim, n_stars)
  real, intent(in) :: r(rdim, n_stars)
  logical, intent(in) :: error_flag
  logical, intent(in) :: dv_default
  real, optional :: verr_in(vdim, n_stars)
  real :: verr(vdim, n_stars)
  real :: dr, dv
  real, intent(out) :: dr_dv(n_pairs, 3)
  real, intent(out) :: max_dr


  !This how fortran does optional arguments. If velocity errors are put in (verr_in)
  !then assign them to the variable verr, otherwise verr just gets assigned zero.
  if(present(verr_in))then
     verr = verr_in
  else
     verr = 0.
  endif


  !State which variables are inputs and/or inputs for f2py which
  !is the package used to make this subroutine compatible with python
  !f2py intent(in) n_stars, n_pairs, bin_width, rdim, r, vdim, v, error_flag, dv_default, verr_in
  !f2py intent(out) dr_dv, max_dr
  !f2py depend(n_stars) r, v, verr, verr_in
  !f2py depend(n_pairs) dr_dv
  !f2py depend(rdim) r
  !f2py depend(vdim) v, verr, verr_in

  !Set some of the currently empty variables that will later be filled to zero
  !A good practice thing
  dr_dv = 0.
  count = 1

  !Loop through each pair and figure out its dr and dv
  do star_a = 2, n_stars
     do star_b = 1, star_a - 1

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        !Get the dv and dr of the pair
        !Do it with a loop through each dimension because the number of dimensions is not known ahead of time.
        !i.e. get dx^2, dy^2 and dz^2 if you have the data, then sqrt to get dr.
        dr = 0.
        do dim = 1, rdim
           dr  = dr + ((r(dim, star_a) - r(dim, star_b))**2) 
        end do
        dr = sqrt(dr)

        !There are two ways to define dv
        !Use the chosen one as indicated by the dv_default flag
        !The default way calculates dv in the same way as dr is
        !calculated: as the magnitude of the difference. In the same
        !Way it's impossible to get a negative distance between stars, by
        !this definition dv is always positive.
       
        if (dv_default) then
           
           dv = 0.
           do dim = 1, vdim
              dv = dv + ((v(dim, star_a) - v(dim, star_b))**2) 
           end do
           dv = sqrt(dv)

        !The alternative definition of dv is the rate at which the
        !distance between the stars is changing, so dr/dt. If the
        !distance between the stars is decreasing (so moving towards
        !one another) this value is negative. If they are moving
        !further apart the distance between the strs increases, so
        !the velue returned is positive

        !In 2D that'd be ((x_a - x_b)*(vx_a - vx_b)) + ((y_a - y_b)*(vy_a - vy_b)) / sqrt(dx^2 + dy^2)
        !So the denominator is just dr (which has already been calculated).
        !Just calculated the numerato and divide by dr.
        else
           
           dv = 0.
           do dim = 1, vdim
              dv = dv + (r(dim, star_a) - r(dim, star_b))*(v(dim, star_a) - v(dim, star_b))
           end do
           dv = dv / dr

        end if


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        
        !If there's errors then need to calculate the error on dv
        if (error_flag) then

           !But first a quick check

           !It is essentially impossible that any two stars will have the *exact* same velocity.
           !However it is impossible to measure an infinite number of decimal places.
           !Therefore due to instrument resulution it is possible two
           !stars may be measured to have identical velocities, so dv = 0.
           !This would result in dividing by zero in the error calculation, which is the next step.
           
           !To avoid this if dv is calculated to be zero we replace it with a value that
           !is just small, 0.1. (assumed km/s).

           !From a mathematical perspective this is fudging things somewhat
           !However it is unlikely two stars will be measured to have identical velocities
           !and even in the case this does occur it would impact one pair of stars
           !among likely 1000s of others.
           
           if (dv .eq. 0.) then
              dv = 0.1
           end if

           !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

           !How the error is calculated depends on how dv is defined
           if (dv_default) then
        
              !Get the error on dv if dv is calculated with the default definition

              !You can work through the maths, but in 2D it's:
              !sqrt(  ((v1x - v2x)^2 * (err_v1x^2 + err_v2x^2)) + ((v1y - v2y)^2 * (err_v1y^2 + err_v2y^2))  )
              !Divided by
              ! sqrt(  (v1x - v2x)^2 + (v1y - v2y)^2   )
              !i.e denominator is dv

              !Get numerator
              dv_err = 0.
              do dim = 1, vdim    
                 dv_err  = dv_err + (((v(dim, star_a) - v(dim, star_b))**2) * ((verr(dim, star_a)**2) + (verr(dim, star_b)**2)))
              end do
              dv_err = sqrt(dv_err)
           
              !Divide by denominator to get the error on dv, then add it to the output array.
              dv_err  = dv_err / dv
              dr_dv(count, 3) = dv_err

           else
              
              !Get the error on dv if dv is calculated with the alternative definition.

              !You can work through the maths, but in 2D it's:
              !sqrt(  (  (x_a - x_b) * sqrt(err_vx_a^2 + err_vx_b^2) )^2 + the same for y   )
              !Divided by
              ! sqrt(  (x_a - x_b)^2 + (y_a - y_b)^2   )
              !i.e denominator is dr

              !Get numerator
              dv_err = 0.
              do dim = 1, vdim    
                 dv_err  = dv_err + ( ((r(dim, star_a) - r(dim, star_b))**2) * ((verr(dim, star_a)**2) + (verr(dim, star_b)**2)))
              end do
              dv_err = sqrt(dv_err)
           
              !Divide by denominator to get the error on dv, then add it to the output array.
              dv_err  = dv_err / dr
              dr_dv(count, 3) = dv_err

           end if

        end if
 
        !Add the data for this pair to the array
        dr_dv(count, 1) = dr
        dr_dv(count, 2) = dv
  
        count = count + 1

     end do

  end do

  !Get the max dist between a pair of stars in this cluster
  !Add a little padding to the end so the last bin with defintly include
  !this dr.
  max_dr = maxval(dr_dv(:,1)) + (2.*bin_width)
  
end subroutine calc_dr_dv




!This subroutine is for the case where there are no observational errors.
!It sorts the pairs of stars into dr bins
!It calculates the mean dv in each bin and returns it in the array mean_dv
!It also calculates out the error on the mean dv, i.e. the standard error, and returns it in the array error.
!It notes how many pairs are sorted into each bin, and returns that in the array n_in_bins
!It also notes how many times each star appears in each bin, and returns that in the array count_stars_bins
subroutine sort_no_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width, mean_dv, error, n_in_bins, count_stars_bins)

  !Declare variables.
  integer, intent(in) :: n_stars, n_pairs, n_bins
  double precision, intent(in) :: dr_dv(n_pairs, 3)
  double precision, intent(in) :: bin_width
  double precision, intent(out) :: mean_dv(n_bins) 
  double precision, intent(out) :: n_in_bins(n_bins)
  integer, intent(out) :: count_stars_bins(n_bins, n_stars)
  !In otder to calculate the error that will be returned (the standard error) first calculate std deviation
  !of the dvs in each bin.
  !Standard deviation demands the sum of (value - mean value)^2 which implies I'd have to sort twice;
  !once to get the mean, then again to get the (value - mean value)^2. However if I 
  !expand square inside the sum and work through the maths it's possible to do it all in one go.
  !That's what I've done here, because sorting twice is ugly and a waste of computer time. The maths is a
  !little more fidly as a result, so as I sort into bins I total the dvs and dv^2s sorted into each bin
  !Those are held in the arrays declared below, and are used to get the standard deviation.
  double precision :: sum_dv(n_bins) !List that holds the total dv in each bin for all the pairs it contains
  double precision :: sum_dv_squ(n_bins) !List that holds the summed (dv^2) in each bin for all the pairs it contains
  double precision :: std_dev(n_bins)!List to hold the std dev of dvs in each bin
  double precision, intent(out) :: error(n_bins) !Standard error due to scatter on dv and number in the bin
                                                 !Equals std dev of dvs in the bin / sqrt(n of pairs in the bin)
  integer :: pair, pair_bin, bin, star_a, star_b
  double precision :: dr

  !State which variables are inputs and/or inputs for f2py which
  !is the package used to make this subroutine compatible with python
  !f2py intent(in) n_stars, n_pairs, n_bins, dr_dv, bin_width
  !f2py intent(out) mean_dv, n_in_bins, error, count_stars_bins
  !f2py depend(n_pairs) dr_dv
  !f2py depend(n_bins) mean_dv, n_in_bins, sum_dv_sq, std_dev, error, count_stars_bins
  !f2py depend(n_stars) count_stars_bins

  !Set to 0 to start with for safety
  mean_dv = 0.
  n_in_bins = 0.
  sum_dv = 0.
  sum_dv_squ = 0.
  std_dev = 0.
  error = 0.
  count_stars_bins = 0
  star_a = 1
  star_b = 2


  !Go through every pair and figure out which bins it belongs to
  !Then add that pair's dv and dv^2 to those bins
  do pair = 1, n_pairs

     !Get what the dr is for this pair
     dr = dr_dv(pair,1)

     !Figure out the bin this dr is in
     pair_bin = ceiling(dr/bin_width)

     !Add the dv for this pair to the correct bins
     n_in_bins(pair_bin) =  n_in_bins(pair_bin) + 1
     sum_dv(pair_bin) =  sum_dv(pair_bin) + dr_dv(pair, 2)
     sum_dv_squ(pair_bin) =  sum_dv_squ(pair_bin) + (dr_dv(pair, 2)**2)

     !Figure out which two stars make up this pair
     star_b = star_b + 1
     if (star_b .ge. star_a) then
        star_b = 1
        star_a = star_a + 1
     end if

     !Note in the count_stars_bins array that these two stars appeared in this dr bin
     count_stars_bins(pair_bin, star_a) = count_stars_bins(pair_bin, star_a) + 1
     count_stars_bins(pair_bin, star_b) = count_stars_bins(pair_bin, star_b) + 1

  end do


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !Now the pairs have been sorted into dr bins use the values recorded to get
  !the mean, standard deviation, and from that standard error in each bin

  !Go through each bin and work them out.
  do bin = 1, n_bins

     !First get the mean dv in the bin
     mean_dv(bin) = sum_dv(bin)/n_in_bins(bin)

     !Next get the standard deviation of the dvs in each bin
     !This is
     std_dev(bin) = sqrt((sum_dv_squ(bin) + (n_in_bins(bin)*(mean_dv(bin)**2)) - (2*mean_dv(bin)*sum_dv(bin))) / n_in_bins(bin))

     !Finally get the standard error in the bin
     error(bin) = std_dev(bin) / sqrt(n_in_bins(bin))
         
  end do


end subroutine sort_no_errors



!This subroutine is for the case where there are observational errors.
!It sorts the pairs of stars into dr bins
!It calculates the weighted mean dv in each bin and returns it in the array mean_dv
!It also calculates out the combined error on the mean dv, and returns it in the array error.
!Combined error is the standard error (std_er) added in quadrature with the error on the weighted mean (err_dv)
!The subroutine also notes how many pairs are sorted into each bin, and returns that in the array n_in_bins
!It also notes how many times each star appears in each bin, and returns that in the array count_stars_bins
subroutine sort_errors(n_stars, n_pairs, n_bins, dr_dv, bin_width, mean_dv, error, n_in_bins, count_stars_bins)

  !Declare variables.
  integer, intent(in) :: n_stars, n_pairs, n_bins
  double precision, intent(in) :: dr_dv(n_pairs, 3)
  double precision, intent(in) :: bin_width
  double precision, intent(out) :: mean_dv(n_bins) ! Weighted mean equals c / b
  double precision, intent(out) :: error(n_bins)
  double precision, intent(out) :: n_in_bins(n_bins)
  integer, intent(out) :: count_stars_bins(n_bins, n_stars)
  !Standard deviation demands the sum of (value - mean value)^2 which implies I'd have to sort twice;
  !once to get the mean, then again to get the (value - mean value). However if I 
  !expand square inside the sum and work through the maths it's possible to do it all in one go.
  !That's what I've done here, because sorting twice is ugly and a waste of computer time.
  !The maths is a little more fidly as a result.
  !To prevent equations from getting absurdly long the terms in the error calculations
  !are broken up into a, b, and c.
  double precision :: a(n_bins) ! a = sum (dv^2) / (error_on_dv^2)   (sum is over all pairs in the bin)
  double precision :: b(n_bins) ! b = sum 1 / (error_on_dv^2)        (sum is over all pairs in the bin)
  double precision :: c(n_bins) ! c = sum ( dv / (error_on_dv^2) )   (sum is over all pairs in the bin)
  double precision :: err_dv(n_bins) !Error on weighted mean due to uncertainties on dvs
                                     !Equals sqrt(1/b)
  double precision :: std_er(n_bins) !Standard error due to scatter on dv and number in the bin
                                     !Equals (1/b)*sqrt( (ab - (c^2)) / n of pairs in the bin )
  integer :: pair, pair_bin, bin, star_a, star_b
  double precision :: dr

  !State which variables are inputs and/or inputs for f2py which
  !is the package used to make this subroutine compatible with python
  !f2py intent(in) n_stars, n_pairs, n_bins, dr_dv, bin_width
  !f2py intent(out) mean_dv, error, n_in_bins, count_stars_bins
  !f2py depend(n_pairs) dr_dv
  !f2py depend(n_bins) mean_dv, error, n_in_bins, a, b, c, err_dv, std_er, count_stars_bins
  !f2py depend(n_stars) count_stars_bins
  
  !Set to 0 to start with for safety
  mean_dv = 0.
  error = 0.
  n_in_bins = 0.
  a = 0.
  b = 0.
  c = 0.
  err_dv = 0.
  std_er = 0.
  count_stars_bins = 0
  star_a = 1
  star_b = 2

  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

  !Go through every pair and figure out which bins it belongs to
  !Then add that pair's dv to those bins
  do pair = 1, n_pairs

     !Get what the dr is for this pair
     dr = dr_dv(pair,1)

     !Figure out the dr bin this pair is in
     pair_bin = ceiling(dr/bin_width)

     !Add the relevent values for this pair to the correct bins
     n_in_bins(pair_bin) =  n_in_bins(pair_bin) + 1
     a(pair_bin) = a(pair_bin) + ((dr_dv(pair, 2) / dr_dv(pair, 3))**2)
     b(pair_bin) = b(pair_bin) + (1. / dr_dv(pair, 3)**2)
     c(pair_bin) = c(pair_bin) + (dr_dv(pair, 2) / (dr_dv(pair, 3)**2))

     !Figure out which two stars were compared for this pair
     star_b = star_b + 1
     if (star_b .ge. star_a) then
        star_b = 1
        star_a = star_a + 1
     end if

     !Note in the count_stars_bin array that these two stars appeared in this dr bin
     count_stars_bins(pair_bin, star_a) = count_stars_bins(pair_bin, star_a) + 1
     count_stars_bins(pair_bin, star_b) = count_stars_bins(pair_bin, star_b) + 1      


  end do
  
     
  
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  !Use a, b, c and the n_in_bins to calculate the weighted mean, the error on it, and the standard error.
  !Finaly combine the errors.
  !Do this in each bin
  !The weighted mean is c/b
  !The error on it (err_dv) is sqrt(1/b)
  !The standard error is  (1/b)*sqrt( (ab - (c^2)) / n_in_this_bin )
  !The combined_errors are sqrt(standard error^2 + error on weighted mean^2 + error on no struct average^2)
  !That's what gets put in the errors array to be returned.
  do bin = 1, n_bins

     mean_dv(bin) =  (c(bin)/b(bin))
     err_dv(bin) = sqrt(1./b(bin))
     std_er(bin) = (1./b(bin)) * sqrt( ((a(bin)*b(bin)) - (c(bin)**2)) / n_in_bins(bin) )
     error(bin) = sqrt( (err_dv(bin)**2) + (std_er(bin)**2) )
   
  end do


end subroutine sort_errors


!This subroutiuine gets the velocity diffrence for every possible
!pair of stars and averages them
subroutine calc_av_dv(n_stars, v, vdim, av_dv)

  implicit none

  integer, intent(in) :: n_stars, vdim
  integer :: i, j, k
  double precision, intent(in) :: v(1:vdim, n_stars)
  double precision, intent(out) :: av_dv
  double precision :: dv, total_dv, n_pairs

  !f2py intent(in) n_stars, vdim, v
  !f2py intent(out) av_dv
  !f2py depend(n_stars) v
  !f2py depend(vdim) v

  !Initialise total velocity difference as 0
  total_dv = 0.

  !Go through every pair of stars
  do i = 1,n_stars

     do j = 1,i

        !Find the v difference between this pair of stars
        !Have to go through all dimensions to do that
        dv = 0.
        do k = 1,vdim
           dv =  dv + (( v(k,i) - v(k,j) )**2)
        end do
        dv = sqrt(dv)

        !Add it to the total dv
        total_dv = total_dv + dv

     end do

  end do

  !For n stars the number of possible pairs is
  n_pairs = ((n_stars*n_stars) - n_stars)/2.

  !Change the total dv to average dv 
  av_dv = total_dv / n_pairs 

  return

end subroutine calc_av_dv



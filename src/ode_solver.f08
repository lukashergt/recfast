module ode_solver
    use precision, only : dp
    implicit none

    contains

!##############################################################################
      subroutine dverk (n, fcn, x, y, xend, tol, ind, c, nw, w)
        implicit none

        integer n, ind, nw, k
        real(dp) :: x, y(n), xend, tol, c(*), w(nw,9), temp

        external fcn

!       ******************************************************************
!       * begin initialization, parameter checking, interrupt re-entries *
!       ******************************************************************

!  ......abort if ind out of range 1 to 6
        if (ind < 1 .or. ind > 6) go to 500

        ! cases - initial entry, normal re-entry, interrupt re-entries
        !go to (5, 5, 45, 1111, 2222, 2222), ind
        if (ind == 3) goto 45
        if (ind == 4) goto 1111
        if (ind == 5 .or. ind == 6) goto 2222

        ! case 1 - initial entry (ind == 1 or 2)
!  .........abort if n > nw or tol <= 0
        if (n > nw .or. tol <= 0._dp) go to 500
        if (ind == 1) then
            ! initial entry without options (ind == 1)
            ! set c(1) to c(9) equal to 0
            do k = 1, 9
                c(k) = 0._dp
            end do
        else if (ind == 2) then
            ! initial entry with options (ind == 2)
            ! make c(1) to c(9) non-negative
            do k = 1, 9
                c(k) = abs(c(k))
            end do
            ! make floor values non-negative if they are to be used
            if (abs(c(1) - 4._dp) < tiny(1._dp) .or. abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 4 .or. c(1) == 5
                do k = 1, n
                    c(k + 30) = abs(c(k + 30))
                end do
            end if
        end if
        ! initialize rreb, dwarf, prev xend, flag, counts
        c(10) = 2._dp**(-56)
        c(11) = 1.e-35_dp
        ! set previous xend initially to initial value of x
        c(20) = x
        do k = 21, 24
            c(k) = 0._dp
        end do
        go to 50
        ! case 2 - normal re-entry (ind == 3)
!  .........abort if xend reached, and either x changed or xend not
        ! if c(21) /= 0  .and.  ( x /= c(20) .or. xend == c(20) )
   45   if (abs(c(21)) >= tiny(1._dp) .and. (abs(x - c(20)) >= tiny(1._dp) .or. abs(xend - c(20)) < tiny(1._dp))) go to 500
!           re-initialize flag
            c(21) = 0._dp
            ! case 3 - re-entry following an interrupt (ind == 4 to 6)
            ! transfer control to the appropriate re-entry point..........
            ! this has already been handled by the computed go to        .
            ! end cases                                                     v
   50   continue

!       end initialization, etc.

!       ******************************************************************
!       * loop through the following 4 stages, once for each trial  step *
!       * until the occurrence of one of the following                   *
!       *    (a) the normal return (with ind == 3) on reaching xend in *
!       *        stage 4                                                 *
!       *    (b) an error return (with ind < 0) in stage 1 or stage 4 *
!       *    (c) an interrupt return (with ind  ==  4,  5  or  6),  if *
!       *        requested, in stage 1 or stage 4                        *
!       ******************************************************************

99999   continue

!       ***************************************************************
!       * stage 1 - prepare - do calculations of  hmin,  hmax,  etc., *
!       * and some parameter  checking,  and  end  up  with  suitable *
!       * values of hmag, xtrial and htrial in preparation for taking *
!       * an integration step.                                        *
!       ***************************************************************

!***********error return (with ind=-1) if no of fcn evals too great
        if (abs(c(7)) >= tiny(1._dp) .and. c(24) >= c(7)) then  ! if c(7) /= 0 .and. c(24) >= c(7)
            ind = -1
            return
        end if

        ! calculate slope (adding 1 to no of fcn evals) if ind /= 6
        if (ind /= 6) then
            call fcn(n, x, y, w(1,1))
            c(24) = c(24) + 1._dp
        end if

        ! calculate hmin - use default unless value prescribed
        c(13) = c(3)
        if (abs(c(3)) < tiny(1._dp)) then  ! if c(3) == 0
            ! calculate default value of hmin
            ! first calculate weighted norm y - c(12) - as specified
            ! by the error control indicator c(1)
            temp = 0._dp
            if (abs(c(1) - 1._dp) < tiny(1._dp)) then  ! if c(1) == 1
                ! absolute error control - weights are 1
                do k = 1, n
                    temp = max(temp, abs(y(k)))
                end do
                c(12) = temp
            else if (abs(c(1) - 2._dp) < tiny(1._dp)) then  ! if c(1) == 2
                ! relative error control - weights are 1/abs(y(k)) so
                ! weighted norm y is 1
                c(12) = 1._dp
            else if (abs(c(1) - 3._dp) < tiny(1._dp)) then  ! if c(1) == 3
                ! weights are 1/max(c(2), abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(2))
                end do
                c(12) = min(temp, 1._dp)
            else if (abs(c(1) - 4._dp) < tiny(1._dp)) then  ! if c(1) == 4
                ! weights are 1/max(c(k + 30), abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(k + 30))
                end do
                c(12) = min(temp, 1._dp)
            else if (abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 5
                ! weights are 1 / c(k + 30)
                do k = 1, n
                    temp = max(temp, abs(y(k)) / c(k + 30))
                end do
                c(12) = temp
            else
                ! default case - weights are 1/max(1, abs(y(k)))
                do k = 1, n
                    temp = max(temp, abs(y(k)))
                end do
                c(12) = min(temp, 1._dp)
            end if
            c(13) = 10._dp * max(c(11), c(10) * max(c(12) / tol, abs(x)))
        end if

        ! calculate scale - use default unless value prescribed
        c(15) = c(5)
        if (abs(c(5)) < tiny(1._dp)) then  ! if c(5) == 0
            c(15) = 1._dp
        end if

        ! calculate hmax - consider 4 cases
        ! case 1 both hmax and scale prescribed
        if (abs(c(6)) >= tiny(1._dp) .and. abs(c(5)) >= tiny(1._dp)) then  ! if c(6) /= 0 .and. c(5) /= 0
            c(16) = min(c(6), 2._dp / c(5))
        ! case 2 - hmax prescribed, but scale not
        else if (abs(c(6)) >= tiny(1._dp) .and. abs(c(5)) < tiny(1._dp)) then  ! if c(6) /= 0 .and. c(5) == 0
            c(16) = c(6)
        ! case 3 - hmax not prescribed, but scale is
        else if (abs(c(6)) < tiny(1._dp) .and. abs(c(5)) >= tiny(1._dp)) then  ! if c(6) == 0 .and. c(5) /= 0
            c(16) = 2._dp / c(5)
        ! case 4 - neither hmax nor scale is provided
        else if (abs(c(6)) < tiny(1._dp) .and. abs(c(5)) < tiny(1._dp)) then  ! if c(6) == 0 .and. c(5) == 0
            c(16) = 2._dp
        end if

!***********error return (with ind=-2) if hmin > hmax
        if (c(13) > c(16)) then
            ind = -2
            return
        end if

        ! calculate preliminary hmag - consider 3 cases
        if (ind <= 2) then
            ! case 1 - initial entry - use prescribed value of hstart, if
            ! any, else default
            c(14) = c(4)
            if (abs(c(4)) < tiny(1._dp)) c(14) = c(16) * tol**(1._dp / 6._dp)  ! if c(4) == 0
        else if (c(23) <= 1._dp) then
            ! case 2 - after a successful step, or at most  one  failure,
            ! use min(2, .9 * (tol/est)**(1/6)) * hmag, but avoid possible
            ! overflow. then avoid reduction by more than half.
            temp = 2._dp * c(14)
                if (tol < (2._dp / .9_dp)**6 * c(19)) then
                    temp = .9_dp * (tol / c(19))**(1._dp / 6._dp) * c(14)
                end if
            c(14) = max(temp, .5_dp * c(14))
        else
            ! case 3 - after two or more successive failures
            c(14) = .5_dp * c(14)
        end if

        ! check against hmax
        c(14) = min(c(14), c(16))

        ! check against hmin
        c(14) = max(c(14), c(13))

!***********interrupt no 1 (with ind=4) if requested
        if (abs(c(8)) >= tiny(1._dp)) then  ! if c(8) /= 0
            ind = 4
            return
        end if
!       resume here on re-entry with ind == 4   ........re-entry..
 1111   continue

        ! calculate hmag, xtrial - depending on preliminary hmag, xend
        if (c(14) < abs(xend - x)) then
            ! do not step more than half way to xend
            c(14) = min(c(14), .5_dp * abs(xend - x))
            c(17) = x + dsign(c(14), xend - x)
        else
            ! hit xend exactly
            c(14) = abs(xend - x)
            c(17) = xend
        end if

        ! calculate htrial
        c(18) = c(17) - x

!       end stage 1

!       ***************************************************************
!       * stage 2 - calculate ytrial (adding 7 to no of  fcn  evals). *
!       * w(*,2), ... w(*,8)  hold  intermediate  results  needed  in *
!       * stage 3. w(*,9) is temporary storage until finally it holds *
!       * ytrial.                                                     *
!       ***************************************************************

        temp = c(18) / 1398169080000._dp

        do k = 1, n
            w(k,9) = y(k) + temp * w(k,1) * 233028180000._dp
        end do
        call fcn(n, x + c(18) / 6._dp, w(1,9), w(1,2))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 74569017600._dp &
                                      + w(k,2) * 298276070400._dp  )
        end do
        call fcn(n, x + c(18) * (4._dp / 15._dp), w(1,9), w(1,3))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 1165140900000._dp &
                - w(k,2) * 3728450880000._dp &
                + w(k,3) * 3495422700000._dp )
        end do
        call fcn(n, x + c(18) * (2._dp / 3._dp), w(1,9), w(1,4))

        do k = 1, n
            w(k,9) = y(k) + temp * ( - w(k,1) * 3604654659375._dp &
                + w(k,2) * 12816549900000._dp &
                - w(k,3) * 9284716546875._dp &
                + w(k,4) * 1237962206250._dp )
        end do
        call fcn(n, x + c(18) * (5._dp / 6._dp), w(1,9), w(1,5))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 3355605792000._dp &
                - w(k,2) * 11185352640000._dp &
                + w(k,3) * 9172628850000._dp &
                - w(k,4) * 427218330000._dp &
                + w(k,5) * 482505408000._dp  )
        end do
        call fcn(n, x + c(18), w(1,9), w(1,6))

        do k = 1, n
            w(k,9) = y(k) + temp * ( - w(k,1) * 770204740536._dp &
                + w(k,2) * 2311639545600._dp &
                - w(k,3) * 1322092233000._dp &
                - w(k,4) * 453006781920._dp &
                + w(k,5) * 326875481856._dp  )
        end do
        call fcn(n, x + c(18) / 15._dp, w(1,9), w(1,7))

        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 2845924389000._dp &
                - w(k,2) * 9754668000000._dp &
                + w(k,3) * 7897110375000._dp &
                - w(k,4) * 192082660000._dp &
                + w(k,5) * 400298976000._dp &
                + w(k,7) * 201586000000._dp  )
        end do
        call fcn(n, x + c(18), w(1,9), w(1,8))

        ! calculate ytrial, the extrapolated approximation and store
        ! in w(*,9)
        do k = 1, n
            w(k,9) = y(k) + temp * (   w(k,1) * 104862681000._dp &
                + w(k,3) * 545186250000._dp &
                + w(k,4) * 446637345000._dp &
                + w(k,5) * 188806464000._dp &
                + w(k,7) * 15076875000._dp &
                + w(k,8) * 97599465000._dp   )
        end do

        ! add 7 to the no of fcn evals
        c(24) = c(24) + 7._dp

!       end stage 2

!       ***************************************************************
!       * stage 3 - calculate the error estimate est. first calculate *
!       * the  unweighted  absolute  error  estimate vector (per unit *
!       * step) for the unextrapolated approximation and store it  in *
!       * w(*,2).  then  calculate the weighted max norm of w(*,2) as *
!       * specified by the error  control  indicator  c(1).  finally, *
!       * modify  this result to produce est, the error estimate (per *
!       * unit step) for the extrapolated approximation ytrial.       *
!       ***************************************************************

        ! calculate the unweighted absolute error estimate vector
        do k = 1, n
            w(k,2) = (   w(k,1) * 8738556750._dp &
                + w(k,3) * 9735468750._dp &
                - w(k,4) * 9709507500._dp &
                + w(k,5) * 8582112000._dp &
                + w(k,6) * 95329710000._dp &
                - w(k,7) * 15076875000._dp &
                - w(k,8) * 97599465000._dp) / 1398169080000._dp
        end do

        ! calculate the weighted max norm of w(*,2) as specified by
        ! the error control indicator c(1)
        temp = 0._dp
        if (abs(c(1) - 1._dp) < tiny(1._dp)) then  ! if c(1) == 1
            ! absolute error control
            do k = 1, n
                temp = max(temp, abs(w(k,2)))
            end do
        else if (abs(c(1) - 2._dp) < tiny(1._dp)) then  ! if c(1) == 2
            ! relative error control
            do k = 1, n
                temp = max(temp, abs(w(k,2) / y(k)))
            end do
        else if (abs(c(1) - 3._dp) < tiny(1._dp)) then  ! if c(1) == 3
            ! weights are 1/max(c(2), abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) &
                    / max(c(2), abs(y(k))) )
            end do
        else if (abs(c(1) - 4._dp) < tiny(1._dp)) then  ! if c(1) == 4
            ! weights are 1/max(c(k + 30), abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) / max(c(k + 30), abs(y(k))) )
            end do
        else if (abs(c(1) - 5._dp) < tiny(1._dp)) then  ! if c(1) == 5
            ! weights are 1/c(k + 30)
            do k = 1, n
                temp = max(temp, abs(w(k,2) / c(k + 30)))
            end do
        else
            ! default case - weights are 1/max(1, abs(y(k)))
            do k = 1, n
                temp = max(temp, abs(w(k,2)) &
                    / max(1._dp, abs(y(k))) )
            end do
        end if

        ! calculate est - (the weighted max norm of w(*,2)) * hmag * scale
        ! - est is intended to be a measure of the error  per  unit
        ! step in ytrial
        c(19) = temp * c(14) * c(15)

!       end stage 3

!       ***************************************************************
!       * stage 4 - make decisions.                                   *
!       ***************************************************************

        ! set ind=5 if step acceptable, else set ind=6
        ind = 5
        if (c(19) > tol) ind = 6

!***********interrupt no 2 if requested
        if (abs(c(9)) >= tiny(1._dp)) then  ! if c(9) /= 0
            return
        end if
!       resume here on re-entry with ind == 5 or 6   ...re-entry..
 2222   continue

        if (ind /= 6) then
            ! step accepted (ind == 5), so update x, y from xtrial,
            ! ytrial, add 1 to the no of successful steps, and set
            ! the no of successive failures to zero
            x = c(17)
            do k = 1, n
                y(k) = w(k,9)
            end do
            c(22) = c(22) + 1._dp
            c(23) = 0._dp
!**************return(with ind=3, xend saved, flag set) if x == xend
            if (abs(x - xend) < tiny(1._dp)) then  ! if x == xend
                ind = 3
                c(20) = xend
                c(21) = 1._dp
                return
            end if
        else
            ! step not accepted (ind == 6), so add 1 to the no of
            ! successive failures
            c(23) = c(23) + 1._dp
!**************error return (with ind=-3) if hmag <= hmin
            if (c(14) <= c(13)) then
                ind = -3
                return
            end if
        end if

!       end stage 4

        go to 99999
!       end loop

!       begin abort action
  500   continue

        write(*,505) ind, tol, x, &
                     n, c(13), xend, &
                     nw, c(16), c(20), &
                     c(22), &
                     c(23), &
                     c(24), &
                     (y(k), k = 1, n)
  505   format( /// 1h0, 58hcomputation stopped in DVERK with the following values    &
                  / 1h0, 5hind =, i4, 5x, 6htol  =, 1pd13.6, 5x, 11hx         =, 1pd22.15 &
                  / 1h , 5hn   =, i4, 5x, 6hhmin =, 1pd13.6, 5x, 11hxend      =, 1pd22.15 &
                  / 1h , 5hnw  =, i4, 5x, 6hhmax =, 1pd13.6, 5x, 11hprev xend =, 1pd22.15 &
                  / 1h0, 14x, 27hno of successful steps    =, 0pf8.0 &
                  / 1h , 14x, 27hno of successive failures =, 0pf8.0 &
                  / 1h , 14x, 27hno of function evals      =, 0pf8.0 &
                  / 1h0, 23hthe components of y are &
                 // (1h , 1p5d24.15) )

!       end abort action

      end subroutine dverk
end module ode_solver

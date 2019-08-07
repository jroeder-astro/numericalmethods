module functions
implicit none
contains
real function e_temp(T) result(y)
      implicit none
    real,intent(in)::T
    real,parameter::lc=1.64603517E-1,T_p=1E+9,T_g=1E+7,q=2.7E-10
    y = q*sqrt(T)*(T-T_g)-lc*(T_p/T-1)
end function e_temp
end module functions

module method
implicit none
contains
logical function crit(one,two) result(answer)
      implicit none
    real,intent(in)::one,two
    if (one*two<0) then
      answer=.true.
    else
      answer=.false.
    endif
end function crit

real function bisection(dynfunc,lower,upper,eps) result(root)
    implicit none
  real::lower,upper,eps
  real,external::dynfunc
  real::error=5.0,root_prev,low,up,mid
  logical::criterium

  do while (error>eps)
    root_prev = root
    low = dynfunc(lower)
    up = dynfunc(upper)
    mid = dynfunc((lower+upper)/2.0)
    
    criterium = crit(low,mid)
    if (criterium) then
      upper = (lower+upper)/2.0
    else
      lower = (lower+upper)/2.0
    endif

    root = (lower+upper)/2.0
    error = abs((root/root_prev)-1.0)
  enddo
end function bisection
end module method

program roots
use method
use functions
  implicit none
real::lower_bracket=1.6E+7,upper_bracket=2.0E+7,prec=1E-7
real::root_bisection

root_bisection = bisection(e_temp,lower_bracket,upper_bracket,prec)
print*,'Result:',root_bisection
end program roots

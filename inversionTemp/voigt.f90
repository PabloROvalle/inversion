function voigt(x,y)

!-------------------------------------------------------------------------
!----------DESCRIPTION FUNCTION VOIGT-------------------------------------
!     Calcule la fonction de Voigt a partir de l'algorithme decrit
!     par Humlicek JQSRT, 27, 437 (1982).
!     Calule la fonction complexe W(Z)=Exp(-Z*Z)*Erfc(-Z*Z)
!     dans le plan complexe superieur (cad pour y>=0., z=x+iy)
!     La partie reelle est la fonction de Voigt
!     L'erreur relative maximale sur les parties imaginaire et reelle
!     est <1e-4
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  implicit none
  real, dimension(:), intent(in) :: x
  real, intent(in) :: y
  real, dimension(size(x)) :: voigt, s
  complex, dimension(size(x)) :: t,u,voigt1
  real, parameter :: sqrtpi=0.5641896
!
  t = cmplx(reshape((/y/),(/size(x)/),pad=(/y/)),-x)
  s = Abs(x) + y
  where (s >= 15)
     voigt1 = t*0.5641896/(0.5+t*t)
     voigt = real(voigt1*sqrtpi)
  elsewhere (s >= 5.5)
     u=t*t
     voigt1 = t*(1.410474 + u*.5641896)/(.75 + u*(3.+u))
     voigt = real(voigt1*sqrtpi)
  elsewhere (0.195*Abs(x)-0.176 <= y)
     voigt1 = (16.4955+t*(20.20933+t*(11.96482+t*(3.778987+t*.5642236))))/&
          (16.4955+t*(38.82363+t*(39.27121+t*(21.69274+t*(6.699398+t)))))
     voigt = Real(voigt1*sqrtpi)
  elsewhere
     u = t*t
     voigt1 = Cexp(u)-t*(36183.31-u*(3321.9905-u*(1540.787-u*(219.0313-u &
          *(35.76683-u*(1.320522-u*.56419))))))/(32066.6-u*(24322.84-u &
          * (9022.228-u*(2186.181-u*(364.2191-u*(61.57037-u*(1.841439-u)))))))
     voigt = Real(voigt1*sqrtpi)
  end where

  return
End function voigt

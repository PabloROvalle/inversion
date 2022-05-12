program extract_opa

  use declaration
  implicit none

  integer :: k, n, rec

  real, dimension(nfreq*nrep) :: total
!--------------------------------------------------------------------------
  open(8,file='/data/opacite/jupiter/jupiter_ch4.opa',access='direct',&
       &recl=nfreq,status='old',form='unformatted')
  open(10,file='/data/opacite/jupiter/jupiter_ch6.opa',access='direct',&
       &recl=nfreq,status='new',form='unformatted')

  do k=1, 1
     rec = (k-1) * nrep
     do n = 1, nrep
        read(8,rec=rec+n) total(1+(n-1)*nfreq:n*nfreq)
        write(10,rec=rec+n) total(1+(n-1)*nfreq:n*nfreq)
     enddo
  enddo

  close(8)
  close(10)
 
end program extract_opa

!! Calculate diagonal elements of the kinetic operator matrix.
!! This subroutine is meant for testing purpose only.
!!
!! author: Fadjar Fathurrahman

SUBROUTINE init_K_diag()

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     D2jl_x => LF3d_D2jl_x, &
                     D2jl_y => LF3d_D2jl_y, &
                     D2jl_z => LF3d_D2jl_z, &
                     lin2xyz => LF3d_lin2xyz
  USE m_hamiltonian, ONLY : K_diag

  IMPLICIT NONE 
  INTEGER :: ip
  INTEGER :: i, j, k

  ALLOCATE( K_diag(Npoints) )

  DO ip = 1, Npoints
    i = lin2xyz(1,ip)
    j = lin2xyz(2,ip)
    k = lin2xyz(3,ip)
    K_diag(ip) = -0.5d0 * ( D2jl_x(i,i) + D2jl_y(j,j) + D2jl_z(k,k) )
  ENDDO

END SUBROUTINE 


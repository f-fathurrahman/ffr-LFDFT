!!>
!!> \section{Subroutine \texttt{init\_strfact\_shifted()}}
!!>
!!> The structure factor $S_{f}(\mathbf{G})$ for atom $i_a$ is calculated as
!!> \begin{equation}
!!> S_{f,i_a}(\mathbf{G}) = \exp\left[ -\imath \mathbf{G}\cdot\mathbf{R}_{i_a}\right]
!!> \label{eq:strfact}
!!> \end{equation}
!!> where $\mathbf{G}$ is reciprocal vector and $R_{i_a}$ is position of $i_a$-th atom.
!!>
!!>  This equation can be simplified into
!!> \begin{equation}
!!> S_{f,i_s} = \sum_{i_a}^{N_{a,s}} S_{f,i_a}
!!> \end{equation}
!!> where the sum is done over number of atoms of species $i_s$, which is denoted
!!> by $N_{a,s}$
!!>
!!> Note that, in this subroutine the position is shifted to fit the usual FFT grid.
SUBROUTINE init_strfact_shifted()

  USE m_atoms, ONLY : Na => Natoms, &
                      Xpos => AtomicCoords, &
                      Nspecies, &
                      atm2species, &
                      strf => StructureFactor
  USE m_LF3d, ONLY : Ng => LF3d_Npoints, &
                     Gv => LF3d_Gv, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z
  IMPLICIT NONE
  !
  INTEGER :: ia, isp, ig
  REAL(8) :: GX, shiftx, shifty, shiftz

  WRITE(*,*)
  WRITE(*,*) 'Calculating structure factor: shifted to Fourier grid'

!!> \begin{itemize}
!!>
!!> \item
!!> Calculate the shift with respect to the usual FFT grid.
  shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )

!!> \item
!!> Allocate memory for structure factor
  ALLOCATE( strf(Ng,Nspecies) )

!!> \item
!!> Initialize to zero
  strf(:,:) = cmplx(0.d0,0.d0,kind=8)

!!>
!!> \item
!!> Implementation of equation \ref{eq:strfact} is done here.
  DO ia = 1,Na
    isp = atm2species(ia)
    DO ig = 1,Ng
      GX = (Xpos(1,ia)-shiftx)*Gv(1,ig) + (Xpos(2,ia)-shifty)*Gv(2,ig) + &
           (Xpos(3,ia)-shiftz)*Gv(3,ig)
      strf(ig,isp) = strf(ig,isp) + cmplx( cos(GX), -sin(GX), kind=8 )
    ENDDO
  ENDDO

!!> \item
!!> Finish
  CALL flush(6)

END SUBROUTINE
!!> \end{itemize}

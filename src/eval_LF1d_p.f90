!!>
!!> \section{Function \texttt{eval\_LF1d\_p}}
!!> 
!!> Evaluate periodic Lagrange function:
!!> \begin{equation}
!!> \phi_{\alpha}(x) = \frac{1}{\sqrt{NL}}
!!> \sum_{i=1}^{N} \cos\left[k_{i}(x - x_{\alpha})\right]
!!> \end{equation}
!!> with $i = 1,2,\ldots,N$ and
!!> \begin{equation}
!!> k_{i} = \frac{2\pi(i - N' - 1)}{L}
!!> \end{equation}
!!> $N' = (N-1)/2$ and $N$ must be an odd number.
!!>
FUNCTION eval_LF1d_p(N, L, grid, ibf, x) RESULT(ff)
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  ! arguments
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: grid(N)
  INTEGER :: ibf
  REAL(8) :: x, ff
  !
  REAL(8) :: pre1
  INTEGER :: ii
  !
  pre1 = 1.d0/sqrt(N*L)
  ff = 0.d0
  DO ii = 1,N
    ff = ff + cos( PI*( 2*ii - N - 1 )*( x - grid(ibf) )/ L )
  ENDDO
  ff = ff*pre1
END FUNCTION

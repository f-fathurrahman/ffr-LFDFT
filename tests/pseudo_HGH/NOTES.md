Initializing HGH pseudopotential

```fortran
  CALL hgh_init( ps, filename )
  CALL hgh_process(ps)

  ! do something ...

  CALL hgh_end( ps )
```

Evaluating local pseudopotential:

```fortran
FUNCTION vlocalr_scalar(r, p)
  TYPE(hgh_t), INTENT(in) :: p
  REAL(8), INTENT(in)     :: r
  REAL(8)                 :: vlocalr_scalar
```

Evaluating projectors for nonlocal pseudopotential

```fortran
FUNCTION projectorr_scalar(r, p, i, l)
  TYPE(hgh_t), INTENT(IN) :: p
  REAL(8), INTENT(IN)     :: r
  INTEGER, INTENT(IN)     :: i, l
  REAL(8)                 :: projectorr_scalar
```



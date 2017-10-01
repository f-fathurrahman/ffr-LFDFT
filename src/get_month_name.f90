FUNCTION get_month_name( imonth ) result( name_ )
  IMPLICIT NONE
  INTEGER :: imonth
  CHARACTER(9) :: name_

  SELECT CASE( imonth )
  CASE( 1 )
    name_ = 'January'
  CASE( 2 )
    name_ = 'February'
  CASE( 3 )
    name_ = 'March'
  CASE( 4 )
    name_ = 'April'
  CASE( 5 )
    name_ = 'May'
  CASE( 6 )
    name_ = 'June'
  CASE( 7 )
    name_ = 'July'
  CASE( 8 )
    name_ = 'August'
  CASE( 9 ) 
    name_ = 'September'
  CASE( 10 )
    name_ = 'October'
  CASE( 11 )
    name_ = 'November'
  CASE( 12 )
    name_ = 'December'
  CASE DEFAULT
    WRITE(*,*) 'Invalid input to get_month_names: ', imonth
    STOP 
  END SELECT 

END FUNCTION 

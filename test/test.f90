program test

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64, &
        stdout => OUTPUT_UNIT, &
        compiler_version, &
        compiler_options

    use type_FortranFFTW, only: &
        FortranFFTW

    ! Explicity typing only
    implicit none

    !---------------------------------------------------------------------------------
    ! Dictionary: calling arguments
    !---------------------------------------------------------------------------------
    integer (ip), parameter :: RANGE_N(*) = [ 8, 16, 32, 64, 128 ]
    integer (ip)             :: n !! Counter
    !---------------------------------------------------------------------------------

    do n = 1, size( RANGE_N )
        associate( sample_size => RANGE_N(n))
            write( stdout, '(A)') ''
            write( stdout, '(A,I3)') 'For sample size n = ', sample_size
            call test_derivatives( sample_size )
        end associate
        write( stdout, '(A)') ''
        write( stdout, '(A)') '****************************************'
        write( stdout, '(A)') ''
    end do

    ! Print compiler info
    write( stdout, '(A)' ) ' '
    write( stdout, '(4A)' ) 'This result was compiled by ', &
        compiler_version(), ' using the options ', &
        compiler_options()
    write( stdout, '(A)' ) ' '

contains
    !
    !*****************************************************************************************
    !
    subroutine get_first_derivative( periodic_array, derivative )
        !
        ! Purpose:
        !
        ! Implements a spectral approximation for the first derivative
        ! of a periodic arrays using the discrete fast fourier transform (FFT)
        !
        !---------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !---------------------------------------------------------------------------------
        real (wp), intent (in)  :: periodic_array(:)
        real (wp), intent (out) :: derivative(:)
        !---------------------------------------------------------------------------------
        ! Dictionary: local variables
        !---------------------------------------------------------------------------------
        type (FortranFFTW)         :: fftw
        integer (ip)               :: i !! Counter
        real (wp),    allocatable :: x(:)
        complex (wp), allocatable :: input(:)
        complex (wp), allocatable :: output(:)
        !---------------------------------------------------------------------------------

        associate( &
            n => size( periodic_array ), &
            f => periodic_array, &
            df => derivative &
            )

            ! Allocate memory
            allocate( x(n), input(n), output(n) )

            ! Initialize indices
            x = [ ((i-1), i = 1, n/2), ((-n+i-1), i = n/2+1, n) ]

            ! Set the wave number for the nyquist frequency to zero
            x( n/2+1 ) = 0.0_wp

            ! Set input
            input = cmplx( f, kind = wp )

            ! Perform discrete fast fourier transfrom
            call fftw%fft( n, input, output )

            ! Compute the derivative in spectral space and normalize the FFT
            do i = 1, n
                associate( normalization_constant => cmplx( 0.0_wp, x(i), kind = wp ) / n )
                    output(i) = normalization_constant * output (i)
                end associate
            end do

            ! Perform inverse discrete fast fourier transfrom
            call fftw%ifft( n, output, input )

            ! Set derivative
            df = real( input, kind = wp )

            ! Free memory
            deallocate( x, input, output )

            ! Destroy type
            call fftw%destroy()

        end associate

    end subroutine get_first_derivative
    !
    !*****************************************************************************************
    !
    subroutine get_second_derivative( periodic_array, derivative )
        !
        ! Purpose:
        !
        ! Implements a spectral approximation for the first derivative
        ! of a periodic arrays using the discrete fast fourier transform (FFT)
        !
        !---------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !---------------------------------------------------------------------------------
        real (wp), intent (in)  :: periodic_array(:)
        real (wp), intent (out) :: derivative(:)
        !---------------------------------------------------------------------------------
        ! Dictionary: local variables
        !---------------------------------------------------------------------------------
        type (FortranFFTW)         :: fftw
        integer (ip)               :: i !! Counter
        real (wp),    allocatable :: x(:)
        complex (wp), allocatable :: input(:)
        complex (wp), allocatable :: output(:)
        !---------------------------------------------------------------------------------

        associate( &
            n => size( periodic_array ), &
            f => periodic_array, &
            df => derivative &
            )

            ! Allocate memory
            allocate( x(n), input(n), output(n) )

            ! Initialize indices
            x = [ (-(i-1)**2, i = 1, n/2), (-(-n+i-1)**2, i = n/2+1, n) ]

            ! Set the wave number for the nyquist frequency to zero
            x( n/2+1 ) = 0.0_wp

            ! Set input
            input = cmplx( f, kind = wp )

            ! Perform discrete fast fourier transfrom
            call fftw%fft( n, input, output )

            ! Compute the derivative in spectral space and normalize the FFT
            do i = 1, n
                associate( normalization_constant => cmplx( x(i), 0.0_wp, kind = wp) / n )
                    output(i) = normalization_constant * output (i)
                end associate
            end do

            ! Perform inverse discrete fast fourier transfrom
            call fftw%ifft( n, output, input )

            ! Set derivative
            df = real( input, kind = wp )

            ! Free memory
            deallocate( x, input, output )

            ! Destroy type
            call fftw%destroy()

        end associate

    end subroutine get_second_derivative
    !
    !*****************************************************************************************
    !
    subroutine test_derivatives ( n )
        !
        ! Purpose:
        !
        !---------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !---------------------------------------------------------------------------------
        integer (ip), intent (in) :: n
        !---------------------------------------------------------------------------------
        ! Dictionary: local variables
        !---------------------------------------------------------------------------------
        real (wp), allocatable :: f(:)
        real (wp), allocatable :: df_exact(:)
        real (wp), allocatable :: df(:)
        real (wp), allocatable :: ddf_exact(:)
        real (wp), allocatable :: ddf(:)
        real (wp), parameter   :: TWO_PI = 2.0_wp * acos( -1.0_wp )
        integer (ip)            :: i !! Counter
        !---------------------------------------------------------------------------------

        ! Allocate arrays
        allocate( f(n), df_exact(n), df(n), ddf_exact(n), ddf(n) )

        ! Set periodic array and corresponding derivatives
        associate( MESH => TWO_PI / n )
            do i = 1,n
                associate( x => real(i - 1, kind = wp) * MESH )

                    f(i) = &
                        cos( 2.0_wp * x) + sin( 3.0_wp * x)


                    df_exact(i) = &
                        -2.0_wp * sin( 2.0_wp * x) + 3.0_wp * cos ( 3.0_wp * x)


                    ddf_exact(i) = &
                        -4.0_wp * cos( 2.0_wp * x) - 9.0_wp * sin ( 3.0_wp * x)
                end associate
            end do
        end associate

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) 'For f(x) = cos(2x) + sin(3x)'

        ! For the first derivative
        call get_first_derivative( f, df )

        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) 'First derivative from periodic periodic_array'
        write( stdout, * ) 'Discretization error = ', maxval( abs( df - df_exact) )

        ! For the second derivative
        call get_second_derivative( f, ddf )


        write( stdout, '(A)' ) ''
        write( stdout, '(A)' ) 'Second derivative from the periodic periodic_array'
        write( stdout, * ) 'Discretization error = ', maxval( abs( ddf - ddf_exact) )

        ! Free memory
        deallocate( f, df, ddf, df_exact, ddf_exact )

    end subroutine test_derivatives
    !
    !*****************************************************************************************
    !
end program test

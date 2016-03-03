module type_FortranFFTW

    use, intrinsic :: iso_fortran_env, only: &
        ip => INT32, &
        wp => REAL64

    use, intrinsic :: iso_c_binding

    ! Explicit typing only
    implicit none

    include 'fftw3.f03'

    ! Everything is private unless stated otherwise
    private

    type, public :: FortranFFTW
        !---------------------------------------------------------------------------------
        ! Class variables
        !---------------------------------------------------------------------------------
        logical,      private :: forward_plan_usable = .false.
        logical,      private :: backward_plan_usable = .false.
        type (C_PTR), private :: forward_plan
        type (C_PTR), private :: backward_plan
        !---------------------------------------------------------------------------------
    contains
        !---------------------------------------------------------------------------------
        ! Class methods
        !---------------------------------------------------------------------------------
        procedure, public  :: fft => compute_one_dimensional_discrete_fourier_transform
        procedure, public  :: ifft => compute_one_dimensional_inverse_discrete_fourier_transform
        procedure, private :: create_forward_plan
        procedure, private :: destroy_forward_plan
        procedure, private :: create_backward_plan
        procedure, private :: destroy_backward_plan
        procedure, public  :: destroy => destroy_FortranFFTW
        final              :: finalize_FortranFFTW
        !---------------------------------------------------------------------------------
    end type FortranFFTW

contains
    !
    !*****************************************************************************************
    !
    subroutine compute_one_dimensional_discrete_fourier_transform( this, n, input, output )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in out)  :: this
        integer (ip),        intent (in)      :: n
        complex (wp),        intent (in out)  :: input(:)
        complex (wp),        intent (out)     :: output(:)
        !--------------------------------------------------------------------------------

        ! Ensure that the forward plan is usable
        if ( .not. this%forward_plan_usable ) then
            call this%create_forward_plan( n, input , output)
        end if

        ! Perform the forward transform
        call fftw_execute_dft( this%forward_plan, input, output )

    end subroutine compute_one_dimensional_discrete_fourier_transform
    !
    !*****************************************************************************************
    !
    subroutine compute_one_dimensional_inverse_discrete_fourier_transform( this, n, output, input )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in out)  :: this
        integer (ip),        intent (in)      :: n
        complex (wp),        intent (in out)  :: output(:)
        complex (wp),        intent (in out)  :: input(:)
        !--------------------------------------------------------------------------------

        ! Ensure that the backward plan is usable
        if ( .not. this%backward_plan_usable ) then
            call this%create_backward_plan( n, output, input )
        end if

        ! Perform the backward transform
        call fftw_execute_dft( this%backward_plan, output, input )

    end subroutine compute_one_dimensional_inverse_discrete_fourier_transform
    !
    !*****************************************************************************************
    !
    subroutine create_forward_plan( this, n, input, output )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in out) :: this
        integer (ip),        intent (in)     :: n
        complex (wp),        intent (in out) :: input(:)
        complex (wp),        intent (in out) :: output(:)
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy_forward_plan()

        ! Create forward plan
        this%forward_plan = fftw_plan_dft_1d( n, input, output, FFTW_FORWARD, FFTW_ESTIMATE )

    end subroutine create_forward_plan
    !
    !*****************************************************************************************
    !
    subroutine destroy_forward_plan( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in) :: this
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( .not. this%forward_plan_usable ) return

        ! Destroy the plan
        call fftw_destroy_plan( this%forward_plan )

    end subroutine destroy_forward_plan
    !
    !*****************************************************************************************
    !
    subroutine create_backward_plan( this, n, output, input )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in out)  :: this
        integer (ip),        intent (in)      :: n
        complex (wp),        intent (in out)  :: output(:)
        complex (wp),        intent (in out)  :: input(:)
        !--------------------------------------------------------------------------------

        ! Ensure that object is usable
        call this%destroy_backward_plan()

        ! Create backward plan
        this%backward_plan = fftw_plan_dft_1d( n, output, input, FFTW_BACKWARD, FFTW_ESTIMATE )

    end subroutine create_backward_plan
    !
    !*****************************************************************************************
    !
    subroutine destroy_backward_plan( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in) :: this
        !--------------------------------------------------------------------------------

        ! Check if object is usable
        if ( .not. this%backward_plan_usable ) return

        ! Destroy the plan
        call fftw_destroy_plan( this%backward_plan )

    end subroutine destroy_backward_plan
    !
    !*****************************************************************************************
    !
    subroutine destroy_FortranFFTW( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        class (FortranFFTW), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy_forward_plan()

        call this%destroy_backward_plan()

    end subroutine destroy_FortranFFTW
    !
    !*****************************************************************************************
    !
    subroutine finalize_FortranFFTW( this )
        !--------------------------------------------------------------------------------
        ! Dictionary: calling arguments
        !--------------------------------------------------------------------------------
        type (FortranFFTW), intent (in out) :: this
        !--------------------------------------------------------------------------------

        call this%destroy()

    end subroutine finalize_FortranFFTW
    !
    !*****************************************************************************************
    !
end module type_FortranFFTW

    interface
        subroutine cg_method(A,X,B,n)
            double precision, dimension(n,n) :: A
            double precision, dimension(n) :: B, X 
        end subroutine

        subroutine disp_vec(v)
            double precision, dimension(:), allocatable  :: v
        end subroutine

        subroutine disp_mat(t)
            double precision, dimension(:,:), allocatable  :: t
        end subroutine
    end interface
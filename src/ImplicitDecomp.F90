module implicit_decomp

    use param
    use decomp_2d

    implicit none

    type(DECOMP_INFO)   :: decomp_diff
    integer             :: dxst(3)
    integer             :: dxen(3)
    integer             :: dyst(3)
    integer             :: dyen(3)
    integer             :: dzst(3)
    integer             :: dzen(3)

    contains

    subroutine InitImplicitDecomp

        use param
        use decomp_2d
        
        implicit none

        character*4         :: rankname
        character*50        :: filename

        call decomp_info_init(nx,nym,nzm,decomp_diff)
        dxst(:) = decomp_diff%xst(:)
        dxen(:) = decomp_diff%xen(:)
        dyst(:) = decomp_diff%yst(:)
        dyen(:) = decomp_diff%yen(:)
        dzst(:) = decomp_diff%zst(:)
        dzen(:) = decomp_diff%zen(:)

        write(rankname,'(I4.4)') nrank
        filename = trim('Debug/'//rankname//'.out')

        open(1000+nrank,file=filename,status='unknown',position='append',access='sequential')
        write(1000+nrank,'(1(A10,I4))') ' nrank = ',nrank!,'prow = ',dims(1),'pcol = ',dims(2)
        write(1000+nrank,'(6(A12,I4))') &
            ' pxst(1) = ',xstart(1),&
            ' pxen(1) = ',xend(1),&
            ' pxst(2) = ',xstart(2),&
            ' pxst(2) = ',xend(2),&
            ' pxst(3) = ',xstart(3),&
            ' pxst(3) = ',xend(3)
        write(1000+nrank,'(6(A12,I4))') &
            ' pyst(1) = ',ystart(1),&
            ' pyen(1) = ',yend(1),&
            ' pyst(2) = ',ystart(2),&
            ' pyst(2) = ',yend(2),&
            ' pyst(3) = ',ystart(3),&
            ' pyst(3) = ',yend(3)
        write(1000+nrank,'(6(A12,I4))') &
            ' pzst(1) = ',zstart(1),&
            ' pzen(1) = ',zend(1),&
            ' pzst(2) = ',zstart(2),&
            ' pzst(2) = ',zend(2),&
            ' pzst(3) = ',zstart(3),&
            ' pzst(3) = ',zend(3)
        write(1000+nrank,'(6(A12,I4))') &
            ' dxst(1) = ',decomp_diff%xst(1),&
            ' dxen(1) = ',decomp_diff%xen(1),&
            ' dxst(2) = ',decomp_diff%xst(2),&
            ' dxst(2) = ',decomp_diff%xen(2),&
            ' dxst(3) = ',decomp_diff%xst(3),&
            ' dxst(3) = ',decomp_diff%xen(3)
        write(1000+nrank,'(6(A12,I4))') &
            ' dyst(1) = ',decomp_diff%yst(1),&
            ' dyen(1) = ',decomp_diff%yen(1),&
            ' dyst(2) = ',decomp_diff%yst(2),&
            ' dyst(2) = ',decomp_diff%yen(2),&
            ' dyst(3) = ',decomp_diff%yst(3),&
            ' dyst(3) = ',decomp_diff%yen(3)
        write(1000+nrank,'(6(A12,I4))') &
            ' dzst(1) = ',decomp_diff%zst(1),&
            ' dzen(1) = ',decomp_diff%zen(1),&
            ' dzst(2) = ',decomp_diff%zst(2),&
            ' dzst(2) = ',decomp_diff%zen(2),&
            ' dzst(3) = ',decomp_diff%zst(3),&
            ' dzst(3) = ',decomp_diff%zen(3)
        close(1000+nrank)

    end subroutine InitImplicitDecomp

end module implicit_decomp


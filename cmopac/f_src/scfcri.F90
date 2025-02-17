subroutine scfcri (selcon, isgood)
    use molkst_C, only: numcal, keywrd, efield
    use chanel_C, only: iw
    use reada_I
 !   use common_convrg, only: scfref
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    double precision, intent (inout) :: selcon
    integer, intent (in) :: isgood
   !
   !.. Local Scalars ..
    logical :: precis
    integer :: i
    integer :: icalcn = 0
    double precision, save  :: scfcrt, scfref 
    
    if (isgood == 1) then
        icalcn = 0
        return
    endif

    if (icalcn /= numcal) then
      icalcn = numcal
      !
      !   DETERMINE THE SELF-CONSISTENCY CRITERION
      !
      ! SCFCRT IS MACHINE-PRECISION DEPENDENT
      !
      !   IF SCFCRT IS CHANGED, THEN ALSO CHANGE DEFAULT IN WRTKEY
      !
      scfcrt = 1.d-2
      !
      !  THE USER CAN STATE THE SCF CRITERION, IF DESIRED.
      !
      i = Index (keywrd, " TS") + Index (keywrd, " FORCETS") + Index (keywrd, " IRC=")
      if (i /= 0) scfcrt = 1.d-3
      precis = (Index (keywrd, " PRECIS") /= 0)
      i = Index (keywrd, " RELSCF")
      if (i /= 0) then
        scfcrt = reada (keywrd, i) * scfcrt
        scfref = scfcrt
        write (iw, "('  SCF CRITERION =',G14.4)") scfcrt
        if (scfcrt < 1.d-5) then
10000     format (//2 x, " THERE IS A RISK OF INFINITE LOOPING WITH", &
         & " THE SCFCRT LESS THAN 1.D-5")
          write (iw, 10000)
        end if
      end if
      i = Index (keywrd, " SCFCRT")
      if (i /= 0) then
        scfcrt = reada (keywrd, i)
        scfref = scfcrt
        write (iw, "('  SCF CRITERION =',G14.4)") scfcrt
        if (scfcrt < 1.d-5) then
          write (iw, "(//2x,' THERE IS A RISK OF INFINITE LOOPING WITH', &
         & ' THE SCFCRT LESS THAN 1.D-5')")
        end if
      end if
      !
      !   SELF-CONSISTENCY CRITERIA: SELCON IS IN KCAL/MOL
      !   IF GNORM IS LARGE, MAKE SELCON BIGGER
      !
      if (precis) then
        scfcrt = scfcrt * 0.01d0
        scfref = scfcrt
      end if
      !
      ! For polarizability calculations, the default convergence should be 
      ! tightened
      !
      if (Index (keywrd, " POLAR") /= 0 .and. scfref == 0.0d0) then
        scfcrt = 1.d-4
      end if
      selcon = scfcrt
    else
      !
      !  IF POLARIZATION IS BEING CALCULATED, TIGHTEN SCF CRITERION
      !
      if (Abs (efield(1))+Abs(efield(2))+Abs(efield(3)) > 1.d-6) then
        selcon = 1.d-4
      end if
      if (scfref /= 0.d0) then
        selcon = scfref
      end if
    end if
end subroutine scfcri

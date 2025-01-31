 subroutine init_filenames
    use molkst_C, only: jobnam, gui, line
    use chanel_C, only : output_fn, restart_fn, brillouin_fn, &
     & density_fn, log_fn, end_fn, archive_fn, esp_fn, ump_fn, &
     mep_fn, pol_fn, gpt_fn, esr_fn, input_fn, xyz_fn, syb_fn, &
     cosmo_fn
    implicit none
    integer :: text_length
    text_length = len_trim (jobnam)
    line = trim(jobnam)
    call upcase(line, len_trim(jobnam))
    if (text_length > 3) then
      if (jobnam(text_length - 3:text_length - 3) == ".") then
        if ((line(text_length - 3:text_length) == ".MOP") .or. &
            (line(text_length - 3:text_length) == ".DAT") .or. &
            (line(text_length - 3:text_length) == ".ARC") .or. &
            (line(text_length - 3:text_length) == ".PDB") .or. &
            (line(text_length - 3:text_length) == ".ENT") .or. &
            (line(text_length - 3:text_length) == ".NEW"))      &
        text_length = text_length - 4
      end if
    end if
!
!  Set up the name of the files that are related to the input file name
!

    density_fn    = jobnam(1:text_length) // ".den"

  end subroutine init_filenames

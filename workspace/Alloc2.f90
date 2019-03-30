
    allocate( gstif(ntotv,ntotv),       stat = alloc_stat(compt_alloc+29))    !Tableau 29
!     allocate( gstif_backup(ntotv,ntotv),       stat = alloc_stat(compt_alloc+29))    !Tableau 29
    allocate( rgstif(nvlib,nvlib),       stat = alloc_stat(compt_alloc+31))    !Tableau 29
    allocate( gload(ntotv),       stat = alloc_stat(compt_alloc+32))    !Tableau 28
    allocate( gload_backup(ntotv),       stat = alloc_stat(compt_alloc+32))    !Tableau 28
    allocate( rgload(nvlib),       stat = alloc_stat(compt_alloc+35))    !Tableau 29
!     allocate( rgload_backup(nvlib),       stat = alloc_stat(compt_alloc+35))    !Tableau 29
    allocate( rasdis(nvlib),       stat = alloc_stat(compt_alloc+36))    !Tableau 28
    allocate( vec_residu(nvlib),       stat = alloc_stat(compt_alloc+37))    !Tableau 29
    allocate( vec_residu_pd(nvlib),       stat = alloc_stat(compt_alloc+37))    !Tableau 29
!         allocate( vec_residu_backup(nvlib),       stat = alloc_stat(compt_alloc+37))    !Tableau 29

    allocate(rstfor(nvlib), stat = alloc_stat(compt_alloc+37))
    allocate( riffix(nvlib),       stat = alloc_stat(compt_alloc+35))    !Tableau 29
    allocate( rfixed(nvlib),       stat = alloc_stat(compt_alloc+35))    !Tableau 29

    !For lapack
    allocate( ipiv(nvlib),       stat = alloc_stat(compt_alloc+35))    !Tableau 29
    
    allocate( rgstif_lapack(nvlib,nvlib),       stat = alloc_stat(compt_alloc+31))    !Tableau 29
    allocate( rasdis_lapack(nvlib),       stat = alloc_stat(compt_alloc+36))    !Tableau 28
     allocate( vec_residu_lapack(nvlib),       stat = alloc_stat(compt_alloc+37))    !Tableau 29

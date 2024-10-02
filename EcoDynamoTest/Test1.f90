        program Test1
  
        USE ecodynamocpp_mod
        
        implicit none
        character(30) :: fname
        INTEGER(kind=8) :: EulerSteps,HollingsType
        DOUBLE PRECISION :: PARtop,KValue,Depth
        DOUBLE PRECISION :: PARbottom, PARopt,PlattPARopt
        DOUBLE PRECISION :: Pmax,beta,slope
        DOUBLE PRECISION :: a, b, c
        DOUBLE PRECISION :: LightLim,TempLim
        DOUBLE PRECISION :: PAR
        DOUBLE PRECISION :: Temp, TempAugRate, Tmin     
        DOUBLE PRECISION :: NutLim, CellQuota
        DOUBLE PRECISION :: MinCellQuota, HalfSaturation
        DOUBLE PRECISION :: NO3, KNO3,NH4, KNH4 
        DOUBLE PRECISION :: NO2, N, KN
        DOUBLE PRECISION :: GrazingLim,K,Prey,Lambda
        PARtop = 239.789527279837
        PAR = PARtop
        KValue = 0.01
        Depth = 100.0
        Pmax = 1.0
        beta = 0.001
        slope = 0.01
        EulerSteps = 10
        PARopt = 250
        PlattPARopt = (Pmax*log((beta+slope)/beta))/slope
        a = 1/(slope*PARopt*PARopt)
        b = 1/Pmax-2/(slope*PARopt)
        c = 1/slope
        Temp = 1.0
        TempAugRate = 0.069
        Tmin = 0.0
        CellQuota = 0.18
        MinCellQuota = 0.05
        HalfSaturation = 0.09
        NO3 = 1.0
        NO2 = 0.1
        NH4 = 0.1
        N = NO3
        KNO3 = 2
        KNH4 = 0.2
        KN = KNO3
        K = 1
        Prey = 2
        HollingsType = 0
        Lambda = 1.0
        write(*,*) 'Test1'
      ! Calling a C++ function through the ecodynamo_cpp
      ! interface.
      ! Arguments passed by reference, i.e.
      ! the C++ function expects memory addresses of the
      ! variables and not their values, which
      ! the function must retrieve through the
      ! memory addresses  
        write(*,*) 'PARtop =',PARtop
      ! Platt beginning  
        LightLim =  Platt1(PARtop,KValue,Depth, &
        & Pmax,beta,slope,EulerSteps) 

        write(*,*) 'Platt LightLim1 =',LightLim

        LightLim =  Platt2(PAR, &
        & Pmax,beta,slope)

        write(*,*) 'Platt LightLim2 =',LightLim
      ! Platt done 

      ! Steele beginning 
        LightLim =  Steele1(PARtop,             &
        &                  KValue,Depth,PARopt)

        write(*,*) 'Steele LightLim1 =',LightLim

        LightLim =  Steele2(PAR,PARopt)
  
        write(*,*) 'Steele LightLim2 =',LightLim
      ! Steele done

      ! Eilers & Peeters beginning

        LightLim =  EilersAndPeeters1(PARtop,   &
        &                 KValue,Depth,a,b,c)

        write(*,*) 'EilersAndPeeters LightLim1 =',LightLim

        LightLim =  EilersAndPeeters2(PAR,   &
        &                 a,b,c)

        write(*,*) 'EilersAndPeetersX LightLim2 =',LightLim


      ! Eiler & Peeters done


      ! Temperature limitation
          
        TempLim = TemperatureExponentialLimitation( &
        & Temp,TempAugRate, Tmin)
        write(*,*) 'Temperature limitation =',TempLim

      ! Temperature limitation done

      ! Internal nutrient limitation

        NutLim = InternalNutrientLimitation(CellQuota, &
        &  MinCellQuota, HalfSaturation)

        write(*,*) 'Inter nut limitation =',NutLim

      ! Internal nutrient limitation done

      ! External nutrient limitation
        NutLim = NitrateAndAmmoniumLimitation(NH4,KNH4,&
        & NO3, KNO3, NO2)

        write(*,*) 'NO3+NH4 limitation =',NutLim

      ! External nutrient limitation done

      ! Michaelis-Menten

        NutLim = MichaelisMentenLimitation(N,KN)

        write(*,*) 'Michaelis-Menten lim =',NutLim
      ! Michaelis-Menten done

      ! Hollings Type 0, 1, 2....

        write(*,*) 'Prey = ',Prey
        write(*,*) 'K =',K
        write(*,*) 'HollingsType =',HollingsType
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 0 =',GrazingLim

        HollingsType = 1
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 1 =',GrazingLim

        HollingsType = 2
        GrazingLim = Hollings(K,Prey,HollingsType)
        write(*,*) 'Grazing Hollings type 2 =',GrazingLim

        GrazingLim = Ivlev(Lambda,Prey)
        write(*,*) 'Grazing Ivlev =',GrazingLim

        end program Test1

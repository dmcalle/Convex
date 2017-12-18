program main

    implicit none
    integer,parameter::ndmax=10000,ncmax=ndmax,dimnd=ndmax+2
    real,parameter::rpi=4*atan(1.),c_final=0.65,rd=1d0,rd2=rd**2,&
    r_domaine=sqrt(ndmax/c_final)*rd,sqrx=sqrt(ndmax*rpi/c_final)*rd,sqry=sqrx!(sqrx,sqry)=size of square box // ndmax=total number of discs // rd = discs radius
    real,dimension(2),parameter::sqr=[sqrx,sqry]
    real,dimension(dimnd)::xd,yd ! positions : first  2 rows = positions x & y // number of the cluste/ third row = number of the next disc in the cluster
    integer,dimension(dimnd)::nextCH,listDiscsCH
     ! nbInClust(ic) nb of discs in cluster ic // first
    integer,dimension(ncmax)::nbCH,firstCH,lastCH
    real,dimension(ncmax)::coordinateStrip,thicknessStrip
    !common xd,yd,cluster,nextDisc,nbInClust,first,last
    logical::foundCluster,exterior,keepOn,afterPrevious,contact,switch,&
    alreadyIn,contactStrip,isInsideStrip,autoconvex,contactICbeforeIOC,contactIOCbeforeIC,&
    contactThisImage,contactOrNotStrip,ICinsideIOC,IOCinsideIC
    !integer,dimension(ndmax)::  ! 6=theoretical max : max number of discs  to be in contact with a given disc at the same time
    integer,dimension(dimnd,2)::nextThrough,limitDomain
    logical,dimension(ncmax)::isStrip
    !logical,dimension(ndmax)::limitDomain
    real::angle,maxAngle,minAngle,vect_product,ux,uy,angleBetweenTwoVectors,areaOfAClusterPeriodic&
    ,normOrientation,orthogonalDistance,orthogonalDistanceMin,orthogonalDistanceMax&
    ,coordinate,maxStrip,minStrip,proj,relativePosition,slope,area_tot,c,frac_area,&
    coordinateIC,coordinateIOC,coordinateStripOut,thicknessStripOut,newThickness,x1,y1&
    ,x2,y2,areaOfaStrip
    integer::i,ic,icmax,icount,idch,iMinAngle,iMaxAngle,iod,iseed,nbImages,nd,iMinOrthDistance,iMaxOrthDistance,&
    iBigLoop,idoch,ioc,nbInThisCH,idochAltern,&
    discContactStrip,discIOCContact,iContact,nbImagesContact,period,whichStripContact&
    ,whichStripContactIC,ich,ifile,nextTemp
    integer,dimension(2)::goingThrough,nextThroughTemp
    integer,dimension(3)::imagesContact
    integer,dimension(10,2)::whereImages
    logical,dimension(10)::exteriors,extLimits
    integer,dimension(2)::wherePeriodic,wherePeriodicNext,whereMinAngle,whereMaxAngle,throughFromEnd2End&
    ,wherePeriodicIOC,wherePeriodicIOCNext,wherePeriodicTemp
    integer,dimension(ncmax,2)::orientationStrip
    integer::mainOrientation
    logical::bool=.true.
    switch=.false.
    !write(*,*) sqrx
    !iseed=3546883
    !iseed=64247561
    !iseed=6448235 seed gros cluster
    !iseed=48847
    iseed=166555
    call srand(iseed)
    icmax=0
    do nd=1,ndmax
        write(*,*) nd,'density',c!2000!ndmax
!        x=r_domaine; y=r_domaine
!        do while (x**2+y**2>r_domaine**2)
!            x=r_domaine*(-1+2*rand(0))
!            y=r_domaine*(-1+2*rand(0))
!        end do
!        xd(nd)=x
!        yd(nd)=y
        xd(nd)=sqrx*rand(0)
        yd(nd)=sqry*rand(0)
!        if (xd(nd)<2*rd) then
!            limitDomain(nd)=1
!        else if (yd(nd)<2*rd) then
!            limitDomain(nd)=2
!        else if (sqrx-xd(nd)<2*rd) then
!            limitDomain(nd)=3
!        else if (sqry-yd(nd)<2*rd) then
!            limitDomain(nd)=4
!        else
!            limitDomain(nd)=0
!        end if
        limitDomain(nd,:)=[0,0]!value=100 for not contact with limit
        if (xd(nd)<2*rd) then
            limitDomain(nd,1)=-1
        else if (sqrx-xd(nd)<2*rd) then
            limitDomain(nd,1)=1
        end if
        if (yd(nd)<2*rd) then
            limitDomain(nd,2)=-1
        else if (sqry-yd(nd)<2*rd) then
            limitDomain(nd,2)=1
        end if
        !limitDomain(nd) = (xd(nd)<2*rd .or. sqrx-xd(nd)<2*rd .or. yd(nd)<2*rd .or. sqry-yd(nd)<2*rd)
        foundCluster=.false.
        ic=0
        do while((.not. foundCluster) .and. ic<icmax)
            !write(*,*) nd,ic
            ic=ic+1
!            write(*,*) ic,nbCH(ic)
            if (nbCH(ic)/=0) then
                if (nbCH(ic)==1) then
                    goingThrough=[0,0]
                    exterior=.false.
                    iod=firstCH(ic)
                    exterior = ((xd(nd)-xd(iod))**2+(yd(nd)-yd(iod))**2>4*rd2)
                    if (exterior .and. any(limitDomain(nd,:)/=[0,0]) .and. any(limitDomain(nd,:)==-limitDomain(iod,:))) then
                        call periodicContact(xd,yd,rd2,nd,iod,sqrx,sqry,dimnd,goingThrough)
                        if (any(goingThrough/=[0,0])) exterior=.false.
                    end if
                    if (.not. exterior) then
                        !write(*,*) 'test'
                        foundCluster=.true.
                        nbCH(ic)=nbCH(ic)+1
                        nextCH(lastCH(ic))=nd
                        nextCH(nd)=firstCH(ic)

                        lastCH(ic)=nd
                        if (any(goingThrough/=[0,0])) then
                            nextThrough(nd,:)=goingThrough
                            nextThrough(firstCH(ic),:)=-goingThrough
                        end if
                    end if
                else if (isStrip(ic)) then
                    if (orientationStrip(ic,2)>=orientationStrip(ic,1)) then
                        if (orientationStrip(ic,1)==0) then
                            slope=0.
                        else
                            slope=sign(1,orientationStrip(ic,1))/float(orientationStrip(ic,2))
                        end if
                        proj=orientationStrip(ic,2)/sqrt(float(orientationStrip(ic,2)**2+orientationStrip(ic,1)**2))
                        relativePosition=mod(xd(nd)-yd(nd)*slope-coordinateStrip(ic)&
                        +2*rd*proj+2*sqrx,sqrx/float(orientationStrip(ic,2)))
                        contactStrip=(relativePosition<=thicknessStrip(ic)+4*rd*proj)
                        if (contactStrip) then
                            relativePosition=mod(xd(nd)-yd(nd)*slope-coordinateStrip(ic)+2*sqrx,sqrx/float(orientationStrip(ic,2)))
                            isInsideStrip=(relativePosition<=thicknessStrip(ic))
                        end if
                    else
                        if (orientationStrip(ic,2)==0) then
                            slope=0
                        else
                            slope=sign(1,orientationStrip(ic,2))/float(orientationStrip(ic,1))
                        end if
                        proj=orientationStrip(ic,1)/sqrt(float(orientationStrip(ic,2)**2+orientationStrip(ic,1)**2))
                        relativePosition=mod(yd(nd)-xd(nd)*slope-coordinateStrip(ic)&
                        +2*rd*proj+2*sqry,sqry/float(orientationStrip(ic,1)))
                        contactStrip=(relativePosition<=thicknessStrip(ic)+4*rd*proj)
                        if (contactStrip) then
                            relativePosition=mod(yd(nd)-xd(nd)*slope-coordinateStrip(ic)+2*sqry,sqry/float(orientationStrip(ic,1)))
                            isInsideStrip=(relativePosition<=thicknessStrip(ic))
                        end if
                    end if
                    if (contactStrip) then
                        foundCluster=.true.
                        if (.not. isInsideStrip) then
                            if (relativePosition>thicknessStrip(ic)) then
                                thicknessStrip(ic)=relativePosition
                            else
                                if (relativePosition>=0) then
                                    write(*,*) "error, disc is inside band and shouldn't"
                                    call exit(0)
                                else
                                    coordinateStrip(ic)=coordinateStrip(ic)+relativePosition
                                end if
                            end if
                        end if
                    end if
                else

                    !!! this bloc determines the number and positions of the cluster with witch to verify contact of the disc
                    idch=firstCH(ic)
                    icount=1
                    wherePeriodic=[0,0]
                    do i=1,10
                        whereImages(i,:)=[0,0]
                    end do
                    nbImages=1

                    do while(icount<=nbCH(ic))
                        if (any(limitDomain(nd,:)/=[0,0]) .and. any(limitDomain(nd,:)==-limitDomain(idch,:))) then
                            alreadyIn=.false.
                            do i=1,nbImages
                                if (all(wherePeriodic+limitDomain(idch,:)==whereImages(i,:))) alreadyIn=.true.
                            end do
                            if (.not. alreadyIn) then
                                nbImages=nbImages+1
                                whereImages(nbImages,:)=wherePeriodic+limitDomain(idch,:)
                            end if
                        end if

                        if (any(nextThrough(idch,:)/=[0,0])) then
                            if (nextThrough(idch,1)==0 .or. nextThrough(idch,2)==0) then
                                wherePeriodicTemp=wherePeriodic+nextThrough(idch,:)
                                alreadyIn=.false.
                                do i=1,nbImages
                                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                                end do
                                if (.not. alreadyIn) then
                                    nbImages=nbImages+1
                                    whereImages(nbImages,:)=wherePeriodicTemp
                                end if
                            else
                                wherePeriodicTemp=wherePeriodic+[nextThrough(idch,1),0]
                                alreadyIn=.false.
                                do i=1,nbImages
                                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                                end do
                                if (.not. alreadyIn) then
                                    nbImages=nbImages+1
                                    whereImages(nbImages,:)=wherePeriodicTemp
                                end if

                                wherePeriodicTemp=wherePeriodic+[0,nextThrough(idch,2)]
                                alreadyIn=.false.
                                do i=1,nbImages
                                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                                end do
                                if (.not. alreadyIn) then
                                    nbImages=nbImages+1
                                    whereImages(nbImages,:)=wherePeriodicTemp
                                end if

                                wherePeriodicTemp=wherePeriodic+nextThrough(idch,:)
                                alreadyIn=.false.
                                do i=1,nbImages
                                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                                end do
                                if (.not. alreadyIn) then
                                    nbImages=nbImages+1
                                    whereImages(nbImages,:)=wherePeriodicTemp
                                end if
                            end if
                        end if
                        wherePeriodic=wherePeriodic+nextThrough(idch,:)
                        idch=nextCH(idch)
                        icount=icount+1
                    end do
                    !write(*,*) whereImages
                    !write(*,*) 'nbImages', nbImages

                    !!!
                    icount=1
                    idch=firstCH(ic)
                    !possibleContactThrough=.false.
                    do while(icount<=nbCH(ic))
                        if (any(limitDomain(nd,:)==-limitDomain(idch,:))) then
                            alreadyIn=.false.
                            do i=1,nbImages
                                if (all(limitDomain(idch,:)==whereImages(i,:))) alreadyIn=.true.
                            end do
                            if (.not. alreadyIn) then
                                nbImages=nbImages+1
                                whereImages(nbImages,:)=limitDomain(idch,:)
                            end if
                        end if
                        icount=icount+1
                        idch=nextCH(idch)
                    end do
                    !!!
                    do i=1,10
                        exteriors(i)=.false.
                        extLimits(i)=.false.
                    end do

                    do i=1,nbImages
                        idch=firstCH(ic)
                        icount=1
                        afterPrevious=.false.
                        wherePeriodic=-whereImages(i,:)!!!minus sign very important
                        do while((.not. exteriors(i)) .and. icount<=nbCH(ic)+1)
                            icount=icount+1
                            !write(*,*) wherePeriodic
                            wherePeriodicNext=wherePeriodic+nextThrough(idch,:)
                            call inOrOutLimit(xd(idch)+wherePeriodic(1)*sqrx,yd(idch)&
                            +wherePeriodic(2)*sqry,xd(nextCH(idch))+wherePeriodicNext(1)*sqrx,yd(nextCH(idch))&
                            +wherePeriodicNext(2)*sqry, &
                            xd(nd),yd(nd),rd,exteriors(i),afterPrevious,extLimits(i))
                            idch=nextCH(idch)
                            wherePeriodic=wherePeriodicNext
                        end do
                    end do
                    exterior=.true.
                    imagesContact=[0,0,0]
                    do i=1,nbImages
                        if (.not. exteriors(i)) then
                            if (exterior) then
                                imagesContact(1)=i
                                exterior=.false.
                            else if (imagesContact(2)==0) then
                                imagesContact(2)=i
                            else if (imagesContact(3)==0) then
                                imagesContact(3)=i
                                !write(*,*) 'koala'
!                                if (nd==29020) then
!                                    write(*,*) 'nbCH(ic) = ',nbCH(ic)
!                                    idch=firstCH(ic)
!                                    do
!                                        write(389,*) xd(idch),yd(idch)
!                                        idch=nextCH(idch)
!                                        if (idch==firstCH(ic)) exit
!                                    end do
!                                    write(389,*) ' '
!                                end if
                                write(*,*) 'contact with 3 periodic images of the same cluster : TOTAL DETACHMENT'
                                call exit(0)
                            else
                                write(*,*) 'contact with more than 3 images of the same cluster : TOTAL DETACHMENT'
                                call exit(0)
                            end if
                        end if
                    end do

                    if (.not. exterior) then
!                        if (nd==324) then
!                            !write(*,*) ic,nbCH(ic)
!                            idch=firstCh(ic)
!                            do i=1,nbCH(ic)
!                                write(*,*) nextThrough(idch,:)
!                                idch=nextCH(idch)
!                            end do
!                        end if
                        foundCluster=.true.
                        if (extLimits(imagesContact(1))) then
                            idch=firstCH(ic)
                            wherePeriodic=-whereImages(imagesContact(1),:)
                            wherePeriodicNext=wherePeriodic+nextThrough(idch,:)
                            angle=0.
                            minAngle=angle; maxAngle=angle
                            iMinAngle=idch; iMaxAngle=idch
                            whereMinAngle=wherePeriodic;whereMaxAngle=wherePeriodic
                            ux=xd(firstCH(ic))+wherePeriodic(1)*sqrx-xd(nd)
                            uy=yd(firstCH(ic))+wherePeriodic(2)*sqry-yd(nd)
                            idch=nextCH(idch)
                            wherePeriodic=wherePeriodicNext
                            whereperiodicNext=wherePeriodicNext+nextThrough(idch,:)

                            do while(idch/=firstCH(ic))
                                angle=angleBetweenTwoVectors(ux,uy,xd(idch)+wherePeriodic(1)&
                                *sqrx-xd(nd),yd(idch)+wherePeriodic(2)*sqry-yd(nd))
                                if (angle<=minAngle) then
                                    minAngle=angle
                                    iMinAngle=idch
                                    whereMinAngle=wherePeriodic
                                    !write(*,*) 'tsoin tsoin', whereMinAngle,wherePeriodic
                                else if (angle>maxAngle) then
                                    maxAngle=angle
                                    iMaxAngle=idch
                                    whereMaxAngle=wherePeriodic
                                end if
                                idch=nextCH(idch)
                                wherePeriodic=wherePeriodicNext
                                wherePeriodicNext=wherePeriodicNext+nextThrough(idch,:)

                            end do
                            vect_product=(xd(nd)-xd(iMinAngle)-whereMinAngle(1)*sqrx)*(yd(iMaxAngle)-whereMaxAngle(2)*sqry-yd(nd))&
                            -(yd(nd)-yd(iMinAngle)-whereMinAngle(2)*sqry)* &
                            (xd(iMaxAngle)-whereMaxAngle(1)*sqrx-xd(nd))
                            if (vect_product>0.) then
                                nextCH(iMinAngle)=nd
                                nextCH(nd)=iMaxAngle
                                lastCH(ic)=iMinAngle
                                firstCH(ic)=nd
                                nextThrough(nd,:)=whereMaxAngle
                                nextThrough(iMinAngle,:)=-whereMinAngle
                                !write(*,*) 'cucu',whereMinAngle
!                                if (nd==15898) then
!                                    write(*,*) 'signal a'
!                                    write(*,*) nextThrough(nd,:)
!                                    write(*,*) nextThrough(iMinAngle,:)
!                                end if
                            else
                                nextCH(iMaxAngle)=nd
                                nextCH(nd)=iMinAngle
                                lastCH(ic)=iMaxAngle
                                firstCH(ic)=nd
                                nextThrough(nd,:)=whereMinAngle
                                !write(*,*) 'cucu2',nextThrough(nd,:)
                                nextThrough(iMaxAngle,:)=-whereMaxAngle
!                                if (nd==15898) then
!                                    write(*,*) 'signal b'
!                                    write(*,*) nextThrough(nd,:)
!                                    write(*,*) whereMinAngle
!                                    write(*,*) whereImages(imagesContact(1),:)
!                                    !write(*,*) nextThrough(iMaxAngle,:)
!                                end if
                            end if
                            icount=1
                            idch=nextCH(nd)
                            do while(idch/=nd)
                                icount=icount+1
                                idch=nextCH(idch)
                            end do
                            nbCH(ic)=icount
!                            if (nd==324) then
!                            write(*,*) ic,nbCH(ic)
!                            idch=firstCh(ic)
!                            do i=1,nbCH(ic)
!                                write(*,*) nextThrough(idch,:)
!                                idch=nextCH(idch)
!                            end do
!                        end if

                            if (imagesContact(2)/=0) then
                                !programmer l'auto convexification : cluster en contact avec lui même
                                if (imagesContact(3)/=0) then
                                    write(*,*) '3 periodic images of the same cluster in contact : total detachment'
                                    call exit(0)
                                    !!!cluster invades whole domain
                                else
                                    isStrip(ic)=.true.
                                    print*, 'whereImages',nbImages
                                    print*, whereImages(:,:)
                                    throughFromEnd2End=whereImages(imagesContact(2),:)-whereImages(imagesContact(1),:)
                                    if (throughFromEnd2End(2)>=throughFromEnd2End(1)) then
                                        mainOrientation=2
                                    else
                                        mainOrientation=1
                                    end if

                                    if (throughFromEnd2End(mainOrientation)<0) throughFromEnd2End=-throughfromEnd2End
                                    orientationStrip(ic,:)=throughFromEnd2End
                                    if (all(throughFromEnd2End==[0,0])) then
                                        !!not supposed to happen
                                        !write(*,*) 'test'
                                        write(*,*) 'error : strip with no orientation : aborting program'
                                        call exit(0)
                                    else if (throughFromEnd2End(1)==0) then
                                        if(throughFromEnd2End(2)>1) then
                                            write(*,*) 'erreur : cluster 2 fois enroulé sans inclinaison !'
                                            call exit(0)
                                        end if
                                        idch=firstCH(ic)
                                        minStrip=xd(idch)
                                        maxStrip=xd(idch)
                                        wherePeriodic=nextThrough(idch,:)
                                        idch=nextCH(idch)
                                        do while (idch/=firstCH(ic))
                                            coordinate=xd(idch)+wherePeriodic(1)*sqrx
                                            if (coordinate<minStrip) then
                                                minStrip=coordinate
                                            else if (coordinate>maxStrip) then
                                                maxStrip=coordinate
                                            end if
                                            wherePeriodic=wherePeriodic+nextThrough(idch,:)
                                            idch=nextCH(idch)
                                        end do
                                        coordinateStrip(ic)=minStrip
                                        thicknessStrip(ic)=maxStrip-minStrip
                                    else if (throughFromEnd2End(2)==0) then
                                        if(throughFromEnd2End(1)>1) then
                                            write(*,*) 'erreur : cluster 2 fois enroulé sans inclinaison !'
                                            call exit(0)
                                        end if
                                        idch=firstCH(ic)
                                        minStrip=yd(idch)
                                        maxStrip=yd(idch)
                                        wherePeriodic=nextThrough(idch,:)
                                        idch=nextCH(idch)
                                        do while (idch/=firstCH(ic))
                                            coordinate=yd(idch)+wherePeriodic(2)*sqry
                                            if (coordinate<minStrip) then
                                                minStrip=coordinate
                                            else if (coordinate>maxStrip) then
                                                maxStrip=coordinate
                                            end if
                                            wherePeriodic=wherePeriodic+nextThrough(idch,:)
                                            idch=nextCH(idch)
                                        end do
                                        coordinateStrip(ic)=minStrip
                                        thicknessStrip(ic)=maxStrip-minStrip
                                    else if (throughFromEnd2End(1)>1 .and. throughFromEnd2End(2)>1) then
                                        !!cluster invades whole domain : to discuss
                                    else




                                        !! to calculate : 2 discs borders of the strip
                                        !! "orthogonal distance" between them
                                        !!! this bloc calculates the two discs at the border of the strip
                                        idch=firstCH(ic)
                                        orthogonalDistanceMax=0.
                                        orthogonalDistanceMin=0.
                                        iMinOrthDistance=idch
                                        iMaxOrthDistance=idch
                                        wherePeriodic=nextThrough(idch,:)
                                        idch=nextCH(idch)
                                        normOrientation=sqrt(float(throughFromEnd2End(1)**2+throughFromEnd2End(2)**2))
                                        do while(idch/=firstCH(ic))
                                            orthogonalDistance=(throughFromEnd2End(1)*&
                                            (yd(idch)+wherePeriodic(2)*sqry-yd(firstCH(ic)))&
                                            -throughFromEnd2End(2)*(xd(idch)+wherePeriodic(1)*sqrx-xd(firstCH(ic))))/normOrientation
                                            if (orthogonalDistance>orthogonalDistanceMax) then
                                                orthogonalDistanceMax=orthogonalDistance
                                                iMaxOrthDistance=idch
                                            else if (orthogonalDistance<orthogonalDistanceMin) then
                                                orthogonalDistanceMin=orthogonalDistance
                                                iMinOrthDistance=idch
                                            end if
                                            wherePeriodic=wherePeriodic+nextThrough(idch,:)
                                            idch=nextCH(idch)
                                        end do

                                        if (mainOrientation==2) then
                                            coordinateStrip(ic)=mod(xd(iMaxOrthDistance)&
                                            -float(throughfromEnd2End(1))/throughFromEnd2End(2)*yd(iMaxOrthDistance),sqrx)
                                            thicknessStrip(ic)=(orthogonalDistanceMax-orthogonalDistanceMin)&
                                            *normOrientation/throughFromEnd2End(2)
                                        else
                                            coordinateStrip(ic)=mod(yd(iMaxOrthDistance)&
                                            -float(throughfromEnd2End(2))/throughFromEnd2End(1)*xd(iMaxOrthDistance),sqry)
                                            thicknessStrip(ic)=(orthogonalDistanceMax-orthogonalDistanceMin)&
                                            *normOrientation/throughFromEnd2End(1)
                                        end if

                                    end if

                                end if
                            end if
                        end if

!                        if (nbCH(ic)==1) then
!                            nbCH(ic)=nbCH(ic)+1
!                            nextCH(lastCH(ic))=nd
!                            nextCH(nd)=firstCH(ic)
!                            lastCH(ic)=nd
!                            if (goingThrough) then
!                                nextThrough(nd,:)=goingThrough
!                                nextThrough(firstCH(ic),:)=-goingThrough
!                            end if
!                        else
!                            call convexHull(xd,yd,firstCH(ic),nextCH,ndmax,nbInThisCH,listDiscsCH)
!                            nbCH(ic)=nbInThisCH
!                            firstCH(ic)=listDiscsCH(1)
!                            lastCH(ic)=listDiscsCH(nbInThisCH)
!                            do i=1,nbInThisCH
!                                nextCH(listDiscsCH(i))=listDiscsCH(i+1)
!                            end do
!                            nextCH(listDiscsCH(nbInThisCH))=listDiscsCH(1)







                    end if
                end if
            end if
        end do

        if (foundCluster) then
            keepOn=.true.
            ifile=0
            do while(keepOn)
                bool=.false.
                idch=firstCH(ic)
                do
                    if (any(nextThrough(idch,:)/=[0,0])) bool=.true.
                    idch=nextCH(idch)
                    if (idch==firstCH(ic)) exit
                end do

                if (isStrip(ic)) then
                    if (orientationStrip(ic,2)>=orientationStrip(ic,1)) then
                        i=1
                    else
                        i=2
                    end if
                    if (thicknessStrip(ic)+2*rd>sqr(i)/orientationStrip(ic,3-i)) then
                        print*, 'strip in contact with itself, invades the whole domain'
                        call exit(0)
                    end if
                    ioc=0
                    foundCluster=.false.
                    do while(ioc<icmax .and. (.not. foundCluster))
                        ioc=ioc+1
                        if (ioc/=ic) then
                            if (isStrip(ioc)) then
                                if (any(orientationStrip(ic,:)/=orientationStrip(ioc,:))) then
                                    !!!!event : decollement total
                                else
                                    !!Here : specific procedure to test contact between and merge two strips
                                    if (orientationStrip(ic,2)>=orientationStrip(ic,1)) then
                                        i=2
                                    else
                                        i=1
                                    end if
                                    period=orientationStrip(ic,i)
                                    coordinateIC=mod(coordinateStrip(ic)+sqr(3-i),sqr(3-i)/period)
                                    coordinateIOC=mod(coordinateStrip(ioc)+sqr(3-i),sqr(3-i)/period)
                                    if (coordinateIC<coordinateIOC) then
                                        if (coordinateIOC<coordinateIC+thicknessStrip(ic)) then
                                            if (coordinateIOC+thicknessStrip(ioc)<coordinateIC+thicknessStrip(ic)) then
                                                foundCluster=.true.
                                                IOCinsideIC=.true.
                                            else
                                                if (mod(coordinateIOC+thicknessStrip(ioc),sqr(3-i)/period)>coordinateIC) then
                                                    write(*,*) 'union of two strips cover the whole domain : TOTAL DETACHMENT'
                                                    call exit(0)
                                                    !!!! union of two strips cover the whole domain -----> total detachment
                                                else
                                                    foundCluster=.true.
                                                    contactICbeforeIOC=.true.
                                                    newThickness=coordinateIOC-coordinateIC+thicknessStrip(ioc)

                                                end if
                                            end if
                                        else
                                            if (mod(coordinateIOC+thicknessStrip(ioc),sqr(3-i)/period)>coordinateIC) then
                                                foundCluster=.true.
                                                contactIOCbeforeIC=.true.
                                                newThickness=mod(coordinateIC+thicknessStrip(ic)&
                                                -coordinateIOC+sqr(3-i),sqr(3-i)/period)
                                            end if
                                        end if
                                    else
                                        if (coordinateIC<coordinateIOC+thicknessStrip(ioc)) then
                                            if (coordinateIC+thicknessStrip(ic)<coordinateIOC+thicknessStrip(ioc)) then
                                                foundCluster=.true.
                                                ICinsideIOC=.true.
                                            else
                                                if (mod(coordinateIC+thicknessStrip(ic),sqr(3-i)/period)>coordinateIOC) then
                                                    !!!! union of two strips cover the whole domain -----> total detachment
                                                    write(*,*) 'union of two strips cover the whole domain : TOTAL DETACHMENT'
                                                    call exit(0)
                                                else
                                                    foundCluster=.true.
                                                    contactIOCbeforeIC=.true.
                                                    newThickness=coordinateIC-coordinateIOC+thicknessStrip(ic)
                                                end if
                                            end if
                                        else
                                            if (mod(coordinateIC+thicknessStrip(ic),sqr(3-i)/period)>coordinateIOC) then
                                                foundCluster=.true.
                                                contactICbeforeIOC=.true.
                                                newThickness=mod(coordinateIOC+thicknessStrip(iOc)&
                                                -coordinateIC+sqr(3-i),sqr(3-i)/period)
                                            end if
                                        end if
                                    end if
                                end if

                            else if (nbCH(ioc)/=0) then
                                !!test contact strip+cluster
                                idoch=firstCH(ioc)
    !                            iWhichStrip=0
                                whichStripContactIC=-1
                                do i=1,nbCH(ioc)
                                    call contactWhereStrip(sqrx,sqry,xd(idch),yd(idoch),rd,orientationStrip(ic,:),&
                                    coordinateStrip(ic),thicknessStrip(ic),contactOrNotStrip,whichStripContact)
                                    if (contactOrNotStrip) then
                                        if (contact) then
                                            if (whichStripContact/=whichStripContactIC) then
                                                !!! contact with 2 strips of the same cluster : cluster invades whole domain
                                            end if
                                        else
                                            whichStripContactIC=whichStripContact
                                            discContactStrip=idoch
                                        end if
                                        foundCluster=.true.
    !                                    do i =1,iWhichStrip
    !                                        if (whichStripContact=whichStripContactTab(i)) then
    !                                            alreadyIn=.true.
    !                                        end if
    !                                    end do
    !                                    if (.not. alreadyIn) then
    !                                        iWhichStrip=iWhichStrip+1
    !                                        whichStripContactTab(iWhichStrip)=whichStripContact
    !                                    end if
                                    end if
                                    idoch=nextCH(idoch)
                                end do
                            end if
                        end if

                    end do
                    if (foundCluster) then
                        if (isStrip(ioc)) then
                            if (IOCinsideIC) then
                                isStrip(ioc)=.false.
                                nbCH(ioc)=0
                            else if (ICinsideIOC) then
                                coordinateStrip(ic)=coordinateStrip(ioc)
                                thicknessStrip(ic)=thicknessStrip(ioc)
                                isStrip(ioc)=.false.
                                nbCH(ioc)=0
                            else if (contactICbeforeIOC) then
                                thicknessStrip(ic)=newThickness
                                isStrip(ioc)=.false.
                                nbCH(ioc)=0
                            else if (contactIOCbeforeIC) then
                                thicknessStrip(ic)=newThickness
                                coordinatestrip(ic)=coordinateStrip(ioc)
                                isStrip(ioc)=.false.
                                nbCH(ioc)=0
                            else
                                write(*,*) '2 strips merging : case not recognized, aborting program'
                                call exit(0)
                            end if
                        else
                            !!merging of strip ic with cluster ioc
                            call twoPointsStrip(coordinateStrip(ic),thicknessStrip(ic),orientationStrip(ic,:),&
                            whichStripContactIC,sqrx,sqry,x1,y1,x2,y2)
                            xd(ndmax+1)=x1;yd(ndmax+1)=y1;xd(ndmax+2)=x2;yd(ndmax+2)=y2
                            nextCH(ndmax+2)=nextCH(discContactStrip)
                            nextCH(discContactStrip)=ndmax+1
                            nextCH(ndmax+1)=ndmax+2
                            nextThrough(ndmax+2,:)=nextThrough(discContactStrip,:)
                            nextThrough(discContactStrip,:)=[0,0]
                            nextThrough(ndmax+1,:)=[0,0]
                            !write(*,*) 'test'
                            !write(*,*) 'makeStip strip2Cluster'
                            call makeStrip(xd,yd,dimnd,orientationStrip(ic,:),firstCH(ioc),&
                            nextThrough,nextCH,sqrx,sqry,coordinateStripOut,thicknessStripOut)
                            !write(*,*) 'fin test'
                            !isStrip(ic)=.true.
                            !orientationStrip(ic,:)=orientationStrip(ioc,:)
                            coordinateStrip(ic)=coordinateStripOut
                            thicknessStrip(ic)=thicknessStripOut
                            !isStrip(ioc)=.false.
                            nbCH(ioc)=0
                        end if
                    else
                        keepOn=.false.
                    end if
                else
                    autoconvex=.false.

                    !stripInList=.false.
                    call describeImages(firstCH(ic),dimnd,nbCH(ic),limitDomain,nextCH,nextThrough,whereImages,nbImages)
!                    write(*,*) whereImages
!                    write(*,*) 'nbImages',nbImages
                    ioc=0
                    foundCluster=.false.
                    do while(ioc<icmax .and. (.not. foundCluster))!ioc=1,icmax
                        ioc=ioc+1
                        if (isStrip(ioc)) then
                            if (ioc/=ic) then
                                contact=.false.
                                idch=firstCH(ic)
    !                            iWhichStrip=0
                                whichStripContactIC=-1
                                do i=1,nbCH(ic)
                                    call contactWhereStrip(sqrx,sqry,xd(idch),yd(idch),rd,orientationStrip(ioc,:),&
                                    coordinateStrip(ioc),thicknessStrip(ioc),contactOrNotStrip,whichStripContact)
                                    if (contactOrNotStrip) then
                                        if (contact) then
                                            if (whichStripContact/=whichStripContactIC) then
                                                !!! contact with 2 strips of the same cluster : cluster invades whole domain
                                            end if
                                        else
                                            whichStripContactIC=whichStripContact
                                            discContactStrip=idch
                                        end if
                                        contact=.true.
    !                                    do i =1,iWhichStrip
    !                                        if (whichStripContact=whichStripContactTab(i)) then
    !                                            alreadyIn=.true.
    !                                        end if
    !                                    end do
    !                                    if (.not. alreadyIn) then
    !                                        iWhichStrip=iWhichStrip+1
    !                                        whichStripContactTab(iWhichStrip)=whichStripContact
    !                                    end if
                                    end if
                                    idch=nextCH(ic)
                                end do
                                if (contact) then
    !                                nbOfClusToMerge=nbOfClusToMerge+1
    !                                listOfClusToMerge(nbOfClusToMerge)=ioc
    !                                stripInList=.true.
                                    foundCluster=.true.
                                end if
                            end if


                        else if (ioc==ic) then
                            if (nbImages>1) then
                                contact=.false.
                                nbImagesContact=0
                                do i=2,nbImages
                                    !! firstly, test contacts of discs of image 1 with other images
                                    contactThisImage=.false.
                                    idoch=firstCH(ioc)
                                    wherePeriodicIOC=-whereImages(1,:)
                                    do while ((.not. contactThisImage) .and. iBigLoop<=nbCH(ioc))
                                        exterior=.false.
                                        idch=firstCH(ic)
                                        icount=1
                                        afterPrevious=.false.
                                        wherePeriodic=-whereImages(i,:)
                                        do while ((.not. exterior).and. icount <=nbCH(ic)+1)
                                            icount = icount+1
                                            wherePeriodicNext=wherePeriodic+nextThrough(idch,:)
                                            call inOrOut(xd(idch)+wherePeriodic(1)*sqrx,yd(idch)&
                                            +wherePeriodic(2)*sqry,xd(nextCH(idch))+wherePeriodicNext(1)*sqrx,&
                                            yd(nextCH(idch))+wherePeriodicNext(2)*sqry,&
                                            xd(idoch)+wherePeriodicIOC(1)*sqrx,yd(idoch)+&
                                            wherePeriodicIOC(2)*sqry,rd,exterior,afterPrevious)
                                            wherePeriodic=wherePeriodicNext
                                            idch=nextCH(idch)
                                        end do
                                        if (.not. exterior) then
                                            contactThisImage=.true.
                                        end if
                                        iBigLoop=iBigLoop+1
                                        wherePeriodicIOC=wherePeriodicIOC+nextThrough(idoch,:)
                                        idoch=nextCH(idoch)
                                    end do
                                    !!!
                                    !!!then, test contacts of discs of other images with image 1 !!!
                                    idch=firstCH(ic)
                                    wherePeriodic=-whereImages(i,:)
                                    do while ((.not. contactThisImage) .and. iBigLoop<=nbCH(ic))
                                        exterior=.false.
                                        idoch=firstCH(ioc)
                                        icount=1
                                        afterPrevious=.false.
                                        wherePeriodicIOC=-whereImages(1,:)
                                        do while ((.not. exterior).and. icount <=nbCH(ioc)+1)
                                            icount = icount+1
                                            wherePeriodicIOCNext=wherePeriodicIOC+nextThrough(idoch,:)
                                            call inOrOut(xd(idoch)+wherePeriodicIOC(1)*sqrx,yd(idoch)&
                                            +wherePeriodicIOC(2)*sqry,xd(nextCH(idoch))+wherePeriodicIOCNext(1)*sqrx,&
                                            yd(nextCH(idoch))+wherePeriodicIOCNext(2)*sqry,&
                                            xd(idch)+wherePeriodic(1)*sqrx,yd(idch)+&
                                            wherePeriodic(2)*sqry,rd,exterior,afterPrevious)
                                            wherePeriodicIOC=wherePeriodicIOCNext
                                            idoch=nextCH(idoch)
                                        end do
                                        if (.not. exterior) then
                                            contactThisImage=.true.
                                        end if
                                        iBigLoop=iBigLoop+1
                                        wherePeriodic=whereperiodic+nextThrough(idch,:)
                                        idch=nextCH(idch)
                                    end do
                                    if (contactThisImage) then
                                        nbImagesContact=nbImagesContact+1
                                        throughFromEnd2End=whereImages(i,:)-whereImages(1,:)
                                    end if
                                end do
                                if (nbImagesContact==1) then
                                    autoconvex=.true.
                                    foundCluster=.true.
                                else if (nbImagesContact>1) then
!                                    if (nd==15984) then
!                                        write(*,*) nbCH(ic)
!                                        idch=firstCH(ic)
!                                        do
!                                            write(175,*) xd(idch),yd(idch)
!                                            idch=nextCH(idch)
!                                            if (idch==firstCH(ic)) exit
!                                        end do
!                                        write(175,*) ' '
!                                    end if
!                                    write(*,*) 'koala'
                                    write(*,*) '3 periodic images of the same cluster in contact : TOTAL DETACHMENT'
                                    call exit(0)
                                    !!!!!!event : total detachment
                                end if
                            end if

                        else
                            if (nbCH(ioc)/=0) then
                                idoch=firstCH(ioc)
                                !print*, 'idoch',idochAltern
                                iBigLoop=1
                                contact=.false.
                                do while ((.not. contact) .and. iBigLoop<=nbCH(ioc))
                                    do i=1,10
                                        exteriors(i)=.false.
                                    end do
                                    !write(*,*) 'test gamma',nbImages
                                    do i=1,nbImages
                                        idch=firstCH(ic)
                                        icount=1
                                        afterPrevious=.false.
                                        wherePeriodic=-whereImages(i,:)!!!minus sign very important
                                        do while((.not. exteriors(i)) .and. icount<=nbCH(ic)+1)
                                            icount=icount+1
                                            !write(*,*) wherePeriodic
                                            wherePeriodicNext=wherePeriodic+nextThrough(idch,:)
                                            call inOrOut(xd(idch)+wherePeriodic(1)*sqrx,yd(idch)&
                                            +wherePeriodic(2)*sqry,xd(nextCH(idch))+wherePeriodicNext(1)*sqrx,&
                                            yd(nextCH(idch))+wherePeriodicNext(2)*sqry,&
                                            xd(idoch),yd(idoch),rd,exteriors(i),afterPrevious)
                                            idch=nextCH(idch)
                                            wherePeriodic=wherePeriodicNext
                                        end do
                                        if (.not. exteriors(i)) then
                                            iContact=i
                                            discIOCContact=idoch
                                            !write(*,*) 'idoch',idochAltern
                                            !write(*,*) 'discIOCCOntact',discIOCContact
                                        end if
                                    end do
                                    !write(*,*) 'exteriors',exteriors(1:nbImages)
    !                                exterior=.false.
    !                                idch=firstCH(ic)
    !                                icount=1
    !                                idchPrevious=lastCH(ic)
    !                                afterPrevious=.false.
    !                                do while((.not. exterior) .and. icount<=nbCH(ic)+1)
    !                                    icount=icount+1
    !                                    call inOrOut(xd(idch),yd(idch),xd(nextCH(idch)),yd(nextCH(idch)) &
    !                                    ,xd(idoch),yd(idoch),rd,exterior,afterPrevious)
    !                                    idchPrevious=idch
    !                                    idch=nextCH(idch)
    !                                end do
                                    if (.not. all(exteriors(1:nbImages))) then
                                        contact=.true.
                                        !write(*,*) 'test a'
                                    end if

                                    iBigLoop=iBigLoop+1
                                    idoch=nextCH(idoch)
                                end do
                                if (contact) then
    !                                nbOfClusToMerge=nbOfClusToMerge+1
    !                                listOfClusToMerge(nbOfClusToMerge)=ioc
                                    foundCluster=.true.
                                end if
                            end if

                        end if

                    end do
                    if (foundCluster) then
                        if (isStrip(ioc)) then
                            !write(*,*) 'contact with a strip'
                            call twoPointsStrip(coordinateStrip(ioc),thicknessStrip(ioc),orientationStrip(ioc,:),&
                            whichStripContactIC,sqrx,sqry,x1,y1,x2,y2)
                            xd(ndmax+1)=x1;yd(ndmax+1)=y1;xd(ndmax+2)=x2;yd(ndmax+2)=y2
                            nextCH(ndmax+2)=nextCH(discContactStrip)
                            nextCH(discContactStrip)=ndmax+1
                            nextCH(ndmax+1)=ndmax+2
                            nextThrough(ndmax+2,:)=nextThrough(discContactStrip,:)
                            nextThrough(discContactStrip,:)=[0,0]
                            nextThrough(ndmax+1,:)=[0,0]
                            !write(*,*) 'pangolin'
                            write(*,*) 'makeStrip clusterToStrip'
                            call makeStrip(xd,yd,dimnd,orientationStrip(ioc,:),firstCH(ic),nextThrough,nextCH,sqrx&
                            ,sqry,coordinateStripOut,thicknessStripOut)
                            isStrip(ic)=.true.
                            orientationStrip(ic,:)=orientationStrip(ioc,:)
                            coordinateStrip(ic)=coordinateStripOut
                            thicknessStrip(ic)=thicknessStripOut
                            isStrip(ioc)=.false.
                            nbCH(ioc)=0

                        else if (autoconvex) then
                            !write(*,*) 'autoconvex'
!                            write(*,*) 'bebop'
!                            write(*,*) throughFromEnd2End
                            write(*,*) 'makeStrip autoconvex'
                            call makeStrip(xd,yd,dimnd,throughFromEnd2End,firstCH(ic),nextThrough,nextCH,sqrx,sqry,&
                            coordinateStripOut,thicknessStripOut)
                            isstrip(ic)=.true.
                            orientationStrip(ic,:)=throughFromEnd2End
                            coordinateStrip(ic)=coordinateStripOut
                            thicknessStrip(ic)=thicknessStripOut
                        else


                            nextTemp=nextCH(firstCH(ic))
                            nextCH(firstCH(ic))=nextCH(discIOCContact)
                            nextCH(discIOCContact)=nextTemp
                            !!!!!!especially this following
                            nextThroughTemp=nextThrough(firstCH(ic),:)
                            nextThrough(firstCH(ic),:)=+whereImages(iContact,:)+&
                            nextThrough(discIOCContact,:)
                            nextThrough(discIOCContact,:)=-whereImages(iContact,:)+nextThroughTemp

                            call periodicConvexification(xd,yd,firstCH(ic),nextCH,dimnd,sqrx,sqry,&
                            nextThrough,nbInThisCH,listDiscsCH)

                            nbCH(ic)=nbInThisCH
                            firstCH(ic)=listDiscsCH(1)
                            lastCH(ic)=listDiscsCH(nbInThisCH)
                            do i=1,nbInThisCH
                                nextCH(listDiscsCH(i))=listDiscsCH(i+1)
                            end do
                            nextCH(listDiscsCH(nbInThisCH))=listDiscsCH(1)
                            nbCH(ioc)=0
                        end if

                    else
                        keepOn=.false.
                    end if
    !                if (nbOfClusToMerge>1) then
    !                    do iMerge=2,nbOfClusToMerge
    !                        icHead=ic; icTail=listOfClusToMerge(iMerge)!!!
    !                        idcLastOfHead=lastCH(icHead)                   !!!
    !                        idcFirstofTail=firstCH(icTail)                 !!!
    !                        nextCH(idcLastOfHead)=idcFirstofTail       !!!
    !                        !cluster(idcFirstOfTail)=icHead               !!!
    !                        idcPrev=idcFirstOfTail                       !!!
    !                        idcNext=nextCH(idcFirstOfTail)             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                        nbToAdd=1                                       !!!         Merging of        !!!
    !                        do while(idcNext/=idcFirstOfTail)            !!!          2 Clusters       !!!
    !                            nbToAdd=nbToAdd+1                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            idcPrev=idcNext                          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                            idcNext=nextCH(idcNext)                !!!
    !                            !cluster(idcPrev)=icHead                  !!!
    !                        end do                                       !!!
    !                        nbCH(icHead)=nbCH(icHead)+nbToAdd  !!!
    !                        nbCH(icTail)=0                          !!!
    !                        lastCH(icHead)=lastCH(icTail)
    !                        nextCH(lastCH(icTail))=firstCH(icHead)
    !                    end do
    !                    call convexHull(xd,yd,firstCH(ic),nextCH,ndmax,nbInThisCH,listDiscsCH)
    !                    nbCH(ic)=nbInThisCH
    !                        firstCH(ic)=listDiscsCH(1)
    !                    lastCH(ic)=listDiscsCH(nbInThisCH)
    !                    do i=1,nbInThisCH
    !                        nextCH(listDiscsCH(i))=listDiscsCH(i+1)
    !                    end do
    !                    nextCH(listDiscsCH(nbInThisCH))=listDiscsCH(1)
    !                else
    !                    keepty(fv;On=.false.
    !                end if
                end if
!                if (nd==943) then
!                    print*, foundCluster,isStrip(ic),nbCH(ic),firstCH(ic),nextCH(firstCH(ic)),nextCH(nextCH(firstCH(ic)))
!                end if
            end do
        else
            icmax=icmax+1
            !cluster(nd)=icmax
            firstCH(icmax)=nd
            lastCH(icmax)=nd
            nbCH(icmax)=1
            nextCH(nd)=nd
            nextThrough(nd,:)=[0,0]
        end if


        if (mod(nd,100)==0) then!mod(nd,1)==0) then !.and. nd>10030) then! .or.(nd>=59300 .and. nd<=59400 .and. mod(nd,1)==0)) then! .and. nd>1500 .and. nd<1701) then

            c=rpi*rd**2*nd/(sqrx*sqry)

            area_tot=0
            do ic=1,icmax
                if (nbCH(ic)/=0) then
                    if (isStrip(ic)) then
                        area_tot=area_tot+areaOfaStrip(orientationStrip(ic,:),thicknessStrip(ic),sqrx,sqry,rd)
                    else if (nbCH(ic)==1) then
                        area_tot=area_tot+rpi*rd**2
                    else
                        area_tot=area_tot+areaOfaClusterPeriodic(xd,yd,firstCH(ic),nextCH,nextThrough,rd,dimnd,sqrx,sqry)
                    end if
                end if
            end do
            frac_area=1-area_tot/(sqrx*sqry)
            write(666,*) c,frac_area
            write(*,*) nd,frac_area
            if (frac_area<0.3 .and. (.not. switch)) then
                write(*,*) 'switch'
                switch = .true.
            end if



!            do ix=1,nMaillageTest
!                xi=sqrx*float(ix)/nMaillageTest
!                do jy=1,nMaillageTest
!                    yj=sqry*float(jy)/nMaillageTest
!                    contact=.false.
!                    do ic=1,icmax
!                        exterior=.false.
!                        idch=firstCH(ic)
!                        icount=1
!                        idchPrevious=lastCH(ic)
!                        afterPrevious=.false.
!                        do while((.not. exterior) .and. icount<=nbCH(ic)+1)
!                            icount=icount+1
!                            call inOrOut(xd(idch),yd(idch),xd(nextCH(idch)),yd(nextCH(idch)) &
!                            ,xi,yj,rd,exterior,afterPrevious)
!                            idchPrevious=idch
!                            idch=nextCH(idch)
!                        end do
!                        if (.not.exterior) contact=.true.
!                    end do
!                    if (contact) write(421,*) xi,yj
!                    if (.not. contact) write(422,*) xi,yj
!                end do
!            end do
        end if

    end do






end program
!    do ic=1,icmax
!        do idc=1,ndmax
!            if (cluster==ic) then
!                do iangle=1,30
!                    x=xd(idc)+cos(2*rpi*iangle/30)
!                    y=yd(idc)+sin(2*rpi*iangle/30)
!                    write(13,*) x,y
!                end do
!                write(12,*) xd(idc),yd(idc),cluster(idc)
!            end if
!        end do
!        write(13,*) ' '
!        write (12,*) ' '
!    end do

subroutine inOrOut(xa,ya,xb,yb,xd,yd,rd,outCH,afterPrevious)
    implicit none
    real, intent(in)::xa,ya,xb,yb,xd,yd,rd
    logical, intent(inout)::outCH,afterPrevious
    real::normal_distance,longit_distance,distanceFromDiscToNext
    logical::afterPreviousNext
    afterPreviousNext=.false.
    distanceFromDiscToNext=sqrt((xb-xa)**2+(yb-ya)**2)
    normal_distance=((xb-xa)*(yd-ya)-(yb-ya)*(xd-xa))/distanceFromDiscToNext
    longit_distance=((xb-xa)*(xd-xa)+(yb-ya)*(yd-ya))/distanceFromDiscToNext
    if (normal_distance<-2*rd) then
        outCH=.true.
    elseif ((longit_distance<0.).and.((xd-xa)**2+(yd-ya)**2>4*rd**2).and.(afterPrevious)) then
       outCH=.true.
    end if
    if(longit_distance>distanceFromDiscToNext) then
        afterPreviousNext=.true.
    end if
    afterPrevious=afterPreviousNext

    return
end subroutine

subroutine inOrOutLimit(xa,ya,xb,yb,xd,yd,rd,outCH, &
    afterPrevious,extLimit)
    implicit none
    real, intent(in)::xa,ya,xb,yb,xd,yd,rd
    logical, intent(inout)::outCH,extLimit,afterPrevious
    logical::afterPreviousNext
    real::normal_distance,longit_distance,distanceFromDiscToNext
    afterPreviousNext=.false.
    distanceFromDiscToNext=sqrt((xb-xa)**2+(yb-ya)**2)
    normal_distance=((xb-xa)*(yd-ya)-(yb-ya)*(xd-xa))/distanceFromDiscToNext
    longit_distance=((xb-xa)*(xd-xa)+(yb-ya)*(yd-ya))/distanceFromDiscToNext
    if (normal_distance<-2*rd) then
        outCH=.true.
    elseif ((longit_distance<0.).and.(afterPrevious).and.((xd-xa)**2+(yd-ya)**2>4*rd**2)) then
            outCH=.true.
    end if
    if ((.not. outCH).and.(normal_distance<0)) then
        extLimit=.true.
    end if

    if(longit_distance>distanceFromDiscToNext) then
        afterPreviousNext=.true.
    end if
    afterPrevious=afterPreviousNext

    return
end subroutine


subroutine convexHull(xd,yd,firstCluster,nextDisc,dimnd,nbInThisCH,listDiscsCH)
    implicit none
!    implicit real (a-h,o-z)
    integer,intent(in):: dimnd,firstCluster
    real, dimension(dimnd), intent(in)::xd,yd
    integer, dimension(dimnd), intent(in)::nextDisc
    integer, intent(out)::nbInThisCh
    integer, dimension(dimnd), intent(out)::listDiscsCH
    real::maxCos,cosAngleIter,cosAngle,ux,uy
    integer::iBigLoop,icount,iminy,iloop,imaxCos

    !!find point of minimal y
    iminy=firstCluster
    iloop=firstCluster
    do while (nextDisc(iloop)/=firstCluster)
        iloop=nextDisc(iloop)
        if (yd(iloop)<yd(iminy)) then
            iminy=iloop
        end if
    end do



    ux=1.
    uy=0.
    iBigLoop=iminy
    listDiscsCH(1)=iminy
    icount=1
    do
        imaxCos=firstCluster
        iloop=firstCluster
        maxCos=-10.
        do
            if (iloop/=iBigLoop) then
                cosAngleIter=cosAngle(ux,uy,xd(iloop)-xd(iBigLoop),yd(iloop)-yd(iBigLoop))
                if (cosAngleIter>maxCos) then
                    maxCos=cosAngleIter
                    imaxCos=iloop
                end if
            end if
            iloop=nextDisc(iloop)
            if (iloop==firstCluster) exit
        end do
        if (imaxCos==iminy) exit
        icount=icount+1
        listDiscsCH(icount)=imaxCos
        ux=xd(imaxCos)-xd(iBigLoop)
        uy=yd(imaxCos)-yd(iBigLoop)
        iBigLoop=imaxCos
        if (icount>dimnd) then
            write(*,*) 'error during construction of convex hull : '
            write(*,*) "l'algorithm doesn't find initial exterior disc !!"
            call exit(0)
        end if
    end do
    nbInThisCH=icount


    return
end subroutine

function cosAngle(ux,uy,vx,vy)
    implicit none
    real,intent(in)::ux,uy,vx,vy
    real::cosAngle
    cosAngle=(ux*vx+uy*vy)/sqrt(ux**2+uy**2)/sqrt(vx**2+vy**2)
    return
end function

function angleBetweenTwoVectors(ux,uy,vx,vy)
    implicit none
    real,intent(in)::ux,uy,vx,vy
    real::cosAngle,vect_product
    real::angleBetweenTwoVectors
    cosAngle=(ux*vx+uy*vy)/sqrt(ux**2+uy**2)/sqrt(vx**2+vy**2)
    vect_product=ux*vy-uy*vx
    if (vect_product>0) then
        angleBetweenTwoVectors=acos(cosAngle)
    else
        angleBetweenTwoVectors=-acos(cosAngle)
    end if

    return
end function

function areaOfaClusterPeriodic(xd,yd,first,nextCH,nextThrough,rd,ndmax,sqrx,sqry)
    !!! attention : cette fonction ne marche que pour des clusters implémentés dans le sens géomètrique !!!
    implicit none
    integer,intent(in)::first,ndmax
    real,dimension(ndmax),intent(in)::xd,yd
    integer,dimension(ndmax),intent(in)::nextCH
    integer,dimension(ndmax,2),intent(in)::nextThrough
    real,intent(in)::rd,sqrx,sqry
    real::areaOfaClusterPeriodic,rpi
    integer,dimension(2)::wherePeriodic,wherePeriodicNext
    integer::idch
    rpi=4*atan(1.)
    areaOfaClusterPeriodic=rpi*rd**2
    idch=first
    wherePeriodic=[0,0]
    wherePeriodicNext=nextThrough(first,:)
    do
        areaOfaClusterPeriodic=areaOfaClusterPeriodic+sqrt((xd(nextCH(idch))+wherePeriodicNext(1)&
        *sqrx-xd(idch)-wherePeriodic(1)*sqrx)**2&
        +(yd(nextCH(idch))+wherePeriodicNext(2)*sqry-yd(idch)-wherePeriodic(2)*sqry)**2)*rd &
        +((xd(idch)+wherePeriodic(1)*sqrx)*(yd(nextCH(idch))+wherePeriodicNext(2)*sqry)&
        -(yd(idch)+wherePeriodic(2)*sqry)*(xd(nextCH(idch))+wherePeriodicNext(1)*sqrx))/2.
        idch=nextCH(idch)
        if(idch==first) exit
        wherePeriodic=wherePeriodicNext
        wherePeriodicNext=wherePeriodicNext+nextThrough(idch,:)
    end do
    return

end function

subroutine periodicContact(xd,yd,rd2,nd,iod,sqrx,sqry,ndmax,goingThrough)
    implicit none
    integer,intent(in)::nd,iod,ndmax
    real,intent(in)::sqrx,sqry,rd2
    real,dimension(ndmax),intent(in)::xd,yd
    integer,dimension(2),intent(out)::goingThrough
    if ((xd(nd)+sqrx-xd(iod))**2+(yd(nd)-yd(iod))**2<4*rd2) then
        goingThrough = [-1,0]
    else if ((xd(nd)-sqrx-xd(iod))**2+(yd(nd)-yd(iod))**2<4*rd2) then
        goingThrough = [1,0]
    else if ((xd(nd)-xd(iod))**2+(yd(nd)+sqry-yd(iod))**2<4*rd2) then
        goingThrough = [0,-1]
    else if ((xd(nd)-xd(iod))**2+(yd(nd)-sqry-yd(iod))**2<4*rd2) then
        goingThrough = [0,1]
    else if ((xd(nd)+sqrx-xd(iod))**2+(yd(nd)+sqry-yd(iod))**2<4*rd2) then
        goingThrough = [-1,-1]
    else if ((xd(nd)+sqrx-xd(iod))**2+(yd(nd)-sqry-yd(iod))**2<4*rd2) then
        goingThrough = [-1,1]
    else if ((xd(nd)-sqrx-xd(iod))**2+(yd(nd)+sqry-yd(iod))**2<4*rd2) then
        goingThrough = [1,-1]
    else if ((xd(nd)-sqrx-xd(iod))**2+(yd(nd)-sqry-yd(iod))**2<4*rd2) then
        goingThrough = [1,1]
    else
        goingThrough = [0,0]
    end if

    return

end subroutine

subroutine makeStrip(xd,yd,ndmax,throughFromEnd2End,first,nextThrough,nextCH,sqrx,sqry,coordinateStripOut,thicknessStripOut)
    implicit none
    integer,intent(in)::ndmax,first
    real,intent(in)::sqrx,sqry
    real,dimension(ndmax),intent(in)::xd,yd
    integer,dimension(2),intent(inout)::throughFromEnd2End
    integer,dimension(ndmax,2),intent(in)::nextThrough
    integer,dimension(ndmax)::nextCH
    integer::mainOrientation,idch,iMinOrthDistance,iMaxOrthDistance
    real::minStrip,maxStrip,coordinate,orthogonalDistanceMax,orthogonalDistanceMin,orthogonalDistance,&
    normOrientation
    integer,dimension(2)::wherePeriodic
    real,intent(out)::coordinateStripOut,thicknessStripOut
    !throughFromEnd2End=whereImages(imagesContact(2),:)-whereImages(imagesContact(1),:)
    if (throughFromEnd2End(2)>=throughFromEnd2End(1)) then
        mainOrientation=2
    else
        mainOrientation=1
    end if

    if (throughFromEnd2End(mainOrientation)<0) throughFromEnd2End=-throughfromEnd2End
    !orientationStrip(ic,:)=throughFromEnd2End
    if (all(throughFromEnd2End==[0,0])) then
        !!not supposed to happen
        write(*,*) 'test omicron'
        write(*,*) 'error : strip with no orientation : aborting program'
        call exit(0)
    else if (throughFromEnd2End(1)==0) then
        if(throughFromEnd2End(2)>1) then
            write(*,*) 'erreur : cluster 2 fois enroulé sans inclinaison !'
            call exit(0)
        end if
        idch=first
        minStrip=xd(idch)
        maxStrip=xd(idch)
        wherePeriodic=nextThrough(idch,:)
        idch=nextCH(idch)
        do while (idch/=first)
            coordinate=xd(idch)+wherePeriodic(1)*sqrx
            if (coordinate<minStrip) then
                minStrip=coordinate
            else if (coordinate>maxStrip) then
                maxStrip=coordinate
            end if
            wherePeriodic=wherePeriodic+nextThrough(idch,:)
            idch=nextCH(idch)
        end do
        coordinateStripOut=minStrip
        thicknessStripOut=maxStrip-minStrip
    else if (throughFromEnd2End(2)==0) then
        if(throughFromEnd2End(1)>1) then
            write(*,*) 'erreur : cluster 2 fois enroulé sans inclinaison !'
            call exit(0)
        end if
        idch=first
        minStrip=yd(idch)
        maxStrip=yd(idch)
        wherePeriodic=nextThrough(idch,:)
        idch=nextCH(idch)
        do while (idch/=first)
            coordinate=yd(idch)+wherePeriodic(2)*sqry
            if (coordinate<minStrip) then
                minStrip=coordinate
            else if (coordinate>maxStrip) then
                maxStrip=coordinate
            end if
            wherePeriodic=wherePeriodic+nextThrough(idch,:)
            idch=nextCH(idch)
        end do
        coordinateStripOut=minStrip
        thicknessStripOut=maxStrip-minStrip
    else if (throughFromEnd2End(1)>1 .and. throughFromEnd2End(2)>1) then
        !!cluster invades whole domain : to discuss
    else




        !! to calculate : 2 discs borders of the strip
        !! "orthogonal distance" between them
        !!! this bloc calculates the two discs at the border of the strip
        idch=first
        orthogonalDistanceMax=0.
        orthogonalDistanceMin=0.
        iMinOrthDistance=idch
        iMaxOrthDistance=idch
        wherePeriodic=nextThrough(idch,:)
        idch=nextCH(idch)
        normOrientation=sqrt(float(throughFromEnd2End(1)**2+throughFromEnd2End(2)**2))
        do while(idch/=first)
            orthogonalDistance=(throughFromEnd2End(1)*&
            (yd(idch)+wherePeriodic(2)*sqry-yd(first))&
            -throughFromEnd2End(2)*(xd(idch)+wherePeriodic(1)*sqrx-xd(first)))/normOrientation
            if (orthogonalDistance>orthogonalDistanceMax) then
                orthogonalDistanceMax=orthogonalDistance
                iMaxOrthDistance=idch
            else if (orthogonalDistance<orthogonalDistanceMin) then
                orthogonalDistanceMin=orthogonalDistance
                iMinOrthDistance=idch
            end if
            wherePeriodic=wherePeriodic+nextThrough(idch,:)
            idch=nextCH(idch)
        end do

        if (mainOrientation==2) then
            coordinateStripOut=mod(xd(iMaxOrthDistance)&
            -float(throughfromEnd2End(1))/throughFromEnd2End(2)*yd(iMaxOrthDistance),sqrx)
            thicknessStripOut=(orthogonalDistanceMax-orthogonalDistanceMin)&
            *normOrientation/throughFromEnd2End(2)
        else
            coordinateStripOut=mod(yd(iMaxOrthDistance)&
            -float(throughfromEnd2End(2))/throughFromEnd2End(1)*xd(iMaxOrthDistance),sqry)
            thicknessStripOut=(orthogonalDistanceMax-orthogonalDistanceMin)&
            *normOrientation/throughFromEnd2End(1)
        end if

    end if
    return
end subroutine


subroutine describeImages(first,ndmax,nInCH,limitDomain,nextCH,nextThrough,whereImages,nbImages)
    implicit none
    integer,intent(in)::ndmax
    integer,intent(in)::first,nInCH
    integer,dimension(ndmax,2),intent(in)::limitDomain,nextThrough
    integer,dimension(ndmax),intent(in)::nextCH
    integer,dimension(10,2),intent(out)::whereImages
    integer,intent(out)::nbImages
    integer::idch,icount,i
    integer,dimension(2)::wherePeriodic,wherePeriodicTemp
    logical::alreadyIn


    idch=first
    icount=1
    wherePeriodic=[0,0]
    do i=1,10
        whereImages(i,:)=[0,0]
    end do
    nbImages=1

    do while(icount<=nInCH)
        if (any(limitDomain(idch,:)/=[0,0])) then
            alreadyIn=.false.
            do i=1,nbImages
                if (all(wherePeriodic+limitDomain(idch,:)==whereImages(i,:))) alreadyIn=.true.
            end do
            if (.not. alreadyIn) then
                nbImages=nbImages+1
                whereImages(nbImages,:)=wherePeriodic+limitDomain(idch,:)
            end if
        end if

        if (any(nextThrough(idch,:)/=[0,0])) then
            if (nextThrough(idch,1)==0 .or. nextThrough(idch,2)==0) then
                wherePeriodicTemp=wherePeriodic+nextThrough(idch,:)
                alreadyIn=.false.
                do i=1,nbImages
                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                end do
                if (.not. alreadyIn) then
                    nbImages=nbImages+1
                    whereImages(nbImages,:)=wherePeriodicTemp
                end if
            else
                wherePeriodicTemp=wherePeriodic+[nextThrough(idch,1),0]
                alreadyIn=.false.
                do i=1,nbImages
                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                end do
                if (.not. alreadyIn) then
                    nbImages=nbImages+1
                    whereImages(nbImages,:)=wherePeriodicTemp
                end if

                wherePeriodicTemp=wherePeriodic+[0,nextThrough(idch,2)]
                alreadyIn=.false.
                do i=1,nbImages
                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                end do
                if (.not. alreadyIn) then
                    nbImages=nbImages+1
                    whereImages(nbImages,:)=wherePeriodicTemp
                end if

                wherePeriodicTemp=wherePeriodic+nextThrough(idch,:)
                alreadyIn=.false.
                do i=1,nbImages
                    if (all(wherePeriodicTemp==whereImages(i,:))) alreadyIn=.true.
                end do
                if (.not. alreadyIn) then
                    nbImages=nbImages+1
                    whereImages(nbImages,:)=wherePeriodicTemp
                end if
            end if
        end if
        wherePeriodic=wherePeriodic+nextThrough(idch,:)
        idch=nextCH(idch)
        icount=icount+1
    end do

!    !!!
!    icount=1
!    idch=first
!    !possibleContactThrough=.false.
!    do while(icount<=nInCH)
!        if (any(limitDomain(nd,:)==-limitDomain(idch,:))) then
!            alreadyIn=.false.
!            do i=1,nbImages
!                if (all(limitDomain(idch,:)==whereImages(i,:))) alreadyIn=.true.
!            end do
!            if (.not. alreadyIn) then
!                nbImages=nbImages+1
!                whereImages(nbImages,:)=limitDomain(idch,:)
!            end if
!        end if
!        icount=icount+1
!        idch=nextCH(idch)
!    end do
    return
end subroutine

function areaOfaStrip(orientation,thickness,sqrx,sqry,rd)
    implicit none
    integer,dimension(2),intent(in)::orientation
    real,intent(in)::thickness,sqrx,sqry,rd
    real::areaOfaStrip
    real,dimension(2)::domainSize
    integer,dimension(1)::dumb
    integer::mainDirection
    dumb=maxloc(orientation)
    mainDirection=dumb(1)
    domainSize=[sqrx,sqry]
    areaOfaStrip=domainSize(mainDirection)*(thickness+2*rd)*orientation(mainDirection)
end function

subroutine contactOrInsideStrip(sqrx,sqry,xdnd,ydnd,rd,orientationStripIn&
    ,coordinateStripIn,thicknessStripIn,contactStrip,isInsideStrip)
    implicit none
    integer,dimension(2),intent(in)::orientationStripIn
    real,intent(in)::coordinateStripIn,thicknessStripIn,rd,sqrx,sqry,xdnd,ydnd
    logical,intent(out)::contactStrip,isInsideStrip
    real::slope,proj,relativePosition
    if (orientationStripIn(2)>=orientationStripIn(1)) then
        if (orientationStripIn(1)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(1))/float(orientationStripIn(2))
        end if
        proj=orientationStripIn(2)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqrx,sqrx/float(orientationStripIn(2)))
        contactStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
        if (contactStrip) then
            relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn+2*sqrx,sqrx/float(orientationStripIn(2)))
            isInsideStrip=(relativePosition<=thicknessStripIn)
        end if
    else
        if (orientationStripIn(2)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(2))/float(orientationStripIn(1))
        end if
        proj=orientationStripIn(1)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqry,sqry/float(orientationStripIn(1)))
        contactStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
        if (contactStrip) then
            relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn+2*sqry,sqry/float(orientationStripIn(1)))
            isInsideStrip=(relativePosition<=thicknessStripIn)
        end if
    end if
    return
end subroutine

function contactOrNotStrip(sqrx,sqry,xdnd,ydnd,rd,orientationStripIn,coordinateStripIn,thicknessStripIn)
    implicit none
    integer,dimension(2),intent(in)::orientationStripIn
    real,intent(in)::coordinateStripIn,thicknessStripIn,rd,sqrx,sqry,xdnd,ydnd
    logical::contactOrNotStrip
    real::slope,proj,relativePosition
    if (orientationStripIn(2)>=orientationStripIn(1)) then
        if (orientationStripIn(1)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(1))/float(orientationStripIn(2))
        end if
        proj=orientationStripIn(2)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqrx,sqrx/float(orientationStripIn(2)))
        contactOrNotStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
!        if (contactStrip) then
!            relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn+2*sqrx,sqrx/float(orientationStripIn(2)))
!            isInsideStrip=(relativePosition<=thicknessStripIn)
!        end if
    else
        if (orientationStripIn(2)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(2))/float(orientationStripIn(1))
        end if
        proj=orientationStripIn(1)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqry,sqry/float(orientationStripIn(1)))
        contactOrNotStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
!        if (contactStrip) then
!            relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn+2*sqry,sqry/float(orientationStripIn(1)))
!            isInsideStrip=(relativePosition<=thicknessStripIn)
!        end if
    end if
    return
end function

subroutine periodicConvexification(xd,yd,firstCluster,nextDisc,ndmax,sqrx,sqry,nextThrough,nbInThisCH,listDiscsCH)
    implicit none
    integer,intent(in):: ndmax,firstCluster
    real, dimension(ndmax), intent(in)::xd,yd
    integer, dimension(ndmax), intent(in)::nextDisc
    real,intent(in)::sqrx,sqry
    integer, dimension(ndmax,2), intent(inout)::nextThrough
    integer, intent(out)::nbInThisCH
    integer, dimension(ndmax), intent(out)::listDiscsCH
    integer, dimension(ndmax,2)::whereDisc
    integer,dimension(2)::wherePeriodic
    integer::iloop,iminy,iBigLoop,icount,idc,imaxcos,i
    real::maxCos,ux,uy,cosAngleIter,cosAngle

    wherePeriodic=[0,0]
    idc=firstCluster
    whereDisc(idc,:)=wherePeriodic
    wherePeriodic=wherePeriodic+nextThrough(idc,:)
    idc=nextDisc(idc)
    do while(idc/=firstCluster)
        whereDisc(idc,:)=wherePeriodic
        wherePeriodic=wherePeriodic+nextThrough(idc,:)
        idc=nextDisc(idc)
    end do


    !!find point of minimal y
    iminy=firstCluster
    iloop=firstCluster
    do while (nextDisc(iloop)/=firstCluster)
        iloop=nextDisc(iloop)
        if (yd(iloop)+whereDisc(iloop,2)*sqry<yd(iminy)+whereDisc(iminy,2)*sqry) then
            iminy=iloop
        end if
    end do



    ux=1.
    uy=0.
    iBigLoop=iminy
    listDiscsCH(1)=iminy
    icount=1
    do
        imaxCos=firstCluster
        iloop=firstCluster
        maxCos=-2.
        do
            if (iloop/=iBigLoop) then
                cosAngleIter=cosAngle(ux,uy,xd(iloop)-xd(iBigLoop)+(whereDisc(iloop,1)-whereDisc(iBigLoop,1))*sqrx,&
                yd(iloop)-yd(iBigLoop)+(whereDisc(iloop,2)-whereDisc(iBigLoop,2))*sqry)
                if (cosAngleIter>maxCos) then
                    maxCos=cosAngleIter
                    imaxCos=iloop
                end if
            end if
            iloop=nextDisc(iloop)
            if (iloop==firstCluster) exit
        end do
        if (imaxCos==iminy) exit
        icount=icount+1
        listDiscsCH(icount)=imaxCos
        ux=xd(imaxCos)-xd(iBigLoop)+(whereDisc(imaxCos,1)-whereDisc(iBigLoop,1))*sqrx
        uy=yd(imaxCos)-yd(iBigLoop)+(whereDisc(imaxCos,2)-whereDisc(iBigLoop,2))*sqry
        iBigLoop=imaxCos
        if (icount>ndmax) then
            write(*,*) "error during construction of convex hull : "
            write(*,*) "the algorithm doesn't find initial exterior disc !!"
            call exit()
        end if
    end do
    nbInThisCH=icount

    do i=1,nbInThisCH-1
        nextThrough(listDiscsCH(i),:)=whereDisc(listDiscsCH(i+1),:)-whereDisc(listDiscsCH(i),:)
    end do
    nextThrough(listDiscsCH(nbInThisCH),:)=whereDisc(listDiscsCH(1),:)-whereDisc(listDiscsCH(nbInThisCH),:)


    return
end subroutine

subroutine contactWhereStrip(sqrx,sqry,xdnd,ydnd,rd,orientationStripIn&
    ,coordinateStripIn,thicknessStripIn,contactOrNotStrip,whichStripContact)
    implicit none
    integer,dimension(2),intent(in)::orientationStripIn
    real,intent(in)::coordinateStripIn,thicknessStripIn,rd,sqrx,sqry,xdnd,ydnd
    logical,intent(out)::contactOrNotStrip
    integer,intent(out)::whichStripContact
    real::slope,proj,relativePosition
    whichStripContact=-1
    if (orientationStripIn(2)>=orientationStripIn(1)) then
        if (orientationStripIn(1)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(1))/float(orientationStripIn(2))
        end if
        proj=orientationStripIn(2)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqrx,sqrx/float(orientationStripIn(2)))
        contactOrNotStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
        if (contactOrNotStrip) then
            whichStripContact=floor((xdnd-ydnd*slope-coordinateStripIn&
            +2*rd*proj+2*sqrx)/(sqrx/float(orientationStripIn(2))))
        end if
        !xd=whichStrip*sqrx/flo
!        if (contactStrip) then
!            relativePosition=mod(xdnd-ydnd*slope-coordinateStripIn+2*sqrx,sqrx/float(orientationStripIn(2)))
!            isInsideStrip=(relativePosition<=thicknessStripIn)
!        end if
    else
        if (orientationStripIn(2)==0) then
            slope=0.
        else
            slope=sign(1,orientationStripIn(2))/float(orientationStripIn(1))
        end if
        proj=orientationStripIn(1)/sqrt(float(orientationStripIn(2)**2+orientationStripIn(1)**2))
        relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqry,sqry/float(orientationStripIn(1)))
        contactOrNotStrip=(relativePosition<=thicknessStripIn+4*rd*proj)
        if (contactOrNotStrip) then
            whichStripContact=floor((ydnd-xdnd*slope-coordinateStripIn&
        +2*rd*proj+2*sqry)/(sqry/float(orientationStripIn(1))))
        end if
!        if (contactStrip) then
!            relativePosition=mod(ydnd-xdnd*slope-coordinateStripIn+2*sqry,sqry/float(orientationStripIn(1)))
!            isInsideStrip=(relativePosition<=thicknessStripIn)
!        end if
    end if
    return
end subroutine
!!think of unite the previous and the following subroutines to gain time
subroutine twoPointsStrip(coordinateStrip,thicknessStrip,orientationStripIn,whichStrip,sqrx,sqry,x1,y1,x2,y2)
    implicit none
    real,intent(in)::coordinateStrip,thicknessStrip,sqrx,sqry
    integer,intent(in),dimension(2)::orientationStripIn
    integer,intent(in)::whichStrip
    real,intent(out)::x1,y1,x2,y2
    if (orientationStripIn(2)>=orientationStripIn(1)) then
        y1=0.; y2=0.
        x1=whichStrip*sqrx/float(orientationStripIn(2))+coordinateStrip-2*sqrx
        x2=x1+thicknessStrip
    else
        x1=0.;x2=0.
        y1=whichStrip*sqry/float(orientationStripIn(1))+coordinateStrip-2*sqry
        y2=y1+thicknessStrip
    end if

end subroutine



!function cosAngle(ux,uy,vx,vy)
!    real,intent(in)::ux,uy,vx,vy
!    real::cosAngle
!    cosAngle=(ux*vx+uy*vy)/sqrt(ux**2+uy**2)/sqrt(vx**2+vy**2)
!    return
!end function

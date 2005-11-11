      subroutine circes (xx1m, xx2m, xroots, xacc, xver, xrev, xchat)
      implicit none
      double precision xx1m, xx2m, xroots
      integer xacc, xver, xrev, xchat
      integer SBAND, TESLA, XBAND
      parameter (SBAND  =  1, TESLA  =  2, XBAND  =  3)
      integer SBNDEE, TESLEE, XBNDEE
      parameter (SBNDEE =  4, TESLEE =  5, XBNDEE =  6)
      integer NACC
      parameter (NACC = 6)
      integer MAGIC0
      parameter (MAGIC0 = 1904 06 16)
      double precision x1m, x2m, roots
      common /circom/ x1m, x2m, roots
      double precision lumi
      common /circom/ lumi
      double precision a1(0:7)
      common /circom/ a1
      integer acc, ver, rev, chat
      common /circom/ acc, ver, rev, chat
      integer magic
      common /circom/ magic
      integer e, r
      common /circom/ e, r
      save /circom/
      character*60 msgbuf
      character*6 accnam(NACC)
      integer zver
      save zver
      integer ver34
      integer GEV350, GEV500, GEV800, TEV1, TEV16
      parameter (GEV350 = 1, GEV500 = 2, GEV800 = 3,
     $           TEV1 =   4, TEV16  = 5)
      integer A1NEGY, A1NREV
      parameter (A1NEGY = 5, A1NREV = 5)
      integer i
      real xa1lum(A1NEGY,NACC,0:A1NREV)
      real xa1(0:7,A1NEGY,NACC,0:A1NREV)
      integer A3NEGY, A3NREV
      parameter (A3NEGY = 5, A3NREV = 5)
      real xa3lum(A3NEGY,NACC,0:A3NREV)
      real xa3(0:7,A3NEGY,NACC,0:A3NREV)
      integer A5NEGY, A5NREV
      parameter (A5NEGY = 5, A5NREV = 1)
      real xa5lum(A5NEGY,NACC,0:A5NREV)
      real xa5(0:7,A5NEGY,NACC,0:A5NREV)
      data accnam(SBAND)  /'SBAND'/
      data accnam(TESLA)  /'TESLA'/
      data accnam(XBAND)  /'XBAND'/
      data accnam(SBNDEE) /'SBNDEE'/
      data accnam(TESLEE) /'TESLEE'/
      data accnam(XBNDEE) /'XBNDEE'/
      data zver / -1 /
      data xa1lum(GEV500,SBAND,1) /  5.212299E+01 /
      data (xa1(i,GEV500,SBAND,1),i=0,7) /
     $    .39192E+00,   .66026E+00,   .11828E+02,  -.62543E+00, 
     $    .52292E+00,  -.69245E+00,   .14983E+02,   .65421E+00 /
      data xa1lum(GEV500,TESLA,1) /  6.066178E+01 /
      data (xa1(i,GEV500,TESLA,1),i=0,7) /
     $    .30196E+00,   .12249E+01,   .21423E+02,  -.57848E+00, 
     $    .68766E+00,  -.69788E+00,   .23121E+02,   .78399E+00 /
      data xa1lum(GEV500,XBAND,1) /  5.884699E+01 /
      data (xa1(i,GEV500,XBAND,1),i=0,7) /
     $    .48594E+00,   .52435E+00,   .83585E+01,  -.61347E+00, 
     $    .30703E+00,  -.68804E+00,   .84109E+01,   .44312E+00 /
      data xa1lum(TEV1,SBAND,1)   /  1.534650E+02 /
      data (xa1(i,TEV1,SBAND,1),i=0,7) /
     $    .24399E+00,   .87464E+00,   .66751E+01,  -.56808E+00, 
     $    .59295E+00,  -.68921E+00,   .94232E+01,   .83351E+00 /
      data xa1lum(TEV1,TESLA,1)   /  1.253381E+03 /
      data (xa1(i,TEV1,TESLA,1),i=0,7) /
     $    .39843E+00,   .70097E+00,   .11602E+02,  -.61061E+00, 
     $    .40737E+00,  -.69319E+00,   .14800E+02,   .51382E+00 /
      data xa1lum(TEV1,XBAND,1)   /  1.901783E+02 /
      data (xa1(i,TEV1,XBAND,1),i=0,7) /
     $    .32211E+00,   .61798E+00,   .28298E+01,  -.54644E+00, 
     $    .45674E+00,  -.67301E+00,   .41703E+01,   .74536E+00 /
      data (xa1lum(GEV350,i,1),i=1,NACC) / NACC*-1d0 /
      data (xa1lum(GEV800,i,1),i=1,NACC) / NACC*-1d0 /
      data (xa1lum(GEV500,i,1),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,1),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,1),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND,2) /   .31057E+02 /
      data (xa1(i,GEV500,SBAND,2),i=0,7) /
     $    .38504E+00,   .79723E+00,   .14191E+02,  -.60456E+00, 
     $    .53411E+00,  -.68873E+00,   .15105E+02,   .65151E+00 /
      data xa1lum(TEV1,SBAND,2) /   .24297E+03 /
      data (xa1(i,TEV1,SBAND,2),i=0,7) /
     $    .24374E+00,   .89466E+00,   .70242E+01,  -.56754E+00, 
     $    .60910E+00,  -.68682E+00,   .96083E+01,   .83985E+00 /
      data xa1lum(GEV350,TESLA,2) /   .73369E+02 /
      data (xa1(i,GEV350,TESLA,2),i=0,7) /
     $    .36083E+00,   .12819E+01,   .37880E+02,  -.59492E+00, 
     $    .69109E+00,  -.69379E+00,   .40061E+02,   .65036E+00 /
      data xa1lum(GEV500,TESLA,2) /   .10493E+03 /
      data (xa1(i,GEV500,TESLA,2),i=0,7) /
     $    .29569E+00,   .11854E+01,   .21282E+02,  -.58553E+00, 
     $    .71341E+00,  -.69279E+00,   .24061E+02,   .77709E+00 /
      data xa1lum(GEV800,TESLA,2) /   .28010E+03 /
      data (xa1(i,GEV800,TESLA,2),i=0,7) /
     $    .22745E+00,   .11265E+01,   .10483E+02,  -.55711E+00, 
     $    .69579E+00,  -.69068E+00,   .13093E+02,   .89605E+00 /
      data xa1lum(TEV1,TESLA,2) /   .10992E+03 /
      data (xa1(i,TEV1,TESLA,2),i=0,7) /
     $    .40969E+00,   .66105E+00,   .11972E+02,  -.62041E+00, 
     $    .40463E+00,  -.69354E+00,   .14669E+02,   .51281E+00 /
      data xa1lum(GEV500,XBAND,2) /   .35689E+02 /
      data (xa1(i,GEV500,XBAND,2),i=0,7) /
     $    .48960E+00,   .46815E+00,   .75249E+01,  -.62769E+00, 
     $    .30341E+00,  -.68754E+00,   .85545E+01,   .43453E+00 /
      data xa1lum(TEV1,XBAND,2) /   .11724E+03 /
      data (xa1(i,TEV1,XBAND,2),i=0,7) /
     $    .31939E+00,   .62415E+00,   .30763E+01,  -.55314E+00, 
     $    .45634E+00,  -.67089E+00,   .41529E+01,   .73807E+00 /
      data xa1lum(GEV350,SBAND,2) / -1d0 /
      data xa1lum(GEV350,XBAND,2) / -1d0 /
      data xa1lum(GEV800,SBAND,2) / -1d0 /
      data xa1lum(GEV800,XBAND,2) / -1d0 /
      data (xa1lum(GEV350,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,2),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,2),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND, 3) /   .31469E+02 /
      data (xa1(i,GEV500,SBAND, 3),i=0,7) /
     $    .38299E+00,   .72035E+00,   .12618E+02,  -.61611E+00, 
     $    .51971E+00,  -.68960E+00,   .15066E+02,   .63784E+00 /
      data xa1lum(TEV1,  SBAND, 3) /   .24566E+03 /
      data (xa1(i,TEV1,  SBAND, 3),i=0,7) /
     $    .24013E+00,   .95763E+00,   .69085E+01,  -.55151E+00, 
     $    .59497E+00,  -.68622E+00,   .94494E+01,   .82158E+00 /
      data xa1lum(GEV350,TESLA, 3) /   .74700E+02 /
      data (xa1(i,GEV350,TESLA, 3),i=0,7) /
     $    .34689E+00,   .12484E+01,   .33720E+02,  -.59523E+00, 
     $    .66266E+00,  -.69524E+00,   .38488E+02,   .63775E+00 /
      data xa1lum(GEV500,TESLA, 3) /   .10608E+03 /
      data (xa1(i,GEV500,TESLA, 3),i=0,7) /
     $    .28282E+00,   .11700E+01,   .19258E+02,  -.58390E+00, 
     $    .68777E+00,  -.69402E+00,   .23638E+02,   .75929E+00 /
      data xa1lum(GEV800,TESLA, 3) /   .28911E+03 /
      data (xa1(i,GEV800,TESLA, 3),i=0,7) /
     $    .21018E+00,   .12039E+01,   .96763E+01,  -.54024E+00, 
     $    .67220E+00,  -.69083E+00,   .12733E+02,   .87355E+00 /
      data xa1lum(TEV1,  TESLA, 3) /   .10936E+03 /
      data (xa1(i,TEV1,  TESLA, 3),i=0,7) /
     $    .41040E+00,   .68099E+00,   .11610E+02,  -.61237E+00, 
     $    .40155E+00,  -.69073E+00,   .14698E+02,   .49989E+00 /
      data xa1lum(GEV500,XBAND, 3) /   .36145E+02 /
      data (xa1(i,GEV500,XBAND, 3),i=0,7) /
     $    .51285E+00,   .45812E+00,   .75135E+01,  -.62247E+00, 
     $    .30444E+00,  -.68530E+00,   .85519E+01,   .43062E+00 /
      data xa1lum(TEV1,  XBAND, 3) /   .11799E+03 /
      data (xa1(i,TEV1,  XBAND, 3),i=0,7) /
     $    .31241E+00,   .61241E+00,   .29938E+01,  -.55848E+00, 
     $    .44801E+00,  -.67116E+00,   .41119E+01,   .72753E+00 /
      data xa1lum(GEV350,SBAND,3) / -1d0 /
      data xa1lum(GEV350,XBAND,3) / -1d0 /
      data xa1lum(GEV800,SBAND,3) / -1d0 /
      data xa1lum(GEV800,XBAND,3) / -1d0 /
      data (xa1lum(GEV350,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,3),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,3),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV500,SBAND, 4) /   .31528E+02 /
      data (xa1(i,GEV500,SBAND, 4),i=0,7) /
     $    .38169E+00,   .73949E+00,   .12543E+02,  -.61112E+00, 
     $    .51256E+00,  -.69009E+00,   .14892E+02,   .63314E+00 /
      data xa1lum(TEV1,  SBAND, 4) /   .24613E+03 /
      data (xa1(i,TEV1,  SBAND, 4),i=0,7) /
     $    .24256E+00,   .94117E+00,   .66775E+01,  -.55160E+00, 
     $    .57484E+00,  -.68891E+00,   .92271E+01,   .81162E+00 /
      data xa1lum(GEV350,TESLA, 4) /   .74549E+02 /
      data (xa1(i,GEV350,TESLA, 4),i=0,7) /
     $    .34120E+00,   .12230E+01,   .32932E+02,  -.59850E+00, 
     $    .65947E+00,  -.69574E+00,   .38116E+02,   .63879E+00 /
      data xa1lum(GEV500,TESLA, 4) /   .10668E+03 /
      data (xa1(i,GEV500,TESLA, 4),i=0,7) /
     $    .28082E+00,   .11074E+01,   .18399E+02,  -.59118E+00, 
     $    .68880E+00,  -.69375E+00,   .23463E+02,   .76073E+00 /
      data xa1lum(GEV800,TESLA, 4) /   .29006E+03 /
      data (xa1(i,GEV800,TESLA, 4),i=0,7) /
     $    .21272E+00,   .11443E+01,   .92564E+01,  -.54657E+00, 
     $    .66799E+00,  -.69137E+00,   .12498E+02,   .87571E+00 /
      data xa1lum(TEV1,  TESLA, 4) /   .11009E+03 /
      data (xa1(i,TEV1,  TESLA, 4),i=0,7) /
     $    .41058E+00,   .64745E+00,   .11271E+02,  -.61996E+00, 
     $    .39801E+00,  -.69150E+00,   .14560E+02,   .49924E+00 /
      data xa1lum(GEV500,XBAND, 4) /   .36179E+02 /
      data (xa1(i,GEV500,XBAND, 4),i=0,7) /
     $    .51155E+00,   .43313E+00,   .70446E+01,  -.63003E+00, 
     $    .29449E+00,  -.68747E+00,   .83489E+01,   .42458E+00 /
      data xa1lum(TEV1,  XBAND, 4) /   .11748E+03 /
      data (xa1(i,TEV1,  XBAND, 4),i=0,7) /
     $    .32917E+00,   .54322E+00,   .28493E+01,  -.57959E+00, 
     $    .39266E+00,  -.68217E+00,   .38475E+01,   .68478E+00 /
      data xa1lum(GEV350,SBAND,4) / -1d0 /
      data xa1lum(GEV350,XBAND,4) / -1d0 /
      data xa1lum(GEV800,SBAND,4) / -1d0 /
      data xa1lum(GEV800,XBAND,4) / -1d0 /
      data (xa1lum(GEV350,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV500,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(GEV800,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV1,i,4),i=SBNDEE,NACC) / 3*-1d0 /
      data (xa1lum(TEV16,i,4),i=1,NACC) / 6*-1d0 /
      data xa1lum(GEV350,SBAND, 5) /  0.21897E+02 /
      data (xa1(i,GEV350,SBAND, 5),i=0,7) /
     $   0.57183E+00,  0.53877E+00,  0.19422E+02, -0.63064E+00, 
     $   0.49112E+00, -0.69109E+00,  0.24331E+02,  0.52718E+00 /
      data xa1lum(GEV500,SBAND, 5) /  0.31383E+02 /
      data (xa1(i,GEV500,SBAND, 5),i=0,7) /
     $   0.51882E+00,  0.49915E+00,  0.11153E+02, -0.63017E+00, 
     $   0.50217E+00, -0.69113E+00,  0.14935E+02,  0.62373E+00 /
      data xa1lum(GEV800,SBAND, 5) /  0.95091E+02 /
      data (xa1(i,GEV800,SBAND, 5),i=0,7) /
     $   0.47137E+00,  0.46150E+00,  0.56562E+01, -0.61758E+00, 
     $   0.46863E+00, -0.68897E+00,  0.85876E+01,  0.67577E+00 /
      data xa1lum(TEV1,SBAND, 5) /  0.11900E+03 /
      data (xa1(i,TEV1,SBAND, 5),i=0,7) /
     $   0.43956E+00,  0.45471E+00,  0.42170E+01, -0.61180E+00, 
     $   0.48711E+00, -0.68696E+00,  0.67145E+01,  0.74551E+00 /
      data xa1lum(TEV16,SBAND, 5) /  0.11900E+03 /
      data (xa1(i,TEV16,SBAND, 5),i=0,7) /
     $   0.43956E+00,  0.45471E+00,  0.42170E+01, -0.61180E+00, 
     $   0.48711E+00, -0.68696E+00,  0.67145E+01,  0.74551E+00 /
      data xa1lum(GEV350,TESLA, 5) /  0.97452E+02 /
      data (xa1(i,GEV350,TESLA, 5),i=0,7) /
     $   0.39071E+00,  0.84996E+00,  0.17614E+02, -0.60609E+00, 
     $   0.73920E+00, -0.69490E+00,  0.28940E+02,  0.77286E+00 /
      data xa1lum(GEV500,TESLA, 5) /  0.10625E+03 /
      data (xa1(i,GEV500,TESLA, 5),i=0,7) /
     $   0.42770E+00,  0.71457E+00,  0.15284E+02, -0.61664E+00, 
     $   0.68166E+00, -0.69208E+00,  0.24165E+02,  0.73806E+00 /
      data xa1lum(GEV800,TESLA, 5) /  0.17086E+03 /
      data (xa1(i,GEV800,TESLA, 5),i=0,7) /
     $   0.36025E+00,  0.69118E+00,  0.76221E+01, -0.59440E+00, 
     $   0.71269E+00, -0.69077E+00,  0.13117E+02,  0.91780E+00 /
      data xa1lum(TEV1,TESLA, 5) /  0.21433E+03 /
      data (xa1(i,TEV1,TESLA, 5),i=0,7) /
     $   0.33145E+00,  0.67075E+00,  0.55438E+01, -0.58468E+00, 
     $   0.72503E+00, -0.69084E+00,  0.99992E+01,  0.10112E+01 /
      data xa1lum(TEV16,TESLA, 5) /  0.34086E+03 /
      data (xa1(i,TEV16,TESLA, 5),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa1lum(GEV350,XBAND, 5) /  0.31901E+02 /
      data (xa1(i,GEV350,XBAND, 5),i=0,7) /
     $   0.65349E+00,  0.31752E+00,  0.94342E+01, -0.64291E+00, 
     $   0.30364E+00, -0.68989E+00,  0.11446E+02,  0.40486E+00 /
      data xa1lum(GEV500,XBAND, 5) /  0.36386E+02 /
      data (xa1(i,GEV500,XBAND, 5),i=0,7) /
     $   0.65132E+00,  0.28728E+00,  0.69853E+01, -0.64440E+00, 
     $   0.28736E+00, -0.68758E+00,  0.83227E+01,  0.41492E+00 /
      data xa1lum(GEV800,XBAND, 5) /  0.10854E+03 /
      data (xa1(i,GEV800,XBAND, 5),i=0,7) /
     $   0.49478E+00,  0.36221E+00,  0.30116E+01, -0.61548E+00, 
     $   0.39890E+00, -0.68418E+00,  0.45183E+01,  0.67243E+00 /
      data xa1lum(TEV1,XBAND, 5) /  0.11899E+03 /
      data (xa1(i,TEV1,XBAND, 5),i=0,7) /
     $   0.49992E+00,  0.34299E+00,  0.26184E+01, -0.61584E+00, 
     $   0.38450E+00, -0.68342E+00,  0.38589E+01,  0.67408E+00 /
      data xa1lum(TEV16,XBAND, 5) /  0.13675E+03 /
      data (xa1(i,TEV16,XBAND, 5),i=0,7) /
     $   0.50580E+00,  0.30760E+00,  0.18339E+01, -0.61421E+00, 
     $   0.35233E+00, -0.68315E+00,  0.26708E+01,  0.67918E+00 /
      data xa1lum(GEV500,SBNDEE, 0) /   .92914E+01 /
      data (xa1(i,GEV500,SBNDEE, 0),i=0,7) /
     $    .34866E+00,   .78710E+00,   .10304E+02,  -.59464E+00, 
     $    .40234E+00,  -.69741E+00,   .20645E+02,   .47274E+00 /
      data xa1lum(TEV1,  SBNDEE, 0) /   .45586E+02 /
      data (xa1(i,TEV1,  SBNDEE, 0),i=0,7) /
     $    .21084E+00,   .99168E+00,   .54407E+01,  -.52851E+00, 
     $    .47493E+00,  -.69595E+00,   .12480E+02,   .64027E+00 /
      data xa1lum(GEV350,TESLEE, 0) /   .15175E+02 /
      data (xa1(i,GEV350,TESLEE, 0),i=0,7) /
     $    .33093E+00,   .11137E+01,   .25275E+02,  -.59942E+00, 
     $    .49623E+00,  -.70403E+00,   .60188E+02,   .44637E+00 /
      data xa1lum(GEV500,TESLEE, 0) /   .21622E+02 /
      data (xa1(i,GEV500,TESLEE, 0),i=0,7) /
     $    .27175E+00,   .10697E+01,   .14858E+02,  -.58418E+00, 
     $    .50824E+00,  -.70387E+00,   .36129E+02,   .53002E+00 /
      data xa1lum(GEV800,TESLEE, 0) /   .43979E+02 /
      data (xa1(i,GEV800,TESLEE, 0),i=0,7) /
     $    .22994E+00,   .10129E+01,   .81905E+01,  -.55751E+00, 
     $    .46551E+00,  -.70461E+00,   .19394E+02,   .58387E+00 /
      data xa1lum(TEV1,  TESLEE, 0) /   .25465E+02 /
      data (xa1(i,TEV1,  TESLEE, 0),i=0,7) /
     $    .37294E+00,   .67522E+00,   .87504E+01,  -.60576E+00, 
     $    .35095E+00,  -.69821E+00,   .18526E+02,   .42784E+00 /
      data xa1lum(GEV500,XBNDEE, 0) /   .13970E+02 /
      data (xa1(i,GEV500,XBNDEE, 0),i=0,7) /
     $    .47296E+00,   .46800E+00,   .58897E+01,  -.61689E+00, 
     $    .27181E+00,  -.68923E+00,   .10087E+02,   .37462E+00 /
      data xa1lum(TEV1,  XBNDEE, 0) /   .41056E+02 /
      data (xa1(i,TEV1,  XBNDEE, 0),i=0,7) /
     $    .27965E+00,   .74816E+00,   .27415E+01,  -.50491E+00, 
     $    .38320E+00,  -.67945E+00,   .47506E+01,   .62218E+00 /
      data xa1lum(GEV350,SBNDEE,0) / -1d0 /
      data xa1lum(GEV350,XBNDEE,0) / -1d0 /
      data xa1lum(GEV800,SBNDEE,0) / -1d0 /
      data xa1lum(GEV800,XBNDEE,0) / -1d0 /
      data xa1lum(GEV500,SBAND, 0) /   .31528E+02 /
      data (xa1(i,GEV500,SBAND, 0),i=0,7) /
     $    .38169E+00,   .73949E+00,   .12543E+02,  -.61112E+00, 
     $    .51256E+00,  -.69009E+00,   .14892E+02,   .63314E+00 /
      data xa1lum(TEV1,  SBAND, 0) /   .24613E+03 /
      data (xa1(i,TEV1,  SBAND, 0),i=0,7) /
     $    .24256E+00,   .94117E+00,   .66775E+01,  -.55160E+00, 
     $    .57484E+00,  -.68891E+00,   .92271E+01,   .81162E+00 /
      data xa1lum(GEV350,TESLA, 0) /   .74549E+02 /
      data (xa1(i,GEV350,TESLA, 0),i=0,7) /
     $    .34120E+00,   .12230E+01,   .32932E+02,  -.59850E+00, 
     $    .65947E+00,  -.69574E+00,   .38116E+02,   .63879E+00 /
      data xa1lum(GEV500,TESLA, 0) /   .10668E+03 /
      data (xa1(i,GEV500,TESLA, 0),i=0,7) /
     $    .28082E+00,   .11074E+01,   .18399E+02,  -.59118E+00, 
     $    .68880E+00,  -.69375E+00,   .23463E+02,   .76073E+00 /
      data xa1lum(GEV800,TESLA, 0) /   .29006E+03 /
      data (xa1(i,GEV800,TESLA, 0),i=0,7) /
     $    .21272E+00,   .11443E+01,   .92564E+01,  -.54657E+00, 
     $    .66799E+00,  -.69137E+00,   .12498E+02,   .87571E+00 /
      data xa1lum(TEV1,  TESLA, 0) /   .11009E+03 /
      data (xa1(i,TEV1,  TESLA, 0),i=0,7) /
     $    .41058E+00,   .64745E+00,   .11271E+02,  -.61996E+00, 
     $    .39801E+00,  -.69150E+00,   .14560E+02,   .49924E+00 /
      data xa1lum(GEV500,XBAND, 0) /   .36179E+02 /
      data (xa1(i,GEV500,XBAND, 0),i=0,7) /
     $    .51155E+00,   .43313E+00,   .70446E+01,  -.63003E+00, 
     $    .29449E+00,  -.68747E+00,   .83489E+01,   .42458E+00 /
      data xa1lum(TEV1,  XBAND, 0) /   .11748E+03 /
      data (xa1(i,TEV1,  XBAND, 0),i=0,7) /
     $    .32917E+00,   .54322E+00,   .28493E+01,  -.57959E+00, 
     $    .39266E+00,  -.68217E+00,   .38475E+01,   .68478E+00 /
      data xa1lum(GEV350,SBAND,0) / -1d0 /
      data xa1lum(GEV350,XBAND,0) / -1d0 /
      data xa1lum(GEV800,SBAND,0) / -1d0 /
      data xa1lum(GEV800,XBAND,0) / -1d0 /
      data xa3lum(GEV800,TESLA, 3) /   .17196E+03 /
      data (xa3(i,GEV800,TESLA, 3),i=0,7) /
     $    .21633E+00,   .11333E+01,   .95928E+01,  -.55095E+00, 
     $    .73044E+00,  -.69101E+00,   .12868E+02,   .94737E+00 /
      data xa3lum(GEV800,TESLA, 4) /   .16408E+03 /
      data (xa3(i,GEV800,TESLA, 4),i=0,7) /
     $    .41828E+00,   .72418E+00,   .14137E+02,  -.61189E+00, 
     $    .36697E+00,  -.69205E+00,   .17713E+02,   .43583E+00 /
      data xa3lum(GEV350,TESLA, 5) /  0.66447E+02 /
      data (xa3(i,GEV350,TESLA, 5),i=0,7) /
     $   0.69418E+00,  0.50553E+00,  0.48430E+02, -0.63911E+00, 
     $   0.34074E+00, -0.69533E+00,  0.55502E+02,  0.29397E+00 /
      data xa3lum(GEV500,TESLA, 5) /  0.95241E+02 /
      data (xa3(i,GEV500,TESLA, 5),i=0,7) /
     $   0.64882E+00,  0.45462E+00,  0.27103E+02, -0.64535E+00, 
     $   0.35101E+00, -0.69467E+00,  0.33658E+02,  0.35024E+00 /
      data xa3lum(GEV800,TESLA, 5) /  0.16974E+03 /
      data (xa3(i,GEV800,TESLA, 5),i=0,7) /
     $   0.58706E+00,  0.43771E+00,  0.13422E+02, -0.63804E+00, 
     $   0.35541E+00, -0.69467E+00,  0.17528E+02,  0.43051E+00 /
      data xa3lum(TEV1,TESLA, 5) /  0.21222E+03 /
      data (xa3(i,TEV1,TESLA, 5),i=0,7) /
     $   0.55525E+00,  0.42577E+00,  0.96341E+01, -0.63587E+00, 
     $   0.36448E+00, -0.69365E+00,  0.13161E+02,  0.47715E+00 /
      data xa3lum(TEV16,TESLA, 5) /  0.34086E+03 /
      data (xa3(i,TEV16,TESLA, 5),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa3lum(GEV350,TESLA, 0) /  0.66447E+02 /
      data (xa3(i,GEV350,TESLA, 0),i=0,7) /
     $   0.69418E+00,  0.50553E+00,  0.48430E+02, -0.63911E+00, 
     $   0.34074E+00, -0.69533E+00,  0.55502E+02,  0.29397E+00 /
      data xa3lum(GEV500,TESLA, 0) /  0.95241E+02 /
      data (xa3(i,GEV500,TESLA, 0),i=0,7) /
     $   0.64882E+00,  0.45462E+00,  0.27103E+02, -0.64535E+00, 
     $   0.35101E+00, -0.69467E+00,  0.33658E+02,  0.35024E+00 /
      data xa3lum(GEV800,TESLA, 0) /  0.16974E+03 /
      data (xa3(i,GEV800,TESLA, 0),i=0,7) /
     $   0.58706E+00,  0.43771E+00,  0.13422E+02, -0.63804E+00, 
     $   0.35541E+00, -0.69467E+00,  0.17528E+02,  0.43051E+00 /
      data xa3lum(TEV1,TESLA, 0) /  0.21222E+03 /
      data (xa3(i,TEV1,TESLA, 0),i=0,7) /
     $   0.55525E+00,  0.42577E+00,  0.96341E+01, -0.63587E+00, 
     $   0.36448E+00, -0.69365E+00,  0.13161E+02,  0.47715E+00 /
      data xa3lum(TEV16,TESLA, 0) /  0.34086E+03 /
      data (xa3(i,TEV16,TESLA, 0),i=0,7) /
     $   0.49058E+00,  0.42609E+00,  0.50550E+01, -0.61867E+00, 
     $   0.39225E+00, -0.68916E+00,  0.75514E+01,  0.58754E+00 /
      data xa5lum(GEV350,TESLA, 1) /  -1.0 /
      data xa5lum(GEV500,TESLA, 1) /  0.33980E+03 /
      data (xa5(i,GEV500,TESLA, 1),i=0,7) /
     $   0.49808E+00,  0.54613E+00,  0.12287E+02, -0.62756E+00, 
     $   0.42817E+00, -0.69120E+00,  0.17067E+02,  0.51143E+00 /
      data xa5lum(GEV800,TESLA, 1) /  0.35936E+03 /
      data (xa5(i,GEV800,TESLA, 1),i=0,7) /
     $   0.58751E+00,  0.43128E+00,  0.13324E+02, -0.64006E+00, 
     $   0.30682E+00, -0.69235E+00,  0.16815E+02,  0.37078E+00 /
      data xa5lum(TEV1,  TESLA, 1) /  -1.0 /
      data xa5lum(TEV16, TESLA, 1) /  -1.0 /
      data xa5lum(GEV350,TESLA, 0) /  -1.0 /
      data xa5lum(GEV500,TESLA, 0) /  0.33980E+03 /
      data (xa5(i,GEV500,TESLA, 0),i=0,7) /
     $   0.49808E+00,  0.54613E+00,  0.12287E+02, -0.62756E+00, 
     $   0.42817E+00, -0.69120E+00,  0.17067E+02,  0.51143E+00 /
      data xa5lum(GEV800,TESLA, 0) /  0.35936E+03 /
      data (xa5(i,GEV800,TESLA, 0),i=0,7) /
     $   0.58751E+00,  0.43128E+00,  0.13324E+02, -0.64006E+00, 
     $   0.30682E+00, -0.69235E+00,  0.16815E+02,  0.37078E+00 /
      data xa5lum(TEV1,  TESLA, 0) /  -1.0 /
      data xa5lum(TEV16, TESLA, 0) /  -1.0 /

      if (magic .ne. 1904 06 16) then
         magic = 1904 06 16
      x1m = 0d0
      x2m = 0d0
      roots = 500D0
      acc = TESLA
      ver = 0
      rev = 0
      chat = 1
      if (xchat .ne. 0) then
         call circem ('MESSAGE', 'starting up ...')
         call circem ('MESSAGE',
     $        '$Id: circes.F,v 1.1 2004/02/11 18:47:02 stavrop Exp $')
      endif
      endif
      if ((xchat .ge. 0) .and. (xchat .ne. chat)) then
         chat = xchat
         if (chat .ge. 1) then
            write (msgbuf, 1000) 'chat', chat
 1000       format ('updating `', A, ''' to ', I2)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1100) 'chat', chat
 1100       format ('keeping `', A, ''' at ', I2)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xx1m .ge. 0d0) .and. (xx1m .ne. x1m)) then
         x1m = xx1m
         if (chat .ge. 1) then
            write (msgbuf, 1001) 'x1min', x1m
 1001       format ('updating `', A, ''' to ', E12.4)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1101) 'x1min', x1m
 1101       format ('keeping `', A, ''' at ', E12.4)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xx2m .ge. 0d0) .and. (xx2m .ne. x2m)) then
         x2m = xx2m
         if (chat .ge. 1) then
            write (msgbuf, 1001) 'x2min', x2m
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1101) 'x2min', x2m
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xroots .ge. 0d0) .and.(xroots .ne. roots)) then
         roots = xroots
         if (chat .ge. 1) then
            write (msgbuf, 1002) 'roots', roots
 1002       format ('updating `', A, ''' to ', F6.1)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1102) 'roots', roots
 1102       format ('keeping `', A, ''' at ', F6.1)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xacc .ge. 0) .and.(xacc .ne. acc)) then
         if ((xacc .ge. 1) .and. (xacc .le. NACC)) then
            acc = xacc
            if (chat .ge. 1) then
               write (msgbuf, 1003) 'acc', accnam(acc)
 1003          format ('updating `', A, ''' to ', A)
               call circem ('MESSAGE', msgbuf)
            endif
         else
            write (msgbuf, 1203) xacc
 1203       format ('invalid `acc'': ', I8)
            call circem ('ERROR', msgbuf)
            write (msgbuf, 1103) 'acc', accnam(acc)
 1103       format ('keeping `', A, ''' at ', A)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1003) 'acc', accnam(acc)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((acc .eq. SBNDEE) .or. (acc .eq. TESLEE)
     $    .or. (acc .eq. XBNDEE)) then
      call circem ('WARNING', '***********************************')
      call circem ('WARNING', '* The accelerator parameters have *')
      call circem ('WARNING', '* not been endorsed for use in    *')
      call circem ('WARNING', '* an e-e- collider yet!!!         *')
      call circem ('WARNING', '***********************************')
      endif
      if ((xver .ge. 0) .and. (xver .ne. zver)) then
         ver = xver
         zver = xver
         if (chat .ge. 1) then
            write (msgbuf, 1000) 'ver', ver
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1100) 'ver', ver
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      if ((xrev .ge. 0) .and.(xrev .ne. rev)) then
         rev = xrev
         if (chat .ge. 1) then
            write (msgbuf, 1004) 'rev', rev
 1004       format ('updating `', A, ''' to ', I8)
            call circem ('MESSAGE', msgbuf)
         endif
      else
         if (chat .ge. 2) then
            write (msgbuf, 1104) 'rev', rev
 1104       format ('keeping `', A, ''' at ', I8)
            call circem ('MESSAGE', msgbuf)
         endif
      endif
      ver34 = 0
      if ((ver .eq. 1) .or. (ver .eq. 0)) then
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1997 04 17) then
         r = 5
      elseif (rev .ge. 1996 09 02) then
         r = 4
      elseif (rev .ge. 1996 07 29) then
         r = 3
      elseif (rev .ge. 1996 07 11) then
         r = 2
      elseif (rev .ge. 1996 04 01) then
         r = 1
      elseif (rev .lt. 1996 04 01) then
         call circem ('ERROR',
     $        'no revision of version 1 before 96/04/01 available')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
 2000    format ('mapping date ', I8, ' to revision index ', I2)
         call circem ('MESSAGE', msgbuf)
      endif
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa1lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2003) roots, e
         call circem ('MESSAGE', msgbuf)
      endif
      lumi = xa1lum (e,acc,r)
      do 10 i = 0, 7
         a1(i) = xa1(i,e,acc,r)
 10   continue
      elseif ((ver .eq. 3) .or. (ver .eq. 4)) then
         ver34 = ver
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1997 04 17) then
         r = 5
         if (ver34 .eq. 3) then
            call circem ('WARNING', 'version 3 retired after 97/04/17')
            call circem ('MESSAGE', 'falling back to version 4')
         endif
      elseif (rev .ge. 1996 10 22) then
         r = ver34
         if ((roots .ne. 800d0) .or. (acc .ne. TESLA)) then
            call circem ('ERROR', 'versions 3 and 4 before 97/04/17')
            call circem ('ERROR', 'apply to TESLA at 800 GeV only')
            call circem ('MESSAGE', 'falling back to TESLA at 800GeV')
            acc = TESLA
            e = GEV800
         endif
      elseif (rev .lt. 1996 10 22) then
         call circem ('ERROR',
     $     'no revision of versions 3 and 4 available before 96/10/22')
         call circem ('MESSAGE', 'falling back to default')
         r = 5
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa3lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2003) roots, e
         call circem ('MESSAGE', msgbuf)
      endif
      lumi = xa3lum (e,acc,r)
      do 20 i = 0, 7
         a1(i) = xa3(i,e,acc,r)
 20   continue
      elseif (ver .eq. 5) then
         ver = 1
      if (rev .eq. 0) then
         r = 0
      call circem ('WARNING', '*************************************')
      call circem ('WARNING', '* This release is not official yet, *')
      call circem ('WARNING', '* do not use it in publications!    *')
      call circem ('WARNING', '*************************************')
      elseif (rev .ge. 1998 05 05) then
         r = 1
      elseif (rev .lt. 1998 05 05) then
         call circem ('ERROR',
     $     'no revision of version 5 available before 98/05/05')
         call circem ('MESSAGE', 'falling back to default')
         r = 1
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2000) rev, r
         call circem ('MESSAGE', msgbuf)
      endif
      if (acc .ne. TESLA) then
         call circem ('ERROR', 'versions 5 applies to TESLA only')
         acc = TESLA
      end if
      if (roots .eq. 350d0) then
         e = GEV350
      elseif ((roots .ge. 340d0) .and. (roots .le. 370d0)) then
         write (msgbuf, 2001) roots, 350d0
         call circem ('MESSAGE', msgbuf)
         e = GEV350
      elseif (roots .eq. 500d0) then
         e = GEV500
      elseif ((roots .ge. 480d0) .and. (roots .le. 520d0)) then
         write (msgbuf, 2001) roots, 500d0
         call circem ('MESSAGE', msgbuf)
         e = GEV500
      elseif (roots .eq. 800d0) then
         e = GEV800
      elseif ((roots .ge. 750d0) .and. (roots .le. 850d0)) then
         write (msgbuf, 2001) roots, 800d0
         call circem ('MESSAGE', msgbuf)
         e = GEV800
      elseif (roots .eq. 1000d0) then
         e = TEV1
      elseif ((roots .ge. 900d0) .and. (roots .le. 1100d0)) then
         write (msgbuf, 2001) roots, 1000d0
         call circem ('MESSAGE', msgbuf)
         e = TEV1
      elseif (roots .eq. 1600d0) then
         e = TEV16
      elseif ((roots .ge. 1500d0) .and. (roots .le. 1700d0)) then
         write (msgbuf, 2001) roots, 1600d0
         call circem ('MESSAGE', msgbuf)
         e = TEV16
      else
         call circem ('ERROR',
     $        'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (xa5lum(e,acc,r) .lt. 0d0) then
         write (msgbuf, 2002) roots, accnam(acc), r
         call circem ('ERROR', msgbuf)
         call circem ('MESSAGE', 'falling back to 500GeV')
         e = GEV500
      endif
      if (chat .ge. 2) then
         write (msgbuf, 2003) roots, e
         call circem ('MESSAGE', msgbuf)
      endif
      lumi = xa5lum (e,acc,r)
      do 30 i = 0, 7
         a1(i) = xa5(i,e,acc,r)
 30   continue
      elseif (ver .eq. 2) then
      call circem ('PANIC', '*********************************')
      call circem ('PANIC', '* version 2 has been retired,   *')
      call circem ('PANIC', '* please use version 1 instead! *')
      call circem ('PANIC', '*********************************')
      return
      elseif (ver .gt. 5) then
         call circem ('PANIC', 'versions >5 not available yet')
         return
      else
         call circem ('PANIC', 'version must be positive')
         return
      endif
 2001 format ('treating energy ', F6.1, 'GeV as ',  F6.1, 'GeV')
 2002 format ('energy ', F6.1, ' not available for ', A6,
     $        ' in revison ', I2)
 2003 format ('mapping energy ', F6.1, ' to energy index ', I2)
      end

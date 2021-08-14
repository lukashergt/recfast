      SUBROUTINE DVERK(N, FCN, X, Y, XEND, TOL, IND, C, NW, W)
      use precision, only : dl => dp
      implicit none
      INTEGER N, IND, NW, K
      real(dl) X, Y(N), XEND, TOL, C(*), W(NW,9), TEMP
      ! Added external statement for fcn to avoid a warning message.
      external fcn
      !
      !***********************************************************************
      !                                                                      *
      !     PURPOSE - THIS IS A RUNGE-KUTTA  SUBROUTINE  BASED  ON  VERNER'S *
      ! FIFTH AND SIXTH ORDER PAIR OF FORMULAS FOR FINDING APPROXIMATIONS TO *
      ! THE SOLUTION OF  A  SYSTEM  OF  FIRST  ORDER  ORDINARY  DIFFERENTIAL *
      ! EQUATIONS  WITH  INITIAL  CONDITIONS. IT ATTEMPTS TO KEEP THE GLOBAL *
      ! ERROR PROPORTIONAL TO  A  TOLERANCE  SPECIFIED  BY  THE  USER.  (THE *
      ! PROPORTIONALITY  DEPENDS  ON THE KIND OF ERROR CONTROL THAT IS USED, *
      ! AS WELL AS THE DIFFERENTIAL EQUATION AND THE RANGE OF INTEGRATION.)  *
      !                                                                      *
      !     VARIOUS OPTIONS ARE AVAILABLE TO THE USER,  INCLUDING  DIFFERENT *
      ! KINDS  OF  ERROR CONTROL, RESTRICTIONS ON STEP SIZES, AND INTERRUPTS *
      ! WHICH PERMIT THE USER TO EXAMINE THE STATE OF THE  CALCULATION  (AND *
      ! PERHAPS MAKE MODIFICATIONS) DURING INTERMEDIATE STAGES.              *
      !                                                                      *
      !     THE PROGRAM IS EFFICIENT FOR NON-STIFF SYSTEMS.  HOWEVER, A GOOD *
      ! VARIABLE-ORDER-ADAMS  METHOD  WILL PROBABLY BE MORE EFFICIENT IF THE *
      ! FUNCTION EVALUATIONS ARE VERY COSTLY.  SUCH A METHOD WOULD  ALSO  BE *
      ! MORE SUITABLE IF ONE WANTED TO OBTAIN A LARGE NUMBER OF INTERMEDIATE *
      ! SOLUTION VALUES BY INTERPOLATION, AS MIGHT BE THE CASE  FOR  EXAMPLE *
      ! WITH GRAPHICAL OUTPUT.                                               *
      !                                                                      *
      !                                    HULL-ENRIGHT-JACKSON   1/10/76    *
      !                                                                      *
      !***********************************************************************
      !                                                                      *
      !     USE - THE USER MUST SPECIFY EACH OF THE FOLLOWING                *
      !                                                                      *
      !     N  NUMBER OF EQUATIONS                                           *
      !                                                                      *
      !   FCN  NAME OF SUBROUTINE FOR EVALUATING FUNCTIONS - THE  SUBROUTINE *
      !           ITSELF MUST ALSO BE PROVIDED BY THE USER - IT SHOULD BE OF *
      !           THE FOLLOWING FORM                                         *
      !              SUBROUTINE FCN(N, X, Y, YPRIME)                         *
      !              INTEGER N                                               *
      !              REAL(DL) X, Y(N), YPRIME(N)                             *
      !                      *** ETC ***                                     *
      !           AND IT SHOULD EVALUATE YPRIME, GIVEN N, X AND Y            *
      !                                                                      *
      !     X  INDEPENDENT VARIABLE - INITIAL VALUE SUPPLIED BY USER         *
      !                                                                      *
      !     Y  DEPENDENT VARIABLE - INITIAL VALUES OF COMPONENTS Y(1), Y(2), *
      !           ..., Y(N) SUPPLIED BY USER                                 *
      !                                                                      *
      !  XEND  VALUE OF X TO WHICH INTEGRATION IS TO BE CARRIED OUT - IT MAY *
      !           BE LESS THAN THE INITIAL VALUE OF X                        *
      !                                                                      *
      !   TOL  TOLERANCE - THE SUBROUTINE ATTEMPTS TO CONTROL A NORM OF  THE *
      !           LOCAL  ERROR  IN  SUCH  A  WAY  THAT  THE  GLOBAL ERROR IS *
      !           PROPORTIONAL TO TOL. IN SOME PROBLEMS THERE WILL BE ENOUGH *
      !           DAMPING  OF  ERRORS, AS WELL AS SOME CANCELLATION, SO THAT *
      !           THE GLOBAL ERROR WILL BE LESS THAN TOL. ALTERNATIVELY, THE *
      !           CONTROL   CAN   BE  VIEWED  AS  ATTEMPTING  TO  PROVIDE  A *
      !           CALCULATED VALUE OF Y AT XEND WHICH IS THE EXACT  SOLUTION *
      !           TO  THE  PROBLEM Y' = F(X,Y) + E(X) WHERE THE NORM OF E(X) *
      !           IS PROPORTIONAL TO TOL.  (THE NORM  IS  A  MAX  NORM  WITH *
      !           WEIGHTS  THAT  DEPEND ON THE ERROR CONTROL STRATEGY CHOSEN *
      !           BY THE USER.  THE DEFAULT WEIGHT FOR THE K-TH COMPONENT IS *
      !           1/MAX(1,ABS(Y(K))),  WHICH THEREFORE PROVIDES A MIXTURE OF *
      !           ABSOLUTE AND RELATIVE ERROR CONTROL.)                      *
      !                                                                      *
      !   IND  INDICATOR - ON INITIAL ENTRY IND MUST BE SET EQUAL TO  EITHER *
      !           1  OR  2. IF THE USER DOES NOT WISH TO USE ANY OPTIONS, HE *
      !           SHOULD SET IND TO 1 - ALL THAT REMAINS FOR THE USER TO  DO *
      !           THEN  IS  TO  DECLARE C AND W, AND TO SPECIFY NW. THE USER *
      !           MAY ALSO  SELECT  VARIOUS  OPTIONS  ON  INITIAL  ENTRY  BY *
      !           SETTING IND = 2 AND INITIALIZING THE FIRST 9 COMPONENTS OF *
      !           C AS DESCRIBED IN THE NEXT SECTION.  HE MAY ALSO  RE-ENTER *
      !           THE  SUBROUTINE  WITH IND = 3 AS MENTIONED AGAIN BELOW. IN *
      !           ANY EVENT, THE SUBROUTINE RETURNS WITH IND EQUAL TO        *
      !              3 AFTER A NORMAL RETURN                                 *
      !              4, 5, OR 6 AFTER AN INTERRUPT (SEE OPTIONS C(8), C(9))  *
      !              -1, -2, OR -3 AFTER AN ERROR CONDITION (SEE BELOW)      *
      !                                                                      *
      !     C  COMMUNICATIONS VECTOR - THE DIMENSION MUST BE GREATER THAN OR *
      !           EQUAL TO 24, UNLESS OPTION C(1) = 4 OR 5 IS USED, IN WHICH *
      !           CASE THE DIMENSION MUST BE GREATER THAN OR EQUAL TO N+30   *
      !                                                                      *
      !    NW  FIRST DIMENSION OF WORKSPACE W -  MUST  BE  GREATER  THAN  OR *
      !           EQUAL TO N                                                 *
      !                                                                      *
      !     W  WORKSPACE MATRIX - FIRST DIMENSION MUST BE NW AND SECOND MUST *
      !           BE GREATER THAN OR EQUAL TO 9                              *
      !                                                                      *
      !     THE SUBROUTINE  WILL  NORMALLY  RETURN  WITH  IND  =  3,  HAVING *
      ! REPLACED THE INITIAL VALUES OF X AND Y WITH, RESPECTIVELY, THE VALUE *
      ! OF XEND AND AN APPROXIMATION TO Y AT XEND.  THE  SUBROUTINE  CAN  BE *
      ! CALLED  REPEATEDLY  WITH NEW VALUES OF XEND WITHOUT HAVING TO CHANGE *
      ! ANY OTHER ARGUMENT.  HOWEVER, CHANGES IN TOL, OR ANY OF THE  OPTIONS *
      ! DESCRIBED BELOW, MAY ALSO BE MADE ON SUCH A RE-ENTRY IF DESIRED.     *
      !                                                                      *
      !     THREE ERROR RETURNS ARE ALSO POSSIBLE, IN WHICH  CASE  X  AND  Y *
      ! WILL BE THE MOST RECENTLY ACCEPTED VALUES -                          *
      !     WITH IND = -3 THE SUBROUTINE WAS UNABLE  TO  SATISFY  THE  ERROR *
      !        REQUIREMENT  WITH A PARTICULAR STEP-SIZE THAT IS LESS THAN OR *
      !        EQUAL TO HMIN, WHICH MAY MEAN THAT TOL IS TOO SMALL           *
      !     WITH IND = -2 THE VALUE OF HMIN  IS  GREATER  THAN  HMAX,  WHICH *
      !        PROBABLY  MEANS  THAT THE REQUESTED TOL (WHICH IS USED IN THE *
      !        CALCULATION OF HMIN) IS TOO SMALL                             *
      !     WITH IND = -1 THE ALLOWED MAXIMUM NUMBER OF FCN EVALUATIONS  HAS *
      !        BEEN  EXCEEDED,  BUT  THIS  CAN ONLY OCCUR IF OPTION C(7), AS *
      !        DESCRIBED IN THE NEXT SECTION, HAS BEEN USED                  *
      !                                                                      *
      !     THERE ARE SEVERAL CIRCUMSTANCES THAT WILL CAUSE THE CALCULATIONS *
      ! TO  BE  TERMINATED,  ALONG WITH OUTPUT OF INFORMATION THAT WILL HELP *
      ! THE USER DETERMINE THE CAUSE OF  THE  TROUBLE.  THESE  CIRCUMSTANCES *
      ! INVOLVE  ENTRY WITH ILLEGAL OR INCONSISTENT VALUES OF THE ARGUMENTS, *
      ! SUCH AS ATTEMPTING A NORMAL  RE-ENTRY  WITHOUT  FIRST  CHANGING  THE *
      ! VALUE OF XEND, OR ATTEMPTING TO RE-ENTER WITH IND LESS THAN ZERO.    *
      !                                                                      *
      !***********************************************************************
      !                                                                      *
      !     OPTIONS - IF THE SUBROUTINE IS ENTERED WITH IND = 1, THE FIRST 9 *
      ! COMPONENTS OF THE COMMUNICATIONS VECTOR ARE INITIALIZED TO ZERO, AND *
      ! THE SUBROUTINE USES ONLY DEFAULT VALUES  FOR  EACH  OPTION.  IF  THE *
      ! SUBROUTINE  IS  ENTERED  WITH IND = 2, THE USER MUST SPECIFY EACH OF *
      ! THESE 9 COMPONENTS - NORMALLY HE WOULD FIRST SET THEM ALL  TO  ZERO, *
      ! AND  THEN  MAKE  NON-ZERO  THOSE  THAT  CORRESPOND TO THE PARTICULAR *
      ! OPTIONS HE WISHES TO SELECT. IN ANY EVENT, OPTIONS MAY BE CHANGED ON *
      ! RE-ENTRY  TO  THE  SUBROUTINE  -  BUT IF THE USER CHANGES ANY OF THE *
      ! OPTIONS, OR TOL, IN THE COURSE OF A CALCULATION HE SHOULD BE CAREFUL *
      ! ABOUT  HOW  SUCH CHANGES AFFECT THE SUBROUTINE - IT MAY BE BETTER TO *
      ! RESTART WITH IND = 1 OR 2. (COMPONENTS 10 TO 24 OF C ARE USED BY THE *
      ! PROGRAM  -  THE INFORMATION IS AVAILABLE TO THE USER, BUT SHOULD NOT *
      ! NORMALLY BE CHANGED BY HIM.)                                         *
      !                                                                      *
      !  C(1)  ERROR CONTROL INDICATOR - THE NORM OF THE LOCAL ERROR IS  THE *
      !           MAX  NORM  OF  THE  WEIGHTED  ERROR  ESTIMATE  VECTOR, THE *
      !           WEIGHTS BEING DETERMINED ACCORDING TO THE VALUE OF C(1) -  *
      !              IF C(1)=1 THE WEIGHTS ARE 1 (ABSOLUTE ERROR CONTROL)    *
      !              IF C(1)=2 THE WEIGHTS ARE 1/ABS(Y(K))  (RELATIVE  ERROR *
      !                 CONTROL)                                             *
      !              IF C(1)=3 THE  WEIGHTS  ARE  1/MAX(ABS(C(2)),ABS(Y(K))) *
      !                 (RELATIVE  ERROR  CONTROL,  UNLESS ABS(Y(K)) IS LESS *
      !                 THAN THE FLOOR VALUE, ABS(C(2)) )                    *
      !              IF C(1)=4 THE WEIGHTS ARE 1/MAX(ABS(C(K+30)),ABS(Y(K))) *
      !                 (HERE INDIVIDUAL FLOOR VALUES ARE USED)              *
      !              IF C(1)=5 THE WEIGHTS ARE 1/ABS(C(K+30))                *
      !              FOR ALL OTHER VALUES OF C(1), INCLUDING  C(1) = 0,  THE *
      !                 DEFAULT  VALUES  OF  THE  WEIGHTS  ARE  TAKEN  TO BE *
      !                 1/MAX(1,ABS(Y(K))), AS MENTIONED EARLIER             *
      !           (IN THE TWO CASES C(1) = 4 OR 5 THE USER MUST DECLARE  THE *
      !           DIMENSION OF C TO BE AT LEAST N+30 AND MUST INITIALIZE THE *
      !           COMPONENTS C(31), C(32), ..., C(N+30).)                    *
      !                                                                      *
      !  C(2)  FLOOR VALUE - USED WHEN THE INDICATOR C(1) HAS THE VALUE 3    *
      !                                                                      *
      !  C(3)  HMIN SPECIFICATION - IF NOT ZERO, THE SUBROUTINE CHOOSES HMIN *
      !           TO BE ABS(C(3)) - OTHERWISE IT USES THE DEFAULT VALUE      *
      !              10*MAX(DWARF,RREB*MAX(WEIGHTED NORM Y/TOL,ABS(X))),     *
      !           WHERE DWARF IS A VERY SMALL POSITIVE  MACHINE  NUMBER  AND *
      !           RREB IS THE RELATIVE ROUNDOFF ERROR BOUND                  *
      !                                                                      *
      !  C(4)  HSTART SPECIFICATION - IF NOT ZERO, THE SUBROUTINE  WILL  USE *
      !           AN  INITIAL  HMAG EQUAL TO ABS(C(4)), EXCEPT OF COURSE FOR *
      !           THE RESTRICTIONS IMPOSED BY HMIN AND HMAX  -  OTHERWISE IT *
      !           USES THE DEFAULT VALUE OF HMAX*(TOL)**(1/6)                *
      !                                                                      *
      !  C(5)  SCALE SPECIFICATION - THIS IS INTENDED TO BE A MEASURE OF THE *
      !           SCALE OF THE PROBLEM - LARGER VALUES OF SCALE TEND TO MAKE *
      !           THE METHOD MORE RELIABLE, FIRST  BY  POSSIBLY  RESTRICTING *
      !           HMAX  (AS  DESCRIBED  BELOW) AND SECOND, BY TIGHTENING THE *
      !           ACCEPTANCE REQUIREMENT - IF C(5) IS ZERO, A DEFAULT  VALUE *
      !           OF  1  IS  USED.  FOR  LINEAR  HOMOGENEOUS  PROBLEMS  WITH *
      !           CONSTANT COEFFICIENTS, AN APPROPRIATE VALUE FOR SCALE IS A *
      !           NORM  OF  THE  ASSOCIATED  MATRIX.  FOR OTHER PROBLEMS, AN *
      !           APPROXIMATION TO  AN  AVERAGE  VALUE  OF  A  NORM  OF  THE *
      !           JACOBIAN ALONG THE TRAJECTORY MAY BE APPROPRIATE           *
      !                                                                      *
      !  C(6)  HMAX SPECIFICATION - FOUR CASES ARE POSSIBLE                  *
      !           IF C(6).NE.0 AND C(5).NE.0, HMAX IS TAKEN TO BE            *
      !              MIN(ABS(C(6)),2/ABS(C(5)))                              *
      !           IF C(6).NE.0 AND C(5).EQ.0, HMAX IS TAKEN TO BE  ABS(C(6)) *
      !           IF C(6).EQ.0 AND C(5).NE.0, HMAX IS TAKEN TO BE            *
      !              2/ABS(C(5))                                             *
      !           IF C(6).EQ.0 AND C(5).EQ.0, HMAX IS GIVEN A DEFAULT  VALUE *
      !              OF 2                                                    *
      !                                                                      *
      !  C(7)  MAXIMUM NUMBER OF FUNCTION EVALUATIONS  -  IF  NOT  ZERO,  AN *
      !           ERROR  RETURN WITH IND = -1 WILL BE CAUSED WHEN THE NUMBER *
      !           OF FUNCTION EVALUATIONS EXCEEDS ABS(C(7))                  *
      !                                                                      *
      !  C(8)  INTERRUPT NUMBER  1  -  IF  NOT  ZERO,  THE  SUBROUTINE  WILL *
      !           INTERRUPT   THE  CALCULATIONS  AFTER  IT  HAS  CHOSEN  ITS *
      !           PRELIMINARY VALUE OF HMAG, AND JUST BEFORE CHOOSING HTRIAL *
      !           AND  XTRIAL  IN  PREPARATION FOR TAKING A STEP (HTRIAL MAY *
      !           DIFFER FROM HMAG IN SIGN, AND MAY  REQUIRE  ADJUSTMENT  IF *
      !           XEND  IS  NEAR) - THE SUBROUTINE RETURNS WITH IND = 4, AND *
      !           WILL RESUME CALCULATION AT THE POINT  OF  INTERRUPTION  IF *
      !           RE-ENTERED WITH IND = 4                                    *
      !                                                                      *
      !  C(9)  INTERRUPT NUMBER  2  -  IF  NOT  ZERO,  THE  SUBROUTINE  WILL *
      !           INTERRUPT   THE  CALCULATIONS  IMMEDIATELY  AFTER  IT  HAS *
      !           DECIDED WHETHER OR NOT TO ACCEPT THE RESULT  OF  THE  MOST *
      !           RECENT  TRIAL STEP, WITH IND = 5 IF IT PLANS TO ACCEPT, OR *
      !           IND = 6 IF IT PLANS TO REJECT -  Y(*)  IS  THE  PREVIOUSLY *
      !           ACCEPTED  RESULT, WHILE W(*,9) IS THE NEWLY COMPUTED TRIAL *
      !           VALUE, AND W(*,2) IS THE UNWEIGHTED ERROR ESTIMATE VECTOR. *
      !           THE  SUBROUTINE  WILL  RESUME CALCULATIONS AT THE POINT OF *
      !           INTERRUPTION ON RE-ENTRY WITH IND = 5 OR 6. (THE USER  MAY *
      !           CHANGE IND IN THIS CASE IF HE WISHES, FOR EXAMPLE TO FORCE *
      !           ACCEPTANCE OF A STEP THAT WOULD OTHERWISE BE REJECTED,  OR *
      !           VICE VERSA. HE CAN ALSO RESTART WITH IND = 1 OR 2.)        *
      !                                                                      *
      !***********************************************************************
      !                                                                      *
      !  SUMMARY OF THE COMPONENTS OF THE COMMUNICATIONS VECTOR              *
      !                                                                      *
      !     PRESCRIBED AT THE OPTION       DETERMINED BY THE PROGRAM         *
      !           OF THE USER                                                *
      !                                                                      *
      !                                    C(10) RREB(REL ROUNDOFF ERR BND)  *
      !     C(1) ERROR CONTROL INDICATOR   C(11) DWARF (VERY SMALL MACH NO)  *
      !     C(2) FLOOR VALUE               C(12) WEIGHTED NORM Y             *
      !     C(3) HMIN SPECIFICATION        C(13) HMIN                        *
      !     C(4) HSTART SPECIFICATION      C(14) HMAG                        *
      !     C(5) SCALE SPECIFICATION       C(15) SCALE                       *
      !     C(6) HMAX SPECIFICATION        C(16) HMAX                        *
      !     C(7) MAX NO OF FCN EVALS       C(17) XTRIAL                      *
      !     C(8) INTERRUPT NO 1            C(18) HTRIAL                      *
      !     C(9) INTERRUPT NO 2            C(19) EST                         *
      !                                    C(20) PREVIOUS XEND               *
      !                                    C(21) FLAG FOR XEND               *
      !                                    C(22) NO OF SUCCESSFUL STEPS      *
      !                                    C(23) NO OF SUCCESSIVE FAILURES   *
      !                                    C(24) NO OF FCN EVALS             *
      !                                                                      *
      !  IF C(1) = 4 OR 5, C(31), C(32), ... C(N+30) ARE FLOOR VALUES        *
      !                                                                      *
      !  RREB and DWARF are machine dependent constants currently set so     *
      !  that they should be appropriate for most machines.  However, it may *
      !  be appropriate to change them when this program is installed on a   *
      !  new machine.                            K.R.J.   3 Oct 1991.        *
      !                                                                      *
      !***********************************************************************
      !                                                                      *
      !  AN OVERVIEW OF THE PROGRAM                                          *
      !                                                                      *
      !     BEGIN INITIALIZATION, PARAMETER CHECKING, INTERRUPT RE-ENTRIES   *
      !  ......ABORT IF IND OUT OF RANGE 1 TO 6                              *
      !  .     CASES - INITIAL ENTRY, NORMAL RE-ENTRY, INTERRUPT RE-ENTRIES  *
      !  .     CASE 1 - INITIAL ENTRY (IND .EQ. 1 OR 2)                      *
      !  V........ABORT IF N.GT.NW OR TOL.LE.0                               *
      !  .        IF INITIAL ENTRY WITHOUT OPTIONS (IND .EQ. 1)              *
      !  .           SET C(1) TO C(9) EQUAL TO ZERO                          *
      !  .        ELSE INITIAL ENTRY WITH OPTIONS (IND .EQ. 2)               *
      !  .           MAKE C(1) TO C(9) NON-NEGATIVE                          *
      !  .           MAKE FLOOR VALUES NON-NEGATIVE IF THEY ARE TO BE USED   *
      !  .        END IF                                                     *
      !  .        INITIALIZE RREB, DWARF, PREV XEND, FLAG, COUNTS            *
      !  .     CASE 2 - NORMAL RE-ENTRY (IND .EQ. 3)                         *
      !  .........ABORT IF XEND REACHED, AND EITHER X CHANGED OR XEND NOT    *
      !  .        RE-INITIALIZE FLAG                                         *
      !  .     CASE 3 - RE-ENTRY FOLLOWING AN INTERRUPT (IND .EQ. 4 TO 6)    *
      !  V        TRANSFER CONTROL TO THE APPROPRIATE RE-ENTRY POINT.......  *
      !  .     END CASES                                                  .  *
      !  .  END INITIALIZATION, ETC.                                      .  *
      !  .                                                                V  *
      !  .  LOOP THROUGH THE FOLLOWING 4 STAGES, ONCE FOR EACH TRIAL STEP .  *
      !  .     STAGE 1 - PREPARE                                          .  *
      !***********ERROR RETURN (WITH IND=-1) IF NO OF FCN EVALS TOO GREAT .  *
      !  .        CALC SLOPE (ADDING 1 TO NO OF FCN EVALS) IF IND .NE. 6  .  *
      !  .        CALC HMIN, SCALE, HMAX                                  .  *
      !***********ERROR RETURN (WITH IND=-2) IF HMIN .GT. HMAX            .  *
      !  .        CALC PRELIMINARY HMAG                                   .  *
      !***********INTERRUPT NO 1 (WITH IND=4) IF REQUESTED.......RE-ENTRY.V  *
      !  .        CALC HMAG, XTRIAL AND HTRIAL                            .  *
      !  .     END STAGE 1                                                .  *
      !  V     STAGE 2 - CALC YTRIAL (ADDING 7 TO NO OF FCN EVALS)        .  *
      !  .     STAGE 3 - CALC THE ERROR ESTIMATE                          .  *
      !  .     STAGE 4 - MAKE DECISIONS                                   .  *
      !  .        SET IND=5 IF STEP ACCEPTABLE, ELSE SET IND=6            .  *
      !***********INTERRUPT NO 2 IF REQUESTED....................RE-ENTRY.V  *
      !  .        IF STEP ACCEPTED (IND .EQ. 5)                              *
      !  .           UPDATE X, Y FROM XTRIAL, YTRIAL                         *
      !  .           ADD 1 TO NO OF SUCCESSFUL STEPS                         *
      !  .           SET NO OF SUCCESSIVE FAILURES TO ZERO                   *
      !**************RETURN(WITH IND=3, XEND SAVED, FLAG SET) IF X .EQ. XEND *
      !  .        ELSE STEP NOT ACCEPTED (IND .EQ. 6)                        *
      !  .           ADD 1 TO NO OF SUCCESSIVE FAILURES                      *
      !**************ERROR RETURN (WITH IND=-3) IF HMAG .LE. HMIN            *
      !  .        END IF                                                     *
      !  .     END STAGE 4                                                   *
      !  .  END LOOP                                                         *
      !  .                                                                   *
      !  BEGIN ABORT ACTION                                                  *
      !     OUTPUT APPROPRIATE  MESSAGE  ABOUT  STOPPING  THE  CALCULATIONS, *
      !        ALONG WITH VALUES OF IND, N, NW, TOL, HMIN,  HMAX,  X,  XEND, *
      !        PREVIOUS XEND,  NO OF  SUCCESSFUL  STEPS,  NO  OF  SUCCESSIVE *
      !        FAILURES, NO OF FCN EVALS, AND THE COMPONENTS OF Y            *
      !     STOP                                                             *
      !  END ABORT ACTION                                                    *
      !                                                                      *
      !***********************************************************************
      !
      !     ******************************************************************
      !     * BEGIN INITIALIZATION, PARAMETER CHECKING, INTERRUPT RE-ENTRIES *
      !     ******************************************************************
      !
      !  ......ABORT IF IND OUT OF RANGE 1 TO 6
         IF (IND.LT.1 .OR. IND.GT.6) GO TO 500
      !
      !        CASES - INITIAL ENTRY, NORMAL RE-ENTRY, INTERRUPT RE-ENTRIES
      !  GO TO (5, 5, 45, 1111, 2222, 2222), IND
         if (ind==3) goto 45
         if (ind==4) goto 1111
         if (ind==5 .or. ind==6) goto 2222

      !        CASE 1 - INITIAL ENTRY (IND .EQ. 1 OR 2)
      !  .........ABORT IF N.GT.NW OR TOL.LE.0
            IF (N.GT.NW .OR. TOL.LE.0._dl) GO TO 500
            IF (IND.EQ. 2) GO TO 15
      !              INITIAL ENTRY WITHOUT OPTIONS (IND .EQ. 1)
      !              SET C(1) TO C(9) EQUAL TO 0
               DO K = 1, 9
                  C(K) = 0._dl
               end do
               GO TO 35
   15       CONTINUE
      !              INITIAL ENTRY WITH OPTIONS (IND .EQ. 2)
      !              MAKE C(1) TO C(9) NON-NEGATIVE
               DO K = 1, 9
                  C(K) = DABS(C(K))
               end do
      !              MAKE FLOOR VALUES NON-NEGATIVE IF THEY ARE TO BE USED
               IF (C(1).NE.4._dl .AND. C(1).NE.5._dl) GO TO 30
                  DO K = 1, N
                     C(K+30) = DABS(C(K+30))
                  end do
   30          CONTINUE
   35       CONTINUE
      !           INITIALIZE RREB, DWARF, PREV XEND, FLAG, COUNTS
            C(10) = 2._dl**(-56)
            C(11) = 1.D-35
      !           SET PREVIOUS XEND INITIALLY TO INITIAL VALUE OF X
            C(20) = X
            DO K = 21, 24
               C(K) = 0._dl
            end do
            GO TO 50
      !        CASE 2 - NORMAL RE-ENTRY (IND .EQ. 3)
      !  .........ABORT IF XEND REACHED, AND EITHER X CHANGED OR XEND NOT
   45       IF (C(21).NE.0._dl .AND. &
                              (X.NE.C(20) .OR. XEND.EQ.C(20))) GO TO 500
      !           RE-INITIALIZE FLAG
            C(21) = 0._dl
            GO TO 50
      !        CASE 3 - RE-ENTRY FOLLOWING AN INTERRUPT (IND .EQ. 4 TO 6)
      !           TRANSFER CONTROL TO THE APPROPRIATE RE-ENTRY POINT..........
      !           THIS HAS ALREADY BEEN HANDLED BY THE COMPUTED GO TO        .
      !        END CASES                                                     V
   50    CONTINUE
      !
      !     END INITIALIZATION, ETC.
      !
      !     ******************************************************************
      !     * LOOP THROUGH THE FOLLOWING 4 STAGES, ONCE FOR EACH TRIAL  STEP *
      !     * UNTIL THE OCCURRENCE OF ONE OF THE FOLLOWING                   *
      !     *    (A) THE NORMAL RETURN (WITH IND .EQ. 3) ON REACHING XEND IN *
      !     *        STAGE 4                                                 *
      !     *    (B) AN ERROR RETURN (WITH IND .LT. 0) IN STAGE 1 OR STAGE 4 *
      !     *    (C) AN INTERRUPT RETURN (WITH IND  .EQ.  4,  5  OR  6),  IF *
      !     *        REQUESTED, IN STAGE 1 OR STAGE 4                        *
      !     ******************************************************************
      !
99999 CONTINUE
      !
      !        ***************************************************************
      !        * STAGE 1 - PREPARE - DO CALCULATIONS OF  HMIN,  HMAX,  ETC., *
      !        * AND SOME PARAMETER  CHECKING,  AND  END  UP  WITH  SUITABLE *
      !        * VALUES OF HMAG, XTRIAL AND HTRIAL IN PREPARATION FOR TAKING *
      !        * AN INTEGRATION STEP.                                        *
      !        ***************************************************************
      !
      !***********ERROR RETURN (WITH IND=-1) IF NO OF FCN EVALS TOO GREAT
            IF (C(7).EQ.0._dl .OR. C(24).LT.C(7)) GO TO 100
               IND = -1
               RETURN
  100       CONTINUE
      !
      !           CALCULATE SLOPE (ADDING 1 TO NO OF FCN EVALS) IF IND .NE. 6
            IF (IND .EQ. 6) GO TO 105
               CALL FCN(N, X, Y, W(1,1))
               C(24) = C(24) + 1._dl
  105       CONTINUE
      !
      !           CALCULATE HMIN - USE DEFAULT UNLESS VALUE PRESCRIBED
            C(13) = C(3)
            IF (C(3) .NE. 0._dl) GO TO 165
      !              CALCULATE DEFAULT VALUE OF HMIN
      !              FIRST CALCULATE WEIGHTED NORM Y - C(12) - AS SPECIFIED
      !              BY THE ERROR CONTROL INDICATOR C(1)
               TEMP = 0._dl
               IF (C(1) .NE. 1) GO TO 115
      !                 ABSOLUTE ERROR CONTROL - WEIGHTS ARE 1
                  DO 110 K = 1, N
                     TEMP = DMAX1(TEMP, DABS(Y(K)))
  110             CONTINUE
                  C(12) = TEMP
                  GO TO 160
  115          IF (C(1) .NE. 2._dl) GO TO 120
      !                 RELATIVE ERROR CONTROL - WEIGHTS ARE 1/DABS(Y(K)) SO
      !                 WEIGHTED NORM Y IS 1
                  C(12) = 1._dl
                  GO TO 160
  120          IF (C(1) .NE. 3._dl) GO TO 130
      !                 WEIGHTS ARE 1/MAX(C(2),ABS(Y(K)))
                  DO 125 K = 1, N
                     TEMP = DMAX1(TEMP, DABS(Y(K))/C(2))
  125             CONTINUE
                  C(12) = DMIN1(TEMP, 1._dl)
                  GO TO 160
  130          IF (C(1) .NE. 4._dl) GO TO 140
      !                 WEIGHTS ARE 1/MAX(C(K+30),ABS(Y(K)))
                  DO 135 K = 1, N
                     TEMP = DMAX1(TEMP, DABS(Y(K))/C(K+30))
  135             CONTINUE
                  C(12) = DMIN1(TEMP, 1._dl)
                  GO TO 160
  140          IF (C(1) .NE. 5._dl) GO TO 150
      !                 WEIGHTS ARE 1/C(K+30)
                  DO 145 K = 1, N
                     TEMP = DMAX1(TEMP, DABS(Y(K))/C(K+30))
  145             CONTINUE
                  C(12) = TEMP
                  GO TO 160
  150          CONTINUE
      !                 DEFAULT CASE - WEIGHTS ARE 1/MAX(1,ABS(Y(K)))
                  DO 155 K = 1, N
                     TEMP = DMAX1(TEMP, DABS(Y(K)))
  155             CONTINUE
                  C(12) = DMIN1(TEMP, 1._dl)
  160          CONTINUE
               C(13) = 10._dl*DMAX1(C(11),C(10)*DMAX1(C(12)/TOL,DABS(X)))
  165       CONTINUE
      !
      !           CALCULATE SCALE - USE DEFAULT UNLESS VALUE PRESCRIBED
            C(15) = C(5)
            IF (C(5) .EQ. 0._dl) C(15) = 1._dl
      !
      !           CALCULATE HMAX - CONSIDER 4 CASES
      !           CASE 1 BOTH HMAX AND SCALE PRESCRIBED
               IF (C(6).NE.0._dl .AND. C(5).NE.0._dl) &
                                          C(16) = DMIN1(C(6), 2._dl/C(5))
      !           CASE 2 - HMAX PRESCRIBED, BUT SCALE NOT
               IF (C(6).NE.0._dl .AND. C(5).EQ.0._dl) C(16) = C(6)
      !           CASE 3 - HMAX NOT PRESCRIBED, BUT SCALE IS
               IF (C(6).EQ.0._dl .AND. C(5).NE.0._dl) C(16) = 2._dl/C(5)
      !           CASE 4 - NEITHER HMAX NOR SCALE IS PROVIDED
               IF (C(6).EQ.0._dl .AND. C(5).EQ.0._dl) C(16) = 2._dl
      !
      !***********ERROR RETURN (WITH IND=-2) IF HMIN .GT. HMAX
            IF (C(13) .LE. C(16)) GO TO 170
               IND = -2
               RETURN
  170       CONTINUE
      !
      !           CALCULATE PRELIMINARY HMAG - CONSIDER 3 CASES
            IF (IND .GT. 2) GO TO 175
      !           CASE 1 - INITIAL ENTRY - USE PRESCRIBED VALUE OF HSTART, IF
      !              ANY, ELSE DEFAULT
               C(14) = C(4)
               IF (C(4) .EQ. 0._dl) C(14) = C(16)*TOL**(1._dl/6._dl)
               GO TO 185
  175       IF (C(23) .GT. 1._dl) GO TO 180
      !           CASE 2 - AFTER A SUCCESSFUL STEP, OR AT MOST  ONE  FAILURE,
      !              USE MIN(2, .9*(TOL/EST)**(1/6))*HMAG, BUT AVOID POSSIBLE
      !              OVERFLOW. THEN AVOID REDUCTION BY MORE THAN HALF.
               TEMP = 2._dl*C(14)
               IF (TOL .LT. (2._dl/.9_dl)**6*C(19)) &
                                  TEMP = .9_dl*(TOL/C(19))**(1._dl/6._dl)*C(14)
               C(14) = DMAX1(TEMP, .5_dl*C(14))
               GO TO 185
  180       CONTINUE
      !           CASE 3 - AFTER TWO OR MORE SUCCESSIVE FAILURES
               C(14) = .5_dl*C(14)
  185       CONTINUE
      !
      !           CHECK AGAINST HMAX
            C(14) = DMIN1(C(14), C(16))
      !
      !           CHECK AGAINST HMIN
            C(14) = DMAX1(C(14), C(13))
      !
      !***********INTERRUPT NO 1 (WITH IND=4) IF REQUESTED
            IF (C(8) .EQ. 0._dl) GO TO 1111
               IND = 4
               RETURN
      !           RESUME HERE ON RE-ENTRY WITH IND .EQ. 4   ........RE-ENTRY..
 1111       CONTINUE
      !
      !           CALCULATE HMAG, XTRIAL - DEPENDING ON PRELIMINARY HMAG, XEND
            IF (C(14) .GE. DABS(XEND - X)) GO TO 190
      !              DO NOT STEP MORE THAN HALF WAY TO XEND
               C(14) = DMIN1(C(14), .5_dl*DABS(XEND - X))
               C(17) = X + DSIGN(C(14), XEND - X)
               GO TO 195
  190       CONTINUE
      !              HIT XEND EXACTLY
               C(14) = DABS(XEND - X)
               C(17) = XEND
  195       CONTINUE
      !
      !           CALCULATE HTRIAL
            C(18) = C(17) - X
      !
      !        END STAGE 1
      !
      !        ***************************************************************
      !        * STAGE 2 - CALCULATE YTRIAL (ADDING 7 TO NO OF  FCN  EVALS). *
      !        * W(*,2), ... W(*,8)  HOLD  INTERMEDIATE  RESULTS  NEEDED  IN *
      !        * STAGE 3. W(*,9) IS TEMPORARY STORAGE UNTIL FINALLY IT HOLDS *
      !        * YTRIAL.                                                     *
      !        ***************************************************************
      !
            TEMP = C(18)/1398169080000._dl
      !
            DO 200 K = 1, N
               W(K,9) = Y(K) + TEMP*W(K,1)*233028180000._dl
  200       CONTINUE
            CALL FCN(N, X + C(18)/6._dl, W(1,9), W(1,2))
      !
            DO 205 K = 1, N
               W(K,9) = Y(K) + TEMP*(   W(K,1)*74569017600._dl &
                                      + W(K,2)*298276070400._dl  )
  205       CONTINUE
            CALL FCN(N, X + C(18)*(4._dl/15._dl), W(1,9), W(1,3))
      !
            DO 210 K = 1, N
               W(K,9) = Y(K) + TEMP*(   W(K,1)*1165140900000._dl &
                                      - W(K,2)*3728450880000._dl &
                                      + W(K,3)*3495422700000._dl )
  210       CONTINUE
            CALL FCN(N, X + C(18)*(2._dl/3._dl), W(1,9), W(1,4))
      !
            DO 215 K = 1, N
               W(K,9) = Y(K) + TEMP*( - W(K,1)*3604654659375._dl &
                                      + W(K,2)*12816549900000._dl &
                                      - W(K,3)*9284716546875._dl &
                                      + W(K,4)*1237962206250._dl )
  215       CONTINUE
            CALL FCN(N, X + C(18)*(5._dl/6._dl), W(1,9), W(1,5))
      !
            DO 220 K = 1, N
               W(K,9) = Y(K) + TEMP*(   W(K,1)*3355605792000._dl &
                                      - W(K,2)*11185352640000._dl &
                                      + W(K,3)*9172628850000._dl &
                                      - W(K,4)*427218330000._dl &
                                      + W(K,5)*482505408000._dl  )
  220       CONTINUE
            CALL FCN(N, X + C(18), W(1,9), W(1,6))
      !
            DO 225 K = 1, N
               W(K,9) = Y(K) + TEMP*( - W(K,1)*770204740536._dl &
                                      + W(K,2)*2311639545600._dl &
                                      - W(K,3)*1322092233000._dl &
                                      - W(K,4)*453006781920._dl &
                                      + W(K,5)*326875481856._dl  )
  225       CONTINUE
            CALL FCN(N, X + C(18)/15._dl, W(1,9), W(1,7))
      !
            DO 230 K = 1, N
               W(K,9) = Y(K) + TEMP*(   W(K,1)*2845924389000._dl &
                                      - W(K,2)*9754668000000._dl &
                                      + W(K,3)*7897110375000._dl &
                                      - W(K,4)*192082660000._dl &
                                      + W(K,5)*400298976000._dl &
                                      + W(K,7)*201586000000._dl  )
  230       CONTINUE
            CALL FCN(N, X + C(18), W(1,9), W(1,8))
      !
      !           CALCULATE YTRIAL, THE EXTRAPOLATED APPROXIMATION AND STORE
      !              IN W(*,9)
            DO 235 K = 1, N
               W(K,9) = Y(K) + TEMP*(   W(K,1)*104862681000._dl &
                                      + W(K,3)*545186250000._dl &
                                      + W(K,4)*446637345000._dl &
                                      + W(K,5)*188806464000._dl &
                                      + W(K,7)*15076875000._dl &
                                      + W(K,8)*97599465000._dl   )
  235       CONTINUE
      !
      !           ADD 7 TO THE NO OF FCN EVALS
            C(24) = C(24) + 7._dl
      !
      !        END STAGE 2
      !
      !        ***************************************************************
      !        * STAGE 3 - CALCULATE THE ERROR ESTIMATE EST. FIRST CALCULATE *
      !        * THE  UNWEIGHTED  ABSOLUTE  ERROR  ESTIMATE VECTOR (PER UNIT *
      !        * STEP) FOR THE UNEXTRAPOLATED APPROXIMATION AND STORE IT  IN *
      !        * W(*,2).  THEN  CALCULATE THE WEIGHTED MAX NORM OF W(*,2) AS *
      !        * SPECIFIED BY THE ERROR  CONTROL  INDICATOR  C(1).  FINALLY, *
      !        * MODIFY  THIS RESULT TO PRODUCE EST, THE ERROR ESTIMATE (PER *
      !        * UNIT STEP) FOR THE EXTRAPOLATED APPROXIMATION YTRIAL.       *
      !        ***************************************************************
      !
      !           CALCULATE THE UNWEIGHTED ABSOLUTE ERROR ESTIMATE VECTOR
            DO 300 K = 1, N
               W(K,2) = (   W(K,1)*8738556750._dl &
                          + W(K,3)*9735468750._dl &
                          - W(K,4)*9709507500._dl &
                          + W(K,5)*8582112000._dl &
                          + W(K,6)*95329710000._dl &
                          - W(K,7)*15076875000._dl &
                          - W(K,8)*97599465000._dl)/1398169080000._dl
  300       CONTINUE
      !
      !           CALCULATE THE WEIGHTED MAX NORM OF W(*,2) AS SPECIFIED BY
      !           THE ERROR CONTROL INDICATOR C(1)
            TEMP = 0._dl
            IF (C(1) .NE. 1) GO TO 310
      !              ABSOLUTE ERROR CONTROL
               DO 305 K = 1, N
                  TEMP = DMAX1(TEMP,DABS(W(K,2)))
  305          CONTINUE
               GO TO 360
  310       IF (C(1) .NE. 2._dl) GO TO 320
      !              RELATIVE ERROR CONTROL
               DO 315 K = 1, N
                  TEMP = DMAX1(TEMP, DABS(W(K,2)/Y(K)))
  315          CONTINUE
               GO TO 360
  320       IF (C(1) .NE. 3._dl) GO TO 330
      !              WEIGHTS ARE 1/MAX(C(2),ABS(Y(K)))
               DO 325 K = 1, N
                  TEMP = DMAX1(TEMP, DABS(W(K,2)) &
                                   / DMAX1(C(2), DABS(Y(K))) )
  325          CONTINUE
               GO TO 360
  330       IF (C(1) .NE. 4._dl) GO TO 340
      !              WEIGHTS ARE 1/MAX(C(K+30),ABS(Y(K)))
               DO 335 K = 1, N
                  TEMP = DMAX1(TEMP, DABS(W(K,2)) &
                                   / DMAX1(C(K+30), DABS(Y(K))) )
  335          CONTINUE
               GO TO 360
  340       IF (C(1) .NE. 5._dl) GO TO 350
      !              WEIGHTS ARE 1/C(K+30)
               DO 345 K = 1, N
                  TEMP = DMAX1(TEMP, DABS(W(K,2)/C(K+30)))
  345          CONTINUE
               GO TO 360
  350       CONTINUE
      !              DEFAULT CASE - WEIGHTS ARE 1/MAX(1,ABS(Y(K)))
               DO 355 K = 1, N
                  TEMP = DMAX1(TEMP, DABS(W(K,2)) &
                                   / DMAX1(1._dl, DABS(Y(K))) )
  355          CONTINUE
  360       CONTINUE
      !
      !           CALCULATE EST - (THE WEIGHTED MAX NORM OF W(*,2))*HMAG*SCALE
      !              - EST IS INTENDED TO BE A MEASURE OF THE ERROR  PER  UNIT
      !              STEP IN YTRIAL
            C(19) = TEMP*C(14)*C(15)
      !
      !        END STAGE 3
      !
      !        ***************************************************************
      !        * STAGE 4 - MAKE DECISIONS.                                   *
      !        ***************************************************************
      !
      !           SET IND=5 IF STEP ACCEPTABLE, ELSE SET IND=6
            IND = 5
            IF (C(19) .GT. TOL) IND = 6
      !
      !***********INTERRUPT NO 2 IF REQUESTED
            IF (C(9) .EQ. 0._dl) GO TO 2222
               RETURN
      !           RESUME HERE ON RE-ENTRY WITH IND .EQ. 5 OR 6   ...RE-ENTRY..
 2222       CONTINUE
      !
            IF (IND .EQ. 6) GO TO 410
      !              STEP ACCEPTED (IND .EQ. 5), SO UPDATE X, Y FROM XTRIAL,
      !                 YTRIAL, ADD 1 TO THE NO OF SUCCESSFUL STEPS, AND SET
      !                 THE NO OF SUCCESSIVE FAILURES TO ZERO
               X = C(17)
               DO 400 K = 1, N
                  Y(K) = W(K,9)
  400          CONTINUE
               C(22) = C(22) + 1._dl
               C(23) = 0._dl
      !**************RETURN(WITH IND=3, XEND SAVED, FLAG SET) IF X .EQ. XEND
               IF (X .NE. XEND) GO TO 405
                  IND = 3
                  C(20) = XEND
                  C(21) = 1._dl
                  RETURN
  405          CONTINUE
               GO TO 420
  410       CONTINUE
      !              STEP NOT ACCEPTED (IND .EQ. 6), SO ADD 1 TO THE NO OF
      !                 SUCCESSIVE FAILURES
               C(23) = C(23) + 1._dl
      !**************ERROR RETURN (WITH IND=-3) IF HMAG .LE. HMIN
               IF (C(14) .GT. C(13)) GO TO 415
                  IND = -3
                  RETURN
  415          CONTINUE
  420       CONTINUE
      !
      !        END STAGE 4
      !
      GO TO 99999
      !     END LOOP
      !
      !  BEGIN ABORT ACTION
  500 CONTINUE
      !
      WRITE(6,505) IND, TOL, X, &
                   N, C(13), XEND, &
                   NW, C(16), C(20), &
                   C(22), &
                   C(23), &
                   C(24), &
                   (Y(K), K = 1, N)
  505 FORMAT( /// 1H0, 58HCOMPUTATION STOPPED IN DVERK WITH THE FOLLOWING VALUES    &
                / 1H0, 5HIND =, I4, 5X, 6HTOL  =, 1PD13.6, 5X, 11HX         =, 1PD22.15 &
                / 1H , 5HN   =, I4, 5X, 6HHMIN =, 1PD13.6, 5X, 11HXEND      =, 1PD22.15 &
                / 1H , 5HNW  =, I4, 5X, 6HHMAX =, 1PD13.6, 5X, 11HPREV XEND =, 1PD22.15 &
                / 1H0, 14X, 27HNO OF SUCCESSFUL STEPS    =, 0PF8.0 &
                / 1H , 14X, 27HNO OF SUCCESSIVE FAILURES =, 0PF8.0 &
                / 1H , 14X, 27HNO OF FUNCTION EVALS      =, 0PF8.0 &
                / 1H0, 23HTHE COMPONENTS OF Y ARE &
                // (1H , 1P5D24.15)                                           )
      !
      STOP
      !
      !  END ABORT ACTION
      !
      END

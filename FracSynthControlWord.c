#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef GLOBAL_RTN
#define GLOBAL_RTN
#endif
#ifdef DEBUG_PRINT
#define DEBUGPRINT(interest, globalFlag, args) \
    {                                          \
        if (interest <= globalFlag)            \
            printf args;                       \
    }
#else
#define DEBUGPRINT(interest, globalFlag, args)
#endif

#ifndef LOCAL
#define LOCAL static
#endif

#define DP_INFO 4

/*---------------------
 * Define the size of each field
 */

#define CONTROL_MDIV_BITS 3     /* Denominator of correction term                       */
#define CONTROL_NDIV_BITS 3     /* Numerator of correction term                         */
#define CONTROL_POSTDIV_BITS 5  /* Post-Divider                                         */
#define CONTROL_MFG_BITS 3      /* Must Be Zero                                         */
#define CONTROL_P_BITS 4        /* Integer part of fractional frequency                 */
#define CONTROL_QPM1_BITS 5     /* Value for Q(p-1) term of fractional frequency        */
#define CONTROL_QP_BITS 5       /* Value for Q(p) term of fractional frequency          */
#define CONTROL_PREAMBLE_BITS 4 /* Must Be Zero                                         */

#define MAX_CORRECTION_RATIO (17. / 14.) /* Maximum value for correction term.                   */
#define MAX_VCO_FREQ 729.0               /* Maximum frequency for voltage-controlled oscillator  */
#define MIN_VCO_FREQ 540.0               /* Minimum frequency for voltage-controlled oscillator  */
#define MIN_P_VALUE 17                   /* Minimum val for integer part of fractional frequency */
#define MAX_FRAC_DIVISOR 31              /* Maximum divisor for fractional frequency value       */

#define NUM_POST_DIVIDES 31     /* Number of unique post-divider values                 */
#define NUM_POST_DIVIDE_VALS 32 /* Number of post-divider codes                         */
#define NUM_CORRECTIONS 23      /* Number of valid correction values (plus one)         */
#define NUM_CORRECTION_VALS 8   /* Number of correction value codes                     */

#define MAX_ERROR 100.0       /* Artifical error maximum                              */
#define ZERO_THRESHOLD 1.0e-9 /* Floating point threshold for zero detection          */
#define NSEC_PER_SEC 1.e9

typedef struct
{                      /* PostDivideStruct                                      */
    double Divisor;    /*   Divisor value                                       */
    unsigned int Code; /*   Field code for post-divide value                    */
} PostDivideStruct;
/*---------------------
 * Correction Factor Structure
 */
typedef struct
{                        /* CorrectionStruct                                      */
    double Ratio;        /*   Correction factor ratio value                       */
    unsigned short nDiv; /*   Field code for dividend (N)                         */
    unsigned short mDiv; /*   Field code for divisor (M)                          */
} CorrectionStruct;
/*---------------------
 * Correction Factor Component Value Structure
 */
typedef struct
{              /* CorrectionValStruct                                   */
    int Value; /*   Field code value                                    */
    int Class; /*   Component class (1 or 2)                            */
} CorrectionValStruct;
/*---------------------
 * Set up a bitmask for each field
 */

#define CONTROL_MDIV_MASK (((1 << CONTROL_MDIV_BITS) - 1) << CONTROL_MDIV_SHIFT)
#define CONTROL_NDIV_MASK (((1 << CONTROL_NDIV_BITS) - 1) << CONTROL_NDIV_SHIFT)
#define CONTROL_POSTDIV_MASK (((1 << CONTROL_POSTDIV_BITS) - 1) << CONTROL_POSTDIV_SHIFT)
#define CONTROL_MFG_MASK (((1 << CONTROL_MFG_BITS) - 1) << CONTROL_MFG_SHIFT)
#define CONTROL_P_MASK (((1 << CONTROL_P_BITS) - 1) << CONTROL_P_SHIFT)
#define CONTROL_QPM1_MASK (((1 << CONTROL_QPM1_BITS) - 1) << CONTROL_QPM1_SHIFT)
#define CONTROL_QP_MASK (((1 << CONTROL_QP_BITS) - 1) << CONTROL_QP_SHIFT)
#define CONTROL_PREAMBLE_MASK (((1 << CONTROL_PREAMBLE_BITS) - 1) << CONTROL_PREAMBLE_SHIFT)

/*---------------------
 * Define the field values for the correction factor components
 */
#define CORRECTION_DIV_14 5 /* Numerator or denominator of 14                        */
#define CORRECTION_DIV_15 7 /* Numerator or denominator of 15                        */
#define CORRECTION_DIV_16 1 /* Numerator or denominator of 16                        */
#define CORRECTION_DIV_17 3 /* Numerator or denominator of 17                        */
#define CORRECTION_DIV_18 2 /* Numerator or denominator of 18                        */
#define CORRECTION_DIV_31 4 /* Numerator or denominator of 31                        */
#define CORRECTION_DIV_32 6 /* Numerator or denominator of 32                        */

/*---------------------
 * Define the offset of each field
 */

#define CONTROL_MDIV_SHIFT 0      /* Denominator of correction term                       */
#define CONTROL_NDIV_SHIFT 3      /* Numerator of correction term                         */
#define CONTROL_POSTDIV_SHIFT 6   /* Post-Divider                                         */
#define CONTROL_MFG_SHIFT 11      /* Must Be Zero                                         */
#define CONTROL_P_SHIFT 14        /* Integer part of fractional frequency                 */
#define CONTROL_QPM1_SHIFT 18     /* Value for Q(p-1) term of fractional frequency        */
#define CONTROL_QP_SHIFT 23       /* Value for Q(p) term of fractional frequency          */
#define CONTROL_PREAMBLE_SHIFT 28 /* Must Be Zero                                         */

/*---------------------
 * Correction Factor Table:
 *---------------------
 *   Contains legal correction factor ratios and their control-word field codes.
 */

LOCAL const CorrectionStruct CorrectionList[NUM_CORRECTIONS] = {
    { (1.0), CORRECTION_DIV_14, CORRECTION_DIV_14 },
    { (14. / 18.), CORRECTION_DIV_14, CORRECTION_DIV_18 },
    { (14. / 17.), CORRECTION_DIV_14, CORRECTION_DIV_17 },
    { (15. / 18.), CORRECTION_DIV_15, CORRECTION_DIV_18 },
    { (14. / 16.), CORRECTION_DIV_14, CORRECTION_DIV_16 },
    { (15. / 17.), CORRECTION_DIV_15, CORRECTION_DIV_17 },
    { (16. / 18.), CORRECTION_DIV_18, CORRECTION_DIV_18 },
    { (14. / 15.), CORRECTION_DIV_14, CORRECTION_DIV_15 },
    { (15. / 16.), CORRECTION_DIV_15, CORRECTION_DIV_16 },
    { (16. / 17.), CORRECTION_DIV_16, CORRECTION_DIV_17 },
    { (17. / 18.), CORRECTION_DIV_17, CORRECTION_DIV_18 },
    { (31. / 32.), CORRECTION_DIV_31, CORRECTION_DIV_32 },
    { (1.0), CORRECTION_DIV_14, CORRECTION_DIV_14 },
    { (32. / 31.), CORRECTION_DIV_32, CORRECTION_DIV_31 },
    { (18. / 17.), CORRECTION_DIV_18, CORRECTION_DIV_17 },
    { (17. / 16.), CORRECTION_DIV_17, CORRECTION_DIV_16 },
    { (16. / 15.), CORRECTION_DIV_16, CORRECTION_DIV_15 },
    { (15. / 14.), CORRECTION_DIV_15, CORRECTION_DIV_14 },
    { (18. / 16.), CORRECTION_DIV_18, CORRECTION_DIV_16 },
    { (17. / 15.), CORRECTION_DIV_17, CORRECTION_DIV_15 },
    { (16. / 14.), CORRECTION_DIV_16, CORRECTION_DIV_14 },
    { (18. / 15.), CORRECTION_DIV_18, CORRECTION_DIV_15 },
    { (17. / 14.), CORRECTION_DIV_17, CORRECTION_DIV_14 }
}; /*Correction List*/

/*---------------------
 * Fractional Synthesizer Component Structure
 */
typedef struct
{                         /* FracSynthComponents                                   */
    double Error;         /*   Deviation from desired output frequency             */
    double EffectiveFreq; /*   Actual frequency generated by these components      */
    int PostDivIndex;     /*   Index into post-divider table                       */
    int CorrectionIndex;  /*   Index into correction factor table                  */
    int P;                /*   Integer part of fractional frequency                */
    int Qp;               /*   Value of Q(p) term of fractional frequency          */
    int Qpm1;             /*   Value of Q(p-1) term of fractional frequency        */
} FracSynthComponents;

LOCAL const PostDivideStruct PostDivideList[NUM_POST_DIVIDES] = {
    { 1.0, 0x00 }, { 2.0, 0x02 }, { 3.0, 0x03 }, { 4.0, 0x04 }, { 5.0, 0x05 }, { 6.0, 0x06 }, { 7.0, 0x07 }, { 8.0, 0x08 }, { 9.0, 0x09 }, { 10.0, 0x0A }, { 11.0, 0x0B }, { 12.0, 0x0C }, { 13.0, 0x0D }, { 14.0, 0x0E }, { 15.0, 0x0F }, { 16.0, 0x10 }, { 18.0, 0x11 }, { 20.0, 0x12 }, { 22.0, 0x13 }, { 24.0, 0x14 }, { 26.0, 0x15 }, { 28.0, 0x16 }, { 30.0, 0x17 }, { 32.0, 0x18 }, { 36.0, 0x19 }, { 40.0, 0x1A }, { 44.0, 0x1B }, { 48.0, 0x1C }, { 52.0, 0x1D }, { 56.0, 0x1E }, { 60.0, 0x1F }
};

/*---------------------
 * Post-Divider Value Table:
 *---------------------
 *   Translates the post-divider field code into the actual post-divider value.
 */

LOCAL int PostDivideValList[NUM_POST_DIVIDE_VALS] = {
    1, 3, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
    16, 18, 20, 22, 24, 26, 28, 30, 32, 36, 40, 44, 48, 52, 56, 60
};

/*---------------------
 * Correction Factor Value Table:
 *---------------------
 *   Translates correction factor field codes into the actual correction value term values
 */

LOCAL const CorrectionValStruct CorrectionValList[NUM_CORRECTION_VALS] = {
    { 16, 1 }, { 16, 1 }, { 18, 1 }, { 17, 1 }, { 31, 2 }, { 14, 1 }, { 32, 2 }, { 15, 1 }
}; /* CorrectionValList*/

unsigned int CalcControlWord(FracSynthComponents Best, double DesiredFreq, double ReferenceFreq);

GLOBAL_RTN
unsigned int FracSynthControlWord(

    /**********************************************************************************************/
    /*  Parameter Declarations                                                                    */
    /**********************************************************************************************/

    double DesiredFreq,   /* Desired output frequency                         */
    double ReferenceFreq) /* SY87739L input reference frequency               */

{
    /**********************************************************************************************/
    /*  Local Variables                                                                           */
    /**********************************************************************************************/

    FracSynthComponents Best;         /* Best overall parameters seen so far              */
    FracSynthComponents BestFracFreq; /* Best fractional frequency parameters seen so far */
    unsigned int ControlWord;         /* Computed control word value                      */
    double CorrectionErr;             /* Error value for current correction factor        */
    double CorrectionFreq;            /* Frequency after applying the correction factor   */
    int CorrectionIndex;              /* Index to correction factor list                  */
    double d;                         /* FP representation of the current denominator     */
    double EffectiveFreq = 0.0;       /* Effective output frequency                       */
    double FracFreqErr;               /* Best fractional frequency error seen so far      */
    double FractionalFreq;            /* Desired fractional frequency                     */
    double FreqErr;                   /* Error value for the desired frequency            */
    int i, j, k, n;                   /* Loop indicies                                    */
    int Numerator = 1;                /* Best numerator found for current denominator     */
    double OldCorrectionErr;          /* Error value for previous correction factor       */
    int p;                            /* Integer part of fractional frequency             */
    int p1;                           /* Numerator of fractional frequency ratio          */
    double PostDivide;                /* Current post divider value                       */
    int Qp;                           /* Number of "P" pulses in the frac. div. circuit   */
    int Qpm1;                         /* Number of "P-1" pulses in the frac. div. circuit */
    double TestFreq;                  /* Fractional frequency to try in correction loop   */
    double VcoFreq;                   /* Desired VCO frequency                            */

    /*-------------------------------------------------------------------------------------------------
    |* NOTES:
    |* Debug printing uses the SLAC "debugPrint" utility. Debug levels are defined as:
    |*    0: DP_NONE -  No debug output is produced
    |*    1: DP_FATAL - Print messages only for fatal errors
    |*    2: DP_ERROR - Print fatal and non-fatal error messages
    |*    3: DP_WARN  - Print warning, fatal and non-fatal error messages
    |*    4: DP_INFO  - Print informational messages, warnings, and errors
    |*    5: DP_DEBUG - Print all of the above plus special messages inserted for code debugging.
    |*
    \**************************************************************************************************/
    /**********************************************************************************************/
    /*  Code                                                                                      */
    /**********************************************************************************************/

    /*---------------------
     * Initialize the "Best Fractional Frequency So Far" parameters
     */
    BestFracFreq.Error = MAX_ERROR;   /* Deviation from desired output frequency        */
    BestFracFreq.EffectiveFreq = 0.0; /* Actual frequency generated by these components */
    BestFracFreq.PostDivIndex = 0;    /* Index into post-divider table                  */
    BestFracFreq.CorrectionIndex = 0; /* Index into correction factor table             */
    BestFracFreq.P = 0;               /* Integer part of fractional frequency           */
    BestFracFreq.Qp = 1;              /* Value of Q(p) term of fractional frequency     */
    BestFracFreq.Qpm1 = 1;            /* Value of Q(p-1) term of fractional frequency   */
    Best = BestFracFreq;              /* Best overall parameters                        */

    /*---------------------
     * Post-Divider Loop:
     *---------------------
     *  Check all post-divider values that would put the VCO frequency within
     *  the allowable range (540 - 729 MHz).
     */
    Best.Error = MAX_ERROR;

    printf("Desired Frequency\tEffective Frequency\tControl Word\tError\t\tQp\tQpm1\tP\tn\tm\n");

    for (i = 0; i < NUM_POST_DIVIDES; i++) {
        PostDivide = PostDivideList[i].Divisor;

        /*---------------------
         * Compute the VCO frequency and make sure it is within the allowable range.
         */
        VcoFreq = DesiredFreq * PostDivide;

        if (VcoFreq >= MAX_VCO_FREQ)
            break;

        if (VcoFreq >= MIN_VCO_FREQ) {

            /*---------------------
             * We have a VCO frequency inside the allowable range.
             * Now compute the desired fractional frequency.
             *
             * The desired fractional frequency is derived by dividing the VCO frequency by the
             * reference frequency.  The actual fractional-N frequency is created by a fractional-N
             * P/P-1 divider which attempts to represent the desired fractional frequency as a
             * rational fraction (with a divisor less than 32) of the reference frequency.
             *
             * Note that the actual fractional-N frequency may not be exactly identical to the
             * desired fractional frequency.
             */
            FractionalFreq = VcoFreq / ReferenceFreq;
            BestFracFreq.Error = MAX_ERROR;

            /*---------------------
             * Divisor Loop:
             *---------------------
             *  Search for a divisor that will produce the closest fractional-N P/P-1 frequency to
             *  the desired fractional frequency.  If we are lucky, we will find a number that
             *  exactly divides the desired fractional frequency and we can stop the search here.
             *  If not, we will have to try different numerator and correction factor combinations
             *  in order to find the closest fit.
             */
            for (j = 1; j <= MAX_FRAC_DIVISOR; j++) {
                d = j;               /* Floating point representation of denominator  */
                CorrectionIndex = 0; /* No correction factor                          */

                /*---------------------
                 * Compute the integer and fractional parts of the fractional-N frequency
                 */
                p1 = (FractionalFreq * j);
                p = (p1 / j) + 1;
                Qpm1 = j - (p1 % j);
                Qp = j - Qpm1;

                /*---------------------
                 * Compute the actual frequency generated by the fractional-N P/P-1 divider
                 * for these values.  Also compute the error between the actual fractional
                 * frequency and the desired fractional frequency.
                 */
                EffectiveFreq = (double)p - ((double)Qpm1 / (double)(Qp + Qpm1));
                FracFreqErr = fabs(FractionalFreq - EffectiveFreq);

                /*---------------------
                 * If the current divisor does not exactly divide the desired fractional frequency,
                 * search for a numerator and correction-factor pair that will come closest to
                 * generating the desired fractional frequency.
                 */
                if (FracFreqErr >= ZERO_THRESHOLD) {
                    FracFreqErr = MAX_ERROR;

                    /*---------------------
                     * Numerator Loop:
                     *---------------------
                     *  Search for the numerator/correction-factor pair with the lowest error.
                     */
                    for (n = 1; n <= j; n++) {
                        TestFreq = (double)p - ((double)n / d);
                        OldCorrectionErr = MAX_ERROR;

                        /*---------------------
                         * Correction Factor Loop:
                         *---------------------
                         *  Search for the correction factor that produces the lowest error
                         *  with the current numerator/denominator pair.  Note that since the
                         *  correction factor list is arranged in ascending order, we can
                         *  stop the search as soon as the error value stops decreasing.
                         */
                        for (k = 1; k < NUM_CORRECTIONS; k++) {
                            CorrectionFreq = CorrectionList[k].Ratio * TestFreq;
                            CorrectionErr = fabs(FractionalFreq - CorrectionFreq);
                            if (CorrectionErr > OldCorrectionErr)
                                break;
                            OldCorrectionErr = CorrectionErr;

                            /*---------------------
                             * If we found a numerator/correction-factor pair with a lower error
                             * than the previous error for this denominator, replace the previous
                             * values.
                             */
                            if (CorrectionErr < FracFreqErr) {
                                EffectiveFreq = CorrectionFreq;
                                FracFreqErr = CorrectionErr;
                                Numerator = n;
                                CorrectionIndex = k;
                            } /*end if we found a better numerator/correction pair*/

                        } /*end correction factor loop*/
                    }     /*end numerator loop*/

                    /*---------------------
                     * At this point, we now have the best numerator/correction-factor pair
                     * for the current denominator.  Recompute the Q(p) and Q(p-1) values
                     * based on the new numerator value.
                     *
                     * It could happen that the numerator/denominator/correction set we found
                     * produces the desired fractional frequency exactly (FracFreqErr == 0).
                     * If this occurs, we will set FracFreqErr to the zero-threshold value so that
                     * the algorithm gives preference to solutions that do not require a correction
                     * factor.
                     */
                    Qpm1 = Numerator;
                    Qp = j - Numerator;
                    if (FracFreqErr < ZERO_THRESHOLD)
                        FracFreqErr = ZERO_THRESHOLD;

                } /*end denominator does not exactly divide the desired fractional frequency*/

                /*---------------------
                 * Store the parameters for the numerator/correction pair that produces the lowest
                 * fractional frequency error for this denominator.  If it turns out that the
                 * denominator exactly divides the desired fractional frequency (FracFreqErr == 0),
                 * exit the denominator loop now.
                 */
                if (FracFreqErr <= 1.e-4 * DesiredFreq) {
                    BestFracFreq.Error = FracFreqErr;
                    BestFracFreq.EffectiveFreq = EffectiveFreq;
                    BestFracFreq.P = p;
                    BestFracFreq.Qp = Qp;
                    BestFracFreq.Qpm1 = Qpm1;
                    BestFracFreq.CorrectionIndex = CorrectionIndex;

                    FreqErr = (BestFracFreq.Error * ReferenceFreq) / PostDivide;
                    if (FreqErr <= 1.e-4 * DesiredFreq) {
                        Best = BestFracFreq;
                        Best.PostDivIndex = i;
                        Best.Error = FreqErr;
                        Best.EffectiveFreq = (BestFracFreq.EffectiveFreq * ReferenceFreq) / PostDivide;

                        CalcControlWord(Best, DesiredFreq, ReferenceFreq);
                        if (FreqErr < ZERO_THRESHOLD)
                            break;

                    } /*end if this post-divide solution is the best so far*/
                }     /*end if fractional frequency is best so far*/

            } /*end denominator loop*/

            /*---------------------
             * Adjust the fractional frequency error for the current post divider.
             * If the adjusted frequency error is less than the current best, make
             * this the current best.  If the new parameters produce an exact match
             * (FreqErr == 0), exit the post-divider loop now.
             */

        } /*end if VCO frequency is within range*/
    }     /*end for each post divider*/

    /*---------------------
     * Abort if we could not come up with a set of parameters that would produce the
     * desired frequency.
     */
    if (MAX_ERROR == Best.Error)
        return 0;

    /*---------------------
     * Construct the control word from the parameters found in the search above.
     */

    return 1;
} /*end FracSynthControlWord()*/

unsigned int CalcControlWord(FracSynthComponents Best, double DesiredFreq, double ReferenceFreq)
{
    unsigned int ControlWord;
    double Error;  /* Error value (in parts per million)               */
    int debugFlag; /* Flag for debug/informational output              */
    double EffectiveFreq;
    int CorrectionIndex = Best.CorrectionIndex;

    ControlWord = ((CorrectionList[CorrectionIndex].mDiv << CONTROL_MDIV_SHIFT) & CONTROL_MDIV_MASK) | ((CorrectionList[CorrectionIndex].nDiv << CONTROL_NDIV_SHIFT) & CONTROL_NDIV_MASK) | ((PostDivideList[Best.PostDivIndex].Code << CONTROL_POSTDIV_SHIFT) & CONTROL_POSTDIV_MASK) | (((Best.P - 1) << CONTROL_P_SHIFT) & CONTROL_P_MASK) | ((Best.Qpm1 << CONTROL_QPM1_SHIFT) & CONTROL_QPM1_MASK) | ((Best.Qp << CONTROL_QP_SHIFT) & CONTROL_QP_MASK);

    /*---------------------
     * Compute the frequency that will actually be generated from these parameters.
     */
    EffectiveFreq = CorrectionList[CorrectionIndex].Ratio * ReferenceFreq * ((double)Best.P - ((double)Best.Qpm1 / (double)(Best.Qp + Best.Qpm1))) / PostDivideList[Best.PostDivIndex].Divisor;

    /*---------------------
     * Return the error (in parts-per-million) between the desired and actual frequencies.
     */
    Error = 1.e6 * (EffectiveFreq - DesiredFreq) / DesiredFreq;

    /*---------------------
     * Output debug information about the results of this call
     */
    DEBUGPRINT(DP_INFO, debugFlag,
        ("Desired Frequency = %f,  Control Word = %08X\n", DesiredFreq, ControlWord));
    DEBUGPRINT(DP_INFO, debugFlag,
        ("Effective Frequency = %15.12f, Error = %5.3f ppm\n", EffectiveFreq, Error));

    /*---------------------
     * Return the computed control word.
     */

    int nDivValue = CorrectionValList[CorrectionList[CorrectionIndex].nDiv].Value;
    int mDivValue = CorrectionValList[CorrectionList[CorrectionIndex].mDiv].Value;
    if (DesiredFreq >= 100.0) {
        printf("%15.12f\t%15.12f\t%08X\t%5.3f ppm\t%d\t%d\t%d\t%d\t%d\n", DesiredFreq, EffectiveFreq, ControlWord, Error, Best.Qp, Best.Qpm1, Best.P, nDivValue, mDivValue);
    } else {
        printf("%15.12f\t\t%15.12f\t\t%08X\t%5.3f ppm\t%d\t%d\t%d\t%d\t%d\n", DesiredFreq, EffectiveFreq, ControlWord, Error, Best.Qp, Best.Qpm1, Best.P, nDivValue, mDivValue);
    }
    // printf("Desired Frequency = %f,  Control Word = %08X\n", DesiredFreq, ControlWord);
    // printf("Effective Frequency = %15.12f, Error = %5.3f ppm\n", EffectiveFreq, Error);
    // printf("Best.P = %d, Best.Qp = %d, Best.Qpm1 = %d\n", Best.P, Best.Qp, Best.Qpm1);
    return ControlWord;
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        fprintf(stderr, "\nUsage:\t%s { Desired Frequency } { Input Frequency }\n" ,argv[0]);
        exit(1);
    }
    FracSynthControlWord(atof(argv[1]), atof(argv[2]));
    return 0;
}

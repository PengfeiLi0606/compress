      real*8 function phisc( r )
c  Copyright (C) 2014, Jerry Sellwood
      implicit none
c gravitational potential in the plane of a full mass SC disc
c
c calling argument
      real*8 r
c
c externals
      real*8 aitken2
c
c local arrays
      integer ntab
      parameter ( ntab = 100 )
      real*8 x( ntab + 1 ), y( ntab + 1 )
c
c local variable
      integer i
c
      data x /  0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
     +          0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17,
     +          0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26,
     +          0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35,
     +          0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44,
     +          0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53,
     +          0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62,
     +          0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71,
     +          0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80,
     +          0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89,
     +          0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98,
     +          0.99, 1.00 /
      data ( y( i ), i = 1, 51 ) /
     +    -2.381048602689929,  -2.404907367842266,  -2.428868385633761,
     +    -2.452940196978267,  -2.477131351318171,  -2.501450418294868,
     +    -2.525905999325122,  -2.550506739118950,  -2.575261337173992,
     +    -2.600178559280778,  -2.625267249072980,  -2.650536339656575,
     +    -2.675994865351856,  -2.701651973582442,  -2.727516936945797,
     +    -2.753599165500324,  -2.779908219304827,  -2.806453821247045,
     +    -2.833245870199068,  -2.860294454538746,  -2.887609866077692,
     +    -2.915202614438221,  -2.943083441923483,  -2.971263338927265,
     +    -2.999753559932356,  -3.028565640149104,  -3.057711412848801,
     +    -3.087203027449869,  -3.117052968418496,  -3.147274075049428,
     +    -3.177879562197094,  -3.208883042032123,  -3.240298546903744,
     +    -3.272140553394450,  -3.304424007659841,  -3.337164352153713,
     +    -3.370377553846318,  -3.404080134052373,  -3.438289199994894,
     +    -3.473022478241414,  -3.508298350160644,  -3.544135889560368,
     +    -3.580554902681331,  -3.617575970737384,  -3.655220495209197,
     +    -3.693510746117772,  -3.732469913524888,  -3.772122162530804,
     +    -3.812492692065274,  -3.853607797796490,  -3.895494939514364 /
      data ( y( i ), i = 52, 101 ) /
     +    -3.938182813379963,  -3.981701429472359,  -4.026082195108252,
     +    -4.071358004459007,  -4.117563335044915,  -4.164734351748405,
     +    -4.212909019057447,  -4.262127222328589,  -4.312430898947207,
     +    -4.363864180362038,  -4.416473546083520,  -4.470307990862905,
     +    -4.525419206413675,  -4.581861779201215,  -4.639693406013972,
     +    -4.698975129243152,  -4.759771594042582,  -4.822151329820760,
     +    -4.886187058839332,  -4.951956035063237,  -5.019540416836171,
     +    -5.089027677450723,  -5.160511058257727,  -5.234090069628564,
     +    -5.309871045864749,  -5.387967761062434,  -5.468502114011121,
     +    -5.551604891467242,  -5.637416620632773,  -5.726088523433720,
     +    -5.817783587291424,  -5.912677769583120,  -6.010961355986768,
     +    -6.112840496511033,  -6.218538947365210,  -6.328300052104039,
     +    -6.442389001914672,  -6.561095422786339,  -6.684736346988433,
     +    -6.813659638259390,  -6.948247955000883,  -7.088923354397226,
     +    -7.236152663817550,  -7.390453775543647,  -7.552403058723671,
     +    -7.722644131083569,  -7.901898295885868,  -8.090977031803339,
     +    -8.290797031587220,  -9. /
c look up value in table
      phisc = aitken2( x, y, ntab + 1, r / ( 1. + r ) ) / ( 1. + r )
      return
      end

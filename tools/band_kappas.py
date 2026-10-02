#!/usr/bin/env python3
"""
Band-averaged number- (kappa_N) and energy-weighted (kappa_E) opacities for the
AREPO star_radiation waveband scheme.

Dust bands (IR/OPT/UV/LW-dust):
  Beam depletion and momentum transfer are NOT the same opacity and are emitted
  separately; see the ctr_t / catt_t / cabs_t block below for the derivation.
  Tabulated is kappa_att = C_ext/H * (1 - albedo*max(<cos>,0)), plus the
  dimensionless ratios f_mom = kappa_tr/kappa_att and f_abs = kappa_abs/kappa_att.
  All from the Draine (2003)-renormalized WD01 Milky Way R_V=3.1 model,
  kext_albedo_WD_MW_3.1_60_D03.all (astro.princeton.edu/~draine).
  Converted to per gram of GAS with M_gas/H = 2.311e-24 g
  (table header: M_dust/H = 1.398e-26 g, M_gas/M_dust = 165.3, He/H = 0.096).

Ionizing bands (HI/H2/HeI/HeII):
  Every ionizing band is averaged against ALL FOUR absorbing species, giving a
  4x4 (band, species) table.  A species only absorbs in bands lying above its
  threshold, so e.g. the HeII band (>54.4 eV) contributes to HI, HeI, HeII and
  H2, while the softest band (13.6-15.43 eV) contributes to HI only.

  The old single ionizing band 13.6-24.6 eV is SPLIT at 15.4259 eV, the H2
  ionization threshold, into IONIZING_HI_SOFT and IONIZING_H2.  Both halves are
  renamed on purpose: dropping the name "IONIZING_HI" makes every stale
  designated initializer in the C a compile error rather than a silently wrong
  number averaged over the old, wider band.  EVERY ionizing-band number in this
  file's output changes as a result, not just the new H2 column.

  H2 photoionization (H2 + gamma -> H2+ + e-) is new here.  Its cross section
  lives in sigma_h2_ion(); read the long comment there before trusting it --
  the near-threshold branch is a documented interpolation, not the published
  fit, and the script reports how much of the answer rests on it.

  The run starts with legacy_regression(), which re-derives the original
  pre-split numbers from the explicit old band edges and aborts if they no
  longer come out.  That is what licenses trusting the new numbers.

  Cross sections are the Verner, Ferland, Korista & Yakovlev (1996) analytic
  fits, NOT a sigma0*(nu/nu0)^-3 power law.  The power law is a decent
  approximation for the hydrogenic ions (HI, HeII) but is badly wrong for
  HeI, which falls off roughly as nu^-2 near threshold; nu^-3 underestimates
  HeI absorption of hard photons by a large factor.

Weights: blackbody at T_SOURCE.
  number-weighted  w_N ~ B_nu / (h nu)      (photon counts)  -> use for rates
  energy-weighted  w_E ~ B_nu               (energy flux)    -> use for kappa

Also tabulated: <E>_sigma, the sigma-and-photon-weighted mean photon energy
per (band, species).  <E>_sigma - E_th is the mean excess energy carried into
photoheating per photoionization of that species by photons of that band.

Replace the blackbody with per-age stellar SEDs later; only OPT is sensitive
among the dust bands, but ALL the ionizing ratios move with the SED shape.
"""
import numpy as np

T_SOURCE = 4.0e4          # K
M_GAS_PER_H = 2.311e-24   # g

h, c, kB = 6.626e-27, 2.998e10, 1.381e-16
EV_UM = 1.23984           # lambda[um] = EV_UM / E[eV]
EV_ERG = 1.602177e-12     # erg per eV

trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz  # renamed in numpy 2.0

# --- D03 table subsample: lambda[um], albedo, g, C_ext/H [cm^2/H] -------------
D03 = np.array([
[0.0912011,0.2386,0.6619,2.406e-21],[0.0933254,0.2457,0.6562,2.356e-21],
[0.0954993,0.2536,0.6516,2.292e-21],[0.0977237,0.2658,0.6496,2.181e-21],
[0.100000,0.2761,0.6501,2.083e-21],[0.102329,0.2850,0.6494,2.010e-21],
[0.104713,0.2963,0.6503,1.920e-21],[0.107152,0.3091,0.6532,1.818e-21],
[0.109648,0.3216,0.6584,1.715e-21],[0.110000,0.3232,0.6591,1.703e-21],
[0.112202,0.3302,0.6636,1.645e-21],[0.114815,0.3357,0.6698,1.593e-21],
[0.117490,0.3424,0.6730,1.545e-21],[0.120226,0.3494,0.6755,1.500e-21],
[0.123027,0.3565,0.6774,1.458e-21],[0.125893,0.3632,0.6789,1.420e-21],
[0.128825,0.3692,0.6801,1.385e-21],[0.131826,0.3748,0.6808,1.350e-21],
[0.134896,0.3789,0.6807,1.321e-21],[0.138038,0.3825,0.6795,1.294e-21],
[0.141254,0.3863,0.6766,1.268e-21],[0.144544,0.3916,0.6719,1.244e-21],
[0.147911,0.3986,0.6673,1.218e-21],[0.151356,0.4068,0.6633,1.193e-21],
[0.154882,0.4162,0.6593,1.171e-21],[0.158489,0.4262,0.6550,1.153e-21],
[0.162181,0.4377,0.6497,1.139e-21],[0.165959,0.4505,0.6433,1.129e-21],
[0.169824,0.4643,0.6358,1.122e-21],[0.173780,0.4784,0.6278,1.119e-21],
[0.177828,0.4924,0.6190,1.120e-21],[0.181970,0.5050,0.6100,1.124e-21],
[0.186209,0.5150,0.6014,1.133e-21],[0.190546,0.5211,0.5932,1.149e-21],
[0.194984,0.5234,0.5850,1.172e-21],[0.199526,0.5212,0.5774,1.203e-21],
[0.204174,0.5155,0.5705,1.240e-21],[0.208930,0.5082,0.5644,1.275e-21],
[0.213796,0.5035,0.5586,1.302e-21],[0.217500,0.5046,0.5544,1.308e-21],
[0.220000,0.5077,0.5520,1.302e-21],[0.223872,0.5151,0.5492,1.281e-21],
[0.229087,0.5279,0.5463,1.242e-21],[0.234423,0.5426,0.5439,1.196e-21],
[0.239883,0.5567,0.5424,1.151e-21],[0.245471,0.5675,0.5424,1.108e-21],
[0.251189,0.5762,0.5428,1.072e-21],[0.257039,0.5834,0.5436,1.040e-21],
[0.263027,0.5888,0.5447,1.013e-21],[0.269153,0.5933,0.5462,9.884e-22],
[0.275423,0.5973,0.5480,9.662e-22],[0.281838,0.6008,0.5503,9.452e-22],
[0.288403,0.6042,0.5528,9.257e-22],[0.295121,0.6075,0.5550,9.078e-22],
[0.301995,0.6108,0.5571,8.909e-22],[0.309029,0.6150,0.5591,8.732e-22],
[0.316228,0.6195,0.5610,8.555e-22],[0.323594,0.6239,0.5629,8.383e-22],
[0.331131,0.6283,0.5647,8.215e-22],[0.338844,0.6324,0.5662,8.051e-22],
[0.346737,0.6362,0.5676,7.892e-22],[0.354814,0.6396,0.5687,7.736e-22],
[0.363078,0.6428,0.5695,7.583e-22],[0.371535,0.6459,0.5700,7.431e-22],
[0.380189,0.6489,0.5703,7.279e-22],[0.389045,0.6518,0.5702,7.127e-22],
[0.398107,0.6546,0.5699,6.973e-22],[0.407380,0.6573,0.5693,6.818e-22],
[0.416869,0.6601,0.5685,6.662e-22],[0.426580,0.6627,0.5674,6.505e-22],
[0.436516,0.6653,0.5661,6.348e-22],[0.446684,0.6676,0.5645,6.192e-22],
[0.457088,0.6697,0.5626,6.036e-22],[0.467735,0.6716,0.5604,5.882e-22],
[0.478630,0.6731,0.5579,5.729e-22],[0.489779,0.6744,0.5552,5.577e-22],
[0.501187,0.6754,0.5522,5.427e-22],[0.512861,0.6762,0.5489,5.277e-22],
[0.524807,0.6768,0.5454,5.129e-22],[0.537032,0.6772,0.5415,4.983e-22],
[0.549541,0.6774,0.5374,4.839e-22],[0.562341,0.6774,0.5330,4.697e-22],
[0.575440,0.6772,0.5283,4.558e-22],[0.588844,0.6769,0.5233,4.421e-22],
[0.602560,0.6763,0.5181,4.287e-22],[0.616595,0.6755,0.5126,4.155e-22],
[0.630958,0.6745,0.5069,4.025e-22],[0.645654,0.6734,0.5010,3.898e-22],
[0.660694,0.6721,0.4949,3.773e-22],[0.676083,0.6707,0.4886,3.651e-22],
[0.691831,0.6691,0.4821,3.532e-22],[0.707946,0.6672,0.4754,3.416e-22],
[0.724436,0.6652,0.4687,3.303e-22],[0.741310,0.6629,0.4619,3.192e-22],
[0.758578,0.6606,0.4550,3.085e-22],[0.776247,0.6580,0.4481,2.980e-22],
[0.794328,0.6553,0.4411,2.878e-22],[0.812830,0.6525,0.4340,2.780e-22],
[0.831764,0.6495,0.4268,2.683e-22],[0.851138,0.6463,0.4194,2.590e-22],
[0.870964,0.6429,0.4119,2.499e-22],[0.891251,0.6393,0.4043,2.411e-22],
[0.912011,0.6355,0.3964,2.326e-22],[0.933254,0.6314,0.3883,2.243e-22],
[0.954993,0.6270,0.3801,2.164e-22],[0.977237,0.6219,0.3716,2.089e-22],
[1.00000,0.6153,0.3630,2.021e-22],[1.02329,0.6047,0.3542,1.967e-22],
[1.04713,0.5903,0.3454,1.928e-22],[1.07152,0.5927,0.3366,1.837e-22],
[1.09648,0.5960,0.3278,1.748e-22],[1.12202,0.5945,0.3192,1.676e-22],
[1.14815,0.5909,0.3107,1.613e-22],[1.17490,0.5863,0.3024,1.555e-22],
[1.20226,0.5809,0.2943,1.501e-22],[1.23027,0.5748,0.2865,1.451e-22],
[1.25893,0.5690,0.2790,1.402e-22],[1.31826,0.5622,0.2647,1.298e-22],
[1.41254,0.5503,0.2451,1.157e-22],[1.51356,0.5354,0.2271,1.035e-22],
[1.62181,0.5192,0.2097,9.253e-23],[1.73780,0.5022,0.1923,8.258e-23],
[1.86209,0.4842,0.1745,7.355e-23],[1.99526,0.4648,0.1561,6.543e-23],
[2.13796,0.4445,0.1370,5.809e-23],[2.29087,0.4232,0.1170,5.143e-23],
[2.45471,0.4009,0.0968,4.543e-23],[2.63027,0.3779,0.0763,4.003e-23],
[2.81838,0.3539,0.0560,3.519e-23],[3.01995,0.3290,0.0367,3.088e-23],
[3.23594,0.3032,0.0189,2.709e-23],[3.46737,0.2773,0.0031,2.372e-23],
[3.71535,0.2511,-0.0105,2.076e-23],[3.98107,0.2249,-0.0219,1.818e-23],
[4.26580,0.1993,-0.0313,1.593e-23],[4.57088,0.1744,-0.0388,1.400e-23],
[4.89779,0.1505,-0.0445,1.238e-23],[5.24808,0.1273,-0.0489,1.107e-23],
[5.62341,0.1060,-0.0519,1.000e-23],[6.16595,0.0757,-0.0538,9.503e-24],
[6.60693,0.0644,-0.0540,8.319e-24],[7.07946,0.0511,-0.0535,7.771e-24],
[7.58578,0.0318,-0.0525,9.209e-24],[8.12831,0.0135,-0.0508,1.610e-23],
[8.70964,0.0062,-0.0470,2.701e-23],[9.50000,0.0032,-0.0385,4.129e-23],
[10.0000,0.0030,-0.0354,3.736e-23],[10.7152,0.0029,-0.0324,2.956e-23],
[11.7490,0.0028,-0.0301,2.017e-23],[12.3027,0.0028,-0.0294,1.653e-23]])

# --- D03 EUV/soft-X-ray extension: 13.6 - 206 eV -----------------------------
# Same file, rows short-ward of the Lyman edge.  NOTE: the published table has
# two corrupt rows in the X-ray section (commented "156 eV" / "157 eV" but
# carrying lambda = 7.94770E-02 / 7.89710E-02 um, i.e. a wrong exponent and the
# 15.6/15.7 eV data).  They are excluded here; if you re-parse the file
# yourself, drop them or sorting by lambda will silently duplicate FUV points.
# NOTE on <cos>: g -> 1 as lambda -> 0 in this model.  That is physical, not a
# tabulation artefact: EUV/X-ray scattering off grains is strongly forward
# (Rayleigh-Gans, Draine 2003c).  All seven bands therefore use the same
# transport opacity C_ext*(1-a<g>), which tends to C_ext*(1-a) = C_abs in this
# limit -- 8%/7%/3% above C_abs in the HI/HeI/HeII bands.  These g values are
# NOT usable as a Henyey-Greenstein phase function for direction sampling.
D03_EUV = np.array([
[0.00603,0.5837,0.9975,5.544e-22],[0.00617,0.5793,0.9974,5.627e-22],
[0.00631,0.5749,0.9972,5.709e-22],[0.00646,0.5707,0.9970,5.791e-22],
[0.00661,0.5665,0.9968,5.871e-22],[0.00676,0.5620,0.9966,5.953e-22],
[0.00692,0.5571,0.9964,6.037e-22],[0.00708,0.5501,0.9962,6.121e-22],
[0.00724,0.5456,0.9960,6.180e-22],[0.00741,0.5437,0.9957,6.256e-22],
[0.00759,0.5358,0.9954,6.336e-22],[0.00776,0.5339,0.9952,6.397e-22],
[0.00794,0.5213,0.9948,6.516e-22],[0.00796,0.5189,0.9948,6.513e-22],
[0.00797,0.5174,0.9948,6.505e-22],[0.00798,0.5164,0.9949,6.491e-22],
[0.00813,0.5265,0.9948,6.467e-22],[0.00832,0.5251,0.9944,6.557e-22],
[0.00851,0.5225,0.9941,6.633e-22],[0.00871,0.5192,0.9937,6.711e-22],
[0.00891,0.5164,0.9933,6.790e-22],[0.00912,0.5140,0.9929,6.874e-22],
[0.00933,0.5103,0.9924,6.976e-22],[0.00955,0.5052,0.9918,7.093e-22],
[0.00977,0.4961,0.9910,7.242e-22],[0.01000,0.4814,0.9908,7.283e-22],
[0.01023,0.4781,0.9906,7.259e-22],[0.01047,0.4831,0.9901,7.295e-22],
[0.01072,0.4765,0.9891,7.503e-22],[0.01096,0.4597,0.9887,7.564e-22],
[0.01122,0.4636,0.9885,7.512e-22],[0.01148,0.4640,0.9872,7.726e-22],
[0.01170,0.4323,0.9853,8.206e-22],[0.01172,0.4267,0.9855,8.202e-22],
[0.01174,0.4213,0.9857,8.181e-22],[0.01175,0.4199,0.9860,8.138e-22],
[0.01176,0.4177,0.9864,8.062e-22],[0.01181,0.4128,0.9875,7.808e-22],
[0.01192,0.4199,0.9894,7.153e-22],[0.01202,0.4340,0.9894,7.077e-22],
[0.01230,0.4496,0.9884,7.215e-22],[0.01240,0.4523,0.9881,7.238e-22],
[0.01259,0.4553,0.9874,7.328e-22],[0.01288,0.4561,0.9866,7.415e-22],
[0.01318,0.4533,0.9855,7.572e-22],[0.01349,0.4532,0.9851,7.539e-22],
[0.01380,0.4551,0.9840,7.645e-22],[0.01413,0.4537,0.9829,7.747e-22],
[0.01445,0.4513,0.9818,7.847e-22],[0.01479,0.4486,0.9807,7.936e-22],
[0.01514,0.4459,0.9796,8.027e-22],[0.01549,0.4428,0.9784,8.117e-22],
[0.01585,0.4397,0.9772,8.203e-22],[0.01622,0.4367,0.9759,8.289e-22],
[0.01660,0.4341,0.9744,8.389e-22],[0.01698,0.4284,0.9732,8.479e-22],
[0.01738,0.4262,0.9718,8.550e-22],[0.01778,0.4217,0.9705,8.628e-22],
[0.01820,0.4216,0.9691,8.675e-22],[0.01862,0.4197,0.9670,8.815e-22],
[0.01905,0.4153,0.9649,8.962e-22],[0.01950,0.4093,0.9631,9.095e-22],
[0.01995,0.4015,0.9618,9.176e-22],[0.02042,0.3995,0.9608,9.166e-22],
[0.02089,0.3992,0.9588,9.248e-22],[0.02138,0.3974,0.9564,9.382e-22],
[0.02188,0.3967,0.9530,9.584e-22],[0.02239,0.3854,0.9487,9.982e-22],
[0.02279,0.3689,0.9489,1.010e-21],[0.02291,0.3644,0.9502,9.996e-22],
[0.02344,0.3666,0.9529,9.471e-22],[0.02399,0.3721,0.9507,9.496e-22],
[0.02455,0.3741,0.9476,9.640e-22],[0.02512,0.3734,0.9446,9.791e-22],
[0.02570,0.3720,0.9416,9.942e-22],[0.02630,0.3695,0.9385,1.010e-21],
[0.02692,0.3663,0.9355,1.026e-21],[0.02754,0.3632,0.9327,1.040e-21],
[0.02818,0.3608,0.9294,1.056e-21],[0.02884,0.3574,0.9261,1.074e-21],
[0.02951,0.3537,0.9227,1.092e-21],[0.03020,0.3498,0.9193,1.110e-21],
[0.03090,0.3458,0.9159,1.128e-21],[0.03162,0.3416,0.9125,1.147e-21],
[0.03236,0.3373,0.9092,1.167e-21],[0.03311,0.3329,0.9057,1.186e-21],
[0.03388,0.3269,0.9031,1.206e-21],[0.03467,0.3241,0.9010,1.213e-21],
[0.03548,0.3285,0.8981,1.208e-21],[0.03631,0.3282,0.8928,1.228e-21],
[0.03715,0.3260,0.8878,1.249e-21],[0.03802,0.3237,0.8829,1.269e-21],
[0.03890,0.3215,0.8777,1.291e-21],[0.03981,0.3191,0.8720,1.315e-21],
[0.04074,0.3167,0.8662,1.339e-21],[0.04169,0.3148,0.8594,1.367e-21],
[0.04266,0.3123,0.8513,1.401e-21],[0.04365,0.3070,0.8442,1.441e-21],
[0.04467,0.3016,0.8377,1.479e-21],[0.04571,0.2968,0.8313,1.514e-21],
[0.04677,0.2914,0.8250,1.552e-21],[0.04786,0.2860,0.8189,1.592e-21],
[0.04881,0.2817,0.8138,1.624e-21],[0.04898,0.2810,0.8129,1.630e-21],
[0.05012,0.2762,0.8070,1.668e-21],[0.05129,0.2716,0.8009,1.708e-21],
[0.05248,0.2672,0.7946,1.749e-21],[0.05370,0.2634,0.7880,1.791e-21],
[0.05495,0.2601,0.7810,1.834e-21],[0.05623,0.2579,0.7735,1.875e-21],
[0.05754,0.2580,0.7655,1.906e-21],[0.05888,0.2478,0.7637,1.999e-21],
[0.06026,0.2416,0.7543,2.094e-21],[0.06166,0.2351,0.7445,2.199e-21],
[0.06310,0.2281,0.7342,2.320e-21],[0.06457,0.2206,0.7238,2.454e-21],
[0.06607,0.2127,0.7132,2.606e-21],[0.06761,0.2047,0.7033,2.766e-21],
[0.06918,0.1965,0.6945,2.935e-21],[0.07079,0.1908,0.6880,3.060e-21],
[0.07244,0.1902,0.6829,3.096e-21],[0.07413,0.1956,0.6791,3.028e-21],
[0.07586,0.2041,0.6764,2.908e-21],[0.07762,0.2135,0.6750,2.773e-21],
[0.07943,0.2211,0.6759,2.653e-21],[0.08128,0.2272,0.6765,2.558e-21],
[0.08318,0.2313,0.6769,2.488e-21],[0.08511,0.2330,0.6763,2.447e-21],
[0.08710,0.2329,0.6751,2.427e-21],[0.08913,0.2329,0.6715,2.420e-21],
[0.09120,0.2386,0.6619,2.406e-21]])


# --- D03 far-IR extension: 12.6 um - 3.0 mm ----------------------------------
# Needed only for the IR trapping opacity: a 100 K Planck peaks at 29 um and a
# 20 K one at 145 um, so the original 12.3 um cutoff sampled the Wien tail only.
# Columns are [lambda_um, C_ext/H]; albedo and <cos> are taken as zero because
# the tabulated albedo stays below 0.003 across this whole range, so absorption,
# extinction and transport opacity agree to better than a third of a percent.
D03_FIR = np.array([
[12.589,1.505e-23],[12.883,1.362e-23],[13.183,1.246e-23],
[13.49,1.146e-23],[13.804,1.067e-23],[14.125,1.037e-23],
[14.454,1.052e-23],[14.791,1.098e-23],[15.136,1.153e-23],
[15.488,1.216e-23],[15.849,1.286e-23],[16.218,1.359e-23],
[16.596,1.434e-23],[16.982,1.500e-23],[17.378,1.556e-23],
[17.783,1.596e-23],[18.197,1.614e-23],[18.621,1.611e-23],
[19.055,1.585e-23],[19.498,1.538e-23],[19.953,1.476e-23],
[20.417,1.411e-23],[20.893,1.358e-23],[21.38,1.304e-23],
[21.878,1.257e-23],[22.387,1.207e-23],[22.909,1.158e-23],
[23.442,1.117e-23],[23.988,1.074e-23],[24.547,1.035e-23],
[25.119,9.952e-24],[25.704,9.562e-24],[26.303,9.185e-24],
[26.915,8.820e-24],[27.542,8.467e-24],[28.184,8.124e-24],
[28.84,7.792e-24],[29.512,7.470e-24],[30.2,7.156e-24],
[31.623,6.555e-24],[33.113,5.999e-24],[34.674,5.493e-24],
[36.308,5.028e-24],[38.019,4.596e-24],[39.811,4.186e-24],
[41.687,3.806e-24],[43.652,3.457e-24],[45.709,3.138e-24],
[47.863,2.847e-24],[50.119,2.582e-24],[52.481,2.340e-24],
[54.954,2.119e-24],[57.544,1.917e-24],[60.256,1.734e-24],
[63.096,1.568e-24],[66.069,1.417e-24],[69.183,1.280e-24],
[72.444,1.157e-24],[75.858,1.045e-24],[79.433,9.444e-25],
[83.176,8.534e-25],[87.096,7.715e-25],[91.201,6.981e-25],
[95.499,6.320e-25],[100,5.726e-25],[104.71,5.192e-25],
[109.65,4.714e-25],[114.81,4.288e-25],[120.23,3.911e-25],
[125.89,3.585e-25],[131.83,3.304e-25],[138.04,3.033e-25],
[144.54,2.731e-25],[151.36,2.455e-25],[158.49,2.216e-25],
[165.96,2.007e-25],[173.78,1.822e-25],[181.97,1.655e-25],
[190.55,1.505e-25],[199.53,1.369e-25],[208.93,1.246e-25],
[218.78,1.134e-25],[229.09,1.033e-25],[239.88,9.407e-26],
[251.19,8.533e-26],[263.03,7.690e-26],[275.42,6.935e-26],
[288.4,6.264e-26],[302,5.665e-26],[316.23,5.130e-26],
[331.13,4.651e-26],[346.74,4.223e-26],[363.08,3.838e-26],
[380.19,3.493e-26],[398.11,3.184e-26],[416.87,2.905e-26],
[436.52,2.654e-26],[457.09,2.428e-26],[478.63,2.223e-26],
[501.19,2.038e-26],[524.81,1.871e-26],[549.54,1.719e-26],
[575.44,1.582e-26],[602.56,1.457e-26],[630.96,1.343e-26],
[660.69,1.240e-26],[691.83,1.145e-26],[724.44,1.059e-26],
[758.58,9.801e-27],[794.33,9.074e-27],[831.76,8.402e-27],
[870.96,7.780e-27],[912.01,7.203e-27],[954.99,6.669e-27],
[1000,6.174e-27],[1047.1,5.717e-27],[1096.5,5.295e-27],
[1148.2,4.904e-27],[1202.3,4.543e-27],[1258.9,4.206e-27],
[1318.3,3.896e-27],[1380.4,3.608e-27],[1445.4,3.342e-27],
[1513.6,3.095e-27],[1584.9,2.866e-27],[1659.6,2.654e-27],
[1737.8,2.458e-27],[1819.7,2.276e-27],[1905.5,2.108e-27],
[1995.3,1.952e-27],[2089.3,1.808e-27],[2187.8,1.674e-27],
[2290.9,1.550e-27],[2398.8,1.435e-27],[2511.9,1.329e-27],
[2630.3,1.230e-27],[2754.2,1.138e-27],[2884,1.053e-27],
[3019.9,9.750e-28]])

_fir4 = np.column_stack([D03_FIR[:, 0], np.zeros(len(D03_FIR)),
                         np.zeros(len(D03_FIR)), D03_FIR[:, 1]])
D03_ALL = np.vstack([D03_EUV, D03, _fir4])
lam_t, alb_t, g_t, cext_t = D03_ALL.T
# Three distinct opacities.  They are NOT interchangeable.
#
#   ctr_t   transport / momentum-transfer, kappa_abs + (1-g) kappa_sca.  This is
#           the correct coefficient for the radiation FORCE at any g in [-1,1],
#           because momentum transfer per scattering goes as (1 - <cos t>).
#
#   catt_t  beam depletion.  The transport (delta-scaling) approximation writes
#           Phi(mu) ~ 2g*delta(1-mu) + (1-g)*Phi_iso and drops the forward delta
#           as indistinguishable from "did not scatter".  That decomposition has
#           non-negative weights ONLY for g >= 0.  For g < 0 there is no forward
#           delta to drop and the beam loses exactly kappa_ext: removal from the
#           beam is a zeroth-moment statement and is blind to direction.  Hence
#           the clamp max(g,0), equivalently catt = min(ctr, cext).  Using ctr
#           here would delete photons that do not exist.
#
#   cabs_t  true absorption; the only part that heats grains.
#
# Ordering: cabs <= catt <= cext always, and ctr > cext exactly where g < 0.
# For D03 MW R_V=3.1 that is 3.715-12.3 um only, where the albedo has already
# collapsed to <=0.25, so ctr/cext peaks at just 1.00677 (4.57 um).  Band-
# averaged the floor is therefore worth -0.004% in INFRARED and nothing at all
# in the other seven bands, i.e. f_mom rounds to 1.0000 everywhere and the
# emitted tables do not move.  It is in here for correctness of definition, so
# that the three quantities stay separable when the weighting changes -- a cool
# reprocessed-IR spectrum, or a grain model with a deeper backscattering dip,
# puts real weight where g < 0.  gfloor_check() prints the current margin.
ctr_t  = cext_t * (1.0 - alb_t * g_t)
catt_t = cext_t * (1.0 - alb_t * np.maximum(g_t, 0.0))
cabs_t = cext_t * (1.0 - alb_t)

def planck_weight(l_um, T, kind):
    """kind='N': photon number spectrum ~ B_lambda*lambda; 'E': ~ B_lambda"""
    l = l_um * 1e-4
    x = np.minimum(h * c / (l * kB * T), 500.0)
    B = 1.0 / (l**5 * np.expm1(x))
    return B * l if kind == 'N' else B

def band_avg(y_tab, l0, l1, T, kind, n=4000):
    lg = np.linspace(l0, l1, n)
    y = np.interp(lg, lam_t, y_tab)
    w = planck_weight(lg, T, kind)
    return trapz(w * y, lg) / trapz(w, lg)

# --- Dust bands: [name, lambda_lo, lambda_hi] in um ---------------------------
DUST_BANDS = [("INFRARED",    EV_UM/1.0, 12.3984),   # 0-1 eV, integrated to 12.4um only!
              ("OPTICAL",     EV_UM/6.0, EV_UM/1.0),
              ("ULTRAVIOLET", EV_UM/11.2, EV_UM/6.0),
              ("LYMAN_WERNER",EV_UM/13.6, EV_UM/11.2)]

# --- Ionizing bands: [name, E_lo eV, E_hi eV] --------------------------------
# Band edges sit exactly on the thresholds so the (band, species) table comes
# out lower-triangular.  The 200 eV top cutoff is arbitrary: a 40 kK Planck has
# nothing there (exp(-58)), but it matters if you swap in an SED with X-rays.
#
# The old single IONIZING_HI band (13.60-24.59 eV) is split at the H2
# ionization threshold, 15.4259 eV, so that H2 gets a band whose lower edge IS
# its threshold.  Without the split, a band-averaged sigma_H2 over 13.6-24.59
# would be diluted by 1.8 eV of sub-threshold photons AND the photon- vs
# energy-weighted means would no longer bracket the threshold cleanly.
#
# IONIZING_HI KEEPS ITS NAME and now means 13.60-15.4259 eV only.  Because the
# name is unchanged, the compiler will NOT flag the five designated-initializer
# tables in star_radiation.c / star_radiation_statistics.c that still hold
# values averaged over the old, wider 13.60-24.59 eV band.  Every one of them
# must be repasted from this script's output by hand -- see the checklist in
# the plan.  The numbers below are the whole point: sigma_HI in the narrowed
# band is 5.40e-18, not 3.67e-18, a 47% change.
E_TH_H2 = 15.4259            # eV, H2 X(1Sigma_g+) v=0 -> H2+ X(2Sigma_g+) v=0

ION_BANDS = [("IONIZING_HI",   13.60,   E_TH_H2),
             ("IONIZING_H2",   E_TH_H2, 24.59),
             ("IONIZING_HeI",  24.59,   54.42),
             ("IONIZING_HeII", 54.42,   200.0)]

# --- Verner, Ferland, Korista & Yakovlev (1996) ApJ 465, 487 -----------------
#   x = E/E0 - y0,   y = sqrt(x^2 + y1^2)
#   F(y) = ((x-1)^2 + yw^2) * y^(0.5P - 5.5) * (1 + sqrt(y/ya))^(-P)
#   sigma(E) = sigma0 * F(y),  sigma0 in Mb = 1e-18 cm^2,  zero below E_th.
# Threshold values reproduced: HI 6.30e-18, HeI 7.42e-18, HeII 1.58e-18 cm^2.
# (The old 7.83e-18 for HeI is the Osterbrock 1989 value; Verner supersedes it.)
MB = 1.0e-18
# Order matters and is load-bearing: it is the order of the SP_* enum emitted
# for the C, and the C relies on CH_HI + s == the channel for species s.  So
# SPECIES, the CH_HI..CH_HeII block of the channel enum, and the SP_* enum must
# all stay in this same order -- ascending in threshold energy:
#   HI 13.60   H2 15.4259   HeI 24.59   HeII 54.42  eV
SPECIES = ["HI", "H2", "HeI", "HeII"]
VERNER = {
    #        E_th    E_max     E0        sigma0[Mb]  ya       P      yw     y0       y1
    "HI":   (13.60,  5.0e4,  4.298e-1,   5.475e4,  3.288e1, 2.963, 0.000, 0.0000,  0.000),
    "HeI":  (24.59,  5.0e3,  1.361e1,    9.492e2,  1.469,   3.188, 2.039, 4.434e-1, 2.136),
    "HeII": (54.42,  5.0e4,  1.720,      1.369e4,  3.288e1, 2.963, 0.000, 0.0000,  0.000),
}
E_TH = {s: VERNER[s][0] for s in VERNER}
E_TH["H2"] = E_TH_H2

def sigma_verner(sp, E_eV):
    """Photoionization cross section [cm^2]; 0 below threshold."""
    Eth, Emax, E0, s0, ya, P, yw, y0, y1 = VERNER[sp]
    E = np.atleast_1d(np.asarray(E_eV, dtype=float))
    out = np.zeros_like(E)
    m = (E >= Eth) & (E <= Emax)
    if np.any(m):
        x = E[m] / E0 - y0
        y = np.sqrt(x * x + y1 * y1)
        F = ((x - 1.0) ** 2 + yw * yw) * y ** (0.5 * P - 5.5) \
            * (1.0 + np.sqrt(y / ya)) ** (-P)
        out[m] = s0 * MB * F
    return out if np.ndim(E_eV) else out[0]


# --- H2 photoionization: H2 + gamma -> H2+ + e-, threshold 15.4259 eV --------
#
# ### READ THIS BEFORE TRUSTING THE NUMBERS ###
#
# H2 is NOT a Verner species, so it needs its own cross section.  The standard
# astrophysical reference is Yan, Sadeghpour & Dalgarno (1998) ApJ 496, 1044,
# *with* the erratum (2001) ApJ 559, 1194 -- the erratum exists precisely
# because the H2 analytic representation (their eqs. 17-19) was misprinted in
# the original paper.  Neither the paper nor the erratum was reachable from the
# machine this was written on, so the two branches below have different
# confidence levels and are labelled accordingly.  Replace HI_BRANCH when you
# have the erratum in hand; the function is the ONLY place H2 sigma enters.
#
# Branch A, 18-85 eV -- HIGH CONFIDENCE, independently cross-checked.
#   sigma = 2e-17 * (c0 x^-s + c1 x^(-s-1) + c2 x^(-s-2) + c3 x^(-s-3))
#   x = E/15.4,  s = 0.252
# Evaluated at 18 eV this gives 9.79 Mb against the recommended experimental
# value of 9.85 Mb (Samson & Haddad 1994, total photoabsorption, +-2-3%),
# i.e. agreement to 0.6%.  That independent anchor is why this branch is
# trusted.  Above 18 eV photoabsorption is essentially all ionization.
#
# Branch B, 15.4259-18 eV -- LOW CONFIDENCE, documented interpolation.
# A cubic in x circulates for this range but it is unusable as published at
# five significant figures: at 18 eV its four terms are -37.9 + 116.6 - 119.2
# + 40.6, a 0.05% residual from ~100-sized terms, so rounding alone moves the
# answer by >100%, and it evaluates to 5.9 Mb where the measured value is
# 9.85 Mb.  Catastrophic cancellation plus a 40% miss against experiment means
# it is either the misprinted form or truncated beyond use.  Rather than
# propagate that, this branch is a monotone interpolation pinned to two
# sourced endpoints:
#     sigma(15.4259 eV) = 1.0 Mb   (theory, "~1e-18 cm^2 at threshold")
#     sigma(18.0    eV) = 9.85 Mb  (experiment, recommended)
# linear in E between them.  The true shape carries strong autoionization
# resonances, so this is a smooth stand-in for a resonance-averaged rate --
# defensible for a band average, wrong in detail.
#
# THIS BRANCH DOMINATES THE UNCERTAINTY for a 40 kK source: hnu/kT = 4.5 at
# 15.4 eV, so the photon weight falls steeply across the band and the
# 15.4-18 eV sub-range carries much of it.  The script prints that weight
# fraction explicitly (see report_h2_uncertainty) -- read it before quoting
# sigma_H2 to better than a factor ~1.5.
#
# SCOPE NOTE: above 18.08 eV some ionization goes to H + H+ + e (dissociative
# ionization) rather than H2+ + e.  This cross section is the TOTAL, so using
# it for grackle's k29 (H2 + gamma -> H2+ + e) attributes every ionization to
# the H2+ channel and overestimates H2+ production by the dissociative
# branching ratio -- a few per cent near 20 eV, growing with energy.  That is
# a deliberate, agreed simplification; see the plan.

H2_X0 = 15.4                       # the fit's reference energy, NOT E_TH_H2
H2_S = 0.252
H2_C = (0.071, -0.673, 1.977, -0.692)
H2_BRANCH_A_LO = 18.0
H2_BRANCH_A_HI = 85.0

H2_ANCHORS = ((E_TH_H2, 1.00e-18),   # theory, low confidence
              (18.0,    9.85e-18))   # experiment, recommended

def _sigma_h2_branchA(E):
    """YSD98 18-85 eV form.  Extrapolated above 85 eV -- see sigma_h2_ion."""
    x = E / H2_X0
    xs = x ** (-H2_S)
    return 2.0e-17 * (H2_C[0] * xs + H2_C[1] * xs / x
                      + H2_C[2] * xs / x**2 + H2_C[3] * xs / x**3)

def sigma_h2_ion(E_eV):
    """H2 photoionization cross section [cm^2]; 0 below 15.4259 eV.

    Piecewise: interpolated near-threshold (B), YSD98 above 18 eV (A).
    Above 85 eV branch A is extrapolated beyond its stated validity; for a
    40 kK Planck that region holds <0.02% of the IONIZING_HeII band photons,
    but swap in an X-ray SED and it must be replaced with the YSD98 >85 eV
    branch.
    """
    E = np.atleast_1d(np.asarray(E_eV, dtype=float))
    out = np.zeros_like(E)

    (e_lo, s_lo), (e_hi, s_hi) = H2_ANCHORS
    mB = (E >= e_lo) & (E < e_hi)
    if np.any(mB):
        out[mB] = s_lo + (s_hi - s_lo) * (E[mB] - e_lo) / (e_hi - e_lo)

    mA = E >= H2_BRANCH_A_LO
    if np.any(mA):
        out[mA] = np.maximum(_sigma_h2_branchA(E[mA]), 0.0)

    return out if np.ndim(E_eV) else out[0]

def sigma_species(sp, E_eV):
    """Dispatch to the right cross section for species `sp`."""
    return sigma_h2_ion(E_eV) if sp == "H2" else sigma_verner(sp, E_eV)

def planck_E(E_eV, T, kind):
    """kind='N': photons  ~ E^2/(e^x-1);  'E': energy ~ E^3/(e^x-1)  (per dE)"""
    x = np.minimum(E_eV * EV_ERG / (kB * T), 700.0)
    p = 3.0 if kind == 'E' else 2.0
    return E_eV ** p / np.expm1(x)

def _egrid(e0, e1, sp=None, n=6000):
    """Log grid over the FULL band; kinks inserted so they are resolved.

    Breakpoints are the species threshold and, for H2, the 18 eV boundary
    between the interpolated and YSD98 branches.  Without the 18 eV point the
    trapezoid rule smears a ~10% slope discontinuity across one coarse cell.
    """
    g = np.logspace(np.log10(e0), np.log10(e1), n)
    if sp is not None:
        kinks = [E_TH[sp]]
        if sp == "H2":
            kinks += [H2_BRANCH_A_LO]
        for ek in kinks:
            if e0 < ek < e1:
                g = np.concatenate([g, [ek * (1 - 1e-12), ek, ek * (1 + 1e-12)]])
        g = np.sort(g)
    return g

def band_sigma(sp, e0, e1, T, kind):
    """<sigma_sp> over [e0,e1].  Denominator is the FULL band weight: photons
    below threshold are in the band but simply do not ionize this species."""
    if e1 <= E_TH[sp] * (1.0 + 1e-9):
        return 0.0          # band entirely below threshold (edges may coincide)
    E = _egrid(e0, e1, sp)
    w = planck_E(E, T, kind)
    return trapz(w * sigma_species(sp, E), E) / trapz(w, E)

def dust_avg_E(y_tab, e0, e1, T, kind, n=4000):
    """Band-average a dust quantity over [e0,e1] in eV, interpolating the D03
    table in log-energy.  Same weights as band_sigma so the gas and dust
    opacities in a given band are averaged consistently."""
    E = np.logspace(np.log10(e0), np.log10(e1), n)
    y = np.interp(EV_UM / E, lam_t, y_tab)
    w = planck_E(E, T, kind)
    return trapz(w * y, E) / trapz(w, E)

def kappa_rosseland(T_d, n=6000):
    """Rosseland mean dust opacity [cm^2/g gas] at dust temperature T_d.
    1/kappa_R = <1/kappa>_{dB/dT}.  This is the flux-mean in the diffusion
    limit, so it is the correct opacity for the tau_IR trapping term in the
    (1 + tau_IR) momentum boost.  Weighted by T_d, NOT by the source.
    Averages the TRANSPORT opacity: a flux mean is a first-moment quantity, so
    ctr_t is right here even where g < 0."""
    lam = np.logspace(np.log10(lam_t.min()), np.log10(lam_t.max()), n)
    kap = np.interp(lam, lam_t, ctr_t) / M_GAS_PER_H
    l = lam * 1e-4
    x = np.minimum(h * c / (l * kB * T_d), 500.0)
    # x^4 e^x/(e^x-1)^2 = x^4 / (4 sinh^2(x/2)); the sinh form avoids overflow
    dBdT = x**4 / (4.0 * np.sinh(np.minimum(0.5 * x, 350.0))**2) / l**2
    return trapz(dBdT, lam) / trapz(dBdT / kap, lam)

def kappa_planck(T_d, n=6000):
    """Planck mean dust opacity [cm^2/g gas]: the emission/absorption mean.
    Use for the grain thermal balance, NOT for the trapping optical depth.
    Averages cabs_t, not ctr_t: scattered light never entered a grain, so it
    neither heats it nor is re-emitted by it.  (Averaging the transport opacity
    here runs 1.0% high at 300 K, 27% at 800 K and 66% at 1500 K.)"""
    lam = np.logspace(np.log10(lam_t.min()), np.log10(lam_t.max()), n)
    kap = np.interp(lam, lam_t, cabs_t) / M_GAS_PER_H
    B = planck_weight(lam, T_d, 'E')
    return trapz(B * kap, lam) / trapz(B, lam)

def band_mean_energy(e0, e1, T, sp=None):
    """<E> over the band [eV].  With sp: sigma-weighted, i.e. the mean energy of
    the photons that actually ionize sp -> excess heat = <E>_sigma - E_th."""
    if sp is not None and e1 <= E_TH[sp] * (1.0 + 1e-9):
        return 0.0
    E = _egrid(e0, e1, sp)
    w = planck_E(E, T, 'N')
    if sp is not None:
        w = w * sigma_species(sp, E)
        if trapz(w, E) <= 0.0:
            return 0.0
    return trapz(w * E, E) / trapz(w, E)

def legacy_regression():
    """Re-derive the ORIGINAL 7-band numbers committed in star_radiation.c.

    Independent of ION_BANDS/SPECIES: the old band edges are passed explicitly,
    so this keeps working after the split.  It proves the weighting, the dust
    tables and the Verner fits are untouched, which is what lets you trust the
    NEW numbers from the same run.  Expected values are the ones that were in
    star_radiation.c before the split (git show HEAD:src/stars/star_radiation.c).
    """
    HI_OLD = (13.60, 24.59)          # the unsplit ionizing band
    checks = [
        ("Sigma_N[HI][HI]",   band_sigma("HI", *HI_OLD, T_SOURCE, 'N'),   3.6742e-18),
        ("Sigma_E[HI][HI]",   band_sigma("HI", *HI_OLD, T_SOURCE, 'E'),   3.4457e-18),
        ("Sigma_N[HeI][HI]",  band_sigma("HI",   24.59, 54.42, T_SOURCE, 'N'), 8.4911e-19),
        ("Sigma_N[HeI][HeI]", band_sigma("HeI",  24.59, 54.42, T_SOURCE, 'N'), 5.8894e-18),
        ("Sigma_E[HeII][HeII]", band_sigma("HeII", 54.42, 200.0, T_SOURCE, 'E'), 1.3294e-18),
        # The dust checks average ctr_t, the UNFLOORED C_ext*(1-a<g>) the
        # committed tables were built from, NOT the catt_t the run now emits.
        # That is deliberate: it keeps this a test of the weighting and the
        # tables rather than a test of the g floor, which is a change of
        # definition and is MEANT to move the emitted numbers.  gfloor_check()
        # below reports by how much, band by band.
        ("Kappa_N[HI]",  dust_avg_E(ctr_t, *HI_OLD, T_SOURCE, 'N') / M_GAS_PER_H, 917.2),
        ("Kappa_E[HI]",  dust_avg_E(ctr_t, *HI_OLD, T_SOURCE, 'E') / M_GAS_PER_H, 902.4),
        ("Kappa_E[LW]",  band_avg(ctr_t, EV_UM / 13.6, EV_UM / 11.2,
                                  T_SOURCE, 'E') / M_GAS_PER_H, 736.6),
        ("Kappa_E[IR]",  band_avg(ctr_t, EV_UM / 1.0, 12.3984,
                                  T_SOURCE, 'E') / M_GAS_PER_H, 34.9),
    ]
    print("Regression against the committed (pre-split) tables")
    worst_rel, bad = 0.0, 0
    for name, got, want in checks:
        rel = abs(got - want) / abs(want)
        worst_rel = max(worst_rel, rel)
        ok = rel < 5e-4          # the committed values are rounded for printing
        bad += not ok
        print("  %-22s got %12.5g  want %12.5g  %s"
              % (name, got, want, "ok" if ok else "*** MISMATCH ***"))
    print("  %d/%d match, worst relative deviation %.2e"
          % (len(checks) - bad, len(checks), worst_rel))
    if bad:
        raise SystemExit("regression FAILED: the generator no longer reproduces "
                         "the committed tables, so its new output is not "
                         "trustworthy. Fix this before pasting anything.")
    print()


def gfloor_check():
    """Verify the three opacities and show what the g floor actually changes.

    cabs <= catt <= cext is the structural invariant: a beam cannot lose more
    than extinction, and cannot heat grains by more than it absorbs.  ctr is
    the one that is allowed to exceed cext, and only where g < 0.  If any of
    this trips, catt/ctr have been wired up wrong and every Kappa_* and
    MomentumFraction below is suspect.
    """
    eps = 1e-12 * cext_t
    assert np.all(cabs_t <= catt_t + eps), "cabs > catt"
    assert np.all(catt_t <= cext_t + eps), "catt > cext: the g floor is not applied"
    assert np.allclose(catt_t, np.minimum(ctr_t, cext_t), rtol=1e-12), \
        "catt != min(ctr, cext)"
    neg = g_t < 0.0
    assert np.all(ctr_t[~neg] <= cext_t[~neg] + eps[~neg]), "ctr > cext where g >= 0"

    print("g floor: <g> < 0 in %d of %d table rows, %.3f-%.3f um"
          % (neg.sum(), len(g_t), lam_t[neg].min(), lam_t[neg].max()))
    r = ctr_t[neg] / cext_t[neg]
    print("  max kappa_tr/kappa_ext = %.5f at %.3f um"
          % (r.max(), lam_t[neg][np.argmax(r)]))
    print("  effect on the emitted bands (kappa_E, cm^2/g gas):")
    for name, l0, l1 in DUST_BANDS:
        ktr  = band_avg(ctr_t,  l0, l1, T_SOURCE, 'E') / M_GAS_PER_H
        katt = band_avg(catt_t, l0, l1, T_SOURCE, 'E') / M_GAS_PER_H
        print("    %-13s kappa_tr = %8.3f  kappa_att = %8.3f  (%+.3f%%)"
              % (name, ktr, katt, 100.0 * (katt / ktr - 1.0)))
    for name, e0, e1 in ION_BANDS:
        ktr  = dust_avg_E(ctr_t,  e0, e1, T_SOURCE, 'E') / M_GAS_PER_H
        katt = dust_avg_E(catt_t, e0, e1, T_SOURCE, 'E') / M_GAS_PER_H
        print("    %-13s kappa_tr = %8.3f  kappa_att = %8.3f  (%+.3f%%)"
              % (name, ktr, katt, 100.0 * (katt / ktr - 1.0)))
    print()


if __name__ == "__main__":
    legacy_regression()
    gfloor_check()
    print("T_source = %.0f K, kappas in cm^2/g GAS at solar Z, sigmas in cm^2\n" % T_SOURCE)
    kN, kE, fabs, fmom = {}, {}, {}, {}
    for name, l0, l1 in DUST_BANDS:
        kN[name] = band_avg(catt_t, l0, l1, T_SOURCE, 'N') / M_GAS_PER_H
        kE[name] = band_avg(catt_t, l0, l1, T_SOURCE, 'E') / M_GAS_PER_H
        kabs     = band_avg(cabs_t, l0, l1, T_SOURCE, 'E') / M_GAS_PER_H
        ktr      = band_avg(ctr_t,  l0, l1, T_SOURCE, 'E') / M_GAS_PER_H
        fabs[name] = kabs / kE[name]
        fmom[name] = ktr / kE[name]
        print("%-13s kappa_N = %6.1f   kappa_E = %6.1f   f_abs(E) = %.2f"
              "   f_mom(E) = %.4f"
              % (name, kN[name], kE[name], fabs[name], fmom[name]))
    # IR is really weighted by the reprocessed dust spectrum, not the star:
    kIR100 = band_avg(catt_t, EV_UM/1.0, 12.3984, 100.0, 'E') / M_GAS_PER_H
    print("  (note: INFRARED re-weighted with a 100 K Planck gives kappa_E = "
          "%.2f; the table is truncated at 12.4 um so this is a lower bound)"
          % kIR100)

    # --- 3x3 (band, species) cross sections ---------------------------------
    sN = {b: {s: band_sigma(s, e0, e1, T_SOURCE, 'N') for s in SPECIES}
          for b, e0, e1 in ION_BANDS}
    sE = {b: {s: band_sigma(s, e0, e1, T_SOURCE, 'E') for s in SPECIES}
          for b, e0, e1 in ION_BANDS}
    eS = {b: {s: band_mean_energy(e0, e1, T_SOURCE, s) for s in SPECIES}
          for b, e0, e1 in ION_BANDS}
    eB = {b: band_mean_energy(e0, e1, T_SOURCE) for b, e0, e1 in ION_BANDS}

    print("\nBand-averaged sigma [cm^2], rows = band, cols = absorbing species")
    for tag, S in [("photon-number weighted (use for ionization rates)", sN),
                   ("energy weighted (use for band opacity/attenuation)", sE)]:
        print("\n  %s" % tag)
        print("  %-17s %s   <E>_band"
              % ("band", " ".join("%11s" % s for s in SPECIES)))
        for b, _, _ in ION_BANDS:
            print("  %-17s %s   %6.2f eV"
                  % (b, " ".join("%11.4e" % S[b][s] for s in SPECIES), eB[b]))

    # --- dust in the ionizing bands -----------------------------------------
    # Same catt_t / ctr_t / cabs_t split as the four dust bands: no special
    # case.  Here the g floor is inert -- g -> 1 across the whole EUV, so
    # catt == ctr exactly and f_mom comes out 1 -- but it costs nothing to run
    # them through the same path, and an SED with an X-ray tail would reach
    # table rows the floor does touch.  g -> 1 is physically right (Rayleigh-
    # Gans small-angle scattering, Draine 2003c), and the form self-corrects,
    # returning C_abs in that limit.  Do NOT reuse these g values as a Henyey-
    # Greenstein phase function -- g ~ 0.9999 is degenerate for direction
    # sampling.  As a scalar transport weight it is well conditioned; nothing
    # here divides by (1-g).
    for name, e0, e1 in ION_BANDS:
        kN[name] = dust_avg_E(catt_t, e0, e1, T_SOURCE, 'N') / M_GAS_PER_H
        kE[name] = dust_avg_E(catt_t, e0, e1, T_SOURCE, 'E') / M_GAS_PER_H
        kabs     = dust_avg_E(cabs_t, e0, e1, T_SOURCE, 'E') / M_GAS_PER_H
        ktr      = dust_avg_E(ctr_t,  e0, e1, T_SOURCE, 'E') / M_GAS_PER_H
        fabs[name] = kabs / kE[name]
        fmom[name] = ktr / kE[name]
    print("\nDust in the ionizing bands (kappa in cm^2/g gas at solar Z)")
    print("  %-14s %10s %10s %8s %9s   %s"
          % ("band", "kappa_N", "kappa_E", "f_abs(E)", "f_mom(E)", "sigma_d/H"))
    for name, e0, e1 in ION_BANDS:
        print("  %-14s %10.1f %10.1f %8.2f %9.4f   %.3e"
              % (name, kN[name], kE[name], fabs[name], fmom[name],
                 kN[name] * M_GAS_PER_H))

    print("\nMean excess energy <E>_sigma - E_th [eV] deposited as heat")
    print("  %-17s %s" % ("band", " ".join("%11s" % s for s in SPECIES)))
    for b, _, _ in ION_BANDS:
        row = [(eS[b][s] - E_TH[s]) if sN[b][s] > 0.0 else 0.0 for s in SPECIES]
        print("  %-17s %s" % (b, " ".join("%11.3f" % v for v in row)))

    # --- photoheating headroom ------------------------------------------------
    # radiation_feedback() in star_radiation.c forms the photoheating as
    #   E_exc = E_abs - N_abs * E_th
    # and warns (unrate-limited!) if it comes out negative.  The margin that
    # protects it is <E>_sigma / E_th.  Splitting IONIZING_HI shrinks the soft
    # band's margin a lot, because its upper edge is only 1.8 eV above the HI
    # threshold.  Two biases in the RT eat into it: the channel attribution is
    # photon/energy asymmetric (Sigma_N/Sigma_E vs Kappa_N/Kappa_E), and
    # Ch[].Energy is reduced by momentum recoil while Ch[].Photons is not.
    # Print the margin so the risk is visible before anything is pasted.
    print("\nPhotoheating headroom <E>_sigma / E_th  (1.0 = warn() fires)")
    print("  %-17s %s" % ("band", " ".join("%11s" % s for s in SPECIES)))
    worst = (1e30, None, None)
    for b, _, _ in ION_BANDS:
        cells = []
        for s in SPECIES:
            if sN[b][s] > 0.0:
                r = eS[b][s] / E_TH[s]
                cells.append("%11.3f" % r)
                if r < worst[0]:
                    worst = (r, b, s)
            else:
                cells.append("%11s" % "-")
        print("  %-17s %s" % (b, " ".join(cells)))
    print("  tightest: %s / %s = %.3f  (%+.1f%% headroom)"
          % (worst[1], worst[2], worst[0], 100.0 * (worst[0] - 1.0)))
    if worst[0] < 1.10:
        print("  WARNING: under 10%% headroom. Clamp E_exc at 0 and rate-limit")
        print("           the warn() in radiation_feedback(), or switch that")
        print("           channel's heating to the DeltaE_eV table below.")

    # --- how much of the H2 answer rests on the weak branch -------------------
    if "H2" in SPECIES:
        print("\nH2 sigma provenance check")
        print("  sigma_h2_ion(15.0 eV)  = %.4e cm^2  (must be 0)"
              % sigma_h2_ion(15.0))
        for e in (E_TH_H2, 16.0, 17.0, 17.999, 18.0, 20.0, 24.59, 30.0, 50.0, 85.0):
            br = "B interp" if e < H2_BRANCH_A_LO else "A YSD98"
            if e > H2_BRANCH_A_HI:
                br = "A extrap"
            print("  sigma_h2_ion(%8.3f eV) = %6.3f Mb   [%s]"
                  % (e, sigma_h2_ion(e) / MB, br))
        jump = sigma_h2_ion(18.0) / sigma_h2_ion(17.999)
        print("  branch A/B ratio at 18 eV = %.4f  (1.0 = continuous)" % jump)

        # Fraction of the H2 band's ionization rate coming from the weak,
        # interpolated 15.4259-18 eV sub-range.  This is the honest error bar.
        for bname, e0, e1 in ION_BANDS:
            if sN[bname]["H2"] <= 0.0:
                continue
            E = _egrid(e0, e1, "H2")
            w = planck_E(E, T_SOURCE, 'N') * sigma_h2_ion(E)
            tot = trapz(w, E)
            if tot <= 0.0:
                continue
            m = E < H2_BRANCH_A_LO
            weak = trapz(np.where(m, w, 0.0), E)
            print("  %-17s %5.1f%% of H2 ionizations come from the"
                  " interpolated <18 eV branch" % (bname, 100.0 * weak / tot))
        print("  -> quote sigma_H2 to no better than ~a factor 1.5 until the")
        print("     YSD98 erratum replaces branch B.")

    print("\n/* ---- paste into star_radiation.c ---- */")
    print("/* D03 MW R_V=3.1 dust, %.0f kK Planck weights.  Beam-depletion opacity"
          % (T_SOURCE/1e3))
    print("   kappa_att = kappa_ext*(1 - a*max(<g>,0)) = min(kappa_tr, kappa_ext): the")
    print("   delta-scaled transport opacity where <g> >= 0, capped at true extinction")
    print("   where <g> < 0, since a beam cannot lose more than kappa_ext.  Momentum is")
    print("   NOT this number -- see MomentumFraction. */")
    # Emitted E-then-N to match the declaration order in star_radiation.c.
    # These two differ by only 1-3% in the ionizing bands, so a swapped paste
    # produces plausible numbers and no symptom -- it has happened once already.
    # Paste whole blocks, do not hand-edit individual entries.
    for tag, K in [("Kappa_E", kE), ("Kappa_N", kN)]:
        print("double %s[WAVEBANDS] =\n{" % tag)
        for name, _, _ in DUST_BANDS:
            print("  [%s] = %.1f," % (name, K[name]))
        for name, _, _ in ION_BANDS:
            print("  [%s] = %.1f, /* dust only; gas is per-species, see sigma */"
                  % (name, K[name]))
        print("};")

    zeros = "{ " + ", ".join(["0.0"] * len(SPECIES)) + " }"
    print("\n/* Band-averaged photoionization cross sections [cm^2].")
    print("   Row = waveband, column = absorbing species.  Zero where the band")
    print("   lies below that species' threshold.")
    print("   HI/HeI/HeII: Verner, Ferland, Korista & Yakovlev (1996).")
    print("   H2: see sigma_h2_ion() in tools/band_kappas.py -- the 15.43-18 eV")
    print("   branch is an interpolation, NOT the published fit, and it is the")
    print("   dominant uncertainty.  Total ionization cross section, so k29 also")
    print("   absorbs the dissociative-ionization channel above 18.08 eV.")
    print("   %.0f kK Planck weights, same as Kappa_* above. */" % (T_SOURCE / 1e3))
    print("enum { %s, N_ION_SPECIES };"
          % ", ".join("SP_" + s for s in SPECIES))
    for tag, S in [("Sigma_N", sN), ("Sigma_E", sE)]:
        print("double %s[WAVEBANDS][N_ION_SPECIES] =\n{" % tag)
        for name, _, _ in DUST_BANDS:
            print("  [%s] = %s," % (name, zeros))
        for name, _, _ in ION_BANDS:
            print("  [%s] = { %s },"
                  % (name, ", ".join("%.4e" % S[name][s] for s in SPECIES)))
        print("};")
    # --- energy-partition fractions -----------------------------------------
    # eps_pe is a fraction of GRAIN-ABSORBED energy, so it multiplies f_abs
    # rather than being subtracted from it: light that scattered never entered
    # a grain and cannot eject a photoelectron.  Only the FUV bands feed the
    # photoelectric channel (6-13.6 eV, the field Grackle's PE heating expects);
    # EUV-absorbed energy has no PE channel to enter, so it is all reradiated.
    EPS_PE = {"ULTRAVIOLET": 0.05, "LYMAN_WERNER": 0.05}
    ALL_BANDS = [b[0] for b in DUST_BANDS] + [b[0] for b in ION_BANDS]
    print("\n/* f_abs = kappa_abs/kappa_att = (1-a)/(1-a*max(<g>,0)), D03 MW dust,")
    print("   band-averaged.  Remainder is scattered light: removed from the ray, still")
    print("   delivers momentum, but must NOT heat.  In the ionizing bands this is the")
    print("   DUST CHANNEL only; the gas channel is 1.0 absorbed / 0.0 reradiated by")
    print("   construction. */")
    print("double AbsorbedFraction[WAVEBANDS] =\n{")
    for name in ALL_BANDS:
        print("  [%s] = %.2f," % (name, fabs[name]))
    print("};")
    print("\n/* f_mom = kappa_tr/kappa_att: momentum transfer per unit energy removed")
    print("   from the beam.  Exactly 1 wherever <g> >= 0, because there the depletion")
    print("   opacity IS the transport opacity.  Exceeds 1 only where <g> < 0, where")
    print("   back-scattering delivers more than h*nu/c per photon removed; bounded")
    print("   above by 2 (retroreflection).  Only meaningful while <g> keeps one sign")
    print("   across the band. */")
    print("double MomentumFraction[WAVEBANDS] =\n{")
    for name in ALL_BANDS:
        print("  [%s] = %.4f," % (name, fmom[name]))
    print("};")
    print("\n/* f_rerad = f_abs*(1-eps_pe); eps_pe = 0.05 for the two FUV bands only */")
    print("double ReradiatedFraction[WAVEBANDS] =\n{")
    for name in ALL_BANDS:
        print("  [%s] = %.2f," % (name, fabs[name] * (1.0 - EPS_PE.get(name, 0.0))))
    print("};")

    print("\n/* Mean excess energy per photoionization [eV] -> photoheating.")
    print("   Not currently used: radiation_feedback() derives the heating from")
    print("   E_abs - N_abs*E_th instead.  This table is the robust fallback if")
    print("   that difference goes negative (see the headroom check above); it")
    print("   needs a per-band photon accumulator to use. */")
    print("double DeltaE_eV[WAVEBANDS][N_ION_SPECIES] =\n{")
    for name, _, _ in DUST_BANDS:
        print("  [%s] = %s," % (name, zeros))
    for name, _, _ in ION_BANDS:
        row = [(eS[name][s] - E_TH[s]) if sN[name][s] > 0.0 else 0.0 for s in SPECIES]
        print("  [%s] = { %s },"
              % (name, ", ".join("%.4f" % v for v in row)))
    print("};")

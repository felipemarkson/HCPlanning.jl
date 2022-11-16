module SystemData
n_bus = 138
n_branches = 151

load_factor = [0.7, 0.83, 1]

branch = [
    #    (s,r) length  type
    ((1, 136), 0.41, "EFF"),
    ((1, 2), 0.84, "EFF"),
    ((2, 3), 0.29, "EFF"),
    ((3, 4), 0.08, "EFF"),
    ((4, 5), 0.59, "EFF"),
    ((5, 6), 0.15, "EFF"),
    ((6, 8), 0.45, "EFF"),
    ((6, 7), 0.97, "EFF"),
    ((8, 10), 0.17, "EFF"),
    ((8, 9), 0.57, "EFF"),
    ((10, 13), 0.08, "EFF"),
    ((10, 12), 0.33, "EFF"),
    ((10, 11), 0.84, "NAF"),
    ((11, 121), 1.15, "EFF"),
    ((13, 15), 0.00, "EFF"),
    ((13, 14), 0.56, "EFF"),
    ((15, 121), 0.27, "NAF"),
    ((15, 120), 0.33, "NAF"),
    ((15, 16), 0.35, "NAF"),
    ((16, 120), 0.62, "ERF"),
    ((17, 136), 0.21, "ERF"),
    ((17, 18), 0.08, "ERF"),
    ((18, 19), 0.34, "EFF"),
    ((19, 20), 0.41, "EFF"),
    ((20, 22), 0.80, "EFF"),
    ((20, 21), 0.46, "EFF"),
    ((22, 24), 0.34, "EFF"),
    ((22, 23), 0.06, "EFF"),
    ((24, 25), 0.11, "EFF"),
    ((25, 26), 0.05, "EFF"),
    ((26, 27), 0.16, "EFF"),
    ((27, 28), 0.03, "EFF"),
    ((28, 31), 0.21, "EFF"),
    ((28, 29), 0.45, "NAF"),
    ((29, 30), 0.08, "EFF"),
    ((30, 113), 0.61, "NAF"),
    ((31, 123), 0.11, "EFF"),
    ((31, 32), 0.07, "EFF"),
    ((32, 33), 0.47, "ERF"),
    ((33, 34), 1.13, "EFF"),
    ((35, 136), 0.31, "ERF"),
    ((35, 36), 0.19, "ERF"),
    ((36, 38), 0.11, "EFF"),
    ((36, 37), 0.37, "EFF"),
    ((38, 39), 0.26, "EFF"),
    ((39, 46), 0.42, "EFF"),
    ((39, 40), 0.30, "EFF"),
    ((40, 41), 0.08, "EFF"),
    ((41, 42), 0.84, "EFF"),
    ((42, 51), 0.30, "EFF"),
    ((42, 43), 3.13, "EFF"),
    ((43, 52), 0.00, "NAF"),
    ((43, 44), 0.17, "EFF"),
    ((44, 45), 0.92, "NAF"),
    ((45, 108), 0.54, "EFF"),
    ((46, 129), 0.16, "EFF"),
    ((46, 47), 0.97, "EFF"),
    ((47, 48), 0.33, "NAF"),
    ((48, 49), 0.16, "NAF"),
    ((49, 50), 0.29, "EFF"),
    ((50, 101), 0.64, "EFF"),
    ((51, 101), 0.88, "EFF"),
    ((52, 53), 1.40, "EFF"),
    ((53, 57), 0.40, "EFF"),
    ((53, 54), 0.11, "NAF"),
    ((54, 55), 0.07, "EFF"),
    ((55, 56), 0.15, "EFF"),
    ((56, 107), 1.22, "NAF"),
    ((57, 58), 0.08, "EFF"),
    ((58, 59), 0.05, "EFF"),
    ((59, 104), 1.21, "EFF"),
    ((60, 137), 0.27, "EFF"),
    ((60, 61), 0.00, "EFF"),
    ((61, 62), 0.29, "NAF"),
    ((62, 63), 0.47, "NAF"),
    ((63, 64), 0.57, "EFF"),
    ((64, 65), 0.37, "EFF"),
    ((65, 133), 0.54, "EFF"),
    ((65, 131), 0.81, "EFF"),
    ((66, 137), 0.24, "EFF"),
    ((66, 67), 0.56, "EFF"),
    ((67, 68), 0.35, "EFF"),
    ((68, 69), 0.02, "NAF"),
    ((69, 70), 1.68, "EFF"),
    ((70, 71), 0.96, "EFF"),
    ((71, 72), 0.83, "EFF"),
    ((72, 119), 0.83, "EFF"),
    ((72, 73), 0.37, "EFF"),
    ((74, 137), 0.43, "EFF"),
    ((74, 75), 0.25, "EFF"),
    ((75, 76), 0.63, "EFF"),
    ((76, 78), 0.15, "EFF"),
    ((76, 77), 0.79, "EFF"),
    ((78, 80), 1.16, "EFF"),
    ((78, 79), 0.17, "NAF"),
    ((80, 82), 1.40, "ERF"),
    ((80, 81), 0.03, "EFF"),
    ((82, 83), 1.84, "ERF"),
    ((83, 84), 0.57, "ERF"),
    ((84, 85), 0.52, "EFF"),
    ((85, 111), 0.31, "EFF"),
    ((86, 137), 0.34, "EFF"),
    ((86, 87), 0.84, "EFF"),
    ((87, 89), 0.12, "NAF"),
    ((87, 88), 0.02, "NAF"),
    ((89, 90), 0.42, "EFF"),
    ((90, 91), 0.88, "EFF"),
    ((91, 92), 0.48, "EFF"),
    ((92, 93), 0.64, "EFF"),
    ((93, 94), 0.15, "EFF"),
    ((94, 110), 0.03, "NAF"),
    ((94, 106), 1.63, "NAF"),
    ((95, 137), 0.64, "NAF"),
    ((95, 96), 0.39, "NAF"),
    ((96, 97), 0.03, "NAF"),
    ((97, 99), 1.05, "NAF"),
    ((97, 98), 0.17, "NAF"),
    ((99, 100), 1.10, "NAF"),
    ((100, 128), 0.11, "NAF"),
    ((101, 102), 0.56, "NAF"),
    ((102, 103), 1.05, "NAF"),
    ((103, 138), 0.30, "NAF"),
    ((104, 105), 0.60, "NAF"),
    ((105, 122), 0.25, "NAF"),
    ((106, 107), 0.30, "NAF"),
    ((107, 125), 0.35, "NAF"),
    ((107, 122), 1.11, "NAF"),
    ((108, 138), 0.23, "NAF"),
    ((108, 122), 0.41, "NAF"),
    ((108, 109), 2.42, "NAF"),
    ((109, 110), 0.53, "NAF"),
    ((111, 112), 0.56, "NAF"),
    ((112, 113), 0.41, "NAF"),
    ((113, 114), 0.92, "NAF"),
    ((114, 123), 1.70, "NAF"),
    ((115, 129), 0.28, "NAF"),
    ((115, 116), 0.90, "NAF"),
    ((116, 117), 1.25, "NAF"),
    ((117, 118), 0.86, "NAF"),
    ((118, 130), 0.56, "EFF"),
    ((119, 120), 0.64, "ERF"),
    ((121, 132), 0.96, "ERF"),
    ((123, 124), 0.58, "EFF"),
    ((124, 125), 0.36, "EFF"),
    ((126, 138), 1.36, "EFF"),
    ((126, 127), 1.15, "EFF"),
    ((127, 128), 0.36, "ERF"),
    ((130, 138), 1.75, "NAF"),
    ((132, 133), 0.31, "NAF"),
    ((133, 134), 0.45, "NAF"),
    ((134, 135), 0.56, "NAF"),
]

peak_demand = [
    #Stage1 2       3       4       5       6       7       8       8       10
    259.00 271.95 284.90 297.85 310.80 323.75 336.70 349.65 362.60 375.55
    51.42 53.99 56.56 59.14 61.71 64.28 66.85 69.42 71.99 74.56
    45.79 48.08 50.37 52.66 54.95 57.24 59.53 61.82 64.11 66.40
    93.66 98.34 103.02 107.70 112.39 117.07 121.75 126.44 131.12 135.80
    335.05 351.80 368.55 385.30 402.05 418.81 435.56 452.31 469.06 485.82
    160.22 168.23 176.24 184.25 192.26 200.27 208.29 216.30 224.31 232.32
    256.87 269.71 282.55 295.40 308.24 321.08 333.93 346.77 359.61 372.46
    67.05 70.40 73.75 77.11 80.46 83.81 87.16 90.52 93.87 97.22
    134.09 140.79 147.50 154.20 160.91 167.61 174.32 181.02 187.73 194.43
    150.85 158.39 165.93 173.47 181.02 188.56 196.10 203.64 211.19 218.73
    125.72 132.00 138.29 144.57 150.86 157.14 163.43 169.72 176.00 182.29
    268.20 281.61 295.02 308.43 321.84 335.25 348.66 362.07 375.48 388.89
    313.53 329.21 344.89 360.57 376.24 391.92 407.60 423.27 438.95 454.63
    326.88 343.22 359.56 375.91 392.25 408.60 424.94 441.28 457.63 473.97
    231.81 243.40 254.99 266.58 278.17 289.76 301.35 312.95 324.54 336.13
    213.73 224.42 235.10 245.79 256.47 267.16 277.85 288.53 299.22 309.91
    429.00 450.45 471.90 493.35 514.80 536.25 557.70 579.15 600.60 622.05
    158.00 165.90 173.80 181.70 189.60 197.50 205.40 213.30 221.20 229.10
    75.00 78.75 82.50 86.25 90.00 93.75 97.50 101.25 105.00 108.75
    33.53 35.21 36.89 38.56 40.24 41.92 43.60 45.27 46.95 48.63
    257.10 269.95 282.80 295.66 308.51 321.37 334.22 347.08 359.93 372.79
    67.07 70.43 73.78 77.13 80.49 83.84 87.19 90.55 93.90 97.25
    257.10 269.95 282.80 295.66 308.51 321.37 334.22 347.08 359.93 372.79
    134.14 140.85 147.55 154.26 160.97 167.68 174.38 181.09 187.80 194.50
    269.00 282.45 295.90 309.35 322.80 336.25 349.70 363.15 376.60 390.05
    63.43 66.60 69.77 72.94 76.11 79.28 82.45 85.62 88.80 91.97
    405.92 426.21 446.51 466.80 487.10 507.40 527.69 547.99 568.28 588.58
    216.00 226.80 237.60 248.40 259.20 270.00 280.80 291.60 302.40 313.20
    138.75 145.69 152.62 159.56 166.50 173.44 180.37 187.31 194.25 201.19
    63.43 66.60 69.77 72.94 76.11 79.28 82.45 85.62 88.80 91.97
    377.00 395.85 414.70 433.55 452.40 471.25 490.10 508.95 527.80 546.65
    95.14 99.90 104.65 109.41 114.17 118.93 123.68 128.44 133.20 137.95
    171.00 179.55 188.10 196.65 205.20 213.75 222.30 230.85 239.40 247.95
    441.61 463.69 485.78 507.86 529.94 552.02 574.10 596.18 618.26 640.34
    373.00 391.65 410.30 428.95 447.60 466.25 484.90 503.55 522.20 540.85
    201.64 211.72 221.80 231.89 241.97 252.05 262.13 272.21 282.30 292.38
    269.56 283.04 296.52 309.99 323.47 336.95 350.43 363.91 377.38 390.86
    83.83 88.03 92.22 96.41 100.60 104.79 108.99 113.18 117.37 121.56
    143.00 150.15 157.30 164.45 171.60 178.75 185.90 193.05 200.20 207.35
    1.36 1.43 1.50 1.57 1.63 1.70 1.77 1.84 1.91 1.97
    6.81 7.16 7.50 7.84 8.18 8.52 8.86 9.20 9.54 9.88
    171.00 179.55 188.10 196.65 205.20 213.75 222.30 230.85 239.40 247.95
    128.03 134.44 140.84 147.24 153.64 160.04 166.44 172.85 179.25 185.65
    67.68 71.07 74.45 77.83 81.22 84.60 87.99 91.37 94.76 98.14
    187.13 196.49 205.84 215.20 224.56 233.91 243.27 252.63 261.98 271.34
    498.06 522.96 547.87 572.77 597.67 622.58 647.48 672.38 697.29 722.19
    285.61 299.89 314.17 328.45 342.73 357.01 371.29 385.57 399.86 414.14
    256.07 268.87 281.68 294.48 307.28 320.09 332.89 345.69 358.50 371.30
    204.00 214.20 224.40 234.60 244.80 255.00 265.20 275.40 285.60 295.80
    118.63 124.56 130.49 136.42 142.35 148.28 154.22 160.15 166.08 172.01
    289.00 303.45 317.90 332.35 346.80 361.25 375.70 390.15 404.60 419.05
    79.08 83.04 86.99 90.94 94.90 98.85 102.81 106.76 110.71 114.67
    280.74 294.77 308.81 322.85 336.88 350.92 364.96 378.99 393.03 407.07
    75.13 78.88 82.64 86.40 90.15 93.91 97.67 101.42 105.18 108.93
    23.72 24.91 26.10 27.28 28.47 29.66 30.84 32.03 33.21 34.40
    442.00 464.10 486.20 508.30 530.40 552.50 574.60 596.70 618.80 640.90
    22.30 23.41 24.52 25.64 26.75 27.87 28.98 30.10 31.21 32.33
    163.52 171.69 179.87 188.05 196.22 204.40 212.57 220.75 228.93 237.10
    239.70 251.68 263.67 275.65 287.64 299.62 311.61 323.59 335.58 347.56
    100.34 105.36 110.38 115.39 120.41 125.43 130.44 135.46 140.48 145.50
    425.00 446.25 467.50 488.75 510.00 531.25 552.50 573.75 595.00 616.25
    246.22 258.53 270.84 283.15 295.46 307.77 320.08 332.39 344.70 357.02
    365.00 383.25 401.50 419.75 438.00 456.25 474.50 492.75 511.00 529.25
    316.43 332.25 348.08 363.90 379.72 395.54 411.36 427.18 443.01 458.83
    89.34 93.81 98.28 102.75 107.21 111.68 116.15 120.61 125.08 129.55
    89.34 93.81 98.28 102.75 107.21 111.68 116.15 120.61 125.08 129.55
    111.68 117.27 122.85 128.43 134.02 139.60 145.19 150.77 156.35 161.94
    189.86 199.35 208.84 218.34 227.83 237.32 246.82 256.31 265.80 275.29
    89.34 93.81 98.28 102.75 107.21 111.68 116.15 120.61 125.08 129.55
    234.53 246.26 257.99 269.71 281.44 293.17 304.89 316.62 328.35 340.07
    25.07 26.32 27.58 28.83 30.08 31.34 32.59 33.84 35.10 36.35
    5.46 5.73 6.01 6.28 6.55 6.83 7.10 7.37 7.65 7.92
    78.18 82.08 85.99 89.90 93.81 97.72 101.63 105.54 109.45 113.35
    436.94 458.79 480.63 502.48 524.33 546.18 568.02 589.87 611.72 633.56
    185.00 194.25 203.50 212.75 222.00 231.25 240.50 249.75 259.00 268.25
    108.81 114.25 119.69 125.13 130.57 136.01 141.45 146.89 152.33 157.77
    154.80 162.54 170.28 178.02 185.76 193.50 201.24 208.98 216.72 224.46
    104.31 109.53 114.75 119.96 125.18 130.39 135.61 140.83 146.04 151.26
    326.33 342.65 358.97 375.28 391.60 407.92 424.23 440.55 456.87 473.18
    153.41 161.08 168.75 176.42 184.09 191.76 199.43 207.10 214.77 222.44
    303.96 319.15 334.35 349.55 364.75 379.94 395.14 410.34 425.54 440.74
    94.83 99.57 104.32 109.06 113.80 118.54 123.28 128.02 132.77 137.51
    264.86 278.10 291.34 304.58 317.83 331.07 344.31 357.55 370.80 384.04
    269.09 282.55 296.00 309.46 322.91 336.37 349.82 363.28 376.73 390.19
    424.00 445.20 466.40 487.60 508.80 530.00 551.20 572.40 593.60 614.80
    97.62 102.50 107.38 112.26 117.14 122.03 126.91 131.79 136.67 141.55
    1235.27 1297.03 1358.79 1420.56 1482.32 1544.08 1605.85 1667.61 1729.37 1791.14
    497.82 522.71 547.61 572.50 597.39 622.28 647.17 672.06 696.95 721.84
    418.38 439.30 460.22 481.14 502.06 522.98 543.90 564.81 585.73 606.65
    67.00 70.35 73.70 77.05 80.40 83.75 87.10 90.45 93.80 97.15
    86.47 90.79 95.11 99.44 103.76 108.08 112.41 116.73 121.05 125.38
    94.83 99.57 104.32 109.06 113.80 118.54 123.28 128.02 132.77 137.51
    260.00 273.00 286.00 299.00 312.00 325.00 338.00 351.00 364.00 377.00
    80.38 84.39 88.41 92.43 96.45 100.47 104.49 108.51 112.53 116.54
    252.04 264.64 277.24 289.84 302.45 315.05 327.65 340.25 352.85 365.46
    154.04 161.74 169.44 177.14 184.84 192.55 200.25 207.95 215.65 223.35
    67.00 70.35 73.70 77.05 80.40 83.75 87.10 90.45 93.80 97.15
    83.03 87.19 91.34 95.49 99.64 103.79 107.94 112.10 116.25 120.40
    363.00 381.15 399.30 417.45 435.60 453.75 471.90 490.05 508.20 526.35
    55.74 58.53 61.32 64.10 66.89 69.68 72.47 75.25 78.04 80.83
    65.03 68.28 71.53 74.79 78.04 81.29 84.54 87.79 91.04 94.30
    130.00 136.50 143.00 149.50 156.00 162.50 169.00 175.50 182.00 188.50
    224.00 235.20 246.40 257.60 268.80 280.00 291.20 302.40 313.60 324.80
    0.00 18.18 19.09 19.99 20.90 21.81 22.72 23.63 24.54 25.45
    0.00 1636.27 1718.09 1799.90 1881.71 1963.53 2045.34 2127.16 2208.97 2290.78
    0.00 339.98 356.98 373.98 390.98 407.98 424.98 441.98 458.98 475.98
    0.00 86.71 91.04 95.38 99.71 104.05 108.38 112.72 117.06 121.39
    0.00 0.00 55.74 58.53 61.32 64.10 66.89 69.68 72.47 75.25
    0.00 0.00 73.00 76.65 80.30 83.95 87.60 91.25 94.90 98.55
    0.00 0.00 219.88 230.87 241.87 252.86 263.85 274.85 285.84 296.83
    0.00 0.00 0.00 66.10 69.40 72.71 76.01 79.32 82.62 85.93
    0.00 0.00 0.00 49.55 52.02 54.50 56.98 59.46 61.93 64.41
    0.00 0.00 0.00 152.00 159.60 167.20 174.80 182.40 190.00 197.60
    0.00 0.00 0.00 170.60 179.13 187.66 196.19 204.72 213.25 221.78
    0.00 0.00 0.00 0.00 117.00 122.85 128.70 134.55 140.40 146.25
    0.00 0.00 0.00 0.00 271.70 285.28 298.87 312.45 326.04 339.62
    0.00 0.00 0.00 0.00 355.00 372.75 390.50 408.25 426.00 443.75
    0.00 0.00 0.00 0.00 74.51 78.24 81.96 85.69 89.42 93.14
    0.00 0.00 0.00 0.00 0.00 34.83 36.58 38.32 40.06 41.80
    0.00 0.00 0.00 0.00 0.00 66.35 69.66 72.98 76.30 79.61
    0.00 0.00 0.00 0.00 0.00 273.00 286.65 300.30 313.95 327.60
    0.00 0.00 0.00 0.00 0.00 0.00 105.32 110.59 115.86 121.12
    0.00 0.00 0.00 0.00 0.00 0.00 55.50 58.27 61.05 63.82
    0.00 0.00 0.00 0.00 0.00 0.00 137.09 143.95 150.80 157.66
    0.00 0.00 0.00 0.00 0.00 0.00 87.21 91.57 95.93 100.29
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 161.93 170.03 178.13
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 23.79 24.98 26.16
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 83.25 87.41 91.57
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 253.71 266.40
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 39.64 41.62
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 277.49 291.36
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 352.54
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 371.58
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 277.49
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 317.00
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 #add eq14 problem
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 #add eq14 problem
    0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 0.00 #add eq14 problem
]/1000

node_zone = [
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    2,
    2,
    2,
    3,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    2,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    2,
    2,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    1,
    1,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    2,
    3,
    3,
    3,
    3,
    3,
    3,
    2,
    2,
    2,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    3,
    3,
    3,
    3,
    2,
    2,
    3,
    1,
    1,
    1,
    1,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    3,
    2,
    3,
    3,
]

wind_speed = [
    #Load Level
    #1   2    3
    8.53 9.12 10.04 #Zone A
    6.13 7.26 7.11 #Zone B
    4.13 5.10 5.56 #Zone C
]
end  # moduleSystemData

#
# module AssetsData
#
#     Ωᵗʳ = Dict(
#     "ET" => [21, 22]
#     "NT" => [24, 23]
#     )
#
#     Ωᵖ = Dict(
#     "C" => [2 ,3 ,7, 13, 15, 16, 17, 20],
#     "W" => [1, 4, 5, 9 ,15, 17, 18, 19]
#     )
# end
#

#
# n_bus = 24
# n_brach = 33
# pf = 0.9
#
# i = 0.071
#
# load_level_perc = [0.7, 0.83, 1]
#
# load_hours_peryear = [2000, 5760, 1000]
#
# cost_non_server = 2000
#
# cost_generation_per_loadlevel = [57.7, 70, 85.3]

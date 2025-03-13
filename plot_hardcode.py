import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

esno_ldpc_20_10_classic = [-7, -6.50000000000000, -6, -5.50000000000000, -5, -4.50000000000000, -4, -3.50000000000000, -3, -2.50000000000000, -2, -1.50000000000000, -1, -0.500000000000000, 0, 0.500000000000000, 1, 1.50000000000000, 2, 2.50000000000000, 3, 3.50000000000000, 4, 4.50000000000000, 5, 5.50000000000000]
fer_ldpc_20_10_classic = [1, 0.909090909090909, 0.909090909090909, 1, 0.833333333333333, 0.666666666666667, 0.625000000000000, 0.625000000000000, 0.666666666666667, 0.588235294117647, 0.312500000000000, 0.333333333333333, 0.400000000000000, 0.163934426229508, 0.222222222222222, 0.172413793103448, 0.0833333333333333, 0.0492610837438424, 0.0342465753424658, 0.0200000000000000, 0.00733137829912024, 0.00788643533123028, 0.00282645562464669, 0.00149186931224825, 0.000486570650058389, 0.000344174840819136]
ber_ldpc_20_10_classic = [0.210000000000000, 0.231818181818182, 0.240909090909091, 0.265000000000000, 0.191666666666667, 0.106666666666667, 0.146875000000000, 0.131250000000000, 0.116666666666667, 0.102941176470588, 0.0640625000000000, 0.0566666666666667, 0.0700000000000000, 0.0262295081967213, 0.0344444444444444, 0.0206896551724138, 0.00958333333333333, 0.00615763546798030, 0.00547945205479452, 0.00280000000000000, 0.000989736070381232, 0.000867507886435331, 0.000423968343697004, 0.000216321050275996, 7.78513040093422e-05, 4.47427293064877e-05]

esno_ldpc_20_10_bcjr = [-6, -5.5, -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
fer_ldpc_20_10_bcjr = [1.0, 1.0, 1.0, 1.0, 0.9090909090909091, 0.9090909090909091, 0.9090909090909091, 0.7142857142857143, 0.8333333333333334, 0.625, 0.625, 0.5263157894736842, 0.3225806451612903, 0.4166666666666667, 0.24390243902439024, 0.23809523809523808, 0.10309278350515463, 0.056179775280898875, 0.05025125628140704, 0.036101083032490974, 0.02824858757062147, 0.019569471624266144, 0.008319467554076539, 0.0033647375504710633]
ber_ldpc_20_10_bcjr = [0.05, 0.05, 0.05, 0.05, 0.045454545454545456, 0.045454545454545456, 0.045454545454545456, 0.03571428571428571, 0.041666666666666664, 0.03125, 0.03125, 0.02631578947368421, 0.016129032258064516, 0.020833333333333332, 0.012195121951219513, 0.011904761904761904, 0.005154639175257732, 0.0028089887640449437, 0.002512562814070352, 0.0018050541516245488, 0.0014124293785310734, 0.0009784735812133072, 0.00041597337770382697, 0.00016823687752355316]


esno_ldpc_16_10_classic = [-7, -6.90000000000000, -6.80000000000000, -6.70000000000000, -6.60000000000000, -6.50000000000000, -6.40000000000000, -6.30000000000000, -6.20000000000000, -6.10000000000000, -6, -5.90000000000000, -5.80000000000000, -5.70000000000000, -5.60000000000000, -5.50000000000000, -5.40000000000000, -5.30000000000000, -5.20000000000000, -5.10000000000000, -5, -4.90000000000000, -4.80000000000000, -4.70000000000000, -4.60000000000000, -4.50000000000000, -4.40000000000000, -4.30000000000000, -4.20000000000000, -4.10000000000000, -4, -3.90000000000000, -3.80000000000000, -3.70000000000000, -3.60000000000000, -3.50000000000000, -3.40000000000000, -3.30000000000000, -3.20000000000000, -3.10000000000000, -3, -2.90000000000000, -2.80000000000000, -2.70000000000000, -2.60000000000000, -2.50000000000000, -2.40000000000000, -2.30000000000000, -2.20000000000000, -2.10000000000000, -2, -1.90000000000000, -1.80000000000000, -1.70000000000000, -1.60000000000000, -1.50000000000000, -1.40000000000000, -1.30000000000000, -1.20000000000000, -1.10000000000000, -1, -0.900000000000000, -0.800000000000000, -0.699999999999999, -0.600000000000000, -0.500000000000000, -0.399999999999999, -0.300000000000000, -0.199999999999999, -0.100000000000001, 0, 0.0999999999999996, 0.200000000000000, 0.300000000000000, 0.399999999999999, 0.500000000000000, 0.600000000000000, 0.700000000000000, 0.800000000000000, 0.900000000000000, 1, 1.10000000000000, 1.20000000000000, 1.30000000000000, 1.40000000000000, 1.50000000000000, 1.60000000000000, 1.70000000000000, 1.80000000000000, 1.90000000000000, 2, 2.10000000000000, 2.20000000000000, 2.30000000000000, 2.40000000000000, 2.50000000000000, 2.60000000000000, 2.70000000000000, 2.80000000000000, 2.90000000000000, 3, 3.10000000000000, 3.20000000000000, 3.30000000000000, 3.40000000000000, 3.50000000000000, 3.60000000000000, 3.70000000000000, 3.80000000000000, 3.90000000000000, 4, 4.10000000000000, 4.20000000000000, 4.30000000000000, 4.40000000000000, 4.50000000000000, 4.60000000000000, 4.70000000000000, 4.80000000000000, 4.90000000000000, 5, 5.10000000000000, 5.20000000000000, 5.30000000000000, 5.40000000000000, 5.50000000000000, 5.60000000000000, 5.70000000000000, 5.80000000000000, 5.90000000000000, 6, 6.10000000000000, 6.20000000000000, 6.30000000000000, 6.40000000000000, 6.50000000000000, 6.60000000000000, 6.70000000000000]
fer_ldpc_16_10_classic = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.993333333333333, 1, 0.979090909090909, 0.962380952380952, 0.952380952380952, 0.869565217391304, 1, 0.909090909090909, 0.952380952380952, 1, 0.869565217391304, 0.952380952380952, 0.909090909090909, 0.952380952380952, 1, 0.833333333333333, 0.909090909090909, 0.833333333333333, 0.909090909090909, 0.909090909090909, 1, 0.869565217391304, 0.952380952380952, 0.800000000000000, 0.833333333333333, 0.740740740740741, 0.769230769230769, 0.740740740740741, 0.952380952380952, 0.833333333333333, 0.800000000000000, 0.800000000000000, 0.740740740740741, 0.869565217391304, 0.769230769230769, 0.952380952380952, 0.740740740740741, 0.625000000000000, 0.555555555555556, 0.526315789473684, 0.588235294117647, 0.645161290322581, 0.454545454545455, 0.487804878048781, 0.526315789473684, 0.434782608695652, 0.465116279069767, 0.416666666666667, 0.714285714285714, 0.344827586206897, 0.357142857142857, 0.400000000000000, 0.285714285714286, 0.322580645161290, 0.370370370370370, 0.266666666666667, 0.235294117647059, 0.322580645161290, 0.263157894736842, 0.219780219780220, 0.263157894736842, 0.194174757281553, 0.198019801980198, 0.161290322580645, 0.243902439024390, 0.165289256198347, 0.121951219512195, 0.103626943005181, 0.105820105820106, 0.127388535031847, 0.109289617486339, 0.0956937799043062, 0.132450331125828, 0.0869565217391304, 0.0851063829787234, 0.0985221674876847, 0.0836820083682008, 0.0576368876080692, 0.0619195046439629, 0.0518134715025907, 0.0475059382422803, 0.0617283950617284, 0.0365630712979890, 0.0298062593144560, 0.0284090909090909, 0.0242718446601942, 0.0372439478584730, 0.0241254523522316, 0.0180018001800180, 0.0131406044678055, 0.0129785853341986, 0.0125313283208020, 0.0125786163522013, 0.0117164616285882, 0.0169779286926995, 0.00729660707770887, 0.0110436223081171, 0.00594530321046373, 0.00695652173913044, 0.00573230151905990, 0.00913659205116492, 0.00453206435531385, 0.00321905681635281, 0.00247005063603804, 0.00214132762312634, 0.00437254044599913, 0.00201085863663784, 0.00283889283179560, 0.00173958423936679, 0.00108926529056152, 0.00108559952233621, 0.000996412913511359, 0.000650026001040042, 0.000821557673348669, 0.000416553849998959, 0.000404440759539746, 0.000288475407471513, 0.000314035831488373, 0.000255007713983348, 0.000237619998099040, 0.000193894269454866, 0.000151257326526754, 0.000138225597998493, 0.000116274337672304, 8.27444685322786e-05]
ber_ldpc_16_10_classic = [0.246875000000000, 0.241071428571429, 0.256250000000000, 0.238095238095238, 0.240625000000000, 0.271875000000000, 0.262500000000000, 0.247159090909091, 0.237500000000000, 0.253125000000000, 0.241071428571429, 0.256250000000000, 0.278125000000000, 0.205729166666667, 0.268750000000000, 0.198863636363636, 0.287500000000000, 0.199404761904762, 0.179347826086957, 0.209375000000000, 0.173295454545455, 0.193452380952381, 0.200000000000000, 0.206521739130435, 0.193452380952381, 0.178977272727273, 0.187500000000000, 0.193750000000000, 0.195312500000000, 0.167613636363636, 0.138020833333333, 0.167613636363636, 0.227272727272727, 0.206250000000000, 0.138586956521739, 0.160714285714286, 0.137500000000000, 0.153645833333333, 0.127314814814815, 0.134615384615385, 0.134259259259259, 0.172619047619048, 0.161458333333333, 0.130000000000000, 0.112500000000000, 0.136574074074074, 0.146739130434783, 0.144230769230769, 0.172619047619048, 0.134259259259259, 0.0996093750000000, 0.0989583333333333, 0.0970394736842105, 0.0919117647058824, 0.0846774193548387, 0.0681818181818182, 0.0899390243902439, 0.0838815789473684, 0.0733695652173913, 0.0843023255813954, 0.0755208333333333, 0.111607142857143, 0.0614224137931034, 0.0580357142857143, 0.0625000000000000, 0.0428571428571429, 0.0463709677419355, 0.0613425925925926, 0.0433333333333333, 0.0338235294117647, 0.0463709677419355, 0.0411184210526316, 0.0357142857142857, 0.0476973684210526, 0.0260922330097087, 0.0303217821782178, 0.0262096774193548, 0.0327743902439024, 0.0284090909090909, 0.0228658536585366, 0.0171632124352332, 0.0158730158730159, 0.0218949044585987, 0.0160519125683060, 0.0119617224880383, 0.0215231788079470, 0.0152173913043478, 0.0130319148936170, 0.0135467980295567, 0.0125523012552301, 0.00828530259365994, 0.0104489164086687, 0.00825777202072539, 0.00653206650831354, 0.00887345679012346, 0.00479890310786106, 0.00428464977645306, 0.00417258522727273, 0.00402002427184466, 0.00488826815642458, 0.00361881785283474, 0.00253150315031503, 0.00201215505913272, 0.00206846203763790, 0.00191885964912281, 0.00184748427672956, 0.00201376684241359, 0.00233446519524618, 0.00107168916453849, 0.00148398674765323, 0.000798900118906064, 0.00104347826086957, 0.000859845227858985, 0.00114207400639561, 0.000665646952186721, 0.000462739417350716, 0.000401383228356181, 0.000247591006423983, 0.000642216878006122, 0.000257641262819224, 0.000425833924769340, 0.000244629033660955, 0.000153177931485213, 0.000132307441784726, 0.000155689517736150, 9.95352314092564e-05, 0.000100127341439369, 7.02934621873243e-05, 6.57216234252088e-05, 4.14683398240300e-05, 4.61240127498548e-05, 3.34697624603144e-05, 2.67322497861420e-05, 3.09018991943693e-05, 1.79618075250520e-05, 2.07338396997740e-05, 1.41709349038121e-05, 1.34459761364953e-05]

esno_ldpc_16_10_bcjr = [-6, -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 6.4]
fer_ldpc_16_10_bcjr = [1.0, 1.0, 1.0, 1.0, 0.9090909090909091, 1.0, 0.7142857142857143, 0.7692307692307693, 0.9090909090909091, 1.0, 0.7692307692307693, 0.8333333333333334, 0.7692307692307693, 0.7692307692307693, 0.8333333333333334, 0.37037037037037035, 0.7692307692307693, 0.45454545454545453, 0.3225806451612903, 0.25, 0.3225806451612903, 0.30303030303030304, 0.2857142857142857, 0.17543859649122806, 0.07575757575757576, 0.04484304932735426, 0.062111801242236024, 0.021231422505307854, 0.019417475728155338, 0.01422475106685633, 0.009074410163339383, 0.010040160642570281]
ber_ldpc_16_10_bcjr = [0.0625, 0.0625, 0.0625, 0.0625, 0.056818181818181816, 0.0625, 0.044642857142857144, 0.04807692307692308, 0.056818181818181816, 0.0625, 0.04807692307692308, 0.052083333333333336, 0.04807692307692308, 0.04807692307692308, 0.052083333333333336, 0.023148148148148147, 0.04807692307692308, 0.028409090909090908, 0.020161290322580645, 0.015625, 0.020161290322580645, 0.01893939393939394, 0.017857142857142856, 0.010964912280701754, 0.004734848484848485, 0.002802690582959641, 0.0038819875776397515, 0.001326963906581741, 0.0012135922330097086, 0.0008890469416785206, 0.0005671506352087115, 0.0006275100401606426]

esno_ldpc_16_10_bcjr_v2 = [-4, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, -0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6]
fer_ldpc_16_10_bcjr_v2 = [0.8333333333333334, 0.8333333333333334, 0.8333333333333334, 1.0, 0.8333333333333334, 0.7142857142857143, 1.0, 0.625, 0.7142857142857143, 0.8333333333333334, 0.2777777777777778, 0.8333333333333334, 0.5, 0.3333333333333333, 0.2631578947368421, 0.21739130434782608, 0.1724137931034483, 0.14285714285714285, 0.18518518518518517, 0.08771929824561403, 0.05263157894736842, 0.00847457627118644, 0.006738544474393531, 0.0035739814152966403, 0.0003619254433586681]
ber_ldpc_16_10_bcjr_v2 = [0.052083333333333336, 0.052083333333333336, 0.052083333333333336, 0.0625, 0.052083333333333336, 0.044642857142857144, 0.0625, 0.0390625, 0.044642857142857144, 0.052083333333333336, 0.017361111111111112, 0.052083333333333336, 0.03125, 0.020833333333333332, 0.01644736842105263, 0.01358695652173913, 0.010775862068965518, 0.008928571428571428, 0.011574074074074073, 0.005482456140350877, 0.003289473684210526, 0.005296610169491525, 0.0004211590296495957, 0.00022337383845604002, 2.2620340209916757e-05]

esno_ldpc_16_10_bcjr_v2 = [-4.9, -4.5, -4.1, -3.7, -3.3, -2.9, -2.5, -2.1, -1.7, -1.3, -0.9, -0.5, -0.1, 0.3, 0.7, 1.1, 1.5, 1.9, 2.3, 2.7, 3.1, 3.5, 3.9, 4.3, 4.7, 5.1, 5.5, 5.9]
fer_ldpc_16_10_bcjr_v2 = [1.0, 1.0, 1.0, 1.0, 0.8333333333333334, 0.8333333333333334, 0.8333333333333334, 0.8333333333333334, 0.8333333333333334, 0.7142857142857143, 0.8333333333333334, 0.8333333333333334, 0.5, 0.7142857142857143, 0.4166666666666667, 0.3125, 0.2631578947368421, 0.4166666666666667, 0.1724137931034483, 0.20833333333333334, 0.09259259259259259, 0.07747126436781609, 0.06329113924050633, 0.050240963855421686, 0.012468827930174564, 0.004240882103477523, 0.0006351626016260162, 0.0001050626016260143]
ber_ldpc_16_10_bcjr_v2 = [0.0625, 0.0625, 0.0625, 0.0625, 0.052083333333333336, 0.052083333333333336, 0.052083333333333336, 0.052083333333333336, 0.052083333333333336, 0.044642857142857144, 0.052083333333333336, 0.052083333333333336, 0.03125, 0.044642857142857144, 0.026041666666666668, 0.01953125, 0.01644736842105263, 0.026041666666666668, 0.010775862068965518, 0.013020833333333334, 0.005787037037037037, 0.0035919540229885057, 0.003955696202531646, 0.0037650602409638554, 0.0007793017456359103, 0.0002650551314673452, 3.9697662601626015e-05, 1.0697662601626015e-05]



esno_bch_15_7_super_iter_2 = [-6, -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2]
fer_bch_15_7_super_iter_2 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.7692307692307693, 0.9090909090909091, 0.625, 0.8333333333333334, 0.7142857142857143, 0.5263157894736842, 0.37037037037037035, 0.4, 0.35714285714285715, 0.43478260869565216, 0.2777777777777778, 0.18181818181818182, 0.11363636363636363, 0.12658227848101267, 0.10752688172043011, 0.04608294930875576, 0.036101083032490974, 0.029850746268656716, 0.02079002079002079, 0.013422818791946308, 0.007455746268656716]

esno_bch_15_7_super_iter_5 = [-6, -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2]
fer_bch_15_7_super_iter_5 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.9375, 0.7142857142857143, 0.6818181818181818, 0.5357142857142857, 0.6818181818181818, 0.28846153846153844, 0.5, 0.26785714285714285, 0.2542372881355932, 0.20833333333333334, 0.15789473684210525, 0.1595744680851064, 0.12195121951219512, 0.04807692307692308, 0.06437768240343347, 0.04, 0.02094972067039106, 0.018867924528301886, 0.008710801393728223, 0.005048805116122518, 0.0027752081406105457, 0.001010801393728513]

esno_bch_15_7_super_iter_8 = [-6, -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8]
fer_bch_15_7_super_iter_8 = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5714285714285714, 0.5714285714285714, 0.27586206896551724, 0.5, 0.21621621621621623, 0.47058823529411764, 0.32, 0.18181818181818182, 0.18181818181818182, 0.1951219512195122, 0.1568627450980392, 0.125, 0.061068702290076333, 0.10526315789473684, 0.044444444444444446, 0.023323615160349854, 0.013582342954159592, 0.008456659619450317, 0.006477732793522267, 0.002191780821917808]

esno_bch_15_7 = [-6, -5.6, -5.2, -4.8, -4.4, -4.0, -3.6, -3.2, -2.8, -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2]
fer_bch_15_7 = [0.967741935483871, 0.967741935483871, 0.967741935483871, 0.8823529411764706, 0.8571428571428571, 0.7142857142857143, 0.8823529411764706, 0.6976744186046512, 0.6, 0.625, 0.410958904109589, 0.46153846153846156, 0.33707865168539325, 0.37037037037037035, 0.234375, 0.18292682926829268, 0.16304347826086957, 0.10600706713780919, 0.08310249307479224, 0.06465517241379311, 0.06593406593406594, 0.03787878787878788, 0.017899761336515514, 0.01430615164520744, 0.007234145165179648, 0.003906758692538091, 0.0018176310209027568, 0.0011545566502463055, 0.0005786811203927681]


# plt.plot(esno_bch_15_7_super_iter_5, fer_bch_15_7_super_iter_5, label="super iter 5", alpha=0.5, linewidth=2)
plt.plot(esno_bch_15_7_super_iter_2, gaussian_filter1d(fer_bch_15_7_super_iter_2, sigma=0.6).tolist(), label="super iter=2", alpha=1, linewidth=2)
plt.plot(esno_bch_15_7_super_iter_5, gaussian_filter1d(fer_bch_15_7_super_iter_5, sigma=0.6).tolist(), label="super iter=5", alpha=1, linewidth=2)
# plt.plot(esno_bch_15_7, fer_bch_15_7, label="classic bcjr", alpha=0.5, linewidth=2)
plt.plot(esno_bch_15_7, gaussian_filter1d(fer_bch_15_7, sigma=0.6).tolist(), label="bcjr", alpha=1, linewidth=2)


plt.yscale("log")  # Логарифмическая шкала по Y
plt.xlabel("EsNo")
plt.ylabel("FER")
plt.title("BCH(15, 7)")
plt.legend()
plt.grid(True, which="both", linestyle="--")
plt.show()
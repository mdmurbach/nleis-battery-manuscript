import os
import numpy as np
import math
from scipy.optimize import curve_fit


def phase_correction(i, v, cyc):
    v_shifted = np.zeros(v.size, dtype=np.complex128)
    i_shifted = np.zeros(v.size, dtype=np.complex128)
    index = cyc
    N = len(v)
    i_phase_reference = 0
    phase_correction = i_phase_reference - np.arctan2(np.imag(i[index]),
                                                      np.real(i[index]))
    for k in range(1, int(N/2)):
        v_phase_corrected = math.fmod((k/float(index))*phase_correction +
                                      np.arctan2(np.imag(v[k]), np.real(v[k])),
                                      2*np.pi)

        v_shifted[k] = np.abs(v[k])*(np.cos(v_phase_corrected) +
                                     1j*np.sin(v_phase_corrected))
        v_shifted[N - k] = np.real(v_shifted[k]) - 1j*np.imag(v_shifted[k])

        i_phase_corrected = (k/float(index))*phase_correction + \
            np.arctan2(np.imag(i[k]), np.real(i[k]))

        i_shifted[k] = np.abs(i[k])*(np.cos(i_phase_corrected) +
                                     1j*np.sin(i_phase_corrected))
        i_shifted[N - k] = np.real(i_shifted[k]) - 1j*np.imag(i_shifted[k])

    return i_shifted, v_shifted


def fit_amp_dependence(I, V, fit_type):
    def fit_Z1(i1, Z1):
        return Z1*i1

    def fit_Z2(i1, Z2):
        return Z2*i1**2

    if fit_type == 'linear':
        fit_real, tmp = curve_fit(fit_Z1, I, np.real(V))
        fit_imag, tmp = curve_fit(fit_Z1, I, np.imag(V))
    elif fit_type == 'quadratic':
        fit_real, tmp = curve_fit(fit_Z2, I, np.real(V))
        fit_imag, tmp = curve_fit(fit_Z2, I, np.imag(V))
    else:
        raise ValueError('{} is not an allowable fit_type, use "linear"' +
                         'or "quadratic"'.format(fit_type))

    return fit_real[0] + 1j*fit_imag[0]


def get_impedances(directories, SoCs, correction=None):
    Z1s, Z2s = [], []
    for directory in directories:
        d = {}

        for test in os.listdir(directory):
            if 'SOC' in test and not test.startswith('OCP'):
                SoC = test.split('SOC')[0]
                if SoC not in d.keys():
                    d[SoC] = {}

                data = np.genfromtxt(directory + test, delimiter=',',
                                     skip_header=1)
                frequency, amplitude = data[0, 3], data[0, 4]
                time, current, voltage = data[:, 0], data[:, 1], data[:, 2]

                if frequency not in d[SoC].keys():
                    d[SoC][frequency] = {}
                if amplitude not in d[SoC][frequency].keys():
                    d[SoC][frequency][amplitude] = {}

                cyc = (time[-1] + time[1])*frequency
                index = int(round(cyc))
                N = len(voltage)
                M = (np.concatenate([np.arange(0, N/2),
                                     np.arange(-N/2, 0)])/(cyc)).tolist()

                i_hat = np.fft.fft(current)/(N/2)  # Normalize FFT results
                v_hat = np.fft.fft(voltage)/(N/2)  # using N/2 (in mA/mV)

                i_hat_shft, v_hat_shft = phase_correction(i_hat, v_hat, index)
                i_baselined = fft_baseline_subt(M, i_hat_shft, index)
                v_baselined = fft_baseline_subt(M, v_hat_shft, index)

                d[SoC][frequency][amplitude]['i1'] = i_baselined[0]
                d[SoC][frequency][amplitude]['i2'] = i_baselined[1]
                d[SoC][frequency][amplitude]['i3'] = i_baselined[2]

                d[SoC][frequency][amplitude]['v1'] = v_baselined[0]
                d[SoC][frequency][amplitude]['v2'] = v_baselined[1]
                d[SoC][frequency][amplitude]['v3'] = v_baselined[2]

        for SoC in SoCs:
            if SoC in d.keys():
                Z1 = []
                Z2 = []

                frequencies = sorted(d[SoC].keys())
                for frequency in frequencies:
                    I1s = []
                    V1s = []
                    V2s = []
                    for amplitude in d[SoC][frequency].keys():
                        I1s.append(np.real(d[SoC][frequency][amplitude]['i1']))
                        V1s.append(d[SoC][frequency][amplitude]['v1'])
                        V2s.append(d[SoC][frequency][amplitude]['v2'])

                    Z1.append(fit_amp_dependence(I1s, V1s, 'linear'))

                    if correction is not None:
                        Z2.append(fit_amp_dependence(I1s, V2s, 'quadratic') -
                                  correction)
                    else:
                        Z2.append(fit_amp_dependence(I1s, V2s, 'quadratic'))

                Z1s.append(np.array(Z1))
                Z2s.append(np.array(Z2))
            else:
                print('{} is not in the dataset'.format(SoC))

    return frequencies, Z1s, Z2s


def fft_baseline_subt(harmonics_idx, x_hat_shift, index):
    k_idx = np.array(harmonics_idx[0:4*index])
    x_idx = np.array(x_hat_shift[0:4*index])
    k_ids = np.array([np.argwhere(np.isclose(k_idx, 1.0, rtol=5e-4)),
                      np.argwhere(np.isclose(k_idx, 2.0, rtol=5e-4)),
                      np.argwhere(np.isclose(k_idx, 3.0, rtol=5e-4))]).flatten()  # noqa E501

    x = np.zeros(3, dtype='complex128')
    count = 0
    for n in k_ids:
        k_baseline = np.concatenate((np.array(k_idx[n-5:n-1]),
                                     np.array(k_idx[n+2:n+6])), axis=0)
        xr_baseline = np.concatenate((np.real(np.array(x_idx[n-5:n-1])),
                                      np.real(np.array(x_idx[n+2:n+6]))),
                                     axis=0)
        xi_baseline = np.concatenate((np.imag(np.array(x_idx[n-5:n-1])),
                                      np.imag(np.array(x_idx[n+2:n+6]))),
                                     axis=0)

        # Quadratic Fits
        baseline1_real = np.polyfit(k_baseline, xr_baseline, 2)
        baseline1_imag = np.polyfit(k_baseline, xi_baseline, 2)

        x[count] = x_idx[n] - np.complex(np.polyval(baseline1_real, k_idx[n]),
                                         np.polyval(baseline1_imag, k_idx[n]))
        count += 1

    return x

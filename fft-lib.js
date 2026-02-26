/**
 * FFT Library - Browser compatible version of fft-js
 * Based on: https://github.com/vail-systems/node-fft
 * Using Cooley-Tukey algorithm
 */

const FFTLib = (function() {
    'use strict';

    // Complex number operations
    const complex = {
        add: function(a, b) {
            return [a[0] + b[0], a[1] + b[1]];
        },
        subtract: function(a, b) {
            return [a[0] - b[0], a[1] - b[1]];
        },
        multiply: function(a, b) {
            return [(a[0] * b[0] - a[1] * b[1]), (a[0] * b[1] + a[1] * b[0])];
        },
        magnitude: function(c) {
            return Math.sqrt(c[0] * c[0] + c[1] * c[1]);
        }
    };

    // Calculate exponent for FFT
    function exponent(k, N) {
        var x = -2 * Math.PI * (k / N);
        return [Math.cos(x), Math.sin(x)];
    }

    // FFT implementation using Cooley-Tukey algorithm
    function fft(vector) {
        var X = [],
            N = vector.length;

        // Base case
        if (N === 1) return [[vector[0], 0]];

        // Split into even and odd indices
        var evens = [], odds = [];
        for (var i = 0; i < N; i++) {
            if (i % 2 === 0) evens.push(vector[i]);
            else odds.push(vector[i]);
        }

        // Recurse
        var X_evens = fft(evens);
        var X_odds = fft(odds);

        // Combine
        for (var k = 0; k < N / 2; k++) {
            var t = X_evens[k];
            var e = complex.multiply(exponent(k, N), X_odds[k]);

            X[k] = complex.add(t, e);
            X[k + (N / 2)] = complex.subtract(t, e);
        }

        return X;
    }

    // Calculate FFT magnitudes
    function fftMag(fftBins) {
        return fftBins.map(complex.magnitude);
    }

    // Calculate frequency bins
    function fftFreq(fftBins, sampleRate) {
        var stepFreq = sampleRate / (fftBins.length * 2);
        return fftBins.map(function(__, ix) {
            return ix * stepFreq;
        });
    }

    // Apply window function to signal
    function applyWindow(signal, windowType) {
        windowType = windowType || 'hanning';
        var N = signal.length;
        var windowed = new Array(N);

        for (var i = 0; i < N; i++) {
            var w = 1;
            if (windowType === 'hanning') {
                // Hanning window: 0.5 * (1 - cos(2*pi*n/N))
                w = 0.5 * (1 - Math.cos(2 * Math.PI * i / N));
            } else if (windowType === 'hamming') {
                // Hamming window: 0.54 - 0.46 * cos(2*pi*n/N)
                w = 0.54 - 0.46 * Math.cos(2 * Math.PI * i / N);
            } else if (windowType === 'blackman') {
                // Blackman window
                w = 0.42 - 0.5 * Math.cos(2 * Math.PI * i / N) + 0.08 * Math.cos(4 * Math.PI * i / N);
            }
            windowed[i] = signal[i] * w;
        }

        return windowed;
    }

    // STFT (Short-Time Fourier Transform)
    function stft(signal, fftSize, hopSize, windowType) {
        fftSize = fftSize || 2048;
        hopSize = hopSize || 512;
        windowType = windowType || 'hanning';

        var numFrames = Math.floor((signal.length - fftSize) / hopSize) + 1;
        var spectrogram = [];

        for (var i = 0; i < numFrames; i++) {
            var start = i * hopSize;
            var frame = signal.slice(start, start + fftSize);

            // Apply window
            var windowed = applyWindow(frame, windowType);

            // Pad to power of 2 if needed
            var padded = padToPowerOf2(windowed);

            // Compute FFT
            var fftResult = fft(padded);

            // Get magnitudes (only first half - positive frequencies)
            var mags = fftMag(fftResult.slice(0, fftResult.length / 2));

            spectrogram.push(mags);
        }

        return spectrogram;
    }

    // Pad array to power of 2
    function padToPowerOf2(arr) {
        var N = arr.length;
        var nextPow2 = Math.pow(2, Math.ceil(Math.log2(N)));

        if (N === nextPow2) return arr;

        var padded = new Array(nextPow2);
        for (var i = 0; i < nextPow2; i++) {
            padded[i] = (i < N) ? arr[i] : 0;
        }
        return padded;
    }

    // Public API
    return {
        fft: fft,
        fftMag: fftMag,
        fftFreq: fftFreq,
        stft: stft,
        applyWindow: applyWindow
    };
})();

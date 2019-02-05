clc; close all; clear all;
[data,Fs] = audioread('Sound files/LDC93S1.wav');
plot(data)
sound(data)
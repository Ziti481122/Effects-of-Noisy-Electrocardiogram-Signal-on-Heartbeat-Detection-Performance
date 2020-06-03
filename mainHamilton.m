
%main hamilton algorithm
%input: signal, fs
%output:amplitude and index of QRS

clear all;
clc;

load sig;
fs = 360;

val = sig;
[amp,ind]=hamilton_algo(val,fs);
csvwrite('hamiltonQRS.csv',ind(:));

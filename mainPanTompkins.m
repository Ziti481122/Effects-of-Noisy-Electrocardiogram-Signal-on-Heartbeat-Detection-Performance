
%main pan tompkins algorithm
%input: signal, fs
%output:amplitude and index of QRS

clear all;
clc;

load sig;
fs = 360;

val = sig;
[amp,ind]=pan_tompkins_algo(val,fs);
csvwrite('pantompkinsQRS.csv',ind(:));

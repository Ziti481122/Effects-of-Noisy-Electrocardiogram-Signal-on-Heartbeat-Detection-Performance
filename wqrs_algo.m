
%wqrs algorithm
%input: signal, fs
%output:amplitude and index of QRS
%requirement= wfdb toolbox from physionet website
%
%Notes: 
%The WQRS algorithm is written in C and included in the WFDB Software Package 
%and can be downloaded as source code (C/C++) from PhysioNet website.
%The WFDB Toolbox provides access from MATLAB to WQRS applications 
%in the WFDB Software Package.
%This file is used to access the WQRS algorithm.

clear all;
clc;

load sig
%fs = 360;
sigV = [sig];
mat2wfdb(sigV,'Ex1',360,[],[]);
wqrs('Ex1',[],[],[]);
[ann]=rdann('Ex1','wqrs');
csvwrite('wqrsQRS.csv',ann)


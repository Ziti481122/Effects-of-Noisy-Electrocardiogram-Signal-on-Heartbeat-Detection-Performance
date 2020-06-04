
%wqrs algorithm
%input: signal, fs
%output:amplitude and index of QRS
%requirement= wfdb toolbox from physionet website
%The WFDB Toolbox provides access from MATLAB to WQRS applications 
%in the WFDB Software Package in Physionet website.

clear all;
clc;

load sig
%fs = 360;
sigV = [sig];
mat2wfdb(sigV,'Ex1',360,[],[]);
wqrs('Ex1',[],[],[]);
[ann]=rdann('Ex1','wqrs');
csvwrite('wqrsQRS.csv',ann)


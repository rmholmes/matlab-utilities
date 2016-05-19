function [PSD,FRQ,PSD_C] = WelchPSD(x,y,varargin)
%
% This is a wrapper for the welch psd estimator that sets options
% that work well.
%
% x = independent variable (if length 1, then assumed dt = x)
%
% y = dependent variable, where the dimension to estimate the psd
% is the final dimension.
%
% varargin : can set the
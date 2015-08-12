function [Params] = fitThermoModel(tf, sites, expr, act, rep)
%Use simulated annealing to fit the thermodynamic model based on Reinitz'
%papers to a given set of expression data from a set of binding sites and
%transcription factor concentrations
%INPUTS:
%tf   - dictionary of protein concentrations, keys are identifier strings
%sites- struct with fields "Name" as in tf, "m", start base, "n", end base,
%       and "Score" for pwm score by that particular TF
%expr - a tx100 matrix of observed expression values
%
%OUTPUTS:
%Params - struct with fields "K", a dictionary of binding affinities, "E",
%       a dictionary of effectivenesses (A or Q), "R_max", for maximum
%       transcription rate and "THETA" for activation energy of 
%       transcription initiation 
Params = struct('E', [], 'K', [], 'R_max', [], 'THETA', []);

%Pick set of parameters based on temperature and old set of parameters

%test parameters
%stuff = useThermoModel(stuff)

%compare test with experimental data and use some metric to determine
%similarity of prediction and experiment

%compare similarity of newly chosen parameters to old parameters 

%either pick new set of parameters or don't


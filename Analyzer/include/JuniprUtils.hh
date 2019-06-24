#include "fastjet/ClusterSequence.hh"
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <algorithm>
using namespace fastjet;
using namespace std;

vector<double> PJ_3mom(PseudoJet PJ);
double dot3(vector<double> v1, vector<double> v2);
vector<double> hat(vector<double> v);
double get_theta(vector<double> v1, vector<double> v2);
double get_thetaPJ(PseudoJet PJ1, PseudoJet PJ2);
vector<double> cross(vector<double> v1, vector<double> v2);
double get_phiPJ(PseudoJet p, PseudoJet ref);
vector<double> PJ_to_ETPM(PseudoJet PJ, PseudoJet ref);
